//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file jet_blast.cpp
//! \brief Problem generator for a jet blast wave inside a star.  Works in Cartesian,
//!        cylindrical, and spherical coordinates.  Contains both non-relativistic and relativistic
//!        implementations. 
//!
//! REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
//!   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <vector>

struct PolytropeData {
  std::vector<double> r, rho, m, P;
};

// Load a CSV with header "r,rho,m,P"
static PolytropeData LoadPolytropeCSV(const std::string &filename) {
  std::ifstream in(filename);
  if (!in.is_open())
    throw std::runtime_error("Failed to open " + filename);
  std::string line;
  // skip header
  if (!std::getline(in, line))
    throw std::runtime_error("Empty file: " + filename);
  PolytropeData out;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    double rv, rv_rho, rv_m, rv_P;
    char comma;
    if (!(ss >> rv >> comma >> rv_rho >> comma >> rv_m >> comma >> rv_P))
      throw std::runtime_error("Parse error: " + line);
    out.r.push_back(rv);
    out.rho.push_back(rv_rho);
    out.m.push_back(rv_m);
    out.P.push_back(rv_P);
  }
  return out;
}

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

#include "../parameter_input.hpp"


int RefinementCondition(MeshBlock *pmb);

Real threshold;
// Optional AMR velocity masking / trigger: configured in InitUserMeshData
// amr_v_cut: below this speed, ignore rho/P gradients (treat as static star)
// amr_v_ref: reference speed to scale the velocity-based indicator (eta_v = |v|/amr_v_ref)
static Real amr_v_cut = 0.0;
static Real amr_v_ref = 0.0;

// --- breakout tracking (configured in InitUserMeshData) ---
static Real breakout_rout, breakout_x1_0, breakout_x2_0, breakout_x3_0;
static int  breakout_nbins;
static Real breakout_factor, breakout_vmin;
// If <=0, we auto-compute ~2 * finest dx on first call
static Real breakout_ringw;
static bool breakout_params_inited = false;
static Real breakout_phi0_deg;
// --- jet driving parameters (configured in InitUserMeshData) ---
static Real jet_t_stop = 0.0;      // stop time for driving
static Real jet_rinj   = 0.0;      // injection radius (nozzle size)
static Real jet_Gam    = 1.0;      // Lorentz factor of injected flow
static Real jet_rho    = 0.0;      // comoving rest-mass density in jet
static Real jet_p      = 0.0;      // gas pressure in jet
static Real gate_dir_min = 0.0;    // angle gate min (same meaning as problem/dir_angle_min)
static Real gate_dir_max = M_PI;   // angle gate max (same as problem/dir_angle_max)
static Real gate_phi0   = 0.0;     // center direction for 2D Cartesian wedge (radians)
static bool gate_cone_bipolar = true; // also include opposite wedge
static bool jet_enabled = false;   // enable jet driving when inputs provided
// ----------------------------------------------------------


  
// ---- User history outputs: global integrals ----
static Real Hst_Etot(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += u(IEN,k,j,i) * pco->GetCellVolume(k,j,i); // Etot density * dV
  return sum;
}


static Real Hst_Eexcess(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += (u(IEN,k,j,i) - u(IDN,k,j,i)) * pco->GetCellVolume(k,j,i); // Etot - D = tau
  return sum;
}

// Integrate over only positive x and then reflect over x=0 plane

static Real Hst_Px_pos(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real Mx;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Mx = u(IM1,k,j,i);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          Real ph = pco->x2v(j);
          Real MR   = u(IM1,k,j,i);
          Real Mphi = u(IM2,k,j,i);
          Mx = MR*std::cos(ph) - Mphi*std::sin(ph);
        } else { // spherical_polar
          Real th = pco->x2v(j);
          Real ph = pco->x3v(k);
          Real Mr  = u(IM1,k,j,i);
          Real Mth = u(IM2,k,j,i);
          Real Mph = u(IM3,k,j,i);
          Real erx  = std::sin(th)*std::cos(ph);
          Real etx  = std::cos(th)*std::cos(ph);
          Real ephx = -std::sin(ph);
          Mx = Mr*erx + Mth*etx + Mph*ephx;
        }
        if (Mx > 0.0)
          sum += Mx * pco->GetCellVolume(k,j,i);
      }
  return sum;
}

static Real Hst_Px_neg(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real Mx;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Mx = u(IM1,k,j,i);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          Real ph = pco->x2v(j);
          Real MR   = u(IM1,k,j,i);
          Real Mphi = u(IM2,k,j,i);
          Mx = MR*std::cos(ph) - Mphi*std::sin(ph);
        } else { // spherical_polar
          Real th = pco->x2v(j);
          Real ph = pco->x3v(k);
          Real Mr  = u(IM1,k,j,i);
          Real Mth = u(IM2,k,j,i);
          Real Mph = u(IM3,k,j,i);
          Real erx  = std::sin(th)*std::cos(ph);
          Real etx  = std::cos(th)*std::cos(ph);
          Real ephx = -std::sin(ph);
          Mx = Mr*erx + Mth*etx + Mph*ephx;
        }
        if (Mx < 0.0)
          sum += Mx * pco->GetCellVolume(k,j,i); // remains negative
      }
  return sum;
}



static Real Hst_Py(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real My;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          My = u(IM2,k,j,i);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          // native basis: (R, phi, z). Transform to Cartesian y: My = MR*sin(phi) + Mphi*cos(phi)
          Real ph = pco->x2v(j);
          Real MR   = u(IM1,k,j,i);
          Real Mphi = u(IM2,k,j,i);
          My = MR*std::sin(ph) + Mphi*std::cos(ph);
        } else { // spherical_polar: native basis (r, theta, phi)
          // Cartesian y components of (e_r, e_theta, e_phi):
          // e_r_y =  sin(theta) sin(phi)
          // e_th_y = cos(theta) sin(phi)
          // e_ph_y =  cos(phi)
          Real th = pco->x2v(j);
          Real ph = pco->x3v(k);
          Real Mr  = u(IM1,k,j,i);
          Real Mth = u(IM2,k,j,i);
          Real Mph = u(IM3,k,j,i);
          Real ery  = std::sin(th)*std::sin(ph);
          Real ety  = std::cos(th)*std::sin(ph);
          Real ephy =  std::cos(ph);
          My = Mr*ery + Mth*ety + Mph*ephy;
        }
        sum += My * pco->GetCellVolume(k,j,i);
      }
  return sum;
}

static Real Hst_Pz(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real Mz;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Mz = u(IM3,k,j,i);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          // native z is Cartesian z in cylindrical coords
          Mz = u(IM3,k,j,i);
        } else { // spherical_polar
          // Cartesian z components: e_r_z = cos(theta), e_th_z = -sin(theta), e_ph_z = 0
          Real th = pco->x2v(j);
          Real Mr  = u(IM1,k,j,i);
          Real Mth = u(IM2,k,j,i);
          Mz = Mr*std::cos(th) - Mth*std::sin(th);
        }
        sum += Mz * pco->GetCellVolume(k,j,i);
      }
  return sum;
}

static Real Hst_PgasInt(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &w = pmb->phydro->w;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += w(IPR,k,j,i) * pco->GetCellVolume(k,j,i); //  p_gas*dV
  return sum;
}

// ---- User-defined per-cell outputs for VTK/HDF5 ----
// Constant (or cell-dependent, if EOS supports it) adiabatic index Γ
static Real Out_Gamma(MeshBlock *pmb, int iout, int k, int j, int i) {
  return pmb->peos->GetGamma();
}

// Gas enthalpy w_gas = rho + Γ/(Γ-1) * p_gas (SR/NR gas-only contribution)
static Real Out_wgas(MeshBlock *pmb, int iout, int k, int j, int i) {
  AthenaArray<Real> &w = pmb->phydro->w;
  const Real gamma = pmb->peos->GetGamma();
  const Real rho   = w(IDN,k,j,i);
  const Real pgas  = w(IPR,k,j,i);
  return rho + (gamma/(gamma - 1.0))*pgas;
}

// ---- History helpers (global integrals) ----
// ∫ Γ dV (use with V_tot to form volume-avg Γ = Gamma_int / V_tot)
static Real Hst_GammaInt(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  const Real gamma = pmb->peos->GetGamma();
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += gamma * pco->GetCellVolume(k,j,i);
  return sum;
}

// ∫ dV (total volume)
static Real Hst_Volume(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += pco->GetCellVolume(k,j,i);
  return sum;
}




void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem","thr");
    amr_v_cut = pin->GetOrAddReal("problem", "amr_v_cut", 0.0);
    amr_v_ref = pin->GetOrAddReal("problem", "amr_v_ref", 0.0);
  }


  // Breakout tracker inputs (available to UserWorkInLoop without a pin)
  breakout_rout   = pin->GetReal("problem", "star_radius");
  breakout_x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
  breakout_x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
  breakout_x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
  breakout_nbins  = pin->GetOrAddInteger("problem", "breakout_nbins", 90);   // bins over [0, 90 deg]
  breakout_factor = pin->GetOrAddReal("problem", "breakout_factor", 10.0);   // p > factor * pamb
  breakout_vmin   = pin->GetOrAddReal("problem", "breakout_vmin", 0.01);     // outward vr threshold
  breakout_ringw  = pin->GetOrAddReal("problem", "breakout_ring_width", -1.0); // auto if <= 0
  breakout_phi0_deg = pin->GetOrAddReal("problem", "breakout_phi0_deg", 0.0); // rotate ring coords
  // Jet driving inputs (optional)
  jet_t_stop = pin->GetOrAddReal("problem", "t_stop", 0.0);
  jet_rinj   = pin->GetOrAddReal("problem", "jet_rinj", 0.0);
  jet_Gam    = pin->GetOrAddReal("problem", "jet_Gam", 1.0);
  jet_rho    = pin->GetOrAddReal("problem", "jet_rho", 0.0);
  jet_p      = pin->GetOrAddReal("problem", "jet_p", 0.0);
  gate_dir_min = pin->GetOrAddReal("problem", "dir_angle_min", 0.0);
  gate_dir_max = pin->GetOrAddReal("problem", "dir_angle_max", M_PI);
  gate_phi0    = pin->GetOrAddReal("problem", "phi0", 0.0);
  gate_cone_bipolar = pin->GetOrAddBoolean("problem", "cone_bipolar", true);
  // Enable jet only if a positive stop time and radius are provided
  jet_enabled = (jet_t_stop > 0.0) && (jet_rinj > 0.0) && (jet_Gam >= 1.0);

  // Register global integrals in history output (allocate N slots, then enroll by index)
  AllocateUserHistoryOutput(9);
  EnrollUserHistoryOutput(0, Hst_Etot,     "E_tot");
  EnrollUserHistoryOutput(1, Hst_Eexcess,  "E_excess"); // Etot - D (tau)
  EnrollUserHistoryOutput(2, Hst_Px_pos,   "Px_pos");
  EnrollUserHistoryOutput(3, Hst_Px_neg,   "Px_neg");
  EnrollUserHistoryOutput(4, Hst_Py,       "Py_tot");
  EnrollUserHistoryOutput(5, Hst_Pz,       "Pz_tot");
  EnrollUserHistoryOutput(6, Hst_PgasInt,  "Pgas_int");
  EnrollUserHistoryOutput(7, Hst_GammaInt, "Gamma_int");
  EnrollUserHistoryOutput(8, Hst_Volume,   "V_tot");

  if (Globals::my_rank == 0) {
    std::fprintf(stderr,
      "[jet:init] enabled=%d t_stop=%g rinj=%g Gam=%g rho=%g p=%g gate=[%g,%g] phi0=%g bipolar=%d\n",
      (int)jet_enabled, (double)jet_t_stop, (double)jet_rinj, (double)jet_Gam,
      (double)jet_rho, (double)jet_p,
      (double)gate_dir_min, (double)gate_dir_max, (double)gate_phi0,
      (int)gate_cone_bipolar);
  }
  breakout_params_inited = true;
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rout =  pin->GetReal("problem", "star_radius");
  Real poly_idx = pin->GetReal("problem", "poly_index");
  std::string poly_csv_path = pin->GetString("problem", "poly_csv_path");
  Real cs_factor = pin->GetReal("problem","cs_factor");

  // Compile-time detector for dynamics mode
#if defined(RELATIVISTIC_DYNAMICS) && (RELATIVISTIC_DYNAMICS != 0)
  const char* dyn_mode = "SR";
#else
  const char* dyn_mode = "NR";
#endif

  // --- diagnostics ---
  long long tot_cells = 0;
  long long nan_rho = 0, nan_p = 0, nan_v2 = 0;
  Real rho_min=1e99, rho_max=-1e99, p_min=1e99, p_max=-1e99;
  // --- end diagnostics ---

    // Load precomputed polytrope data from CSV
  static PolytropeData poly;
  if (poly.r.empty()) {
    poly = LoadPolytropeCSV(poly_csv_path);
  }

  Real b0, angle;
  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = pin->GetReal("problem", "b0");
    angle = (PI/180.0)*pin->GetReal("problem", "angle");
  }
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;



  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }


  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        // Cartesian components of the cell center relative to origin (filled per coord system)
        Real x, y, z;
        Real rad;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          x = pcoord->x1v(i);
          y = pcoord->x2v(j);
          z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
          x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }

        Real rr     = rad;
        Real r_star = rout;
        Real rho, pgas;


        // Use precomputed Polytropic structure
        Real rho0, pgas0;
        {
          const auto &r_arr   = poly.r;
          const auto &rho_arr = poly.rho;
          const auto &P_arr   = poly.P;
          int N = static_cast<int>(r_arr.size());
          if (rr <= r_star) {
            auto it = std::lower_bound(r_arr.begin(), r_arr.end(), rr);
            int idx = static_cast<int>(std::distance(r_arr.begin(), it));
            if (idx <= 0) {
              rho0  = rho_arr[0];
              pgas0 = P_arr[0];
            } else if (idx >= N) {
              rho0  = rho_arr[N-1];
              pgas0 = P_arr[N-1];
            } else {
              double t = (rr - r_arr[idx-1]) / (r_arr[idx] - r_arr[idx-1]);
              rho0  = rho_arr[idx-1] + t * (rho_arr[idx] - rho_arr[idx-1]);
              pgas0 = P_arr[idx-1]   + t * (P_arr[idx]   - P_arr[idx-1]);
            }
          } else {
            rho0  = 1.0e-19;
            pgas0 = 1.0e-19;
          }
        }


        // Initialize velocity components
        Real vx = 0.0, vy = 0.0, vz = 0.0;

        // Get rid of Pressure
          pgas0 *= cs_factor;
          
        

        // Final cell primitives (BM-modified inside shell, star/ambient elsewhere)
        // floors to avoid NaNs in conserved mapping
        if (!std::isfinite(rho0) || rho0 <= 0.0) rho0 = 1e-30;
        if (!std::isfinite(pgas0) || pgas0 <= 0.0) pgas0 = 1e-30;
        rho  = rho0;
        pgas = pgas0;
        // diagnostics
        tot_cells++;
        if (!std::isfinite(rho)) nan_rho++;
        if (!std::isfinite(pgas)) nan_p++;
        if (rho < rho_min) rho_min = rho;
        if (rho > rho_max) rho_max = rho;
        if (pgas < p_min) p_min = pgas;
        if (pgas > p_max) p_max = pgas;

        // Write primitives
        phydro->w(IDN,k,j,i) = rho;     // comoving rest-mass density
        phydro->w(IPR,k,j,i) = pgas;    // gas pressure
        phydro->w(IVX,k,j,i) = vx;
        phydro->w(IVY,k,j,i) = vy;
        phydro->w(IVZ,k,j,i) = vz;

        // Compute conserved from primitives to keep arrays consistent
        #if defined(RELATIVISTIC_DYNAMICS) && (RELATIVISTIC_DYNAMICS != 0)
        // Special-relativistic conserved variables (c=1)
        Real v2      = vx*vx + vy*vy + vz*vz;
        v2 = std::min(v2, 1.0 - 1e-12);     // keep subluminal
        Real gL = 1.0/std::sqrt(1.0 - v2);
        Real rho_safe = std::max(rho, 1e-30);
        Real h_spec  = 1.0 + (gamma/(gamma - 1.0)) * (pgas / rho_safe); // ideal-gas EOS
        Real D       = rho_safe * gL;
        Real momx    = rho_safe * h_spec * gL*gL * vx;
        Real momy    = rho_safe * h_spec * gL*gL * vy;
        Real momz    = rho_safe * h_spec * gL*gL * vz;
        Real Etot    = rho_safe * h_spec * gL*gL - pgas;  // includes rest-mass energy

        phydro->u(IDN,k,j,i) = D;
        phydro->u(IM1,k,j,i) = momx;
        phydro->u(IM2,k,j,i) = momy;
        phydro->u(IM3,k,j,i) = momz;
        // store Etot (include rest mass)
        phydro->u(IEN,k,j,i) = Etot;
        // diagnostic: v^2 sanity
        if (!std::isfinite(v2)) nan_v2++;
        // catch obviously bad initial energies
        if (!std::isfinite(phydro->u(IEN,k,j,i)) || phydro->u(IEN,k,j,i) <= 0.0) {
          // keep a tiny positive energy to avoid inversion failure
          phydro->u(IEN,k,j,i) = std::max(phydro->u(IEN,k,j,i), (Real)1e-30); 
        }
        #else
        // Non-relativistic conserved variables
        Real v2 = vx*vx + vy*vy + vz*vz;  // for diagnostics
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho * vx;
        phydro->u(IM2,k,j,i) = rho * vy;
        phydro->u(IM3,k,j,i) = rho * vz;
        phydro->u(IEN,k,j,i) = pgas/gm1 + 0.5*rho*v2;
        // diagnostic: v^2 sanity
        if (!std::isfinite(v2)) nan_v2++;
        #endif

      }
    }
  }

  // diagnostics summary
  if (Globals::my_rank == 0) {
    std::fprintf(stderr,
      "[init] dyn=%s  rho[min,max]=[%g,%g]  p[min,max]=[%g,%g]  NaNs: rho=%lld p=%lld v2=%lld\n",
      dyn_mode, (double)rho_min, (double)rho_max, (double)p_min,
      (double)p_max, nan_rho, nan_p, nan_v2);
  }


  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            pfield->b.x1f(k,j,i) = b0 * std::cos(angle);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            Real phi = pcoord->x2v(j);
            pfield->b.x1f(k,j,i) =
                b0 * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x1f(k,j,i) = b0 * std::abs(std::sin(theta))
                                   * (std::cos(angle) * std::cos(phi)
                                      + std::sin(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            pfield->b.x2f(k,j,i) = b0 * std::sin(angle);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            Real phi = pcoord->x2v(j);
            pfield->b.x2f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x2f(k,j,i) = b0 * std::cos(theta)
                                   * (std::cos(angle) * std::cos(phi)
                                      + std::sin(angle) * std::sin(phi));
            if (std::sin(theta) < 0.0)
              pfield->b.x2f(k,j,i) *= -1.0;
          }
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0
              || std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            pfield->b.x3f(k,j,i) = 0.0;
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            Real phi = pcoord->x3v(k);
            pfield->b.x3f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IEN,k,j,i) += 0.5*b0*b0;
        }
      }
    }
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkInLoop()
//! \brief JET DRIVING inside cells add a radial velocity.
//!       Track shock breakout times at r≈star_radius for angles in [0, 90°] and log CSV.
//========================================================================================
void Mesh::UserWorkInLoop() {
  // Ensure parameters were cached
  if (!breakout_params_inited) return;

  // Determine a default ring width ~ 2 * local finest dx, lazily on first call
  static Real ringw_default = -1.0;
  if (ringw_default < 0.0 && nblocal > 0) {
    MeshBlock *pmb0 = my_blocks(0);
    Real dx1 = pmb0->pcoord->dx1f(pmb0->is);
    Real dx2 = pmb0->pcoord->dx2f(pmb0->js);
    Real dx3 = pmb0->pcoord->dx3f(pmb0->ks);
    ringw_default = 2.0 * std::max(dx1, std::max(dx2, dx3));
  }
  const Real rout   = breakout_rout;
  const Real x1_0   = breakout_x1_0;
  const Real x2_0   = breakout_x2_0;
  const Real x3_0   = breakout_x3_0;
  const int  nbins  = breakout_nbins;
  const Real factor = breakout_factor;
  const Real vmin   = breakout_vmin;
  const Real ringw  = (breakout_ringw > 0.0) ? breakout_ringw
                                             : (ringw_default > 0.0 ? ringw_default : 1e-2);

  // Compute breakout center in Cartesian coordinates for coordinate-agnostic math
  Real x0c, y0c, z0c;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0c = breakout_x1_0; y0c = breakout_x2_0; z0c = breakout_x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    Real R0 = breakout_x1_0, ph0 = breakout_x2_0, z0 = breakout_x3_0;
    x0c = R0*std::cos(ph0); y0c = R0*std::sin(ph0); z0c = z0;
  } else { // spherical_polar
    Real r0 = breakout_x1_0, th0 = breakout_x2_0, ph0 = breakout_x3_0;
    x0c = r0*std::sin(th0)*std::cos(ph0);
    y0c = r0*std::sin(th0)*std::sin(ph0);
    z0c = r0*std::cos(th0);
  }

  // Static state across calls: per-bin breakout time and ambient reference
  static bool inited = false;
  static std::vector<Real> t_break;     // breakout time per bin (or <0 if unknown)
  static std::vector<Real> pamb_bin;    // ambient reference pressure per bin (sampled at t=0)
  static std::vector<int>  logged;      // whether we've already written the line for a bin

  // Initialize on first entry or if nbins changed
  if (!inited || (int)t_break.size() != nbins) {
    t_break.assign(nbins, -1.0);
    pamb_bin.assign(nbins, -1.0);
    logged.assign(nbins, 0);
    inited = true;
    if (Globals::my_rank == 0) {
      // Create/overwrite CSV with a header
      FILE* f = std::fopen("shock_breakout.csv", "w");
      if (f) {
        std::fprintf(f, "angle_deg,radius,breakout_time\n");
        std::fclose(f);
      }
    }
  }

  // Helper to convert cell center to Cartesian
  auto cell_to_cart = [&](MeshBlock* pmb, int k, int j, int i, Real &x, Real &y, Real &z){
    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
      x = pmb->pcoord->x1v(i);
      y = pmb->pcoord->x2v(j);
      z = pmb->pcoord->x3v(k);
    } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
      Real R = pmb->pcoord->x1v(i);
      Real ph = pmb->pcoord->x2v(j);
      x = R*std::cos(ph); y = R*std::sin(ph); z = pmb->pcoord->x3v(k);
    } else { // spherical_polar
      Real r  = pmb->pcoord->x1v(i);
      Real th = pmb->pcoord->x2v(j);
      Real ph = pmb->pcoord->x3v(k);
      x = r*std::sin(th)*std::cos(ph);
      y = r*std::sin(th)*std::sin(ph);
      z = r*std::cos(th);
    }
  };

  // --- Jet driving: enforce jet state inside a small nozzle while time <= t_stop ---
  if (jet_enabled && (time <= jet_t_stop)) {
    // ---- debug counters (per-step, local to this rank) ----
    int dbg_in_rad = 0;        // cells with rad <= jet_rinj (geometric nozzle)
    int dbg_angle_ok = 0;      // cells that also pass the angle/wedge gate
    int dbg_written = 0;       // cells actually overwritten with jet state

    for (int nb=0; nb<nblocal; ++nb) {
      MeshBlock* pmb = my_blocks(nb);
      auto &coord = *pmb->pcoord;
      AthenaArray<Real> &w = pmb->phydro->w;
      AthenaArray<Real> &u = pmb->phydro->u;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            // Cartesian position of cell center (independent of native mesh coordinates)
            Real x, y, z; cell_to_cart(pmb, k, j, i, x, y, z);
            Real dx = x - x0c, dy = y - y0c, dz = z - z0c;
            bool is2d_local = (pmb->ks == pmb->ke);
            // In 2D runs (single zone in x3), treat the nozzle radius as in-plane (ignore dz)
            if (is2d_local) {
              dz = 0.0;           // collapse the slab
              z  = z0c;           // set to center for downstream angle calcs
            }
            Real rad = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (rad > jet_rinj) {
              continue; // outside geometric nozzle
            } else {
              ++dbg_in_rad; // count geometric nozzle cells
            }

            // Compute gate angles analogous to ProblemGenerator
            bool is2d = is2d_local;
            bool angle_ok = true;
            // Angles from Cartesian position relative to breakout center (coordinate-agnostic)
            Real zloc = z - z0c;
            Real ct_dir;
            if (is2d) {
              // In 2D slices, define theta_dir = π/2 so gating falls back to φ only (as intended)
              ct_dir = 0.0;
            } else {
              ct_dir = (rad > 0.0) ? (zloc / rad) : 1.0; // cos(theta)
            }
            ct_dir = std::max((Real)-1.0, std::min((Real)1.0, ct_dir));
            Real theta_dir = std::acos(ct_dir);
            Real xloc = x - x0c;
            Real yloc = y - y0c;
            Real phi_dir = std::atan2(yloc, xloc);
            if (phi_dir < 0.0) phi_dir += 2.0*M_PI;

            if (is2d && std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
              auto wrap_pm_pi = [](Real a)->Real { a = std::fmod(a + M_PI, 2.0*M_PI); if (a < 0.0) a += 2.0*M_PI; return a - M_PI; };
              Real dphi  = wrap_pm_pi(phi_dir - gate_phi0);
              bool angle_ok_primary = (dphi >= gate_dir_min && dphi <= gate_dir_max);
              bool angle_ok_bip = false;
              if (gate_cone_bipolar) {
                Real dphi2 = wrap_pm_pi(phi_dir - (gate_phi0 + M_PI));
                angle_ok_bip = (dphi2 >= gate_dir_min && dphi2 <= gate_dir_max);
              }
              angle_ok = angle_ok_primary || angle_ok_bip;
            }
            else if (is2d && std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
              auto wrap_pm_pi = [](Real a)->Real { a = std::fmod(a + M_PI, 2.0*M_PI); if (a < 0.0) a += 2.0*M_PI; return a - M_PI; };
              Real dphi  = wrap_pm_pi(phi_dir - gate_phi0);
              bool angle_ok_primary = (dphi >= gate_dir_min && dphi <= gate_dir_max);
              bool angle_ok_bip = false;
              if (gate_cone_bipolar) {
                Real dphi2 = wrap_pm_pi(phi_dir - (gate_phi0 + M_PI));
                angle_ok_bip = (dphi2 >= gate_dir_min && dphi2 <= gate_dir_max);
              }
              angle_ok = angle_ok_primary || angle_ok_bip;
            }
            else {
              angle_ok = (theta_dir >= gate_dir_min && theta_dir <= gate_dir_max);
            }

            if (!angle_ok) {
              continue;
            } else {
              ++dbg_angle_ok; // inside nozzle AND angle gate
            }

            // Enforce jet state: density, pressure, and velocity (purely radial)
            Real beta = 0.0;
            if (jet_Gam > 1.0) {
              Real invG2 = 1.0/(jet_Gam*jet_Gam);
              beta = std::sqrt(std::max(0.0, 1.0 - invG2));
            }
            // Spherical angles from Cartesian position (works for all coordinate systems)
            Real ct;
            if (is2d) {
              ct = 0.0; // theta = π/2 in 2D slice → purely in-plane
            } else {
              ct = (rad>0.0) ? ((z - z0c) / rad) : 1.0;  // cos(theta)
            }
            ct = std::max((Real)-1.0, std::min((Real)1.0, ct));
            Real th = std::acos(ct);
            Real ph = std::atan2(y - y0c, x - x0c);
            // Unit vector e_r in Cartesian
            Real sth = std::sin(th), cth = std::cos(th);
            Real cph = std::cos(ph),  sph = std::sin(ph);
            Real erx = sth*cph, ery = sth*sph, erz = cth;
            // Purely radial velocity in Cartesian
            Real vx = beta * erx;
            Real vy = beta * ery;
            Real vz = beta * erz;

            // Map Cartesian (vx,vy,vz) to native mesh-basis velocity components (v1,v2,v3)
            Real v1, v2c, v3;
            if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
              v1 = vx;
              v2c = vy;
              v3 = vz;
            } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
              // native basis is (v_R, v_phi, v_z); phi is the azimuth of this cell
              // Use the same φ we computed above for position
              Real vR   =  vx*std::cos(ph) + vy*std::sin(ph);
              Real vphi = -vx*std::sin(ph) + vy*std::cos(ph);
              v1 = vR; v2c = vphi; v3 = vz;
            } else { // spherical_polar: native basis is (v_r, v_theta, v_phi)
              // Orthonormal basis at (th, ph):
              // e_r     = (sinθ cosφ, sinθ sinφ, cosθ)  == (erx, ery, erz)
              // e_theta = (cosθ cosφ, cosθ sinφ, -sinθ)
              // e_phi   = (-sinφ, cosφ, 0)
              Real etx = cth*cph, ety = cth*sph, etz = -sth;
              Real ephx = -sph,   ephy =  cph,   ephz = 0.0;
              Real vr   = vx*erx + vy*ery + vz*erz;
              Real vth  = vx*etx + vy*ety + vz*etz;
              Real vph  = vx*ephx + vy*ephy + vz*ephz;
              v1 = vr; v2c = vth; v3 = vph;
            }

            // Write primitives in native basis
            Real rho = std::max(jet_rho, (Real)1e-30); //define a 1/r^2 rho profile? Define the density using luminosity and lorentz factor
            Real pgas = std::max(jet_p, (Real)1e-30);
            w(IDN,k,j,i) = rho;
            w(IPR,k,j,i) = pgas;
            w(IVX,k,j,i) = v1;   // v^1
            w(IVY,k,j,i) = v2c;  // v^2
            w(IVZ,k,j,i) = v3;   // v^3

#if defined(RELATIVISTIC_DYNAMICS) && (RELATIVISTIC_DYNAMICS != 0)
            // Consistent conserved variables (still scalar v^2 from physical speed)
            Real v2 = vx*vx + vy*vy + vz*vz; v2 = std::min(v2, 1.0 - 1e-12);
            Real gL = 1.0/std::sqrt(1.0 - v2);
            Real gamma = pmb->peos->GetGamma();
            Real h_spec  = 1.0 + (gamma/(gamma - 1.0)) * (pgas / rho);
            u(IDN,k,j,i) = rho * gL;
            u(IM1,k,j,i) = rho * h_spec * gL*gL * v1;
            u(IM2,k,j,i) = rho * h_spec * gL*gL * v2c;
            u(IM3,k,j,i) = rho * h_spec * gL*gL * v3;
            Real tau  = rho * h_spec * gL*gL - pgas - rho * gL;
            Real Etot = tau + rho * gL; // store Etot (include rest mass)
            u(IEN,k,j,i) = Etot;
#else
            Real v2mag = vx*vx + vy*vy + vz*vz;
            Real gm1 = pmb->peos->GetGamma() - 1.0;
            u(IDN,k,j,i) = rho;
            u(IM1,k,j,i) = rho * v1;
            u(IM2,k,j,i) = rho * v2c;
            u(IM3,k,j,i) = rho * v3;
            u(IEN,k,j,i) = pgas/gm1 + 0.5*rho*v2mag;
#endif
            // Count cells actually written with jet state this step
            ++dbg_written;
          }
        }
      }
    }
  }

  // Scan only cells within a thin spherical shell around r = rout, and angles in [0, 90 deg]
  for (int nb=0; nb<nblocal; ++nb) {
    MeshBlock* pmb = my_blocks(nb);
    AthenaArray<Real> &w = pmb->phydro->w;
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real x, y, z; cell_to_cart(pmb, k, j, i, x, y, z);
          Real dx = x - x0c, dy = y - y0c, dz = z - z0c;
          Real rad = std::sqrt(dx*dx + dy*dy + dz*dz);
          if (rad <= rout || rad > rout + ringw) continue;  // only near the stellar surface

          // in-plane azimuth measured from +x, wrap to [0, 2pi)
          Real phi = std::atan2(dy, dx);
          if (phi < 0.0) phi += 2.0*M_PI;
          const Real phi0 = breakout_phi0_deg * (M_PI/180.0); // center direction in radians
          Real phi_rel = phi - phi0;
          if (phi_rel < 0.0) phi_rel += 2.0*M_PI;
          if (phi_rel >= 2.0*M_PI) phi_rel -= 2.0*M_PI;
          if (phi_rel < 0.0 || phi_rel > (M_PI/2.0)) continue; // only angles in [0, 90 deg]
          int b = static_cast<int>( (phi_rel / (0.5*M_PI)) * nbins );
          if (b >= nbins) b = nbins-1;

          // Establish ambient reference at t<=0
          if (pamb_bin[b] < 0.0 && time <= 0.0) {
            pamb_bin[b] = w(IPR, k, j, i);
          }

          // If breakout not yet recorded, test condition on this cell (outside surface)
          if (t_break[b] < 0.0 && !logged[b]) {
            Real pamb = (pamb_bin[b] > 0.0 ? pamb_bin[b] : 1.0e-30);
            Real p  = w(IPR, k, j, i);
            Real vx = w(IVX, k, j, i), vy = w(IVY, k, j, i), vz = w(IVZ, k, j, i);
            Real vr = (rad > 0.0) ? (dx*vx + dy*vy + dz*vz)/rad : 0.0; // outward radial speed
            if ( (p > factor * pamb) && (vr > vmin) ) {
              t_break[b] = time;  // record first time this angle breaks out
              logged[b]  = 1;     // ensure we only log once per bin
              if (Globals::my_rank == 0) {
                FILE* f = std::fopen("shock_breakout.csv", "a");
                if (f) {
                  // bin center angle in degrees [0,90]
                  Real phi_center_deg = ((static_cast<Real>(b) + 0.5) * 90.0 / static_cast<Real>(nbins));
                  std::fprintf(f, "%g,%g,%g\n", phi_center_deg, rad, t_break[b]);
                  std::fclose(f);
                }
              }
            }
          }
        }
      }
    }
  }
}

//========================================================================================
//! \fn int RefinementCondition(MeshBlock *pmb)
//! \brief AMR refinement condition based on velocity shear
//========================================================================================


// refinement condition: check the maximum pressure gradient
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Coordinates &coord = *pmb->pcoord;

  // Optional radial cap (kept from existing code)
  Real rcen = coord.x1v(pmb->is) + coord.x1v(pmb->ie);
  rcen *= 0.5;
  if (rcen > 4.0) return 0;   // no AMR out beyond 4.0R

  // Dimensionality flags
  bool has_y = (pmb->block_size.nx2 > 1) || pmb->pmy_mesh->f2;
  bool has_z = (pmb->block_size.nx3 > 1) || pmb->pmy_mesh->f3;

  // Floors for safe normalization
  const Real pfloor  = (Real)1e-20;
  const Real rhofloor= (Real)1e-30;
  const Real tiny    = (Real)1e-99;

  // Coordinate-aware spacing helpers (cell-centered, handles edges)
  auto dx1_at = [&](int i)->Real {
    if (i == pmb->is) return coord.x1v(i+1) - coord.x1v(i);
    if (i == pmb->ie) return coord.x1v(i)   - coord.x1v(i-1);
    Real hp = coord.x1v(i+1) - coord.x1v(i);
    Real hm = coord.x1v(i)   - coord.x1v(i-1);
    return (Real)0.5*(hp + hm);
  };
  auto dx2_at = [&](int j)->Real {
    if (!has_y) return dx1_at(pmb->is); // unused, but return something finite
    if (j == pmb->js) return coord.x2v(j+1) - coord.x2v(j);
    if (j == pmb->je) return coord.x2v(j)   - coord.x2v(j-1);
    Real hp = coord.x2v(j+1) - coord.x2v(j);
    Real hm = coord.x2v(j)   - coord.x2v(j-1);
    return (Real)0.5*(hp + hm);
  };
  auto dx3_at = [&](int k)->Real {
    if (!has_z) return dx1_at(pmb->is);
    if (k == pmb->ks) return coord.x3v(k+1) - coord.x3v(k);
    if (k == pmb->ke) return coord.x3v(k)   - coord.x3v(k-1);
    Real hp = coord.x3v(k+1) - coord.x3v(k);
    Real hm = coord.x3v(k)   - coord.x3v(k-1);
    return (Real)0.5*(hp + hm);
  };

  // One-sided at edges, centered otherwise
  auto dA_dx1 = [&](int k,int j,int i,int comp)->Real {
    if (i == pmb->is)
      return (w(comp,k,j,i+1) - w(comp,k,j,i)) / (coord.x1v(i+1) - coord.x1v(i) + tiny);
    if (i == pmb->ie)
      return (w(comp,k,j,i)   - w(comp,k,j,i-1)) / (coord.x1v(i)   - coord.x1v(i-1) + tiny);
    return (w(comp,k,j,i+1) - w(comp,k,j,i-1)) / (coord.x1v(i+1) - coord.x1v(i-1) + tiny);
  };
  auto dA_dx2 = [&](int k,int j,int i,int comp)->Real {
    if (!has_y) return 0.0;
    if (j == pmb->js)
      return (w(comp,k,j+1,i) - w(comp,k,j,i)) / (coord.x2v(j+1) - coord.x2v(j) + tiny);
    if (j == pmb->je)
      return (w(comp,k,j,i)   - w(comp,k,j-1,i)) / (coord.x2v(j)   - coord.x2v(j-1) + tiny);
    return (w(comp,k,j+1,i) - w(comp,k,j-1,i)) / (coord.x2v(j+1) - coord.x2v(j-1) + tiny);
  };
  auto dA_dx3 = [&](int k,int j,int i,int comp)->Real {
    if (!has_z) return 0.0;
    if (k == pmb->ks)
      return (w(comp,k+1,j,i) - w(comp,k,j,i)) / (coord.x3v(k+1) - coord.x3v(k) + tiny);
    if (k == pmb->ke)
      return (w(comp,k,j,i)   - w(comp,k-1,j,i)) / (coord.x3v(k)   - coord.x3v(k-1) + tiny);
    return (w(comp,k+1,j,i) - w(comp,k-1,j,i)) / (coord.x3v(k+1) - coord.x3v(k-1) + tiny);
  };

  Real max_eta = 0.0; // composite indicator across the block

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // Local cell size for normalization (level-independent)
        Real dx1 = dx1_at(i);
        Real dx2 = has_y ? dx2_at(j) : dx1;
        Real dx3 = has_z ? dx3_at(k) : dx1;
        Real dmin = std::min(dx1, std::min(dx2, dx3));

        // Gradients
        Real gpx = dA_dx1(k,j,i, IPR);
        Real gpy = dA_dx2(k,j,i, IPR);
        Real gpz = dA_dx3(k,j,i, IPR);
        Real grx = dA_dx1(k,j,i, IDN);
        Real gry = dA_dx2(k,j,i, IDN);
        Real grz = dA_dx3(k,j,i, IDN);

        Real gp = std::sqrt(gpx*gpx + gpy*gpy + gpz*gpz);
        Real gr = std::sqrt(grx*grx + gry*gry + grz*grz);

        // Dimensionless indicators from pressure/density
        Real p   = std::max(w(IPR,k,j,i), pfloor);
        Real rho = std::max(w(IDN,k,j,i), rhofloor);
        Real eta_p = gp * dmin / p;
        Real eta_r = gr * dmin / rho;

        // Velocity magnitude (native basis components)
        Real vx = w(IVX,k,j,i);
        Real vy = has_y ? w(IVY,k,j,i) : 0.0;
        Real vz = has_z ? w(IVZ,k,j,i) : 0.0;
        Real v2 = vx*vx + vy*vy + vz*vz;
        Real vmag = std::sqrt(v2);

        // If below a configured speed, treat as static star: ignore rho/P gradients there
        if (vmag < amr_v_cut) {
          eta_p = 0.0;
          eta_r = 0.0;
        }

        // Optional velocity-based refinement indicator (e.g. for jets / blast wave)
        Real eta_v = 0.0;
        if (amr_v_ref > 0.0) {
          eta_v = vmag / amr_v_ref;
        }

        // Combine all indicators
        Real eta = eta_p;
        if (eta_r > eta) eta = eta_r;
        if (eta_v > eta) eta = eta_v;

        if (eta > max_eta) max_eta = eta;
      }
    }
  }

  // Hysteresis using existing 'threshold' input as tau_ref
  Real tau_ref   = threshold;
  Real tau_deref = (Real)0.9 * threshold; // prevent thrashing

  if (max_eta > tau_ref)   return 1;   // refine
  if (max_eta < tau_deref) return -1;  // derefine
  return 0;                            // keep
}