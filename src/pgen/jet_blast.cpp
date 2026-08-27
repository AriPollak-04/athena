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
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../scalars/scalars.hpp"

#include "../parameter_input.hpp"


int RefinementCondition(MeshBlock *pmb);

Real threshold;
// Optional AMR velocity masking / trigger: configured in InitUserMeshData
// amr_v_cut: below this speed, ignore rho/P gradients (treat as static star)
// amr_v_ref: reference speed to scale the velocity-based indicator (eta_v = |v|/amr_v_ref)
static Real amr_v_cut = 0.0;
static Real amr_v_ref = 0.0;
// AMR indicator selector:
// 0: max |v| (beta) in block
// 1: max conserved energy density u(IEN)
// 2: max normalized gradients of ln(rho) and ln(p)
// 3: hybrid (use gradients in quasi-static star, velocity/energy in jet)
static int  amr_mode = 0;
// Optional energy scale used in hybrid mode (if <=0, energy is ignored)
static Real amr_e_ref = 0.0;

// --- continuous shock front tracking (configured in InitUserMeshData) ---
static Real shock_x1_0 = 0.0, shock_x2_0 = 0.0, shock_x3_0 = 0.0;
static int  shock_nbins_g       = 180;   // angular bins over [0, 2pi]
static Real shock_dt_log_g      = 0.1;   // log cadence (same units as time)
static Real shock_r_min_g       = 0.0;
static Real shock_r_max_g       = 2.0;
static int  shock_nr_g          = 512;   // radial bins
static Real shock_peak_height_g = 100.0; // min |d log10(S)/dr| for a peak
static int  shock_peak_dist_g   = 10;    // min separation between peaks in bins
static Real shock_max_dist_g    = 0.3;   // max distance from outermost peak to keep
static Real shock_t_last        = -1e99;
static bool shock_params_inited = false;
// --- jet driving parameters ---
static Real jet_t_stop = 0.0;      // stop time for driving
static Real jet_rinj   = 0.0;      // injection radius (nozzle size)
static Real jet_Gam    = 1.0;      // Lorentz factor of injected flow
static Real jet_rho    = 0.0;      // comoving rest-mass density in jet
static Real jet_p      = 0.0;      // gas pressure in jet
static Real gate_theta0 = M_PI;    // half-opening angle; inject within theta_0 of each pole
static Real gate_phi0   = 0.0;     // center direction for 2D Cartesian wedge (radians)
static Real jet_t_ramp  = 0.0;     // spin-down duration ending at t_stop (0 = off)
static Real jet_Gam_end = 1.0;     // Lorentz factor reached at the end of the spin-down
static bool jet_enabled = false;   // enable jet driving when inputs provided

//----------------------------------------------------------------------------------------
//! \fn static Real JetGammaOfTime(Real t)
//! \brief Injected Lorentz factor at time t.
//!
//! Full strength until t_stop - t_ramp, then a raised-cosine taper down to jet_Gam_end
//! at t_stop.  A hard cutoff (t_ramp = 0) releases the last-stamped Gamma = jet_Gam
//! material as a free-coasting slug that overtakes the decelerating jet head; tapering
//! Gamma instead makes every parcel slower than the one ahead of it, so the outflow
//! stretches rather than colliding.
//!
//! Energy bookkeeping: one second of spin-down is worth I = <Gam^2 v>/(Gam^2 v) of a
//! full-strength second (I = 0.38311 for 31 -> 1), so the effective drive duration is
//! t_stop - (1 - I) * t_ramp.  See jet_energy.py, which inverts this for jet_rho.

static Real JetGammaOfTime(Real t) {
  if (jet_t_ramp <= 0.0) return jet_Gam;
  Real t_on = jet_t_stop - jet_t_ramp;
  if (t <= t_on) return jet_Gam;
  Real s = std::min(std::max((t - t_on)/jet_t_ramp, (Real)0.0), (Real)1.0);
  return jet_Gam_end + (jet_Gam - jet_Gam_end)*0.5*(1.0 + std::cos((Real)M_PI*s));
}
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
// ∫ Γ dV (use with V_tot to form volume-avg Γ = Gamma_int / V_tot).
// Γ is taken as D/rho from the conserved and primitive densities, which is exact and
// needs no metric factors.  This previously integrated the adiabatic index by mistake,
// which made the column a constant 4/3 * V_tot rather than a Lorentz-factor diagnostic.
static Real Hst_GammaInt(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real gL = 1.0;
#if defined(RELATIVISTIC_DYNAMICS) && (RELATIVISTIC_DYNAMICS != 0)
        const Real rho = pmb->phydro->w(IDN,k,j,i);
        if (rho > 0.0) gL = pmb->phydro->u(IDN,k,j,i) / rho;   // D / rho
#endif
        sum += gL * pco->GetCellVolume(k,j,i);
      }
    }
  }
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
    amr_mode  = pin->GetOrAddInteger("problem", "amr_mode", 0);
    amr_v_cut = pin->GetOrAddReal("problem", "amr_v_cut", 0.0);
    amr_v_ref = pin->GetOrAddReal("problem", "amr_v_ref", 0.0);
    amr_e_ref = pin->GetOrAddReal("problem", "amr_e_ref", 0.0);
    if (Globals::my_rank == 0) {
      std::fprintf(stderr,
        "[AMR:init] thr=%g  amr_mode=%d  amr_v_cut=%g  amr_v_ref=%g  amr_e_ref=%g\n",
        (double)threshold, (int)amr_mode, (double)amr_v_cut, (double)amr_v_ref, (double)amr_e_ref);
    }
  }


  // Shock front tracker inputs
  shock_x1_0          = pin->GetOrAddReal("problem", "x1_0", 0.0);
  shock_x2_0          = pin->GetOrAddReal("problem", "x2_0", 0.0);
  shock_x3_0          = pin->GetOrAddReal("problem", "x3_0", 0.0);
  shock_nbins_g       = pin->GetOrAddInteger("problem", "shock_nbins", 180);
  shock_dt_log_g      = pin->GetOrAddReal("problem", "shock_dt_log", 0.1);
  shock_r_min_g       = pin->GetOrAddReal("problem", "shock_r_min", 0.0);
  shock_r_max_g       = pin->GetOrAddReal("problem", "shock_r_max", 2.0);
  shock_nr_g          = pin->GetOrAddInteger("problem", "shock_nr", 512);
  shock_peak_height_g = pin->GetOrAddReal("problem", "shock_peak_height", 100.0);
  shock_peak_dist_g   = pin->GetOrAddInteger("problem", "shock_peak_dist", 10);
  shock_max_dist_g    = pin->GetOrAddReal("problem", "shock_max_dist", 0.3);
  if (Globals::my_rank == 0) {
    FILE* fcsv = std::fopen("shock_front.csv", "w");
    if (fcsv) {
      std::fprintf(fcsv, "time,angle_deg,r_shock,vr,v_tang,mushroom\n");
      std::fclose(fcsv);
    }
  }
  // Jet driving inputs (optional)
  jet_t_stop = pin->GetOrAddReal("problem", "t_stop", 0.0);
  jet_rinj   = pin->GetOrAddReal("problem", "jet_rinj", 0.0);
  jet_Gam    = pin->GetOrAddReal("problem", "jet_Gam", 1.0);
  jet_rho    = pin->GetOrAddReal("problem", "jet_rho", 0.0);
  jet_p      = pin->GetOrAddReal("problem", "jet_p", 0.0);
  gate_theta0 = pin->GetOrAddReal("problem", "theta_0", M_PI);
  gate_phi0   = pin->GetOrAddReal("problem", "phi0", 0.0);
  jet_t_ramp  = pin->GetOrAddReal("problem", "t_ramp", 0.0);
  jet_Gam_end = pin->GetOrAddReal("problem", "jet_Gam_end", 1.0);
  // Spin-down must fit inside the drive window and must not speed the jet back up
  jet_t_ramp  = std::min(std::max(jet_t_ramp, (Real)0.0), jet_t_stop);
  jet_Gam_end = std::min(std::max(jet_Gam_end, (Real)1.0), jet_Gam);
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
      "[jet:init] enabled=%d t_stop=%g rinj=%g Gam=%g rho=%g p=%g theta_0=%g phi0=%g "
      "t_ramp=%g Gam_end=%g (always bipolar)\n",
      (int)jet_enabled, (double)jet_t_stop, (double)jet_rinj, (double)jet_Gam,
      (double)jet_rho, (double)jet_p,
      (double)gate_theta0, (double)gate_phi0,
      (double)jet_t_ramp, (double)jet_Gam_end);
  }
  shock_params_inited = true;
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


  bool is2d = (ks == ke); // true for 2D runs (single zone in x3)

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
          z = is2d ? z0 : pcoord->x3v(k);  // 2D slab: collapse to equatorial plane
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
            rho0  = 1.0e-25;
            pgas0 = 1.0e-25;
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
        // Create angle_ok gate for passive scalar injection (mirrors UserWorkInLoop)
        bool angle_ok = true;
        if (jet_enabled) {
          Real xloc = x - x0;
          Real yloc = y - y0;
          Real zloc = z - z0;
          Real phi_dir = std::atan2(yloc, xloc);
          if (phi_dir < 0.0) phi_dir += 2.0*M_PI;

          if (is2d && std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            auto wrap_pm_pi = [](Real a)->Real { a = std::fmod(a + M_PI, 2.0*M_PI); if (a < 0.0) a += 2.0*M_PI; return a - M_PI; };
            Real dphi = wrap_pm_pi(phi_dir - gate_phi0);
            angle_ok = (std::fabs(dphi) <= gate_theta0) || (std::fabs(dphi) >= M_PI - gate_theta0);
          } else if (is2d && std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            auto wrap_pm_pi = [](Real a)->Real { a = std::fmod(a + M_PI, 2.0*M_PI); if (a < 0.0) a += 2.0*M_PI; return a - M_PI; };
            Real dphi = wrap_pm_pi(phi_dir - gate_phi0);
            angle_ok = (std::fabs(dphi) <= gate_theta0) || (std::fabs(dphi) >= M_PI - gate_theta0);
          } else {
            // 3D: bipolar cone check using position relative to center
            Real ct_dir = (rad > 0.0) ? (zloc / rad) : 1.0;
            ct_dir = std::max((Real)-1.0, std::min((Real)1.0, ct_dir));
            Real theta_dir = std::acos(ct_dir);
            angle_ok = (theta_dir <= gate_theta0) || (theta_dir >= M_PI - gate_theta0);
          }
        }



        // Initialize passive scalar: r=1 inside jet nozzle (rad <= jet_rinj AND angle gate), r=0 elsewhere
        if (NSCALARS > 0) {
          for (int n = 0; n < NSCALARS; ++n) {
            pscalars->s(n,k,j,i) = 0.0;
            pscalars->r(n,k,j,i) = 0.0;
            if (angle_ok && rad <= jet_rinj) {
              pscalars->s(n,k,j,i) = rho; // s = rho * r, r=1
              pscalars->r(n,k,j,i) = 1.0;
            }
          }
        }

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
  if (!shock_params_inited) return;

  // Origin in Cartesian
  Real x0c, y0c, z0c;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0c = shock_x1_0; y0c = shock_x2_0; z0c = shock_x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0c = shock_x1_0*std::cos(shock_x2_0);
    y0c = shock_x1_0*std::sin(shock_x2_0);
    z0c = shock_x3_0;
  } else {
    x0c = shock_x1_0*std::sin(shock_x2_0)*std::cos(shock_x3_0);
    y0c = shock_x1_0*std::sin(shock_x2_0)*std::sin(shock_x3_0);
    z0c = shock_x1_0*std::cos(shock_x2_0);
  }

  auto cell_to_cart = [&](MeshBlock* pmb, int k, int j, int i, Real &x, Real &y, Real &z) {
    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
      x = pmb->pcoord->x1v(i); y = pmb->pcoord->x2v(j); z = pmb->pcoord->x3v(k);
    } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
      Real R = pmb->pcoord->x1v(i), ph = pmb->pcoord->x2v(j);
      x = R*std::cos(ph); y = R*std::sin(ph); z = pmb->pcoord->x3v(k);
    } else {
      Real r = pmb->pcoord->x1v(i), th = pmb->pcoord->x2v(j), ph = pmb->pcoord->x3v(k);
      x = r*std::sin(th)*std::cos(ph); y = r*std::sin(th)*std::sin(ph); z = r*std::cos(th);
    }
  };

  auto native_vel_to_cart = [&](MeshBlock* pmb, int k, int j, int i,
                                Real, Real, Real, Real &vx, Real &vy, Real &vz) {
    Real v1 = pmb->phydro->w(IVX,k,j,i);
    Real v2 = pmb->phydro->w(IVY,k,j,i);
    Real v3 = pmb->phydro->w(IVZ,k,j,i);
    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
      vx = v1; vy = v2; vz = v3;
    } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
      Real ph = pmb->pcoord->x2v(j);
      vx = v1*std::cos(ph) - v2*std::sin(ph);
      vy = v1*std::sin(ph) + v2*std::cos(ph);
      vz = v3;
    } else {
      Real th = pmb->pcoord->x2v(j), ph = pmb->pcoord->x3v(k);
      Real sth = std::sin(th), cth = std::cos(th), sph = std::sin(ph), cph = std::cos(ph);
      vx = v1*sth*cph + v2*cth*cph - v3*sph;
      vy = v1*sth*sph + v2*cth*sph + v3*cph;
      vz = v1*cth - v2*sth;
    }
  };

  // ---- Jet driving: enforce jet state inside nozzle while time <= t_stop ----
  if (jet_enabled && (time <= jet_t_stop)) {
    const Real Gam_t = JetGammaOfTime(time);   // full strength, or tapering to Gam_end
    for (int nb=0; nb<nblocal; ++nb) {
      MeshBlock* pmb = my_blocks(nb);
      AthenaArray<Real> &w = pmb->phydro->w;
      AthenaArray<Real> &u = pmb->phydro->u;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real x, y, z; cell_to_cart(pmb, k, j, i, x, y, z);
            Real dx = x - x0c, dy = y - y0c, dz = z - z0c;
            bool is2d = (pmb->ks == pmb->ke);
            if (is2d) { dz = 0.0; z = z0c; }
            Real rad = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (rad > jet_rinj) continue;

            Real zloc = z - z0c, xloc = x - x0c, yloc = y - y0c;
            Real ct_dir = is2d ? 0.0 : ((rad > 0.0) ? (zloc/rad) : 1.0);
            ct_dir = std::max((Real)-1.0, std::min((Real)1.0, ct_dir));
            Real theta_dir = std::acos(ct_dir);
            Real phi_dir = std::atan2(yloc, xloc);
            if (phi_dir < 0.0) phi_dir += 2.0*M_PI;

            bool angle_ok = true;
            if (is2d && (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0 ||
                         std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0)) {
              auto wrap_pm_pi = [](Real a)->Real {
                a = std::fmod(a + M_PI, 2.0*M_PI); if (a < 0.0) a += 2.0*M_PI; return a - M_PI;
              };
              Real dphi = wrap_pm_pi(phi_dir - gate_phi0);
              angle_ok = (std::fabs(dphi) <= gate_theta0) || (std::fabs(dphi) >= M_PI - gate_theta0);
            } else {
              angle_ok = (theta_dir <= gate_theta0) || (theta_dir >= M_PI - gate_theta0);
            }
            if (!angle_ok) continue;

            Real beta = 0.0;
            if (Gam_t > 1.0) beta = std::sqrt(std::max(0.0, 1.0 - 1.0/(Gam_t*Gam_t)));
            Real ct = is2d ? 0.0 : ((rad > 0.0) ? ((z - z0c)/rad) : 1.0);
            ct = std::max((Real)-1.0, std::min((Real)1.0, ct));
            Real th = std::acos(ct);
            Real ph = std::atan2(y - y0c, x - x0c);
            Real sth = std::sin(th), cth = std::cos(th), cph = std::cos(ph), sph = std::sin(ph);
            Real erx = sth*cph, ery = sth*sph, erz = cth;
            Real vx = beta*erx, vy = beta*ery, vz = beta*erz;

            Real v1, v2c, v3;
            if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
              v1 = vx; v2c = vy; v3 = vz;
            } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
              v1  =  vx*std::cos(ph) + vy*std::sin(ph);
              v2c = -vx*std::sin(ph) + vy*std::cos(ph);
              v3  = vz;
            } else {
              Real etx = cth*cph, ety = cth*sph, etz = -sth;
              Real ephx = -sph, ephy = cph, ephz = 0.0;
              v1  = vx*erx  + vy*ery  + vz*erz;
              v2c = vx*etx  + vy*ety  + vz*etz;
              v3  = vx*ephx + vy*ephy + vz*ephz;
            }

            Real rho  = std::max(jet_rho, (Real)1e-30);
            Real pgas = std::max(jet_p,   (Real)1e-30);
            w(IDN,k,j,i) = rho; w(IPR,k,j,i) = pgas;
            w(IVX,k,j,i) = v1; w(IVY,k,j,i) = v2c; w(IVZ,k,j,i) = v3;

#if defined(RELATIVISTIC_DYNAMICS) && (RELATIVISTIC_DYNAMICS != 0)
            Real v2 = vx*vx + vy*vy + vz*vz; v2 = std::min(v2, 1.0 - 1e-12);
            Real gL = 1.0/std::sqrt(1.0 - v2);
            Real gamma_eos = pmb->peos->GetGamma();
            Real h_spec = 1.0 + (gamma_eos/(gamma_eos - 1.0)) * (pgas/rho);
            u(IDN,k,j,i) = rho*gL;
            u(IM1,k,j,i) = rho*h_spec*gL*gL*v1;
            u(IM2,k,j,i) = rho*h_spec*gL*gL*v2c;
            u(IM3,k,j,i) = rho*h_spec*gL*gL*v3;
            Real tau  = rho*h_spec*gL*gL - pgas - rho*gL;
            u(IEN,k,j,i) = tau + rho*gL;
#else
            Real v2mag = vx*vx + vy*vy + vz*vz;
            Real gm1 = pmb->peos->GetGamma() - 1.0;
            u(IDN,k,j,i) = rho;
            u(IM1,k,j,i) = rho*v1; u(IM2,k,j,i) = rho*v2c; u(IM3,k,j,i) = rho*v3;
            u(IEN,k,j,i) = pgas/gm1 + 0.5*rho*v2mag;
#endif
          }
        }
      }
    }
  }

  // ---- Continuous shock front tracking at shock_dt_log cadence ----
  if (time - shock_t_last < shock_dt_log_g) return;

  const int  nb_ang = shock_nbins_g;
  const int  nr     = shock_nr_g;
  const Real rmin   = shock_r_min_g;
  const Real rmax   = shock_r_max_g;
  const Real dr     = (rmax - rmin) / static_cast<Real>(nr);
  const int  ps     = nb_ang * nr;   // total profile size

  // Accumulation arrays: [angle_bin * nr + r_bin]
  std::vector<Real> logS_sum(ps, 0.0), vr_sum(ps, 0.0), vt_sum(ps, 0.0), cnt(ps, 0.0);

  for (int inb = 0; inb < nblocal; ++inb) {
    MeshBlock* pmb = my_blocks(inb);
    AthenaArray<Real> &w = pmb->phydro->w;
    const Real gamma = pmb->peos->GetGamma();
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          Real x, y, z; cell_to_cart(pmb, k, j, i, x, y, z);
          Real dx_c = x - x0c, dy_c = y - y0c, dz_c = z - z0c;
          if (pmb->ks == pmb->ke) dz_c = 0.0;
          Real rad = std::sqrt(dx_c*dx_c + dy_c*dy_c + dz_c*dz_c);
          if (rad < rmin || rad >= rmax) continue;
          int ri = static_cast<int>((rad - rmin) / dr);
          if (ri < 0 || ri >= nr) continue;

          // Full [0, 2pi) azimuthal binning
          Real phi = std::atan2(dy_c, dx_c);
          if (phi < 0.0) phi += 2.0*M_PI;
          int bi = static_cast<int>((phi / (2.0*M_PI)) * nb_ang);
          if (bi < 0) bi = 0;
          if (bi >= nb_ang) bi = nb_ang - 1;

          Real rho = w(IDN,k,j,i), p = w(IPR,k,j,i);
          if (rho <= 0.0 || p <= 0.0 || !std::isfinite(rho) || !std::isfinite(p)) continue;
          Real S = p / std::pow(rho, gamma);
          if (S <= 0.0 || !std::isfinite(S)) continue;

          Real cvx, cvy, cvz;
          native_vel_to_cart(pmb, k, j, i, x, y, z, cvx, cvy, cvz);
          Real vr_c = (rad > 0.0) ? (dx_c*cvx + dy_c*cvy + dz_c*cvz) / rad : 0.0;
          Real vt_c = (rad > 0.0) ? (-dy_c*cvx + dx_c*cvy) / rad : 0.0;

          int idx = bi * nr + ri;
          logS_sum[idx] += std::log10(S);
          vr_sum  [idx] += vr_c;
          vt_sum  [idx] += vt_c;
          cnt     [idx] += 1.0;
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, logS_sum.data(), ps, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, vr_sum.data(),   ps, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, vt_sum.data(),   ps, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, cnt.data(),       ps, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (Globals::my_rank == 0) {
    FILE* f = std::fopen("shock_front.csv", "a");
    if (f) {
      for (int bi = 0; bi < nb_ang; ++bi) {
        // |d(log10 S)/dr| via central differences; skip bins with no data on either side
        std::vector<Real> abs_grad(nr, 0.0);
        for (int ri = 1; ri < nr - 1; ++ri) {
          int im = bi*nr + ri - 1, ip = bi*nr + ri + 1;
          if (cnt[im] <= 0.0 || cnt[ip] <= 0.0) continue;
          abs_grad[ri] = std::fabs(
            (logS_sum[ip]/cnt[ip] - logS_sum[im]/cnt[im]) / (2.0*dr));
        }

        // Local maxima above height threshold with minimum bin separation (mirrors find_peaks)
        const Real ht  = shock_peak_height_g;
        const int  dst = shock_peak_dist_g;
        std::vector<int> peaks;
        for (int ri = 1; ri < nr - 1; ++ri) {
          if (abs_grad[ri] < ht) continue;
          if (abs_grad[ri] <= abs_grad[ri-1] || abs_grad[ri] <= abs_grad[ri+1]) continue;
          if (!peaks.empty() && ri - peaks.back() < dst) {
            if (abs_grad[ri] > abs_grad[peaks.back()]) peaks.back() = ri;
          } else {
            peaks.push_back(ri);
          }
        }

        // r at center of bin ri
        auto r_of = [&](int ri_) { return rmin + (ri_ + 0.5)*dr; };

        // Keep peaks with r > 0.2, sort descending by r
        std::vector<int> valid;
        for (int pi : peaks)
          if (r_of(pi) > 0.2) valid.push_back(pi);
        std::sort(valid.begin(), valid.end(),
                  [&](int a, int b_) { return r_of(a) > r_of(b_); });

        // Outermost peak + any within max_dist of it (mushroom-shaped secondary fronts)
        std::vector<int> selected;
        if (!valid.empty()) {
          Real r_outer = r_of(valid[0]);
          for (int pi : valid)
            if (r_outer - r_of(pi) <= shock_max_dist_g) selected.push_back(pi);
        }

        // Also keep stellar surface peak near R=1.0 if not already captured
        std::vector<int> surf;
        for (int pi : valid)
          if (r_of(pi) > 0.85 && r_of(pi) < 1.15) surf.push_back(pi);
        if (!surf.empty()) {
          int best = *std::min_element(surf.begin(), surf.end(),
            [&](int a, int b_) { return std::fabs(r_of(a)-1.0) < std::fabs(r_of(b_)-1.0); });
          bool already = false;
          for (int pi : selected)
            if (std::fabs(r_of(pi) - r_of(best)) < 0.05) { already = true; break; }
          if (!already) selected.push_back(best);
        }

        Real angle_deg = ((Real)bi + 0.5) * 360.0 / (Real)nb_ang;
        bool has_surface = false, has_outer = false;
        for (int pi : selected) {
          Real rv = r_of(pi);
          if (rv > 0.95 && rv < 1.05) has_surface = true;
          if (rv > 1.05)              has_outer   = true;
        }
        int mushroom = (has_surface && has_outer) ? 1 : 0;
        for (int pi : selected) {
          int idx = bi*nr + pi;
          std::fprintf(f, "%.6g,%.4g,%.6g,%.6g,%.6g,%d\n",
            (double)time, (double)angle_deg,
            (double)r_of(pi),
            (double)(vr_sum[idx]/cnt[idx]),
            (double)(vt_sum[idx]/cnt[idx]),
            mushroom);
        }
      }
      std::fclose(f);
    }
  }

  shock_t_last = time;
}

//========================================================================================
//! \fn int RefinementCondition(MeshBlock *pmb)
//! \brief AMR refinement condition based on velocity shear
//========================================================================================


// refinement condition: check the maximum pressure gradient
int RefinementCondition(MeshBlock *pmb) {
  // ---------------------------------------------------------------------------
  // AMR INDICATORS
  // We compute several candidate indicators and then select/combine them based
  // on `amr_mode`.
  // ---------------------------------------------------------------------------
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
    if (!has_y) return dx1_at(pmb->is);
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

  // Max physical speed (beta) in this block.
  // NOTE: In this problem generator + jet driver, primitives store 3-velocity.
  Real max_beta = 0.0;

  // Max conserved total energy density in this block (u(IEN)).
  Real max_e = 0.0;

  // Max normalized gradients of ln(rho) and ln(p)
  Real max_grad = 0.0;

  AthenaArray<Real> &u = pmb->phydro->u;

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // --- speed indicator (beta) ---
        Real vx = w(IVX,k,j,i);
        Real vy = has_y ? w(IVY,k,j,i) : 0.0;
        Real vz = has_z ? w(IVZ,k,j,i) : 0.0;
        Real v2  = vx*vx + vy*vy + vz*vz;
#if defined(RELATIVISTIC_DYNAMICS) && (RELATIVISTIC_DYNAMICS != 0)
        // Keep subluminal
        v2 = std::min(v2, (Real)(1.0 - 1e-12));
#endif
        Real beta = std::sqrt(std::max((Real)0.0, v2));
        if (beta > max_beta) max_beta = beta;

        // --- energy indicator ---
        Real ecell = u(IEN,k,j,i);
        if (std::isfinite(ecell) && ecell > max_e) max_e = ecell;

        // --- gradient indicator (dimensionless) ---
        // Use ln(rho) and ln(p) to weight relative changes. Add floors to avoid NaNs.
        Real rho = std::max(w(IDN,k,j,i), rhofloor);
        Real p   = std::max(w(IPR,k,j,i), pfloor);

        // Compute derivatives of ln(rho) and ln(p)
        // d ln A / dx = (1/A) dA/dx
        Real drho_dx1 = dA_dx1(k,j,i,IDN) / rho;
        Real dp_dx1   = dA_dx1(k,j,i,IPR) / p;
        Real gr1 = std::fabs(drho_dx1)*dx1_at(i) + std::fabs(dp_dx1)*dx1_at(i);

        Real gr2 = 0.0;
        if (has_y) {
          Real drho_dx2 = dA_dx2(k,j,i,IDN) / rho;
          Real dp_dx2   = dA_dx2(k,j,i,IPR) / p;
          gr2 = std::fabs(drho_dx2)*dx2_at(j) + std::fabs(dp_dx2)*dx2_at(j);
        }

        Real gr3 = 0.0;
        if (has_z) {
          Real drho_dx3 = dA_dx3(k,j,i,IDN) / rho;
          Real dp_dx3   = dA_dx3(k,j,i,IPR) / p;
          gr3 = std::fabs(drho_dx3)*dx3_at(k) + std::fabs(dp_dx3)*dx3_at(k);
        }

        // Dimensionless gradient magnitude proxy
        Real gtot = std::sqrt(gr1*gr1 + gr2*gr2 + gr3*gr3);
        if (std::isfinite(gtot) && gtot > max_grad) max_grad = gtot;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // SELECT / COMBINE INDICATOR
  // ---------------------------------------------------------------------------
  Real indicator = 0.0;

  if (amr_mode == 0) {
    // pure speed
    indicator = max_beta;
  } else if (amr_mode == 1) {
    // pure energy density
    indicator = max_e;
  } else if (amr_mode == 2) {
    // pure gradients
    indicator = max_grad;
  } else {
    // hybrid: in quasi-static regions, prioritize gradients; in fast jet regions,
    // prioritize speed, plus optionally energy (scaled by amr_e_ref).
    Real eta_v = 0.0;
    if (amr_v_ref > 0.0) eta_v = max_beta / amr_v_ref;

    Real eta_e = 0.0;
    if (amr_e_ref > 0.0) eta_e = max_e / amr_e_ref;

    if (max_beta < amr_v_cut) {
      indicator = max_grad;
    } else {
      indicator = std::max(max_grad, std::max(eta_v, eta_e));
    }
  }

  // Hysteresis using `threshold` against the selected indicator.
  const Real thr_ref   = threshold;
  const Real thr_deref = (Real)0.9 * threshold; // prevent thrashing

  if (Globals::my_rank == 0 && pmb->gid == 0) {
    std::fprintf(stderr,
      "[AMR:step] t=%g  mode=%d  ind=%g  beta=%g  grad=%g  E=%g  thr=%g\n",
      (double)pmb->pmy_mesh->time, (int)amr_mode,
      (double)indicator, (double)max_beta, (double)max_grad, (double)max_e, (double)threshold);
  }

  if (indicator > thr_ref)   return 1;   // refine
  if (indicator < thr_deref) return -1;  // derefine
  return 0;                               // keep
}