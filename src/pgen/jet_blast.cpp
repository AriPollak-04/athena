//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file blast.cpp
//! \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//!        cylindrical, and spherical coordinates.  Contains post-processing code
//!        to check whether blast is spherical for regression tests
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

Real threshold;

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
// --- physical jet normalization (optional): use total energy instead of fixed (rho,p)
static Real jet_E_tot = 0.0;     // total injected energy (excluding rest mass) over [0, t_stop]
static Real jet_mu     = 1.0;     // mean molecular weight (dimensionless, ~1 for baryons)
static Real jet_fth    = 0.01;    // fractional thermal content (0=cold; small numbers recommended)
static Real jet_L      = 0.0;     // computed luminosity = E_tot / t_stop (if enabled)
static Real gate_omega = -1.0;    // solid angle of the gate (3D only); -1 means compute on init
// driving mode: 0 = overwrite primitives in nozzle; 1 = conservative source (add ΔU)
static int  jet_mode = 0;
// optional convenience input: half-opening angle; if >=0 it sets dir_angle_min/max automatically
static Real theta_j = -1.0;
// ----------------------------------------------------------

int RefinementCondition(MeshBlock *pmb);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem","thr");
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
  // Physical (energy-normalized) jet options
  jet_E_tot = pin->GetOrAddReal("problem", "E_jet", 0.0);        // total energy over [0,t_stop]
  jet_mu    = pin->GetOrAddReal("problem", "jet_mu", 1.0);        // composition via mean molecular weight
  jet_fth   = pin->GetOrAddReal("problem", "jet_fth", 0.01);      // small thermal fraction (0..0.5 typical)
  jet_L     = (jet_t_stop > 0.0 ? jet_E_tot / jet_t_stop : 0.0);    // constant over the driving window

  // Mode & convenience angle input
  std::string jet_mode_s = pin->GetOrAddString("problem", "jet_mode", "overwrite");
  if (jet_mode_s == "source" || jet_mode_s == "conservative") jet_mode = 1; else jet_mode = 0;
  theta_j = pin->GetOrAddReal("problem", "theta_j", -1.0);
  if (theta_j >= 0.0) {
    bool is2d = (nblocal>0 ? (my_blocks(0)->ks == my_blocks(0)->ke) : true);
    if (is2d) { gate_dir_min = -theta_j; gate_dir_max = +theta_j; }
    else      { gate_dir_min = 0.0;      gate_dir_max =  theta_j; }
    // recompute Ω if 3D
    if (!is2d) {
      Real th_min = std::max((Real)0.0, gate_dir_min);
      Real th_max = std::min((Real)M_PI, gate_dir_max);
      th_min = std::min(th_min, th_max);
      gate_omega = 2.0 * M_PI * (std::cos(th_min) - std::cos(th_max));
    }
  }

  // For 3D runs, precompute the gate solid angle Ω; for 2D we will treat normalization per-unit-length
  if (nblocal > 0) {
    MeshBlock *pmb0 = my_blocks(0);
    bool is2d = (pmb0->ks == pmb0->ke);
    if (!is2d) {
      // θ gate: dir_angle_min..dir_angle_max, possibly bipolar in 2D only; for 3D we use θ only
      Real th_min = std::max((Real)0.0, gate_dir_min);
      Real th_max = std::min((Real)M_PI, gate_dir_max);
      th_min = std::min(th_min, th_max);
      // Solid angle of a band in θ: Ω = 2π (cos th_min - cos th_max)
      gate_omega = 2.0 * M_PI * (std::cos(th_min) - std::cos(th_max));
    } else {
      gate_omega = -1.0; // mark as 2D; we will use per-unit-length normalization
    }
  }

  if (Globals::my_rank == 0) {
    std::fprintf(stderr,
      "[jet:init] mode=%s enabled=%d  t_stop=%g  rinj=%g  Gam=%g  (rho=%g p=%g)  E_jet=%g  L=%g  mu=%g  fth=%g  Omega=%g  theta_j=%g  gate=[%g,%g] phi0=%g bipolar=%d\n",
      (jet_mode==1?"source":"overwrite"), (int)jet_enabled, (double)jet_t_stop, (double)jet_rinj, (double)jet_Gam,
      (double)jet_rho, (double)jet_p, (double)jet_E_tot, (double)jet_L,
      (double)jet_mu, (double)jet_fth, (double)gate_omega, (double)theta_j,
      (double)gate_dir_min, (double)gate_dir_max, (double)gate_phi0, (int)gate_cone_bipolar);
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
  Real dir_angle_min    = pin->GetReal("problem","dir_angle_min");
  Real dir_angle_max    = pin->GetReal("problem","dir_angle_max");

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
        // compute local azimuthal angle
        Real rr     = rad;
        Real r_star = rout;
        Real rho, pgas;
        // Interpolate from precomputed polytrope CSV
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



        // Default ambient velocities (will be overwritten inside BM shell)
        Real vx = 0.0, vy = 0.0, vz = 0.0;

        Real theta_dir;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          // polar angle θ = arccos((z - z0)/rad)
          Real zloc = pcoord->x3v(k) - z0;
          if (rad > 0.0) {
            Real ct = zloc / rad;
            if (ct > 1.0) ct = 1.0;
            if (ct < -1.0) ct = -1.0;
            theta_dir = std::acos(ct);
          } else {
            theta_dir = 0.0;
          }
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          // in cylindrical axisym, x2v(j) is the azimuthal angle φ, so θ = atan2(|z|, r)
          Real zloc = pcoord->x3v(k) - z0;
          Real rloc = pcoord->x1v(i);
          theta_dir = (rloc > 0.0 ? std::atan2(std::abs(zloc), rloc) : 0.0);
        } else {  // spherical_polar
          theta_dir = pcoord->x2v(j) - x2_0;
        }

        Real phi_dir;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          // azimuthal angle φ = atan2(y - y0, x - x0)
          // compute local x,y relative to blast center
          Real xloc = pcoord->x1v(i) - x0;
          Real yloc = pcoord->x2v(j) - y0;
          phi_dir = std::atan2(yloc, xloc);
          if (phi_dir < 0.0) phi_dir += 2.0*M_PI;
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          // φ coordinate in cylindrical, offset by x2_0
            phi_dir = pcoord->x2v(j) - x2_0;
            phi_dir = std::fmod(phi_dir + 2.0*M_PI, 2.0*M_PI);
        } else {  // spherical_polar
          // φ coordinate in spherical_polar, offset by x3_0
          phi_dir = pcoord->x3v(k) - x3_0;
          phi_dir = std::fmod(phi_dir + 2.0*M_PI, 2.0*M_PI);
        }

        // angle gating
        bool angle_ok = true;
        bool is2d = (ks == ke);

        // In 2D Cartesian, gate by in-plane azimuth φ around phi0 with width [dir_angle_min, dir_angle_max]
        if (is2d && std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Real phi0 = pin->GetOrAddReal("problem","phi0", 0.0);                    // center direction (radians), default +x
          bool cone_bipolar = pin->GetOrAddBoolean("problem","cone_bipolar", true); // also include opposite wedge by default

          auto wrap_pm_pi = [](Real a)->Real {
            a = std::fmod(a + M_PI, 2.0*M_PI);
            if (a < 0.0) a += 2.0*M_PI;
            return a - M_PI;
          };

          Real dphi = wrap_pm_pi(phi_dir - phi0);
          angle_ok = (dphi >= dir_angle_min && dphi <= dir_angle_max);

          if (cone_bipolar && !angle_ok) {
            Real dphi2 = wrap_pm_pi(phi_dir - (phi0 + M_PI));
            angle_ok = (dphi2 >= dir_angle_min && dphi2 <= dir_angle_max);
          }
        } else {
          // In 3D, keep θ gate relative to +z
          angle_ok = (theta_dir >= dir_angle_min && theta_dir <= dir_angle_max);
        }

        
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
        phydro->u(IEN,k,j,i) = Etot;
        // diagnostic: v^2 sanity
        if (!std::isfinite(v2)) nan_v2++;
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
//! \brief Track shock breakout times at r≈star_radius for angles in [0, 90°] and log CSV.
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
    for (int nb=0; nb<nblocal; ++nb) {
      MeshBlock* pmb = my_blocks(nb);
      auto &coord = *pmb->pcoord;
      AthenaArray<Real> &w = pmb->phydro->w;
      AthenaArray<Real> &u = pmb->phydro->u;

      // --- Physical normalization: if E_jet>0 use power L=E_jet/t_stop to set rho (cold-ish jet)
      // We assume a uniform energy flux across the spherical patch at r=rinj within the gate.
      bool use_energy_norm = (jet_L > 0.0 && jet_Gam > 1.0 && jet_rinj > 0.0);
      Real rho_from_power = jet_rho; // fallback to user rho if energy mode not active
      Real p_from_power   = jet_p;   // fallback
      if (use_energy_norm) {
        // Speed/lorentz
        Real invG2 = 1.0/(jet_Gam*jet_Gam);
        Real beta  = std::sqrt(std::max(0.0, 1.0 - invG2));
        // Effective area for 3D: A_eff = Ω * r_inj^2 (flux through spherical surface)
        // For 2D, we treat power per unit length along the out-of-plane axis: A_2D = (Δφ) * r_inj (units of length)
        Real A_eff = 0.0;
        bool is2d_block = (pmb->ks == pmb->ke);
        if (!is2d_block && gate_omega > 0.0) {
          A_eff = gate_omega * (jet_rinj * jet_rinj);
        } else {
          // 2D Cartesian wedge: approximate using azimuthal span around phi0
          Real dphi = gate_dir_max - gate_dir_min; if (dphi < 0.0) dphi = -dphi;
          if (gate_cone_bipolar) dphi *= 2.0;
          // Per-unit-length normalization (out-of-plane): use A2D = dphi * jet_rinj
          A_eff = std::max((Real)1e-30, dphi * jet_rinj);
        }
        // Cold-jet kinetic power density: L = (Gamma-1) * (rho * Gamma * beta) * A_eff
        // => rho = L / [Gamma * beta * (Gamma-1) * A_eff]
        rho_from_power = jet_L / std::max((Real)1e-30, (jet_Gam*beta*(jet_Gam - 1.0) * A_eff));
        // Add a small thermal fraction so h > 1; p = fth * rho (scaled by EOS)
        Real gamma = pmb->peos->GetGamma();
        Real eps = jet_fth;                    // fractional internal energy per unit rest-mass ~ p/(rho*(gamma-1))
        p_from_power = std::max((Real)1e-30, eps * (gamma - 1.0) * rho_from_power);
      }

      // Pre-pass: total injection volume for conservative source mode
      Real V_tot = 0.0;
      if (jet_mode == 1) {
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
          for (int j=pmb->js; j<=pmb->je; ++j) {
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              Real x,y,z; cell_to_cart(pmb, k, j, i, x, y, z);
              Real dx = x - breakout_x1_0, dy = y - breakout_x2_0, dz = z - breakout_x3_0;
              Real rad = std::sqrt(dx*dx + dy*dy + dz*dz);
              if (rad > jet_rinj) continue;
              // Minimal angle gate (same as main logic)
              bool is2d = (pmb->ks == pmb->ke);
              bool angle_ok = true;
              Real theta_dir, phi_dir;
              if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
                Real zloc = pmb->pcoord->x3v(k) - breakout_x3_0;
                if (rad > 0.0) { Real ct = zloc / rad; ct = std::max(-1.0, std::min(1.0, ct)); theta_dir = std::acos(ct); } else { theta_dir = 0.0; }
                Real xloc = pmb->pcoord->x1v(i) - breakout_x1_0;
                Real yloc = pmb->pcoord->x2v(j) - breakout_x2_0;
                phi_dir = std::atan2(yloc, xloc); if (phi_dir < 0.0) phi_dir += 2.0*M_PI;
              } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
                Real zloc = pmb->pcoord->x3v(k) - breakout_x3_0;
                Real rloc = pmb->pcoord->x1v(i);
                theta_dir = (rloc > 0.0 ? std::atan2(std::abs(zloc), rloc) : 0.0);
                phi_dir   = pmb->pcoord->x2v(j) - breakout_x2_0; phi_dir = std::fmod(phi_dir + 2.0*M_PI, 2.0*M_PI);
              } else { theta_dir = pmb->pcoord->x2v(j) - breakout_x2_0; phi_dir = pmb->pcoord->x3v(k) - breakout_x3_0; phi_dir = std::fmod(phi_dir + 2.0*M_PI, 2.0*M_PI); }
              if (is2d && std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
                auto wrap_pm_pi = [](Real a)->Real { a = std::fmod(a + M_PI, 2.0*M_PI); if (a < 0.0) a += 2.0*M_PI; return a - M_PI; };
                Real dphi = wrap_pm_pi(phi_dir - gate_phi0);
                angle_ok = (dphi >= gate_dir_min && dphi <= gate_dir_max);
                if (gate_cone_bipolar && !angle_ok) { Real dphi2 = wrap_pm_pi(phi_dir - (gate_phi0 + M_PI)); angle_ok = (dphi2 >= gate_dir_min && dphi2 <= gate_dir_max); }
              } else { angle_ok = (theta_dir >= gate_dir_min && theta_dir <= gate_dir_max); }
              if (!angle_ok) continue;
              V_tot += pmb->pcoord->GetCellVolume(k,j,i);
            }
          }
        }
        if (V_tot <= 0.0) V_tot = 1e-30;
      }

      long long inj_cells = 0;    // number of cells injected this step
      double    vmag_acc  = 0.0;  // accumulate |v|
      double    vmag_max  = 0.0;  // track max |v|

      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            // Cartesian position of cell center
            Real x, y, z; cell_to_cart(pmb, k, j, i, x, y, z);
            Real dx = x - breakout_x1_0, dy = y - breakout_x2_0, dz = z - breakout_x3_0;
            Real rad = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (rad > jet_rinj) continue; // nozzle region only

            // Compute gate angles analogous to ProblemGenerator
            bool is2d = (pmb->ks == pmb->ke);
            bool angle_ok = true;
            Real theta_dir, phi_dir;
            if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
              // theta from +z
              Real zloc = coord.x3v(k) - breakout_x3_0;
              if (rad > 0.0) {
                Real ct = zloc / rad; ct = std::max(-1.0, std::min(1.0, ct));
                theta_dir = std::acos(ct);
              } else { theta_dir = 0.0; }
              // in-plane azimuth
              Real xloc = coord.x1v(i) - breakout_x1_0;
              Real yloc = coord.x2v(j) - breakout_x2_0;
              phi_dir = std::atan2(yloc, xloc); if (phi_dir < 0.0) phi_dir += 2.0*M_PI;
            } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
              Real zloc = coord.x3v(k) - breakout_x3_0;
              Real rloc = coord.x1v(i);
              theta_dir = (rloc > 0.0 ? std::atan2(std::abs(zloc), rloc) : 0.0);
              phi_dir = coord.x2v(j) - breakout_x2_0; phi_dir = std::fmod(phi_dir + 2.0*M_PI, 2.0*M_PI);
            } else { // spherical_polar
              theta_dir = coord.x2v(j) - breakout_x2_0;
              phi_dir   = coord.x3v(k) - breakout_x3_0; phi_dir = std::fmod(phi_dir + 2.0*M_PI, 2.0*M_PI);
            }

            if (is2d && std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
              auto wrap_pm_pi = [](Real a)->Real { a = std::fmod(a + M_PI, 2.0*M_PI); if (a < 0.0) a += 2.0*M_PI; return a - M_PI; };
              Real dphi = wrap_pm_pi(phi_dir - gate_phi0);
              angle_ok = (dphi >= gate_dir_min && dphi <= gate_dir_max);
              if (gate_cone_bipolar && !angle_ok) {
                Real dphi2 = wrap_pm_pi(phi_dir - (gate_phi0 + M_PI));
                angle_ok = (dphi2 >= gate_dir_min && dphi2 <= gate_dir_max);
              }
            } else {
              angle_ok = (theta_dir >= gate_dir_min && theta_dir <= gate_dir_max);
            }

            if (!angle_ok) continue;

            // Enforce jet state: density, pressure, and radial velocity with Lorentz factor jet_Gam
            Real beta = 0.0;
            if (jet_Gam > 1.0) {
              Real invG2 = 1.0/(jet_Gam*jet_Gam);
              beta = std::sqrt(std::max(0.0, 1.0 - invG2));
            }
            // radial unit vector
            Real invr = (rad>0.0) ? 1.0/rad : 0.0;
            Real erx = invr * dx, ery = invr * dy, erz = invr * dz;
            Real vx = beta * erx;
            Real vy = beta * ery;
            Real vz = beta * erz;

            // debug: accumulate statistics
            Real vmag = std::sqrt(vx*vx + vy*vy + vz*vz);
            vmag_acc += vmag;
            if (vmag > vmag_max) vmag_max = vmag;
            ++inj_cells;

            if (jet_mode == 0) {
              // --- OVERWRITE MODE (existing behavior) ---
              Real rho, pgas;
              if (use_energy_norm) { rho = rho_from_power; pgas = p_from_power; }
              else { rho = std::max(jet_rho, (Real)1e-30); pgas = std::max(jet_p, (Real)1e-30); }
              w(IDN,k,j,i) = rho;
              w(IPR,k,j,i) = pgas;
              w(IVX,k,j,i) = vx;
              w(IVY,k,j,i) = vy;
              w(IVZ,k,j,i) = vz;

#if defined(RELATIVISTIC_DYNAMICS) && (RELATIVISTIC_DYNAMICS != 0)
              Real v2 = vx*vx + vy*vy + vz*vz; v2 = std::min(v2, 1.0 - 1e-12);
              Real gL = 1.0/std::sqrt(1.0 - v2);
              Real gamma = pmb->peos->GetGamma();
              Real h_spec  = 1.0 + (gamma/(gamma - 1.0)) * (pgas / rho);
              u(IDN,k,j,i) = rho * gL;
              u(IM1,k,j,i) = rho * h_spec * gL*gL * vx;
              u(IM2,k,j,i) = rho * h_spec * gL*gL * vy;
              u(IM3,k,j,i) = rho * h_spec * gL*gL * vz;
              // total energy density excluding rest mass: tau = rho*h*gamma^2 - p - rho*gamma
              u(IEN,k,j,i) = rho * h_spec * gL*gL - pgas - rho * gL;
#else
              Real v2 = vx*vx + vy*vy + vz*vz;
              Real gm1 = pmb->peos->GetGamma() - 1.0;
              u(IDN,k,j,i) = rho;
              u(IM1,k,j,i) = rho * vx;
              u(IM2,k,j,i) = rho * vy;
              u(IM3,k,j,i) = rho * vz;
              u(IEN,k,j,i) = pgas/gm1 + 0.5*rho*v2;
#endif
            } else {
              // --- CONSERVATIVE SOURCE MODE (add ΔU; do not overwrite primitives) ---
              Real invG2 = 1.0/(jet_Gam*jet_Gam);
              Real beta  = std::sqrt(std::max(0.0, 1.0 - invG2));
              Real dV    = pmb->pcoord->GetCellVolume(k,j,i);
              Real wgt   = dV / V_tot; // fraction of total source to this cell
              Real dtloc = pmb->pmy_mesh->dt;
              // Energy increment (tau excludes rest mass)
              Real dE = jet_L * dtloc * wgt;
              // Cold-jet relation: L = (Gamma-1) * Mdot  => Mdot = L/(Gamma-1)
              Real Mdot = jet_L / std::max((Real)1e-30, (jet_Gam - 1.0));
              Real dD   = jet_Gam * Mdot * dtloc * wgt;         // ΔD
              Real dS   = jet_Gam * beta * Mdot * dtloc * wgt;  // |ΔS|
              u(IDN,k,j,i) += dD;
              u(IM1,k,j,i) += dS * erx;
              u(IM2,k,j,i) += dS * ery;
              u(IM3,k,j,i) += dS * erz;
              u(IEN,k,j,i) += dE;  // tau increment
            }
          }
        }
      }
      if (Globals::my_rank == 0) {
        double vmag_avg = (inj_cells > 0 ? vmag_acc / (double)inj_cells : 0.0);
        std::fprintf(stderr,
          "[jet:step] t=%g nb=%d inj_cells=%lld vmag_avg=%g vmag_max=%g (rinj=%g, Gam=%g)\n",
          (double)time, nb, inj_cells, vmag_avg, vmag_max, (double)jet_rinj, (double)jet_Gam);
      }
      if (Globals::my_rank == 0 && jet_mode == 1) {
        std::fprintf(stderr, "[jet:src] t=%g nb=%d added_E=%g (L=%g dt=%g)\n",
          (double)time, nb, (double)(jet_L * pmb->pmy_mesh->dt), (double)jet_L, (double)pmb->pmy_mesh->dt);
      }
    }
  }
  else if (time <= jet_t_stop && Globals::my_rank == 0) {
    std::fprintf(stderr, "[jet:step] t=%g jet_enabled=0 (check t_stop, rinj, Gam, gating)\n", (double)time);
  }

  // Scan only cells within a thin spherical shell around r = rout, and angles in [0, 90 deg]
  for (int nb=0; nb<nblocal; ++nb) {
    MeshBlock* pmb = my_blocks(nb);
    AthenaArray<Real> &w = pmb->phydro->w;
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real x, y, z; cell_to_cart(pmb, k, j, i, x, y, z);
          Real dx = x - x1_0, dy = y - x2_0, dz = z - x3_0;
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
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Check radius of sphere to make sure it is round
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // delay AMR enrollment until after first timestep
  static bool enrolled_amr = false;
  if (adaptive && !enrolled_amr && time > 0.01) {
    EnrollUserRefinementCondition(RefinementCondition);
    enrolled_amr = true;
  }
  
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;
  MeshBlock *pmb = my_blocks(0);

  // analysis - check shape of the spherical blast wave
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> pr;
  pr.InitWithShallowSlice(pmb->phydro->w, 4, IPR, 1);

  // get coordinate location of the center, convert to Cartesian
  Real x1_0 = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0 = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0 = pin->GetOrAddReal("problem", "x3_0", 0.0);
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
    msg << "### FATAL ERROR in blast.cpp ParameterInput" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // find indices of the center
  int ic, jc, kc;
  for (ic=is; ic<=ie; ic++)
    if (pmb->pcoord->x1f(ic) > x1_0) break;
  ic--;
  for (jc=pmb->js; jc<=pmb->je; jc++)
    if (pmb->pcoord->x2f(jc) > x2_0) break;
  jc--;
  for (kc=pmb->ks; kc<=pmb->ke; kc++)
    if (pmb->pcoord->x3f(kc) > x3_0) break;
  kc--;

  // search pressure maximum in each direction
  Real rmax = 0.0, rmin = 100.0, rave = 0.0;
  int nr = 0;
  for (int o=0; o<=6; o++) {
    int ios = 0, jos = 0, kos = 0;
    if (o == 1) ios=-10;
    else if (o == 2) ios =  10;
    else if (o == 3) jos = -10;
    else if (o == 4) jos =  10;
    else if (o == 5) kos = -10;
    else if (o == 6) kos =  10;
    for (int d=0; d<6; d++) {
      Real pmax = 0.0;
      int imax(0), jmax(0), kmax(0);
      if (d == 0) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i>=is; i--) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 1) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i<=ie; i++) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 2) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j>=js; j--) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 3) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j<=je; j++) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 4) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k>=ks; k--) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      } else { // if (d == 5) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k<=ke; k++) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      }

      Real xm, ym, zm;
      Real x1m = pmb->pcoord->x1v(imax);
      Real x2m = pmb->pcoord->x2v(jmax);
      Real x3m = pmb->pcoord->x3v(kmax);
      if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
        xm = x1m;
        ym = x2m;
        zm = x3m;
      } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
        xm = x1m*std::cos(x2m);
        ym = x1m*std::sin(x2m);
        zm = x3m;
      } else {  // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
        xm = x1m*std::sin(x2m)*std::cos(x3m);
        ym = x1m*std::sin(x2m)*std::sin(x3m);
        zm = x1m*std::cos(x2m);
      }
      Real rad = std::sqrt(SQR(xm-x0)+SQR(ym-y0)+SQR(zm-z0));
      if (rad > rmax) rmax = rad;
      if (rad < rmin) rmin = rad;
      rave += rad;
      nr++;
    }
  }
  rave /= static_cast<Real>(nr);

  // use physical grid spacing at center of blast
  Real dr_max;
  Real  x1c = pmb->pcoord->x1v(ic);
  Real dx1c = pmb->pcoord->dx1f(ic);
  Real  x2c = pmb->pcoord->x2v(jc);
  Real dx2c = pmb->pcoord->dx2f(jc);
  Real dx3c = pmb->pcoord->dx3f(kc);
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    dr_max = std::max(std::max(dx1c, dx2c), dx3c);
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), dx3c);
  } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), x1c*std::sin(x2c)*dx3c);
  }
  Real deform=(rmax-rmin)/dr_max;

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign("blastwave-shape.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(),"r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(),"a",pfile)) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }

      // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(),"w")) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }
    }
    std::fprintf(pfile,"# Offset blast wave test in %s coordinates:\n",COORDINATE_SYSTEM);
    std::fprintf(pfile,"# Rmax       Rmin       Rave        Deformation\n");
    std::fprintf(pfile,"%e  %e  %e  %e \n",rmax,rmin,rave,deform);
    std::fclose(pfile);
  }
  return;
}


// refinement condition: check the maximum pressure gradient
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  // compute block center
  auto &coord = *pmb->pcoord;
  Real rcen = coord.x1v(pmb->is) + coord.x1v(pmb->ie);
  rcen *= 0.5;
  if (rcen > 4.0) return 0;   // no AMR out beyond 4.0R

  Real maxcurv = 0.0;
  if (pmb->pmy_mesh->f3) {
    for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
          Real p_im1 = w(IPR,k,j,i-1);
          Real p_i   = w(IPR,k,j,i);
          Real p_ip1 = w(IPR,k,j,i+1);
          Real d2p1 = (p_ip1 - 2.0*p_i + p_im1)
                      / SQR(pmb->pcoord->dx1f(i));
          Real p_jm1 = w(IPR,k,j-1,i);
          Real p_jp1 = w(IPR,k,j+1,i);
          Real d2p2 = (p_jp1 - 2.0*p_i + p_jm1)
                      / SQR(pmb->pcoord->dx2f(j));
          Real p_km1 = w(IPR,k-1,j,i);
          Real p_kp1 = w(IPR,k+1,j,i);
          Real d2p3 = (p_kp1 - 2.0*p_i + p_km1)
                      / SQR(pmb->pcoord->dx3f(k));
          Real curv = std::sqrt(d2p1*d2p1 + d2p2*d2p2 + d2p3*d2p3)
                      / p_i;
          maxcurv = std::max(maxcurv, curv);
        }
      }
    }
  }
  else if (pmb->pmy_mesh->f2) {
    int k = pmb->ks;
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real p_im1 = w(IPR,k,j,i-1);
        Real p_i   = w(IPR,k,j,i);
        Real p_ip1 = w(IPR,k,j,i+1);
        Real d2p1 = (p_ip1 - 2.0*p_i + p_im1)
                    / SQR(pmb->pcoord->dx1f(i));
        Real p_jm1 = w(IPR,k,j-1,i);
        Real p_jp1 = w(IPR,k,j+1,i);
        Real d2p2 = (p_jp1 - 2.0*p_i + p_jm1)
                    / SQR(pmb->pcoord->dx2f(j));
        Real curv = std::sqrt(d2p1*d2p1 + d2p2*d2p2)
                    / p_i;
        maxcurv = std::max(maxcurv, curv);
      }
    }
  }
  else {
    return 0;
  }

  if (maxcurv > threshold) return 1;
  if (maxcurv < 0.25*threshold) return -1;
  return 0;

}