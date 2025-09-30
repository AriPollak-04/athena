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
  const bool use_BM_R0 = pin->GetBoolean("problem", "use_BM_R0");
  Real r_blast      = pin->GetReal("problem","r_blast");
  Real cs_factor = pin->GetReal("problem","cs_factor");
  Real dir_angle_min    = pin->GetReal("problem","dir_angle_min");
  Real dir_angle_max    = pin->GetReal("problem","dir_angle_max");

  Real E_iso = pin->GetReal("problem","E_iso");
  Real mp = pin->GetReal("hydro","mp");
  Real Gam0 = pin->GetReal("problem","Gam0");

  // Compile-time detector for dynamics mode
#if defined(RELATIVISTIC_DYNAMICS) && (RELATIVISTIC_DYNAMICS != 0)
  const char* dyn_mode = "SR";
#else
  const char* dyn_mode = "NR";
#endif

  // --- diagnostics ---
  long long tot_cells = 0, bm_cells = 0;
  long long nan_rho = 0, nan_p = 0, nan_v2 = 0;
  Real rho_min=1e99, rho_max=-1e99, p_min=1e99, p_max=-1e99;
  Real chi_min_seen = 1e99, chi_max_seen = -1e99;
  // --- end diagnostics ---

    // Load precomputed polytrope data from CSV
  static PolytropeData poly;
  if (poly.r.empty()) {
    poly = LoadPolytropeCSV(poly_csv_path);
  }

  Real R0 = 0.0;
  const Real chi_max = 1.0 + 8.0*Gam0*Gam0; // ensure full shell coverage
  if (use_BM_R0) {
    // choose upstream density for BM normalization (use polytrope at r = rout or an input-specified r_up)
    Real r_up = pin->GetOrAddReal("problem","r_up_for_BM", rout/3);
    Real rho_up_for_R0 = 0.0, p_up_for_R0 = 0.0;
    {
      const auto &r_arr   = poly.r;
      const auto &rho_arr = poly.rho;
      const auto &P_arr   = poly.P;
      int N = static_cast<int>(r_arr.size());
      auto it = std::lower_bound(r_arr.begin(), r_arr.end(), r_up);
      int idx = static_cast<int>(std::distance(r_arr.begin(), it));
      if (idx <= 0) {
        rho_up_for_R0 = rho_arr[0];
        p_up_for_R0   = P_arr[0];
      } else if (idx >= N) {
        rho_up_for_R0 = rho_arr[N-1];
        p_up_for_R0   = P_arr[N-1];
      } else {
        double t = (r_up - r_arr[idx-1]) / (r_arr[idx] - r_arr[idx-1]);
        rho_up_for_R0 = rho_arr[idx-1] + t * (rho_arr[idx] - rho_arr[idx-1]);
        p_up_for_R0   = P_arr[idx-1]   + t * (P_arr[idx]   - P_arr[idx-1]);
      }
    }
    // tiny floors
    rho_up_for_R0 = std::max(rho_up_for_R0, 1e-20);
    if (!(E_iso > 0.0) || !(Gam0 > 0.0) || !(rho_up_for_R0 > 0.0)) {
      std::stringstream msg;
      msg << "Bad BM inputs for R0: E_iso=" << E_iso
          << " Gam0=" << Gam0
          << " rho_up_for_R0=" << rho_up_for_R0
          << " (r_up_for_BM=" << r_up << ")";
      ATHENA_ERROR(msg);
    }
    R0 = cbrt( (17.0*E_iso)/(16.0*M_PI*rho_up_for_R0*Gam0*Gam0) );
  } else {
    // Use direct blast radius specified in code units
    if (!(r_blast > 0.0)) {
      std::stringstream msg;
      msg << "Invalid r_blast (must be >0) when use_BM_R0=false: r_blast=" << r_blast;
      ATHENA_ERROR(msg);
    }
    R0 = r_blast;
  }
  if (!std::isfinite(R0) || R0 <= 0.0) {
    std::stringstream msg;
    msg << "Computed/selected R0 is invalid: " << R0
        << (use_BM_R0 ? " (BM mode)" : " (direct r_blast)");
    ATHENA_ERROR(msg);
  }
  // Upstream (unshocked) adiabatic index for enthalpy normalization (defaults to EOS gamma)
  Real Gamma_up = pin->GetOrAddReal("problem", "Gamma_up", 1.333333333);




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

        // Save upstream (unshocked) state for BM enthalpy normalization
        Real rho_up = rho0;
        Real p_up   = pgas0;

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


        Real chi = 1.0 + 8.0*Gam0*Gam0*(1.0 - rad/R0);
        // diagnostics: track chi range
        if (std::isfinite(chi)) {
          if (chi < chi_min_seen) chi_min_seen = chi;
          if (chi > chi_max_seen) chi_max_seen = chi;
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

        if (chi >= 1.0 && chi <= chi_max && angle_ok) {
          bm_cells++;
          // BM profiles
          Real lorentz = (Gam0 / sqrt(2.0)) * pow(chi, -0.5);
          Real beta    = sqrt(1.0 - 1.0/(lorentz*lorentz));



          Real n_up   = rho_up / mp; // local upstream number density
          Real nprime = 4.0 * Gam0 * n_up * pow(chi, -7.0/4.0);     // χ^{-7/4}
          // Upstream enthalpy density w1 = e1 + p1 + rho1 (with c=1); for ideal gas upstream, e1 = p1/(Gamma_up-1)
          Real w1 = rho_up + (Gamma_up/(Gamma_up - 1.0)) * p_up;  // cold upstream -> p_up≈0, so w1≈rho_up
          Real pprime = (2.0/3.0) * Gam0*Gam0 * w1 * pow(chi, -17.0/12.0);

          rho0 = nprime * mp;
          pgas0 = pprime;

          // radial unit vector (relative to blast center)
          Real dx = x - x0, dy = y - y0, dz = z - z0;
          Real invr = (rad>0.0) ? 1.0/rad : 0.0;
          Real erx = invr * dx, ery = invr * dy, erz = invr * dz;

          vx = beta * erx;
          vy = beta * ery;
          vz = beta * erz;
        }
        if (!angle_ok) {
          pgas0 *= cs_factor;
          //rho0 *= cs_factor;
        }
        if (angle_ok && rad > R0) {
          pgas0 *= cs_factor;
          //rho0 *= cs_factor;
        }

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
      "[BM init] R0=%g  mode=%s  dyn=%s  chi_range=[%g,%g]  bm_cells=%lld / %lld  rho[min,max]=[%g,%g]  p[min,max]=[%g,%g]  NaNs: rho=%lld p=%lld v2=%lld\n",
      (double)R0, use_BM_R0 ? "BM" : "r_blast", dyn_mode,
      (double)chi_min_seen, (double)chi_max_seen,
      bm_cells, tot_cells, (double)rho_min, (double)rho_max, (double)p_min,
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