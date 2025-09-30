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
// cache ParameterInput pointer for use outside ProblemGenerator
static ParameterInput *g_pin = nullptr;

int RefinementCondition(MeshBlock *pmb);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem","thr");
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  g_pin = pin;
  // reference to primitive array for velocity assignments
  auto &w = phydro->w;
  Real rout =  pin->GetReal("problem", "star_radius");
  Real poly_idx = pin->GetReal("problem", "poly_index");
  std::string poly_csv_path = pin->GetString("problem", "poly_csv_path");
  Real r_blast      = pin->GetReal("problem","r_blast");
  Real blast_factor = pin->GetReal("problem","blast_factor");
  Real dir_angle_min    = pin->GetReal("problem","dir_angle_min");
  Real dir_angle_max    = pin->GetReal("problem","dir_angle_max");
  Real dir_blast_factor = pin->GetReal("problem","dir_blast_factor");
  Real vel_perturb = pin->GetOrAddReal("problem","vel_perturb",0.1);
  Real b0, angle;
  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = pin->GetReal("problem", "b0");
    angle = (PI/180.0)*pin->GetReal("problem", "angle");
  }
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // --- debug tables setup ---
  Real r_fixed      = pin->GetOrAddReal("problem","r_fixed",      0.5*r_blast);
  Real theta_fixed  = pin->GetOrAddReal("problem","theta_fixed", 3.14159265 );
  Real phi_fixed    = pin->GetOrAddReal("problem","phi_fixed",    3.14159265);
  std::vector<std::pair<Real,Real>> phi_table;   // (phi, pgas0)
  std::vector<std::pair<Real,Real>> theta_table; // (theta, pgas0)
  // --- end debug tables setup ---

  // Load precomputed polytrope data from CSV
  static PolytropeData poly;
  if (poly.r.empty()) {
    poly = LoadPolytropeCSV(poly_csv_path);
  }

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

  // --- debug: generic sampling points (rr,theta,phi) ---
  struct DebugPoint { Real rr, theta, phi; bool printed; };
  static std::vector<DebugPoint> dbg_pts;
  static bool dbg_init = false;
  const Real tol_rr    = r_blast/20.0;
  const Real tol_angle = M_PI/20.0;
  const int  NRR       = 5, NANG = 5;
  if (!dbg_init) {
    // sample 5 radii between 0 and r_blast
    for (int n=0; n<NRR; n++) {
      dbg_pts.push_back({ (n+1)*r_blast/Real(NRR+1), 0.0, 0.0, false });
    }
    // sample NANG angles at fixed theta, varying phi
    Real theta_fixed = pin->GetOrAddReal("problem","theta_fixed", M_PI/2.0);
    for (int n=0; n<NANG; n++) {
      Real phi = n*2.0*M_PI/Real(NANG);
      dbg_pts.push_back({ 0.0, theta_fixed, phi, false });
    }
    dbg_init = true;
  }
  // --- end debug init ---

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real rad;
        // compute Cartesian offsets about blast center once
        Real xloc = 0.0, yloc = 0.0, zloc = 0.0;
        if (std::strcmp(COORDINATE_SYSTEM,"cartesian")==0) {
          xloc = pcoord->x1v(i) - x0;
          yloc = pcoord->x2v(j) - y0;
          zloc = pcoord->x3v(k) - z0;
        } else if (std::strcmp(COORDINATE_SYSTEM,"cylindrical")==0) {
          xloc = pcoord->x1v(i)*std::cos(pcoord->x2v(j)) - x0;
          yloc = pcoord->x1v(i)*std::sin(pcoord->x2v(j)) - y0;
          zloc = pcoord->x3v(k)              - z0;
        } else { // spherical_polar
          xloc = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k)) - x0;
          yloc = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k)) - y0;
          zloc = pcoord->x1v(i)*std::cos(pcoord->x2v(j))                         - z0;
        }
        Real rloc = std::sqrt(xloc*xloc + yloc*yloc + zloc*zloc);
        rad = rloc;
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
                // collect debug entries at fixed radius
        Real theta_dir;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          // polar angle θ = arccos(zloc/rloc)
          theta_dir = (rloc > 0.0 ? std::acos(zloc/ rloc) : 0.0);
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          // in cylindrical axisym, x2v(j) is the azimuthal angle φ, so θ = atan2(|zloc|, r)
          Real rr_cyl = pcoord->x1v(i);
          theta_dir = (rr_cyl > 0.0 ? std::atan2(std::abs(zloc), rr_cyl) : 0.0);
        } else {  // spherical_polar
          theta_dir = pcoord->x2v(j) - x2_0;
        }

        Real phi_dir;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          // azimuthal angle φ = atan2(yloc, xloc)
          phi_dir = std::atan2(yloc, xloc);
          if (phi_dir < 0.0) phi_dir += 2.0*M_PI;
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          // φ coordinate in cylindrical, offset by x3_0
          phi_dir = pcoord->x2v(j) - x3_0;
          phi_dir = std::fmod(phi_dir + 2.0*M_PI, 2.0*M_PI);
        } else {  // spherical_polar
          // φ coordinate in spherical_polar, offset by x3_0
          phi_dir = pcoord->x3v(k) - x3_0;
          phi_dir = std::fmod(phi_dir + 2.0*M_PI, 2.0*M_PI);
        }

        // apply inner blast: set initial radial velocity inside r_blast
        if (rr <= r_blast) {
          // compute in_lobe once for both Newtonian and relativistic
          bool in_lobe = false;
          {
            Real rxy = std::sqrt(xloc*xloc + yloc*yloc);
            if (rxy > 0.0) {
              Real cosphi = xloc / rxy;
              if (std::fabs(cosphi) >= std::cos(dir_angle_max)) {
                in_lobe = true;
              }
            }
          }
#ifndef RELATIVISTIC_DYNAMICS
          // assign initial radial velocity to all cells inside r_blast
          {
            Real v0 = pin->GetOrAddReal("problem","v_blast", 0.1);
            if (rloc > 0.0) {
              w(IVX,k,j,i) = v0 * xloc/rloc;
              w(IVY,k,j,i) = v0 * yloc/rloc;
              w(IVZ,k,j,i) = v0 * zloc/rloc;
            } else {
              w(IVX,k,j,i) = 0.0;
              w(IVY,k,j,i) = 0.0;
              w(IVZ,k,j,i) = 0.0;
            }
          }
#else
          // Relativistic: initial radial velocity injection inside r_star in lobes
          {
            // read desired blast velocity (default 0.1c)
            Real v0 = pin->GetOrAddReal("problem","v_blast", 0.1);
            if (rr <= r_star) {
              if (in_lobe && rloc > 0.0) {
                // set radial velocity in direction of (xloc,yloc,zloc)
                Real vr = v0;
                w(IVX,k,j,i) = vr * xloc/rloc;
                w(IVY,k,j,i) = vr * yloc/rloc;
                w(IVZ,k,j,i) = vr * zloc/rloc;
              } else {
                // zero velocity outside lobes
                w(IVX,k,j,i) = 0.0;
                w(IVY,k,j,i) = 0.0;
                w(IVZ,k,j,i) = 0.0;
              }
            }
          }
#endif
        }

        // unified debug: print once for each (rr,theta,phi) target
        if (Globals::my_rank == 0) {
          for (auto &dp : dbg_pts) {
            bool match_rr    = (dp.rr > 0.0 && std::fabs(rr - dp.rr) < tol_rr) ||
                               (dp.rr == 0.0 && rr < r_blast);
            bool match_theta = std::fabs(theta_dir - dp.theta) < tol_angle;
            bool match_phi   = std::fabs(phi_dir   - dp.phi)   < tol_angle;
            if (!dp.printed && match_rr && match_theta && match_phi) {
              printf("DBG_PT: rr=%.6f θ=%.6f φ=%.6f -> pgas0=%e\n",
                     dp.rr, dp.theta, dp.phi, pgas0);
              fflush(stdout);
              dp.printed = true;
            }
          }
        }

        rho  = rho0;
        pgas = pgas0;
        phydro->u(IDN,k,j,i) = rho;
#ifndef RELATIVISTIC_DYNAMICS
        // zero momentum in Newtonian mode only
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
#endif
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = pgas/(gamma-1.0);
          if (RELATIVISTIC_DYNAMICS) phydro->u(IEN,k,j,i) += rho;
        }
      }
    }
  }

  // --- debug tables printout ---
  static bool tables_printed = false;
  if (Globals::my_rank == 0 && !tables_printed) {
    printf("DEBUG TABLE: phi vs pgas at r=%.6f, theta=%.6f\n", r_fixed, theta_fixed);
    printf("   phi          pgas\n");
    for (auto &p : phi_table) {
      printf(" % .6f   % .6e\n", p.first, p.second);
    }
    printf("DEBUG TABLE: theta vs pgas at r=%.6f, phi=%.6f\n", r_fixed, phi_fixed);
    printf("  theta         pgas\n");
    for (auto &p : theta_table) {
      printf(" % .6f   % .6e\n", p.first, p.second);
    }
    fflush(stdout);
    tables_printed = true;
  }
  // --- end debug tables printout ---

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


// Damp velocities outside the angular lobe only for relativistic runs
#include <cstdio>
void MeshBlock::UserWorkInLoop() {
#ifndef RELATIVISTIC_DYNAMICS
  return; // do nothing in Newtonian runs
#else
  // read (once) the parameters we need
  static bool inited = false;
  static Real ang_min, ang_max;
  static Real clamp_alpha;
  static Real x1_0, x2_0, x3_0;
  static Real x0, y0, z0;
  static Real time_zero, time_end, r_zero;
  if (!inited) {
    ParameterInput *pin = g_pin;
    if (pin == nullptr) return;  // safety
    ang_min     = pin->GetReal("problem","dir_angle_min");
    ang_max     = pin->GetReal("problem","dir_angle_max");
    clamp_alpha = pin->GetOrAddReal("problem","clamp_alpha", 0.0); // 0 => kill, 1 => keep
    time_zero   = pin->GetOrAddReal("problem","time_zero", 0.1);
    r_zero      = pin->GetOrAddReal("problem","r_zero",    0.5);
    time_end    = pin->GetOrAddReal("problem","time_end", time_zero + 0.1);
    x1_0 = pin->GetOrAddReal("problem","x1_0",0.0);
    x2_0 = pin->GetOrAddReal("problem","x2_0",0.0);
    x3_0 = pin->GetOrAddReal("problem","x3_0",0.0);
    // convert center to Cartesian once
    if (std::strcmp(COORDINATE_SYSTEM,"cartesian")==0) {
      x0 = x1_0; y0 = x2_0; z0 = x3_0;
    } else if (std::strcmp(COORDINATE_SYSTEM,"cylindrical")==0) {
      x0 = x1_0*std::cos(x2_0);
      y0 = x1_0*std::sin(x2_0);
      z0 = x3_0;
    } else { // spherical_polar
      x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
      y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
      z0 = x1_0*std::cos(x2_0);
    }
    inited = true;
  }

  auto &w = phydro->w;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // Cartesian offsets
        Real xloc, yloc, zloc;
        if (std::strcmp(COORDINATE_SYSTEM,"cartesian")==0) {
          xloc = pcoord->x1v(i)-x0;
          yloc = pcoord->x2v(j)-y0;
          zloc = pcoord->x3v(k)-z0;
        } else if (std::strcmp(COORDINATE_SYSTEM,"cylindrical")==0) {
          Real r = pcoord->x1v(i), phi = pcoord->x2v(j);
          xloc = r*std::cos(phi)-x0;
          yloc = r*std::sin(phi)-y0;
          zloc = pcoord->x3v(k)-z0;
        } else { // spherical_polar
          Real r = pcoord->x1v(i);
          Real th = pcoord->x2v(j);
          Real ph = pcoord->x3v(k);
          xloc = r*std::sin(th)*std::cos(ph)-x0;
          yloc = r*std::sin(th)*std::sin(ph)-y0;
          zloc = r*std::cos(th)                -z0;
        }

        Real rloc = std::sqrt(xloc*xloc + yloc*yloc + zloc*zloc);

        // wedge test (same as in ProblemGenerator)
        Real rxy = std::sqrt(xloc*xloc + yloc*yloc);
        bool in_lobe = false;
        if (rxy > 0.0) {
          Real cosphi = xloc / rxy;
          if (std::fabs(cosphi) >= std::cos(ang_max)) in_lobe = true;
        }

        // after time_zero and before time_end and within r_zero, zero velocities outside lobes
        if (pmy_mesh->time > time_zero && rloc < r_zero && !in_lobe) {
          w(IVX,k,j,i) = 0.0;
          w(IVY,k,j,i) = 0.0;
          w(IVZ,k,j,i) = 0.0;
        }
      }
    }
  }
#endif
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
  if (rcen > 6.0) return 0;   // no AMR out beyond 6.0R

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
    // decide refine/derefine based on curvature threshold
  if (maxcurv > threshold) return 1;
  if (maxcurv < 0.25*threshold) return -1;
  return 0;
}
//  if (maxcurv > threshold) return 1;
//  if (maxcurv < 0.25*threshold) return -1;
//  return 0;
// }
//  if (pmb->pmy_mesh->f3) {
//    for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
//      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
//        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
//          Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
//                               +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
//                               +SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
//          maxeps = std::max(maxeps, eps);
//        }
//      }
//    }
//  } else if (pmb->pmy_mesh->f2) {
//    int k = pmb->ks;
//    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
//      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
//        Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
//                             + SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i))))/w(IPR,k,j,i);
//        maxeps = std::max(maxeps, eps);
//      }
//    }
// grad(logarithm of density/pressure)