//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_blast.cpp
//! \brief Problem generator for GRMHD spherical blast wave in flat spacetime.

// C headers

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // sqrt()
#include <cstring>    // strcmp()
#include <cstdio>     // fopen(), fprintf(), freopen()   // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <vector>

//Loading in CSV
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
#include "../athena.hpp"                   // macros, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"          // ParameterInput

// Configuration checking
#if not GENERAL_RELATIVITY
#error "This problem generator must be used with general relativity"
#endif

Real threshold;

int RefinementCondition(MeshBlock *pmb);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem","thr");
  }
  return;
}

// Declarations
namespace {
void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
                             Real *px, Real *py, Real *pz);
void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
                     Real *pa0, Real *pa1, Real *pa2, Real *pa3);
Real DistanceBetweenPoints(Real x1, Real x2, Real x3, Real y1, Real y2, Real y3);
} // namespace

//----------------------------------------------------------------------------------------
// Function for setting initial conditions
// Inputs:
//   pin: parameters

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  // Get ratio of specific heats
  // Real gamma_adi = peos->GetGamma();
  // Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read problem parameters
  Real rout =  pin->GetReal("problem", "star_radius");
  Real poly_idx = pin->GetReal("problem", "poly_index");
  std::string poly_csv_path = pin->GetString("problem", "poly_csv_path");
  Real r_blast      = pin->GetReal("problem","r_blast");
  Real blast_factor = pin->GetReal("problem","blast_factor");
  Real dir_angle_min    = pin->GetReal("problem","dir_angle_min");
  Real dir_angle_max    = pin->GetReal("problem","dir_angle_max");
  Real dir_blast_factor = pin->GetReal("problem","dir_blast_factor");
  Real num_x = pin->GetReal("problem", "num_x");
  Real num_y = pin->GetReal("problem", "num_y");
  Real x_spacing = pin->GetReal("problem", "x_spacing");
  Real y_spacing = pin->GetReal("problem", "y_spacing");
  // Real radius = pin->GetReal("problem", "radius");
  // Real rho_inner = pin->GetReal("problem", "rho_inner");
  // Real pgas_inner = pin->GetReal("problem", "pgas_inner");
  // Real rho_outer = pin->GetReal("problem", "rho_outer");
  // Real pgas_outer = pin->GetReal("problem", "pgas_outer");
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // Load precomputed polytrope data from CSV
  static PolytropeData poly;
  if (poly.r.empty()) {
    poly = LoadPolytropeCSV(poly_csv_path);
  }

  Real bx = 0.0, by = 0.0, bz = 0.0;
  if (MAGNETIC_FIELDS_ENABLED) {
    bx = pin->GetReal("problem", "bx");
    by = pin->GetReal("problem", "by");
    bz = pin->GetReal("problem", "bz");
  }

  // Prepare auxiliary arrays
  AthenaArray<Real> b, g, gi;
  b.NewAthenaArray(3, ncells3, ncells2, ncells1);
  g.NewAthenaArray(NMETRIC, ncells1);
  gi.NewAthenaArray(NMETRIC, ncells1);

  // Initialize hydro variables
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      pcoord->CellMetric(k, j, il, iu, g, gi);
      for (int i=il; i<=iu; ++i) {
        // Calculate distance to nearest blast center
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real min_separation = DistanceBetweenPoints(x1, x2, x3, 0.0, 0.0, 0.0);
        for (int x_index = -num_x; x_index <= num_x; ++x_index) {
          Real center_x = x_index * x_spacing;
          for (int y_index = -num_y; y_index <= num_y; ++y_index) {
            if (x_index == 0 && y_index == 0) {
              continue;
            }
            Real center_y = y_index * y_spacing;
            Real separation = DistanceBetweenPoints(x1, x2, x3, center_x, center_y, 0.0);
            min_separation = std::min(min_separation, separation);
          }
        }
        
          Real x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
          Real x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
          Real x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
          Real x0, y0, z0;
          if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
            x0 = x1_0;
            y0 = x2_0;
            z0 = x3_0;
          } else if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0) {
            x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
            y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
            z0 = x1_0*std::cos(x2_0);
          } else {
            // Only check legality of COORDINATE_SYSTEM once in this function
            std::stringstream msg;
           msg << "### FATAL ERROR in gr_ap_blast.cpp ProblemGenerator" << std::endl
                << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
            ATHENA_ERROR(msg);
          }
        Real rad;
        if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0) {
          Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }
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
            rho0  = 1.0e-9;
            pgas0 = 1.0e-9;
          }
        }
        Real theta_dir;
        if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
          // polar angle θ = arccos((z - z0)/rad)
          Real zloc = pcoord->x3v(k) - z0;
          theta_dir = (rad > 0.0 ? std::acos(zloc/ rad) : 0.0);
        } else {  // spherical_polar
          theta_dir = pcoord->x2v(j) - x2_0;
        }

        Real phi_dir;
        if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
          // azimuthal angle φ = atan2(y - y0, x - x0)
          // compute local x,y relative to blast center
          Real xloc = pcoord->x1v(i) - x0;
          Real yloc = pcoord->x2v(j) - y0;
          phi_dir = std::atan2(yloc, xloc);
          if (phi_dir < 0.0) phi_dir += 2.0*M_PI;
        } else {  // spherical_polar
          // φ coordinate in spherical_polar, offset by x3_0
          phi_dir = pcoord->x3v(k) - x3_0;
          phi_dir = std::fmod(phi_dir + 2.0*M_PI, 2.0*M_PI);
        }

        // apply inner blast multiplier: symmetric lobes along x-axis in Cartesian
        if (rr <= r_blast) {
          bool in_lobe = false;
          if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
            // compute local x-y position relative to blast center
            Real xloc = pcoord->x1v(i) - x0;
            Real yloc = pcoord->x2v(j) - y0;
            Real rxy  = std::sqrt(xloc*xloc + yloc*yloc);
            if (rxy > 0.0) {
              Real cosphi = xloc / rxy;  // cos of azimuthal angle
              // lobe if angle between position and +x or -x is within dir_angle_max
              if (std::fabs(cosphi) >= std::cos(dir_angle_max)) {
                in_lobe = true;
              }
            }
          } else {
            // fallback to existing φ-based wedge
            Real cosphi = std::cos(phi_dir);
            if (std::fabs(cosphi) >= std::cos(dir_angle_max)) {
              in_lobe = true;
            }
          }
          if (in_lobe) {
            pgas0 *= dir_blast_factor;
          } else {
            pgas0 *= blast_factor;
          }
        }
        rho  = rho0;
        pgas = pgas0;
        
        // Get Minkowski coordinates of point
        Real t, x, y, z;
        GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);

        // Set velocity
        Real ut = 1.0;
        Real ux = 0.0;
        Real uy = 0.0;
        Real uz = 0.0;
        Real u0, u1, u2, u3;
        TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = u1 - gi(I01,i)/gi(I00,i) * u0;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = u2 - gi(I02,i)/gi(I00,i) * u0;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = u3 - gi(I03,i)/gi(I00,i) * u0;

        // Calculate cell-centered magnetic fields given Minkowski values
        Real bcont = 0.0;
        Real bconx = bx;
        Real bcony = by;
        Real bconz = bz;
        Real bcon0, bcon1, bcon2, bcon3;
        TransformVector(bcont, bconx, bcony, bconz, x, y, z, &bcon0, &bcon1, &bcon2,
                        &bcon3);
        b(IB1,k,j,i) = bcon1 * u0 - bcon0 * u1;
        b(IB2,k,j,i) = bcon2 * u0 - bcon0 * u2;
        b(IB3,k,j,i) = bcon3 * u0 - bcon0 * u3;
      }
    }
  }
  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  // Delete auxiliary array

  // Initialize magnetic field
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
        for (int i=il; i<=iu+1; ++i) {
          Real ut = 1.0;
          Real ux = 0.0;
          Real uy = 0.0;
          Real uz = 0.0;
          Real bcont = 0.0;
          Real bconx = bx;
          Real bcony = by;
          Real bconz = bz;
          Real u0, u1, u2, u3;
          Real bcon0, bcon1, bcon2, bcon3;
          if (j != ju+1 && k != ku+1) {
            Real x1 = pcoord->x1f(i);
            Real x2 = pcoord->x2v(j);
            Real x3 = pcoord->x3v(k);
            Real t, x, y, z;
            GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
            TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
            TransformVector(bcont, bconx, bcony, bconz, x, y, z, &bcon0, &bcon1, &bcon2,
                            &bcon3);
            pfield->b.x1f(k,j,i) = bcon1 * u0 - bcon0 * u1;
          }
          if (i != iu+1 && k != ku+1) {
            Real x1 = pcoord->x1v(i);
            Real x2 = pcoord->x2f(j);
            Real x3 = pcoord->x3v(k);
            Real t, x, y, z;
            GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
            TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
            TransformVector(bcont, bconx, bcony, bconz, x, y, z, &bcon0, &bcon1, &bcon2,
                            &bcon3);
            pfield->b.x2f(k,j,i) = bcon2 * u0 - bcon0 * u2;
          }
          if (i != iu+1 && j != ju+1) {
            Real x1 = pcoord->x1v(i);
            Real x2 = pcoord->x2v(j);
            Real x3 = pcoord->x3f(k);
            Real t, x, y, z;
            GetMinkowskiCoordinates(0.0, x1, x2, x3, &t, &x, &y, &z);
            TransformVector(ut, ux, uy, uz, x, y, z, &u0, &u1, &u2, &u3);
            TransformVector(bcont, bconx, bcony, bconz, x, y, z, &bcon0, &bcon1, &bcon2,
                            &bcon3);
            pfield->b.x3f(k,j,i) = bcon3 * u0 - bcon0 * u3;
          }
        }
      }
    }
  }
  return;
}

namespace {
//----------------------------------------------------------------------------------------
// Function for returning corresponding Minkowski coordinates of point
// Inputs:
//   x0,x1,x2,x3: global coordinates to be converted
// Outputs:
//   pt,px,py,pz: variables pointed to set to Minkowski coordinates
// Notes:
//   conversion is trivial
//   useful to have if other coordinate systems for Minkowski space are developed

void GetMinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
                             Real *px, Real *py, Real *pz) {
  if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
    *pt = x0;
    *px = x1;
    *py = x2;
    *pz = x3;
  } else if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0) {
    // For m=0 Schwarzschild: treat (x1, x2, x3) as spherical (r,θ,φ)
    *pt = x0;
    *px = x1 * std::sin(x2) * std::cos(x3);
    *py = x1 * std::sin(x2) * std::sin(x3);
    *pz = x1 * std::cos(x2);
  }
  return;
}


//----------------------------------------------------------------------------------------
// Function for transforming 4-vector from Minkowski to desired coordinates
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   x,y,z: Minkowski coordinates of point
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in desired coordinates
// Notes:
//   conversion is trivial
//   useful to have if other coordinate systems for Minkowski space are developed

void TransformVector(Real at, Real ax, Real ay, Real az, Real x, Real y, Real z,
                     Real *pa0, Real *pa1, Real *pa2, Real *pa3) {
  if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
    *pa0 = at;
    *pa1 = ax;
    *pa2 = ay;
    *pa3 = az;
  } else if (std::strcmp(COORDINATE_SYSTEM, "schwarzschild") == 0) {
    // No additional rotation for flat Schwarzschild (m=0)
    *pa0 = at;
    *pa1 = ax;
    *pa2 = ay;
    *pa3 = az;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for returning spatial separation between points at same time
// Inputs:
//   x1,x2,x3: spatial coordinates of one point
//   y1,y2,y3: spatial coordinates of other point
// Outputs:
//   returned value: spatial separation between x and y
// Notes:
//   distance function is Euclidean in Minkowski coordinates

Real DistanceBetweenPoints(Real x1, Real x2, Real x3, Real y1, Real y2, Real y3) {
  Real distance = 0.0;
  if (std::strcmp(COORDINATE_SYSTEM, "minkowski") == 0) {
    distance = std::sqrt(SQR(x1-y1) + SQR(x2-y2) + SQR(x3-y3));
  }
  return distance;
}
} // namespace
void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // delay AMR enrollment until after first timestep
  static bool enrolled_amr = false;
  if (adaptive && !enrolled_amr && time > 0.01) {
    EnrollUserRefinementCondition(RefinementCondition);
    enrolled_amr = true;
  }
}
// refinement condition: check the maximum pressure gradient
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  // compute block center
  auto &coord = *pmb->pcoord;
  Real rcen = coord.x1v(pmb->is) + coord.x1v(pmb->ie);
  rcen *= 0.5;
  if (rcen > 3.0) return 0;   // no AMR out beyond 3.0R

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
  // decide refine/derefine based on curvature threshold
  if (maxcurv > threshold) return 1;
  if (maxcurv < 0.25*threshold) return -1;
  return 0;
}