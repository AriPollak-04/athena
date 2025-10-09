//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file jet.cpp
//! \brief Sets up a nonrelativistic jet introduced through L-x1 boundary (left edge)
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>


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
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

// BCs on L-x1 (left edge) of grid with jet inflow conditions
void JetInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
// Make radius of jet and jet variables global so they can be accessed by BC functions
// Real r_amb,
Real d_amb, p_amb, vx_amb, vy_amb, vz_amb, bx_amb, by_amb, bz_amb;
Real r_jet, d_jet, p_jet, vx_jet, vy_jet, vz_jet, bx_jet, by_jet, bz_jet;
Real gm1, x2_0, x3_0;
Real t_stop;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
  d_amb  = pin->GetReal("problem", "d");
  p_amb  = pin->GetReal("problem", "p");
  vx_amb = pin->GetReal("problem", "vx");
  vy_amb = pin->GetReal("problem", "vy");
  vz_amb = pin->GetReal("problem", "vz");
  if (MAGNETIC_FIELDS_ENABLED) {
    bx_amb = pin->GetReal("problem", "bx");
    by_amb = pin->GetReal("problem", "by");
    bz_amb = pin->GetReal("problem", "bz");
  }
  d_jet  = pin->GetReal("problem", "djet");
  p_jet  = pin->GetReal("problem", "pjet");
  vx_jet = pin->GetReal("problem", "vxjet");
  vy_jet = pin->GetReal("problem", "vyjet");
  vz_jet = pin->GetReal("problem", "vzjet");
  t_stop = pin->GetReal("problem", "t_stop");
  if (MAGNETIC_FIELDS_ENABLED) {
    bx_jet = pin->GetReal("problem", "bxjet");
    by_jet = pin->GetReal("problem", "byjet");
    bz_jet = pin->GetReal("problem", "bzjet");
  }
  r_jet = pin->GetReal("problem", "rjet");
  x2_0 = 0.5*(mesh_size.x2max + mesh_size.x2min);
  x3_0 = 0.5*(mesh_size.x3max + mesh_size.x3min);

  // enroll boundary value function pointers
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, JetInnerX1);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Jet problem
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gm1 = peos->GetGamma() - 1.0;
  std::string poly_csv_path = pin->GetString("problem", "poly_csv_path");
  Real rout = pin->GetReal("problem", "rout");
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

  // initialize conserved variables
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
                Real rad;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
          Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
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
            rho0  = d_amb;
            pgas0 = p_amb;
          }
        }
        phydro->u(IDN,k,j,i) = rho0;
        phydro->u(IM1,k,j,i) = rho0*vx_amb;
        phydro->u(IM2,k,j,i) = rho0*vy_amb;
        phydro->u(IM3,k,j,i) = rho0*vz_amb;
        if (NON_BAROTROPIC_EOS) {
          // Use local (possibly stellar) gas pressure and density for internal + kinetic energy
          phydro->u(IEN,k,j,i) = pgas0/gm1
                                 + 0.5*rho0*(SQR(vx_amb)+SQR(vy_amb)+SQR(vz_amb));
        }
      }
    }
  }

  // initialize interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = bx_amb;
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = by_amb;
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = bz_amb;
        }
      }
    }
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            phydro->u(IEN,k,j,i) += 0.5*(SQR(bx_amb) + SQR(by_amb) + SQR(bz_amb));
          }
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void JetInnerX1()
//  \brief Sets boundary condition on left X boundary (iib) for jet problem

void JetInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  const bool jet_active = (t_stop <= 0.0) || (time <= t_stop);
  // set primitive variables in inlet ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        Real rad = std::sqrt(SQR(pco->x2v(j)-x2_0) + SQR(pco->x3v(k)-x3_0));
        const bool in_nozzle = (rad<=r_jet);
        if (in_nozzle && jet_active) {
          prim(IDN,k,j,il-i) = d_jet;
          prim(IVX,k,j,il-i) = vx_jet;
          prim(IVY,k,j,il-i) = vy_jet;
          prim(IVZ,k,j,il-i) = vz_jet;
          prim(IPR,k,j,il-i) = p_jet;
        } else {
          prim(IDN,k,j,il-i) = prim(IDN,k,j,il);
          prim(IVX,k,j,il-i) = prim(IVX,k,j,il);
          prim(IVY,k,j,il-i) = prim(IVY,k,j,il);
          prim(IVZ,k,j,il-i) = prim(IVZ,k,j,il);
          prim(IPR,k,j,il-i) = prim(IPR,k,j,il);
        }
      }
    }
  }

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          Real rad = std::sqrt(SQR(pco->x2v(j)-x2_0) + SQR(pco->x3v(k)-x3_0));
          if (rad <= r_jet) {
            b.x1f(k,j,il-i) = bx_jet;
          } else {
            b.x1f(k,j,il-i) = b.x1f(k,j,il);
          }
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          Real rad = std::sqrt(SQR(pco->x2v(j)-x2_0) + SQR(pco->x3v(k)-x3_0));
          if (rad <= r_jet) {
            b.x2f(k,j,il-i) = by_jet;
          } else {
            b.x2f(k,j,il-i) = b.x2f(k,j,il);
          }
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          Real rad = std::sqrt(SQR(pco->x2v(j)-x2_0) + SQR(pco->x3v(k)-x3_0));
          if (rad <= r_jet) {
            b.x3f(k,j,il-i) = bz_jet;
          } else {
            b.x3f(k,j,il-i) = b.x3f(k,j,il);
          }
        }
      }
    }
  }
}
