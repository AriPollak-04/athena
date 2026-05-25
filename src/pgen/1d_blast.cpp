//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file 1d_blast.cpp
//! \brief 1D blast wave in a polytropic stellar profile.
//!        Cartesian 1D (nx2=nx3=1); supports SR and Newtonian dynamics.
//!        A pressure bomb at x < bomb_radius drives a rightward shock through a
//!        polytropic star (0 <= x <= star_radius) into an ambient medium.
//!        An optional uniform transverse velocity vy_bulk (collapsed y-direction)
//!        is imposed on all cells and contributes to the Lorentz factor.
//!
//! Build: python configure.py --prob=1d_blast --coord=cartesian --eos=adiabatic \
//!                             --flux=hlle -s --mpi --hdf5

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct PolytropeData {
  std::vector<double> r, rho, m, P;
};

static PolytropeData LoadPolytropeCSV(const std::string &filename) {
  std::ifstream in(filename);
  if (!in.is_open())
    throw std::runtime_error("Failed to open " + filename);
  std::string line;
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
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

int RefinementCondition(MeshBlock *pmb);

static Real threshold;
static Real vy_bulk_g      = 0.0;
static Real bomb_radius_g  = 0.05;
static Real bomb_p_ratio_g = 1.0e4;
static Real ambient_rho_g  = 1.0e-25;
static Real ambient_p_g    = 1.0e-25;

// ---- User history outputs ----

static Real Hst_Etot(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += u(IEN,k,j,i) * pco->GetCellVolume(k,j,i);
  return sum;
}

static Real Hst_Eexcess(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += (u(IEN,k,j,i) - u(IDN,k,j,i)) * pco->GetCellVolume(k,j,i);
  return sum;
}

static Real Hst_Px(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += u(IM1,k,j,i) * pco->GetCellVolume(k,j,i);
  return sum;
}

static Real Hst_Py(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &u = pmb->phydro->u;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += u(IM2,k,j,i) * pco->GetCellVolume(k,j,i);
  return sum;
}

static Real Hst_PgasInt(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  AthenaArray<Real> &w = pmb->phydro->w;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += w(IPR,k,j,i) * pco->GetCellVolume(k,j,i);
  return sum;
}

static Real Hst_Volume(MeshBlock *pmb, int iout) {
  Coordinates *pco = pmb->pcoord;
  Real sum = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        sum += pco->GetCellVolume(k,j,i);
  return sum;
}

//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem", "thr");
  }

  bomb_radius_g  = pin->GetOrAddReal("problem", "bomb_radius",  0.05);
  bomb_p_ratio_g = pin->GetOrAddReal("problem", "bomb_p_ratio", 1.0e4);
  ambient_rho_g  = pin->GetOrAddReal("problem", "ambient_rho",  1.0e-25);
  ambient_p_g    = pin->GetOrAddReal("problem", "ambient_p",    1.0e-25);

  // Accept either vy_bulk directly, or Gam_bulk (Lorentz factor) from which vy is derived.
  // vy_bulk takes priority if non-zero; otherwise Gam_bulk > 1 is used.
  vy_bulk_g = pin->GetOrAddReal("problem", "vy_bulk", 0.0);
  if (vy_bulk_g == 0.0) {
    Real Gam_bulk = pin->GetOrAddReal("problem", "Gam_bulk", 1.0);
    if (Gam_bulk > 1.0)
      vy_bulk_g = std::sqrt(std::max(0.0, 1.0 - 1.0/(Gam_bulk*Gam_bulk)));
  }

  if (Globals::my_rank == 0) {
    std::fprintf(stderr,
      "[1d_blast:init] vy_bulk=%g  bomb_radius=%g  bomb_p_ratio=%g\n",
      (double)vy_bulk_g, (double)bomb_radius_g, (double)bomb_p_ratio_g);
  }

  AllocateUserHistoryOutput(6);
  EnrollUserHistoryOutput(0, Hst_Etot,    "E_tot");
  EnrollUserHistoryOutput(1, Hst_Eexcess, "E_excess");
  EnrollUserHistoryOutput(2, Hst_Px,      "Px_tot");
  EnrollUserHistoryOutput(3, Hst_Py,      "Py_tot");
  EnrollUserHistoryOutput(4, Hst_PgasInt, "Pgas_int");
  EnrollUserHistoryOutput(5, Hst_Volume,  "V_tot");
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief 1D blast wave in a polytropic star with optional transverse velocity
//========================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const Real rout          = pin->GetReal("problem", "star_radius");
  const std::string csv    = pin->GetString("problem", "poly_csv_path");
  const Real cs_factor     = pin->GetReal("problem", "cs_factor");
  const Real bomb_radius   = bomb_radius_g;
  const Real bomb_p_ratio  = bomb_p_ratio_g;
  const Real vy_bulk       = vy_bulk_g;
  const Real ambient_rho   = ambient_rho_g;
  const Real ambient_p     = ambient_p_g;

  static PolytropeData poly;
  if (poly.r.empty())
    poly = LoadPolytropeCSV(csv);

  const Real gamma = peos->GetGamma();
  const Real gm1   = gamma - 1.0;

  long long tot_cells = 0, nan_rho = 0, nan_p = 0;
  Real rho_min = 1e99, rho_max = -1e99, p_min = 1e99, p_max = -1e99;

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        const Real x = pcoord->x1v(i);

        Real rho, pgas;
        if (x <= rout) {
          const auto &r_arr   = poly.r;
          const auto &rho_arr = poly.rho;
          const auto &P_arr   = poly.P;
          const int N = static_cast<int>(r_arr.size());
          auto it  = std::lower_bound(r_arr.begin(), r_arr.end(), x);
          int idx  = static_cast<int>(std::distance(r_arr.begin(), it));
          if (idx <= 0) {
            rho  = rho_arr[0];
            pgas = P_arr[0];
          } else if (idx >= N) {
            rho  = rho_arr[N-1];
            pgas = P_arr[N-1];
          } else {
            double t = (x - r_arr[idx-1]) / (r_arr[idx] - r_arr[idx-1]);
            rho  = rho_arr[idx-1] + t*(rho_arr[idx]  - rho_arr[idx-1]);
            pgas = P_arr[idx-1]   + t*(P_arr[idx]    - P_arr[idx-1]);
          }
        } else {
          rho  = ambient_rho;
          pgas = ambient_p;
        }

        pgas *= cs_factor;

        // Pressure bomb: over-pressurize a region near x=0 to seed the blast
        if (x <= bomb_radius)
          pgas *= bomb_p_ratio;

        if (!std::isfinite(rho)  || rho  <= 0.0) rho  = 1e-30;
        if (!std::isfinite(pgas) || pgas <= 0.0) pgas = 1e-30;

        // vx=0 everywhere initially; transverse vy_bulk applied uniformly
        const Real vx = 0.0, vy = vy_bulk, vz = 0.0;

        tot_cells++;
        if (!std::isfinite(rho))  nan_rho++;
        if (!std::isfinite(pgas)) nan_p++;
        rho_min = std::min(rho_min, rho);
        rho_max = std::max(rho_max, rho);
        p_min   = std::min(p_min,   pgas);
        p_max   = std::max(p_max,   pgas);

        phydro->w(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = pgas;
        phydro->w(IVX,k,j,i) = vx;
        phydro->w(IVY,k,j,i) = vy;
        phydro->w(IVZ,k,j,i) = vz;

#if defined(RELATIVISTIC_DYNAMICS) && (RELATIVISTIC_DYNAMICS != 0)
        Real v2 = vx*vx + vy*vy + vz*vz;
        v2 = std::min(v2, 1.0 - 1e-12);
        const Real gL    = 1.0 / std::sqrt(1.0 - v2);
        const Real rho_s = std::max(rho, (Real)1e-30);
        const Real h     = 1.0 + (gamma / gm1) * (pgas / rho_s);
        phydro->u(IDN,k,j,i) = rho_s * gL;
        phydro->u(IM1,k,j,i) = rho_s * h * gL*gL * vx;
        phydro->u(IM2,k,j,i) = rho_s * h * gL*gL * vy;
        phydro->u(IM3,k,j,i) = rho_s * h * gL*gL * vz;
        Real Etot = rho_s * h * gL*gL - pgas;
        phydro->u(IEN,k,j,i) = std::max(Etot, (Real)1e-30);
#else
        Real v2 = vx*vx + vy*vy + vz*vz;
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho * vx;
        phydro->u(IM2,k,j,i) = rho * vy;
        phydro->u(IM3,k,j,i) = rho * vz;
        phydro->u(IEN,k,j,i) = pgas / gm1 + 0.5*rho*v2;
#endif
      }
    }
  }

  if (Globals::my_rank == 0) {
    std::fprintf(stderr,
      "[1d_blast:PG] rho[min,max]=[%g,%g]  p[min,max]=[%g,%g]  NaNs: rho=%lld p=%lld\n",
      (double)rho_min, (double)rho_max,
      (double)p_min,   (double)p_max, nan_rho, nan_p);
  }
}

//========================================================================================
//! \fn int RefinementCondition(MeshBlock *pmb)
//! \brief AMR based on normalized log-gradient of density and pressure
//========================================================================================
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  const Real tiny    = 1e-99;
  const Real pfloor  = 1e-20;
  const Real rhofloor= 1e-30;

  Real max_grad = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        const Real rho = std::max(w(IDN,k,j,i), rhofloor);
        const Real p   = std::max(w(IPR,k,j,i), pfloor);

        // Central (or one-sided at edges) finite difference in x1 only
        Real drho_dx, dp_dx, dx;
        if (i == pmb->is) {
          const Real h = pmb->pcoord->x1v(i+1) - pmb->pcoord->x1v(i) + tiny;
          drho_dx = (w(IDN,k,j,i+1) - w(IDN,k,j,i)) / h;
          dp_dx   = (w(IPR,k,j,i+1) - w(IPR,k,j,i)) / h;
          dx = h;
        } else if (i == pmb->ie) {
          const Real h = pmb->pcoord->x1v(i) - pmb->pcoord->x1v(i-1) + tiny;
          drho_dx = (w(IDN,k,j,i) - w(IDN,k,j,i-1)) / h;
          dp_dx   = (w(IPR,k,j,i) - w(IPR,k,j,i-1)) / h;
          dx = h;
        } else {
          const Real h = pmb->pcoord->x1v(i+1) - pmb->pcoord->x1v(i-1) + tiny;
          drho_dx = (w(IDN,k,j,i+1) - w(IDN,k,j,i-1)) / h;
          dp_dx   = (w(IPR,k,j,i+1) - w(IPR,k,j,i-1)) / h;
          dx = 0.5*h;
        }

        // Dimensionless: |d ln rho / dx| * dx  +  |d ln p / dx| * dx
        const Real gr = (std::fabs(drho_dx)/rho + std::fabs(dp_dx)/p) * dx;
        if (std::isfinite(gr) && gr > max_grad) max_grad = gr;
      }
    }
  }

  if (max_grad > threshold)        return  1;
  if (max_grad < 0.9*threshold)    return -1;
  return 0;
}
