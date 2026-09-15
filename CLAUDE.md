# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Working Style

Do not create planning, analysis, or scratch markdown files during tasks. Keep all notes and findings inline in the conversation.

## What This Is

Athena++ is an astrophysical MHD and radiation code with adaptive mesh refinement (AMR). The current local work focuses on jet-blast simulations using cylindrical coordinates with special relativity. The codebase is C++11 with Python tooling for configuration and analysis.

## Build System

Configuration generates a custom `Makefile` and `src/defs.hpp` from templates:

```bash
python configure.py --prob=jet_blast --coord=cylindrical --eos=adiabatic --flux=hlle -s --mpi --hdf5
make -j 192
```

The current build is already configured (see `configure.log`). After changing problem generators or physics modules, re-run `configure.py` with the appropriate flags before `make`.

Key `configure.py` flags:
- `--prob=<name>` — problem generator from `src/pgen/`
- `--coord=<system>` — cartesian, cylindrical, spherical_polar, schwarzschild, kerr-schild
- `--eos=<type>` — adiabatic, isothermal
- `-b` — enable magnetic fields; `-s` — special relativity; `-g` — general relativity
- `-mpi`, `-hdf5`, `-omp` — enable MPI, HDF5, OpenMP

Run a simulation with MPI, switch first into output directory, then copy input file in, then run executable ~/athena/bin/athena:
```bash
cd /scratch/aripoll/athena_out/outputs 
cp /scratch/aripoll/athena/inputs/mhd/athinput.jet_blast .
mpiexec -n 192 /scratch/aripoll/athena/bin/athena -i athinput.jet_blast
```

Output goes to `/scratch/aripoll/athena_out/` by default (set via input file).

## Architecture

### Physics Pipeline

The simulation loop lives in `src/main.cpp`. Each timestep, the `TaskList` system (`src/task_list/`) schedules work across `MeshBlock`s. The core physics chain per block is:

1. **Reconstruction** (`src/reconstruct/`) — interpolates cell-centered to face-centered values (PPM, linear, etc.)
2. **Riemann solver** (`src/hydro/`) — computes fluxes at cell faces
3. **Integration** — updates conserved variables; boundary communication via `src/bvals/`
4. **EOS** (`src/eos/`) — converts between conserved and primitive variables

### Key Abstractions

- **`Mesh`** / **`MeshBlock`** (`src/mesh/`) — spatial domain decomposition; AMR logic lives here. Each `MeshBlock` owns its data arrays.
- **`AthenaArray`** (`src/athena_arrays.hpp`) — custom ND array used everywhere for mesh data.
- **Problem generators** (`src/pgen/`) — define initial and boundary conditions; `jet_blast.cpp` reads a polytrope CSV and injects a jet.
- **Input files** (`.athinput` format, `src/parameter_input.cpp`) — control all runtime parameters: mesh size, boundary conditions, physics module settings, problem-specific values.

### Current Problem: `jet_blast`

`src/pgen/jet_blast.cpp` loads a polytropic stellar profile (CSV from `getpoly.ipynb`) and injects a relativistic jet. Coordinate system is cylindrical with SR. Key input parameters are in `inputs/mhd/athinput.jet_blast` (or equivalent).

## Analysis Notebooks

- **`getpoly.ipynb`** — Generates dimensionless polytrope profiles and scales to physical units; outputs CSV consumed by the `jet_blast` problem generator.
- **`yt.ipynb`** — Post-processing with the `yt` library: loads HDF5 outputs, defines custom fields (Lorentz factor, enthalpy, cylindrical→Cartesian velocity), produces slice plots, tracks shock fronts, and renders frames for video via `ffmpeg`.
- **`analysis.py`** - Analysis of the hst file output `/scratch/aripoll/athena_out/outputs/jet_blast.hst` which watches all of the conserved quantities, plus some more. And csv file output `/scratch/aripoll/athena_out/outputs/shock_breakout.csv` which tracks when the surface of the star is perturbed to get the shock breakout structure. 

Custom fields defined in `yt.ipynb`:
- Lorentz gamma from 4-velocity components
- Enthalpy: `h = 1 + 4P/ρ`
- `u_cyl_x` / `u_cyl_y` — Cartesian four-velocity from cylindrical decomposition

## SR Velocity Convention (important)

In an SR build (`-s`), Athena++ stores the **spatial four-velocity** `u^i = γ v^i` in the
primitive velocity slots `w(IVX/IVY/IVZ)`, **not** the 3-velocity. See
`src/eos/adiabatic_hydro_sr.cpp` (`prim(IVX) = gamma*v1`, inverted as
`u0 = sqrt(1 + u1² + u2² + u3²)`). Consequences:

- `γ = sqrt(1 + |u|²)`. Never compute `γ = 1/sqrt(1 - |u|²)` from these fields.
- `v = u / sqrt(1 + |u|²)`. `|u|` is unbounded and is not a speed; `|u| = 1` is `v = 0.707c`.
- Components in curvilinear coordinates are **physical (orthonormal)**, so magnitudes are
  the plain Pythagorean sum and rotations use a plain rotation matrix, with no metric
  factors. The `1/r` factors belong to the gradient,
  `|∇f|² = (∂f/∂r)² + (1/r²)(∂f/∂φ)²`, not to physical-component vectors.
- The `ur` / `u_tang` columns in `shock_front.csv` are four-velocity components (older
  files spell them `vr` / `v_tang`).

## Regression Tests

```bash
cd tst/regression
python run_tests.py [test_name]
```

Python style: `flake8` with config in `setup.cfg` (max line length 90).
C++ style: `CPPLINT` with config in `CPPLINT.cfg`.
