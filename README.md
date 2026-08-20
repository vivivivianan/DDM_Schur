# SIPG global heat solver snapshot

This directory is a self-contained C++17 comparison solver for three-dimensional
steady/transient heat conduction on COMSOL-exported `.mphtxt` subdomain meshes.
It assembles a global P2 tetrahedral FEM system with SIPG/NIPG interface terms and can
solve it with MKL PARDISO. It is separate from the MFEM Schwarz executable.

## Dependencies

- CMake 3.15 or newer.
- Visual Studio 2019 C++ (validated) or another C++17 compiler (not verified).
- Intel oneAPI MKL for the validated PARDISO path.
- OpenMP is optional but recommended.

Set `MKLROOT` and `CMPLR_ROOT` before configuring. The build deliberately does not
contain machine-specific oneAPI installation paths.

## Build

From the repository root:

```powershell
cmake -S external_sipg_snapshot -B external_sipg_snapshot/build -G "Visual Studio 16 2019" -A x64
cmake --build external_sipg_snapshot/build --config Release -j 4
```

## Run

```powershell
$env:MFEM_HEAT_MODEL_ROOT = "X:\shared\heat-models"
.\external_sipg_snapshot\build\Release\SIPGHeatDDM3D.exe `
  --steady `
  --config cpp_mfem_heat\configs\eight_layer_comsol_sipg_global_pardiso.txt `
  --output-dir results\sipg_eight_layer `
  --solvers direct --direct-mode general
```

Config path fields support `${VARIABLE_NAME}` expansion. Large meshes and generated
results remain outside Git; see the repository-level `README.md` for the data policy.

Use `SIPGHeatDDM3D.exe --help` for all solver, diagnostics and transient options.

## Matrix-free two-level Schur solve

The `schur` solver keeps the condensed interface operator matrix-free and applies a
balanced two-level preconditioner. Its default local coarse space is `{1,x,y,z}`:
one constant plus orthonormalized linear coordinate modes per subdomain. Disable
the x/y enrichment with `--no-schur-linear-xy-coarse` or the z mode with
`--no-schur-linear-z-coarse`. For uniform volumetric heating along a long z-chain,
`--schur-global-quadratic-z-coarse` adds one global curvature mode without adding a
quadratic mode per subdomain.

```powershell
.\build\Release\SIPGHeatDDM3D.exe `
  --steady `
  --config .\configs\ten_cube_schur.txt `
  --output-dir .\results\ten_cube_two_level_verified `
  --solvers schur,direct `
  --schur-linear-xy-coarse `
  --schur-global-quadratic-z-coarse `
  --direct-mode spd `
  --pcg-tolerance 1e-10 `
  --max-pcg-iterations 500 `
  --gmres-restart 400
```

Each immutable local PARDISO matrix is prepared with phase 11 (symbolic analysis)
once and phase 22 (numerical factorization) once. Every condensed Schur matvec,
preconditioner application, condensed right-hand side, and recovery then reuses the
stored factors through phase 33 only. `schur_solver_summary.csv` records the FGMRES
iteration count, coarse dimension, phase-11 and phase-22 call counts and times,
repeated local solve time, coarse solve time, and total FGMRES time.

### Deterministic interface-patch coarse space

`--schur-interface-patch-coarse` replaces the volume/subdomain polynomial basis
with a basis built directly on algebraic Schur-interface DOFs. A patch is keyed
by `(subdomain, neighbor, local boundary entity, neighbor boundary entity)`, so
the two sides of a physical interface remain separate and can represent both
mean and jump errors. Every patch contributes one normalized constant mode.
Shared stencil DOFs are assigned to the first sorted patch, which makes the
supports disjoint and avoids dependent constant modes.

SIPG normal-flux terms couple every P2 DOF in a boundary-tetrahedron stencil,
including DOFs whose shape functions vanish on the geometric face. Patch
support therefore uses the complete boundary-tet stencil and then retains only
DOFs classified as algebraic Schur-interface DOFs. Add centered, locally
orthonormalized x/y modes with `--schur-interface-patch-linear-xy`; this
enrichment is disabled by default.

```powershell
.\build\Release\SIPGHeatDDM3D.exe `
  --steady `
  --config D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt `
  --output-dir .\results\rram26_interface_boundary_stencil_patch_constant_tol1e5 `
  --solvers schur `
  --schur-interface-patch-coarse `
  --no-schur-interface-patch-linear-xy `
  --no-schur-global-quadratic-z-coarse `
  --pcg-tolerance 1e-5 `
  --max-pcg-iterations 500 `
  --gmres-restart 100 `
  --fast-run
```

The RRAM26 deterministic benchmark reduced 381 iterations to 297 with 2500
patch constants. Patch-local x/y enrichment used 6997 modes and reduced the
count to 276, but its 386.45 s setup made total time worse. The enrichment is
therefore intentionally optional. Full results are in
`outputs/interface_patch_benchmark_report.md` and
`outputs/interface_patch_benchmark_summary.csv`.

The validated matrix-free stopping tolerance is `1e-10`. At `1e-12`, the Arnoldi
estimate can reach machine precision while a freshly recomputed condensed residual
stagnates at a few `1e-12` because of cancellation in repeated local elimination.
FGMRES therefore verifies its estimated stopping condition with a fresh matrix-free
Schur residual before declaring convergence.
