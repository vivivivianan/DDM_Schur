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

### Global Schur harmonic coarse space

`--schur-global-slow-coarse` augments the volume `{1,x,y,z}` basis with
RHS-seeded global harmonic Ritz modes of the balanced, right-preconditioned
interface operator. The pilot Arnoldi process uses only the existing
matrix-free Schur apply and the already factored PARDISO subdomain solves. It
does not assemble the Schur complement, collect snapshots, or change FGMRES.
Retained modes are S-orthogonalized against the physics basis and one another,
with a setup-time overlap check.
The default selects six modes from a 100-dimensional candidate space; adjust
these with `--schur-global-slow-modes` and
`--schur-global-slow-subspace-dimension`.

The RRAM26 size sweep gives 284, 278, 275, and 275 iterations for
`m=6,10,20,30` (coarse dimensions 110, 114, 124, and 134). The target of fewer
than 100 iterations was not reached; `m=10` has the best measured total time,
and retaining 20 or 30 modes does not justify the added cost.
The ten-cube regression converged in 97 iterations at `1e-10`. Full commands,
timings, and negative tuning results are recorded in
`outputs/global_slow_mode_benchmark_report.md` and
`outputs/global_slow_mode_benchmark_summary.csv`; the detailed size-sweep
metrics, including S-orthogonality checks, are in
`outputs/harmonic_enrichment_size_benchmark_summary.csv`.

The validated matrix-free stopping tolerance is `1e-10`. At `1e-12`, the Arnoldi
estimate can reach machine precision while a freshly recomputed condensed residual
stagnates at a few `1e-12` because of cancellation in repeated local elimination.
FGMRES therefore verifies its estimated stopping condition with a fresh matrix-free
Schur residual before declaring convergence.

### Physics-informed sparse Schur proxy

The recommended Schur defaults are the validated 1-ring proxy, volume
`{1,x,y,z}` coarse basis, balanced correction, FGMRES, reusable local/proxy
PARDISO factors, and block-64 phase-33 probing. No proxy/coarse flags are
needed after selecting `--solvers schur`; the older coarse/preconditioner
variants remain available only as explicit benchmark options.

`--schur-proxy` replaces only the Level-1 preconditioner action with an
independently factorized sparse approximation of the interface Schur matrix.
The exact matrix-free Schur operator is still used for FGMRES, all residual
checks, harmonic Arnoldi, and final recovery. The existing balanced coarse
correction and reusable local PARDISO factors are unchanged.

The proxy pattern is built from the interface FEM/SIPG graph and optional
high-conductivity connectivity. Use
`--schur-proxy-use-material-connectivity`,
`--schur-proxy-high-k-threshold`, and `--schur-proxy-ring` to control it.
`--schur-proxy-diagnostics` samples exact columns and writes distance-decay and
1/2/3-ring error/fill estimates before construction. `--schur-proxy-only`
disables the coarse correction for an ablation run.

Independent local phase-33 solves can be parallelized with
`--schur-local-solve-threads N`; each local PARDISO instance remains
single-threaded by default (`--schur-local-pardiso-threads 1`) to avoid nested
OpenMP/MKL oversubscription. The same MKL limit is applied to the ring-1 proxy
PARDISO analysis, factorization, and phase-33 solve. Per-subdomain analysis,
factorization, solve, and
factor-memory statistics are written to `schur_subdomain_performance.csv`.
With the ring-1 proxy and non-energy coarse space, unused full local-block
preconditioner factors are not constructed.

`--schur-proxy-cache path` enables a persistent, fingerprint-validated cache
of the interface graph, coloring, proxy sparsity, and proxy matrix. A cache
hit skips graph construction and exact color probing. Symbolic and numerical
proxy factors are reused for identical operators within the same process;
across processes the cached matrix is validated and refactorized. Any mesh,
partition, matrix, or proxy-configuration fingerprint mismatch is treated as
a cache miss rather than silently reused.

On RRAM26, the accepted 1-ring proxy plus volume `{1,x,y,z}` coarse space
reduces 278 strict-harmonic iterations to 31 and total time from 167.74 s to
110.60 s. The 2/3-ring variants are intentionally not factorized because the
locality diagnostic predicts 3.39/14.13 GiB of CSR storage before PARDISO fill.
Implementation details, A-E results, direct-solution accuracy, commands, and
regressions are recorded in `outputs/schur_proxy_benchmark_report.md` and
`outputs/schur_proxy_benchmark_summary.csv`.

The next-stage RRAM26 tolerance, block-probing, equivalence, and fixed-matrix
multi-RHS benchmarks are reproducible with
`scripts/run_schur_stage2_benchmarks.ps1`. Multi-RHS runs construct the
interface graph, exact coloring, proxy sparsity, local factors, and proxy
symbolic/numerical factors once; every online heat-source case then performs
only RHS assembly/condensation, FGMRES, and temperature recovery.
Measured RRAM26 and 595k-DOF Chiplet results are summarized in
`outputs/schur_stage2_benchmark_report.md`.

## Local MOR + physical Schur mainline

The production MOR direction is now domain decomposition with one independently
trained interior basis per physical subdomain, coupled through the original
SIPG/FEM Schur interface. It is deliberately distinct from the preserved
global-basis methods:

- `Global-Reduced-Schur-ROM`: the Stage 2A/2B reduced global interface model;
- `Global-Block-Arnoldi-ROM`: the Stage 2C.1 global-volume transient model;
- `Local-POD-Schur-ROM`: the new steady mainline;
- `Local-Block-Arnoldi-Schur-ROM`: the future transient mainline, started only
  after the steady local-ROM milestones pass.

The steady Local-ROM mainline now includes Milestone 1 (two-cube) and
Milestone 2 (ten-cube). It reduces only each subdomain interior. Interface DOFs
remain full order, every reduced block is projected from the original assembled
matrix, and the reduced local Schur contributions are assembled in the original
global interface ordering. Ten-cube uses ten independent local bases and an
exact sparse PARDISO factorization of the 10,593-DOF full interface. The
canonical validation is reproduced with:

```powershell
.\scripts\run_local_rom_schur.ps1
```

The main CLI is `--solvers local-rom` with `--local-mor-generate`,
`--local-mor-load/save`, `--local-mor-method pod|rb`,
`--local-mor-mode pure|corrected`, per-subdomain rank controls, and
`--local-interface-mode full`. Corrected mode uses the unchanged Stage 1 exact
Schur-FGMRES and true-residual gate. Detailed formulation, current milestone
boundary, diagnostics, and reproducibility contract are in
`docs/local_mor_schur_mainline.md` and `outputs/local_rom_schur_report.md`.
The optional `--local-rom-reuse-identical-subdomains` path is guarded by mesh,
material, boundary, and interface fingerprints and retains instance-specific
maps. Local transient Block Arnoldi remains out of scope until the steady
Milestone 3 cases are addressed.

## Stage 2A: fixed-matrix Reduced Schur MOR (global benchmark)

`--solvers reduced-schur` is an independent optional path for many-query,
steady thermal analysis with a fixed matrix and parameterized heat-source RHS,

\[
 A T(p)=b_0+Bp.
\]

Every configured `heat_source` entry is one independent power channel.  Model
generation solves a reference case and every physical unit channel, forms
interface temperature-rise snapshots, and computes an orthonormal POD basis by
the snapshot Gram matrix.  The reduced operator is the exact Galerkin
projection `V^T S V`: every `S*V` column uses the existing exact matrix-free
Schur apply.  The sparse proxy is never substituted for `S` in the ROM.

The first version deliberately retains exact local PARDISO interior recovery.
`pure` mode performs only the reduced interface solve and exact recovery.
`corrected` mode uses the ROM interface field as the initial guess for the
unchanged Stage 1 exact Schur-FGMRES and true-residual gate.

Generate and benchmark a model:

```powershell
.\build\Release\SIPGHeatDDM3D.exe `
  --steady --config .\configs\ten_cube_schur.txt `
  --output-dir .\results\mor_ten_cube `
  --solvers reduced-schur --mor-generate-model `
  --mor-save-model .\results\mor_ten_cube\mor_model `
  --mor-rank 10 --mor-rank-sweep 5,10,20,40,60,80,100 `
  --mor-training-count 20 --mor-validation-count 20 --mor-test-count 20 `
  --mor-random-seed 20260721 --mor-compare-fom `
  --mor-snapshot-solver auto `
  --reduced-schur-mode pure --mor-exact-interior-recovery `
  --pcg-tolerance 1e-10 --max-pcg-iterations 500 --gmres-restart 100
```

Load the immutable model and request exact correction:

```powershell
.\build\Release\SIPGHeatDDM3D.exe `
  --steady --config .\configs\ten_cube_schur.txt `
  --output-dir .\results\mor_ten_cube_corrected `
  --solvers reduced-schur `
  --mor-load-model .\results\mor_ten_cube\mor_model `
  --reduced-schur-mode corrected `
  --mor-training-count 0 --mor-validation-count 0 --mor-test-count 0 `
  --pcg-tolerance 1e-10 --max-pcg-iterations 500 --gmres-restart 100
```

The loader rejects a model if the mesh, assembled matrix, interface ordering,
boundary/source definition, or number of source channels differs.  Stage 2A is
not applicable to varying conductivity, TIM resistance, convection
coefficients, geometry, mesh, nonlinear materials, or transient systems; those
require the affine/operator and dynamic extensions reserved for Stage 2B/2C.
With MKL available, snapshot solver `auto` selects a single reusable global
PARDISO factorization and phase-33 solves; `schur` remains available as an
explicit audit backend.  Both generate snapshots on the identical assembled
matrix and DOF ordering, while `S_r` is still built only with exact Schur
applies.  Pure ROM with direct snapshots does not build the unused sparse
proxy; corrected mode and explicit Schur snapshots retain the complete Stage 1
proxy/coarse/FGMRES path.
All power vectors and train/validation/test splits are written to CSV with a
fixed seed.  `mor_summary.json`, `mor_rank_sweep.csv`,
`mor_rank_split_summary.csv`, and `mor_case_results.csv` contain the complete
offline/online timing, residual, temperature, interface-jump, heat-flux, memory,
and break-even diagnostics.  The implementation and measured results are
documented in `outputs/reduced_schur_stage2_report.md`; the four-case driver is
`scripts/run_reduced_schur_stage2.ps1`.

Measured Stage 2A decision: RRAM26 requires rank 125 to cover its 125
independent localized source directions.  That accepted model has maximum
validation relative L2 error `1.93e-13`, averages `0.194 s` online versus
`0.604 s` for reused monolithic PARDISO, and amortizes its `108.18 s` offline
cost after about 264 queries.  Rank 100 is rejected (`1.81e-2` validation/test
relative L2 error).  The four-channel Chiplet ROM is accurate to `2.70e-12`
but slower than reused monolithic PARDISO, so direct remains recommended there.

## Stage 2A.1: response-map compression (global benchmark)

Stage 2A.1 leaves the Stage 2A interface basis and reduced Schur operator
unchanged and replaces only the online local-interior recovery.  The supported
strategies are:

- `pardiso`: the Stage 2A reference path with reusable local factors;
- `exact-response`: a float64 merged power-to-temperature response, exact over
  the configured fixed-matrix power-channel space;
- `compressed-rb`: a per-subdomain thin-SVD factorization of that response.

The offline process evaluates the reference and every physical unit power
channel using the existing Reduced Schur solver, records the direct local
power-channel support, and merges the reduced interface coordinates into
`T = T_ref + H_full p`.  The deployment-only path consequently performs only
dense response matrix-vector products and global-DOF scatter.  It does not
read a FEM configuration and does not initialize assembly, Schur/proxy/coarse,
FGMRES, local PARDISO, or monolithic PARDISO state.

Generate a fixed-matrix exact-response model:

```powershell
.\build\Release\SIPGHeatDDM3D.exe `
  --steady --config .\configs\ten_cube_schur.txt `
  --output-dir .\results\stage2a1_ten_exact `
  --solvers reduced-schur --mor-generate-model `
  --mor-save-model .\results\stage2a1_ten_exact\mor_model `
  --mor-rank 10 --mor-rank-sweep 10 `
  --mor-training-count 20 --mor-validation-count 20 --mor-test-count 20 `
  --mor-snapshot-solver direct --mor-compare-fom `
  --reduced-schur-mode pure --mor-interior-mode exact-response `
  --mor-precompute-power-response --mor-storage-precision float64
```

Run 100 online queries without a mesh/configuration file:

```powershell
.\build\Release\SIPGHeatDDM3D.exe `
  --mor-deployment-only `
  --mor-load-model .\results\stage2a1_ten_exact\mor_model `
  --mor-deployment-rhs-count 100 --mor-report-io-time `
  --output-dir .\results\stage2a1_ten_deployment
```

Convert the stored exact map to compressed local bases and audit fixed ranks:

```powershell
.\build\Release\SIPGHeatDDM3D.exe `
  --mor-deployment-only `
  --mor-load-model .\results\stage2a1_rram26_exact\mor_model `
  --mor-save-model .\results\stage2a1_rram26_compressed\mor_model `
  --mor-interior-mode compressed-rb --mor-interior-rank 125 `
  --mor-interior-rank-sweep 5,10,20,40,60,80,100,125 `
  --mor-deployment-rhs-count 100 --mor-compare-interior-modes `
  --output-dir .\results\stage2a1_rram26_compression_sweep
```

For automatic local ranks, omit `--mor-interior-rank` and select discarded
energy with `--mor-interior-energy-tolerance`; the relative singular-value
threshold is controlled by `--mor-interior-singular-value-tolerance`.  Energy
is a selection heuristic only: acceptance must use the reported worst-case
full-field, absolute-temperature, hotspot, and unit-channel errors.  The
default storage is `float64`; `float32` is experimental.

All deployment models are immutable fixed-matrix artifacts.  A mesh,
conductivity, material, boundary condition, reference state, interface basis,
or power-channel-definition change invalidates the model.  Corrected mode
continues to use the unchanged Stage 1 exact correction and defaults to exact
PARDISO recovery when a true global-residual guarantee is required.  Complete
four-case reproduction is available in
`scripts/run_local_interior_rom_stage2a1.ps1`; measured results and CSV tables
are in `outputs/local_interior_rom_stage2a1_report.md` and the adjacent
`local_interior_*.csv` files.

Despite its historical `local-interior` label, Stage 2A.1 is a compressed
global deployment response map. It is not the new Local-ROM-Schur mainline and
does not replace independently trained subdomain bases coupled on a full
physical interface.

## Stage 2B.1: one matrix parameter and arbitrary power channels (global benchmark)

Stage 2B.1 generalizes the fixed-matrix deployment model to

\[
 A(\mu)T(\mu,p)=b(\mu,p),
\]

with exactly one active physical matrix parameter and every configured heat
source retained as an independent power channel.  The supported first-version
parameters are a convection coefficient selected by subdomain/boundary tag and
an isotropic material conductivity selected by subdomain/domain entity.  No
coordinate-box selection or finite-difference operator approximation is used.
Convection uses exact `A0 + h Ah` FEM boundary components.  A material region
that touches an SIPG interface additionally uses exact analytic harmonic-mean
penalty coefficient groups, assembled with the unchanged triangle-overlap
quadrature.

Offline generation factors the FOM once per training parameter and applies
PARDISO phase 33 to the reference plus every unit power channel in one
multi-RHS call.  Two-pass orthonormalized POD bases are built separately on the
global interface and each subdomain interior.  All affine local blocks and the
global interface block are then Galerkin projected.  The pure online path loads
only the self-contained dense model, evaluates affine coefficients, factors
small dense local blocks, solves the reduced Schur system, and reconstructs the
temperature in global DOF order.  It never loads a config or initializes a
sparse matrix, FEM/SIPG assembly, exact Schur, proxy, FGMRES, or PARDISO.

Generate a five-point convection model and validate unseen parameters:

```powershell
.\build\Release\SIPGHeatDDM3D.exe `
  --steady --config .\configs\ten_cube_parametric_h.txt `
  --solvers reduced-schur --mor-parametric-generate `
  --mor-parametric-save .\results\stage2b1_ten_cube\model `
  --mor-matrix-parameter convection-h `
  --mor-parameter-subdomain 9 --mor-parameter-region-id 5 `
  --mor-parameter-min 50 --mor-parameter-max 200 `
  --mor-parameter-reference 100 --mor-parameter-training-count 5 `
  --mor-parameter-validation-count 2 --mor-parameter-test-count 2 `
  --mor-interface-rank 32 --mor-local-rank 32 `
  --mor-parametric-mode pure --output-dir .\results\stage2b1_ten_cube
```

Deploy it without a mesh or config:

```powershell
.\build\Release\SIPGHeatDDM3D.exe `
  --mor-parametric-deployment-only `
  --mor-parametric-load .\results\stage2b1_ten_cube\model `
  --mor-parameter-value 125 `
  --mor-power-values 1,1,1,1,1,1,1,1,1,1 `
  --output-dir .\results\stage2b1_ten_cube_deployment
```

For a deterministic pure-deployment batch over random in-range parameter and
power-channel combinations, omit `--mor-parameter-value` and
`--mor-power-values`, then add
`--mor-deployment-rhs-count 100 --mor-seed 20260721`. The model is loaded once,
the timing CSV contains one row per dense ROM solve, and only the first case's
temperature field is written.

After choosing a rank from validation, `--mor-parametric-save <new-directory>`
may be combined with deployment-only loading and rank options to serialize the
truncated production model. Subsequent deployment should load that rank-specific
model directly so temporary rank truncation is not included in online memory.

Use `--mor-parametric-mode corrected` in a full configured run to initialize
the unchanged Stage 1 exact Schur-FGMRES correction and residual gate.  Pure
deployment rejects values outside the training interval unless
`--mor-allow-extrapolation` is explicit; extrapolated runs are flagged.  Model
loading rejects changed mesh/system/interface/source fingerprints, affine
component hashes, parameter definition, or source-channel count.

The reproducibility driver is
`scripts/run_parametric_reduced_schur_stage2b1.ps1`.  The Release tests cover
affine exactness, fixed Stage 2A.1 regression, parametric accuracy,
serialization/reload, range rejection, fingerprint rejection, and corrected
residual gating:

```powershell
ctest --test-dir build -C Release --output-on-failure
```

## Numerical safety gates

Heterogeneous SIPG interfaces use conductivity/mesh-weighted consistency
fluxes and the same `penalty_scaling` in both `max` and `harmonic` penalty
modes.  The current Schur, PARDISO-SPD, balanced coarse-correction, and local
dynamic Schur implementations require a symmetric system, so
`interface_scheme=nipg` is rejected before assembly; use
`interface_scheme=sipg`.

Local dynamic Schur time stepping applies a full discrete residual gate after
every ROM step.  Set the threshold with
`--mor-full-residual-tolerance` (default `1e-5`).  A failed ROM step is replaced
using the exact full operator for that time step and is accepted only if the
recomputed residual satisfies the gate.  Validation mode retains the eager FOM
reference solve.  Add `--mor-transient-production` to skip that reference and
create the FOM factor lazily only when fallback is required.  Production
`max-temperature` output also omits the full-field and FOM/ROM detailed-flux
CSVs.  The per-step CSV records the residual before the gate, whether fallback
was used, and the residual after fallback.

Full-interface matrix-free Dynamic Schur uses the previous time-step interface
solution as its default FGMRES initial guess.  Select the controlled benchmark
variants with `--mor-interface-initial-guess zero|previous|extrapolated`; no
choice changes the operator, tolerance, or true-residual acceptance gate.
The default interface Krylov method remains FGMRES.  The optional
`--mor-interface-krylov pcg` path is enabled only for the symmetric reduced
Schur/balanced-preconditioner configuration; it checks operator and
preconditioner curvature, verifies the exact true residual, and explicitly
falls back to FGMRES on any PCG breakdown or iteration failure.

Two exact reduced-system direct backends are available for single-node
production studies. `--mor-interface-krylov port-core` explicitly eliminates
independent physical-port leaves and can reuse an operator artifact through
`--mor-port-core-cache`. `--mor-interface-krylov augmented-direct` instead
factors the sparse interface-plus-local-ROM block system directly, avoiding
explicit `D_p^-1 E_p` construction. The latter is the fastest measured cold
path for Package15, but remains opt-in because fill and memory must be checked
on each new decomposition.

For full-interface Local Dynamic Schur, `--mor-transient-save <directory>`
writes both `local_dynamic_interior_model.bin` and the small
`local_dynamic_reference.bin` sidecar.  A later `--mor-transient-load` accepts
them only after mesh/operator/boundary fingerprints and partition ordering
match.  The reference sidecar avoids rebuilding the global steady reference
factor; an older model cache without the sidecar remains compatible and is
upgraded after one validated reference solve.

Cold local-model generation supports two construction-trace modes. The
production default, `--mor-construction-traces global-fom`, factors the steady
conductivity operator once and solves all source channels as two multi-RHS
blocks. A compatible uniform Dirichlet/convection equilibrium is now detected
analytically, so it does not require a separate reference solve. Independent
local Block-Arnoldi prototypes and their projections are built concurrently;
`--schur-local-pardiso-threads` remains the per-subdomain MKL limit, while an
`augmented-direct` interface factor uses the full
`--schur-local-solve-threads` budget.

`--mor-construction-traces operator-coarse` is an operator-only alternative.
It forms an exact matrix-free Schur action from parallel local `K_II` factors,
uses a geometric coarse correction, and converges source particular traces by
PCG without a global FOM factor or training data. It is retained for memory or
distributed-memory experiments; Package15 timing currently makes
`global-fom` the faster cold-start choice. Select the trace rank with
`--mor-interface-rank` and verify provenance/timing through the
`construction_*` columns in `local_dynamic_schur_summary.csv`.

The local-model cache header is version 4. Legacy global-trace caches remain
readable, including version 3 `global-fom` models. Version 3
`operator-coarse` models are deliberately rejected because the source-trace
semantics changed; rebuild those experimental models once with the current
executable.

For the measured fastest single-node Package15 configuration, use the fixed
production preset instead of assembling the CLI manually:

```powershell
.\scripts\run_package15_fastest_dynamic_schur.ps1 `
  -Profile large -Stage Cold -Force
.\scripts\run_package15_fastest_dynamic_schur.ps1 `
  -Profile large -Stage Warm -WarmSteps 100
```

The preset selects `global-fom` construction traces, one local Arnoldi moment,
the full physical interface, and `augmented-direct`. It disables FOM comparison
and residual replacement by default. See
[`docs/package15_fastest_dynamic_schur.md`](docs/package15_fastest_dynamic_schur.md)
for the algorithm, code map, cache contract, timing definition, and flowcharts.

Use the validated cached production wrapper for repeated transient runs.  It
refuses a nonempty output directory and reports only the principal gate,
iteration, timing, cache and memory fields:

```powershell
.\scripts\run_cached_local_dynamic_schur.ps1 `
  -Config D:\models\case.txt `
  -ModelCache E:\cache\local_dynamic_model `
  -ProxyCache E:\cache\case.proxycache `
  -OutputDirectory E:\runs\case_100step `
  -Steps 100 -Dt 1 -OuterThreads 5 -MklThreads 2
```

Both production wrappers cap the independent general worker pool at eight,
the measured Package15 optimum. `OuterThreads` is still passed in full to a
single `augmented-direct` factor, while `MklThreads` limits concurrent local
factors.

The thread defaults above are validated for RRAM5; other models still require
their own thread sweep.  All setup and every residual-gated time step remain
included in `total_seconds`.

This full residual gate is a discrete-equation consistency check;
without an independently established stability constant, it is not a strict
upper bound on temperature error.
