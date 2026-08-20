param(
    [string]$Exe = '.\build\Release\SIPGHeatDDM3D.exe',
    [string]$Root = 'results\milestone8_flux_aware_small',
    [int]$Threads = 8
)
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
Push-Location $project
try {
    foreach ($type in @('numerical','physical','both')) {
        & $Exe --transient --config configs\ten_cube_parametric_h.txt `
            --mor-transient-generate `
            --mor-transient-method local-port-block-arnoldi `
            --port-basis-method hybrid-randomized --mor-arnoldi-moments 2 `
            --randomized-port-rank 16 --randomized-oversampling 5 `
            --randomized-power-iterations 1 --randomized-seed 12345 `
            --residual-krylov-max-rank 4 --residual-krylov-max-sweeps 2 `
            --residual-krylov-tol 1e-4 --residual-krylov-block-size 4 `
            --residual-krylov-probe-mode operator-geometry `
            --optimal-port-source-mode trace-only `
            --optimal-port-inner-solver woodbury-exact `
            --history-compression-method deterministic-rrqr `
            --history-compression-rank 64 `
            --flux-aware-flux-type $type --mor-transient-dt 0.1 `
            --mor-transient-t-end 0.1 --mor-transient-waveform single_step `
            --mor-transient-initial-mode ambient `
            --output-dir (Join-Path $Root $type) --fast-run
        if ($LASTEXITCODE -ne 0) { throw "Small case failed: $type" }
    }
} finally { Pop-Location }
