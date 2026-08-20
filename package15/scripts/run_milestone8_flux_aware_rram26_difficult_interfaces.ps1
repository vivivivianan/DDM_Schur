param(
    [Parameter(Mandatory=$true)][string]$Config,
    [Parameter(Mandatory=$true)][string]$AdaptiveModel,
    [string]$OutputRoot = 'results\milestone8_flux_aware_rram26',
    [int]$Threads = 8
)
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
Push-Location $project
try {
    foreach ($type in @('numerical','both')) {
        $output = Join-Path $OutputRoot $type
        & $exe --transient --config $Config --mor-transient-generate `
            --mor-transient-method local-port-block-arnoldi `
            --port-basis-method hybrid-randomized --mor-arnoldi-moments 1 `
            --randomized-port-rank 16 --randomized-oversampling 5 `
            --randomized-power-iterations 1 --randomized-seed 12345 `
            --residual-krylov-max-rank 4 --residual-krylov-max-sweeps 2 `
            --residual-krylov-tol 1e-4 --residual-krylov-block-size 4 `
            --residual-krylov-probe-mode operator-geometry `
            --optimal-port-source-mode trace-only `
            --optimal-port-inner-solver woodbury-exact `
            --history-compression-method deterministic-rrqr `
            --history-compression-rank 64 `
            --flux-aware-flux-type $type `
            --milestone8-flux-operator-audit `
            --projection-interface-ids 4,10,13,18,23 `
            --mor-transient-load $AdaptiveModel --mor-transient-dt 0.01 `
            --mor-transient-t-end 0.01 `
            --mor-transient-waveform mixed_frequency `
            --mor-transient-initial-mode ambient `
            --output-dir $output --fast-run
        if ($LASTEXITCODE -ne 0) { throw "RRAM prototype failed: $type" }
    }
    throw 'Projection gate evaluation is required. Do not continue to a full 25-interface or transient run automatically.'
} finally { Pop-Location }
