param(
    [Parameter(Mandatory=$true)][string]$Config,
    [Parameter(Mandatory=$true)][string]$Model,
    [string]$InterfaceIds = '4,10,13,18,23',
    [string]$Output = 'results\milestone8_flux_operator_audit',
    [int]$Threads = 8
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $root 'build\Release\SIPGHeatDDM3D.exe'
$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
Push-Location $root
try {
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
        --milestone8-adaptive-production `
        --milestone8-flux-operator-audit `
        --projection-interface-ids $InterfaceIds `
        --mor-transient-load $Model --mor-transient-dt 0.01 `
        --mor-transient-t-end 0.01 `
        --mor-transient-waveform mixed_frequency `
        --mor-transient-initial-mode ambient --output-dir $Output --fast-run
    if ($LASTEXITCODE -ne 0) { throw 'M8.13 flux audit failed.' }
} finally { Pop-Location }
