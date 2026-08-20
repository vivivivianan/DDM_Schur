param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('rram26','chiplet')]
    [string]$Case,
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [string]$ResultsDirectory = '',
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$configPath = (Resolve-Path -LiteralPath $Config).Path
$expectedInterfaces = if ($Case -eq 'rram26') { 25 } else { 1 }
$sourceMode = if ($Case -eq 'rram26') {
    'trace-only'
} else {
    'generalized-dynamic'
}
$dt = if ($Case -eq 'rram26') { 0.01 } else { 0.1 }
$moments = if ($Case -eq 'rram26') { 1 } else { 4 }
if ([string]::IsNullOrWhiteSpace($ResultsDirectory)) {
    $ResultsDirectory = "results\milestone8_production_basis_$Case"
}
$result = [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite production basis directory: $result"
}
if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
    throw "Release executable not found: $exe"
}

$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
New-Item -ItemType Directory -Path $result | Out-Null

Push-Location $project
try {
    & $exe --transient --config $configPath `
        --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --port-basis-method hybrid-randomized `
        --mor-arnoldi-moments $moments `
        --randomized-port-rank 16 `
        --randomized-oversampling 5 `
        --randomized-power-iterations 1 `
        --randomized-seed 12345 `
        --residual-krylov-max-rank 4 `
        --residual-krylov-max-sweeps 2 `
        --residual-krylov-tol 1e-4 `
        --residual-krylov-block-size 4 `
        --residual-krylov-probe-mode operator-geometry `
        --residual-krylov-inner-solver woodbury-exact `
        --optimal-port-source-mode $sourceMode `
        --optimal-port-inner-solver woodbury-exact `
        --optimal-port-inner-tol 1e-10 `
        --optimal-port-inner-refinement-max-iters 3 `
        --optimal-port-inner-refinement-tol 1e-10 `
        --history-compression-method deterministic-rrqr `
        --history-compression-rank 64 `
        --history-compression-tolerance 1e-12 `
        --milestone8-production-basis-only `
        --mor-transient-dt $dt `
        --mor-transient-t-end $dt `
        --mor-transient-waveform single_step `
        --mor-transient-initial-mode ambient `
        --mor-transient-save $result `
        --output-dir $result --fast-run
    if ($LASTEXITCODE -ne 0) {
        throw "M8 production basis build failed for $Case."
    }

    $summaryPath = Join-Path $result `
        'milestone8_production_basis_summary.csv'
    $summary = Import-Csv -LiteralPath $summaryPath |
        Select-Object -First 1
    if ($summary.status -ne 'success' -or
        [int]$summary.physical_interfaces -ne $expectedInterfaces -or
        [int]$summary.history_rank -ne 64 -or
        [int]$summary.randomized_rank -ne 16 -or
        [int]$summary.residual_rank -ne 4 -or
        [int]$summary.snapshot_used -ne 0 -or
        [int]$summary.fom_used_for_basis -ne 0 -or
        [int]$summary.pod_used -ne 0 -or
        [int]$summary.svd_used -ne 0 -or
        [int]$summary.full_field_read -ne 0 -or
        [int]$summary.transient_advanced -ne 0) {
        throw "M8 production basis provenance/gate failed for $Case."
    }
    if (-not (Test-Path -LiteralPath (
            Join-Path $result 'local_port_basis.bin') -PathType Leaf)) {
        throw "Serialized M8 basis is missing for $Case."
    }
    Copy-Item -LiteralPath $summaryPath -Destination (
        Join-Path $project "outputs\milestone8_${Case}_production_basis.csv")
} finally {
    Pop-Location
}
