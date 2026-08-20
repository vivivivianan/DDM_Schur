param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('rram26','chiplet')]
    [string]$Case,
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [Parameter(Mandatory = $true)]
    [string]$M8Model,
    [string]$ResultsDirectory = '',
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$configPath = (Resolve-Path -LiteralPath $Config).Path
$modelPath = (Resolve-Path -LiteralPath $M8Model).Path
if ([string]::IsNullOrWhiteSpace($ResultsDirectory)) {
    $ResultsDirectory = "results\milestone8_efficiency_$Case"
}
$result = [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite production efficiency directory: $result"
}
if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
    throw "Release executable not found: $exe"
}

$dt = if ($Case -eq 'rram26') { '0.01' } else { '0.1' }
$endTime = if ($Case -eq 'rram26') { '0.8' } else { '1.0' }
$moments = if ($Case -eq 'rram26') { '1' } else { '4' }
$waveform = if ($Case -eq 'rram26') {
    'mixed_frequency'
} else {
    'asynchronous_hotspots'
}
$sourceMode = if ($Case -eq 'rram26') {
    'trace-only'
} else {
    'generalized-dynamic'
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
        --mor-interface-rank 0 `
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
        --mor-transient-load $modelPath `
        --mor-transient-dt $dt `
        --mor-transient-t-end $endTime `
        --mor-transient-waveform $waveform `
        --mor-transient-initial-mode ambient `
        --mor-deployment-rhs-count 100 `
        --output-dir $result --fast-run
    if ($LASTEXITCODE -ne 0) {
        throw "M8 production full-horizon/100-RHS run failed: $Case"
    }

    $summary = Import-Csv -LiteralPath (
        Join-Path $result 'local_dynamic_schur_summary.csv') |
        Select-Object -First 1
    $deployment = Import-Csv -LiteralPath (
        Join-Path $result 'local_dynamic_schur_deployment_timing.csv') |
        Select-Object -First 1
    if (($summary.status -ne 'success' -and
            $summary.status -ne 'accuracy_failed') -or
        [int]$summary.port_snapshot_used -ne 0 -or
        [int]$summary.port_fom_used_for_basis -ne 0 -or
        [int]$deployment.waveforms -ne 100 -or
        [int]$deployment.setup_reused -ne 1 -or
        [int]$deployment.port_factor_reused -ne 1 -or
        [int]$deployment.local_factors_reused -ne 1) {
        throw "M8 production efficiency/provenance gate failed: $Case"
    }

    [pscustomobject]@{
        case = $Case
        method = 'M8 Hybrid Operator Port'
        scope = 'full_horizon'
        steps = $summary.steps
        port_dimension = $summary.port_dimension
        total_local_rank = $summary.total_local_rank
        temperature_relative_l2 = $summary.space_time_relative_l2
        max_error_k = $summary.maximum_absolute_k
        flux_relative_l2 = $summary.maximum_fom_rom_flux_relative_l2
        interface_residual =
            $summary.maximum_interface_relative_residual
        full_residual = $summary.maximum_full_residual
        online_time_s = $summary.local_online_core_seconds
        total_time_s = $summary.total_seconds
        rhs_count = $deployment.waveforms
        rhs_total_online_time_s = $deployment.total_online_seconds
        average_online_time_per_rhs_s =
            $deployment.average_online_seconds_per_waveform
        peak_memory_bytes = $summary.peak_working_set_bytes
        model_bytes = $summary.model_bytes
        snapshot_used = 0
        fom_used_for_basis = 0
        corrected = 0
        status = $summary.status
    } | Export-Csv -NoTypeInformation -Encoding UTF8 (
        Join-Path $result 'milestone8_efficiency.csv')
} finally {
    Pop-Location
}
