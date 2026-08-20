param(
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [Parameter(Mandatory = $true)]
    [string]$AdaptiveModel,
    [string]$FixedAccuracyDirectory =
        'results\milestone8_accuracy_rram26',
    [string]$ResultsDirectory =
        'results\milestone8_adaptive_port_accuracy_rram26',
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$configPath = (Resolve-Path -LiteralPath $Config).Path
$modelPath = (Resolve-Path -LiteralPath $AdaptiveModel).Path
$fixedRoot = if ([IO.Path]::IsPathRooted($FixedAccuracyDirectory)) {
    [IO.Path]::GetFullPath($FixedAccuracyDirectory)
} else {
    [IO.Path]::GetFullPath(
        (Join-Path $project $FixedAccuracyDirectory))
}
$result = if ([IO.Path]::IsPathRooted($ResultsDirectory)) {
    [IO.Path]::GetFullPath($ResultsDirectory)
} else {
    [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite M8.9 accuracy directory: $result"
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
        --mor-arnoldi-moments 1 `
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
        --optimal-port-source-mode trace-only `
        --optimal-port-inner-solver woodbury-exact `
        --optimal-port-inner-tol 1e-10 `
        --optimal-port-inner-refinement-max-iters 3 `
        --optimal-port-inner-refinement-tol 1e-10 `
        --history-compression-method deterministic-rrqr `
        --history-compression-rank 64 `
        --history-compression-tolerance 1e-12 `
        --milestone8-adaptive-production `
        --mor-transient-load $modelPath `
        --mor-transient-dt 0.01 `
        --mor-transient-t-end 0.01 `
        --mor-transient-waveform mixed_frequency `
        --mor-transient-initial-mode ambient `
        --output-dir $result --fast-run
    if ($LASTEXITCODE -ne 0) {
        throw 'M8.9 adaptive one-step accuracy run failed.'
    }

    $adaptiveSummary = Import-Csv -LiteralPath (
        Join-Path $result 'local_dynamic_schur_summary.csv') |
        Select-Object -First 1
    if ([int]$adaptiveSummary.port_snapshot_used -ne 0 -or
        [int]$adaptiveSummary.port_fom_used_for_basis -ne 0) {
        throw 'M8.9 adaptive one-step basis provenance gate failed.'
    }
    $fixedCsv = Join-Path $fixedRoot 'milestone8_accuracy.csv'
    $fixedRows = @(Import-Csv -LiteralPath $fixedCsv)
    if ($fixedRows.Count -ne 3) {
        throw 'Expected the validated M8.7 three-method accuracy table.'
    }
    $adaptiveRow = [pscustomobject]@{
        case = 'rram26'
        method = 'M8.9 Adaptive Operator Port'
        scope = 'one_step'
        port_dimension = $adaptiveSummary.port_dimension
        total_local_rank = $adaptiveSummary.total_local_rank
        temperature_relative_l2 =
            $adaptiveSummary.space_time_relative_l2
        max_error_k = $adaptiveSummary.maximum_absolute_k
        max_temperature_error_k =
            $adaptiveSummary.maximum_temperature_error_k
        flux_relative_l2 =
            $adaptiveSummary.maximum_fom_rom_flux_relative_l2
        interface_residual =
            $adaptiveSummary.maximum_interface_relative_residual
        full_residual = $adaptiveSummary.maximum_full_residual
        online_time_s = $adaptiveSummary.local_online_core_seconds
        total_time_s = $adaptiveSummary.total_seconds
        peak_memory_bytes = $adaptiveSummary.peak_working_set_bytes
        snapshot_used = 0
        fom_used_for_basis = 0
        corrected = 0
        status = $adaptiveSummary.status
    }
    @($fixedRows) + @($adaptiveRow) |
        Export-Csv -NoTypeInformation -Encoding UTF8 (
            Join-Path $result 'milestone8_adaptive_port_accuracy.csv')
    Copy-Item -LiteralPath (
        Join-Path $result 'milestone8_adaptive_port_accuracy.csv'
    ) -Destination (
        Join-Path $project `
            'outputs\milestone8_adaptive_port_accuracy.csv')
} finally {
    Pop-Location
}
