param(
    [Parameter(Mandatory = $true)][string]$Config,
    [Parameter(Mandatory = $true)][string]$ModelCache,
    [Parameter(Mandatory = $true)][string]$ProxyCache,
    [Parameter(Mandatory = $true)][string]$OutputDirectory,
    [string]$Exe = '',
    [ValidateSet('local-block-arnoldi', 'local-port-block-arnoldi')]
    [string]$Method = 'local-block-arnoldi',
    [ValidateRange(1, 1000000)][int]$Steps = 100,
    [ValidateRange(1.0e-12, 1.0e12)][double]$Dt = 1.0,
    [ValidateRange(1, 1024)][int]$OuterThreads = 5,
    [ValidateRange(1, 1024)][int]$MklThreads = 2,
    [ValidateSet('fgmres', 'pcg', 'port-core', 'augmented-direct')]
    [string]$Krylov = 'fgmres',
    [ValidateRange(1, 1000000)][int]$Restart = 100,
    [ValidateRange(1, 1000000)][int]$MaxIterations = 500,
    [ValidateRange(0.0, 1.0)][double]$AdaptiveTolerance = 1.0e-9,
    [ValidateRange(0, 512)][int]$OperatorCoarseRank = 0,
    [ValidateRange(1, 8)][int]$OperatorCoarseSweeps = 2,
    [string]$OperatorCoarseCache = '',
    [switch]$DisableOperatorCoarsePredictor,
    [ValidateRange(1, 1000000)][int]$ArnoldiMoments = 2,
    [ValidateRange(1.0e-14, 1.0e-2)][double]$RankTolerance = 1.0e-10,
    [ValidateSet('global-fom', 'operator-coarse')]
    [string]$ConstructionTraces = 'global-fom',
    [ValidateRange(0, 1000)][int]$ConstructionTraceRank = 0,
    [ValidateRange(1.0e-12, 1.0)][double]$SecondMomentEnergy = 1.0,
    [ValidateRange(0, 1000000)][int]$SecondMomentMaxColumns = 0,
    [ValidateRange(1, 1000000)][int]$PortRank = 128,
    [ValidateRange(0.0, 1.0)][double]$PortEnergyTolerance = 1.0e-20,
    [ValidateRange(0, 1000000)][int]$PortEnrichmentRounds = 0,
    [ValidateRange(0.0, 1.0e12)][double]$PortTemperatureWeight = 1.0,
    [ValidateRange(0.0, 1.0e12)][double]$PortFluxWeight = 1.0,
    [ValidateRange(0.0, 1.0e12)][double]$PortResidualWeight = 1.0,
    [ValidateSet('single_step', 'multi_step', 'random', 'asynchronous_hotspots',
        'mixed_frequency')]
    [string]$Waveform = 'multi_step',
    [switch]$CompareFomSummaryOnly,
    [switch]$DisableNativeReducedHistory,
    [switch]$ReuseIdenticalSubdomains,
    [switch]$AllowCacheWarmup
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
if ([string]::IsNullOrWhiteSpace($Exe)) {
    $Exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
}
$exeFull = [IO.Path]::GetFullPath($Exe)
$configFull = [IO.Path]::GetFullPath($Config)
$modelCacheFull = [IO.Path]::GetFullPath($ModelCache)
$proxyCacheFull = [IO.Path]::GetFullPath($ProxyCache)
$operatorCoarseCacheFull = if ([string]::IsNullOrWhiteSpace($OperatorCoarseCache)) {
    "$proxyCacheFull.operator_coarse_r${OperatorCoarseRank}_s${OperatorCoarseSweeps}.bin"
} else {
    [IO.Path]::GetFullPath($OperatorCoarseCache)
}
$outputFull = [IO.Path]::GetFullPath($OutputDirectory)
$isPortMethod = $Method -eq 'local-port-block-arnoldi'

if (-not (Test-Path -LiteralPath $exeFull -PathType Leaf)) {
    throw "Release executable not found: $exeFull"
}
if (-not (Test-Path -LiteralPath $configFull -PathType Leaf)) {
    throw "Configuration file not found: $configFull"
}
$modelFile = if ([IO.Path]::GetExtension($modelCacheFull)) {
    $modelCacheFull
} else {
    Join-Path $modelCacheFull 'local_dynamic_interior_model.bin'
}
$modelDirectory = if ([IO.Path]::GetExtension($modelCacheFull)) {
    [IO.Path]::GetDirectoryName($modelCacheFull)
} else {
    $modelCacheFull
}
$referenceFile = if ([IO.Path]::GetExtension($modelCacheFull)) {
    Join-Path $modelDirectory (
        [IO.Path]::GetFileNameWithoutExtension($modelCacheFull) + '_reference.bin')
} else {
    Join-Path $modelCacheFull 'local_dynamic_reference.bin'
}
$descriptorFile = if ([IO.Path]::GetExtension($modelCacheFull)) {
    Join-Path $modelDirectory (
        [IO.Path]::GetFileNameWithoutExtension($modelCacheFull) +
        '_thermal_descriptor.bin')
} else {
    Join-Path $modelCacheFull 'thermal_descriptor.bin'
}
$portBasisFile = if ([IO.Path]::GetExtension($modelCacheFull)) {
    Join-Path $modelDirectory (
        [IO.Path]::GetFileNameWithoutExtension($modelCacheFull) + '_local_port_basis.bin')
} else {
    Join-Path $modelCacheFull 'local_port_basis.bin'
}
if (-not $isPortMethod -and -not (Test-Path -LiteralPath $modelFile -PathType Leaf)) {
    throw "Local Dynamic model cache not found: $modelFile"
}
$requiredSidecarFiles = if ($isPortMethod) {
    @($portBasisFile)
} else {
    @($referenceFile, $descriptorFile)
}
$missingCacheFiles = @()
foreach ($requiredCacheFile in $requiredSidecarFiles) {
    if (-not (Test-Path -LiteralPath $requiredCacheFile -PathType Leaf)) {
        $missingCacheFiles += $requiredCacheFile
    }
}
if ($missingCacheFiles.Count -gt 0) {
    $message = "Required cached model is missing for cached production." +
        " Missing: " + ($missingCacheFiles -join '; ')
    if (-not $AllowCacheWarmup) {
        throw $message + " Rebuild/package the cache first or pass -AllowCacheWarmup."
    }
    Write-Warning ($message + " Proceeding because -AllowCacheWarmup was set.")
}
if (Test-Path -LiteralPath $outputFull) {
    $existing = Get-ChildItem -LiteralPath $outputFull -Force -ErrorAction Stop |
        Select-Object -First 1
    if ($null -ne $existing) {
        throw "Output directory is not empty; choose a new directory: $outputFull"
    }
}

$invariant = [Globalization.CultureInfo]::InvariantCulture
$dtText = $Dt.ToString('R', $invariant)
$endText = ($Steps * $Dt).ToString('R', $invariant)
$arguments = @(
    '--transient', '--config', $configFull,
    '--mor-transient-generate',
    '--mor-transient-method', $Method,
    '--mor-transient-load', $modelCacheFull,
    '--mor-arnoldi-moments', $ArnoldiMoments,
    '--mor-arnoldi-rank-tolerance', $RankTolerance.ToString('R', $invariant),
    '--mor-construction-traces', $ConstructionTraces,
    '--mor-interface-rank', $ConstructionTraceRank,
    '--mor-arnoldi-second-moment-energy', $SecondMomentEnergy.ToString('R', $invariant),
    '--mor-arnoldi-second-moment-max-columns', $SecondMomentMaxColumns,
    '--mor-transient-dt', $dtText,
    '--mor-transient-t-end', $endText,
    '--mor-transient-waveform', $Waveform,
    '--mor-transient-initial-mode', 'ambient',
    '--mor-transient-production',
    '--mor-interface-initial-guess', 'previous',
    '--mor-interface-krylov', $Krylov,
    '--local-mor-matrix-free-threshold', 0,
    '--max-pcg-iterations', $MaxIterations,
    '--gmres-restart', $Restart,
    '--pcg-tolerance', '1e-10',
    '--mor-adaptive-interface-tolerance', $AdaptiveTolerance,
    '--schur-proxy-ring', 1,
    '--schur-proxy-block-size', 64,
    '--schur-local-solve-threads', $OuterThreads,
    '--schur-local-pardiso-threads', $MklThreads,
    '--schur-proxy-cache', $proxyCacheFull,
    '--output-dir', $outputFull,
    '--fast-run'
)
if ($Method -eq 'local-port-block-arnoldi') {
    $arguments += @(
        '--local-port-rank', $PortRank,
        '--local-port-energy-tolerance', $PortEnergyTolerance.ToString('R', $invariant),
        '--local-port-enrichment-rounds', $PortEnrichmentRounds,
        '--local-port-temperature-weight', $PortTemperatureWeight.ToString('R', $invariant),
        '--local-port-flux-weight', $PortFluxWeight.ToString('R', $invariant),
        '--local-port-residual-weight', $PortResidualWeight.ToString('R', $invariant)
    )
}
if ($CompareFomSummaryOnly) {
    # This must follow production mode, which disables FOM comparison by default.
    $arguments += '--mor-transient-compare-fom-summary-only'
}
if ($DisableNativeReducedHistory) {
    $arguments += '--no-mor-native-reduced-history'
}
if ($ReuseIdenticalSubdomains) {
    $arguments += '--mor-local-transient-reuse-identical-subdomains'
}
if ($OperatorCoarseRank -gt 0) {
    $arguments += @(
        '--schur-interface-operator-coarse-rank', $OperatorCoarseRank,
        '--schur-interface-operator-coarse-sweeps', $OperatorCoarseSweeps,
        '--schur-interface-operator-coarse-cache', $operatorCoarseCacheFull
    )
    if ($DisableOperatorCoarsePredictor) {
        $arguments += '--no-schur-interface-operator-coarse-predictor'
    }
}

$oldOmp = $env:OMP_NUM_THREADS
$oldMkl = $env:MKL_NUM_THREADS
$oldDynamic = $env:MKL_DYNAMIC
$oldWorkers = $env:SIPG_SOLVER_WORKERS
try {
    $env:OMP_NUM_THREADS = [string]$OuterThreads
    $env:MKL_NUM_THREADS = [string]$MklThreads
    $env:MKL_DYNAMIC = 'FALSE'
    $env:SIPG_SOLVER_WORKERS = [string][Math]::Min($OuterThreads, 8)
    Push-Location $project
    try {
        & $exeFull @arguments
        if ($LASTEXITCODE -ne 0) {
            throw "Cached Local Dynamic Schur failed with exit code $LASTEXITCODE."
        }
    } finally {
        Pop-Location
    }
} finally {
    $env:OMP_NUM_THREADS = $oldOmp
    $env:MKL_NUM_THREADS = $oldMkl
    $env:MKL_DYNAMIC = $oldDynamic
    $env:SIPG_SOLVER_WORKERS = $oldWorkers
}

$summaryPath = Join-Path $outputFull 'local_dynamic_schur_summary.csv'
if (-not (Test-Path -LiteralPath $summaryPath -PathType Leaf)) {
    throw "Run succeeded but summary is missing: $summaryPath"
}
$summary = Import-Csv -LiteralPath $summaryPath
[pscustomobject]@{
    status = $summary.status
    steps = [int]$summary.steps
    total_local_rank = [int]$summary.total_local_rank
    port_reduction = [int]$summary.port_reduction
    port_dimension = [int]$summary.port_dimension
    port_model_cache_hit = [int]$summary.port_model_cache_hit
    arnoldi_second_moment_energy =
        [double]$summary.arnoldi_second_moment_energy
    arnoldi_second_moment_max_columns =
        [int]$summary.arnoldi_second_moment_max_columns
    native_reduced_history_enabled =
        [int]$summary.native_reduced_history_enabled
    native_reduced_history_steps =
        [int]$summary.native_reduced_history_steps
    native_reduced_rhs_seconds =
        [double]$summary.native_reduced_rhs_seconds
    interface_iterations_total = [int]$summary.interface_iterations_total
    interface_iterations_maximum = [int]$summary.interface_iterations_maximum
    maximum_interface_relative_residual =
        [double]$summary.maximum_interface_relative_residual
    maximum_full_residual = [double]$summary.maximum_full_residual
    residual_fallback_steps = [int]$summary.residual_fallback_steps
    adaptive_interface_tolerance = [double]$summary.adaptive_interface_tolerance
    adaptive_interface_retry_steps = [int]$summary.adaptive_interface_retry_steps
    adaptive_interface_retry_iterations = [int]$summary.adaptive_interface_retry_iterations
    steady_state_detected_step = [int]$summary.steady_state_detected_step
    steady_state_reused_steps = [int]$summary.steady_state_reused_steps
    reference_cache_hit = [int]$summary.reference_cache_hit
    descriptor_cache_hit = [int]$summary.descriptor_cache_hit
    descriptor_cache_load_seconds = [double]$summary.descriptor_cache_load_seconds
    descriptor_assembly_seconds = [double]$summary.descriptor_assembly_seconds
    local_model_cache_hit = [int]$summary.local_model_cache_hit
    construction_trace_mode = [string]$summary.construction_trace_mode
    construction_trace_setup_seconds =
        [double]$summary.construction_trace_setup_seconds
    construction_global_factor_seconds =
        [double]$summary.construction_global_factor_seconds
    construction_global_solve_seconds =
        [double]$summary.construction_global_solve_seconds
    global_construction_factor_used =
        [int]$summary.global_construction_factor_used
    interface_factor_threads = [int]$summary.interface_factor_threads
    interface_solve_seconds = [double]$summary.interface_solve_seconds
    interface_operator_seconds = [double]$summary.interface_operator_seconds
    interface_preconditioner_seconds = [double]$summary.interface_preconditioner_seconds
    interface_orthogonalization_seconds = [double]$summary.interface_orthogonalization_seconds
    full_step_assembly_seconds = [double]$summary.full_step_assembly_seconds
    online_core_seconds = [double]$summary.local_online_core_seconds
    total_seconds = [double]$summary.total_seconds
    peak_working_set_bytes = [double]$summary.peak_working_set_bytes
} | Format-List
