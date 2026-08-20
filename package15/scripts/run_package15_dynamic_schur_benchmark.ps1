param(
    [ValidateSet('smoke', 'medium', 'large')][string]$Profile = 'medium',
    [ValidateSet('Generate', 'Exact', 'Cold', 'Warm', 'All')]
    [string]$Stage = 'All',
    [string]$Exe = '',
    [string]$Root = '',
    [ValidateRange(1, 12)][int]$ArnoldiMoments = 1,
    [ValidateSet('global-fom', 'operator-coarse')]
    [string]$ConstructionTraces = 'global-fom',
    [ValidateRange(0, 1000)][int]$ConstructionTraceRank = 19,
    [ValidateRange(0, 1000000)][int]$SecondMomentMaxColumns = 0,
    [ValidateRange(1.0e-14, 1.0e-2)][double]$RankTolerance = 1.0e-6,
    [ValidateRange(0, 1000)][int]$EnrichmentRounds = 0,
    [ValidateRange(1, 1000)][int]$ValidationSteps = 10,
    [ValidateRange(1, 1000000)][int]$WarmSteps = 100,
    [ValidateRange(1.0e-12, 1.0e6)][double]$Dt = 0.05,
    [ValidateRange(1, 1024)][int]$OuterThreads = 8,
    [ValidateRange(1, 1024)][int]$MklThreads = 1,
    [ValidateSet('fgmres', 'pcg', 'port-core', 'augmented-direct')]
    [string]$InterfaceSolver = 'fgmres',
    [ValidateRange(1, 1000000)][int]$Restart = 100,
    [ValidateRange(1, 1000000)][int]$MaxIterations = 500,
    [switch]$DisableTemplateReuse,
    [switch]$NoFomComparison,
    [switch]$EnableResidualFallback,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$invariant = [Globalization.CultureInfo]::InvariantCulture
$trimSeparators = [char[]]@(
    [IO.Path]::DirectorySeparatorChar, [IO.Path]::AltDirectorySeparatorChar)

function Test-StrictChildPath([string]$Path, [string]$Parent) {
    $childFull = [IO.Path]::GetFullPath($Path).TrimEnd($trimSeparators)
    $parentFull = [IO.Path]::GetFullPath($Parent).TrimEnd($trimSeparators)
    $prefix = $parentFull + [IO.Path]::DirectorySeparatorChar
    return -not $childFull.Equals(
        $parentFull, [StringComparison]::OrdinalIgnoreCase) -and
        $childFull.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)
}

if ([string]::IsNullOrWhiteSpace($Exe)) {
    $Exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
}
if ([string]::IsNullOrWhiteSpace($Root)) {
    $rankToleranceKey = $RankTolerance.ToString('R', $invariant).
        Replace('.', 'p').Replace('+', 'p').Replace('-', 'm')
    $interfaceKey = if ($InterfaceSolver -eq 'fgmres') {
        ''
    } else {
        "_is$($InterfaceSolver.Replace('-', '_'))"
    }
    $fallbackKey = if ($EnableResidualFallback) { '' } else { '_nofb' }
    $constructionKey = if ($ConstructionTraces -eq 'global-fom') {
        '_ctglobal'
    } else {
        "_ctoperator_r$ConstructionTraceRank"
    }
    $Root = Join-Path $project (
        "build\package15_benchmark_${Profile}_m${ArnoldiMoments}_c${SecondMomentMaxColumns}_e${EnrichmentRounds}_rt${rankToleranceKey}${interfaceKey}${constructionKey}${fallbackKey}")
}
$exeFull = [IO.Path]::GetFullPath($Exe)
$rootFull = [IO.Path]::GetFullPath($Root)
$buildRoot = [IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not (Test-StrictChildPath $rootFull $buildRoot)) {
    throw "Benchmark output must remain under build: $rootFull"
}
if (-not (Test-Path -LiteralPath $exeFull -PathType Leaf)) {
    throw "Release executable not found: $exeFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

function Reset-Directory([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    if (-not (Test-StrictChildPath $full $rootFull)) {
        throw "Refusing to reset directory outside benchmark root: $full"
    }
    if (Test-Path -LiteralPath $full) {
        $resolved = [IO.Path]::GetFullPath((Resolve-Path -LiteralPath $full).Path)
        if (-not (Test-StrictChildPath $resolved $rootFull)) {
            throw "Resolved reset target escaped benchmark root: $resolved"
        }
        Remove-Item -LiteralPath $resolved -Recurse -Force
    }
    New-Item -ItemType Directory -Force -Path $full | Out-Null
}

function Invoke-Solver([string[]]$Arguments, [string]$LogPath) {
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
            & $exeFull @Arguments 2>&1 | Tee-Object -FilePath $LogPath
            if ($LASTEXITCODE -ne 0) {
                throw "SIPGHeatDDM3D failed with exit code $LASTEXITCODE."
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
}

function Add-FomComparison([System.Collections.ArrayList]$Arguments) {
    if (-not $NoFomComparison) {
        [void]$Arguments.Add('--mor-transient-compare-fom-summary-only')
    }
}

function Convert-MetricDouble($Value) {
    $text = ([string]$Value).Trim()
    if ($text -match '^(?i:nan)$') { return [double]::NaN }
    if ($text -match '^\+(?i:inf(?:inity)?)$') {
        return [double]::PositiveInfinity
    }
    if ($text -match '^-(?i:inf(?:inity)?)$') {
        return [double]::NegativeInfinity
    }
    return [double]::Parse(
        $text, [Globalization.NumberStyles]::Float, $invariant)
}

function Convert-OptionalMetricDouble(
    [pscustomobject]$Summary,
    [string]$Name) {
    $property = $Summary.PSObject.Properties[$Name]
    if ($null -eq $property) { return 0.0 }
    return Convert-MetricDouble $property.Value
}

function New-ReportRow([string]$Run, [pscustomobject]$Summary) {
    $online = Convert-MetricDouble $Summary.local_online_core_seconds
    $fomSolve = Convert-MetricDouble $Summary.fom_solve_seconds
    $dynamicTotal = Convert-MetricDouble $Summary.total_seconds
    $fomFactor = Convert-MetricDouble $Summary.fom_factor_seconds
    $fomAssembly = Convert-OptionalMetricDouble `
        $Summary 'full_step_assembly_seconds'
    $relativeL2 = Convert-MetricDouble $Summary.space_time_relative_l2
    $temperatureError = Convert-MetricDouble $Summary.maximum_temperature_error_k
    $maximumResidual = Convert-MetricDouble $Summary.maximum_full_residual
    $fallbackSteps = [int]$Summary.residual_fallback_steps
    $fallbackEnabled = [int](Convert-OptionalMetricDouble `
        $Summary 'residual_fallback_enabled')
    $residualViolationSteps = [int](Convert-OptionalMetricDouble `
        $Summary 'residual_tolerance_violation_steps')
    $fomComparisonEnabled = [int]$Summary.fom_comparison_enabled -eq 1
    # Summary-only FOM trajectory solves are validation overhead. The FOM
    # factor remains production cost when residual fallback actually uses it.
    $validationOnlySeconds = if ($fomComparisonEnabled) {
        $fomSolve + $(if ($fallbackSteps -eq 0) {
            $fomFactor + $fomAssembly
        } else { 0.0 })
    } else { 0.0 }
    $productionTotal = [Math]::Max(0.0, $dynamicTotal - $validationOnlySeconds)
    $onlineSpeedup = if ($online -gt 0.0 -and $fomSolve -gt 0.0) {
        $fomSolve / $online
    } else { [double]::NaN }
    $coldSpeedup = if ($productionTotal -gt 0.0 -and ($fomFactor + $fomSolve) -gt 0.0) {
        ($fomFactor + $fomSolve) / $productionTotal
    } else { [double]::NaN }
    $accuracyPassed =
        $fomComparisonEnabled -and
        $fallbackSteps -eq 0 -and
        -not [double]::IsNaN($relativeL2) -and $relativeL2 -le 1.0e-4 -and
        -not [double]::IsNaN($temperatureError) -and $temperatureError -le 0.05 -and
        $maximumResidual -le 1.0e-4
    $averageIterations = if ([int]$Summary.steps -gt 0) {
        [double]$Summary.interface_iterations_total / [int]$Summary.steps
    } else { [double]::NaN }
    return [pscustomobject]@{
        run = $Run
        profile = $Profile
        arnoldi_moments = $ArnoldiMoments
        construction_trace_mode = [string]$Summary.construction_trace_mode
        construction_trace_rank = $ConstructionTraceRank
        rank_tolerance = $RankTolerance
        enrichment_rounds = $EnrichmentRounds
        status = [string]$Summary.status
        interface_solver = [string]$Summary.interface_solver
        steps = [int]$Summary.steps
        global_dofs = [int64]$Summary.global_dofs
        full_interface_dofs = [int64]$Summary.full_interface_dofs
        total_local_rank = [int]$Summary.total_local_rank
        unique_templates = [int]$Summary.unique_templates
        reused_instances = [int]$Summary.reused_instances
        descriptor_cache_hit = [int]$Summary.descriptor_cache_hit
        reference_cache_hit = [int]$Summary.reference_cache_hit
        local_model_cache_hit = [int]$Summary.local_model_cache_hit
        port_core_cache_hit = [int](Convert-OptionalMetricDouble `
            $Summary 'port_core_cache_hit')
        port_core_cache_load_seconds = Convert-OptionalMetricDouble `
            $Summary 'port_core_cache_load_seconds'
        port_core_cache_save_seconds = Convert-OptionalMetricDouble `
            $Summary 'port_core_cache_save_seconds'
        port_core_cache_gib = (Convert-OptionalMetricDouble `
            $Summary 'port_core_cache_bytes') / 1GB
        interface_iterations_total = [int]$Summary.interface_iterations_total
        interface_iterations_maximum = [int]$Summary.interface_iterations_maximum
        interface_iterations_average = $averageIterations
        residual_fallback_steps = $fallbackSteps
        residual_fallback_enabled = $fallbackEnabled
        residual_tolerance_violation_steps = $residualViolationSteps
        residual_fallback_solve_seconds =
            Convert-MetricDouble $Summary.residual_fallback_solve_seconds
        relative_l2 = $relativeL2
        maximum_temperature_error_k = $temperatureError
        maximum_full_residual = $maximumResidual
        maximum_full_residual_before_gate =
            [double]$Summary.maximum_full_residual_before_gate
        local_basis_setup_seconds = [double]$Summary.local_basis_setup_seconds
        construction_trace_setup_seconds = Convert-OptionalMetricDouble `
            $Summary 'construction_trace_setup_seconds'
        construction_pardiso_threads = [int](Convert-OptionalMetricDouble `
            $Summary 'construction_pardiso_threads')
        construction_global_factor_seconds = Convert-OptionalMetricDouble `
            $Summary 'construction_global_factor_seconds'
        construction_global_solve_seconds = Convert-OptionalMetricDouble `
            $Summary 'construction_global_solve_seconds'
        operator_coarse_trace_dimension = [int](Convert-OptionalMetricDouble `
            $Summary 'operator_coarse_trace_dimension')
        operator_trace_krylov_iterations = [int](Convert-OptionalMetricDouble `
            $Summary 'operator_trace_krylov_iterations')
        analytic_reference_used = [int](Convert-OptionalMetricDouble `
            $Summary 'analytic_reference_used')
        global_construction_factor_used = [int](Convert-OptionalMetricDouble `
            $Summary 'global_construction_factor_used')
        interface_factor_threads = [int](Convert-OptionalMetricDouble `
            $Summary 'interface_factor_threads')
        enrichment_setup_seconds = [double]$Summary.enrichment_total_seconds
        interface_solve_seconds = [double]$Summary.interface_solve_seconds
        interface_operator_seconds = [double]$Summary.interface_operator_seconds
        interface_preconditioner_seconds =
            [double]$Summary.interface_preconditioner_seconds
        interface_orthogonalization_seconds =
            [double]$Summary.interface_orthogonalization_seconds
        port_forward_solve_seconds =
            Convert-OptionalMetricDouble $Summary 'port_forward_solve_seconds'
        port_core_solve_seconds =
            Convert-OptionalMetricDouble $Summary 'port_core_solve_seconds'
        port_back_substitution_seconds =
            Convert-OptionalMetricDouble $Summary 'port_back_substitution_seconds'
        fom_factor_seconds = $fomFactor
        fom_solve_seconds = $fomSolve
        fom_step_assembly_seconds = $fomAssembly
        dynamic_online_seconds = $online
        measured_wall_seconds = $dynamicTotal
        fom_validation_seconds_removed = $validationOnlySeconds
        dynamic_production_total_seconds = $productionTotal
        online_speedup_vs_fom = $onlineSpeedup
        end_to_end_speedup_vs_fom = $coldSpeedup
        peak_working_set_gib = [double]$Summary.peak_working_set_bytes / 1GB
        fom_comparison_enabled = [int]$fomComparisonEnabled
        accuracy_gate_passed = [int]$accuracyPassed
    }
}

$generator = Join-Path $project 'tools\generate_package15_mesh.py'
& python $generator --profile $Profile
if ($LASTEXITCODE -ne 0) { throw 'package15 mesh generation failed.' }
$config = Join-Path $project "configs\package15_$Profile.txt"
$reports = @()

if ($Stage -eq 'Generate') {
    Get-Content (Join-Path $project "data\generated\package15\$Profile\mesh_manifest.json")
    return
}

if ($Stage -eq 'Exact' -or ($Stage -eq 'All' -and $Profile -eq 'smoke')) {
    $exactOutput = Join-Path $rootFull 'exact_one_step'
    if ($Force -or -not (Test-Path -LiteralPath $exactOutput)) {
        Reset-Directory $exactOutput
        $exactConfig = Join-Path $exactOutput 'package15_one_step.txt'
        $absoluteData = (Join-Path $project 'data\generated').Replace('\', '/')
        $exactConfigLines = Get-Content -LiteralPath $config | ForEach-Object {
            if ($_ -match '^\s*time_steps\s*=') {
                'time_steps = 1'
            } else {
                $_ -replace '\.\./data/generated', $absoluteData
            }
        }
        $exactConfigLines | Set-Content -LiteralPath $exactConfig -Encoding ascii
        $exactArguments = @(
            '--transient', '--config', $exactConfig,
            '--solvers', 'direct,schur-direct-exact', '--direct-mode', 'general',
            '--schur-direct-verify-operator', 'true',
            '--schur-direct-random-checks', '1',
            '--output-dir', $exactOutput, '--fast-run'
        )
        Invoke-Solver $exactArguments (Join-Path $exactOutput 'run.log')
    }
    $exactSummary = Import-Csv (Join-Path $exactOutput 'schur_direct_exact_summary.csv')
    $exactStatus = ($exactSummary | Where-Object field -eq 'schur_direct_status').value
    $operatorError = [double](($exactSummary |
        Where-Object field -eq 'operator_relative_error').value)
    if ($exactStatus -notin @('ready', 'success') -or $operatorError -gt 1.0e-10) {
        throw "Exact Schur gate failed: status=$exactStatus, operator_error=$operatorError"
    }
    $comparison = Import-Csv (Join-Path $exactOutput 'solver_comparison.csv')
    $exactRow = $comparison | Where-Object solver -eq 'Exact-Schur-Direct'
    if ($null -eq $exactRow -or $exactRow.status -ne 'success' -or
        [double]$exactRow.relative_l2_diff_vs_global -gt 1.0e-10) {
        throw 'Exact Schur did not match the one-step monolithic solution.'
    }
    [pscustomobject]@{
        exact_status = $exactStatus
        operator_relative_error = $operatorError
        relative_l2_vs_monolithic = [double]$exactRow.relative_l2_diff_vs_global
    } | Export-Csv (Join-Path $rootFull 'exact_gate.csv') -NoTypeInformation
}

$dynamicRoot = Join-Path $rootFull 'dynamic'
$modelCache = Join-Path $dynamicRoot 'local_dynamic_model'
$proxyCache = Join-Path $dynamicRoot 'dynamic_schur.proxycache'
$portCoreCache = Join-Path $dynamicRoot 'port_core_artifacts.bin'
$coldOutput = Join-Path $dynamicRoot 'cold_validation'

if ($Stage -eq 'Cold' -or $Stage -eq 'All') {
    if ($Force -and (Test-Path -LiteralPath $dynamicRoot)) {
        Reset-Directory $dynamicRoot
    } else {
        New-Item -ItemType Directory -Force -Path $dynamicRoot | Out-Null
    }
    if (-not (Test-Path -LiteralPath $coldOutput)) {
        New-Item -ItemType Directory -Force -Path $coldOutput | Out-Null
        $coldArguments = [Collections.ArrayList]@(
            '--transient', '--config', $config,
            '--mor-transient-generate', '--mor-transient-method', 'local-block-arnoldi',
            '--mor-transient-save', $modelCache,
            '--mor-arnoldi-moments', $ArnoldiMoments,
            '--mor-construction-traces', $ConstructionTraces,
            '--mor-interface-rank', $ConstructionTraceRank,
            '--mor-arnoldi-rank-tolerance', $RankTolerance.ToString('R', $invariant),
            '--mor-arnoldi-second-moment-energy', '1',
            '--mor-arnoldi-second-moment-max-columns', $SecondMomentMaxColumns,
            '--mor-transient-dt', $Dt.ToString('R', $invariant),
            '--mor-transient-t-end', ($ValidationSteps * $Dt).ToString('R', $invariant),
            '--mor-transient-waveform', 'asynchronous_hotspots',
            '--mor-transient-initial-mode', 'ambient',
            '--mor-transient-production',
            '--mor-interface-initial-guess', 'previous',
            '--mor-interface-krylov', $InterfaceSolver,
            '--local-mor-matrix-free-threshold', '0',
            '--max-pcg-iterations', $MaxIterations,
            '--gmres-restart', $Restart,
            '--pcg-tolerance', '1e-10',
            '--mor-adaptive-interface-tolerance', '1e-9',
            '--mor-full-residual-tolerance', '1e-4',
            '--schur-proxy-ring', '1', '--schur-proxy-block-size', '64',
            '--schur-local-solve-threads', $OuterThreads,
            '--schur-local-pardiso-threads', $MklThreads,
            '--schur-proxy-cache', $proxyCache,
            '--output-dir', $coldOutput, '--fast-run'
        )
        if (-not $DisableTemplateReuse) {
            [void]$coldArguments.Add('--mor-local-transient-reuse-identical-subdomains')
        }
        if ($InterfaceSolver -eq 'port-core') {
            [void]$coldArguments.Add('--mor-port-core-cache')
            [void]$coldArguments.Add($portCoreCache)
        }
        if ($EnableResidualFallback) {
            [void]$coldArguments.Add('--mor-full-residual-fallback')
        } else {
            [void]$coldArguments.Add('--no-mor-full-residual-fallback')
        }
        if ($EnrichmentRounds -gt 0) {
            [void]$coldArguments.Add('--local-port-enrichment-rounds')
            [void]$coldArguments.Add($EnrichmentRounds)
        }
        Add-FomComparison $coldArguments
        Invoke-Solver $coldArguments (Join-Path $coldOutput 'run.log')
    }
    $cold = Import-Csv (Join-Path $coldOutput 'local_dynamic_schur_summary.csv')
    $reports += New-ReportRow 'cold_validation' $cold
    if (-not $DisableTemplateReuse -and [int]$cold.reused_instances -lt 7) {
        Write-Warning "Expected at least seven HBM template reuses; observed $($cold.reused_instances)."
    }
}

if ($Stage -eq 'Warm' -or $Stage -eq 'All') {
    $modelFile = Join-Path $modelCache 'local_dynamic_interior_model.bin'
    if (-not (Test-Path -LiteralPath $modelFile -PathType Leaf)) {
        throw "Warm stage requires an existing cold cache: $modelFile"
    }
    $warmOutput = Join-Path $dynamicRoot "warm_${WarmSteps}_steps"
    $portCoreCacheExistedBeforeWarm = Test-Path `
        -LiteralPath $portCoreCache -PathType Leaf
    if ($Force -and (Test-Path -LiteralPath $warmOutput)) {
        Reset-Directory $warmOutput
    } elseif (-not (Test-Path -LiteralPath $warmOutput)) {
        New-Item -ItemType Directory -Force -Path $warmOutput | Out-Null
    }
    $summaryPath = Join-Path $warmOutput 'local_dynamic_schur_summary.csv'
    if (-not (Test-Path -LiteralPath $summaryPath -PathType Leaf)) {
        $warmArguments = [Collections.ArrayList]@(
            '--transient', '--config', $config,
            '--mor-transient-generate', '--mor-transient-method', 'local-block-arnoldi',
            '--mor-transient-load', $modelCache,
            '--mor-arnoldi-moments', $ArnoldiMoments,
            '--mor-construction-traces', $ConstructionTraces,
            '--mor-interface-rank', $ConstructionTraceRank,
            '--mor-arnoldi-rank-tolerance', $RankTolerance.ToString('R', $invariant),
            '--mor-arnoldi-second-moment-energy', '1',
            '--mor-arnoldi-second-moment-max-columns', $SecondMomentMaxColumns,
            '--mor-transient-dt', $Dt.ToString('R', $invariant),
            '--mor-transient-t-end', ($WarmSteps * $Dt).ToString('R', $invariant),
            '--mor-transient-waveform', 'asynchronous_hotspots',
            '--mor-transient-initial-mode', 'ambient',
            '--mor-transient-production',
            '--mor-interface-initial-guess', 'previous',
            '--mor-interface-krylov', $InterfaceSolver,
            '--local-mor-matrix-free-threshold', '0',
            '--max-pcg-iterations', $MaxIterations,
            '--gmres-restart', $Restart,
            '--pcg-tolerance', '1e-10',
            '--mor-adaptive-interface-tolerance', '1e-9',
            '--mor-full-residual-tolerance', '1e-4',
            '--schur-proxy-ring', '1', '--schur-proxy-block-size', '64',
            '--schur-local-solve-threads', $OuterThreads,
            '--schur-local-pardiso-threads', $MklThreads,
            '--schur-proxy-cache', $proxyCache,
            '--output-dir', $warmOutput, '--fast-run'
        )
        if (-not $DisableTemplateReuse) {
            [void]$warmArguments.Add('--mor-local-transient-reuse-identical-subdomains')
        }
        if ($InterfaceSolver -eq 'port-core') {
            [void]$warmArguments.Add('--mor-port-core-cache')
            [void]$warmArguments.Add($portCoreCache)
        }
        if ($EnableResidualFallback) {
            [void]$warmArguments.Add('--mor-full-residual-fallback')
        } else {
            [void]$warmArguments.Add('--no-mor-full-residual-fallback')
        }
        Add-FomComparison $warmArguments
        Invoke-Solver $warmArguments (Join-Path $warmOutput 'run.log')
    }
    $warm = Import-Csv $summaryPath
    if ([int]$warm.descriptor_cache_hit -ne 1 -or
        [int]$warm.reference_cache_hit -ne 1 -or
        [int]$warm.local_model_cache_hit -ne 1) {
        throw 'Warm package15 run missed a descriptor/reference/local-model cache.'
    }
    if ($InterfaceSolver -eq 'port-core' -and
        [int]$warm.port_core_cache_hit -ne 1) {
        $reason = if ($portCoreCacheExistedBeforeWarm) {
            'the existing artifact was incompatible or invalid and was rebuilt'
        } else {
            'no artifact existed, so this run constructed it'
        }
        Write-Warning "Warm package15 port/core cache miss: $reason. Rerun Warm to measure a cache hit."
    }
    $reports += New-ReportRow "warm_${WarmSteps}_steps" $warm
}

if ($reports.Count -gt 0) {
    $reportPath = Join-Path $rootFull 'efficiency_summary.csv'
    $existing = if (Test-Path -LiteralPath $reportPath) {
        @(Import-Csv -LiteralPath $reportPath |
            Where-Object { $_.run -notin @($reports.run) })
    } else { @() }
    @($existing) + @($reports) | Export-Csv -LiteralPath $reportPath -NoTypeInformation
    $reports | Format-Table run,steps,global_dofs,full_interface_dofs,total_local_rank,
        reused_instances,interface_iterations_average,dynamic_online_seconds,
        fom_solve_seconds,online_speedup_vs_fom,accuracy_gate_passed -AutoSize
}
