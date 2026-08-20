param(
    [ValidateSet('smoke', 'medium')][string]$Profile = 'medium',
    [ValidateSet('Generate', 'Exact', 'Cold', 'Warm', 'All')]
    [string]$Stage = 'All',
    [string]$Exe = '',
    [string]$Root = '',
    [ValidateRange(1, 12)][int]$ArnoldiMoments = 1,
    [ValidateRange(0, 1000000)][int]$SecondMomentMaxColumns = 0,
    [ValidateRange(1.0e-14, 1.0e-2)][double]$RankTolerance = 1.0e-10,
    [ValidateRange(1, 1000)][int]$ValidationSteps = 10,
    [ValidateRange(1, 1000000)][int]$WarmSteps = 100,
    [ValidateRange(1.0e-12, 1.0e6)][double]$Dt = 0.05,
    [ValidateRange(1, 1024)][int]$OuterThreads = 8,
    [ValidateRange(1, 1024)][int]$MklThreads = 1,
    [ValidateRange(1, 1000000)][int]$Restart = 100,
    [ValidateRange(1, 1000000)][int]$MaxIterations = 500,
    [switch]$DisableTemplateReuse,
    [switch]$NoFomComparison,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
if ([string]::IsNullOrWhiteSpace($Exe)) {
    $Exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
}
if ([string]::IsNullOrWhiteSpace($Root)) {
    $Root = Join-Path $project (
        "build\package15_benchmark_${Profile}_m${ArnoldiMoments}_c${SecondMomentMaxColumns}")
}
$exeFull = [IO.Path]::GetFullPath($Exe)
$rootFull = [IO.Path]::GetFullPath($Root)
$buildRoot = [IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Benchmark output must remain under build: $rootFull"
}
if (-not (Test-Path -LiteralPath $exeFull -PathType Leaf)) {
    throw "Release executable not found: $exeFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

function Reset-Directory([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    if (-not $full.StartsWith($rootFull, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to reset directory outside benchmark root: $full"
    }
    if (Test-Path -LiteralPath $full) {
        $resolved = [IO.Path]::GetFullPath((Resolve-Path -LiteralPath $full).Path)
        if (-not $resolved.StartsWith($rootFull, [StringComparison]::OrdinalIgnoreCase)) {
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
    try {
        $env:OMP_NUM_THREADS = [string]$OuterThreads
        $env:MKL_NUM_THREADS = [string]$MklThreads
        $env:MKL_DYNAMIC = 'FALSE'
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
    }
}

function Add-FomComparison([System.Collections.ArrayList]$Arguments) {
    if (-not $NoFomComparison) {
        [void]$Arguments.Add('--mor-transient-compare-fom-summary-only')
    }
}

function New-ReportRow([string]$Run, [pscustomobject]$Summary) {
    $online = [double]$Summary.local_online_core_seconds
    $fomSolve = [double]$Summary.fom_solve_seconds
    $dynamicTotal = [double]$Summary.total_seconds
    $fomFactor = [double]$Summary.fom_factor_seconds
    $onlineSpeedup = if ($online -gt 0.0 -and $fomSolve -gt 0.0) {
        $fomSolve / $online
    } else { [double]::NaN }
    $coldSpeedup = if ($dynamicTotal -gt 0.0 -and ($fomFactor + $fomSolve) -gt 0.0) {
        ($fomFactor + $fomSolve) / $dynamicTotal
    } else { [double]::NaN }
    $accuracyPassed =
        [int]$Summary.residual_fallback_steps -eq 0 -and
        [double]$Summary.space_time_relative_l2 -le 1.0e-4 -and
        [double]$Summary.maximum_temperature_error_k -le 0.05 -and
        [double]$Summary.maximum_full_residual -le 1.0e-4
    return [pscustomobject]@{
        run = $Run
        profile = $Profile
        status = [string]$Summary.status
        steps = [int]$Summary.steps
        global_dofs = [int64]$Summary.global_dofs
        full_interface_dofs = [int64]$Summary.full_interface_dofs
        total_local_rank = [int]$Summary.total_local_rank
        unique_templates = [int]$Summary.unique_templates
        reused_instances = [int]$Summary.reused_instances
        descriptor_cache_hit = [int]$Summary.descriptor_cache_hit
        reference_cache_hit = [int]$Summary.reference_cache_hit
        local_model_cache_hit = [int]$Summary.local_model_cache_hit
        interface_iterations_total = [int]$Summary.interface_iterations_total
        interface_iterations_maximum = [int]$Summary.interface_iterations_maximum
        residual_fallback_steps = [int]$Summary.residual_fallback_steps
        relative_l2 = [double]$Summary.space_time_relative_l2
        maximum_temperature_error_k = [double]$Summary.maximum_temperature_error_k
        maximum_full_residual = [double]$Summary.maximum_full_residual
        fom_factor_seconds = $fomFactor
        fom_solve_seconds = $fomSolve
        dynamic_online_seconds = $online
        dynamic_total_seconds = $dynamicTotal
        online_speedup_vs_fom = $onlineSpeedup
        end_to_end_speedup_vs_fom = $coldSpeedup
        peak_working_set_gib = [double]$Summary.peak_working_set_bytes / 1GB
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
$coldOutput = Join-Path $dynamicRoot 'cold_validation'

if ($Stage -eq 'Cold' -or $Stage -eq 'All') {
    if ($Force -and (Test-Path -LiteralPath $dynamicRoot)) {
        Reset-Directory $dynamicRoot
    } else {
        New-Item -ItemType Directory -Force -Path $dynamicRoot | Out-Null
    }
    if (-not (Test-Path -LiteralPath $coldOutput)) {
        New-Item -ItemType Directory -Force -Path $coldOutput | Out-Null
        $invariant = [Globalization.CultureInfo]::InvariantCulture
        $coldArguments = [Collections.ArrayList]@(
            '--transient', '--config', $config,
            '--mor-transient-generate', '--mor-transient-method', 'local-block-arnoldi',
            '--mor-transient-save', $modelCache,
            '--mor-arnoldi-moments', $ArnoldiMoments,
            '--mor-arnoldi-rank-tolerance', $RankTolerance.ToString('R', $invariant),
            '--mor-arnoldi-second-moment-energy', '1',
            '--mor-arnoldi-second-moment-max-columns', $SecondMomentMaxColumns,
            '--mor-transient-dt', $Dt.ToString('R', $invariant),
            '--mor-transient-t-end', ($ValidationSteps * $Dt).ToString('R', $invariant),
            '--mor-transient-waveform', 'asynchronous_hotspots',
            '--mor-transient-initial-mode', 'ambient',
            '--mor-transient-production',
            '--mor-interface-initial-guess', 'previous',
            '--mor-interface-krylov', 'fgmres',
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
    if ($Force -and (Test-Path -LiteralPath $warmOutput)) {
        Reset-Directory $warmOutput
    } elseif (-not (Test-Path -LiteralPath $warmOutput)) {
        New-Item -ItemType Directory -Force -Path $warmOutput | Out-Null
    }
    $summaryPath = Join-Path $warmOutput 'local_dynamic_schur_summary.csv'
    if (-not (Test-Path -LiteralPath $summaryPath -PathType Leaf)) {
        $invariant = [Globalization.CultureInfo]::InvariantCulture
        $warmArguments = [Collections.ArrayList]@(
            '--transient', '--config', $config,
            '--mor-transient-generate', '--mor-transient-method', 'local-block-arnoldi',
            '--mor-transient-load', $modelCache,
            '--mor-arnoldi-moments', $ArnoldiMoments,
            '--mor-arnoldi-rank-tolerance', $RankTolerance.ToString('R', $invariant),
            '--mor-arnoldi-second-moment-energy', '1',
            '--mor-arnoldi-second-moment-max-columns', $SecondMomentMaxColumns,
            '--mor-transient-dt', $Dt.ToString('R', $invariant),
            '--mor-transient-t-end', ($WarmSteps * $Dt).ToString('R', $invariant),
            '--mor-transient-waveform', 'asynchronous_hotspots',
            '--mor-transient-initial-mode', 'ambient',
            '--mor-transient-production',
            '--mor-interface-initial-guess', 'previous',
            '--mor-interface-krylov', 'fgmres',
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
        Add-FomComparison $warmArguments
        Invoke-Solver $warmArguments (Join-Path $warmOutput 'run.log')
    }
    $warm = Import-Csv $summaryPath
    if ([int]$warm.descriptor_cache_hit -ne 1 -or
        [int]$warm.reference_cache_hit -ne 1 -or
        [int]$warm.local_model_cache_hit -ne 1) {
        throw 'Warm package15 run missed a descriptor/reference/local-model cache.'
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
        reused_instances,interface_iterations_total,dynamic_online_seconds,
        fom_solve_seconds,online_speedup_vs_fom,accuracy_gate_passed -AutoSize
}
