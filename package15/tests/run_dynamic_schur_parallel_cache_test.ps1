param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$trimSeparators = [char[]]@(
    [IO.Path]::DirectorySeparatorChar, [IO.Path]::AltDirectorySeparatorChar)
$rootFull = [IO.Path]::GetFullPath($Root).TrimEnd($trimSeparators)
$buildRoot = [IO.Path]::GetFullPath((Join-Path $project 'build')).TrimEnd(
    $trimSeparators)
$buildChildPrefix = $buildRoot + [IO.Path]::DirectorySeparatorChar
if ($rootFull.Equals($buildRoot, [StringComparison]::OrdinalIgnoreCase) -or
    -not $rootFull.StartsWith($buildChildPrefix,
        [StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
if (Test-Path -LiteralPath $rootFull) {
    Remove-Item -LiteralPath $rootFull -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null
$cache = Join-Path $rootFull 'dynamic_schur.proxycache'
$model = Join-Path $rootFull 'local_dynamic_model'
$operatorCache = Join-Path $model 'thermal_descriptor.bin'
if (Test-Path -LiteralPath $cache) { Remove-Item -LiteralPath $cache -Force }
if (Test-Path -LiteralPath $operatorCache) {
    Remove-Item -LiteralPath $operatorCache -Force
}

function Invoke-Case(
    [string]$Name, [int]$OuterThreads, [string]$Krylov,
    [string]$ModelAction, [string]$ModelCacheRoot = $model,
    [string]$SaveCacheRoot = '') {
    $output = Join-Path $rootFull $Name
    $modelArguments = switch ($ModelAction) {
        'save' { @('--mor-transient-save', $ModelCacheRoot) }
        'load' { @('--mor-transient-load', $ModelCacheRoot) }
        'load-save' {
            @('--mor-transient-load', $ModelCacheRoot,
                '--mor-transient-save', $SaveCacheRoot)
        }
        default { throw "Unsupported model cache action: $ModelAction" }
    }
    & $Exe --transient --config configs\two_cube_parametric_h.txt `
        --mor-transient-generate --mor-transient-method local-block-arnoldi `
        --mor-arnoldi-moments 3 --mor-transient-dt 0.1 --mor-transient-t-end 0.2 `
        --mor-transient-waveform single_step --mor-transient-initial-mode ambient `
        --mor-interface-initial-guess previous --local-mor-matrix-free-threshold 0 `
        --max-pcg-iterations 200 --gmres-restart 30 --pcg-tolerance 1e-10 `
        --schur-proxy-ring 1 --schur-proxy-block-size 8 `
        --mor-interface-krylov $Krylov `
        --schur-local-solve-threads $OuterThreads --schur-local-pardiso-threads 1 `
        --schur-proxy-cache $cache --output-dir $output --fast-run `
        @modelArguments | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "$Name failed." }
    return Import-Csv (Join-Path $output 'local_dynamic_schur_summary.csv')
}

function Compare-Fields([string]$LeftName, [string]$RightName) {
    $left = Import-Csv (Join-Path $rootFull `
        "$LeftName\local_dynamic_schur_final_temperature.csv")
    $right = Import-Csv (Join-Path $rootFull `
        "$RightName\local_dynamic_schur_final_temperature.csv")
    if ($left.Count -ne $right.Count) { throw 'Temperature vector sizes differ.' }
    [double]$error2 = 0.0
    [double]$reference2 = 0.0
    [double]$maximum = 0.0
    for ($row = 0; $row -lt $left.Count; ++$row) {
        $a = [double]$left[$row].temperature_k
        $b = [double]$right[$row].temperature_k
        $difference = $a - $b
        $error2 += $difference * $difference
        $reference2 += $a * $a
        $maximum = [Math]::Max($maximum, [Math]::Abs($difference))
    }
    $relative = [Math]::Sqrt($error2 / $reference2)
    if ($relative -ge 1e-12 -or $maximum -ge 1e-9) {
        throw "$LeftName/$RightName mismatch: relative=$relative, max=$maximum"
    }
}

Push-Location $project
try {
    $serial = Invoke-Case 'serial_cache_miss' 1 'fgmres' 'save'
    $parallel = Invoke-Case 'parallel_cache_hit' 2 'fgmres' 'load'
    $pcg = Invoke-Case 'protected_pcg' 2 'pcg' 'load'
    if ($serial.status -ne 'success' -or $parallel.status -ne 'success') {
        throw 'Dynamic Schur serial/parallel run did not report success.'
    }
    if ([int]$serial.proxy_matrix_cache_hit -ne 0 -or
        [int]$parallel.proxy_matrix_cache_hit -ne 1) {
        throw 'Dynamic Schur persistent proxy cache miss/hit sequence was not reported.'
    }
    if ([int]$serial.local_model_cache_hit -ne 0 -or
        [int]$parallel.local_model_cache_hit -ne 1 -or
        [double]$serial.local_model_cache_save_seconds -le 0.0 -or
        [double]$parallel.local_model_cache_load_seconds -le 0.0) {
        throw 'Dynamic local-model cache save/load sequence was not reported.'
    }
    if ([int]$serial.reference_cache_hit -ne 0 -or
        [int]$parallel.reference_cache_hit -ne 1 -or
        [double]$serial.reference_cache_save_seconds -le 0.0 -or
        [double]$parallel.reference_cache_load_seconds -le 0.0) {
        throw 'Dynamic reference cache save/load sequence was not reported.'
    }
    if ([int]$serial.descriptor_cache_hit -ne 0 -or
        [int]$parallel.descriptor_cache_hit -ne 1 -or
        [double]$serial.descriptor_cache_save_seconds -le 0.0 -or
        [double]$parallel.descriptor_cache_load_seconds -le 0.0 -or
        [double]$parallel.descriptor_assembly_seconds -ne 0.0) {
        throw 'Thermal descriptor cache miss/hit sequence was not reported.'
    }
    $packageModel = Join-Path $rootFull 'packaged_model_from_load'
    $package = Invoke-Case 'cache_package_from_hit' 1 'fgmres' `
        'load-save' $model $packageModel
    if ($package.status -ne 'success' -or
        [int]$package.descriptor_cache_hit -ne 1 -or
        [double]$package.descriptor_cache_save_seconds -le 0.0 -or
        [int]$package.reference_cache_hit -ne 1 -or
        [double]$package.reference_cache_save_seconds -le 0.0 -or
        [int]$package.local_model_cache_hit -ne 1 -or
        [double]$package.local_model_cache_save_seconds -le 0.0) {
        throw 'Load/save cache packaging did not mirror all cache components.'
    }
    foreach ($required in @('thermal_descriptor.bin',
            'local_dynamic_reference.bin',
            'local_dynamic_interior_model.bin')) {
        if (-not (Test-Path -LiteralPath (Join-Path $packageModel $required))) {
            throw "Load/save cache package is missing $required."
        }
    }
    $incompleteModel = Join-Path $rootFull 'incomplete_cached_runner_model'
    New-Item -ItemType Directory -Force -Path $incompleteModel | Out-Null
    Copy-Item -LiteralPath (Join-Path $model `
        'local_dynamic_interior_model.bin') -Destination $incompleteModel -Force
    $incompleteOutput = Join-Path $rootFull 'incomplete_cached_runner_output'
    $incompleteRejected = $false
    try {
        & (Join-Path $project 'scripts\run_cached_local_dynamic_schur.ps1') `
            -Exe $Exe -Config (Join-Path $project `
                'configs\two_cube_parametric_h.txt') `
            -ModelCache $incompleteModel -ProxyCache $cache `
            -OutputDirectory $incompleteOutput -Steps 1 -Dt 0.1 `
            -OuterThreads 2 -MklThreads 1 -ArnoldiMoments 3 | Out-Null
    } catch {
        $incompleteRejected = $true
        if ($_.Exception.Message -notmatch
            'Required cached model is missing for cached production') {
            throw
        }
    }
    if (-not $incompleteRejected) {
        throw 'Cached production wrapper accepted an incomplete cache.'
    }
    $replaceLoadModel = Join-Path $rootFull 'descriptor_replace_load_model'
    $replaceSaveModel = Join-Path $rootFull 'descriptor_replace_save_model'
    Copy-Item -LiteralPath $model -Destination $replaceLoadModel -Recurse -Force
    Remove-Item -LiteralPath (Join-Path $replaceLoadModel `
        'thermal_descriptor.bin') -Force
    New-Item -ItemType Directory -Force -Path $replaceSaveModel | Out-Null
    Copy-Item -LiteralPath $operatorCache `
        -Destination (Join-Path $replaceSaveModel 'thermal_descriptor.bin') -Force
    $replace = Invoke-Case 'descriptor_cache_replace_save' 1 'fgmres' `
        'load-save' $replaceLoadModel $replaceSaveModel
    if ($replace.status -ne 'success' -or
        [int]$replace.descriptor_cache_hit -ne 0 -or
        [double]$replace.descriptor_cache_save_seconds -le 0.0 -or
        [double]$replace.descriptor_assembly_seconds -le 0.0 -or
        [int]$replace.local_model_cache_hit -ne 1 -or
        [int]$replace.reference_cache_hit -ne 1) {
        throw 'Thermal descriptor cache replacement save was not reported.'
    }
    if ([int]$serial.local_solve_threads -ne 1 -or
        [int]$parallel.local_solve_threads -ne 2 -or
        [int]$parallel.local_pardiso_threads -ne 1) {
        throw 'Dynamic Schur local thread settings were not propagated.'
    }
    if ([int]$serial.interface_iterations_total -ne
        [int]$parallel.interface_iterations_total) {
        throw 'Parallel local solves changed the Dynamic Schur iteration count.'
    }
    $expectedFgmresMatvecs = [int]$serial.interface_iterations_total + 3
    if ([int]$serial.interface_matvecs -ne $expectedFgmresMatvecs -or
        [int]$parallel.interface_matvecs -ne $expectedFgmresMatvecs) {
        throw 'FGMRES repeated a verified or cached initial operator apply.'
    }
    if ([double]$parallel.maximum_interface_relative_residual -ge 1e-9 -or
        [double]$parallel.maximum_full_residual -ge 1e-5) {
        throw 'Parallel Dynamic Schur residual gate failed.'
    }
    if ($pcg.status -ne 'success' -or
        $pcg.interface_solver -ne 'matrix-free-pcg-protected' -or
        $pcg.interface_krylov_actual -ne 'pcg' -or
        [int]$pcg.interface_krylov_fallback_steps -ne 0 -or
        [double]$pcg.maximum_interface_relative_residual -ge 1e-9) {
        throw 'Protected PCG did not pass the exact residual path.'
    }
    Compare-Fields 'serial_cache_miss' 'parallel_cache_hit'
    Compare-Fields 'serial_cache_miss' 'protected_pcg'
    Compare-Fields 'serial_cache_miss' 'cache_package_from_hit'
    Compare-Fields 'serial_cache_miss' 'descriptor_cache_replace_save'

    $wrapperOutput = Join-Path $rootFull 'cached_wrapper_smoke'
    if (Test-Path -LiteralPath $wrapperOutput) {
        Remove-Item -LiteralPath $wrapperOutput -Recurse -Force
    }
    & (Join-Path $project 'scripts\run_cached_local_dynamic_schur.ps1') `
        -Exe $Exe -Config (Join-Path $project 'configs\two_cube_parametric_h.txt') `
        -ModelCache $model -ProxyCache $cache -OutputDirectory $wrapperOutput `
        -Steps 1 -Dt 0.1 -OuterThreads 2 -MklThreads 1 `
        -ArnoldiMoments 3 | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw 'Cached production wrapper smoke test failed.'
    }
    $wrapperSummary = Import-Csv (
        Join-Path $wrapperOutput 'local_dynamic_schur_summary.csv')
    if ($wrapperSummary.status -ne 'success' -or
        [int]$wrapperSummary.reference_cache_hit -ne 1 -or
        [int]$wrapperSummary.local_model_cache_hit -ne 1 -or
        [int]$wrapperSummary.descriptor_cache_hit -ne 1 -or
        [int]$wrapperSummary.native_reduced_history_enabled -ne 1 -or
        [int]$wrapperSummary.native_reduced_history_steps -ne 1 -or
        [double]$wrapperSummary.native_reduced_rhs_seconds -le 0.0 -or
        [double]$wrapperSummary.adaptive_interface_tolerance -ne 1e-9 -or
        [int]$wrapperSummary.residual_fallback_steps -ne 0) {
        throw 'Cached production wrapper did not preserve cache/gate/native-history behavior.'
    }

    $interfaceCoarseCache = Join-Path $rootFull 'interface_operator_coarse_rank4.bin'
    $operatorMissOutput = Join-Path $rootFull 'interface_operator_coarse_miss'
    $operatorHitOutput = Join-Path $rootFull 'interface_operator_coarse_hit'
    & (Join-Path $project 'scripts\run_cached_local_dynamic_schur.ps1') `
        -Exe $Exe -Config (Join-Path $project 'configs\two_cube_parametric_h.txt') `
        -ModelCache $model -ProxyCache $cache `
        -OutputDirectory $operatorMissOutput -Steps 2 -Dt 0.1 `
        -OuterThreads 2 -MklThreads 1 -ArnoldiMoments 3 `
        -OperatorCoarseRank 4 -OperatorCoarseSweeps 2 `
        -OperatorCoarseCache $interfaceCoarseCache | Out-Null
    & (Join-Path $project 'scripts\run_cached_local_dynamic_schur.ps1') `
        -Exe $Exe -Config (Join-Path $project 'configs\two_cube_parametric_h.txt') `
        -ModelCache $model -ProxyCache $cache `
        -OutputDirectory $operatorHitOutput -Steps 2 -Dt 0.1 `
        -OuterThreads 2 -MklThreads 1 -ArnoldiMoments 3 `
        -OperatorCoarseRank 4 -OperatorCoarseSweeps 2 `
        -OperatorCoarseCache $interfaceCoarseCache | Out-Null
    $operatorMiss = Import-Csv (
        Join-Path $operatorMissOutput 'local_dynamic_schur_summary.csv')
    $operatorHit = Import-Csv (
        Join-Path $operatorHitOutput 'local_dynamic_schur_summary.csv')
    if ($operatorMiss.status -ne 'success' -or
        $operatorHit.status -ne 'success' -or
        [int]$operatorMiss.operator_coarse_cache_hit -ne 0 -or
        [int]$operatorHit.operator_coarse_cache_hit -ne 1 -or
        [int]$operatorHit.geometric_coarse_dimension -ne 8 -or
        [int]$operatorHit.operator_coarse_dimension -ne 4 -or
        [int]$operatorHit.coarse_dimension -ne 12 -or
        [int]$operatorHit.interface_predictor_applied_steps -ne 2 -or
        [int]$operatorHit.interface_predictor_accepted_steps -ne 2 -or
        [int]$operatorHit.interface_iterations_total -ge
            [int]$parallel.interface_iterations_total -or
        [int]$operatorHit.interface_matvecs -ne
            ([int]$operatorHit.interface_iterations_total + 3) -or
        [double]$operatorHit.maximum_interface_relative_residual -ge 1e-9 -or
        [double]$operatorHit.maximum_full_residual -ge 1e-5 -or
        $operatorMiss.interface_iterations_total -ne
            $operatorHit.interface_iterations_total -or
        $operatorMiss.maximum_full_residual -ne
            $operatorHit.maximum_full_residual) {
        throw 'Interface operator coarse cold/hot cache regression failed.'
    }

    $adaptiveGuardOutput = Join-Path $rootFull 'adaptive_retry_guard'
    & (Join-Path $project 'scripts\run_cached_local_dynamic_schur.ps1') `
        -Exe $Exe -Config (Join-Path $project 'configs\two_cube_parametric_h.txt') `
        -ModelCache $model -ProxyCache $cache `
        -OutputDirectory $adaptiveGuardOutput -Steps 2 -Dt 0.1 `
        -OuterThreads 2 -MklThreads 1 -ArnoldiMoments 3 `
        -AdaptiveTolerance 1e-4 | Out-Null
    $adaptiveGuardSummary = Import-Csv (
        Join-Path $adaptiveGuardOutput 'local_dynamic_schur_summary.csv')
    if ($adaptiveGuardSummary.status -ne 'success' -or
        [int]$adaptiveGuardSummary.adaptive_interface_retry_steps -lt 1 -or
        [int]$adaptiveGuardSummary.adaptive_interface_retry_iterations -lt 1 -or
        [int]$adaptiveGuardSummary.residual_fallback_steps -ne 0 -or
        [double]$adaptiveGuardSummary.maximum_full_residual -ge 1e-5) {
        throw 'Adaptive interface tolerance did not retry through the strict residual gate.'
    }

    $steadyOutput = Join-Path $rootFull 'steady_exact_reuse'
    & $Exe --transient --config configs\two_cube_parametric_h.txt `
        --mor-transient-generate --mor-transient-method local-block-arnoldi `
        --mor-arnoldi-moments 3 --mor-transient-dt 0.1 `
        --mor-transient-t-end 1.0 `
        --mor-transient-input tests\two_channel_nominal_waveform.csv `
        --mor-transient-initial-mode steady --mor-transient-production `
        --mor-interface-initial-guess previous `
        --mor-adaptive-interface-tolerance 1e-9 `
        --local-mor-matrix-free-threshold 0 --max-pcg-iterations 200 `
        --gmres-restart 30 --pcg-tolerance 1e-10 --schur-proxy-ring 1 `
        --schur-proxy-block-size 8 --mor-interface-krylov fgmres `
        --schur-local-solve-threads 2 --schur-local-pardiso-threads 1 `
        --schur-proxy-cache $cache --mor-transient-load $model `
        --output-dir $steadyOutput --fast-run | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'Exact steady-state reuse case failed.' }
    $steadySummary = Import-Csv (
        Join-Path $steadyOutput 'local_dynamic_schur_summary.csv')
    $steadyTiming = Import-Csv (
        Join-Path $steadyOutput 'local_dynamic_schur_reduced_timing.csv')
    if ($steadySummary.status -ne 'success' -or
        [int]$steadySummary.steady_state_detected_step -lt 1 -or
        [int]$steadySummary.steady_state_reused_steps -lt 1 -or
        [int]$steadySummary.residual_fallback_steps -ne 0 -or
        [double]$steadySummary.full_step_assembly_seconds -ne 0.0 -or
        $steadyTiming[-1].interface_initial_guess -ne 'steady_exact_reuse') {
        throw 'Exact steady-state reuse did not preserve the residual-safe fast path.'
    }

    $mismatchOutput = Join-Path $rootFull 'fingerprint_mismatch'
    $mismatchArguments = @(
        '--transient', '--config', 'configs\two_cube_parametric_h.txt',
        '--mor-transient-generate', '--mor-transient-method',
        'local-block-arnoldi', '--mor-arnoldi-moments', '2',
        '--mor-transient-dt', '0.1', '--mor-transient-t-end', '0.1',
        '--mor-transient-waveform', 'single_step',
        '--mor-transient-load', ('"' + $model + '"'), '--output-dir',
        ('"' + $mismatchOutput + '"'), '--fast-run')
    $mismatchProcess = Start-Process -FilePath $Exe `
        -ArgumentList $mismatchArguments -WorkingDirectory $project `
        -RedirectStandardOutput (Join-Path $rootFull 'mismatch.stdout.txt') `
        -RedirectStandardError (Join-Path $rootFull 'mismatch.stderr.txt') `
        -WindowStyle Hidden -Wait -PassThru
    $mismatchExitCode = $mismatchProcess.ExitCode
    if ($mismatchExitCode -eq 0) {
        throw 'Dynamic local-model cache accepted a basis-parameter mismatch.'
    }

    $operatorMismatchOutput = Join-Path $rootFull 'operator_fingerprint_mismatch'
    $operatorMismatchStdout = Join-Path $rootFull 'operator_mismatch.stdout.txt'
    $operatorMismatchStderr = Join-Path $rootFull 'operator_mismatch.stderr.txt'
    $operatorMismatchArguments = @(
        '--transient', '--config', 'configs\two_cube_parametric_h.txt',
        '--penalty-factor', '16', '--mor-transient-generate',
        '--mor-transient-method', 'local-block-arnoldi',
        '--mor-arnoldi-moments', '3', '--mor-transient-dt', '0.1',
        '--mor-transient-t-end', '0.1', '--mor-transient-waveform',
        'single_step', '--mor-transient-load', ('"' + $model + '"'),
        '--output-dir', ('"' + $operatorMismatchOutput + '"'), '--fast-run')
    $operatorMismatchProcess = Start-Process -FilePath $Exe `
        -ArgumentList $operatorMismatchArguments -WorkingDirectory $project `
        -RedirectStandardOutput $operatorMismatchStdout `
        -RedirectStandardError $operatorMismatchStderr `
        -WindowStyle Hidden -Wait -PassThru
    if ($operatorMismatchProcess.ExitCode -eq 0) {
        throw 'Thermal descriptor cache accepted a penalty mismatch.'
    }
    $operatorMismatchText = (Get-Content -LiteralPath $operatorMismatchStderr `
        -Raw -ErrorAction SilentlyContinue) + (Get-Content -LiteralPath `
        $operatorMismatchStdout -Raw -ErrorAction SilentlyContinue)
    if ($operatorMismatchText -notmatch 'Assembly-input fingerprint mismatch') {
        throw 'Thermal descriptor mismatch did not report a clear refusal.'
    }
} finally {
    Pop-Location
}
