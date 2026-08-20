param(
    [string]$ResultsDirectory =
        'results\milestone8_global_randomized_small_cases',
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$result = if ([IO.Path]::IsPathRooted($ResultsDirectory)) {
    [IO.Path]::GetFullPath($ResultsDirectory)
} else {
    [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite global-randomized results: $result"
}
if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
    throw "Release executable not found: $exe"
}
$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
New-Item -ItemType Directory -Path $result | Out-Null

function Common-LocalArguments(
    [string]$Config,
    [string]$Output,
    [double]$Dt) {
    return @(
        '--transient', '--config', $Config,
        '--mor-transient-generate',
        '--mor-transient-method', 'local-port-block-arnoldi',
        '--port-basis-method', 'hybrid-randomized',
        '--mor-arnoldi-moments', '1',
        '--mor-interface-rank', '0',
        '--randomized-port-rank', '16',
        '--randomized-oversampling', '5',
        '--randomized-power-iterations', '1',
        '--randomized-seed', '12345',
        '--residual-krylov-max-rank', '4',
        '--residual-krylov-max-sweeps', '2',
        '--residual-krylov-tol', '1e-4',
        '--residual-krylov-block-size', '4',
        '--residual-krylov-probe-mode', 'operator-geometry',
        '--residual-krylov-inner-solver', 'woodbury-exact',
        '--optimal-port-inner-solver', 'woodbury-exact',
        '--optimal-port-inner-tol', '1e-10',
        '--history-compression-method', 'deterministic-rrqr',
        '--history-compression-rank', '64',
        '--history-compression-tolerance', '1e-12',
        '--mor-transient-dt', "$Dt",
        '--mor-transient-t-end', "$Dt",
        '--mor-transient-waveform', 'single_step',
        '--mor-transient-initial-mode', 'ambient',
        '--mor-transient-compare-fom',
        '--output-dir', $Output, '--fast-run'
    )
}

function Accuracy-Row(
    [string]$Case,
    [string]$Method,
    [string]$Composition,
    [int]$Rank,
    [string]$RunDirectory) {
    $summary = Import-Csv -LiteralPath (
        Join-Path $RunDirectory 'local_dynamic_schur_summary.csv') |
        Select-Object -First 1
    $diagnosticsPath = Join-Path $RunDirectory (
        'milestone8_global_randomized_diagnostics.csv')
    $diagnostics = if (Test-Path -LiteralPath $diagnosticsPath) {
        Import-Csv -LiteralPath $diagnosticsPath |
            Select-Object -First 1
    } else {
        $null
    }
    $singleStepPath = Join-Path $RunDirectory (
        'milestone8_global_randomized_single_step.csv')
    $singleStep = if (Test-Path -LiteralPath $singleStepPath) {
        Import-Csv -LiteralPath $singleStepPath |
            Select-Object -First 1
    } else {
        $null
    }
    return [pscustomobject]@{
        case = $Case
        method = $Method
        composition = $Composition
        requested_global_rank = $Rank
        local_port_rank = if ($singleStep) {
            $singleStep.local_port_rank
        } elseif ($Method -eq 'M8.9 Local Adaptive Port') {
            $summary.port_dimension
        } else { 0 }
        global_port_rank = if ($diagnostics) {
            $diagnostics.global_port_rank
        } else { 0 }
        active_port_rank = if ($Method -eq
                'Full Interface Dynamic Schur') {
            $summary.full_interface_dofs
        } else {
            $summary.port_dimension
        }
        full_interface_dofs = $summary.full_interface_dofs
        compression_ratio = if ($diagnostics) {
            $diagnostics.compression_ratio
        } else { 1.0 }
        temperature_relative_l2 = $summary.space_time_relative_l2
        max_error_k = $summary.maximum_absolute_k
        flux_relative_l2 = $summary.maximum_fom_rom_flux_relative_l2
        interface_residual =
            $summary.maximum_interface_relative_residual
        full_residual = $summary.maximum_full_residual
        basis_build_time_s = if ($diagnostics) {
            $diagnostics.basis_build_time_s
        } else { $summary.port_basis_seconds }
        global_schur_apply_count = if ($diagnostics) {
            $diagnostics.global_schur_apply_count
        } else { 0 }
        pardiso_phase33_calls = if ($diagnostics) {
            $diagnostics.pardiso_phase33_calls
        } else { 0 }
        global_rhs_count = if ($diagnostics) {
            $diagnostics.global_rhs_count
        } else { 0 }
        global_inner_iterations = if ($diagnostics) {
            $diagnostics.global_inner_iterations
        } else { 0 }
        mean_solve_time_s = if ($diagnostics) {
            $diagnostics.mean_solve_time_s
        } else { 0.0 }
        orthogonality_error = if ($diagnostics) {
            $diagnostics.orthogonality_error
        } else { 0.0 }
        schur_residual = if ($diagnostics) {
            $diagnostics.schur_residual
        } else { 0.0 }
        target_solve_residual = if ($diagnostics) {
            $diagnostics.target_solve_residual
        } else { 0.0 }
        peak_incremental_memory_bytes = if ($diagnostics) {
            $diagnostics.peak_incremental_memory_bytes
        } else { 0 }
        process_peak_memory_bytes = $summary.peak_working_set_bytes
        corrected = 0
        snapshot_used = 0
        fom_used_for_basis = 0
        pod_used = 0
        svd_used = 0
        training_waveform_used = 0
        algebra_status = if ($diagnostics) {
            $diagnostics.status
        } else { 'baseline' }
        failure_reason = if ($diagnostics) {
            $diagnostics.failure_reason
        } else { '' }
        physical_status = $summary.status
    }
}

$cases = @(
    [pscustomobject]@{
        Name = 'two-cube'
        Config = 'configs\two_cube_parametric_h.txt'
        Dt = 0.1
        SourceMode = 'generalized-dynamic'
        Ranks = @(10, 20, 50)
    },
    [pscustomobject]@{
        Name = 'ten-cube'
        Config = 'configs\ten_cube_parametric_h.txt'
        Dt = 0.1
        SourceMode = 'trace-only'
        Ranks = @(25, 50, 100)
    }
)

$allRows = @()
Push-Location $project
try {
    foreach ($case in $cases) {
        $caseRoot = Join-Path $result $case.Name
        $fullRun = Join-Path $caseRoot 'full_interface'
        $localRun = Join-Path $caseRoot 'local_model'
        $modelPath = Join-Path $localRun 'model'
        New-Item -ItemType Directory -Path $fullRun -Force |
            Out-Null
        & $exe --transient --config $case.Config `
            --mor-transient-generate `
            --mor-transient-method local-block-arnoldi `
            --mor-arnoldi-moments 1 --mor-interface-rank 0 `
            --mor-transient-dt $case.Dt `
            --mor-transient-t-end $case.Dt `
            --mor-transient-waveform single_step `
            --mor-transient-initial-mode ambient `
            --mor-transient-compare-fom `
            --output-dir $fullRun --fast-run
        if ($LASTEXITCODE -ne 0) {
            throw "$($case.Name) full-interface run failed."
        }
        $allRows += Accuracy-Row $case.Name `
            'Full Interface Dynamic Schur' 'full-interface' 0 $fullRun

        New-Item -ItemType Directory -Path $localRun -Force |
            Out-Null
        $localArguments = Common-LocalArguments `
            $case.Config $localRun $case.Dt
        $localArguments += @(
            '--optimal-port-source-mode', $case.SourceMode,
            '--mor-transient-save', $modelPath)
        & $exe @localArguments
        if ($LASTEXITCODE -ne 0) {
            throw "$($case.Name) local model run failed."
        }
        $allRows += Accuracy-Row $case.Name `
            'M8.9 Local Adaptive Port' 'local-only' 0 $localRun

        foreach ($rank in $case.Ranks) {
            foreach ($composition in @('global-only', 'augment-local')) {
                $run = Join-Path $caseRoot (
                    "rank_${rank}_$($composition.Replace('-', '_'))")
                New-Item -ItemType Directory -Path $run -Force |
                    Out-Null
                $arguments = Common-LocalArguments `
                    $case.Config $run $case.Dt
                $arguments += @(
                    '--optimal-port-source-mode', $case.SourceMode,
                    '--global-randomized-schur',
                    '--global-randomized-rank', "$rank",
                    '--global-randomized-composition', $composition,
                    '--global-randomized-seed', '12345',
                    '--global-randomized-inner-max-iters', '1000',
                    '--global-randomized-inner-tol', '1e-10',
                    '--mor-transient-load', $modelPath)
                & $exe @arguments
                if ($LASTEXITCODE -ne 0) {
                    throw (
                        "$($case.Name) rank $rank $composition failed.")
                }
                $diagnosticsPath = Join-Path $run (
                    'milestone8_global_randomized_diagnostics.csv')
                $diagnostics =
                    Import-Csv -LiteralPath $diagnosticsPath |
                    Select-Object -First 1
                if ($diagnostics.status -ne 'passed') {
                    throw (
                        "$($case.Name) rank $rank $composition " +
                        "algebra gate failed: $($diagnostics.status)")
                }
                $allRows += Accuracy-Row $case.Name `
                    'M8.11 Global Randomized Port' `
                    $composition $rank $run
            }
        }
    }
} finally {
    Pop-Location
}

$twoRows = @($allRows | Where-Object { $_.case -eq 'two-cube' })
$tenRows = @($allRows | Where-Object { $_.case -eq 'ten-cube' })
$twoRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $project 'outputs\milestone8_global_randomized_two_cube.csv')
$tenRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $project 'outputs\milestone8_global_randomized_ten_cube.csv')
$allRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $project 'outputs\milestone8_global_vs_local_comparison.csv')
$allRows | Select-Object case,method,composition,
    requested_global_rank,basis_build_time_s,
    global_schur_apply_count,pardiso_phase33_calls,
    global_rhs_count,global_inner_iterations,mean_solve_time_s |
    Export-Csv -NoTypeInformation -Encoding UTF8 (
        Join-Path $project 'outputs\milestone8_global_randomized_timing.csv')
$allRows | Select-Object case,method,composition,
    requested_global_rank,global_port_rank,
    peak_incremental_memory_bytes,process_peak_memory_bytes |
    Export-Csv -NoTypeInformation -Encoding UTF8 (
        Join-Path $project 'outputs\milestone8_global_randomized_memory.csv')
