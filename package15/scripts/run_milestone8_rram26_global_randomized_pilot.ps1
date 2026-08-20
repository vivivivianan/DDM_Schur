param(
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [Parameter(Mandatory = $true)]
    [string]$AdaptiveModel,
    [Parameter(Mandatory = $true)]
    [string]$M89Accuracy,
    [string]$ResultsDirectory =
        'results\milestone8_rram26_global_randomized_pilot',
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$configPath = (Resolve-Path -LiteralPath $Config).Path
$modelPath = (Resolve-Path -LiteralPath $AdaptiveModel).Path
$baselinePath = (Resolve-Path -LiteralPath $M89Accuracy).Path
$smallGate = Join-Path $project (
    'outputs\milestone8_global_randomized_ten_cube.csv')
$result = if ([IO.Path]::IsPathRooted($ResultsDirectory)) {
    [IO.Path]::GetFullPath($ResultsDirectory)
} else {
    [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite RRAM26 pilot: $result"
}
if (-not (Test-Path -LiteralPath $smallGate)) {
    throw 'Small-case global-randomized gate CSV is missing.'
}
$failedSmall = @(Import-Csv -LiteralPath $smallGate |
    Where-Object {
        $_.method -eq 'M8.11 Global Randomized Port' -and
        $_.algebra_status -ne 'passed'
    })
if ($failedSmall.Count -ne 0) {
    throw (
        'RRAM26 global-randomized pilot is not authorized: ' +
        'a ten-cube algebra gate failed.')
}
$baselineRows = @(Import-Csv -LiteralPath $baselinePath)
$baseline = $baselineRows |
    Where-Object {
        $_.method -like '*M8.9*' -or
        $_.method -like '*Adaptive*'
    } |
    Select-Object -Last 1
if (-not $baseline) {
    throw 'M8.9 RRAM26 accuracy baseline is missing.'
}
$baselineL2 = [double]$baseline.temperature_relative_l2
$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
New-Item -ItemType Directory -Path $result | Out-Null

$rows = @()
$stop = $false
Push-Location $project
try {
    foreach ($rank in @(50, 100, 200)) {
        $rankRows = @()
        foreach ($composition in @('global-only', 'augment-local')) {
            $run = Join-Path $result (
                "rank_${rank}_$($composition.Replace('-', '_'))")
            New-Item -ItemType Directory -Path $run | Out-Null
            & $exe --transient --config $configPath `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method hybrid-randomized `
                --mor-arnoldi-moments 1 --mor-interface-rank 0 `
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
                --history-compression-method deterministic-rrqr `
                --history-compression-rank 64 `
                --history-compression-tolerance 1e-12 `
                --milestone8-adaptive-production `
                --global-randomized-schur `
                --global-randomized-rank $rank `
                --global-randomized-composition $composition `
                --global-randomized-seed 12345 `
                --global-randomized-inner-max-iters 1000 `
                --global-randomized-inner-tol 1e-10 `
                --mor-transient-load $modelPath `
                --mor-transient-dt 0.01 `
                --mor-transient-t-end 0.01 `
                --mor-transient-waveform mixed_frequency `
                --mor-transient-initial-mode ambient `
                --mor-transient-compare-fom `
                --output-dir $run --fast-run
            if ($LASTEXITCODE -ne 0) {
                $rows += [pscustomobject]@{
                    case = 'rram26'
                    composition = $composition
                    global_interface_dofs = 331331
                    requested_global_rank = $rank
                    local_port_rank = 3384
                    global_port_rank = 0
                    active_port_rank = 0
                    compression_ratio = [double]::NaN
                    basis_build_time_s = [double]::NaN
                    global_schur_apply_count = [double]::NaN
                    pardiso_phase33_calls = 0
                    global_rhs_count = 0
                    global_inner_iterations = 0
                    mean_solve_time_s = [double]::NaN
                    peak_incremental_memory_bytes = [double]::NaN
                    process_peak_memory_bytes = [double]::NaN
                    orthogonality_error = [double]::NaN
                    schur_residual = [double]::NaN
                    target_solve_residual = [double]::NaN
                    temperature_relative_l2 = [double]::NaN
                    max_error_k = [double]::NaN
                    flux_relative_l2 = [double]::NaN
                    interface_residual = [double]::NaN
                    full_residual = [double]::NaN
                    temperature_improvement_fraction = [double]::NaN
                    snapshot_used = 0
                    fom_used_for_basis = 0
                    pod_used = 0
                    svd_used = 0
                    corrected = 0
                    status = 'global_randomized_port_failed'
                    failure_reason =
                        'process_exit_before_diagnostics'
                }
                Write-Warning (
                    "RRAM26 rank $rank $composition exited before " +
                    "diagnostics; stopping staged pilot.")
                $stop = $true
                break
            }
            $diagnostics = Import-Csv -LiteralPath (
                Join-Path $run (
                    'milestone8_global_randomized_diagnostics.csv')) |
                Select-Object -First 1
            $singleStepPath = Join-Path $run (
                'milestone8_global_randomized_single_step.csv')
            $single = if (Test-Path -LiteralPath $singleStepPath) {
                Import-Csv -LiteralPath $singleStepPath |
                    Select-Object -First 1
            } else { $null }
            $row = [pscustomobject]@{
                case = 'rram26'
                composition = $composition
                global_interface_dofs =
                    $diagnostics.global_interface_dofs
                requested_global_rank = $rank
                local_port_rank = 3384
                global_port_rank = $diagnostics.global_port_rank
                active_port_rank = if ($composition -eq
                        'augment-local') {
                    3384 + [int]$diagnostics.global_port_rank
                } else {
                    [int]$diagnostics.global_port_rank
                }
                compression_ratio = $diagnostics.compression_ratio
                basis_build_time_s =
                    $diagnostics.basis_build_time_s
                global_schur_apply_count =
                    $diagnostics.global_schur_apply_count
                pardiso_phase33_calls =
                    $diagnostics.pardiso_phase33_calls
                global_rhs_count = $diagnostics.global_rhs_count
                global_inner_iterations =
                    $diagnostics.global_inner_iterations
                mean_solve_time_s = $diagnostics.mean_solve_time_s
                peak_incremental_memory_bytes =
                    $diagnostics.peak_incremental_memory_bytes
                process_peak_memory_bytes =
                    $diagnostics.peak_working_set_bytes
                orthogonality_error =
                    $diagnostics.orthogonality_error
                schur_residual = $diagnostics.schur_residual
                target_solve_residual =
                    $diagnostics.target_solve_residual
                temperature_relative_l2 = if ($single) {
                    $single.temperature_relative_l2
                } else { [double]::NaN }
                max_error_k = if ($single) {
                    $single.max_error_k
                } else { [double]::NaN }
                flux_relative_l2 = if ($single) {
                    $single.flux_relative_l2
                } else { [double]::NaN }
                interface_residual = if ($single) {
                    $single.interface_residual
                } else { [double]::NaN }
                full_residual = if ($single) {
                    $single.full_residual
                } else { [double]::NaN }
                temperature_improvement_fraction = if ($single) {
                    1.0 - [double]$single.temperature_relative_l2 /
                        $baselineL2
                } else { [double]::NaN }
                snapshot_used = 0
                fom_used_for_basis = 0
                pod_used = 0
                svd_used = 0
                corrected = 0
                status = $diagnostics.status
                failure_reason = $diagnostics.failure_reason
            }
            $rows += $row
            $rankRows += $row
            $rows | Export-Csv -NoTypeInformation -Encoding UTF8 (
                Join-Path $result (
                    'milestone8_global_randomized_rram26_pilot.csv'))
            if ($diagnostics.status -ne 'passed') {
                Write-Warning (
                    "Rank $rank $composition failed algebra/resource " +
                    "gate; stopping.")
                $stop = $true
                break
            }
        }
        if ($stop) { break }
        $bestImprovement = ($rankRows |
            Measure-Object -Property temperature_improvement_fraction `
                -Maximum).Maximum
        if ([double]$bestImprovement -lt 0.10) {
            Write-Warning (
                "Rank $rank did not improve RRAM26 temperature L2 " +
                "by at least 10%; larger ranks are not authorized.")
            $stop = $true
            break
        }
    }
} finally {
    Pop-Location
}

$pilotOutput = Join-Path $project (
    'outputs\milestone8_global_randomized_rram26_pilot.csv')
$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $pilotOutput
$existingComparison = Join-Path $project (
    'outputs\milestone8_global_vs_local_comparison.csv')
$comparisonRows = @(
    Import-Csv -LiteralPath $existingComparison |
        Where-Object { $_.case -ne 'rram26' }
)
$comparisonRows += [pscustomobject]@{
    case = 'rram26'
    method = 'M8.9 Local Adaptive Port'
    composition = 'local-only'
    requested_global_rank = 0
    local_port_rank = 3384
    global_port_rank = 0
    active_port_rank = 3384
    full_interface_dofs = 331331
    compression_ratio = 331331.0 / 3384.0
    temperature_relative_l2 = $baseline.temperature_relative_l2
    max_error_k = $baseline.max_error_k
    flux_relative_l2 = $baseline.flux_relative_l2
    interface_residual = $baseline.interface_residual
    full_residual = $baseline.full_residual
    basis_build_time_s = 363.74
    global_schur_apply_count = 0
    pardiso_phase33_calls = 0
    global_rhs_count = 0
    global_inner_iterations = 0
    mean_solve_time_s = 0
    orthogonality_error = 0
    schur_residual = 0
    target_solve_residual = 0
    peak_incremental_memory_bytes = 853.5MB
    process_peak_memory_bytes = $baseline.peak_memory_bytes
    corrected = 0
    snapshot_used = 0
    fom_used_for_basis = 0
    pod_used = 0
    svd_used = 0
    training_waveform_used = 0
    algebra_status = 'baseline'
    failure_reason = ''
    physical_status = $baseline.status
}
foreach ($row in $rows) {
    $comparisonRows += [pscustomobject]@{
        case = 'rram26'
        method = 'M8.11 Global Randomized Port'
        composition = $row.composition
        requested_global_rank = $row.requested_global_rank
        local_port_rank = $row.local_port_rank
        global_port_rank = $row.global_port_rank
        active_port_rank = $row.active_port_rank
        full_interface_dofs = $row.global_interface_dofs
        compression_ratio = $row.compression_ratio
        temperature_relative_l2 = $row.temperature_relative_l2
        max_error_k = $row.max_error_k
        flux_relative_l2 = $row.flux_relative_l2
        interface_residual = $row.interface_residual
        full_residual = $row.full_residual
        basis_build_time_s = $row.basis_build_time_s
        global_schur_apply_count = $row.global_schur_apply_count
        pardiso_phase33_calls = $row.pardiso_phase33_calls
        global_rhs_count = $row.global_rhs_count
        global_inner_iterations = $row.global_inner_iterations
        mean_solve_time_s = $row.mean_solve_time_s
        orthogonality_error = $row.orthogonality_error
        schur_residual = $row.schur_residual
        target_solve_residual = $row.target_solve_residual
        peak_incremental_memory_bytes =
            $row.peak_incremental_memory_bytes
        process_peak_memory_bytes = $row.process_peak_memory_bytes
        corrected = 0
        snapshot_used = 0
        fom_used_for_basis = 0
        pod_used = 0
        svd_used = 0
        training_waveform_used = 0
        algebra_status = $row.status
        failure_reason = $row.failure_reason
        physical_status = if ($row.status -eq 'passed') {
            'completed'
        } else {
            'not_run'
        }
    }
}
$comparisonRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
    $existingComparison)

$timingPath = Join-Path $project (
    'outputs\milestone8_global_randomized_timing.csv')
$timingRows = @(
    Import-Csv -LiteralPath $timingPath |
        Where-Object { $_.case -ne 'rram26' }
)
$timingRows += $comparisonRows |
    Where-Object { $_.case -eq 'rram26' } |
    Select-Object case,method,composition,requested_global_rank,
        basis_build_time_s,global_schur_apply_count,
        pardiso_phase33_calls,global_rhs_count,
        global_inner_iterations,mean_solve_time_s
$timingRows | Export-Csv -NoTypeInformation -Encoding UTF8 $timingPath

$memoryPath = Join-Path $project (
    'outputs\milestone8_global_randomized_memory.csv')
$memoryRows = @(
    Import-Csv -LiteralPath $memoryPath |
        Where-Object { $_.case -ne 'rram26' }
)
$memoryRows += $comparisonRows |
    Where-Object { $_.case -eq 'rram26' } |
    Select-Object case,method,composition,requested_global_rank,
        global_port_rank,peak_incremental_memory_bytes,
        process_peak_memory_bytes
$memoryRows | Export-Csv -NoTypeInformation -Encoding UTF8 $memoryPath
