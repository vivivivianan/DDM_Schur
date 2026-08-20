param(
    [string]$Exe = ".\build\Release\SIPGHeatDDM3D.exe",
    [string]$RawRoot = ".\outputs\milestone8_randomized_runs",
    [int]$Threads = 8,
    [switch]$TwoOnly
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $project
try {
    $exePath = (Resolve-Path $Exe).Path
    $rawPath = [System.IO.Path]::GetFullPath($RawRoot)
    New-Item -ItemType Directory -Force -Path $rawPath | Out-Null
    $env:OMP_NUM_THREADS = "$Threads"
    $env:MKL_NUM_THREADS = "$Threads"

    $twoRows = @()
    $tenRows = @()
    $comparisonRows = @()
    $timingRows = @()
    $memoryRows = @()
    $failureRows = @()
    if ($TwoOnly) {
        $tenRows = @(Import-Csv `
            'outputs\milestone8_randomized_transfer_ten_cube.csv')
        $comparisonRows = @(Import-Csv `
            'outputs\milestone8_randomized_vs_optimal_transfer.csv' |
            Where-Object { $_.case -eq 'ten-cube' })
        $timingRows = @(Import-Csv `
            'outputs\milestone8_randomized_timing.csv' |
            Where-Object { $_.case -eq 'ten-cube' })
        $memoryRows = @(Import-Csv `
            'outputs\milestone8_randomized_memory.csv' |
            Where-Object { $_.case -eq 'ten-cube' })
    }

    function Invoke-PortCase(
        [string]$CaseName,
        [string]$Config,
        [string]$Method,
        [int]$Rank,
        [int]$Power,
        [string]$SourceMode,
        [bool]$CompareOptimal) {
        $name = "${CaseName}_${Method}_r${Rank}_q${Power}"
        $output = Join-Path $rawPath $name
        New-Item -ItemType Directory -Force -Path $output | Out-Null
        $moments = if ($CaseName -eq 'two-cube') { 2 } else { 1 }
        $common = @(
            '--transient','--config',$Config,
            '--mor-transient-generate',
            '--mor-transient-method','local-port-block-arnoldi',
            '--mor-arnoldi-moments',"$moments",
            '--optimal-port-inner-solver','woodbury-exact',
            '--optimal-port-inner-tol','1e-10',
            '--optimal-port-inner-refinement-max-iters','3',
            '--optimal-port-inner-refinement-tol','1e-10',
            '--optimal-port-source-mode',$SourceMode,
            '--mor-transient-dt','0.1',
            '--mor-transient-t-end','0.1',
            '--mor-transient-waveform','single_step',
            '--mor-transient-initial-mode','ambient',
            '--output-dir',$output,'--fast-run')
        if ($Method -eq 'randomized-transfer') {
            $arguments = @(
                '--port-basis-method','randomized-transfer',
                '--randomized-port-rank',"$Rank",
                '--randomized-oversampling','5',
                '--randomized-power-iterations',"$Power",
                '--randomized-seed','12345') + $common
            if ($CompareOptimal) {
                $arguments += '--randomized-port-compare-optimal'
            }
        } else {
            $ablation = if ($SourceMode -eq 'trace-only') {
                'trace-transfer-only'
            } else {
                'generalized-transfer-only'
            }
            $arguments = @(
                '--port-basis-method','optimal-transfer',
                '--optimal-port-rank',"$Rank",
                '--optimal-port-rank-mode','fixed',
                '--optimal-port-min-rank',"$Rank",
                '--optimal-port-max-rank',"$Rank",
                '--optimal-port-ablation',$ablation,
                '--optimal-port-eigen-tol','1e-8',
                '--optimal-port-eigen-max-iters','300') + $common
        }
        & $exePath @arguments
        if ($LASTEXITCODE -ne 0) {
            $script:failureRows += [pscustomobject]@{
                case = $CaseName; method = $Method; rank = $Rank
                power_iteration = $Power
                reason = "process_exit_$LASTEXITCODE"
            }
            return
        }
        $summary = Import-Csv (
            Join-Path $output 'local_dynamic_schur_summary.csv')
        if ($Method -eq 'randomized-transfer') {
            $diagnostics = @(Import-Csv (
                Join-Path $output `
                    'randomized_transfer_interface_diagnostics.csv'))
            $basisTiming = Import-Csv (
                Join-Path $output 'randomized_transfer_timing.csv')
            $targetSolves = ($diagnostics |
                Measure-Object target_solve_count -Sum).Sum
            $phaseCalls = ($diagnostics |
                Measure-Object target_solve_phase33_calls -Sum).Sum
            $portRank = ($diagnostics |
                Measure-Object accepted_rank -Sum).Sum
            $basisTime = [double]$basisTiming.total_basis_time_s
            $incrementalMemory = ($diagnostics |
                Measure-Object peak_incremental_memory_bytes -Maximum).Maximum
            $orthogonality = ($diagnostics |
                Measure-Object orthogonality_error -Maximum).Maximum
            $adjoint = ($diagnostics |
                Measure-Object weighted_adjoint_error -Maximum).Maximum
            $targetResidual = ($diagnostics |
                Measure-Object target_residual -Maximum).Maximum
            $projectionError = ($diagnostics |
                Measure-Object basis_error_indicator -Maximum).Maximum
            $script:timingRows += [pscustomobject]@{
                case = $CaseName; method = $Method; rank = $Rank
                power_iteration = $Power
                probe_generation_s =
                    [double]$basisTiming.probe_generation_time_s
                transfer_apply_s =
                    [double]$basisTiming.transfer_apply_time_s
                transpose_apply_s =
                    [double]$basisTiming.transpose_apply_time_s
                qr_s = [double]$basisTiming.qr_time_s
                serialization_s =
                    [double]$basisTiming.basis_serialization_time_s
                offline_total_s = $basisTime
                port_solve_s =
                    [double]$summary.interface_solve_seconds
                recovery_s =
                    [double]$summary.local_recovery_seconds
                full_temperature_reconstruction_s =
                    [double]$summary.local_online_core_seconds
            }
            foreach ($item in $diagnostics) {
                $script:memoryRows += [pscustomobject]@{
                    case = $CaseName; method = $Method
                    rank = $Rank; power_iteration = $Power
                    interface_id = [int]$item.interface_id
                    probe_matrix_bytes =
                        [uint64]$item.probe_matrix_bytes
                    Y_matrix_bytes = [uint64]$item.Y_matrix_bytes
                    qr_workspace_bytes =
                        [uint64]$item.qr_workspace_bytes
                    final_basis_bytes =
                        [uint64]$item.final_basis_bytes
                    peak_incremental_memory_bytes =
                        [uint64]$item.peak_incremental_memory_bytes
                    process_peak_memory_bytes =
                        [uint64]$summary.peak_working_set_bytes
                }
            }
        } else {
            $rankRows = @(Import-Csv (
                Join-Path $output 'optimal_port_rank_by_interface.csv'))
            $inner = @(Import-Csv (
                Join-Path $output 'optimal_port_inner_solver.csv'))
            $operators = @(Import-Csv (
                Join-Path $output 'optimal_port_operator_diagnostics.csv'))
            $basisTiming = Import-Csv (
                Join-Path $output 'optimal_port_timing.csv')
            $targetSolves = ($inner |
                Measure-Object solve_calls -Sum).Sum
            $phaseCalls = $targetSolves
            $portRank = ($rankRows |
                Measure-Object total_port_rank -Sum).Sum
            $basisTime =
                [double]$basisTiming.total_port_offline_seconds
            $incrementalMemory = ($inner |
                Measure-Object peak_incremental_memory_bytes -Maximum).Maximum
            $orthogonality = ($rankRows |
                Measure-Object orthogonality_error -Maximum).Maximum
            $adjoint = ($operators |
                Measure-Object adjoint_relative_error -Maximum).Maximum
            $targetResidual = ($inner |
                Measure-Object max_relative_residual -Maximum).Maximum
            $projectionError = ($rankRows |
                Measure-Object transfer_indicator -Maximum).Maximum
            $script:timingRows += [pscustomobject]@{
                case = $CaseName; method = $Method; rank = $Rank
                power_iteration = -1; probe_generation_s = 0.0
                transfer_apply_s = 0.0; transpose_apply_s = 0.0
                qr_s = 0.0; serialization_s = 0.0
                offline_total_s = $basisTime
                port_solve_s =
                    [double]$summary.interface_solve_seconds
                recovery_s =
                    [double]$summary.local_recovery_seconds
                full_temperature_reconstruction_s =
                    [double]$summary.local_online_core_seconds
            }
            foreach ($item in $inner) {
                $script:memoryRows += [pscustomobject]@{
                    case = $CaseName; method = $Method
                    rank = $Rank; power_iteration = -1
                    interface_id = [int]$item.interface_id
                    probe_matrix_bytes = 0; Y_matrix_bytes = 0
                    qr_workspace_bytes = 0
                    final_basis_bytes = 0
                    peak_incremental_memory_bytes =
                        [uint64]$item.peak_incremental_memory_bytes
                    process_peak_memory_bytes =
                        [uint64]$summary.peak_working_set_bytes
                }
            }
        }
        $row = [pscustomobject]@{
            case = $CaseName; method = $Method; rank = $Rank
            power_iteration = if ($Method -eq 'randomized-transfer') {
                $Power
            } else { -1 }
            target_solve_count = [int]$targetSolves
            target_solve_phase33_calls = [int]$phaseCalls
            basis_build_time_s = $basisTime
            peak_memory_bytes = [uint64]$incrementalMemory
            port_rank = [int]$portRank
            relative_L2 = [double]$summary.space_time_relative_l2
            max_temperature_error =
                [double]$summary.maximum_temperature_error_k
            flux_relative_L2 =
                [double]$summary.maximum_fom_rom_flux_relative_l2
            interface_residual =
                [double]$summary.maximum_interface_relative_residual
            full_residual =
                [double]$summary.maximum_full_residual
            port_projection_error = [double]$projectionError
            orthogonality_error = [double]$orthogonality
            weighted_adjoint_error = [double]$adjoint
            target_residual = [double]$targetResidual
            status = $summary.status
        }
        $script:comparisonRows += $row
        if ($CaseName -eq 'two-cube') {
            $principal = $null
            if ($CompareOptimal) {
                $principal = Import-Csv (
                    Join-Path $output `
                        'randomized_port_subspace_comparison.csv')
            }
            $script:twoRows += [pscustomobject]@{
                case = $CaseName; method = $Method; rank = $Rank
                power_iteration = $Power
                maximum_principal_angle_radians =
                    if ($principal) {
                        [double]$principal.maximum_principal_angle_radians
                    } else { [double]::NaN }
                projector_difference =
                    if ($principal) {
                        [double]$principal.projector_difference
                    } else { [double]::NaN }
                port_projection_error = [double]$projectionError
                relative_L2 = [double]$summary.space_time_relative_l2
                max_temperature_error =
                    [double]$summary.maximum_temperature_error_k
                interface_residual =
                    [double]$summary.maximum_interface_relative_residual
                target_solve_count = [int]$targetSolves
                basis_build_time_s = $basisTime
            }
        } elseif ($Method -eq 'randomized-transfer') {
            $script:tenRows += $row
        }
    }

    foreach ($power in 0,1) {
        Invoke-PortCase 'two-cube' `
            'configs\two_cube_parametric_h.txt' `
            'randomized-transfer' 8 $power `
            'generalized-dynamic' $true
    }
    if (-not $TwoOnly) {
        foreach ($rank in 8,16,32) {
            foreach ($power in 0,1) {
                Invoke-PortCase 'ten-cube' `
                    'configs\ten_cube_parametric_h.txt' `
                    'randomized-transfer' $rank $power `
                    'trace-only' $false
            }
            Invoke-PortCase 'ten-cube' `
                'configs\ten_cube_parametric_h.txt' `
                'optimal-transfer' $rank 0 'trace-only' $false
        }
    }

    $twoRows | Export-Csv -NoTypeInformation -Encoding utf8 `
        'outputs\milestone8_randomized_transfer_two_cube.csv'
    $tenRows | Export-Csv -NoTypeInformation -Encoding utf8 `
        'outputs\milestone8_randomized_transfer_ten_cube.csv'
    $comparisonRows | Export-Csv -NoTypeInformation -Encoding utf8 `
        'outputs\milestone8_randomized_vs_optimal_transfer.csv'
    $comparisonRows | Export-Csv -NoTypeInformation -Encoding utf8 `
        'outputs\m8_randomized_vs_optimal_transfer.csv'
    $timingRows | Export-Csv -NoTypeInformation -Encoding utf8 `
        'outputs\milestone8_randomized_timing.csv'
    $memoryRows | Export-Csv -NoTypeInformation -Encoding utf8 `
        'outputs\milestone8_randomized_memory.csv'
    if ($failureRows.Count -eq 0) {
        @([pscustomobject]@{
            case = ''; method = ''; rank = ''; power_iteration = ''
            reason = 'none'
        }) | Export-Csv -NoTypeInformation -Encoding utf8 `
            'outputs\milestone8_randomized_failures.csv'
    } else {
        $failureRows | Export-Csv -NoTypeInformation -Encoding utf8 `
            'outputs\milestone8_randomized_failures.csv'
        throw "One or more randomized-transfer small-case runs failed."
    }
} finally {
    Pop-Location
}
