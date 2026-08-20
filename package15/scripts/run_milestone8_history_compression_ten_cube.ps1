param(
    [string]$Exe = ".\build\Release\SIPGHeatDDM3D.exe",
    [string]$RawRoot =
        ".\outputs\milestone8_history_compression_ten_cube_runs",
    [string]$Aggregate =
        ".\outputs\milestone8_history_compression.csv",
    [int]$Threads = 8,
    [double]$AccuracyRelativeSlack = 1.0e-3,
    [switch]$AggregateOnly
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $project
try {
    $exePath = (Resolve-Path $Exe).Path
    $rawPath = [System.IO.Path]::GetFullPath($RawRoot)
    $aggregatePath = [System.IO.Path]::GetFullPath($Aggregate)
    New-Item -ItemType Directory -Force -Path $rawPath | Out-Null
    $env:OMP_NUM_THREADS = "$Threads"
    $env:MKL_NUM_THREADS = "$Threads"

    $reference = $null
    $referenceCsv =
        'outputs\milestone8_hybrid_rank_accuracy.csv'
    if (Test-Path $referenceCsv) {
        $reference = Import-Csv $referenceCsv |
            Where-Object {
                $_.case -eq 'ten-cube' -and
                $_.configuration -eq
                    'mandatory+randomized+residual'
            } |
            Select-Object -First 1
    }

    function Test-NoRegression(
        [double]$Value,
        [double]$ReferenceValue) {
        $absoluteSlack = 1.0e-14
        return $Value -le (
            $ReferenceValue * (1.0 + $AccuracyRelativeSlack) +
            $absoluteSlack)
    }

    $rows = @()
    foreach ($historyRank in 16,32,64,128) {
        $output = Join-Path $rawPath "rank_$historyRank"
        if (-not $AggregateOnly) {
            & $exePath --transient `
                --config configs\ten_cube_parametric_h.txt `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method hybrid-randomized `
                --mor-arnoldi-moments 1 `
                --randomized-port-rank 8 `
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
                --history-compression-rank $historyRank `
                --history-compression-tolerance 1e-12 `
                --mor-transient-dt 0.1 `
                --mor-transient-t-end 0.1 `
                --mor-transient-waveform single_step `
                --mor-transient-initial-mode ambient `
                --output-dir $output --fast-run
            if ($LASTEXITCODE -ne 0) {
                throw "ten-cube history rank $historyRank failed."
            }
        }

        $summary = Import-Csv (
            Join-Path $output 'local_dynamic_schur_summary.csv')
        $diagnostics = @(Import-Csv (
            Join-Path $output `
                'residual_krylov_interface_diagnostics.csv'))
        $randomized = @(Import-Csv (
            Join-Path $output `
                'randomized_transfer_interface_diagnostics.csv'))
        $portRanks = @(Import-Csv (
            Join-Path $output 'local_port_rank_by_interface.csv'))
        if ($diagnostics.Count -eq 0 -or
            @($diagnostics | Where-Object {
                $_.history_compression_method -ne
                    'deterministic-rrqr' -or
                [int]$_.requested_history_rank -ne $historyRank
            }).Count -gt 0) {
            throw "ten-cube history rank $historyRank diagnostics invalid."
        }

        $temperature = [double]$summary.space_time_relative_l2
        $flux =
            [double]$summary.maximum_fom_rom_flux_relative_l2
        $interfaceResidual =
            [double]$summary.maximum_interface_relative_residual
        $fullResidual =
            [double]$summary.maximum_full_residual
        $temperatureGate = $null
        $fluxGate = $null
        $interfaceGate = $null
        $fullGate = $null
        $accuracyStatus = 'reference_missing'
        $temperatureChange = $null
        $fluxChange = $null
        $interfaceChange = $null
        $fullChange = $null
        if ($null -ne $reference) {
            $referenceTemperature =
                [double]$reference.temperature_relative_L2
            $referenceFlux =
                [double]$reference.flux_relative_L2
            $referenceInterface =
                [double]$reference.interface_residual
            $referenceFull =
                [double]$reference.full_residual
            $temperatureGate =
                Test-NoRegression $temperature $referenceTemperature
            $fluxGate = Test-NoRegression $flux $referenceFlux
            $interfaceGate =
                Test-NoRegression $interfaceResidual $referenceInterface
            $fullGate =
                Test-NoRegression $fullResidual $referenceFull
            $temperatureChange =
                ($temperature - $referenceTemperature) /
                [Math]::Max([Math]::Abs($referenceTemperature), 1e-300)
            $fluxChange =
                ($flux - $referenceFlux) /
                [Math]::Max([Math]::Abs($referenceFlux), 1e-300)
            $interfaceChange =
                ($interfaceResidual - $referenceInterface) /
                [Math]::Max([Math]::Abs($referenceInterface), 1e-300)
            $fullChange =
                ($fullResidual - $referenceFull) /
                [Math]::Max([Math]::Abs($referenceFull), 1e-300)
            $accuracyStatus =
                if ($temperatureGate -and $fluxGate -and
                    $interfaceGate -and $fullGate) {
                    'passed'
                } else {
                    'accuracy_gate_failed'
                }
        }

        $peakIncremental = (
            $diagnostics |
            Measure-Object peak_incremental_memory_bytes -Maximum
        ).Maximum
        $targetResidual = [Math]::Max(
            [double](($diagnostics |
                Measure-Object target_residual -Maximum).Maximum),
            [double](($randomized |
                Measure-Object target_residual -Maximum).Maximum))
        $adjoint = [Math]::Max(
            [double](($diagnostics |
                Measure-Object weighted_adjoint_error -Maximum).Maximum),
            [double](($randomized |
                Measure-Object weighted_adjoint_error -Maximum).Maximum))
        $builderPassed =
            @($diagnostics | Where-Object {
                $_.status -notlike 'success*'
            }).Count -eq 0 -and
            @($randomized | Where-Object {
                $_.status -notlike 'success*'
            }).Count -eq 0
        $timePassed =
            [double]$summary.port_basis_seconds -lt 120.0
        $memoryPassed =
            [uint64]$peakIncremental -lt [uint64]1073741824
        $status =
            if (-not $builderPassed) {
                'basis_failed'
            } elseif (-not $timePassed) {
                'basis_build_time_gate_failed'
            } elseif (-not $memoryPassed) {
                'incremental_memory_gate_failed'
            } elseif ($accuracyStatus -eq 'accuracy_gate_failed') {
                'accuracy_gate_failed'
            } else {
                'success'
            }

        $fingerprints = (
            $diagnostics |
            Sort-Object {[int]$_.interface_id} |
            ForEach-Object {
                $_.history_compression_fingerprint
            }) -join ';'
        $rows += [pscustomobject][ordered]@{
            case = 'ten-cube'
            selection = 'all-interfaces'
            interface_id = 'all'
            configuration = "history-rank-$historyRank"
            compression_method = 'deterministic-rrqr'
            requested_history_rank = $historyRank
            raw_history_channels =
                [int](($diagnostics |
                    Measure-Object raw_history_channels -Sum).Sum)
            active_history_channels =
                [int](($diagnostics |
                    Measure-Object active_history_channels -Sum).Sum)
            compressed_history_rank =
                [int](($diagnostics |
                    Measure-Object compressed_history_rank -Sum).Sum)
            deflated_history_channels =
                [int](($diagnostics |
                    Measure-Object deflated_history_channels -Sum).Sum)
            history_target_rhs =
                [int](($diagnostics |
                    Measure-Object history_target_rhs -Sum).Sum)
            history_compression_relative_error_max =
                [double](($diagnostics |
                    Measure-Object `
                        history_compression_relative_error -Maximum).Maximum)
            history_compression_time_s =
                [double](($diagnostics |
                    Measure-Object history_compression_time_s -Sum).Sum)
            history_compression_workspace_bytes =
                [uint64](($diagnostics |
                    Measure-Object `
                        history_compression_workspace_bytes -Maximum).Maximum)
            history_compression_fingerprint = $fingerprints
            mandatory_rank =
                [int](($diagnostics |
                    Measure-Object mandatory_rank_total -Sum).Sum)
            requested_randomized_rank = 8
            accepted_randomized_rank =
                [int](($diagnostics |
                    Measure-Object accepted_randomized_rank -Sum).Sum)
            requested_enrichment_rank = 4
            accepted_enrichment_rank =
                [int](($diagnostics |
                    Measure-Object accepted_enrichment_rank -Sum).Sum)
            total_port_rank =
                [int](($portRanks |
                    Measure-Object total_port_rank -Sum).Sum)
            randomized_target_solve_count =
                [int](($randomized |
                    Measure-Object target_solve_count -Sum).Sum)
            residual_target_solve_count =
                [int](($diagnostics |
                    Measure-Object target_solve_count -Sum).Sum)
            total_target_solve_count =
                [int](($randomized |
                    Measure-Object target_solve_count -Sum).Sum) +
                [int](($diagnostics |
                    Measure-Object target_solve_count -Sum).Sum)
            target_solve_phase33_calls =
                [int](($randomized |
                    Measure-Object target_solve_phase33_calls -Sum).Sum) +
                [int](($diagnostics |
                    Measure-Object target_solve_phase33_calls -Sum).Sum)
            basis_build_time_s = [double]$summary.port_basis_seconds
            incremental_memory_bytes = [uint64]$peakIncremental
            process_peak_memory_bytes =
                [uint64]$summary.peak_working_set_bytes
            temperature_relative_L2 = $temperature
            temperature_relative_change_vs_uncompressed =
                $temperatureChange
            flux_relative_L2 = $flux
            flux_relative_change_vs_uncompressed = $fluxChange
            interface_residual = $interfaceResidual
            interface_residual_relative_change_vs_uncompressed =
                $interfaceChange
            full_residual = $fullResidual
            full_residual_relative_change_vs_uncompressed = $fullChange
            target_residual = $targetResidual
            weighted_adjoint_error = $adjoint
            time_gate_passed = $timePassed
            memory_gate_passed = $memoryPassed
            temperature_gate_passed = $temperatureGate
            flux_gate_passed = $fluxGate
            interface_residual_gate_passed = $interfaceGate
            full_residual_gate_passed = $fullGate
            accuracy_status = $accuracyStatus
            status = $status
            transient_advanced = 1
            full_field_read = 0
            snapshot_used = 0
            fom_used_for_basis = 0
            pod_used = 0
            svd_used = 0
        }
    }

    $preserved = @()
    if (Test-Path $aggregatePath) {
        $preserved = @(Import-Csv $aggregatePath |
            Where-Object { $_.case -ne 'ten-cube' })
    }
    $combined = @($preserved) + @($rows)
    $combined | Export-Csv -NoTypeInformation -Encoding utf8 `
        $aggregatePath
} finally {
    Pop-Location
}
