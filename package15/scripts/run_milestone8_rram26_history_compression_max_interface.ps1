param(
    [string]$Exe = ".\build\Release\SIPGHeatDDM3D.exe",
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [string]$Audit =
        ".\outputs\milestone8_large_case_topology_audit.csv",
    [string]$RawRoot =
        ".\outputs\milestone8_history_compression_rram26_max_runs",
    [string]$Aggregate =
        ".\outputs\milestone8_history_compression.csv",
    [int]$Threads = 8,
    [switch]$AggregateOnly
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $project
try {
    $exePath = (Resolve-Path $Exe).Path
    $configPath = (Resolve-Path $Config).Path
    $auditPath = (Resolve-Path $Audit).Path
    $rawPath = [System.IO.Path]::GetFullPath($RawRoot)
    $aggregatePath = [System.IO.Path]::GetFullPath($Aggregate)
    New-Item -ItemType Directory -Force -Path $rawPath | Out-Null
    $env:OMP_NUM_THREADS = "$Threads"
    $env:MKL_NUM_THREADS = "$Threads"

    $rows = @()
    foreach ($historyRank in 16,32,64,128) {
        $output = Join-Path $rawPath "rank_$historyRank"
        if (-not $AggregateOnly) {
            & $exePath --transient --config $configPath `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method hybrid-randomized `
                --mor-arnoldi-moments 1 `
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
                --history-compression-rank $historyRank `
                --history-compression-tolerance 1e-12 `
                --history-compression-max-interface-pilot `
                --optimal-port-topology-audit-csv $auditPath `
                --mor-transient-dt 0.1 `
                --mor-transient-t-end 0.1 `
                --mor-transient-waveform single_step `
                --mor-transient-initial-mode ambient `
                --output-dir $output --fast-run
            if ($LASTEXITCODE -ne 0) {
                throw "RRAM26 maximum-interface history rank " +
                    "$historyRank process failed."
            }
        }

        $pilot = @(Import-Csv (
            Join-Path $output `
                'milestone8_rram26_history_compression_max_interface.csv'))
        $stop = Import-Csv (
            Join-Path $output `
                'milestone8_history_compression_max_pilot_stop.csv')
        if ($pilot.Count -ne 1 -or
            $pilot[0].selection -ne 'maximum' -or
            [int]$pilot[0].interface_id -ne 24 -or
            [int]$pilot[0].requested_history_rank -ne $historyRank -or
            $pilot[0].history_compression_method -ne
                'deterministic-rrqr' -or
            [int]$stop.interfaces_requested -ne 1 -or
            [int]$stop.transient_advanced -ne 0 -or
            [int]$stop.full_field_read -ne 0 -or
            [int]$stop.snapshot_used -ne 0) {
            throw "RRAM26 maximum-interface history rank " +
                "$historyRank scope/provenance failed."
        }
        $item = $pilot[0]
        $timePassed =
            [double]$item.total_basis_build_time_s -lt 120.0
        $memoryPassed =
            [uint64]$item.incremental_memory_bytes -lt
                [uint64]1073741824
        $status =
            if ($item.status -notlike 'success*') {
                $item.status
            } elseif (-not $timePassed) {
                'basis_build_time_gate_failed'
            } elseif (-not $memoryPassed) {
                'incremental_memory_gate_failed'
            } else {
                'success'
            }
        $rows += [pscustomobject][ordered]@{
            case = 'RRAM26'
            selection = 'maximum'
            interface_id = [int]$item.interface_id
            configuration = "history-rank-$historyRank"
            compression_method = $item.history_compression_method
            requested_history_rank =
                [int]$item.requested_history_rank
            raw_history_channels =
                [int]$item.raw_history_channels
            active_history_channels =
                [int]$item.active_history_channels
            compressed_history_rank =
                [int]$item.compressed_history_rank
            deflated_history_channels =
                [int]$item.deflated_history_channels
            history_target_rhs = [int]$item.history_target_rhs
            history_compression_relative_error_max =
                [double]$item.history_compression_relative_error
            history_compression_time_s =
                [double]$item.history_compression_time_s
            history_compression_workspace_bytes =
                [uint64]$item.history_compression_workspace_bytes
            history_compression_fingerprint =
                $item.history_compression_fingerprint
            mandatory_rank = [int]$item.mandatory_rank
            requested_randomized_rank =
                [int]$item.requested_randomized_rank
            accepted_randomized_rank =
                [int]$item.accepted_randomized_rank
            requested_enrichment_rank =
                [int]$item.requested_enrichment_rank
            accepted_enrichment_rank =
                [int]$item.accepted_enrichment_rank
            total_port_rank = [int]$item.total_port_rank
            randomized_target_solve_count =
                [int]$item.randomized_target_solve_count
            residual_target_solve_count =
                [int]$item.residual_target_solve_count
            total_target_solve_count =
                [int]$item.total_target_solve_count
            target_solve_phase33_calls =
                [int]$item.target_solve_phase33_calls
            basis_build_time_s =
                [double]$item.total_basis_build_time_s
            incremental_memory_bytes =
                [uint64]$item.incremental_memory_bytes
            process_peak_memory_bytes =
                [uint64]$item.process_peak_memory_bytes
            temperature_relative_L2 = $null
            temperature_relative_change_vs_uncompressed = $null
            flux_relative_L2 = $null
            flux_relative_change_vs_uncompressed = $null
            interface_residual = $null
            interface_residual_relative_change_vs_uncompressed = $null
            full_residual = $null
            full_residual_relative_change_vs_uncompressed = $null
            target_residual = [double]$item.target_residual
            weighted_adjoint_error =
                [double]$item.weighted_adjoint_error
            time_gate_passed = $timePassed
            memory_gate_passed = $memoryPassed
            temperature_gate_passed = $null
            flux_gate_passed = $null
            interface_residual_gate_passed = $null
            full_residual_gate_passed = $null
            accuracy_status = 'accuracy_not_run_by_scope'
            status = $status
            transient_advanced = 0
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
            Where-Object { $_.case -ne 'RRAM26' })
    }
    $combined = @($preserved) + @($rows)
    $combined | Export-Csv -NoTypeInformation -Encoding utf8 `
        $aggregatePath
} finally {
    Pop-Location
}
