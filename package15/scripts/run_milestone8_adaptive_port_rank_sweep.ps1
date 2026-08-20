param(
    [string]$Exe = ".\build\Release\SIPGHeatDDM3D.exe",
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [string]$DifficultInterfaces =
        ".\outputs\milestone8_difficult_interfaces.csv",
    [string]$TopologyAudit =
        ".\outputs\milestone8_large_case_topology_audit.csv",
    [string]$RawRoot =
        ".\outputs\milestone8_adaptive_port_rank_runs",
    [string]$Aggregate =
        ".\outputs\milestone8_adaptive_port_rank_sweep.csv",
    [int]$Threads = 8,
    [switch]$Resume,
    [switch]$AggregateOnly
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path

function Resolve-ProjectPath([string]$Path) {
    if ([IO.Path]::IsPathRooted($Path)) {
        return [IO.Path]::GetFullPath($Path)
    }
    return [IO.Path]::GetFullPath((Join-Path $project $Path))
}

$exePath = (Resolve-Path -LiteralPath (
    Resolve-ProjectPath $Exe)).Path
$configPath = (Resolve-Path -LiteralPath (
    Resolve-ProjectPath $Config)).Path
$difficultPath = (Resolve-Path -LiteralPath (
    Resolve-ProjectPath $DifficultInterfaces)).Path
$auditPath = (Resolve-Path -LiteralPath (
    Resolve-ProjectPath $TopologyAudit)).Path
$rawPath = Resolve-ProjectPath $RawRoot
$aggregatePath = Resolve-ProjectPath $Aggregate

$difficult = @(
    Import-Csv -LiteralPath $difficultPath |
    Sort-Object {[int]$_.difficulty_rank} |
    Select-Object -First 5
)
if ($difficult.Count -ne 5) {
    throw 'Adaptive rank sweep requires exactly five difficult interfaces.'
}

$rramAudit = @(
    Import-Csv -LiteralPath $auditPath |
    Where-Object {$_.case -like 'rram26*'}
)
if ($rramAudit.Count -ne 25) {
    throw 'Adaptive rank sweep requires the verified 25-interface RRAM26 audit.'
}
foreach ($item in $difficult) {
    $match = @($rramAudit | Where-Object {
        [int]$_.interface_id -eq [int]$item.interface_id
    })
    if ($match.Count -ne 1 -or
        [int]$match[0].left_subdomain -ne [int]$item.left_subdomain -or
        [int]$match[0].right_subdomain -ne [int]$item.right_subdomain -or
        [int]$match[0].target_dofs -ne [int]$item.target_dofs -or
        [int]$match[0].source_dofs -ne [int]$item.source_dofs) {
        throw "Difficult-interface topology mismatch for interface " +
            "$($item.interface_id)."
    }
}

$configurations = @(
    [pscustomobject]@{
        Name = 'A'
        HistoryRank = 64
        TransferRank = 16
        ResidualRank = 4
    },
    [pscustomobject]@{
        Name = 'B'
        HistoryRank = 128
        TransferRank = 16
        ResidualRank = 4
    },
    [pscustomobject]@{
        Name = 'C'
        HistoryRank = 128
        TransferRank = 32
        ResidualRank = 8
    },
    [pscustomobject]@{
        Name = 'D'
        HistoryRank = 256
        TransferRank = 32
        ResidualRank = 8
    }
)

if (-not (Test-Path -LiteralPath $rawPath)) {
    New-Item -ItemType Directory -Path $rawPath | Out-Null
} elseif (-not $Resume -and -not $AggregateOnly) {
    $existing = @(Get-ChildItem -LiteralPath $rawPath -Force)
    if ($existing.Count -ne 0) {
        throw "Refusing to overwrite nonempty adaptive rank run root: $rawPath"
    }
}

$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'

$rows = @()
$interfaceIds = @(
    $difficult | ForEach-Object {[int]$_.interface_id}
)
$interfaceIdArgument = $interfaceIds -join ','
Push-Location $project
try {
    foreach ($configuration in $configurations) {
        $runName = "top5_$($configuration.Name)"
        $output = Join-Path $rawPath $runName
        $pilotPath = Join-Path $output `
            'milestone8_adaptive_port_local_pilot.csv'
        $stopPath = Join-Path $output `
            'milestone8_adaptive_port_local_pilot_stop.csv'
        if (-not $AggregateOnly -and
            -not ($Resume -and
                (Test-Path -LiteralPath $pilotPath -PathType Leaf) -and
                (Test-Path -LiteralPath $stopPath -PathType Leaf))) {
            if (Test-Path -LiteralPath $output) {
                throw "Refusing to overwrite incomplete run: $output"
            }
            New-Item -ItemType Directory -Path $output | Out-Null
            Write-Host (
                "[M8.8] interfaces=$interfaceIdArgument configuration=" +
                "$($configuration.Name) history=" +
                "$($configuration.HistoryRank) transfer=" +
                "$($configuration.TransferRank) residual=" +
                "$($configuration.ResidualRank)")
            & $exePath --transient --config $configPath `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method hybrid-randomized `
                --mor-arnoldi-moments 1 `
                --randomized-port-rank $configuration.TransferRank `
                --randomized-oversampling 5 `
                --randomized-power-iterations 1 `
                --randomized-seed 12345 `
                --residual-krylov-max-rank `
                    $configuration.ResidualRank `
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
                --history-compression-rank `
                    $configuration.HistoryRank `
                --history-compression-tolerance 1e-12 `
                --adaptive-port-local-pilot `
                --adaptive-port-interface-ids $interfaceIdArgument `
                --optimal-port-topology-audit-csv $auditPath `
                --mor-transient-dt 0.01 `
                --mor-transient-t-end 0.01 `
                --mor-transient-waveform single_step `
                --mor-transient-initial-mode ambient `
                --output-dir $output --fast-run
            if ($LASTEXITCODE -ne 0) {
                throw "Adaptive pilot process failed: $runName"
            }
        }

        $pilot = @(Import-Csv -LiteralPath $pilotPath)
        $stop = @(Import-Csv -LiteralPath $stopPath)
        if ($pilot.Count -ne $difficult.Count -or
            $stop.Count -ne 1 -or
            [int]$stop[0].interfaces_requested -ne $difficult.Count -or
            [int]$stop[0].transient_advanced -ne 0 -or
            [int]$stop[0].full_field_read -ne 0 -or
            [int]$stop[0].snapshot_used -ne 0) {
            throw "Adaptive pilot scope/provenance failed: $runName"
        }

        foreach ($interface in $difficult) {
            $interfaceId = [int]$interface.interface_id
            $matches = @($pilot | Where-Object {
                [int]$_.interface_id -eq $interfaceId
            })
            if ($matches.Count -ne 1) {
                throw "Adaptive pilot omitted or duplicated interface " +
                    "$interfaceId in $runName."
            }
            $item = $matches[0]
            if ([int]$item.requested_history_rank -ne
                    $configuration.HistoryRank -or
                [int]$item.requested_randomized_rank -ne
                    $configuration.TransferRank -or
                [int]$item.requested_enrichment_rank -ne
                    $configuration.ResidualRank -or
                [int]$item.transient_advanced -ne 0 -or
                [int]$item.full_field_read -ne 0 -or
                [int]$item.snapshot_used -ne 0 -or
                [int]$item.fom_used_for_basis -ne 0) {
                throw "Adaptive pilot rank/provenance failed: $runName"
            }

            $rows += [pscustomobject][ordered]@{
                case = 'RRAM26'
                difficulty_rank = [int]$interface.difficulty_rank
                interface_id = $interfaceId
                left_subdomain = [int]$interface.left_subdomain
                right_subdomain = [int]$interface.right_subdomain
                target_dofs = [int]$item.target_dofs
                source_dofs = [int]$item.source_dofs
                configuration = $configuration.Name
                requested_history_rank =
                    [int]$item.requested_history_rank
                compressed_history_rank =
                    [int]$item.compressed_history_rank
                requested_transfer_rank =
                    [int]$item.requested_randomized_rank
                accepted_transfer_rank =
                    [int]$item.accepted_randomized_rank
                requested_residual_rank =
                    [int]$item.requested_enrichment_rank
                accepted_residual_rank =
                    [int]$item.accepted_enrichment_rank
                mandatory_rank = [int]$item.mandatory_rank
                total_port_rank = [int]$item.total_port_rank
                temperature_error_proxy =
                    [double]$item.history_compression_relative_error
                flux_error_proxy =
                    [double]$item.basis_error_indicator
                interface_residual_proxy =
                    [double]$item.final_max_probe_residual
                initial_interface_residual_proxy =
                    [double]$item.initial_max_probe_residual
                local_basis_time_s =
                    [double]$item.total_basis_build_time_s
                local_incremental_memory_bytes =
                    [uint64]$item.incremental_memory_bytes
                process_peak_memory_bytes =
                    [uint64]$item.process_peak_memory_bytes
                history_rhs = [int]$item.history_target_rhs
                randomized_target_solve_count =
                    [int]$item.randomized_target_solve_count
                residual_target_solve_count =
                    [int]$item.residual_target_solve_count
                total_target_solve_count =
                    [int]$item.total_target_solve_count
                target_residual = [double]$item.target_residual
                weighted_adjoint_error =
                    [double]$item.weighted_adjoint_error
                status = $item.status
                physical_transient_error_run = 0
                transient_advanced = 0
                full_field_read = 0
                snapshot_used = 0
                fom_used_for_basis = 0
                pod_used = 0
                svd_used = 0
            }
        }
    }
} finally {
    Pop-Location
}

$aggregateDirectory = Split-Path -Parent $aggregatePath
if (-not (Test-Path -LiteralPath $aggregateDirectory)) {
    New-Item -ItemType Directory -Path $aggregateDirectory | Out-Null
}
$rows |
    Sort-Object difficulty_rank, configuration |
    Export-Csv -LiteralPath $aggregatePath -NoTypeInformation

$expectedRows = $difficult.Count * $configurations.Count
if ($rows.Count -ne $expectedRows) {
    throw "Expected $expectedRows adaptive rows, found $($rows.Count)."
}
Write-Host "[M8.8] wrote $aggregatePath ($($rows.Count) rows)."
