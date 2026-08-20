param(
    [string]$RramM7Directory =
        'D:\AI agent\kimi\DDM_Schur\results\local_port_milestone7_rram26_final',
    [string]$ChipletM7Directory =
        'D:\AI agent\kimi\DDM_Schur\results\local_port_milestone7_chiplet_final',
    [string]$RramM6Directory =
        'D:\AI agent\kimi\DDM_Schur\results\local_block_arnoldi_milestone6_rram26',
    [string]$ChipletM6Directory =
        'D:\AI agent\kimi\DDM_Schur\results\local_block_arnoldi_milestone6_chiplet_final',
    [string]$RramM6DeploymentDirectory =
        'D:\AI agent\kimi\DDM_Schur\results\local_block_arnoldi_milestone6_rram26_deployment100',
    [string]$ChipletM6DeploymentDirectory =
        'D:\AI agent\kimi\DDM_Schur\results\local_block_arnoldi_milestone6_chiplet_deployment100'
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$outputs = Join-Path $project 'outputs'

function Read-One([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        throw "Missing M8.7 input: $Path"
    }
    return Import-Csv -LiteralPath $Path | Select-Object -First 1
}

function Number($Object, [string]$Name) {
    $property = $Object.PSObject.Properties[$Name]
    if ($null -eq $property -or
        [string]::IsNullOrWhiteSpace([string]$property.Value)) {
        return 0.0
    }
    return [double]$property.Value
}

function Read-LargeBaseline([string]$Directory) {
    return Read-One (Join-Path $Directory 'local_dynamic_schur_summary.csv')
}

function Read-Deployment([string]$Directory) {
    return Read-One (
        Join-Path $Directory 'local_dynamic_schur_deployment_timing.csv')
}

function Make-AccuracyRow(
    [string]$Case, [string]$Method, [string]$Scope, $Summary,
    [int]$SnapshotUsed, [int]$FomUsedForBasis, [string]$Source) {
    return [pscustomobject]@{
        case = $Case
        method = $Method
        scope = $Scope
        port_dimension = Number $Summary 'port_dimension'
        total_local_rank = Number $Summary 'total_local_rank'
        temperature_relative_l2 =
            Number $Summary 'space_time_relative_l2'
        max_error_k = Number $Summary 'maximum_absolute_k'
        max_temperature_error_k =
            Number $Summary 'maximum_temperature_error_k'
        flux_relative_l2 =
            Number $Summary 'maximum_fom_rom_flux_relative_l2'
        interface_residual =
            Number $Summary 'maximum_interface_relative_residual'
        full_residual = Number $Summary 'maximum_full_residual'
        online_time_s = Number $Summary 'local_online_core_seconds'
        total_time_s = Number $Summary 'total_seconds'
        peak_memory_bytes = Number $Summary 'peak_working_set_bytes'
        corrected = 0
        snapshot_used = $SnapshotUsed
        fom_used_for_basis = $FomUsedForBasis
        data_source = $Source
        validation_status =
            $(if ($Method -eq 'M8 Hybrid Operator Port') {
                if (Is-ValidM8Accuracy $Summary) {
                    'success'
                } elseif (
                    (Number $Summary 'space_time_relative_l2') -ge
                        1.0e-4 -or
                    (Number $Summary 'maximum_absolute_k') -ge 0.1) {
                    'temperature_accuracy_failed'
                } else {
                    'flux_interface_accuracy_failed'
                }
            } else {
                'benchmark_reference'
            })
        status = $Summary.status
    }
}

function BreakEven(
    [double]$CandidateOffline, [double]$BaselineOffline,
    [double]$CandidateOnline, [double]$BaselineOnline) {
    $saving = $BaselineOnline - $CandidateOnline
    if ($saving -le 0.0) {
        return ''
    }
    return [Math]::Max(
        0.0, ($CandidateOffline - $BaselineOffline) / $saving)
}

function Is-ValidM8Accuracy($Summary) {
    return (
        (Number $Summary 'space_time_relative_l2') -lt 1.0e-4 -and
        (Number $Summary 'maximum_absolute_k') -lt 0.1 -and
        (Number $Summary 'maximum_fom_rom_flux_relative_l2') -lt 0.05 -and
        (Number $Summary 'maximum_interface_relative_residual') -lt
            1.0e-4 -and
        (Number $Summary 'maximum_full_residual') -lt 0.05)
}

$cases = @(
    [pscustomobject]@{
        name = 'rram26'
        m7 = $RramM7Directory
        m6 = $RramM6Directory
        m6Deployment = $RramM6DeploymentDirectory
    },
    [pscustomobject]@{
        name = 'chiplet'
        m7 = $ChipletM7Directory
        m6 = $ChipletM6Directory
        m6Deployment = $ChipletM6DeploymentDirectory
    })

$accuracyRows = [System.Collections.Generic.List[object]]::new()
$efficiencyRows = [System.Collections.Generic.List[object]]::new()
$memoryRows = [System.Collections.Generic.List[object]]::new()
$breakEvenRows = [System.Collections.Generic.List[object]]::new()
$masterRows = [System.Collections.Generic.List[object]]::new()

foreach ($case in $cases) {
    $name = $case.name
    $basis = Read-One (
        Join-Path $outputs "milestone8_${name}_production_basis.csv")
    $oneStepPath = Join-Path $project (
        "results\milestone8_accuracy_$name\milestone8_accuracy.csv")
    $oneStep = @(Import-Csv -LiteralPath $oneStepPath)
    if ($oneStep.Count -ne 3) {
        throw "Expected three one-step accuracy rows for $name."
    }
    foreach ($row in $oneStep) {
        $accuracyRows.Add([pscustomobject]@{
            case = $row.case
            method = $row.method
            scope = $row.scope
            port_dimension = $row.port_dimension
            total_local_rank = $row.total_local_rank
            temperature_relative_l2 = $row.temperature_relative_l2
            max_error_k = $row.max_error_k
            max_temperature_error_k = $row.max_temperature_error_k
            flux_relative_l2 = $row.flux_relative_l2
            interface_residual = $row.interface_residual
            full_residual = $row.full_residual
            online_time_s = $row.online_time_s
            total_time_s = $row.total_time_s
            peak_memory_bytes = $row.peak_memory_bytes
            corrected = $row.corrected
            snapshot_used = $row.snapshot_used
            fom_used_for_basis = $row.fom_used_for_basis
            data_source = 'm8.7_fresh'
            validation_status =
                $(if ($row.method -eq 'M8 Hybrid Operator Port') {
                    if ([double]$row.temperature_relative_l2 -ge 1.0e-4 -or
                        [double]$row.max_error_k -ge 0.1) {
                        'temperature_accuracy_failed'
                    } elseif (
                        [double]$row.flux_relative_l2 -ge 0.05 -or
                        [double]$row.interface_residual -ge 1.0e-4 -or
                        [double]$row.full_residual -ge 0.05) {
                        'flux_interface_accuracy_failed'
                    } else {
                        'success'
                    }
                } else {
                    'benchmark_reference'
                })
            status = $row.status
        })
    }

    $m8Directory = Join-Path $project "results\milestone8_efficiency_$name"
    $m8Summary = Read-One (
        Join-Path $m8Directory 'local_dynamic_schur_summary.csv')
    $m8Deployment = Read-Deployment $m8Directory
    $m7Summary = Read-LargeBaseline $case.m7
    $m7Deployment = Read-Deployment $case.m7
    $m6Summary = Read-LargeBaseline $case.m6
    $m6Deployment = Read-Deployment $case.m6Deployment
    $m8AccuracyValid = Is-ValidM8Accuracy $m8Summary

    if ([int]$basis.snapshot_used -ne 0 -or
        [int]$basis.fom_used_for_basis -ne 0 -or
        [int]$m8Summary.port_snapshot_used -ne 0 -or
        [int]$m8Summary.port_fom_used_for_basis -ne 0) {
        throw "M8 provenance gate failed while aggregating $name."
    }

    $m6Accuracy = Make-AccuracyRow $name 'Full Interface' `
        'full_horizon' $m6Summary 0 0 'validated_m6_baseline'
    $m7Accuracy = Make-AccuracyRow $name 'M7 Port POD' `
        'full_horizon' $m7Summary 1 1 'validated_m7_baseline'
    $m8Accuracy = Make-AccuracyRow $name 'M8 Hybrid Operator Port' `
        'full_horizon' $m8Summary 0 0 'm8.7_fresh'
    $accuracyRows.Add($m6Accuracy)
    $accuracyRows.Add($m7Accuracy)
    $accuracyRows.Add($m8Accuracy)

    $methods = @(
        [pscustomobject]@{
            method = 'Full Interface'
            summary = $m6Summary
            deployment = $m6Deployment
            source = 'validated_m6_baseline'
            snapshot = 0
            fom = 0
        },
        [pscustomobject]@{
            method = 'M7 Port POD'
            summary = $m7Summary
            deployment = $m7Deployment
            source = 'validated_m7_baseline'
            snapshot = 1
            fom = 1
        },
        [pscustomobject]@{
            method = 'M8 Hybrid Operator Port'
            summary = $m8Summary
            deployment = $m8Deployment
            source = 'm8.7_fresh'
            snapshot = 0
            fom = 0
        })

    foreach ($method in $methods) {
        $summary = $method.summary
        $deployment = $method.deployment
        $efficiencyRows.Add([pscustomobject]@{
            case = $name
            method = $method.method
            full_horizon_steps = Number $summary 'steps'
            full_horizon_online_time_s =
                Number $summary 'local_online_core_seconds'
            full_horizon_total_time_s = Number $summary 'total_seconds'
            rhs_count = Number $deployment 'waveforms'
            rhs_total_online_time_s =
                Number $deployment 'total_online_seconds'
            average_online_time_per_rhs_s =
                Number $deployment 'average_online_seconds_per_waveform'
            setup_reused = Number $deployment 'setup_reused'
            port_factor_reused = Number $deployment 'port_factor_reused'
            local_factors_reused =
                Number $deployment 'local_factors_reused'
            snapshot_used = $method.snapshot
            fom_used_for_basis = $method.fom
            data_source = $method.source
            accuracy_valid =
                $(if ($method.method -eq 'M8 Hybrid Operator Port') {
                    if ($m8AccuracyValid) { 1 } else { 0 }
                } else { '' })
            status = $summary.status
        })
        $memoryRows.Add([pscustomobject]@{
            case = $name
            method = $method.method
            process_peak_memory_bytes =
                Number $summary 'peak_working_set_bytes'
            factor_memory_bytes = Number $summary 'factor_memory_bytes'
            dynamic_schur_factor_memory_bytes =
                Number $summary 'dynamic_schur_factor_memory_bytes'
            serialized_model_bytes = Number $summary 'model_bytes'
            basis_peak_incremental_memory_bytes =
                $(if ($method.method -eq 'M8 Hybrid Operator Port') {
                    Number $basis 'peak_incremental_memory_bytes'
                } else { 0 })
            basis_build_process_peak_memory_bytes =
                $(if ($method.method -eq 'M8 Hybrid Operator Port') {
                    Number $basis 'process_peak_memory_bytes'
                } else { 0 })
            data_source = $method.source
        })
    }

    $m8Offline = Number $basis 'total_basis_build_time_s'
    $m8Offline += Number $m8Summary 'dynamic_schur_setup_seconds'
    $m7Offline = Number $m7Summary 'local_basis_setup_seconds'
    $m7Offline += Number $m7Summary 'dynamic_schur_setup_seconds'
    $m7Offline += Number $m7Summary 'port_local_full_interface_pilot_seconds'
    $m7Offline += Number $m7Summary 'enrichment_total_seconds'
    $m6Offline = Number $m6Summary 'local_basis_setup_seconds'
    $m6Offline += Number $m6Summary 'dynamic_schur_setup_seconds'
    $m8Online = Number $m8Summary 'local_online_core_seconds'
    $m7Online = Number $m7Summary 'local_online_core_seconds'
    $m6Online = Number $m6Summary 'local_online_core_seconds'
    $m8Rhs = Number $m8Deployment 'average_online_seconds_per_waveform'
    $m7Rhs = Number $m7Deployment 'average_online_seconds_per_waveform'
    $m6Rhs = Number $m6Deployment 'average_online_seconds_per_waveform'
    foreach ($comparison in @(
        [pscustomobject]@{
            baseline = 'Full Interface'
            offline = $m6Offline
            online = $m6Online
            rhs = $m6Rhs
        },
        [pscustomobject]@{
            baseline = 'M7 Port POD'
            offline = $m7Offline
            online = $m7Online
            rhs = $m7Rhs
        })) {
        $breakEvenRows.Add([pscustomobject]@{
            case = $name
            candidate = 'M8 Hybrid Operator Port'
            baseline = $comparison.baseline
            candidate_offline_time_s = $m8Offline
            baseline_offline_time_s = $comparison.offline
            candidate_full_horizon_online_time_s = $m8Online
            baseline_full_horizon_online_time_s = $comparison.online
            full_horizon_break_even_queries = BreakEven $m8Offline `
                $comparison.offline $m8Online $comparison.online
            candidate_online_time_per_rhs_s = $m8Rhs
            baseline_online_time_per_rhs_s = $comparison.rhs
            one_step_break_even_rhs = BreakEven $m8Offline `
                $comparison.offline $m8Rhs $comparison.rhs
            formula =
                'max(0,(candidate_offline-baseline_offline)/(baseline_online-candidate_online))'
            candidate_accuracy_valid =
                $(if ($m8AccuracyValid) { 1 } else { 0 })
            decision =
                $(if ($m8AccuracyValid) {
                    'deployable_break_even'
                } else {
                    'informational_only_not_deployable'
                })
        })
    }

    $masterRows.Add([pscustomobject]@{
        case = $name
        method = 'M8 Hybrid Operator Port'
        configuration =
            'mandatory+history64+randomized16(q1,p5,seed12345)+residual4(sweeps2,tol1e-4)'
        physical_interfaces = $basis.physical_interfaces
        total_port_rank = $basis.total_port_rank
        total_basis_build_time_s = $basis.total_basis_build_time_s
        max_interface_time_s = $basis.max_interface_time_s
        basis_peak_incremental_memory_bytes =
            $basis.peak_incremental_memory_bytes
        basis_process_peak_memory_bytes =
            $basis.process_peak_memory_bytes
        target_solve_count = $basis.target_solve_count
        temperature_relative_l2 =
            $m8Accuracy.temperature_relative_l2
        max_error_k = $m8Accuracy.max_error_k
        flux_relative_l2 = $m8Accuracy.flux_relative_l2
        interface_residual = $m8Accuracy.interface_residual
        full_residual = $m8Accuracy.full_residual
        full_horizon_online_time_s = $m8Online
        average_online_time_per_rhs_s = $m8Rhs
        serialized_model_bytes = $basis.serialized_model_bytes
        snapshot_used = 0
        fom_used_for_basis = 0
        corrected = 0
        accuracy_valid = $(if ($m8AccuracyValid) { 1 } else { 0 })
        status = $m8Summary.status
        validation_status =
            $(if ($m8AccuracyValid) {
                'success'
            } elseif (
                (Number $m8Summary 'space_time_relative_l2') -ge 1.0e-4 -or
                (Number $m8Summary 'maximum_absolute_k') -ge 0.1) {
                'temperature_accuracy_failed'
            } else {
                'flux_interface_accuracy_failed'
            })
    })
}

$accuracyRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $outputs 'milestone8_accuracy.csv')
$efficiencyRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $outputs 'milestone8_efficiency.csv')
$memoryRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $outputs 'milestone8_memory.csv')
$breakEvenRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $outputs 'milestone8_break_even.csv')
$masterRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $outputs 'milestone8_final_master_summary.csv')
