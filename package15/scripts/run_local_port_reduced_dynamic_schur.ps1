param(
    [string]$RramConfig = 'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$ChipletConfig = 'D:\CPP\TEST_CHATGPT\chiplet_model\case_chiplet_config_horizontal.txt',
    [string]$ResultsDirectory = 'results',
    [string]$OutputsDirectory = 'outputs',
    [switch]$SkipBuild,
    [switch]$SkipRuns,
    [switch]$SkipTests
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$results = [System.IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
$outputs = [System.IO.Path]::GetFullPath((Join-Path $project $OutputsDirectory))
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$rramDirectory = Join-Path $results 'local_port_milestone7_rram26_final'
$chipletDirectory = Join-Path $results 'local_port_milestone7_chiplet_final'

function Assert-Less([double]$Value, [double]$Limit, [string]$Message) {
    if (-not ($Value -lt $Limit)) { throw "$Message (value=$Value, limit=$Limit)" }
}

function Read-One([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { throw "Missing result: $Path" }
    return Import-Csv $Path | Select-Object -First 1
}

function Assert-PortPartition([string]$Directory, $Summary) {
    $ports = @(Import-Csv (Join-Path $Directory 'local_port_rank_by_interface.csv'))
    $covered = ($ports | Measure-Object full_interface_rows -Sum).Sum
    if ($ports.Count -le 0 -or [int]$covered -ne [int]$Summary.full_interface_dofs) {
        throw "Physical ports are not a disjoint, complete SIPG trace partition " +
            "(covered=$covered, interface=$($Summary.full_interface_dofs))."
    }
}

function Invoke-PortCase(
    [string]$Config, [string]$Output, [int]$Moments, [int]$PortRank,
    [int]$EnrichmentRounds, [double]$Dt, [double]$EndTime,
    [string]$Waveform, [double]$TemperatureWeight) {
    New-Item -ItemType Directory -Force -Path $Output | Out-Null
    & $exe --transient --config $Config --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --mor-arnoldi-moments $Moments --mor-interface-rank 0 `
        --local-port-rank $PortRank --local-port-energy-tolerance 1e-20 `
        --local-port-temperature-weight $TemperatureWeight `
        --local-port-flux-weight 1 --local-port-residual-weight 1 `
        --local-port-enrichment-rounds $EnrichmentRounds --local-port-corrected `
        --mor-transient-dt $Dt --mor-transient-t-end $EndTime `
        --mor-transient-waveform $Waveform --mor-transient-initial-mode ambient `
        --mor-deployment-rhs-count 100 --output-dir $Output --fast-run
    if ($LASTEXITCODE -ne 0) { throw "Local-port run failed: $Output" }
}

Push-Location $project
try {
    New-Item -ItemType Directory -Force -Path $results,$outputs | Out-Null
    if (-not $SkipBuild) {
        cmake -S . -B build -DBUILD_TESTING=ON
        if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed.' }
        cmake --build build --config Release --parallel
        if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
    }
    if (-not (Test-Path -LiteralPath $exe)) { throw "Release executable missing: $exe" }

    if (-not $SkipRuns) {
        Invoke-PortCase $RramConfig $rramDirectory 1 128 5 0.01 0.8 `
            'mixed_frequency' 100
        Invoke-PortCase $ChipletConfig $chipletDirectory 4 200 5 0.1 1.0 `
            'asynchronous_hotspots' 1
    }

    if (-not $SkipTests) {
        $ctestOutput = @(& ctest --test-dir build -C Release --output-on-failure 2>&1)
        $ctestOutput | Set-Content -Encoding UTF8 `
            (Join-Path $outputs 'local_port_ctest_results.txt')
        if ($LASTEXITCODE -ne 0) { throw 'Release CTest failed.' }
    }

    $rram = Read-One (Join-Path $rramDirectory 'local_dynamic_schur_summary.csv')
    $chiplet = Read-One (Join-Path $chipletDirectory 'local_dynamic_schur_summary.csv')
    $rramDeployment = Read-One (Join-Path $rramDirectory `
        'local_dynamic_schur_deployment_timing.csv')
    $chipletDeployment = Read-One (Join-Path $chipletDirectory `
        'local_dynamic_schur_deployment_timing.csv')

    Assert-PortPartition $rramDirectory $rram
    Assert-PortPartition $chipletDirectory $chiplet

    if ([int]$rram.port_dimension -ge [int]$rram.full_interface_dofs) {
        throw 'RRAM port interface was not reduced.'
    }
    Assert-Less ([double]$rram.space_time_relative_l2) 1e-4 `
        'RRAM pure-ROM relative L2 gate failed'
    Assert-Less ([double]$rram.maximum_absolute_k) 0.1 `
        'RRAM pure-ROM maximum error gate failed'
    Assert-Less ([double]$rramDeployment.average_online_seconds_per_waveform) 5.0 `
        'RRAM first online-time target failed'
    Assert-Less ([double]$chiplet.space_time_relative_l2) 1e-4 `
        'Chiplet pure-ROM relative L2 gate failed'
    Assert-Less ([double]$chiplet.maximum_absolute_k) 0.1 `
        'Chiplet pure-ROM maximum error gate failed'
    Assert-Less ([double]$chiplet.maximum_fom_rom_flux_relative_l2) 0.05 `
        'Chiplet flux gate failed'
    Assert-Less ([double]$chiplet.maximum_full_residual) 0.05 `
        'Chiplet full-residual gate failed'
    Assert-Less ([double]$chipletDeployment.average_online_seconds_per_waveform) 4.47 `
        'Chiplet online-time gate failed'

    $cases = @(
        [pscustomobject]@{Name='rram26';Summary=$rram;Deployment=$rramDeployment;Directory=$rramDirectory},
        [pscustomobject]@{Name='chiplet';Summary=$chiplet;Deployment=$chipletDeployment;Directory=$chipletDirectory})

    $summaryRows = foreach ($case in $cases) {
        $s = $case.Summary; $d = $case.Deployment
        [pscustomobject]@{
            case_name=$case.Name; status='success'; method=$s.method
            subdomains=$s.subdomains; global_dofs=$s.global_dofs
            full_interface_dofs=$s.full_interface_dofs; port_dimension=$s.port_dimension
            interface_compression=([double]$s.full_interface_dofs/[double]$s.port_dimension)
            total_local_rank=$s.total_local_rank; enrichment_added_rank=$s.enrichment_added_rank
            space_time_relative_l2=$s.space_time_relative_l2
            maximum_absolute_k=$s.maximum_absolute_k
            maximum_flux_relative_l2=$s.maximum_fom_rom_flux_relative_l2
            maximum_full_residual=$s.maximum_full_residual
            maximum_port_projection_relative_error=$s.maximum_port_projection_relative_error
            maximum_port_reduced_relative_residual=$s.maximum_port_reduced_relative_residual
            port_schur_relative_asymmetry=$s.port_schur_relative_asymmetry
            corrected_relative_l2=$s.corrected_relative_l2
            corrected_maximum_absolute_k=$s.corrected_maximum_absolute_k
            corrected_maximum_full_residual=$s.corrected_maximum_full_residual
            average_online_seconds_100_rhs=$d.average_online_seconds_per_waveform
            total_seconds=$s.total_seconds
        }
    }
    $summaryRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_port_rom_summary.csv')

    $rankOutput = Join-Path $outputs 'local_port_rank_by_interface.csv'
    if (Test-Path -LiteralPath $rankOutput) {
        Remove-Item -LiteralPath $rankOutput -Force
    }
    foreach ($case in $cases) {
        Import-Csv (Join-Path $case.Directory 'local_port_rank_by_interface.csv') |
            Select-Object @{Name='case_name';Expression={$case.Name}},* |
            Export-Csv -NoTypeInformation -Encoding UTF8 -Append `
                $rankOutput
    }

    $temperatureRows = foreach ($case in $cases) {
        $s = $case.Summary
        [pscustomobject]@{case_name=$case.Name;mode='pure';relative_l2=$s.space_time_relative_l2;maximum_absolute_k=$s.maximum_absolute_k;maximum_temperature_error_k=$s.maximum_temperature_error_k;status=$s.status}
        [pscustomobject]@{case_name=$case.Name;mode='corrected';relative_l2=$s.corrected_relative_l2;maximum_absolute_k=$s.corrected_maximum_absolute_k;maximum_temperature_error_k='';status=$(if ([double]$s.corrected_maximum_full_residual -lt 1e-8) {'success'} else {'accuracy_failed'})}
    }
    $temperatureRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_port_temperature_accuracy.csv')

    $fluxRows = foreach ($case in $cases) {
        $s = $case.Summary
        [pscustomobject]@{case_name=$case.Name;maximum_temperature_jump_rms_k=$s.maximum_temperature_jump_rms_k;maximum_relative_flux_imbalance=$s.maximum_relative_flux_imbalance;fom_rom_flux_relative_l2=$s.maximum_fom_rom_flux_relative_l2}
    }
    $fluxRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_port_flux_accuracy.csv')

    $residualRows = foreach ($case in $cases) {
        $s = $case.Summary
        [pscustomobject]@{case_name=$case.Name;full_residual=$s.maximum_full_residual;reduced_residual=$s.maximum_reduced_residual;port_relative_residual=$s.maximum_interface_relative_residual;port_projected_system_relative_residual=$s.maximum_port_reduced_relative_residual;port_projection_relative_error=$s.maximum_port_projection_relative_error;port_schur_relative_asymmetry=$s.port_schur_relative_asymmetry;corrected_full_residual=$s.corrected_maximum_full_residual}
    }
    $residualRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_port_residual_accuracy.csv')

    $timingRows = foreach ($case in $cases) {
        $s = $case.Summary; $d = $case.Deployment
        [pscustomobject]@{case_name=$case.Name;steps=$s.steps;reference_setup_seconds=$s.reference_setup_seconds;local_basis_setup_seconds=$s.local_basis_setup_seconds;enrichment_seconds=$s.enrichment_total_seconds;port_local_full_interface_pilot_seconds=$s.port_local_full_interface_pilot_seconds;port_snapshot_seconds=$s.port_snapshot_seconds;port_basis_seconds=$s.port_basis_seconds;port_schur_assembly_seconds=$s.port_schur_assembly_seconds;port_schur_factor_seconds=$s.port_schur_factor_seconds;interface_solve_seconds=$s.interface_solve_seconds;local_recovery_seconds=$s.local_recovery_seconds;online_core_seconds=$s.local_online_core_seconds;average_online_seconds_100_rhs=$d.average_online_seconds_per_waveform;total_seconds=$s.total_seconds}
    }
    $timingRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_port_dynamic_schur_timing.csv')

    $memoryRows = foreach ($case in $cases) {
        $s = $case.Summary
        [pscustomobject]@{case_name=$case.Name;model_bytes=$s.model_bytes;port_dynamic_factor_bytes=$s.dynamic_schur_factor_memory_bytes;fom_factor_bytes=$s.fom_factor_memory_bytes;combined_factor_bytes=$s.factor_memory_bytes;peak_working_set_bytes=$s.peak_working_set_bytes}
    }
    $memoryRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_port_memory.csv')

    $enrichmentRows = foreach ($case in $cases) {
        Import-Csv (Join-Path $case.Directory 'local_port_enrichment_history.csv') |
            Select-Object @{Name='case_name';Expression={$case.Name}},*
    }
    $enrichmentRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_port_enrichment_history.csv')

    $milestone6 = Import-Csv (Join-Path $outputs `
        'local_block_arnoldi_milestone6_summary.csv')
    $fullComparison = foreach ($case in $cases) {
        $old = $milestone6 | Where-Object case_name -eq $case.Name | Select-Object -First 1
        $s = $case.Summary; $d = $case.Deployment
        $portOffline = [double]$s.reference_setup_seconds +
            [double]$s.local_basis_setup_seconds +
            [double]$s.enrichment_total_seconds +
            [double]$s.port_local_full_interface_pilot_seconds +
            [double]$s.port_snapshot_seconds + [double]$s.port_basis_seconds +
            [double]$s.dynamic_schur_setup_seconds
        $fullOffline = [double]$old.reference_setup_seconds +
            [double]$old.local_basis_setup_seconds +
            [double]$old.dynamic_schur_setup_seconds
        $fullOnline = if ($case.Name -eq 'rram26') {32.333254485} else {4.466560076}
        $breakEven = if ($fullOnline -gt [double]$d.average_online_seconds_per_waveform) {
            [Math]::Ceiling([Math]::Max(0.0,
                ($portOffline - $fullOffline) /
                ($fullOnline - [double]$d.average_online_seconds_per_waveform)))
        } else { [double]::PositiveInfinity }
        [pscustomobject]@{case_name=$case.Name;method='Milestone 6 full-interface';interface_dimension=$old.full_interface_dofs;local_rank=$old.total_local_rank;offline_seconds=([double]$old.reference_setup_seconds+[double]$old.local_basis_setup_seconds+[double]$old.dynamic_schur_setup_seconds);online_seconds_per_one_step_rhs=$(if ($case.Name -eq 'rram26') {32.333254485} else {4.466560076});relative_l2=$old.space_time_relative_l2;maximum_absolute_k=$old.maximum_absolute_k;flux_relative_l2=$old.maximum_fom_rom_flux_relative_l2;full_residual=$old.maximum_full_residual;model_bytes=$old.model_bytes;peak_working_set_bytes=$old.peak_working_set_bytes;break_even_vs_milestone6_waveforms='baseline'}
        [pscustomobject]@{case_name=$case.Name;method='Milestone 7 port-reduced';interface_dimension=$s.port_dimension;local_rank=$s.total_local_rank;offline_seconds=$portOffline;online_seconds_per_one_step_rhs=$d.average_online_seconds_per_waveform;relative_l2=$s.space_time_relative_l2;maximum_absolute_k=$s.maximum_absolute_k;flux_relative_l2=$s.maximum_fom_rom_flux_relative_l2;full_residual=$s.maximum_full_residual;model_bytes=$s.model_bytes;peak_working_set_bytes=$s.peak_working_set_bytes;break_even_vs_milestone6_waveforms=$breakEven}
        [pscustomobject]@{case_name=$case.Name;method='Milestone 7 corrected';interface_dimension=$s.port_dimension;local_rank=$s.total_local_rank;offline_seconds=$portOffline;online_seconds_per_one_step_rhs=([double]$d.average_online_seconds_per_waveform+[double]$s.corrected_solve_seconds/[double]$s.steps);relative_l2=$s.corrected_relative_l2;maximum_absolute_k=$s.corrected_maximum_absolute_k;flux_relative_l2='';full_residual=$s.corrected_maximum_full_residual;model_bytes=$s.model_bytes;peak_working_set_bytes=$s.peak_working_set_bytes;break_even_vs_milestone6_waveforms='diagnostic-only'}
    }
    $fullComparison | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_port_vs_full_interface.csv')

    $global = Import-Csv (Join-Path $outputs `
        'local_transient_vs_global_block_arnoldi.csv')
    $globalRows = foreach ($case in $cases) {
        $g = $global | Where-Object case_name -eq $case.Name
        $s = $case.Summary; $d = $case.Deployment
        $g
        $portOffline = [double]$s.reference_setup_seconds +
            [double]$s.local_basis_setup_seconds +
            [double]$s.enrichment_total_seconds +
            [double]$s.port_local_full_interface_pilot_seconds +
            [double]$s.port_snapshot_seconds + [double]$s.port_basis_seconds +
            [double]$s.dynamic_schur_setup_seconds
        [pscustomobject]@{case_name=$case.Name;method='Local Block Arnoldi + port Dynamic Schur';basis_scope=$(if ($case.Name -eq 'rram26') {'26 local bases + 25 independent ports'} else {'2 local bases + 1 independent port'});rank=$s.total_local_rank;interface_dofs=$s.port_dimension;offline_seconds=$portOffline;online_seconds=$s.local_online_core_seconds;iterations=0;relative_l2=$s.space_time_relative_l2;maximum_absolute_k=$s.maximum_absolute_k;model_bytes=$s.model_bytes;peak_working_set_bytes=$s.peak_working_set_bytes;status='success'}
        [pscustomobject]@{case_name=$case.Name;method='Local port Dynamic Schur + corrected';basis_scope='same local/port bases + one FOM residual correction';rank=$s.total_local_rank;interface_dofs=$s.port_dimension;offline_seconds=$portOffline;online_seconds=([double]$s.local_online_core_seconds+[double]$s.corrected_solve_seconds);iterations=0;relative_l2=$s.corrected_relative_l2;maximum_absolute_k=$s.corrected_maximum_absolute_k;model_bytes=$s.model_bytes;peak_working_set_bytes=$s.peak_working_set_bytes;status=$(if ([double]$s.corrected_maximum_full_residual -lt 1e-8) {'success'} else {'accuracy_failed'})}
    }
    $globalRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_port_vs_global_block_arnoldi.csv')

    $rramSpeedup = 32.333254485 / [double]$rramDeployment.average_online_seconds_per_waveform
    $chipletSpeedup = 4.466560076 / [double]$chipletDeployment.average_online_seconds_per_waveform
    $rramBreakEven = ($fullComparison | Where-Object {
        $_.case_name -eq 'rram26' -and $_.method -eq 'Milestone 7 port-reduced'
    }).break_even_vs_milestone6_waveforms
    $chipletBreakEven = ($fullComparison | Where-Object {
        $_.case_name -eq 'chiplet' -and $_.method -eq 'Milestone 7 port-reduced'
    }).break_even_vs_milestone6_waveforms
    $report = @"
# Milestone 7: Local port-reduced Dynamic Schur

## Decision

Milestone 7 passes the pure-ROM temperature gates on RRAM26 and the temperature, flux, and sampled full-residual gates on Chiplet. Each physical SIPG interface owns an independent, serialized port basis; no global interface POD or sliced global basis is used. Fixed-dt local reduced factors and the projected port Schur factor are constructed once and reused for all time steps and the measured 100-RHS deployment.

| Case | Subdomains | Full interface | Port dimension | Local rank | Relative L2 | Max error (K) | Flux L2 | Full residual | 100-RHS avg (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| RRAM26 | $($rram.subdomains) | $($rram.full_interface_dofs) | $($rram.port_dimension) | $($rram.total_local_rank) | $($rram.space_time_relative_l2) | $($rram.maximum_absolute_k) | $($rram.maximum_fom_rom_flux_relative_l2) | $($rram.maximum_full_residual) | $($rramDeployment.average_online_seconds_per_waveform) |
| Chiplet | $($chiplet.subdomains) | $($chiplet.full_interface_dofs) | $($chiplet.port_dimension) | $($chiplet.total_local_rank) | $($chiplet.space_time_relative_l2) | $($chiplet.maximum_absolute_k) | $($chiplet.maximum_fom_rom_flux_relative_l2) | $($chiplet.maximum_full_residual) | $($chipletDeployment.average_online_seconds_per_waveform) |

RRAM26 one-step online time is $([Math]::Round($rramSpeedup,2))x faster than the Milestone 6 full-interface result; Chiplet is $([Math]::Round($chipletSpeedup,2))x faster. Including all port snapshot, pilot, enrichment, and projected-Schur setup costs, the break-even against Milestone 6 is $rramBreakEven / $chipletBreakEven waveforms for RRAM26 / Chiplet. These are actual setup-once batches of 100 distinct one-step RHS with full-field recovery, not projections from a single run.

## Construction and enrichment

Port snapshots are kept separate per physical interface and scale-balanced across temperature/increment, weak SIPG flux, and discrete residual families. Each port uses a snapshot-Gram POD after projecting the geometry modes, and mandatory target-waveform traces from both the FOM and the Milestone 6 full-interface Local-ROM pilot are retained before energy/fixed-rank truncation. The two nonmatching sides share one physical port coordinate through the existing SIPG interface ordering. Junction trace DOFs receive a deterministic unique physical-port owner, so port supports are disjoint and their row count exactly equals the full SIPG interface dimension. Exact fingerprint equality gates template reuse.

Physical-normal flux on both sides and the SIPG numerical flux continue to be evaluated with the original overlap-triangle integration for every reported time point; their detailed FOM/ROM values remain in each case result directory. No matching-node assumption or replacement interface quadrature is introduced.

Residual-driven enrichment factors each local full interior step operator once, solves the strongest local residual corrections, applies two-pass MGS, and reprojects the local C/K/input/boundary blocks. Multiple rational expansion points were not activated because the stated pure-ROM targets were met without them.

Corrected mode reaches relative L2 $($rram.corrected_relative_l2) / $($chiplet.corrected_relative_l2) and full residual $($rram.corrected_maximum_full_residual) / $($chiplet.corrected_maximum_full_residual) on RRAM26 / Chiplet. Corrected metrics are reported separately and are not used to pass the pure-ROM gates.

The maximum projected-system relative residual is $($rram.maximum_port_reduced_relative_residual) / $($chiplet.maximum_port_reduced_relative_residual), the maximum port projection error is $($rram.maximum_port_projection_relative_error) / $($chiplet.maximum_port_projection_relative_error), and the projected Schur relative asymmetry is $($rram.port_schur_relative_asymmetry) / $($chiplet.port_schur_relative_asymmetry). The full-interface pilot costs $($rram.port_local_full_interface_pilot_seconds) s / $($chiplet.port_local_full_interface_pilot_seconds) s offline and is never executed during deployment.

## Baselines and reproducibility

Milestone 6 full-interface and Global Block Arnoldi remain benchmark-only paths. FEM/SIPG assembly, triangle-overlap/BVH, Stage 1 exact Schur, the FGMRES residual gate, temperature recovery, and the Milestone 1-6 algorithms were not changed. Detailed ranks, timing, memory, temperature, flux, residual, enrichment, and comparison data are in the companion CSV files. Reproduce with `scripts/run_local_port_reduced_dynamic_schur.ps1`.
"@
    $report | Set-Content -Encoding UTF8 (Join-Path $outputs 'local_port_rom_report.md')
} finally {
    Pop-Location
}
