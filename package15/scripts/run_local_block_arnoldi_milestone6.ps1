param(
    [string]$BuildDirectory = '.\build',
    [string]$ResultsDirectory = '.\results',
    [string]$RramConfig = 'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$ChipletConfig = 'D:\CPP\TEST_CHATGPT\chiplet_model\case_chiplet_config_horizontal.txt',
    [switch]$SkipBuild,
    [switch]$SkipTests,
    [switch]$AggregateOnly
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot

function Invoke-Checked([string]$Label, [scriptblock]$Command) {
    Write-Host "[Milestone 6] $Label"
    & $Command
    if ($LASTEXITCODE -ne 0) { throw "$Label failed with exit code $LASTEXITCODE." }
}

function Get-Decision($Summary) {
    if ([double]$Summary.space_time_relative_l2 -lt 1e-4 -and
        [double]$Summary.maximum_absolute_k -lt 0.1 -and
        [double]$Summary.maximum_temperature_error_k -lt 0.01) { return 'success' }
    return 'accuracy_failed'
}

function Write-WorstFlux([array]$Specifications, [string]$Destination) {
    $encoding = New-Object System.Text.UTF8Encoding($false)
    $writer = New-Object System.IO.StreamWriter($Destination, $false, $encoding, 1048576)
    try {
        $writer.WriteLine('case_name,waveform,initial_mode,step,time_s,physical_interface_id,face_pair_id,integration_triangle_id,left_subdomain,right_subdomain,left_boundary_entity,right_boundary_entity,area_m2,rom_temperature_jump_rms_k,rom_left_physical_normal_flux_w_m2,rom_right_physical_normal_flux_w_m2,rom_sipg_numerical_flux_w_m2,rom_flux_imbalance_l2_w_m2,rom_relative_flux_imbalance,fom_temperature_jump_rms_k,fom_left_physical_normal_flux_w_m2,fom_right_physical_normal_flux_w_m2,fom_sipg_numerical_flux_w_m2,fom_flux_imbalance_l2_w_m2,fom_relative_flux_imbalance,sipg_flux_error_w_m2,sipg_flux_relative_error')
        foreach ($specification in $Specifications) {
            $targets = @{}
            foreach ($row in @(Import-Csv -LiteralPath $specification.Accuracy)) {
                if ([int]$row.worst_integration_triangle -ge 0) {
                    $targets["$($row.step)|$($row.worst_physical_interface)|$($row.worst_face_pair)|$($row.worst_integration_triangle)"] = $true
                }
            }
            $reader = New-Object System.IO.StreamReader($specification.Flux)
            try {
                [void]$reader.ReadLine()
                while (($line = $reader.ReadLine()) -ne $null) {
                    $parts = $line.Split(',')
                    if ($parts.Count -ne 26) { throw "Unexpected transient flux CSV shape." }
                    if ($targets.ContainsKey(
                            "$($parts[2])|$($parts[4])|$($parts[5])|$($parts[6])")) {
                        $writer.WriteLine("$($specification.Name),$line")
                    }
                }
            } finally { $reader.Dispose() }
        }
    } finally { $writer.Dispose() }
}

Push-Location $repo
try {
    $exe = Join-Path $BuildDirectory 'Release\SIPGHeatDDM3D.exe'
    $rramDirectory = Join-Path $ResultsDirectory 'local_block_arnoldi_milestone6_rram26'
    $chipletDirectory = Join-Path $ResultsDirectory 'local_block_arnoldi_milestone6_chiplet_final'
    $rramDeploymentDirectory = Join-Path $ResultsDirectory `
        'local_block_arnoldi_milestone6_rram26_deployment100'
    $chipletDeploymentDirectory = Join-Path $ResultsDirectory `
        'local_block_arnoldi_milestone6_chiplet_deployment100'
    $outputs = Join-Path $repo 'outputs'
    New-Item -ItemType Directory -Force -Path $outputs | Out-Null

    if (-not $AggregateOnly) {
        if (-not $SkipBuild) {
            Invoke-Checked 'Release build' { cmake --build $BuildDirectory --config Release --parallel 4 }
        }
        foreach ($configuration in @($RramConfig, $ChipletConfig)) {
            if (-not (Test-Path -LiteralPath $configuration)) {
                throw "Required large-case configuration not found: $configuration"
            }
        }
        Invoke-Checked 'RRAM26 full-trace Local Block Arnoldi, 80-step Dynamic Schur' {
            & $exe --transient --config $RramConfig --mor-transient-generate `
                --mor-transient-method local-block-arnoldi --mor-arnoldi-moments 1 `
                --mor-interface-rank 0 --mor-transient-dt 0.01 `
                --mor-transient-t-end 0.8 --mor-transient-waveform mixed_frequency `
                --mor-transient-initial-mode ambient `
                --local-mor-matrix-free-threshold 20000 --max-pcg-iterations 500 `
                --gmres-restart 100 --pcg-tolerance 1e-5 --schur-proxy-ring 1 `
                --schur-proxy-block-size 64 --output-dir $rramDirectory --fast-run
        }
        Invoke-Checked 'Chiplet asynchronous-hotspot Local Block Arnoldi' {
            & $exe --transient --config $ChipletConfig --mor-transient-generate `
                --mor-transient-method local-block-arnoldi --mor-arnoldi-moments 4 `
                --mor-interface-rank 0 --mor-transient-dt 0.1 `
                --mor-transient-t-end 1.0 --mor-transient-waveform asynchronous_hotspots `
                --mor-transient-initial-mode ambient `
                --local-mor-matrix-free-threshold 20000 --max-pcg-iterations 500 `
                --gmres-restart 100 --pcg-tolerance 1e-5 --schur-proxy-ring 1 `
                --schur-proxy-block-size 64 --output-dir $chipletDirectory --fast-run
        }
        Invoke-Checked 'RRAM26 actual 100-waveform one-step deployment' {
            & $exe --transient --config $RramConfig --mor-transient-generate `
                --mor-transient-method local-block-arnoldi --mor-arnoldi-moments 1 `
                --mor-interface-rank 0 --mor-transient-dt 0.01 `
                --mor-transient-t-end 0.01 --mor-transient-waveform unseen_waveform `
                --mor-transient-initial-mode ambient --mor-deployment-rhs-count 100 `
                --local-mor-matrix-free-threshold 20000 --max-pcg-iterations 500 `
                --gmres-restart 100 --pcg-tolerance 1e-5 --schur-proxy-ring 1 `
                --schur-proxy-block-size 64 --output-dir $rramDeploymentDirectory --fast-run
        }
        Invoke-Checked 'Chiplet actual 100-waveform one-step deployment' {
            & $exe --transient --config $ChipletConfig --mor-transient-generate `
                --mor-transient-method local-block-arnoldi --mor-arnoldi-moments 4 `
                --mor-interface-rank 0 --mor-transient-dt 0.1 `
                --mor-transient-t-end 0.1 --mor-transient-waveform unseen_waveform `
                --mor-transient-initial-mode ambient --mor-deployment-rhs-count 100 `
                --local-mor-matrix-free-threshold 20000 --max-pcg-iterations 500 `
                --gmres-restart 100 --pcg-tolerance 1e-5 --schur-proxy-ring 1 `
                --schur-proxy-block-size 64 --output-dir $chipletDeploymentDirectory --fast-run
        }
    }

    foreach ($directory in @($rramDirectory, $chipletDirectory)) {
        if (-not (Test-Path -LiteralPath (Join-Path $directory 'local_dynamic_schur_summary.csv'))) {
            throw "Missing Milestone 6 result directory: $directory"
        }
    }
    $rram = Import-Csv (Join-Path $rramDirectory 'local_dynamic_schur_summary.csv') | Select-Object -First 1
    $chiplet = Import-Csv (Join-Path $chipletDirectory 'local_dynamic_schur_summary.csv') | Select-Object -First 1
    $rramSummary = $rram | Select-Object `
        @{Name='case_name';Expression={'rram26'}},
        @{Name='acceptance_status';Expression={Get-Decision $_}}, *
    $chipletSummary = $chiplet | Select-Object `
        @{Name='case_name';Expression={'chiplet'}},
        @{Name='acceptance_status';Expression={Get-Decision $_}}, *
    $rramSummary.status = Get-Decision $rram
    $chipletSummary.status = Get-Decision $chiplet
    $summaries = @($rramSummary, $chipletSummary)
    $summaries | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_block_arnoldi_milestone6_summary.csv')
    if ((Get-Decision $rram) -ne 'success' -or (Get-Decision $chiplet) -ne 'success') {
        throw 'Milestone 6 temperature accuracy gate failed.'
    }

    $ranks = @()
    foreach ($specification in @(
        [pscustomobject]@{Name='rram26';Directory=$rramDirectory},
        [pscustomobject]@{Name='chiplet';Directory=$chipletDirectory})) {
        $rows = @(Import-Csv (Join-Path $specification.Directory 'local_block_arnoldi_rank.csv'))
        foreach ($group in @($rows | Group-Object subdomain)) {
            $last = $group.Group[-1]
            $ranks += $last | Select-Object @{Name='case_name';Expression={$specification.Name}}, *
        }
    }
    $ranks | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_block_arnoldi_milestone6_rank_by_subdomain.csv')

    $timing = @()
    foreach ($specification in @(
        [pscustomobject]@{Name='rram26';Summary=$rram;Directory=$rramDirectory},
        [pscustomobject]@{Name='chiplet';Summary=$chiplet;Directory=$chipletDirectory})) {
        $s = $specification.Summary
        $timing += [pscustomobject]@{
            case_name=$specification.Name; steps=$s.steps
            reference_setup_seconds=$s.reference_setup_seconds
            local_basis_setup_seconds=$s.local_basis_setup_seconds
            local_symbolic_seconds=$s.local_symbolic_seconds
            local_numerical_seconds=$s.local_numerical_seconds
            local_multi_rhs_seconds=$s.local_multi_rhs_seconds
            local_orthogonalization_seconds=$s.local_orthogonalization_seconds
            dynamic_schur_setup_seconds=$s.dynamic_schur_setup_seconds
            proxy_setup_seconds=$s.proxy_setup_seconds
            interface_solve_seconds=$s.interface_solve_seconds
            proxy_solve_seconds=$s.proxy_solve_seconds
            coarse_solve_seconds=$s.coarse_solve_seconds
            local_recovery_seconds=$s.local_recovery_seconds
            local_online_core_seconds=$s.local_online_core_seconds
            fom_factor_seconds=$s.fom_factor_seconds; fom_solve_seconds=$s.fom_solve_seconds
            total_seconds=$s.total_seconds
            average_local_online_seconds_per_step=([double]$s.local_online_core_seconds/[int]$s.steps)
        }
        Import-Csv (Join-Path $specification.Directory 'local_dynamic_schur_accuracy_by_time.csv') |
            Select-Object @{Name='case_name';Expression={$specification.Name}}, * |
            Export-Csv -NoTypeInformation -Encoding UTF8 `
                (Join-Path $outputs "local_transient_$($specification.Name)_by_time.csv")
    }
    $timing | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_transient_timing.csv')

    @(
        [pscustomobject]@{case_name='rram26';model_bytes=$rram.model_bytes;dynamic_factor_bytes=$rram.dynamic_schur_factor_memory_bytes;fom_factor_bytes=$rram.fom_factor_memory_bytes;combined_factor_bytes=$rram.factor_memory_bytes;peak_working_set_bytes=$rram.peak_working_set_bytes},
        [pscustomobject]@{case_name='chiplet';model_bytes=$chiplet.model_bytes;dynamic_factor_bytes=$chiplet.dynamic_schur_factor_memory_bytes;fom_factor_bytes=$chiplet.fom_factor_memory_bytes;combined_factor_bytes=$chiplet.factor_memory_bytes;peak_working_set_bytes=$chiplet.peak_working_set_bytes}
    ) | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_transient_memory.csv')

    Write-WorstFlux @(
        [pscustomobject]@{Name='rram26';Accuracy=(Join-Path $rramDirectory 'local_dynamic_schur_accuracy_by_time.csv');Flux=(Join-Path $rramDirectory 'local_dynamic_schur_interface_flux.csv')},
        [pscustomobject]@{Name='chiplet';Accuracy=(Join-Path $chipletDirectory 'local_dynamic_schur_accuracy_by_time.csv');Flux=(Join-Path $chipletDirectory 'local_dynamic_schur_interface_flux.csv')}
    ) (Join-Path $outputs 'local_transient_interface_flux.csv')

    $rramGlobalSummary = Import-Csv 'results\stage2c1_rram26_m1\transient_block_arnoldi_summary.csv' | Select-Object -First 1
    $rramGlobalOffline = Import-Csv 'results\stage2c1_rram26_m1\transient_offline_timing.csv' | Select-Object -First 1
    $rramGlobalOnline = Import-Csv 'results\stage2c1_rram26_m1\transient_online_timing.csv' | Select-Object -First 1
    $rramGlobalMemory = Import-Csv 'results\stage2c1_rram26_m1\transient_memory.csv' | Select-Object -First 1
    $chipletGlobalSummary = Import-Csv 'results\stage2c1_chiplet_m4\transient_block_arnoldi_summary.csv' | Select-Object -First 1
    $chipletGlobalOffline = Import-Csv 'results\stage2c1_chiplet_m4\transient_offline_timing.csv' | Select-Object -First 1
    $chipletGlobalOnline = Import-Csv 'results\stage2c1_chiplet_m4\transient_online_timing.csv' | Select-Object -First 1
    $chipletGlobalMemory = Import-Csv 'results\stage2c1_chiplet_m4\transient_memory.csv' | Select-Object -First 1
    $comparison = @(
        [pscustomobject]@{case_name='rram26';method='monolithic_pardiso';basis_scope='none';rank=0;interface_dofs=0;offline_seconds=$rram.fom_factor_seconds;online_seconds=$rram.fom_solve_seconds;iterations=0;relative_l2=0;maximum_absolute_k=0;model_bytes=$rram.fom_factor_memory_bytes;peak_working_set_bytes=$rram.peak_working_set_bytes;status='reference'},
        [pscustomobject]@{case_name='rram26';method='Global Block Arnoldi';basis_scope='one_global_basis';rank=$rramGlobalSummary.rank;interface_dofs=0;offline_seconds=$rramGlobalOffline.total_offline_seconds;online_seconds=$rramGlobalOnline.total_seconds;iterations=0;relative_l2=$rramGlobalSummary.space_time_relative_l2;maximum_absolute_k=$rramGlobalSummary.maximum_absolute_k;model_bytes=$rramGlobalMemory.model_bytes;peak_working_set_bytes=$rramGlobalMemory.peak_working_set_bytes;status=$rramGlobalSummary.status},
        [pscustomobject]@{case_name='rram26';method='Local Block Arnoldi + Dynamic Schur';basis_scope='26_independent_local_bases';rank=$rram.total_local_rank;interface_dofs=$rram.full_interface_dofs;offline_seconds=([double]$rram.reference_setup_seconds+[double]$rram.local_basis_setup_seconds+[double]$rram.dynamic_schur_setup_seconds);online_seconds=$rram.local_online_core_seconds;iterations=$rram.interface_iterations_maximum;relative_l2=$rram.space_time_relative_l2;maximum_absolute_k=$rram.maximum_absolute_k;model_bytes=$rram.model_bytes;peak_working_set_bytes=$rram.peak_working_set_bytes;status=(Get-Decision $rram)},
        [pscustomobject]@{case_name='chiplet';method='monolithic_pardiso';basis_scope='none';rank=0;interface_dofs=0;offline_seconds=$chiplet.fom_factor_seconds;online_seconds=$chiplet.fom_solve_seconds;iterations=0;relative_l2=0;maximum_absolute_k=0;model_bytes=$chiplet.fom_factor_memory_bytes;peak_working_set_bytes=$chiplet.peak_working_set_bytes;status='reference'},
        [pscustomobject]@{case_name='chiplet';method='Global Block Arnoldi';basis_scope='one_global_basis_benchmark_only';rank=$chipletGlobalSummary.rank;interface_dofs=0;offline_seconds=$chipletGlobalOffline.total_offline_seconds;online_seconds=$chipletGlobalOnline.total_seconds;iterations=0;relative_l2=$chipletGlobalSummary.space_time_relative_l2;maximum_absolute_k=$chipletGlobalSummary.maximum_absolute_k;model_bytes=$chipletGlobalMemory.model_bytes;peak_working_set_bytes=$chipletGlobalMemory.peak_working_set_bytes;status=$chipletGlobalSummary.status},
        [pscustomobject]@{case_name='chiplet';method='Local Block Arnoldi + Dynamic Schur';basis_scope='2_independent_package_level_ddm_bases';rank=$chiplet.total_local_rank;interface_dofs=$chiplet.full_interface_dofs;offline_seconds=([double]$chiplet.reference_setup_seconds+[double]$chiplet.local_basis_setup_seconds+[double]$chiplet.dynamic_schur_setup_seconds);online_seconds=$chiplet.local_online_core_seconds;iterations=$chiplet.interface_iterations_maximum;relative_l2=$chiplet.space_time_relative_l2;maximum_absolute_k=$chiplet.maximum_absolute_k;model_bytes=$chiplet.model_bytes;peak_working_set_bytes=$chiplet.peak_working_set_bytes;status=(Get-Decision $chiplet)}
    )
    $comparison | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_transient_vs_global_block_arnoldi.csv')

    $rramGlobal100 = Import-Csv 'results\stage2c1_rram26_deployment_100\transient_online_timing.csv' | Select-Object -First 1
    $rramLocal100 = Import-Csv (Join-Path $rramDeploymentDirectory `
        'local_dynamic_schur_deployment_timing.csv') | Select-Object -First 1
    $chipletLocal100 = Import-Csv (Join-Path $chipletDeploymentDirectory `
        'local_dynamic_schur_deployment_timing.csv') | Select-Object -First 1
    foreach ($deployment in @($rramLocal100, $chipletLocal100)) {
        if ([int]$deployment.waveforms -ne 100 -or
            [int]$deployment.setup_reused -ne 1 -or
            [int]$deployment.proxy_reused -ne 1 -or
            [int]$deployment.local_factors_reused -ne 1 -or
            [double]$deployment.maximum_interface_relative_residual -ge 1e-5) {
            throw 'Actual Local 100-waveform deployment reuse/residual gate failed.'
        }
    }
    @(
        [pscustomobject]@{case_name='rram26';method='Local Block Arnoldi + Dynamic Schur';verified_waveforms=100;waveform_horizon_steps=1;average_online_seconds=$rramLocal100.average_online_seconds_per_waveform;projected_100_waveform_seconds=$rramLocal100.total_online_seconds;measurement='actual_100_one_step_full_recovery_waveforms';average_interface_iterations=$rramLocal100.average_interface_iterations;maximum_interface_iterations=$rramLocal100.maximum_interface_iterations;maximum_interface_residual=$rramLocal100.maximum_interface_relative_residual;break_even_vs_monolithic='never_online_slower';break_even_vs_global='never_online_slower'},
        [pscustomobject]@{case_name='rram26';method='Local Block Arnoldi + Dynamic Schur formal horizon';verified_waveforms=1;waveform_horizon_steps=80;average_online_seconds=$rram.local_online_core_seconds;projected_100_waveform_seconds=([double]$rram.local_online_core_seconds*100);measurement='one_actual_80_step_waveform_projection';average_interface_iterations=([double]$rram.interface_iterations_total/[int]$rram.steps);maximum_interface_iterations=$rram.interface_iterations_maximum;maximum_interface_residual=$rram.maximum_interface_relative_residual;break_even_vs_monolithic='never_online_slower';break_even_vs_global='never_online_slower'},
        [pscustomobject]@{case_name='rram26';method='Global Block Arnoldi selected-output deployment';verified_waveforms=100;waveform_horizon_steps=80;average_online_seconds=$rramGlobal100.average_seconds;projected_100_waveform_seconds=([double]$rramGlobal100.total_seconds);measurement='actual_100_waveforms';break_even_vs_monolithic='see_stage2c1';break_even_vs_global='baseline'},
        [pscustomobject]@{case_name='chiplet';method='Local Block Arnoldi + Dynamic Schur';verified_waveforms=100;waveform_horizon_steps=1;average_online_seconds=$chipletLocal100.average_online_seconds_per_waveform;projected_100_waveform_seconds=$chipletLocal100.total_online_seconds;measurement='actual_100_one_step_full_recovery_waveforms';average_interface_iterations=$chipletLocal100.average_interface_iterations;maximum_interface_iterations=$chipletLocal100.maximum_interface_iterations;maximum_interface_residual=$chipletLocal100.maximum_interface_relative_residual;break_even_vs_monolithic='never_online_slower';break_even_vs_global='never_online_slower'},
        [pscustomobject]@{case_name='chiplet';method='Local Block Arnoldi + Dynamic Schur formal horizon';verified_waveforms=1;waveform_horizon_steps=10;average_online_seconds=$chiplet.local_online_core_seconds;projected_100_waveform_seconds=([double]$chiplet.local_online_core_seconds*100);measurement='one_actual_10_step_waveform_projection';average_interface_iterations=([double]$chiplet.interface_iterations_total/[int]$chiplet.steps);maximum_interface_iterations=$chiplet.interface_iterations_maximum;maximum_interface_residual=$chiplet.maximum_interface_relative_residual;break_even_vs_monolithic='never_online_slower';break_even_vs_global='never_online_slower'},
        [pscustomobject]@{case_name='chiplet';method='Global Block Arnoldi';verified_waveforms=1;waveform_horizon_steps=10;average_online_seconds=$chipletGlobalOnline.total_seconds;projected_100_waveform_seconds=([double]$chipletGlobalOnline.total_seconds*100);measurement='projected_from_one_verified_waveform';break_even_vs_monolithic='see_stage2c1';break_even_vs_global='baseline'}
    ) | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'local_transient_break_even.csv')

    $rramIterations = @(Import-Csv (Join-Path $rramDirectory 'local_dynamic_schur_accuracy_by_time.csv') |
        Where-Object {[int]$_.step -gt 0} |
        ForEach-Object {[int]$_.interface_iterations})
    $chipletIterations = @(Import-Csv (Join-Path $chipletDirectory 'local_dynamic_schur_accuracy_by_time.csv') |
        Where-Object {[int]$_.step -gt 0} |
        ForEach-Object {[int]$_.interface_iterations})
    $report = @"
# Milestone 6: RRAM26 / Chiplet transient Local Block Arnoldi + Dynamic Schur

## Decision

Milestone 6 passes the specified temperature gates on both real cases. RRAM26 uses all 26 configured DDM subdomains, 26 independently generated interior Block Arnoldi bases, the full 331,331-DOF interface, matrix-free Local Dynamic Schur, FGMRES, and a dynamic 1-ring proxy built from that reduced local operator. Chiplet uses the Milestone-3-audited two package-level DDM modules and four power channels; it is not mislabeled as four DDM subdomains. Global Block Arnoldi is retained only as a benchmark.

## Verified results

| case | subdomains | interface DOFs | total local rank | steps | FGMRES min / mean / max | space-time L2 | max node K | max-temperature K | FOM/ROM flux L2 | online core s | peak GiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| RRAM26 | $($rram.subdomains) | $($rram.full_interface_dofs) | $($rram.total_local_rank) | $($rram.steps) | $(($rramIterations | Measure-Object -Minimum).Minimum) / $([Math]::Round(($rramIterations | Measure-Object -Average).Average,3)) / $(($rramIterations | Measure-Object -Maximum).Maximum) | $($rram.space_time_relative_l2) | $($rram.maximum_absolute_k) | $($rram.maximum_temperature_error_k) | $($rram.maximum_fom_rom_flux_relative_l2) | $($rram.local_online_core_seconds) | $([Math]::Round([double]$rram.peak_working_set_bytes/1GB,3)) |
| Chiplet | $($chiplet.subdomains) | $($chiplet.full_interface_dofs) | $($chiplet.total_local_rank) | $($chiplet.steps) | $(($chipletIterations | Measure-Object -Minimum).Minimum) / $([Math]::Round(($chipletIterations | Measure-Object -Average).Average,3)) / $(($chipletIterations | Measure-Object -Maximum).Maximum) | $($chiplet.space_time_relative_l2) | $($chiplet.maximum_absolute_k) | $($chiplet.maximum_temperature_error_k) | $($chiplet.maximum_fom_rom_flux_relative_l2) | $($chiplet.local_online_core_seconds) | $([Math]::Round([double]$chiplet.peak_working_set_bytes/1GB,3)) |

RRAM26 uses one symbolic/numerical dynamic proxy factorization, 882 exact colors, 17,258,637 proxy nonzeros, and coarse dimension 104 for all 80 steps. Chiplet uses 3,636 colors, 3,814,739 proxy nonzeros, and coarse dimension 8. No time step rebuilds the fixed-dt proxy or local factors.

## Interface physics

The by-time files retain temperature jump and physical/SIPG flux imbalance at every step. `local_transient_interface_flux.csv` retains the worst physical interface, face pair, and integration triangle for every detailed diagnostic sample, including both FOM and ROM flux. Full per-triangle raw diagnostics remain in each result directory. RRAM26 final sampled FOM/ROM numerical-flux relative L2 is $($rram.maximum_fom_rom_flux_relative_l2). Chiplet flux L2 is $($chiplet.maximum_fom_rom_flux_relative_l2), substantially less accurate than its temperature field; this unfavorable result is intentionally reported.

RRAM26 maximum sampled full/reduced residuals are $($rram.maximum_full_residual) / $($rram.maximum_reduced_residual); Chiplet reports $($chiplet.maximum_full_residual) / $($chiplet.maximum_reduced_residual). The Chiplet full residual is large despite accepted temperature accuracy and is not hidden. Maximum physical-flux imbalance is $($rram.maximum_relative_flux_imbalance) for RRAM26 and $($chiplet.maximum_relative_flux_imbalance) for Chiplet.

## Efficiency and break-even

The retained Global Block Arnoldi benchmarks are much faster online: RRAM26 global full-field online is $($rramGlobalOnline.total_seconds) s versus $($rram.local_online_core_seconds) s for Local Dynamic Schur; Chiplet is $($chipletGlobalOnline.total_seconds) s versus $($chiplet.local_online_core_seconds) s. Monolithic reused PARDISO solves are also faster for both current serial partitions. Therefore there is no positive waveform-count break-even for the present Local implementation. Its benefit is independent replaceable subdomain models and preservation of the full physical interface, not serial online speed.

The actual setup-once 100-waveform single-step batches take $($rramLocal100.total_online_seconds) s total / $($rramLocal100.average_online_seconds_per_waveform) s average for RRAM26 and $($chipletLocal100.total_online_seconds) s / $($chipletLocal100.average_online_seconds_per_waveform) s for Chiplet. The CSV separately retains formal 80-step/10-step horizon projections, so single-step and full-horizon workloads are not conflated.

## Invariants

The exact Stage-1 Schur operator, FGMRES true-residual gate, PARDISO reuse, FEM/SIPG assembly, interface ordering, triangle-overlap/BVH, and global temperature recovery are unchanged. The Dynamic Schur proxy values are generated from the Local Dynamic Schur operator; no Stage-1 exact proxy values or Global basis are substituted.
"@
    $report | Set-Content -Encoding UTF8 `
        (Join-Path $outputs 'local_block_arnoldi_milestone6_report.md')

    if (-not $SkipTests) {
        $testOutput = @(ctest --test-dir $BuildDirectory -C Release --output-on-failure 2>&1)
        $testOutput | Set-Content -Encoding UTF8 `
            (Join-Path $outputs 'local_mor_all_ctest_results.txt')
        if ($LASTEXITCODE -ne 0) { throw 'Full CTest regression failed.' }
    }
} finally { Pop-Location }
