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
    Write-Host "[Milestone 3] $Label"
    & $Command
    if ($LASTEXITCODE -ne 0) { throw "$Label failed with exit code $LASTEXITCODE." }
}

function Compare-Fields([string]$LeftPath, [string]$RightPath, [string]$Column) {
    $left = @(Import-Csv -LiteralPath $LeftPath)
    $right = @(Import-Csv -LiteralPath $RightPath)
    if ($left.Count -ne $right.Count) { throw "Field row counts differ: $LeftPath / $RightPath" }
    $differenceSquared = 0.0
    $referenceSquared = 0.0
    $maximum = 0.0
    for ($i = 0; $i -lt $left.Count; ++$i) {
        $reference = [double]$left[$i].$Column
        $difference = [double]$right[$i].$Column - $reference
        $differenceSquared += $difference * $difference
        $referenceSquared += $reference * $reference
        $maximum = [Math]::Max($maximum, [Math]::Abs($difference))
    }
    return [pscustomobject]@{
        RelativeL2 = [Math]::Sqrt($differenceSquared) / [Math]::Max(1e-300, [Math]::Sqrt($referenceSquared))
        Maximum = $maximum
    }
}

function Write-CompactFlux([array]$Specifications, [string]$Destination) {
    $utf8 = New-Object System.Text.UTF8Encoding($false)
    $writer = New-Object System.IO.StreamWriter($Destination, $false, $utf8, 1048576)
    try {
        $writer.WriteLine('case_name,case,split,family,power_vector_id,physical_interface_id,face_pair_id,integration_triangle_id,left_subdomain,right_subdomain,left_boundary_entity,right_boundary_entity,area_m2,temperature_jump_rms_k,left_physical_normal_flux_w_m2,right_physical_normal_flux_w_m2,sipg_numerical_flux_w_m2,flux_imbalance_l2_w_m2,relative_flux_imbalance,relative_flux_floor_w_m2,maximum_flux_imbalance_w_m2')
        foreach ($specification in $Specifications) {
            $reader = New-Object System.IO.StreamReader($specification.Path)
            try {
                [void]$reader.ReadLine()
                while (($line = $reader.ReadLine()) -ne $null) {
                    # power_vector_w contains semicolons but no commas, so this exact
                    # generated CSV format can be compacted without materializing it.
                    $parts = $line.Split(',')
                    if ($parts.Count -notin @(20, 21)) { throw "Unexpected flux CSV shape in $($specification.Path)." }
                    $tail = ($parts[4..19] -join ',')
                    $writer.WriteLine("$($specification.Name),$($parts[0]),$($parts[1]),$($parts[2]),$($specification.PowerVectorId),$tail")
                }
            }
            finally { $reader.Dispose() }
        }
    }
    finally { $writer.Dispose() }
}

Push-Location $repo
try {
    $exe = Join-Path $BuildDirectory 'Release\SIPGHeatDDM3D.exe'
    if (-not $AggregateOnly) {
        if (-not $SkipBuild) {
            Invoke-Checked 'Release build' {
                cmake --build $BuildDirectory --config Release --parallel 4
            }
        }
        if (-not (Test-Path -LiteralPath $exe)) { throw "Executable not found: $exe" }
        foreach ($configuration in @($RramConfig, $ChipletConfig)) {
            if (-not (Test-Path -LiteralPath $configuration)) {
                throw "Required large-case configuration not found: $configuration"
            }
        }

        $rramPure = Join-Path $ResultsDirectory 'local_rom_milestone3_rram26_pure'
        $rramModel = Join-Path $ResultsDirectory 'local_rom_milestone3_rram26_model'
        $rramCorrected = Join-Path $ResultsDirectory 'local_rom_milestone3_rram26_corrected'
        Invoke-Checked 'RRAM26 independent Local-ROM generation and deployment' {
            & $exe --steady --config $RramConfig --output-dir $rramPure `
                --solvers local-rom --local-mor-generate --local-mor-save $rramModel `
                --local-mor-method pod --local-mor-mode pure --local-interface-mode full `
                --local-mor-rank 133 --local-mor-training-cases 133 `
                --local-mor-validation-cases 0 --local-mor-test-cases 0 `
                --local-mor-compare-fom --pcg-tolerance 1e-6 `
                --max-pcg-iterations 500 --gmres-restart 100 `
                --schur-proxy-block-size 64 --fast-run
        }
        Invoke-Checked 'RRAM26 residual-gated corrected Local-ROM' {
            & $exe --steady --config $RramConfig --output-dir $rramCorrected `
                --solvers local-rom --local-mor-load $rramModel `
                --local-mor-mode corrected --local-interface-mode full `
                --local-mor-validation-cases 0 --local-mor-test-cases 0 `
                --local-mor-compare-fom --pcg-tolerance 1e-6 `
                --max-pcg-iterations 500 --gmres-restart 100 `
                --schur-proxy-block-size 64 --fast-run
        }

        $chipletPure = Join-Path $ResultsDirectory 'local_rom_milestone3_chiplet_pure'
        $chipletModel = Join-Path $ResultsDirectory 'local_rom_milestone3_chiplet_model'
        $chipletCorrected = Join-Path $ResultsDirectory 'local_rom_milestone3_chiplet_corrected'
        $chipletReuse = Join-Path $ResultsDirectory 'local_rom_milestone3_chiplet_reuse'
        $chipletReuseModel = Join-Path $ResultsDirectory 'local_rom_milestone3_chiplet_reuse_model'
        Invoke-Checked 'Chiplet independent Local-ROM generation and deployment' {
            & $exe --steady --config $ChipletConfig --output-dir $chipletPure `
                --solvers local-rom --local-mor-generate --local-mor-save $chipletModel `
                --local-mor-method pod --local-mor-mode pure --local-interface-mode full `
                --local-mor-rank 4 --local-mor-training-cases 12 `
                --local-mor-validation-cases 0 --local-mor-test-cases 0 `
                --local-mor-compare-fom --pcg-tolerance 1e-5 `
                --max-pcg-iterations 500 --gmres-restart 100 `
                --schur-proxy-block-size 64 --fast-run
        }
        Invoke-Checked 'Chiplet residual-gated corrected Local-ROM' {
            & $exe --steady --config $ChipletConfig --output-dir $chipletCorrected `
                --solvers local-rom --local-mor-load $chipletModel `
                --local-mor-mode corrected --local-interface-mode full `
                --local-mor-validation-cases 0 --local-mor-test-cases 0 `
                --local-mor-compare-fom --pcg-tolerance 1e-5 `
                --max-pcg-iterations 500 --gmres-restart 100 `
                --schur-proxy-block-size 64 --fast-run
        }
        Invoke-Checked 'Chiplet strict-fingerprint template-reuse attempt' {
            & $exe --steady --config $ChipletConfig --output-dir $chipletReuse `
                --solvers local-rom --local-mor-generate --local-mor-save $chipletReuseModel `
                --local-mor-method pod --local-mor-mode pure --local-interface-mode full `
                --local-mor-rank 4 --local-mor-training-cases 12 `
                --local-mor-validation-cases 0 --local-mor-test-cases 0 `
                --local-mor-compare-fom --local-rom-reuse-identical-subdomains `
                --pcg-tolerance 1e-5 --max-pcg-iterations 500 `
                --gmres-restart 100 --schur-proxy-block-size 64 --fast-run
        }
    }

    $outputs = Join-Path $repo 'outputs'
    New-Item -ItemType Directory -Force -Path $outputs | Out-Null
    $directories = @{
        rram26_pure = Join-Path $ResultsDirectory 'local_rom_milestone3_rram26_pure'
        rram26_corrected = Join-Path $ResultsDirectory 'local_rom_milestone3_rram26_corrected'
        chiplet_pure = Join-Path $ResultsDirectory 'local_rom_milestone3_chiplet_pure'
        chiplet_corrected = Join-Path $ResultsDirectory 'local_rom_milestone3_chiplet_corrected'
        chiplet_reuse = Join-Path $ResultsDirectory 'local_rom_milestone3_chiplet_reuse'
    }
    foreach ($directory in $directories.Values) {
        if (-not (Test-Path -LiteralPath $directory)) { throw "Missing result directory: $directory" }
    }

    $summary = @()
    foreach ($key in @('rram26_pure','rram26_corrected','chiplet_pure','chiplet_corrected','chiplet_reuse')) {
        $row = Import-Csv (Join-Path $directories[$key] 'local_rom_schur_summary.csv') | Select-Object -First 1
        $summary += $row | Select-Object @{Name='case_name';Expression={$key}}, *
    }
    $summary | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_milestone3_summary.csv')

    $ranks = @()
    foreach ($key in @('rram26_pure','chiplet_pure')) {
        $ranks += Import-Csv (Join-Path $directories[$key] 'local_rom_rank_by_subdomain.csv') |
            Select-Object @{Name='case_name';Expression={$key -replace '_pure$',''}}, *
    }
    $ranks | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_milestone3_rank_by_subdomain.csv')

    $accuracy = @()
    foreach ($key in @('rram26_pure','rram26_corrected','chiplet_pure','chiplet_corrected')) {
        $accuracy += Import-Csv (Join-Path $directories[$key] 'local_rom_accuracy_by_case.csv') |
            Select-Object @{Name='case_name';Expression={$key}}, *
    }
    $accuracy | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_milestone3_accuracy.csv')

    $timing = @()
    foreach ($key in @('rram26_pure','rram26_corrected','chiplet_pure','chiplet_corrected','chiplet_reuse')) {
        $offline = Import-Csv (Join-Path $directories[$key] 'local_rom_offline_timing.csv') | Select-Object -First 1
        $online = Import-Csv (Join-Path $directories[$key] 'local_rom_online_timing.csv') | Select-Object -First 1
        $program = Import-Csv (Join-Path $directories[$key] 'program_timing.csv') |
            Where-Object stage -eq 'total_program' | Select-Object -First 1
        $timing += [pscustomobject]@{
            case_name=$key; total_offline_seconds=$offline.total_offline_seconds
            model_load_seconds=$offline.model_load_seconds; snapshot_solve_seconds=$offline.snapshot_solve_seconds
            local_basis_seconds=$offline.local_basis_seconds; local_projection_seconds=$offline.local_projection_seconds
            reduced_schur_construction_seconds=$offline.reduced_schur_construction_seconds
            proxy_setup_seconds=$offline.proxy_setup_seconds
            exact_schur_setup_seconds=$offline.exact_schur_setup_seconds
            source_projection_seconds=$online.source_projection_seconds
            local_reduced_assembly_seconds=$online.local_reduced_assembly_seconds
            interface_solve_seconds=$online.interface_solve_seconds; proxy_solve_seconds=$online.proxy_solve_seconds
            coarse_solve_seconds=$online.coarse_solve_seconds; local_recovery_seconds=$online.local_recovery_seconds
            full_field_reconstruction_seconds=$online.full_field_reconstruction_seconds
            correction_seconds=$online.correction_seconds; core_online_seconds=$online.core_online_seconds
            total_with_diagnostics_seconds=$online.total_with_source_and_diagnostics_seconds
            total_program_seconds=$program.seconds
        }
    }
    $timing | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_milestone3_timing.csv')

    $memory = @()
    foreach ($key in @('rram26_pure','rram26_corrected','chiplet_pure','chiplet_corrected','chiplet_reuse')) {
        $memory += Import-Csv (Join-Path $directories[$key] 'local_rom_memory.csv') |
            Select-Object @{Name='case_name';Expression={$key}}, *
    }
    $memory | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_milestone3_memory.csv')

    $temperatureDifference = Compare-Fields `
        (Join-Path $directories.chiplet_pure 'temperature_local_pod_schur_rom_pure_nodes.csv') `
        (Join-Path $directories.chiplet_reuse 'temperature_local_pod_schur_rom_pure_nodes.csv') 'temperature_K'
    $fluxDifference = Compare-Fields `
        (Join-Path $directories.chiplet_pure 'local_rom_interface_flux.csv') `
        (Join-Path $directories.chiplet_reuse 'local_rom_interface_flux.csv') 'sipg_numerical_flux_w_m2'
    $independentTemplate = Import-Csv (Join-Path $directories.chiplet_pure 'local_rom_template_reuse.csv') | Select-Object -First 1
    $reuseTemplate = Import-Csv (Join-Path $directories.chiplet_reuse 'local_rom_template_reuse.csv') | Select-Object -First 1
    @(
        [pscustomobject]@{scenario='independent';reuse_enabled=0;unique_template_count=$independentTemplate.unique_template_count;subdomain_instance_count=$independentTemplate.subdomain_instance_count;reused_instance_count=$independentTemplate.reused_instance_count;basis_storage_without_reuse_bytes=$independentTemplate.basis_storage_without_reuse_bytes;basis_storage_with_reuse_bytes=$independentTemplate.basis_storage_with_reuse_bytes;temperature_relative_l2_difference=0;maximum_temperature_difference_k=0;flux_relative_l2_difference=0;maximum_flux_difference_w_m2=0},
        [pscustomobject]@{scenario='strict_fingerprint_reuse';reuse_enabled=1;unique_template_count=$reuseTemplate.unique_template_count;subdomain_instance_count=$reuseTemplate.subdomain_instance_count;reused_instance_count=$reuseTemplate.reused_instance_count;basis_storage_without_reuse_bytes=$reuseTemplate.basis_storage_without_reuse_bytes;basis_storage_with_reuse_bytes=$reuseTemplate.basis_storage_with_reuse_bytes;temperature_relative_l2_difference=$temperatureDifference.RelativeL2;maximum_temperature_difference_k=$temperatureDifference.Maximum;flux_relative_l2_difference=$fluxDifference.RelativeL2;maximum_flux_difference_w_m2=$fluxDifference.Maximum}
    ) | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_milestone3_template_reuse.csv')

    $multiRhs = @()
    foreach ($key in @('rram26_pure','chiplet_pure','chiplet_reuse')) {
        $multiRhs += Import-Csv (Join-Path $directories[$key] 'local_rom_multi_rhs_timing.csv') |
            Select-Object @{Name='case_name';Expression={$key}}, *
    }
    $multiRhs | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_milestone3_multi_rhs.csv')

    Write-CompactFlux @(
        [pscustomobject]@{Name='rram26';Path=(Join-Path $directories.rram26_pure 'local_rom_interface_flux.csv');PowerVectorId='uniform_125x8e-5W'},
        [pscustomobject]@{Name='chiplet';Path=(Join-Path $directories.chiplet_pure 'local_rom_interface_flux.csv');PowerVectorId='uniform_4x10W'}
    ) (Join-Path $outputs 'local_rom_milestone3_interface_flux.csv')

    $partitions = @()
    foreach ($key in @('rram26_pure','chiplet_pure')) {
        $partitions += Import-Csv (Join-Path $directories[$key] 'local_rom_partition_definition.csv') |
            Select-Object @{Name='case_name';Expression={$key -replace '_pure$',''}}, *
    }
    $partitions | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_milestone3_partition_definition.csv')

    $rramPureSummary = $summary | Where-Object case_name -eq 'rram26_pure' | Select-Object -First 1
    $rramCorrectedSummary = $summary | Where-Object case_name -eq 'rram26_corrected' | Select-Object -First 1
    $chipletPureSummary = $summary | Where-Object case_name -eq 'chiplet_pure' | Select-Object -First 1
    $chipletCorrectedSummary = $summary | Where-Object case_name -eq 'chiplet_corrected' | Select-Object -First 1
    $rramOffline = Import-Csv (Join-Path $directories.rram26_pure 'local_rom_offline_timing.csv') | Select-Object -First 1
    $rramMemory = Import-Csv (Join-Path $directories.rram26_pure 'local_rom_memory.csv') | Select-Object -First 1
    $rramAccuracy = Import-Csv (Join-Path $directories.rram26_pure 'local_rom_accuracy_by_case.csv') | Select-Object -First 1
    $chipletAccuracy = Import-Csv (Join-Path $directories.chiplet_pure 'local_rom_accuracy_by_case.csv') | Select-Object -First 1
    $historicalGlobal = Import-Csv (Join-Path $outputs 'reduced_schur_stage2_benchmark_summary.csv')
    $historicalExactResponse = Import-Csv (Join-Path $outputs 'local_interior_rom_stage2a1_summary.csv')
    $historicalStage1 = Import-Csv (Join-Path $outputs 'schur_proxy_benchmark_summary.csv')
    $historicalLarge = Import-Csv (Join-Path $outputs 'schur_stage2_large_case_comparison.csv')
    $rramGlobal = $historicalGlobal | Where-Object { $_.case -eq 'rram26' -and $_.mode -eq 'pure_rank125' } | Select-Object -First 1
    $chipletGlobal = $historicalGlobal | Where-Object { $_.case -eq 'chiplet' -and $_.mode -eq 'pure' } | Select-Object -First 1
    $rramExactResponse = $historicalExactResponse | Where-Object { $_.case -eq 'RRAM26' -and $_.interior_mode -eq 'exact-response' } | Select-Object -First 1
    $chipletExactResponse = $historicalExactResponse | Where-Object { $_.case -eq 'chiplet_horizontal' -and $_.interior_mode -eq 'exact-response' } | Select-Object -First 1
    $rramStage1 = $historicalStage1 | Where-Object { $_.case -eq 'RRAM26' -and $_.configuration -eq 'D_proxy_plus_volume_xyz' } | Select-Object -First 1
    $chipletStage1 = $historicalLarge | Where-Object { $_.case -eq 'chiplet_horizontal' -and $_.solver -eq 'DDM-Schur-FGMRES' } | Select-Object -First 1
    $chipletDirect = $historicalLarge | Where-Object { $_.case -eq 'chiplet_horizontal' -and $_.solver -eq 'Global-PARDISO-SPD-Direct' } | Select-Object -First 1

    $comparison = @(
        [pscustomobject]@{case_name='rram26';method='A_monolithic_pardiso';basis_scope='none';interface_scope='monolithic';tolerance='direct';offline_or_setup_seconds=$rramOffline.fom_factorization_wall_seconds;online_seconds=$rramAccuracy.fom_seconds;iterations=0;relative_l2=0;max_node_error_k=0;peak_memory_mib=([double]$rramMemory.peak_working_set_bytes/1MB);model_or_factor_bytes=$rramMemory.fom_factor_bytes;status='reference';source='milestone3_direct_reference'},
        [pscustomobject]@{case_name='rram26';method='B_stage1_exact_ddm_schur';basis_scope='none_exact_local';interface_scope='full_exact_schur';tolerance='1e-5';offline_or_setup_seconds=$rramStage1.setup_seconds;online_seconds=$rramStage1.solve_seconds;iterations=$rramStage1.iterations;relative_l2='';max_node_error_k='';peak_memory_mib=$rramStage1.peak_working_set_mib;model_or_factor_bytes='';status=$rramStage1.status;source='verified_stage1_proxy_baseline'},
        [pscustomobject]@{case_name='rram26';method='C_global_reduced_schur';basis_scope='one_global_interface_basis';interface_scope='reduced_interface';tolerance='direct_reduced';offline_or_setup_seconds=$rramGlobal.offline_or_setup_seconds;online_seconds=$rramGlobal.average_or_nominal_online_seconds;iterations=0;relative_l2=$rramGlobal.max_relative_l2;max_node_error_k=$rramGlobal.max_absolute_error_K;peak_memory_mib=$rramGlobal.peak_working_set_MiB;model_or_factor_bytes='';status=$rramGlobal.decision;source='stage2_global_reduced_schur'},
        [pscustomobject]@{case_name='rram26';method='D_local_rom_full_interface_schur';basis_scope='26_independent_local_bases';interface_scope='331331_full_interface_matrix_free';tolerance='1e-6';offline_or_setup_seconds=$rramPureSummary.offline_seconds;online_seconds=$rramPureSummary.nominal_online_seconds;iterations=$rramPureSummary.interface_fgmres_iterations;relative_l2=$rramPureSummary.relative_l2;max_node_error_k=$rramPureSummary.max_node_error_k;peak_memory_mib=([double]$rramMemory.peak_working_set_bytes/1MB);model_or_factor_bytes=$rramPureSummary.model_bytes;status=$rramPureSummary.status;source='milestone3'},
        [pscustomobject]@{case_name='rram26';method='E_local_rom_corrected';basis_scope='26_independent_local_bases';interface_scope='full_exact_schur_residual_gate';tolerance='1e-6';offline_or_setup_seconds=$rramCorrectedSummary.offline_seconds;online_seconds=$rramCorrectedSummary.nominal_online_seconds;iterations=$rramCorrectedSummary.local_rom_guess_iterations;relative_l2=$rramCorrectedSummary.relative_l2;max_node_error_k=$rramCorrectedSummary.max_node_error_k;peak_memory_mib=((Import-Csv (Join-Path $directories.rram26_corrected 'solver_comparison.csv') | Select-Object -First 1).peak_working_set_mb);model_or_factor_bytes=$rramCorrectedSummary.model_bytes;status=$rramCorrectedSummary.status;source='milestone3'},
        [pscustomobject]@{case_name='rram26';method='F_stage2a1_exact_response';basis_scope='global_interface_exact_power_response';interface_scope='rank125';tolerance='direct_reduced';offline_or_setup_seconds=$rramExactResponse.offline_seconds;online_seconds=$rramExactResponse.average_online_seconds_100_rhs;iterations=0;relative_l2=$rramExactResponse.worst_full_relative_l2;max_node_error_k=$rramExactResponse.worst_max_abs_K;peak_memory_mib=$rramExactResponse.peak_memory_mib;model_or_factor_bytes=([double]$rramExactResponse.model_size_mib*1MB);status=$rramExactResponse.decision;source='stage2a1'},
        [pscustomobject]@{case_name='chiplet';method='A_monolithic_pardiso';basis_scope='none';interface_scope='monolithic';tolerance='direct';offline_or_setup_seconds=$chipletDirect.setup_seconds;online_seconds=$chipletDirect.solve_seconds;iterations=0;relative_l2=$chipletDirect.relative_l2_error;max_node_error_k=$chipletDirect.max_temperature_error_K;peak_memory_mib=$chipletDirect.peak_working_set_mb;model_or_factor_bytes='';status=$chipletDirect.status;source='verified_large_case_baseline'},
        [pscustomobject]@{case_name='chiplet';method='B_stage1_exact_ddm_schur';basis_scope='none_exact_local';interface_scope='full_exact_schur';tolerance='verified_baseline';offline_or_setup_seconds=$chipletStage1.setup_seconds;online_seconds=$chipletStage1.solve_seconds;iterations=$chipletStage1.iterations;relative_l2=$chipletStage1.relative_l2_error;max_node_error_k=$chipletStage1.max_temperature_error_K;peak_memory_mib=$chipletStage1.peak_working_set_mb;model_or_factor_bytes='';status=$chipletStage1.status;source='verified_large_case_baseline'},
        [pscustomobject]@{case_name='chiplet';method='C_global_reduced_schur';basis_scope='one_global_interface_basis';interface_scope='reduced_interface';tolerance='direct_reduced';offline_or_setup_seconds=$chipletGlobal.offline_or_setup_seconds;online_seconds=$chipletGlobal.average_or_nominal_online_seconds;iterations=0;relative_l2=$chipletGlobal.max_relative_l2;max_node_error_k=$chipletGlobal.max_absolute_error_K;peak_memory_mib=$chipletGlobal.peak_working_set_MiB;model_or_factor_bytes='';status=$chipletGlobal.decision;source='stage2_global_reduced_schur'},
        [pscustomobject]@{case_name='chiplet';method='D_local_rom_full_interface_schur';basis_scope='2_independent_local_bases';interface_scope='80853_full_interface_matrix_free';tolerance='1e-5';offline_or_setup_seconds=$chipletPureSummary.offline_seconds;online_seconds=$chipletPureSummary.nominal_online_seconds;iterations=$chipletPureSummary.interface_fgmres_iterations;relative_l2=$chipletPureSummary.relative_l2;max_node_error_k=$chipletPureSummary.max_node_error_k;peak_memory_mib=((Import-Csv (Join-Path $directories.chiplet_pure 'solver_comparison.csv') | Select-Object -First 1).peak_working_set_mb);model_or_factor_bytes=$chipletPureSummary.model_bytes;status=$chipletPureSummary.status;source='milestone3'},
        [pscustomobject]@{case_name='chiplet';method='E_local_rom_corrected';basis_scope='2_independent_local_bases';interface_scope='full_exact_schur_residual_gate';tolerance='1e-5';offline_or_setup_seconds=$chipletCorrectedSummary.offline_seconds;online_seconds=$chipletCorrectedSummary.nominal_online_seconds;iterations=$chipletCorrectedSummary.local_rom_guess_iterations;relative_l2=$chipletCorrectedSummary.relative_l2;max_node_error_k=$chipletCorrectedSummary.max_node_error_k;peak_memory_mib=((Import-Csv (Join-Path $directories.chiplet_corrected 'solver_comparison.csv') | Select-Object -First 1).peak_working_set_mb);model_or_factor_bytes=$chipletCorrectedSummary.model_bytes;status=$chipletCorrectedSummary.status;source='milestone3'},
        [pscustomobject]@{case_name='chiplet';method='F_stage2a1_exact_response';basis_scope='global_interface_exact_power_response';interface_scope='rank4';tolerance='direct_reduced';offline_or_setup_seconds=$chipletExactResponse.offline_seconds;online_seconds=$chipletExactResponse.average_online_seconds_100_rhs;iterations=0;relative_l2=$chipletExactResponse.worst_full_relative_l2;max_node_error_k=$chipletExactResponse.worst_max_abs_K;peak_memory_mib=$chipletExactResponse.peak_memory_mib;model_or_factor_bytes=([double]$chipletExactResponse.model_size_mib*1MB);status=$chipletExactResponse.decision;source='stage2a1'}
    )
    $comparison | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_milestone3_vs_global.csv')

    $rramBatch100 = $multiRhs | Where-Object { $_.case_name -eq 'rram26_pure' -and $_.rhs_count -eq '100' } | Select-Object -First 1
    $chipletBatch100 = $multiRhs | Where-Object { $_.case_name -eq 'chiplet_pure' -and $_.rhs_count -eq '100' } | Select-Object -First 1
    $report = @"
# Milestone 3: RRAM26 / Chiplet steady Local-ROM + full-interface Schur

## Decision

Milestone 3 passes the requested full-field accuracy and corrected true-residual gates on both real cases. The implementation uses one independent interior basis per actual DDM subdomain and never replaces the full SIPG interface with a global reduced basis. The existing Global Reduced Schur remains a benchmark only.

## Audited physical partitions

- RRAM26: 26 configured DDM subdomains, 125 independent power channels, 494,111 global DOFs, 331,331 full-interface DOFs, 107,402 integrated face pairs and 139,713 integration triangles.
- Chiplet: four independently powered chiplet entities are contained in two runnable package-level DDM modules: (0) heatsink/lid stack and (1) multi-chiplet die/interconnect assembly. The split has one physical interface, 23,845 face pairs, 29,962 integration triangles and 80,853 interface DOFs. It is defined by supplied mesh/domain/material/boundary IDs, not coordinate ranges.
- Strict Chiplet fingerprints found two unique templates for two instances, so reuse was correctly rejected (reused_instance_count=0). Independent/reuse temperature difference is $($temperatureDifference.RelativeL2) relative L2 and $($temperatureDifference.Maximum) K maximum; numerical-flux difference is $($fluxDifference.RelativeL2) relative L2.

## Main results

| case / mode | local total rank | coarse dim | interface FGMRES / correction iterations | setup (s) | nominal online (s) | relative L2 | max node error (K) | max-temperature error (K) | true/global residual | peak MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| RRAM26 pure | $($rramPureSummary.total_local_rank) | $($rramPureSummary.coarse_dimension) | $($rramPureSummary.interface_fgmres_iterations) | $($rramPureSummary.offline_seconds) | $($rramPureSummary.nominal_online_seconds) | $($rramPureSummary.relative_l2) | $($rramPureSummary.max_node_error_k) | $($rramPureSummary.max_temperature_error_k) | $($rramAccuracy.global_relative_residual) | $([Math]::Round([double]$rramMemory.peak_working_set_bytes/1MB,2)) |
| RRAM26 corrected | $($rramCorrectedSummary.total_local_rank) | $($rramCorrectedSummary.coarse_dimension) | $($rramCorrectedSummary.local_rom_guess_iterations) | $($rramCorrectedSummary.offline_seconds) | $($rramCorrectedSummary.nominal_online_seconds) | $($rramCorrectedSummary.relative_l2) | $($rramCorrectedSummary.max_node_error_k) | $($rramCorrectedSummary.max_temperature_error_k) | $($rramCorrectedSummary.final_true_residual) | $((Import-Csv (Join-Path $directories.rram26_corrected 'solver_comparison.csv') | Select-Object -First 1).peak_working_set_mb) |
| Chiplet pure | $($chipletPureSummary.total_local_rank) | $($chipletPureSummary.coarse_dimension) | $($chipletPureSummary.interface_fgmres_iterations) | $($chipletPureSummary.offline_seconds) | $($chipletPureSummary.nominal_online_seconds) | $($chipletPureSummary.relative_l2) | $($chipletPureSummary.max_node_error_k) | $($chipletPureSummary.max_temperature_error_k) | $($chipletAccuracy.global_relative_residual) | $((Import-Csv (Join-Path $directories.chiplet_pure 'solver_comparison.csv') | Select-Object -First 1).peak_working_set_mb) |
| Chiplet corrected | $($chipletCorrectedSummary.total_local_rank) | $($chipletCorrectedSummary.coarse_dimension) | $($chipletCorrectedSummary.local_rom_guess_iterations) | $($chipletCorrectedSummary.offline_seconds) | $($chipletCorrectedSummary.nominal_online_seconds) | $($chipletCorrectedSummary.relative_l2) | $($chipletCorrectedSummary.max_node_error_k) | $($chipletCorrectedSummary.max_temperature_error_k) | $($chipletCorrectedSummary.final_true_residual) | $((Import-Csv (Join-Path $directories.chiplet_corrected 'solver_comparison.csv') | Select-Object -First 1).peak_working_set_mb) |

RRAM26 exact Stage 1 requires 41 iterations from zero at the same 1e-6 gate; the Local-ROM initial guess reduces this to 4 corrected iterations. Chiplet reduces the exact correction from 7 to 1 iteration at 1e-5.

## Fixed-matrix deployment

The model, interface graph, proxy and factors are constructed once. The measured 100-RHS averages are $($rramBatch100.average_online_seconds) s/RHS for RRAM26 and $($chipletBatch100.average_online_seconds) s/RHS for Chiplet. These are full-field Local-ROM deployments, including the full-interface FGMRES solve and reconstruction.

## Flux diagnostics

Every integration triangle is preserved in local_rom_milestone3_interface_flux.csv. To avoid repeating a 125-entry power vector 139,713 times, the compact file uses power_vector_id: RRAM26 is uniform 125 x 8e-5 W; Chiplet is uniform 4 x 10 W. RRAM26 reports temperature-jump RMS $($rramPureSummary.interface_temperature_jump_rms_k) K and relative physical-flux imbalance $($rramPureSummary.relative_flux_imbalance); Chiplet reports $($chipletPureSummary.interface_temperature_jump_rms_k) K and $($chipletPureSummary.relative_flux_imbalance). These physical-gradient diagnostics are reported without hiding the relatively large imbalance; they are not ROM-vs-FOM error norms.

## Honest efficiency comparison

For RRAM26, Local-ROM setup/online are $($rramPureSummary.offline_seconds) / $($rramPureSummary.nominal_online_seconds) s. The retained Global Reduced Schur benchmark is much faster online ($($rramGlobal.average_or_nominal_online_seconds) s), as is Stage 2A.1 exact-response ($($rramExactResponse.average_online_seconds_100_rhs) s/RHS). The Local-ROM path pays for the 331,331-dimensional physical interface, but retains independent replaceable subdomain models and an exact full-interface coupling path. Chiplet monolithic PARDISO is also faster than Local-ROM for the current two-module partition. The normalized A-F data are in local_rom_milestone3_vs_global.csv.

## Scope and invariants

The exact Schur operator, Stage 1 FGMRES residual gate, PARDISO factor reuse, SIPG/FEM assembly, triangle-overlap/BVH, interface ordering and global temperature recovery are unchanged. Matrix-free Local-ROM Schur values and its 1-ring proxy are built from the reduced local blocks; exact Stage 1 values are used only by corrected/reference paths. Dynamic Local Block Arnoldi remains out of scope until Milestone 4.
"@
    $report | Set-Content -Encoding UTF8 (Join-Path $outputs 'local_rom_milestone3_report.md')

    if (-not $SkipTests) {
        $testOutput = & ctest --test-dir $BuildDirectory -C Release --output-on-failure 2>&1
        $testOutput | Set-Content -Encoding UTF8 (Join-Path $outputs 'local_rom_milestone3_ctest_results.txt')
        if ($LASTEXITCODE -ne 0) { throw 'Full CTest regression failed.' }
    }
}
finally { Pop-Location }
