param(
    [string]$OutputDirectory = 'outputs',
    [string]$DocsDirectory = 'docs'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$outputs = [System.IO.Path]::GetFullPath((Join-Path $project $OutputDirectory))
$docs = [System.IO.Path]::GetFullPath((Join-Path $project $DocsDirectory))
$invariant = [System.Globalization.CultureInfo]::InvariantCulture
$na = 'N/A'
$maximumSourceBytes = 2MB

function Relative-Path([string]$Path) {
    $root = $project.TrimEnd('\') + '\'
    $full = [System.IO.Path]::GetFullPath($Path)
    if ($full.StartsWith($root, [System.StringComparison]::OrdinalIgnoreCase)) {
        return $full.Substring($root.Length).Replace('\','/')
    }
    return $full.Replace('\','/')
}

function Read-SmallCsv([string]$RelativePath) {
    $path = Join-Path $project $RelativePath
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing source CSV: $RelativePath" }
    $item = Get-Item -LiteralPath $path
    if ($item.Length -gt $maximumSourceBytes) {
        throw "Refusing to read source larger than $maximumSourceBytes bytes: $RelativePath"
    }
    return @(Import-Csv -LiteralPath $path)
}

function Read-SmallText([string]$RelativePath) {
    $path = Join-Path $project $RelativePath
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing source text: $RelativePath" }
    $item = Get-Item -LiteralPath $path
    if ($item.Length -gt $maximumSourceBytes) {
        throw "Refusing to read source larger than $maximumSourceBytes bytes: $RelativePath"
    }
    return [System.IO.File]::ReadAllText($path)
}

function Read-GitCsv([string]$Commit, [string]$RelativePath) {
    $spec = "${Commit}:$RelativePath"
    $lines = @(& git show $spec 2>&1)
    if ($LASTEXITCODE -ne 0) { throw "Cannot read archived CSV $spec`n$($lines -join [Environment]::NewLine)" }
    return @(($lines -join "`n") | ConvertFrom-Csv)
}

function Assert-GitObject([string]$Commit, [string]$RelativePath) {
    & git cat-file -e "${Commit}:$RelativePath" 2>$null
    if ($LASTEXITCODE -ne 0) { throw "Missing archived source ${Commit}:$RelativePath" }
}

function D($Value) {
    if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value) -or $Value -eq $na) {
        return [double]::NaN
    }
    return [double]::Parse([string]$Value, $invariant)
}

function F($Value) {
    if ($Value -is [string] -and $Value -eq $na) { return $na }
    $number = D $Value
    if ([double]::IsNaN($number)) { return $na }
    return $number.ToString('G17', $invariant)
}

function I($Value) {
    if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value) -or $Value -eq $na) {
        return $na
    }
    return ([int64](D $Value)).ToString($invariant)
}

function Maximum([object[]]$Rows, [string]$Column) {
    $values = @($Rows | ForEach-Object { D $_.$Column } | Where-Object { -not [double]::IsNaN($_) })
    if ($values.Count -eq 0) { return $na }
    return F (($values | Measure-Object -Maximum).Maximum)
}

function Rank-Stats([object[]]$Rows, [string]$RankColumn) {
    $values = @($Rows | ForEach-Object { D $_.$RankColumn })
    return [pscustomobject]@{
        minimum = I (($values | Measure-Object -Minimum).Minimum)
        maximum = I (($values | Measure-Object -Maximum).Maximum)
        total = I (($values | Measure-Object -Sum).Sum)
    }
}

function Final-Arnoldi-Ranks([object[]]$Rows, [string]$CaseName = '') {
    $selected = $Rows
    if ($CaseName) { $selected = @($Rows | Where-Object case_name -eq $CaseName) }
    return @($selected | Group-Object subdomain | ForEach-Object {
        $_.Group | Select-Object -Last 1
    })
}

function Interface-Pair-Count([object[]]$PartitionRows, [string]$CaseName) {
    $pairs = [System.Collections.Generic.HashSet[string]]::new()
    foreach ($row in @($PartitionRows | Where-Object case_name -eq $CaseName)) {
        if ([string]::IsNullOrWhiteSpace($row.neighbor_subdomains)) { continue }
        foreach ($neighborText in $row.neighbor_subdomains.Split(';')) {
            if ([string]::IsNullOrWhiteSpace($neighborText)) { continue }
            $left = [int]$row.subdomain_id
            $right = [int]$neighborText
            $key = if ($left -lt $right) { "$left-$right" } else { "$right-$left" }
            [void]$pairs.Add($key)
        }
    }
    return $pairs.Count
}

function Power-Channel-Count([object[]]$PartitionRows, [string]$CaseName) {
    return [int](($PartitionRows | Where-Object case_name -eq $CaseName |
        Measure-Object power_channel_count -Sum).Sum)
}

function MiB($Bytes) {
    $number = D $Bytes
    if ([double]::IsNaN($number)) { return $na }
    return F ($number / 1MB)
}

function Ratio($Numerator, $Denominator) {
    $n = D $Numerator; $d = D $Denominator
    if ([double]::IsNaN($n) -or [double]::IsNaN($d) -or $d -eq 0.0) { return $na }
    return F ($n / $d)
}

function Ceil-BreakEven($CurrentOffline, $BaselineOffline, $BaselineOnline, $CurrentOnline) {
    $denominator = (D $BaselineOnline) - (D $CurrentOnline)
    if ($denominator -le 0.0) { return $na }
    return I ([Math]::Ceiling([Math]::Max(0.0,
        ((D $CurrentOffline) - (D $BaselineOffline)) / $denominator)))
}

function New-MasterRow([hashtable]$Values) {
    $columns = @(
        'milestone_id','case_name','steady_or_transient','measurement_mode',
        'interior_reduction','interface_treatment','interface_solver','basis_construction',
        'corrected_mode','main_purpose','global_dofs','subdomain_count',
        'physical_interface_count','full_interface_dofs','reduced_interface_dofs',
        'interface_compression_ratio','local_rank_total','local_rank_min','local_rank_max',
        'port_rank_total','power_channel_count','time_step_count','dt_s','waveform_count',
        'relative_l2','max_node_error_k','max_temperature_error_k',
        'max_temperature_curve_rmse_k','interface_temperature_jump','flux_relative_l2',
        'maximum_flux_imbalance_w_m2','sampled_full_residual','projected_system_residual',
        'corrected_relative_l2','corrected_full_residual','offline_total_s',
        'snapshot_generation_s','basis_generation_s','reduced_projection_s','schur_setup_s',
        'model_load_s','online_core_s','full_field_reconstruction_s','online_total_s',
        'rhs_count','average_online_s_per_case','baseline_name','baseline_online_s_per_case',
        'speedup_vs_baseline','offline_break_even_waveforms','offline_peak_memory_bytes',
        'offline_peak_memory_mib','online_peak_memory_bytes','online_peak_memory_mib',
        'serialized_model_size_bytes','basis_storage_bytes','template_reuse_reduction_percent',
        'speedup_vs_monolithic','speedup_vs_stage1_schur','speedup_vs_global_block_arnoldi',
        'speedup_vs_milestone6','commit','tag','branch','commit_message','release_status',
        'ctest_passed','ctest_total','reproduction_script','report_file',
        'workspace_clean_at_milestone','accuracy_conclusion','notes')
    $ordered = [ordered]@{}
    foreach ($column in $columns) {
        $ordered[$column] = if ($Values.ContainsKey($column)) { $Values[$column] } else { $na }
    }
    return [pscustomobject]$ordered
}

Push-Location $project
try {
    New-Item -ItemType Directory -Force -Path $outputs,$docs | Out-Null

    # This script intentionally reads only compact, committed summaries.  It
    # never opens full-field temperature or per-triangle flux result files and
    # never invokes SIPGHeatDDM3D, CTest, or another solver.
    $m1Commit = '852c8f3'; $m2Commit = 'da5f874'; $m3Commit = 'dc2d267'
    $m4Commit = '1c50542'; $m5Commit = 'c60d152'; $m6Commit = '6f2d35f'
    $m7Commit = 'd1755bc'

    $m1Summary = (Read-GitCsv $m1Commit 'outputs/local_rom_schur_summary.csv')[0]
    $m1Accuracy = Read-GitCsv $m1Commit 'outputs/local_rom_accuracy_by_case.csv'
    $m1Offline = (Read-GitCsv $m1Commit 'outputs/local_rom_offline_timing.csv')[0]
    $m1Online = (Read-GitCsv $m1Commit 'outputs/local_rom_online_timing.csv' |
        Where-Object { $_.case -eq '0' -and $_.split -eq 'nominal' })[0]
    $m1Memory = (Read-GitCsv $m1Commit 'outputs/local_rom_memory.csv')[0]
    $m1Ranks = Read-GitCsv $m1Commit 'outputs/local_rom_rank_by_subdomain.csv'
    $m1RankStats = Rank-Stats $m1Ranks 'selected_rank'
    $m1Comparison = Read-GitCsv $m1Commit 'outputs/local_rom_vs_global_rom.csv'
    $m1Corrected = ($m1Comparison | Where-Object method -eq 'Local-POD-Schur-ROM-Corrected')[0]
    $m1Stage1 = ($m1Comparison | Where-Object method -eq 'DDM-Schur-FGMRES')[0]
    $m1Monolithic = ($m1Comparison | Where-Object method -eq 'Global-PARDISO-SPD-Direct')[0]

    $m2Summary = (Read-SmallCsv 'outputs/local_rom_schur_summary.csv')[0]
    $m2Accuracy = Read-SmallCsv 'outputs/local_rom_accuracy_by_case.csv'
    $m2OfflineRows = Read-SmallCsv 'outputs/local_rom_offline_timing.csv'
    $m2Offline = ($m2OfflineRows | Where-Object scenario -eq 'generation')[0]
    $m2Reload = ($m2OfflineRows | Where-Object scenario -eq 'serialized_model_reload')[0]
    $m2Multi = (Read-SmallCsv 'outputs/local_rom_multi_rhs_timing.csv' |
        Where-Object rhs_count -eq '100')[0]
    $m2Memory = (Read-SmallCsv 'outputs/local_rom_memory.csv')[0]
    $m2Ranks = Read-SmallCsv 'outputs/local_rom_rank_by_subdomain.csv'
    $m2RankStats = Rank-Stats $m2Ranks 'selected_rank'
    $m2Comparison = Read-SmallCsv 'outputs/local_rom_vs_global_rom.csv'
    $m2Corrected = ($m2Comparison | Where-Object method -eq 'Local-POD-Schur-ROM-Corrected')[0]

    $m3Summary = Read-SmallCsv 'outputs/local_rom_milestone3_summary.csv'
    $m3Timing = Read-SmallCsv 'outputs/local_rom_milestone3_timing.csv'
    $m3Memory = Read-SmallCsv 'outputs/local_rom_milestone3_memory.csv'
    $m3Multi = Read-SmallCsv 'outputs/local_rom_milestone3_multi_rhs.csv'
    $m3Ranks = Read-SmallCsv 'outputs/local_rom_milestone3_rank_by_subdomain.csv'
    $m3Partition = Read-SmallCsv 'outputs/local_rom_milestone3_partition_definition.csv'

    $m4Summary = Read-SmallCsv 'outputs/local_block_arnoldi_milestone4_summary.csv'
    $m4Local = ($m4Summary | Where-Object method -eq 'Local Block Arnoldi + Dynamic Schur')[0]
    $m4Monolithic = ($m4Summary | Where-Object method -eq 'Transient monolithic PARDISO')[0]
    $m4Dynamic = Read-SmallCsv 'outputs/local_dynamic_schur_milestone4.csv'
    $m4Representative = ($m4Dynamic | Where-Object {
        $_.waveform -eq 'single_step' -and $_.initial_mode -eq 'ambient' })[0]
    $m4Ranks = Final-Arnoldi-Ranks (Read-SmallCsv 'outputs/local_block_arnoldi_milestone4_rank.csv')
    $m4RankStats = Rank-Stats $m4Ranks 'cumulative_rank'

    $m5Summary = Read-SmallCsv 'outputs/local_block_arnoldi_milestone5_summary.csv'
    $m5Local = ($m5Summary | Where-Object method -eq 'Local Block Arnoldi + Dynamic Schur')[0]
    $m5Monolithic = ($m5Summary | Where-Object method -eq 'Transient monolithic PARDISO')[0]
    $m5Dynamic = Read-SmallCsv 'outputs/local_dynamic_schur_milestone5.csv'
    $m5Ranks = Final-Arnoldi-Ranks (Read-SmallCsv 'outputs/local_block_arnoldi_milestone5_rank_by_subdomain.csv')
    $m5RankStats = Rank-Stats $m5Ranks 'cumulative_rank'

    $m6Summary = Read-SmallCsv 'outputs/local_block_arnoldi_milestone6_summary.csv'
    $m6BreakEven = Read-SmallCsv 'outputs/local_transient_break_even.csv'
    $m6Memory = Read-SmallCsv 'outputs/local_transient_memory.csv'
    $m6RanksAll = Read-SmallCsv 'outputs/local_block_arnoldi_milestone6_rank_by_subdomain.csv'
    $m6Comparison = Read-SmallCsv 'outputs/local_transient_vs_global_block_arnoldi.csv'

    $m7Summary = Read-SmallCsv 'outputs/local_port_rom_summary.csv'
    $m7Timing = Read-SmallCsv 'outputs/local_port_dynamic_schur_timing.csv'
    $m7Memory = Read-SmallCsv 'outputs/local_port_memory.csv'
    $m7Comparison = Read-SmallCsv 'outputs/local_port_vs_full_interface.csv'
    $m7Enrichment = Read-SmallCsv 'outputs/local_port_enrichment_history.csv'
    $m7Ports = Read-SmallCsv 'outputs/local_port_rank_by_interface.csv'

    $master = [System.Collections.Generic.List[object]]::new()

    $m1WorstMaxNode = Maximum $m1Accuracy 'max_node_error_k'
    $m1WorstMaxTemperature = Maximum $m1Accuracy 'max_temperature_error_k'
    $m1WorstJump = Maximum $m1Accuracy 'interface_temperature_jump_rms_k'
    $m1WorstFlux = Maximum $m1Accuracy 'maximum_flux_imbalance_w_m2'
    $m1WorstFullResidual = Maximum $m1Accuracy 'global_relative_residual'
    $m1WorstProjectedResidual = Maximum $m1Accuracy 'interface_relative_residual'
    $m1Stage1Speedup = Ratio $m1Stage1.solve_seconds $m1Summary.nominal_online_seconds
    $m1MonolithicSpeedup = Ratio $m1Monolithic.solve_seconds $m1Summary.nominal_online_seconds
    $master.Add((New-MasterRow @{
        milestone_id='M1';case_name='two-cube';steady_or_transient='steady';
        measurement_mode='single_case;pure;full_field_recovery';
        interior_reduction='2 independent Local POD/RB bases';
        interface_treatment='full SIPG interface';interface_solver='sparse PARDISO';
        basis_construction='local snapshot POD';corrected_mode='optional exact Schur correction';
        main_purpose='Establish two-cube steady Local-ROM Schur path';
        global_dofs=I $m1Summary.global_dofs;subdomain_count=I $m1Summary.subdomains;
        physical_interface_count='1';full_interface_dofs=I $m1Summary.interface_dofs;
        reduced_interface_dofs=I $m1Summary.interface_dofs;interface_compression_ratio='1';
        local_rank_total=I $m1Summary.total_local_rank;local_rank_min=$m1RankStats.minimum;
        local_rank_max=$m1RankStats.maximum;port_rank_total='0';power_channel_count=$na;
        time_step_count='1';dt_s=$na;waveform_count=$na;
        relative_l2=F $m1Summary.relative_l2;max_node_error_k=$m1WorstMaxNode;
        max_temperature_error_k=$m1WorstMaxTemperature;interface_temperature_jump=$m1WorstJump;
        flux_relative_l2=$na;maximum_flux_imbalance_w_m2=$m1WorstFlux;
        sampled_full_residual=$m1WorstFullResidual;projected_system_residual=$m1WorstProjectedResidual;
        corrected_relative_l2=F $m1Corrected.relative_l2;
        corrected_full_residual='3.1159430437664524e-12';offline_total_s=F $m1Offline.total_offline_seconds;
        snapshot_generation_s=F $m1Offline.snapshot_solve_seconds;
        basis_generation_s=F $m1Offline.local_basis_seconds;
        reduced_projection_s=F $m1Offline.local_projection_seconds;
        schur_setup_s=F ((D $m1Offline.reduced_schur_construction_seconds)+(D $m1Offline.reduced_factorization_seconds));
        model_load_s=$na;online_core_s=F $m1Online.total_seconds;
        full_field_reconstruction_s=F $m1Online.full_field_reconstruction_seconds;
        online_total_s=F $m1Online.total_seconds;rhs_count='1';
        average_online_s_per_case=F $m1Online.total_seconds;
        baseline_name='Stage 1 exact DDM-Schur-FGMRES';
        baseline_online_s_per_case=F $m1Stage1.solve_seconds;speedup_vs_baseline=$m1Stage1Speedup;
        offline_break_even_waveforms=$na;offline_peak_memory_bytes=I $m1Memory.peak_working_set_bytes;
        offline_peak_memory_mib=MiB $m1Memory.peak_working_set_bytes;
        online_peak_memory_bytes=$na;online_peak_memory_mib=$na;
        serialized_model_size_bytes=I $m1Memory.model_bytes;basis_storage_bytes=$na;
        template_reuse_reduction_percent=$na;speedup_vs_monolithic=$m1MonolithicSpeedup;
        speedup_vs_stage1_schur=$m1Stage1Speedup;speedup_vs_global_block_arnoldi=$na;
        speedup_vs_milestone6=$na;commit=$m1Commit;tag=$na;
        branch='codex/feature/local-mor-schur-mainline';
        commit_message='Implement local steady ROM Schur milestone 1';release_status='success';
        ctest_passed='22';ctest_total='22';reproduction_script='scripts/run_local_rom_schur.ps1';
        report_file='git:852c8f3:outputs/local_rom_schur_report.md';
        workspace_clean_at_milestone=$na;accuracy_conclusion='pass';
        notes='Power-channel count and separate online peak memory were not recorded in the archived M1 outputs.'
    }))

    $m2WorstMaxNode = Maximum $m2Accuracy 'max_node_error_k'
    $m2WorstMaxTemperature = Maximum $m2Accuracy 'max_temperature_error_k'
    $m2WorstJump = Maximum $m2Accuracy 'interface_temperature_jump_rms_k'
    $m2WorstFlux = Maximum $m2Accuracy 'maximum_flux_imbalance_w_m2'
    $m2WorstFullResidual = Maximum $m2Accuracy 'global_relative_residual'
    $m2WorstProjectedResidual = Maximum $m2Accuracy 'interface_relative_residual'
    $master.Add((New-MasterRow @{
        milestone_id='M2';case_name='ten-cube';steady_or_transient='steady';
        measurement_mode='100_rhs_average;pure;full_field_recovery';
        interior_reduction='10 independent Local POD/RB bases';
        interface_treatment='full SIPG interface';interface_solver='sparse PARDISO';
        basis_construction='local snapshot POD';corrected_mode='optional exact Schur correction';
        main_purpose='Scale steady Local-ROM Schur to ten subdomains';
        global_dofs=I $m2Summary.global_dofs;subdomain_count=I $m2Summary.subdomains;
        physical_interface_count='9';full_interface_dofs=I $m2Summary.interface_dofs;
        reduced_interface_dofs=I $m2Summary.interface_dofs;interface_compression_ratio='1';
        local_rank_total=I $m2Summary.total_local_rank;local_rank_min=$m2RankStats.minimum;
        local_rank_max=$m2RankStats.maximum;port_rank_total='0';power_channel_count='10';
        time_step_count='1';dt_s=$na;waveform_count=$na;
        relative_l2=F $m2Summary.relative_l2;max_node_error_k=$m2WorstMaxNode;
        max_temperature_error_k=$m2WorstMaxTemperature;interface_temperature_jump=$m2WorstJump;
        flux_relative_l2=$na;maximum_flux_imbalance_w_m2=$m2WorstFlux;
        sampled_full_residual=$m2WorstFullResidual;projected_system_residual=$m2WorstProjectedResidual;
        corrected_relative_l2=F $m2Corrected.relative_l2;
        corrected_full_residual='1.1697903240601744e-11';offline_total_s=F $m2Offline.total_offline_seconds;
        snapshot_generation_s=F $m2Offline.snapshot_solve_seconds;
        basis_generation_s=F $m2Offline.local_basis_seconds;
        reduced_projection_s=F $m2Offline.local_projection_seconds;
        schur_setup_s=F ((D $m2Offline.reduced_schur_construction_seconds)+(D $m2Offline.reduced_factorization_seconds));
        model_load_s=F $m2Reload.model_load_seconds;online_core_s=F $m2Multi.online_total_seconds;
        full_field_reconstruction_s=F $m2Multi.full_field_reconstruction_total_seconds;
        online_total_s=F $m2Multi.online_total_seconds;rhs_count='100';
        average_online_s_per_case=F $m2Multi.average_online_seconds;
        baseline_name=$na;baseline_online_s_per_case=$na;speedup_vs_baseline=$na;
        offline_break_even_waveforms=$na;offline_peak_memory_bytes=I $m2Memory.peak_working_set_bytes;
        offline_peak_memory_mib=MiB $m2Memory.peak_working_set_bytes;
        online_peak_memory_bytes=I $m2Reload.peak_working_set_bytes;
        online_peak_memory_mib=MiB $m2Reload.peak_working_set_bytes;
        serialized_model_size_bytes=I $m2Memory.model_bytes;
        basis_storage_bytes=I $m2Memory.basis_storage_with_reuse_bytes;
        template_reuse_reduction_percent='31.12';speedup_vs_monolithic=$na;
        speedup_vs_stage1_schur=$na;speedup_vs_global_block_arnoldi=$na;
        speedup_vs_milestone6=$na;commit=$m2Commit;tag='baseline-local-rom-schur-milestone2';
        branch='codex/feature/local-mor-schur-mainline';
        commit_message='Complete ten-cube local ROM Schur milestone 2';release_status='success';
        ctest_passed='32';ctest_total='32';reproduction_script='scripts/run_local_rom_schur.ps1';
        report_file='outputs/local_rom_schur_report.md';workspace_clean_at_milestone=$na;
        accuracy_conclusion='pass';
        notes='No cross-method speedup is reported for the 100-RHS average because archived baselines use a different single-case timing definition.'
    }))

    foreach ($caseName in @('rram26','chiplet')) {
        $pureName = "${caseName}_pure"; $correctedName = "${caseName}_corrected"
        $pure = ($m3Summary | Where-Object case_name -eq $pureName)[0]
        $corrected = ($m3Summary | Where-Object case_name -eq $correctedName)[0]
        $timing = ($m3Timing | Where-Object case_name -eq $pureName)[0]
        $memory = ($m3Memory | Where-Object case_name -eq $pureName)[0]
        $batch = ($m3Multi | Where-Object { $_.case_name -eq $pureName -and $_.rhs_count -eq '100' })[0]
        $rankStats = Rank-Stats @($m3Ranks | Where-Object case_name -eq $caseName) 'selected_rank'
        $interfaces = Interface-Pair-Count $m3Partition $caseName
        $channels = Power-Channel-Count $m3Partition $caseName
        $master.Add((New-MasterRow @{
            milestone_id='M3';case_name=$caseName;steady_or_transient='steady';
            measurement_mode='100_rhs_average;pure;full_field_recovery';
            interior_reduction="$($pure.subdomains) independent Local POD/RB bases";
            interface_treatment='full SIPG interface';interface_solver=$pure.interface_solver;
            basis_construction='local snapshot POD';corrected_mode='optional exact Schur correction';
            main_purpose='Validate steady Local-ROM Schur on real RRAM/Chiplet models';
            global_dofs=I $pure.global_dofs;subdomain_count=I $pure.subdomains;
            physical_interface_count=I $interfaces;full_interface_dofs=I $pure.interface_dofs;
            reduced_interface_dofs=I $pure.interface_dofs;interface_compression_ratio='1';
            local_rank_total=I $pure.total_local_rank;local_rank_min=$rankStats.minimum;
            local_rank_max=$rankStats.maximum;port_rank_total='0';power_channel_count=I $channels;
            time_step_count='1';dt_s=$na;waveform_count=$na;relative_l2=F $pure.relative_l2;
            max_node_error_k=F $pure.max_node_error_k;
            max_temperature_error_k=F $pure.max_temperature_error_k;
            interface_temperature_jump=F $pure.interface_temperature_jump_rms_k;
            flux_relative_l2=$na;maximum_flux_imbalance_w_m2=$na;
            sampled_full_residual=F $pure.initial_true_residual;
            projected_system_residual=F $pure.interface_fgmres_true_residual;
            corrected_relative_l2=F $corrected.relative_l2;
            corrected_full_residual=F $corrected.final_true_residual;
            offline_total_s=F $pure.offline_seconds;snapshot_generation_s=F $timing.snapshot_solve_seconds;
            basis_generation_s=F $timing.local_basis_seconds;
            reduced_projection_s=F $timing.local_projection_seconds;
            schur_setup_s=F $timing.reduced_schur_construction_seconds;model_load_s='0';
            online_core_s=F $batch.online_total_seconds;
            full_field_reconstruction_s=F $batch.full_field_reconstruction_total_seconds;
            online_total_s=F $batch.online_total_seconds;rhs_count='100';
            average_online_s_per_case=F $batch.average_online_seconds;
            baseline_name=$na;baseline_online_s_per_case=$na;speedup_vs_baseline=$na;
            offline_break_even_waveforms=$na;offline_peak_memory_bytes=I $memory.peak_working_set_bytes;
            offline_peak_memory_mib=MiB $memory.peak_working_set_bytes;
            online_peak_memory_bytes=$na;online_peak_memory_mib=$na;
            serialized_model_size_bytes=I $memory.model_bytes;
            basis_storage_bytes=I $memory.basis_storage_with_reuse_bytes;
            template_reuse_reduction_percent='0';speedup_vs_monolithic=$na;
            speedup_vs_stage1_schur=$na;speedup_vs_global_block_arnoldi=$na;
            speedup_vs_milestone6=$na;commit=$m3Commit;tag=$na;
            branch='codex/feature/local-mor-schur-mainline';
            commit_message='Complete RRAM26 and Chiplet steady local ROM milestone 3';
            release_status='success';ctest_passed='33';ctest_total='33';
            reproduction_script='scripts/run_local_rom_milestone3.ps1';
            report_file='outputs/local_rom_milestone3_report.md';workspace_clean_at_milestone=$na;
            accuracy_conclusion='pass';
            notes='Physical-gradient flux imbalance is archived, but a FOM/ROM flux relative-L2 error was not recorded for M3.'
        }))
    }

    $master.Add((New-MasterRow @{
        milestone_id='M4';case_name='two-cube';steady_or_transient='transient';
        measurement_mode='single_case;10_step_horizon;pure;full_field_recovery';
        interior_reduction='2 independent Local Block Arnoldi bases';
        interface_treatment='full SIPG interface';interface_solver='explicit Dynamic Schur PARDISO';
        basis_construction='local Block Arnoldi';corrected_mode='steady corrected regression only';
        main_purpose='Establish transient local descriptor and fixed-dt Dynamic Schur';
        global_dofs=I $m4Representative.global_dofs;subdomain_count=I $m4Representative.subdomains;
        physical_interface_count='1';full_interface_dofs=I $m4Representative.full_interface_dofs;
        reduced_interface_dofs=I $m4Representative.full_interface_dofs;interface_compression_ratio='1';
        local_rank_total=I $m4Local.rank;local_rank_min=$m4RankStats.minimum;
        local_rank_max=$m4RankStats.maximum;port_rank_total='0';power_channel_count='2';
        time_step_count=I $m4Representative.steps;dt_s=F $m4Representative.dt_s;
        waveform_count=I (@($m4Dynamic | Select-Object -ExpandProperty waveform -Unique).Count);
        relative_l2=Maximum $m4Dynamic 'space_time_relative_l2';
        max_node_error_k=Maximum $m4Dynamic 'maximum_absolute_k';
        max_temperature_error_k=Maximum $m4Dynamic 'maximum_temperature_error_k';
        interface_temperature_jump=Maximum $m4Dynamic 'maximum_temperature_jump_rms_k';
        flux_relative_l2=$na;maximum_flux_imbalance_w_m2=$na;
        sampled_full_residual=Maximum $m4Dynamic 'maximum_full_residual';
        projected_system_residual=Maximum $m4Dynamic 'maximum_reduced_residual';
        corrected_relative_l2=$na;corrected_full_residual=$na;
        offline_total_s=F $m4Local.offline_seconds;snapshot_generation_s=$na;
        basis_generation_s=F $m4Representative.local_basis_setup_seconds;
        reduced_projection_s=$na;schur_setup_s=F $m4Representative.dynamic_schur_setup_seconds;
        model_load_s=$na;online_core_s=F $m4Local.online_seconds;
        full_field_reconstruction_s=F $m4Representative.local_recovery_seconds;
        online_total_s=F $m4Local.online_seconds;rhs_count='1';
        average_online_s_per_case=F $m4Local.online_seconds;
        baseline_name='Transient monolithic PARDISO';
        baseline_online_s_per_case=F $m4Monolithic.online_seconds;
        speedup_vs_baseline=Ratio $m4Monolithic.online_seconds $m4Local.online_seconds;
        offline_break_even_waveforms=$na;offline_peak_memory_bytes=$na;
        offline_peak_memory_mib=$na;online_peak_memory_bytes=$na;online_peak_memory_mib=$na;
        serialized_model_size_bytes=$na;basis_storage_bytes=$na;
        template_reuse_reduction_percent=$na;
        speedup_vs_monolithic=Ratio $m4Monolithic.online_seconds $m4Local.online_seconds;
        speedup_vs_stage1_schur=$na;speedup_vs_global_block_arnoldi=$na;
        speedup_vs_milestone6=$na;commit=$m4Commit;tag=$na;
        branch='codex/feature/local-mor-schur-mainline';
        commit_message='Complete two-cube local Block Arnoldi milestone 4';release_status='success';
        ctest_passed='36';ctest_total='36';
        reproduction_script='scripts/run_local_block_arnoldi_milestone4.ps1';
        report_file='outputs/local_block_arnoldi_milestone4_report.md';
        workspace_clean_at_milestone=$na;accuracy_conclusion='pass';
        notes='Worst accuracy is aggregated across the validation suite; timing uses the archived representative 10-step local row.'
    }))

    $master.Add((New-MasterRow @{
        milestone_id='M5';case_name='ten-cube';steady_or_transient='transient';
        measurement_mode='single_case;10_step_horizon;pure;full_field_recovery';
        interior_reduction='10 independent Local Block Arnoldi bases';
        interface_treatment='full SIPG interface';interface_solver='explicit Dynamic Schur PARDISO';
        basis_construction='local Block Arnoldi';corrected_mode='steady corrected regression only';
        main_purpose='Scale transient Local Block Arnoldi to ten subdomains';
        global_dofs='60640';subdomain_count='10';physical_interface_count='9';
        full_interface_dofs='10593';reduced_interface_dofs='10593';interface_compression_ratio='1';
        local_rank_total=I $m5Local.rank;local_rank_min=$m5RankStats.minimum;
        local_rank_max=$m5RankStats.maximum;port_rank_total='0';power_channel_count='10';
        time_step_count='10';dt_s='0.10000000000000001';waveform_count='10';
        relative_l2=Maximum $m5Dynamic 'space_time_relative_l2';
        max_node_error_k=Maximum $m5Dynamic 'maximum_absolute_k';
        max_temperature_error_k=Maximum $m5Dynamic 'maximum_temperature_error_k';
        interface_temperature_jump=Maximum $m5Dynamic 'maximum_temperature_jump_rms_k';
        flux_relative_l2=$na;maximum_flux_imbalance_w_m2=$na;
        sampled_full_residual=Maximum $m5Dynamic 'maximum_full_residual';
        projected_system_residual=Maximum $m5Dynamic 'maximum_reduced_residual';
        corrected_relative_l2=$na;corrected_full_residual=$na;
        offline_total_s=F $m5Local.offline_seconds;snapshot_generation_s=$na;
        basis_generation_s='1.0830991000000001';reduced_projection_s=$na;
        schur_setup_s='1.7842393000000001';model_load_s=$na;
        online_core_s=F $m5Local.online_seconds;full_field_reconstruction_s=$na;
        online_total_s=F $m5Local.online_seconds;rhs_count='1';
        average_online_s_per_case=F $m5Local.online_seconds;
        baseline_name='Transient monolithic PARDISO';
        baseline_online_s_per_case=F $m5Monolithic.online_seconds;
        speedup_vs_baseline=Ratio $m5Monolithic.online_seconds $m5Local.online_seconds;
        offline_break_even_waveforms=$na;offline_peak_memory_bytes=I $m5Local.peak_working_set_bytes;
        offline_peak_memory_mib=MiB $m5Local.peak_working_set_bytes;
        online_peak_memory_bytes=$na;online_peak_memory_mib=$na;
        serialized_model_size_bytes=I $m5Local.model_bytes;basis_storage_bytes=$na;
        template_reuse_reduction_percent='0';
        speedup_vs_monolithic=Ratio $m5Monolithic.online_seconds $m5Local.online_seconds;
        speedup_vs_stage1_schur=$na;speedup_vs_global_block_arnoldi=$na;
        speedup_vs_milestone6=$na;commit=$m5Commit;tag=$na;
        branch='codex/feature/local-mor-schur-mainline';
        commit_message='Complete ten-cube local Block Arnoldi milestone 5';release_status='success';
        ctest_passed='38';ctest_total='38';
        reproduction_script='scripts/run_local_block_arnoldi_milestone5.ps1';
        report_file='outputs/local_block_arnoldi_milestone5_report.md';
        workspace_clean_at_milestone=$na;accuracy_conclusion='pass';
        notes='Strict template fingerprints rejected reuse; the archived M5 report does not separate full-field reconstruction from online core.'
    }))

    foreach ($caseName in @('rram26','chiplet')) {
        $summary = ($m6Summary | Where-Object case_name -eq $caseName)[0]
        $batch = ($m6BreakEven | Where-Object {
            $_.case_name -eq $caseName -and $_.measurement -eq 'actual_100_one_step_full_recovery_waveforms' })[0]
        $memory = ($m6Memory | Where-Object case_name -eq $caseName)[0]
        $rankStats = Rank-Stats (Final-Arnoldi-Ranks $m6RanksAll $caseName) 'cumulative_rank'
        $mono = ($m6Comparison | Where-Object {
            $_.case_name -eq $caseName -and $_.method -eq 'monolithic_pardiso' })[0]
        $global = ($m6Comparison | Where-Object {
            $_.case_name -eq $caseName -and $_.method -eq 'Global Block Arnoldi' })[0]
        $master.Add((New-MasterRow @{
            milestone_id='M6';case_name=$caseName;steady_or_transient='transient';
            measurement_mode='100_rhs_average;one_step;pure;full_field_recovery';
            interior_reduction="$($summary.subdomains) independent Local Block Arnoldi bases";
            interface_treatment='full SIPG interface';interface_solver='matrix-free FGMRES + 1-ring proxy';
            basis_construction='local Block Arnoldi';corrected_mode='not used for M6 acceptance';
            main_purpose='Validate transient Local Block Arnoldi on real models and expose full-interface bottleneck';
            global_dofs=I $summary.global_dofs;subdomain_count=I $summary.subdomains;
            physical_interface_count=I (Interface-Pair-Count $m3Partition $caseName);
            full_interface_dofs=I $summary.full_interface_dofs;
            reduced_interface_dofs=I $summary.full_interface_dofs;interface_compression_ratio='1';
            local_rank_total=I $summary.total_local_rank;local_rank_min=$rankStats.minimum;
            local_rank_max=$rankStats.maximum;port_rank_total='0';
            power_channel_count=I (Power-Channel-Count $m3Partition $caseName);
            time_step_count=I $summary.steps;dt_s=F $summary.dt_s;waveform_count='1';
            relative_l2=F $summary.space_time_relative_l2;
            max_node_error_k=F $summary.maximum_absolute_k;
            max_temperature_error_k=F $summary.maximum_temperature_error_k;
            interface_temperature_jump=F $summary.maximum_temperature_jump_rms_k;
            flux_relative_l2=F $summary.maximum_fom_rom_flux_relative_l2;
            maximum_flux_imbalance_w_m2=$na;sampled_full_residual=F $summary.maximum_full_residual;
            projected_system_residual=F $summary.maximum_reduced_residual;
            corrected_relative_l2=$na;corrected_full_residual=$na;
            offline_total_s=F ((D $summary.reference_setup_seconds)+(D $summary.local_basis_setup_seconds)+(D $summary.dynamic_schur_setup_seconds));
            snapshot_generation_s=$na;basis_generation_s=F $summary.local_basis_setup_seconds;
            reduced_projection_s=$na;schur_setup_s=F $summary.dynamic_schur_setup_seconds;
            model_load_s=$na;online_core_s=$na;full_field_reconstruction_s=$na;
            online_total_s=F $batch.projected_100_waveform_seconds;rhs_count='100';
            average_online_s_per_case=F $batch.average_online_seconds;
            baseline_name=$na;baseline_online_s_per_case=$na;speedup_vs_baseline=$na;
            offline_break_even_waveforms=$na;offline_peak_memory_bytes=I $memory.peak_working_set_bytes;
            offline_peak_memory_mib=MiB $memory.peak_working_set_bytes;
            online_peak_memory_bytes=$na;online_peak_memory_mib=$na;
            serialized_model_size_bytes=I $memory.model_bytes;basis_storage_bytes=$na;
            template_reuse_reduction_percent='0';speedup_vs_monolithic=$na;
            speedup_vs_stage1_schur=$na;speedup_vs_global_block_arnoldi=$na;
            speedup_vs_milestone6='1';commit=$m6Commit;tag='baseline-local-block-arnoldi-milestone6';
            branch='codex/feature/local-mor-schur-mainline';
            commit_message='Complete RRAM26 and Chiplet local Block Arnoldi milestone 6';
            release_status='success';ctest_passed='39';ctest_total='39';
            reproduction_script='scripts/run_local_block_arnoldi_milestone6.ps1';
            report_file='outputs/local_block_arnoldi_milestone6_report.md';
            workspace_clean_at_milestone=$na;
            accuracy_conclusion=$(if ($caseName -eq 'chiplet') {'temperature pass; flux/residual limitation'} else {'pass'});
            notes="Formal full-horizon online core is $($summary.local_online_core_seconds) s; it is not mixed with the 100-RHS one-step timing row. Monolithic ($($mono.online_seconds) s) and Global Block Arnoldi ($($global.online_seconds) s) use full-horizon definitions and are therefore not used for this row's speedup."
        }))
    }

    foreach ($caseName in @('rram26','chiplet')) {
        $summary = ($m7Summary | Where-Object case_name -eq $caseName)[0]
        $timing = ($m7Timing | Where-Object case_name -eq $caseName)[0]
        $memory = ($m7Memory | Where-Object case_name -eq $caseName)[0]
        $pure = ($m7Comparison | Where-Object {
            $_.case_name -eq $caseName -and $_.method -eq 'Milestone 7 port-reduced' })[0]
        $baseline = ($m7Comparison | Where-Object {
            $_.case_name -eq $caseName -and $_.method -eq 'Milestone 6 full-interface' })[0]
        $corrected = ($m7Comparison | Where-Object {
            $_.case_name -eq $caseName -and $_.method -eq 'Milestone 7 corrected' })[0]
        $rankStats = Rank-Stats @($m7Enrichment | Where-Object case_name -eq $caseName) 'cumulative_rank'
        $ports = @($m7Ports | Where-Object case_name -eq $caseName)
        $portTotal = I (($ports | Measure-Object selected_rank -Sum).Sum)
        $breakEven = Ceil-BreakEven $pure.offline_seconds $baseline.offline_seconds `
            $baseline.online_seconds_per_one_step_rhs $pure.online_seconds_per_one_step_rhs
        $master.Add((New-MasterRow @{
            milestone_id='M7';case_name=$caseName;steady_or_transient='transient';
            measurement_mode='100_rhs_average;one_step;pure;full_field_recovery';
            interior_reduction="$($summary.subdomains) independent Local Block Arnoldi bases";
            interface_treatment="$($ports.Count) independent physical-interface Port POD bases";
            interface_solver='dense port-reduced Dynamic Schur LU';
            basis_construction='local Block Arnoldi + scale-balanced port snapshot POD';
            corrected_mode='optional FOM residual correction';
            main_purpose='Remove full-interface bottleneck and improve flux/residual accuracy';
            global_dofs=I $summary.global_dofs;subdomain_count=I $summary.subdomains;
            physical_interface_count=I $ports.Count;full_interface_dofs=I $summary.full_interface_dofs;
            reduced_interface_dofs=I $summary.port_dimension;
            interface_compression_ratio=Ratio $summary.full_interface_dofs $summary.port_dimension;
            local_rank_total=I $summary.total_local_rank;local_rank_min=$rankStats.minimum;
            local_rank_max=$rankStats.maximum;port_rank_total=$portTotal;
            power_channel_count=I (Power-Channel-Count $m3Partition $caseName);
            time_step_count=I $timing.steps;dt_s=$(if ($caseName -eq 'rram26') {'0.01'} else {'0.1'});
            waveform_count='1';relative_l2=F $summary.space_time_relative_l2;
            max_node_error_k=F $summary.maximum_absolute_k;
            max_temperature_error_k=$na;interface_temperature_jump=$na;
            flux_relative_l2=F $summary.maximum_flux_relative_l2;
            maximum_flux_imbalance_w_m2=$na;sampled_full_residual=F $summary.maximum_full_residual;
            projected_system_residual=F $summary.maximum_port_reduced_relative_residual;
            corrected_relative_l2=F $corrected.relative_l2;
            corrected_full_residual=F $corrected.full_residual;
            offline_total_s=F $pure.offline_seconds;snapshot_generation_s=F $timing.port_snapshot_seconds;
            basis_generation_s=F ((D $timing.local_basis_setup_seconds)+(D $timing.port_basis_seconds));
            reduced_projection_s=F $timing.port_schur_assembly_seconds;
            schur_setup_s=F ((D $timing.port_schur_assembly_seconds)+(D $timing.port_schur_factor_seconds));
            model_load_s=$na;online_core_s=$na;full_field_reconstruction_s=$na;
            online_total_s=F ((D $pure.online_seconds_per_one_step_rhs)*100.0);rhs_count='100';
            average_online_s_per_case=F $pure.online_seconds_per_one_step_rhs;
            baseline_name='Milestone 6 full-interface 100-RHS one-step deployment';
            baseline_online_s_per_case=F $baseline.online_seconds_per_one_step_rhs;
            speedup_vs_baseline=Ratio $baseline.online_seconds_per_one_step_rhs $pure.online_seconds_per_one_step_rhs;
            offline_break_even_waveforms=$breakEven;
            offline_peak_memory_bytes=I $memory.peak_working_set_bytes;
            offline_peak_memory_mib=MiB $memory.peak_working_set_bytes;
            online_peak_memory_bytes=$na;online_peak_memory_mib=$na;
            serialized_model_size_bytes=I $memory.model_bytes;basis_storage_bytes=$na;
            template_reuse_reduction_percent='0';speedup_vs_monolithic=$na;
            speedup_vs_stage1_schur=$na;speedup_vs_global_block_arnoldi=$na;
            speedup_vs_milestone6=Ratio $baseline.online_seconds_per_one_step_rhs $pure.online_seconds_per_one_step_rhs;
            commit=$m7Commit;tag='baseline-local-port-reduced-dynamic-schur-milestone7';
            branch='codex/feature/local-port-reduced-dynamic-schur';
            commit_message='Complete local port-reduced Dynamic Schur milestone 7';
            release_status='success';ctest_passed='51';ctest_total='51';
            reproduction_script='scripts/run_local_port_reduced_dynamic_schur.ps1';
            report_file='outputs/local_port_rom_report.md';workspace_clean_at_milestone='yes';
            accuracy_conclusion='pass';
            notes="Full-interface Local-ROM pilot is offline-only ($($timing.port_local_full_interface_pilot_seconds) s). Formal full-horizon online core is $($timing.online_core_seconds) s and is not mixed with this 100-RHS one-step row."
        }))
    }

    $methodDefinitions = @(
        [pscustomobject]@{Milestone='M1';Case='two-cube';SteadyTransient='Steady';InteriorReduction='Local POD/RB';InterfaceTreatment='Full SIPG interface';InterfaceSolver='Sparse PARDISO';BasisConstruction='Per-subdomain snapshots';CorrectedMode='Optional exact Schur';MainPurpose='Establish local steady ROM';Commit=$m1Commit;CTest='22/22'},
        [pscustomobject]@{Milestone='M2';Case='ten-cube';SteadyTransient='Steady';InteriorReduction='Local POD/RB';InterfaceTreatment='Full SIPG interface';InterfaceSolver='Sparse PARDISO';BasisConstruction='10 independent POD bases';CorrectedMode='Optional exact Schur';MainPurpose='Scale steady path';Commit=$m2Commit;CTest='32/32'},
        [pscustomobject]@{Milestone='M3';Case='RRAM26 / Chiplet';SteadyTransient='Steady';InteriorReduction='Local POD/RB';InterfaceTreatment='Full SIPG interface';InterfaceSolver='Matrix-free FGMRES';BasisConstruction='Physical-subdomain POD';CorrectedMode='Residual-gated exact Schur';MainPurpose='Real-model validation';Commit=$m3Commit;CTest='33/33'},
        [pscustomobject]@{Milestone='M4';Case='two-cube';SteadyTransient='Transient';InteriorReduction='Local Block Arnoldi';InterfaceTreatment='Full SIPG interface';InterfaceSolver='Fixed-dt Dynamic Schur';BasisConstruction='Local Krylov moments';CorrectedMode='Regression only';MainPurpose='Establish transient path';Commit=$m4Commit;CTest='36/36'},
        [pscustomobject]@{Milestone='M5';Case='ten-cube';SteadyTransient='Transient';InteriorReduction='Local Block Arnoldi';InterfaceTreatment='Full SIPG interface';InterfaceSolver='Fixed-dt Dynamic Schur';BasisConstruction='10 independent Krylov bases';CorrectedMode='Regression only';MainPurpose='Scale transient path';Commit=$m5Commit;CTest='38/38'},
        [pscustomobject]@{Milestone='M6';Case='RRAM26 / Chiplet';SteadyTransient='Transient';InteriorReduction='Local Block Arnoldi';InterfaceTreatment='Full SIPG interface';InterfaceSolver='FGMRES + 1-ring proxy';BasisConstruction='Physical-subdomain Krylov';CorrectedMode='Not acceptance path';MainPurpose='Expose interface bottleneck';Commit=$m6Commit;CTest='39/39'},
        [pscustomobject]@{Milestone='M7';Case='RRAM26 / Chiplet';SteadyTransient='Transient';InteriorReduction='Local Block Arnoldi';InterfaceTreatment='Independent Port POD';InterfaceSolver='Port-Reduced Dynamic Schur';BasisConstruction='Local Krylov + port POD + enrichment';CorrectedMode='Optional FOM residual';MainPurpose='Final scalable local mainline';Commit=$m7Commit;CTest='51/51'})

    $methodEvolution = $methodDefinitions | Select-Object `
        @{n='milestone';e={$_.Milestone}},@{n='case';e={$_.Case}},
        @{n='steady_or_transient';e={$_.SteadyTransient}},
        @{n='interior_reduction';e={$_.InteriorReduction}},
        @{n='interface_treatment';e={$_.InterfaceTreatment}},
        @{n='interface_solver';e={$_.InterfaceSolver}},
        @{n='basis_construction';e={$_.BasisConstruction}},
        @{n='corrected_mode';e={$_.CorrectedMode}},@{n='main_purpose';e={$_.MainPurpose}},
        @{n='commit';e={$_.Commit}},@{n='ctest';e={$_.CTest}}

    $modelScale = $master | Select-Object milestone_id,case_name,global_dofs,subdomain_count,
        physical_interface_count,full_interface_dofs,reduced_interface_dofs,
        interface_compression_ratio,local_rank_total,local_rank_min,local_rank_max,
        port_rank_total,power_channel_count,serialized_model_size_bytes,
        offline_peak_memory_bytes,offline_peak_memory_mib

    $accuracy = $master | Select-Object milestone_id,case_name,relative_l2,max_node_error_k,
        max_temperature_error_k,max_temperature_curve_rmse_k,interface_temperature_jump,
        flux_relative_l2,maximum_flux_imbalance_w_m2,sampled_full_residual,
        projected_system_residual,corrected_relative_l2,corrected_full_residual,
        accuracy_conclusion

    $timing = $master | Select-Object milestone_id,case_name,measurement_mode,offline_total_s,
        snapshot_generation_s,basis_generation_s,reduced_projection_s,schur_setup_s,
        model_load_s,online_core_s,full_field_reconstruction_s,online_total_s,rhs_count,
        average_online_s_per_case,baseline_name,baseline_online_s_per_case,
        speedup_vs_baseline,offline_break_even_waveforms

    $memory = $master | Select-Object milestone_id,case_name,offline_peak_memory_bytes,
        offline_peak_memory_mib,online_peak_memory_bytes,online_peak_memory_mib,
        serialized_model_size_bytes,basis_storage_bytes,template_reuse_reduction_percent

    $speedup = $master | Select-Object milestone_id,case_name,measurement_mode,
        speedup_vs_monolithic,speedup_vs_stage1_schur,speedup_vs_global_block_arnoldi,
        speedup_vs_milestone6,offline_break_even_waveforms,baseline_name,
        baseline_online_s_per_case,average_online_s_per_case

    $version = $master | Group-Object milestone_id | ForEach-Object {
        $row = $_.Group | Select-Object -First 1
        [pscustomobject]@{
            milestone=$row.milestone_id;commit=$row.commit;tag=$row.tag;branch=$row.branch
            commit_message=$row.commit_message;release_build=$row.release_status
            ctest_passed=$row.ctest_passed;ctest_total=$row.ctest_total
            reproduction_script=$row.reproduction_script;report_file=$row.report_file
            workspace_clean=$row.workspace_clean_at_milestone
        }
    }

    $sourceCatalog = @{
        'M1|two-cube'=@{commit=$m1Commit;summary='git:852c8f3:outputs/local_rom_schur_summary.csv';accuracy='git:852c8f3:outputs/local_rom_accuracy_by_case.csv';timing='git:852c8f3:outputs/local_rom_offline_timing.csv';memory='git:852c8f3:outputs/local_rom_memory.csv';rank='git:852c8f3:outputs/local_rom_rank_by_subdomain.csv';report='git:852c8f3:outputs/local_rom_schur_report.md';comparison='git:852c8f3:outputs/local_rom_vs_global_rom.csv'}
        'M2|ten-cube'=@{commit=$m2Commit;summary='outputs/local_rom_schur_summary.csv';accuracy='outputs/local_rom_accuracy_by_case.csv';timing='outputs/local_rom_multi_rhs_timing.csv';memory='outputs/local_rom_memory.csv';rank='outputs/local_rom_rank_by_subdomain.csv';report='outputs/local_rom_schur_report.md';comparison='outputs/local_rom_vs_global_rom.csv'}
        'M3|rram26'=@{commit=$m3Commit;summary='outputs/local_rom_milestone3_summary.csv';accuracy='outputs/local_rom_milestone3_summary.csv';timing='outputs/local_rom_milestone3_multi_rhs.csv';memory='outputs/local_rom_milestone3_memory.csv';rank='outputs/local_rom_milestone3_rank_by_subdomain.csv';report='outputs/local_rom_milestone3_report.md';comparison='outputs/local_rom_milestone3_vs_global.csv'}
        'M3|chiplet'=@{commit=$m3Commit;summary='outputs/local_rom_milestone3_summary.csv';accuracy='outputs/local_rom_milestone3_summary.csv';timing='outputs/local_rom_milestone3_multi_rhs.csv';memory='outputs/local_rom_milestone3_memory.csv';rank='outputs/local_rom_milestone3_rank_by_subdomain.csv';report='outputs/local_rom_milestone3_report.md';comparison='outputs/local_rom_milestone3_vs_global.csv'}
        'M4|two-cube'=@{commit=$m4Commit;summary='outputs/local_block_arnoldi_milestone4_summary.csv';accuracy='outputs/local_dynamic_schur_milestone4.csv';timing='outputs/local_block_arnoldi_milestone4_summary.csv';memory='outputs/local_dynamic_schur_milestone4.csv';rank='outputs/local_block_arnoldi_milestone4_rank.csv';report='outputs/local_block_arnoldi_milestone4_report.md';comparison='outputs/local_block_arnoldi_milestone4_summary.csv'}
        'M5|ten-cube'=@{commit=$m5Commit;summary='outputs/local_block_arnoldi_milestone5_summary.csv';accuracy='outputs/local_dynamic_schur_milestone5.csv';timing='outputs/local_block_arnoldi_milestone5_summary.csv';memory='outputs/local_block_arnoldi_milestone5_summary.csv';rank='outputs/local_block_arnoldi_milestone5_rank_by_subdomain.csv';report='outputs/local_block_arnoldi_milestone5_report.md';comparison='outputs/local_block_arnoldi_milestone5_summary.csv'}
        'M6|rram26'=@{commit=$m6Commit;summary='outputs/local_block_arnoldi_milestone6_summary.csv';accuracy='outputs/local_block_arnoldi_milestone6_summary.csv';timing='outputs/local_transient_break_even.csv';memory='outputs/local_transient_memory.csv';rank='outputs/local_block_arnoldi_milestone6_rank_by_subdomain.csv';report='outputs/local_block_arnoldi_milestone6_report.md';comparison='outputs/local_transient_vs_global_block_arnoldi.csv'}
        'M6|chiplet'=@{commit=$m6Commit;summary='outputs/local_block_arnoldi_milestone6_summary.csv';accuracy='outputs/local_block_arnoldi_milestone6_summary.csv';timing='outputs/local_transient_break_even.csv';memory='outputs/local_transient_memory.csv';rank='outputs/local_block_arnoldi_milestone6_rank_by_subdomain.csv';report='outputs/local_block_arnoldi_milestone6_report.md';comparison='outputs/local_transient_vs_global_block_arnoldi.csv'}
        'M7|rram26'=@{commit=$m7Commit;summary='outputs/local_port_rom_summary.csv';accuracy='outputs/local_port_rom_summary.csv';timing='outputs/local_port_vs_full_interface.csv';memory='outputs/local_port_memory.csv';rank='outputs/local_port_enrichment_history.csv';report='outputs/local_port_rom_report.md';comparison='outputs/local_port_vs_full_interface.csv'}
        'M7|chiplet'=@{commit=$m7Commit;summary='outputs/local_port_rom_summary.csv';accuracy='outputs/local_port_rom_summary.csv';timing='outputs/local_port_vs_full_interface.csv';memory='outputs/local_port_memory.csv';rank='outputs/local_port_enrichment_history.csv';report='outputs/local_port_rom_report.md';comparison='outputs/local_port_vs_full_interface.csv'}
    }

    # Field-level provenance overrides.  The generic catalog above chooses a
    # source family; these maps name the exact archived CSV column (or formula)
    # used to create each important numeric value.
    $sourceOverrides = @{}
    function Add-SourceMap([string]$RowKey, [string]$File, [hashtable]$Columns) {
        foreach ($metric in $Columns.Keys) {
            $sourceOverrides["$RowKey|$metric"] = [pscustomobject]@{
                file = $File
                column = [string]$Columns[$metric]
            }
        }
    }

    Add-SourceMap 'M1|two-cube' 'git:852c8f3:outputs/local_rom_schur_summary.csv' @{
        global_dofs='global_dofs';subdomain_count='subdomains';full_interface_dofs='interface_dofs';
        reduced_interface_dofs='interface_dofs';interface_compression_ratio='interface_dofs / interface_dofs';
        local_rank_total='total_local_rank';relative_l2='relative_l2'
    }
    Add-SourceMap 'M1|two-cube' 'git:852c8f3:outputs/local_rom_accuracy_by_case.csv' @{
        max_node_error_k='max(max_node_error_k)';max_temperature_error_k='max(max_temperature_error_k)';
        interface_temperature_jump='max(interface_temperature_jump_rms_k)';
        maximum_flux_imbalance_w_m2='max(maximum_flux_imbalance_w_m2)';
        sampled_full_residual='max(global_relative_residual)';
        projected_system_residual='max(interface_relative_residual)'
    }
    Add-SourceMap 'M1|two-cube' 'git:852c8f3:outputs/local_rom_rank_by_subdomain.csv' @{
        local_rank_min='min(selected_rank)';local_rank_max='max(selected_rank)'
    }
    Add-SourceMap 'M1|two-cube' 'git:852c8f3:outputs/local_rom_offline_timing.csv' @{
        offline_total_s='total_offline_seconds';snapshot_generation_s='snapshot_solve_seconds';
        basis_generation_s='local_basis_seconds';reduced_projection_s='local_projection_seconds';
        schur_setup_s='reduced_schur_construction_seconds + reduced_factorization_seconds'
    }
    Add-SourceMap 'M1|two-cube' 'git:852c8f3:outputs/local_rom_online_timing.csv' @{
        online_core_s='total_seconds';full_field_reconstruction_s='full_field_reconstruction_seconds';
        online_total_s='total_seconds';average_online_s_per_case='total_seconds'
    }
    Add-SourceMap 'M1|two-cube' 'git:852c8f3:outputs/local_rom_memory.csv' @{
        offline_peak_memory_bytes='peak_working_set_bytes';offline_peak_memory_mib='peak_working_set_bytes / 1048576';
        serialized_model_size_bytes='model_bytes'
    }
    Add-SourceMap 'M1|two-cube' 'git:852c8f3:outputs/local_rom_vs_global_rom.csv' @{
        corrected_relative_l2='relative_l2 where method=Local-POD-Schur-ROM-Corrected';
        baseline_online_s_per_case='solve_seconds where method=DDM-Schur-FGMRES';
        speedup_vs_baseline='DDM-Schur-FGMRES solve_seconds / Local-ROM online';
        speedup_vs_monolithic='Global-PARDISO-SPD-Direct solve_seconds / Local-ROM online';
        speedup_vs_stage1_schur='DDM-Schur-FGMRES solve_seconds / Local-ROM online'
    }

    Add-SourceMap 'M2|ten-cube' 'outputs/local_rom_schur_summary.csv' @{
        global_dofs='global_dofs';subdomain_count='subdomains';full_interface_dofs='interface_dofs';
        reduced_interface_dofs='interface_dofs';interface_compression_ratio='interface_dofs / interface_dofs';
        local_rank_total='total_local_rank';relative_l2='relative_l2';corrected_full_residual='final_true_residual'
    }
    Add-SourceMap 'M2|ten-cube' 'outputs/local_rom_accuracy_by_case.csv' @{
        max_node_error_k='max(max_node_error_k)';max_temperature_error_k='max(max_temperature_error_k)';
        interface_temperature_jump='max(interface_temperature_jump_rms_k)';
        maximum_flux_imbalance_w_m2='max(maximum_flux_imbalance_w_m2)';
        sampled_full_residual='max(global_relative_residual)';
        projected_system_residual='max(interface_relative_residual)'
    }
    Add-SourceMap 'M2|ten-cube' 'outputs/local_rom_rank_by_subdomain.csv' @{
        local_rank_min='min(selected_rank)';local_rank_max='max(selected_rank)'
    }
    Add-SourceMap 'M2|ten-cube' 'outputs/local_rom_offline_timing.csv' @{
        offline_total_s='total_offline_seconds where scenario=generation';
        snapshot_generation_s='snapshot_solve_seconds where scenario=generation';
        basis_generation_s='local_basis_seconds where scenario=generation';
        reduced_projection_s='local_projection_seconds where scenario=generation';
        schur_setup_s='reduced_schur_construction_seconds + reduced_factorization_seconds';
        model_load_s='model_load_seconds where scenario=serialized_model_reload';
        online_peak_memory_bytes='peak_working_set_bytes where scenario=serialized_model_reload';
        online_peak_memory_mib='peak_working_set_bytes / 1048576 where scenario=serialized_model_reload'
    }
    Add-SourceMap 'M2|ten-cube' 'outputs/local_rom_multi_rhs_timing.csv' @{
        online_core_s='online_total_seconds where rhs_count=100';
        full_field_reconstruction_s='full_field_reconstruction_total_seconds where rhs_count=100';
        online_total_s='online_total_seconds where rhs_count=100';rhs_count='rhs_count';
        average_online_s_per_case='average_online_seconds where rhs_count=100'
    }
    Add-SourceMap 'M2|ten-cube' 'outputs/local_rom_memory.csv' @{
        offline_peak_memory_bytes='peak_working_set_bytes';offline_peak_memory_mib='peak_working_set_bytes / 1048576';
        serialized_model_size_bytes='model_bytes';basis_storage_bytes='basis_storage_with_reuse_bytes';
        template_reuse_reduction_percent='calculated from basis_storage_without_reuse_bytes and basis_storage_with_reuse_bytes'
    }
    Add-SourceMap 'M2|ten-cube' 'outputs/local_rom_vs_global_rom.csv' @{
        corrected_relative_l2='relative_l2 where method=Local-POD-Schur-ROM-Corrected'
    }

    foreach ($caseName in @('rram26','chiplet')) {
        $rowKey = "M3|$caseName"
        Add-SourceMap $rowKey 'outputs/local_rom_milestone3_summary.csv' @{
            global_dofs='global_dofs';subdomain_count='subdomains';full_interface_dofs='interface_dofs';
            reduced_interface_dofs='interface_dofs';interface_compression_ratio='interface_dofs / interface_dofs';
            local_rank_total='total_local_rank';relative_l2='relative_l2 where mode=pure';
            max_node_error_k='max_node_error_k where mode=pure';max_temperature_error_k='max_temperature_error_k where mode=pure';
            interface_temperature_jump='interface_temperature_jump_rms_k where mode=pure';
            sampled_full_residual='initial_true_residual where mode=pure';
            projected_system_residual='interface_fgmres_true_residual where mode=pure';
            corrected_relative_l2='relative_l2 where mode=corrected';
            corrected_full_residual='final_true_residual where mode=corrected';offline_total_s='offline_seconds where mode=pure'
        }
        Add-SourceMap $rowKey 'outputs/local_rom_milestone3_timing.csv' @{
            snapshot_generation_s='snapshot_solve_seconds';basis_generation_s='local_basis_seconds';
            reduced_projection_s='local_projection_seconds';schur_setup_s='reduced_schur_construction_seconds'
        }
        Add-SourceMap $rowKey 'outputs/local_rom_milestone3_multi_rhs.csv' @{
            online_core_s='online_total_seconds where rhs_count=100';
            full_field_reconstruction_s='full_field_reconstruction_total_seconds where rhs_count=100';
            online_total_s='online_total_seconds where rhs_count=100';rhs_count='rhs_count';
            average_online_s_per_case='average_online_seconds where rhs_count=100'
        }
        Add-SourceMap $rowKey 'outputs/local_rom_milestone3_memory.csv' @{
            offline_peak_memory_bytes='peak_working_set_bytes';offline_peak_memory_mib='peak_working_set_bytes / 1048576';
            serialized_model_size_bytes='model_bytes';basis_storage_bytes='basis_storage_with_reuse_bytes'
        }
        Add-SourceMap $rowKey 'outputs/local_rom_milestone3_rank_by_subdomain.csv' @{
            local_rank_min='min(selected_rank)';local_rank_max='max(selected_rank)'
        }
        Add-SourceMap $rowKey 'outputs/local_rom_milestone3_partition_definition.csv' @{
            physical_interface_count='unique pairs from neighbor_subdomains';power_channel_count='sum(power_channel_count)'
        }
    }

    Add-SourceMap 'M4|two-cube' 'outputs/local_dynamic_schur_milestone4.csv' @{
        global_dofs='global_dofs';subdomain_count='subdomains';full_interface_dofs='full_interface_dofs';
        reduced_interface_dofs='full_interface_dofs';interface_compression_ratio='full_interface_dofs / full_interface_dofs';
        time_step_count='steps';dt_s='dt_s';waveform_count='count(unique waveform)';
        relative_l2='max(space_time_relative_l2)';max_node_error_k='max(maximum_absolute_k)';
        max_temperature_error_k='max(maximum_temperature_error_k)';
        interface_temperature_jump='max(maximum_temperature_jump_rms_k)';
        sampled_full_residual='max(maximum_full_residual)';projected_system_residual='max(maximum_reduced_residual)';
        basis_generation_s='local_basis_setup_seconds';schur_setup_s='dynamic_schur_setup_seconds';
        full_field_reconstruction_s='local_recovery_seconds'
    }
    Add-SourceMap 'M4|two-cube' 'outputs/local_block_arnoldi_milestone4_summary.csv' @{
        local_rank_total='rank where method=Local Block Arnoldi + Dynamic Schur';
        offline_total_s='offline_seconds where method=Local Block Arnoldi + Dynamic Schur';
        online_core_s='online_seconds where method=Local Block Arnoldi + Dynamic Schur';
        online_total_s='online_seconds where method=Local Block Arnoldi + Dynamic Schur';
        average_online_s_per_case='online_seconds where method=Local Block Arnoldi + Dynamic Schur';
        baseline_online_s_per_case='online_seconds where method=Transient monolithic PARDISO';
        speedup_vs_baseline='monolithic online_seconds / local online_seconds';
        speedup_vs_monolithic='monolithic online_seconds / local online_seconds'
    }
    Add-SourceMap 'M4|two-cube' 'outputs/local_block_arnoldi_milestone4_rank.csv' @{
        local_rank_min='min(final cumulative_rank by subdomain)';local_rank_max='max(final cumulative_rank by subdomain)'
    }

    Add-SourceMap 'M5|ten-cube' 'outputs/local_dynamic_schur_milestone5.csv' @{
        global_dofs='global_dofs';subdomain_count='subdomains';full_interface_dofs='full_interface_dofs';
        reduced_interface_dofs='full_interface_dofs';interface_compression_ratio='full_interface_dofs / full_interface_dofs';
        time_step_count='steps';dt_s='dt_s';waveform_count='count(unique waveform)';
        relative_l2='max(space_time_relative_l2)';max_node_error_k='max(maximum_absolute_k)';
        max_temperature_error_k='max(maximum_temperature_error_k)';
        interface_temperature_jump='max(maximum_temperature_jump_rms_k)';
        sampled_full_residual='max(maximum_full_residual)';projected_system_residual='max(maximum_reduced_residual)';
        basis_generation_s='local_basis_setup_seconds';schur_setup_s='dynamic_schur_setup_seconds'
    }
    Add-SourceMap 'M5|ten-cube' 'outputs/local_block_arnoldi_milestone5_summary.csv' @{
        local_rank_total='rank where method=Local Block Arnoldi + Dynamic Schur';
        offline_total_s='offline_seconds where method=Local Block Arnoldi + Dynamic Schur';
        online_core_s='online_seconds where method=Local Block Arnoldi + Dynamic Schur';
        online_total_s='online_seconds where method=Local Block Arnoldi + Dynamic Schur';
        average_online_s_per_case='online_seconds where method=Local Block Arnoldi + Dynamic Schur';
        baseline_online_s_per_case='online_seconds where method=Transient monolithic PARDISO';
        offline_peak_memory_bytes='peak_working_set_bytes';offline_peak_memory_mib='peak_working_set_bytes / 1048576';
        serialized_model_size_bytes='model_bytes';speedup_vs_baseline='monolithic online_seconds / local online_seconds';
        speedup_vs_monolithic='monolithic online_seconds / local online_seconds'
    }
    Add-SourceMap 'M5|ten-cube' 'outputs/local_block_arnoldi_milestone5_rank_by_subdomain.csv' @{
        local_rank_min='min(final cumulative_rank by subdomain)';local_rank_max='max(final cumulative_rank by subdomain)'
    }

    foreach ($caseName in @('rram26','chiplet')) {
        $rowKey = "M6|$caseName"
        Add-SourceMap $rowKey 'outputs/local_block_arnoldi_milestone6_summary.csv' @{
            global_dofs='global_dofs';subdomain_count='subdomains';full_interface_dofs='full_interface_dofs';
            reduced_interface_dofs='full_interface_dofs';interface_compression_ratio='full_interface_dofs / full_interface_dofs';
            local_rank_total='total_local_rank';time_step_count='steps';dt_s='dt_s';relative_l2='space_time_relative_l2';
            max_node_error_k='maximum_absolute_k';max_temperature_error_k='maximum_temperature_error_k';
            interface_temperature_jump='maximum_temperature_jump_rms_k';flux_relative_l2='maximum_fom_rom_flux_relative_l2';
            sampled_full_residual='maximum_full_residual';projected_system_residual='maximum_reduced_residual';
            offline_total_s='reference_setup_seconds + local_basis_setup_seconds + dynamic_schur_setup_seconds';
            basis_generation_s='local_basis_setup_seconds';schur_setup_s='dynamic_schur_setup_seconds'
        }
        Add-SourceMap $rowKey 'outputs/local_transient_break_even.csv' @{
            online_total_s='projected_100_waveform_seconds where measurement=actual_100_one_step_full_recovery_waveforms';
            rhs_count='verified_waveforms';average_online_s_per_case='average_online_seconds where measurement=actual_100_one_step_full_recovery_waveforms'
        }
        Add-SourceMap $rowKey 'outputs/local_transient_memory.csv' @{
            offline_peak_memory_bytes='peak_working_set_bytes';offline_peak_memory_mib='peak_working_set_bytes / 1048576';
            serialized_model_size_bytes='model_bytes'
        }
        Add-SourceMap $rowKey 'outputs/local_block_arnoldi_milestone6_rank_by_subdomain.csv' @{
            local_rank_min='min(final cumulative_rank by subdomain)';local_rank_max='max(final cumulative_rank by subdomain)'
        }
        Add-SourceMap $rowKey 'outputs/local_rom_milestone3_partition_definition.csv' @{
            physical_interface_count='unique pairs from neighbor_subdomains';power_channel_count='sum(power_channel_count)'
        }
    }

    foreach ($caseName in @('rram26','chiplet')) {
        $rowKey = "M7|$caseName"
        Add-SourceMap $rowKey 'outputs/local_port_rom_summary.csv' @{
            global_dofs='global_dofs';subdomain_count='subdomains';full_interface_dofs='full_interface_dofs';
            reduced_interface_dofs='port_dimension';interface_compression_ratio='full_interface_dofs / port_dimension';
            local_rank_total='total_local_rank';relative_l2='space_time_relative_l2';
            max_node_error_k='maximum_absolute_k';flux_relative_l2='maximum_flux_relative_l2';
            sampled_full_residual='maximum_full_residual';projected_system_residual='maximum_port_reduced_relative_residual'
        }
        Add-SourceMap $rowKey 'outputs/local_port_rank_by_interface.csv' @{
            physical_interface_count='count(interface_id)';port_rank_total='sum(selected_rank)'
        }
        Add-SourceMap $rowKey 'outputs/local_port_enrichment_history.csv' @{
            local_rank_min='min(final cumulative_rank by subdomain)';local_rank_max='max(final cumulative_rank by subdomain)'
        }
        Add-SourceMap $rowKey 'outputs/local_port_dynamic_schur_timing.csv' @{
            time_step_count='steps';snapshot_generation_s='port_snapshot_seconds';
            basis_generation_s='local_basis_setup_seconds + port_basis_seconds';
            reduced_projection_s='port_schur_assembly_seconds';
            schur_setup_s='port_schur_assembly_seconds + port_schur_factor_seconds'
        }
        Add-SourceMap $rowKey 'outputs/local_port_vs_full_interface.csv' @{
            corrected_relative_l2='relative_l2 where method=Milestone 7 corrected';
            corrected_full_residual='full_residual where method=Milestone 7 corrected';
            offline_total_s='offline_seconds where method=Milestone 7 port-reduced';
            online_total_s='100 * online_seconds_per_one_step_rhs where method=Milestone 7 port-reduced';
            rhs_count='100 deployment queries';
            average_online_s_per_case='online_seconds_per_one_step_rhs where method=Milestone 7 port-reduced';
            baseline_online_s_per_case='online_seconds_per_one_step_rhs where method=Milestone 6 full-interface';
            speedup_vs_baseline='M6 online_seconds_per_one_step_rhs / M7 online_seconds_per_one_step_rhs';
            speedup_vs_milestone6='M6 online_seconds_per_one_step_rhs / M7 online_seconds_per_one_step_rhs';
            offline_break_even_waveforms='calculated from offline_seconds and online_seconds_per_one_step_rhs'
        }
        Add-SourceMap $rowKey 'outputs/local_port_memory.csv' @{
            offline_peak_memory_bytes='peak_working_set_bytes';offline_peak_memory_mib='peak_working_set_bytes / 1048576';
            serialized_model_size_bytes='model_bytes'
        }
        Add-SourceMap $rowKey 'outputs/local_rom_milestone3_partition_definition.csv' @{
            power_channel_count='sum(power_channel_count)'
        }
    }

    $accuracyColumns = @('relative_l2','max_node_error_k','max_temperature_error_k',
        'max_temperature_curve_rmse_k','interface_temperature_jump','flux_relative_l2',
        'maximum_flux_imbalance_w_m2','sampled_full_residual','projected_system_residual',
        'corrected_relative_l2','corrected_full_residual','accuracy_conclusion')
    $timingColumns = @('offline_total_s','snapshot_generation_s','basis_generation_s',
        'reduced_projection_s','schur_setup_s','model_load_s','online_core_s',
        'full_field_reconstruction_s','online_total_s','rhs_count','average_online_s_per_case')
    $memoryColumns = @('offline_peak_memory_bytes','offline_peak_memory_mib',
        'online_peak_memory_bytes','online_peak_memory_mib','serialized_model_size_bytes',
        'basis_storage_bytes','template_reuse_reduction_percent')
    $rankColumns = @('local_rank_total','local_rank_min','local_rank_max','port_rank_total')
    $comparisonColumns = @('baseline_name','baseline_online_s_per_case','speedup_vs_baseline',
        'offline_break_even_waveforms','speedup_vs_monolithic','speedup_vs_stage1_schur',
        'speedup_vs_global_block_arnoldi','speedup_vs_milestone6')
    $unitMap = @{
        global_dofs='DOF';subdomain_count='count';physical_interface_count='count';
        full_interface_dofs='DOF';reduced_interface_dofs='DOF';interface_compression_ratio='ratio';
        local_rank_total='count';local_rank_min='count';local_rank_max='count';port_rank_total='count';
        power_channel_count='count';time_step_count='count';dt_s='s';waveform_count='count';
        relative_l2='relative';max_node_error_k='K';max_temperature_error_k='K';
        max_temperature_curve_rmse_k='K';interface_temperature_jump='K';
        flux_relative_l2='relative';maximum_flux_imbalance_w_m2='W/m^2';
        sampled_full_residual='relative';projected_system_residual='relative';
        corrected_relative_l2='relative';corrected_full_residual='relative';
        offline_total_s='s';snapshot_generation_s='s';basis_generation_s='s';
        reduced_projection_s='s';schur_setup_s='s';model_load_s='s';online_core_s='s';
        full_field_reconstruction_s='s';online_total_s='s';rhs_count='count';
        average_online_s_per_case='s/case';baseline_online_s_per_case='s/case';
        speedup_vs_baseline='x';offline_break_even_waveforms='waveforms';
        offline_peak_memory_bytes='bytes';offline_peak_memory_mib='MiB';
        online_peak_memory_bytes='bytes';online_peak_memory_mib='MiB';
        serialized_model_size_bytes='bytes';basis_storage_bytes='bytes';
        template_reuse_reduction_percent='percent';speedup_vs_monolithic='x';
        speedup_vs_stage1_schur='x';speedup_vs_global_block_arnoldi='x';
        speedup_vs_milestone6='x'
    }
    $provenance = foreach ($row in $master) {
        $catalog = $sourceCatalog["$($row.milestone_id)|$($row.case_name)"]
        foreach ($property in $row.PSObject.Properties) {
            $metric = $property.Name; $value = [string]$property.Value
            $sourceKind = 'summary'; $notes = ''
            if ($accuracyColumns -contains $metric) { $sourceKind = 'accuracy' }
            elseif ($timingColumns -contains $metric) { $sourceKind = 'timing' }
            elseif ($memoryColumns -contains $metric) { $sourceKind = 'memory' }
            elseif ($rankColumns -contains $metric) { $sourceKind = 'rank' }
            elseif ($comparisonColumns -contains $metric) { $sourceKind = 'comparison' }
            elseif ($metric -in @('commit','tag','branch','commit_message','release_status',
                    'ctest_passed','ctest_total','workspace_clean_at_milestone')) {
                $sourceKind = 'report'
            }
            elseif ($metric -in @('interior_reduction','interface_treatment','interface_solver',
                    'basis_construction','corrected_mode','main_purpose','notes')) {
                $sourceKind = 'report'
            }
            $sourceFile = $catalog[$sourceKind]
            $sourceColumn = $metric
            $overrideKey = "$($row.milestone_id)|$($row.case_name)|$metric"
            if ($sourceOverrides.ContainsKey($overrideKey)) {
                $sourceFile = $sourceOverrides[$overrideKey].file
                $sourceColumn = $sourceOverrides[$overrideKey].column
            } elseif ($metric -in @('commit','tag','branch','commit_message')) {
                $sourceFile = 'git'
                $sourceColumn = $(switch ($metric) {
                    'commit' {'git log / rev-parse'}
                    'tag' {'git tag --points-at'}
                    'branch' {'git branch --contains'}
                    'commit_message' {'git log --format=%s'}
                })
            } elseif ($sourceKind -eq 'report') {
                $sourceColumn = 'reported method/version metadata'
            }
            if ($metric -in @('interface_compression_ratio','offline_peak_memory_mib',
                    'online_peak_memory_mib','speedup_vs_baseline','speedup_vs_monolithic',
                    'speedup_vs_stage1_schur','speedup_vs_global_block_arnoldi',
                    'speedup_vs_milestone6','offline_break_even_waveforms')) {
                if (-not $sourceOverrides.ContainsKey($overrideKey)) {
                    $sourceColumn = 'calculated from archived source columns'
                }
                $notes = 'Recomputed by build_milestone1_7_final_summary.ps1.'
            }
            if ($value -eq $na) {
                $sourceFile = $na; $sourceColumn = $na
                $notes = 'not recorded in the archived milestone outputs'
            }
            [pscustomobject]@{
                milestone=$row.milestone_id;case=$row.case_name;metric=$metric;value=$value
                unit=$(if ($unitMap.ContainsKey($metric)) {$unitMap[$metric]} else {'text'})
                source_file=$sourceFile;source_column=$sourceColumn
                source_commit=$catalog.commit;notes=$notes
            }
        }
    }

    $csvOutputs = [ordered]@{
        'milestone1_7_method_evolution.csv'=$methodEvolution
        'milestone1_7_model_scale.csv'=$modelScale
        'milestone1_7_accuracy.csv'=$accuracy
        'milestone1_7_timing.csv'=$timing
        'milestone1_7_memory.csv'=$memory
        'milestone1_7_speedup.csv'=$speedup
        'milestone1_7_version_status.csv'=$version
        'milestone1_7_master_summary.csv'=$master
        'milestone1_7_data_provenance.csv'=$provenance
    }
    foreach ($entry in $csvOutputs.GetEnumerator()) {
        @($entry.Value) | Export-Csv -NoTypeInformation -Encoding UTF8 `
            (Join-Path $outputs $entry.Key)
    }

    function Tex-Escape([string]$Text) {
        if ($null -eq $Text) { return '' }
        return $Text.Replace('\','\textbackslash{}').Replace('_','\_').Replace('%','\%').Replace('&','\&')
    }
    function Tex-Number($Value) {
        if ($Value -eq $na) { return '\text{N/A}' }
        $number = D $Value
        if ([double]::IsNaN($number)) { return Tex-Escape ([string]$Value) }
        if ($number -eq 0.0) { return '0' }
        $exponent = [Math]::Floor([Math]::Log10([Math]::Abs($number)))
        if ([Math]::Abs($exponent) -ge 3) {
            $mantissa = $number / [Math]::Pow(10.0,$exponent)
            return ('{0:F3}\times 10^{{{1}}}' -f $mantissa,[int]$exponent)
        }
        return $number.ToString('0.###',$invariant)
    }
    function Write-Utf8([string]$Path, [string[]]$Lines) {
        $Lines | Set-Content -LiteralPath $Path -Encoding UTF8
    }

    $methodTex = @(
        '\begin{table}[t]','\centering','\caption{Evolution of the proposed local reduced-order domain-decomposition framework}',
        '\label{tab:milestone_evolution}','\small','\begin{tabular}{cllll}','\toprule',
        'Stage & Case & Interior model & Interface model & Purpose \\','\midrule')
    foreach ($row in $methodDefinitions) {
        $methodTex += ('{0} & {1} & {2} & {3} & {4} \\' -f
            (Tex-Escape $row.Milestone),(Tex-Escape $row.Case),(Tex-Escape $row.InteriorReduction),
            (Tex-Escape $row.InterfaceTreatment),(Tex-Escape $row.MainPurpose))
    }
    $methodTex += @('\bottomrule','\end{tabular}',
        '\begin{minipage}{0.96\linewidth}\footnotesize Global POD, Global Reduced Schur, and Global Block Arnoldi are benchmark paths; the proposed mainline is local interior reduction plus physical-interface Port POD and Reduced Dynamic Schur.\end{minipage}',
        '\end{table}')
    Write-Utf8 (Join-Path $outputs 'milestone1_7_method_evolution_table.tex') $methodTex

    $accuracyTex = @(
        '\begin{table}[t]','\centering','\caption{Accuracy of the local ROM and port-reduced Dynamic Schur models}',
        '\label{tab:milestone_accuracy}','\small','\begin{tabular}{clrrrr}','\toprule',
        'Stage & Case & Rel. $L_2$ & Max error (K) & Flux rel. $L_2$ & Full residual \\','\midrule')
    foreach ($row in $master) {
        $accuracyTex += ('{0} & {1} & ${2}$ & ${3}$ & ${4}$ & ${5}$ \\' -f
            $row.milestone_id,(Tex-Escape $row.case_name),(Tex-Number $row.relative_l2),
            (Tex-Number $row.max_node_error_k),(Tex-Number $row.flux_relative_l2),
            (Tex-Number $row.sampled_full_residual))
    }
    $accuracyTex += @('\bottomrule','\end{tabular}',
        '\begin{minipage}{0.96\linewidth}\footnotesize Values in the main rows are pure-ROM errors. Corrected-ROM residuals are reported separately in the accompanying CSV; N/A means not recorded in the archived milestone outputs.\end{minipage}',
        '\end{table}')
    Write-Utf8 (Join-Path $outputs 'milestone1_7_accuracy_table.tex') $accuracyTex

    $efficiencyTex = @(
        '\begin{table}[t]','\centering','\caption{Computational performance and break-even analysis}',
        '\label{tab:milestone_efficiency}','\small','\begin{tabular}{clrrrr}','\toprule',
        'Stage & Case & Offline (s) & Online (s/case) & Speedup & Break-even \\','\midrule')
    foreach ($row in $master) {
        $efficiencyTex += ('{0} & {1} & ${2}$ & ${3}$ & ${4}$ & {5} \\' -f
            $row.milestone_id,(Tex-Escape $row.case_name),(Tex-Number $row.offline_total_s),
            (Tex-Number $row.average_online_s_per_case),(Tex-Number $row.speedup_vs_baseline),
            (Tex-Escape $row.offline_break_even_waveforms))
    }
    $efficiencyTex += @('\bottomrule','\end{tabular}',
        '\begin{minipage}{0.96\linewidth}\footnotesize M2, M3, M6, and M7 deployment times are setup-once 100-RHS averages with full temperature recovery. Offline setup is not repeated during deployment. A speedup is shown only when workload definitions match; the M7 baseline is M6 full-interface one-step deployment.\end{minipage}',
        '\end{table}')
    Write-Utf8 (Join-Path $outputs 'milestone1_7_efficiency_table.tex') $efficiencyTex

    # Use a line array here instead of a here-string. Windows PowerShell 5.1
    # can mis-tokenize a UTF-8/no-BOM here-string that contains CJK text.
    $groupMeeting = @(
        '# Milestone 1–7 组会精简结果',
        '',
        '## 第1页：方法演进',
        '',
        '| 阶段 | 案例 | 子区域内部 | 接口 | 结论 |',
        '|---|---|---|---|---|',
        '| M1–M2 | two-/ten-cube | 独立 Local POD/RB | full-interface Schur | 稳态 Local-ROM 主线跑通 |',
        '| M3 | RRAM26 / Chiplet | 真实子区域 Local POD | full interface | 扩展到真实大模型，接口成本显现 |',
        '| M4–M5 | two-/ten-cube | Local Block Arnoldi | full Dynamic Schur | 建立瞬态、固定 dt 因子复用路径 |',
        '| M6 | RRAM26 / Chiplet | Local Block Arnoldi | 331331 / 80853 全阶接口 | 温度准确，但接口成为在线瓶颈 |',
        '| M7 | RRAM26 / Chiplet | Local Block Arnoldi + enrichment | 独立 physical-interface Port POD | 最终主线：Reduced Dynamic Schur |',
        '',
        '最终架构：DDM + 每子区域独立 Local POD/Block Arnoldi + 每物理接口独立 Port POD + Reduced Dynamic Schur + 全局温度恢复。Global ROM 仅作 benchmark。',
        '',
        '## 第2页：精度',
        '',
        '| Case | M6 relative L2 | M7 relative L2 | M7 max error (K) | M7 flux L2 | M7 full residual |',
        '|---|---:|---:|---:|---:|---:|',
        '| RRAM26 | 1.94e-6 | 4.19e-9 | 2.53e-5 | 3.48e-8 | 1.61e-7 |',
        '| Chiplet | 5.50e-6 | 1.53e-8 | 1.77e-4 | 4.47e-3 | 2.16e-3 |',
        '',
        'M7 corrected full residual：RRAM26 2.61e-12；Chiplet 8.28e-13。Pure ROM 指标用于验收，corrected 结果单列。',
        '',
        '## 第3页：效率与最终结论',
        '',
        '| Case | Interface dimension | 100-RHS average | vs M6 | Break-even |',
        '|---|---:|---:|---:|---:|',
        '| RRAM26 | 331331 → 3200 | 0.554 s/case | 58.41× | 74 waveforms |',
        '| Chiplet | 80853 → 154 | 0.236 s/case | 18.91× | 51 waveforms |',
        '',
        '- RRAM26：331331 → 3200 interface dimensions；0.554 s/case；58.41× faster than Milestone 6；relative L2 4.19e-9。',
        '- Chiplet：80853 → 154 interface dimensions；0.236 s/case；18.91× faster than Milestone 6；relative L2 1.53e-8。',
        '- M7 offline pilot 成本较高，只有多波形部署才能摊销；Global Block Arnoldi 在线仍可能更快，但不具备相同的局部替换、分区扩展和接口物理解释能力。'
    )
    Write-Utf8 (Join-Path $outputs 'milestone1_7_group_meeting_summary.md') $groupMeeting

    $finalReport = @"
# Milestone 1–7 final algorithm and results

## 1. Project objective

The project develops a partition-preserving thermal reduced-order model for repeated steady and transient power queries on nonmatching SIPG meshes. The final objective is not merely a small global state, but independently replaceable subdomain models, physically identifiable interface models, factor reuse, full-temperature recovery, and auditable temperature/flux/residual errors.

## 2. Final mathematical architecture

The final mainline is DDM + one independent Local POD or Local Block Arnoldi basis per physical subdomain + one independent Port POD basis per physical interface + Reduced Dynamic Schur + global temperature recovery. Global POD, Global Reduced Schur, and Global Block Arnoldi are retained only as benchmarks.

For fixed time step, each local reduced interior block is factored once. The full Local-ROM Dynamic Schur operator is projected with the block physical-port map, the port matrix is factored once, and subsequent time steps perform condensed-RHS construction, a small port solve, local recovery, and full-field reconstruction.

## 3. Milestone 1–7 evolution

- M1 (852c8f3): two-cube steady Local POD/RB with a full 766-DOF interface.
- M2 (da5f874): ten-cube steady extension with ten local rank-10 bases and 10,593 full interface DOFs.
- M3 (dc2d267): steady RRAM26/Chiplet validation with full physical interfaces.
- M4 (1c50542): two-cube Local Block Arnoldi and fixed-dt Dynamic Schur.
- M5 (c60d152): ten-cube transient extension and strict template audit.
- M6 (6f2d35f): transient RRAM26/Chiplet with full-interface matrix-free Dynamic Schur.
- M7 (d1755bc): independent physical-interface Port POD, port-reduced Dynamic Schur, and residual/flux-aware enrichment.

## 4. Why full-interface Local ROM was insufficient

M6 retained 331,331 RRAM26 and 80,853 Chiplet interface DOFs. Its temperature fields were accurate, but the measured setup-once 100-RHS one-step averages were 32.33 s/case and 4.47 s/case. The interface Krylov solve, not the reduced interiors, dominated. Chiplet additionally exposed an unfavorable flux relative L2 of 0.285 and sampled full residual of 0.212.

## 5. Why Global ROM was not the final mainline

Global Block Arnoldi can remain faster in serial: the archived full-horizon benchmark is 1.621 s for RRAM26 and 0.054 s for Chiplet, versus 43.10 s and 2.29 s M7 formal local-port online core. It is not selected as the mainline because it couples the entire device into one global basis and does not provide the same independent subdomain replacement, partition scaling, strict template identity, or physical-interface coordinate interpretation. These are architectural trade-offs, not a claim that the local method wins every serial timing comparison.

## 6. Junction-DOF ownership bug and its correction

An intermediate M7 implementation allowed junction interface DOFs to appear in more than one physical-port support. Port lift then summed duplicated contributions while the diagnostic projected each port independently, producing a misleading near-machine projection residual but a large thermal error. The corrected construction assigns every junction trace DOF to one deterministic physical-interface owner. Port supports are mutually exclusive, their row counts exactly cover the full SIPG interface, and the regression suite enforces this invariant.

This correction does not assume matching nodes. The existing triangle-overlap/BVH geometry, master/slave mapping, interface ordering, and original SIPG integration remain unchanged.

## 7. Local Block Arnoldi construction

Each physical subdomain builds its own moment basis from local power inputs and compressed interface excitations. PARDISO symbolic/numerical factors are reused, two-pass orthogonalization controls loss of rank, and projected capacity, conductivity, input, boundary, and coupling blocks retain DDM ownership.

## 8. Independent physical-interface Port POD

Each physical interface owns a separate shared port coordinate over its two nonmatching sides. Temperature traces, increments, weak SIPG flux information, and discrete residual snapshots are scale-balanced. Geometry modes and mandatory FOM/M6 Local-ROM pilot traces are retained before snapshot-Gram POD truncation. No global interface POD or sliced global basis is used.

## 9. Reduced Dynamic Schur

The interface state is represented by the assembled physical-port map. The projected fixed-dt Schur matrix is formed from the unchanged Local-ROM Schur apply, retains measured relative asymmetry near 1e-15, and is LU-factorized once. RRAM26 reduces 331,331 interface coordinates to 3,200; Chiplet reduces 80,853 to 154.

## 10. Residual/flux enrichment

Local full discrete residuals are evaluated at difficult waveform/time samples. Reused local interior factors solve residual corrections, two-pass MGS appends independent directions, and all local reduced blocks are reprojected. Original physical-normal and SIPG numerical flux integration remains the validation reference.

## 11. Accuracy comparison

RRAM26 M7 pure relative L2 is 4.190e-9, maximum nodal error 2.526e-5 K, flux relative L2 3.477e-8, and sampled full residual 1.610e-7. Chiplet reports 1.533e-8, 1.769e-4 K, 4.472e-3, and 2.156e-3. Corrected full residuals are 2.614e-12 and 8.275e-13. Pure and corrected results are never conflated in the generated tables.

## 12. Timing and memory comparison

The matching 100-RHS one-step full-recovery averages improve from 32.333 to 0.554 s/case on RRAM26 and from 4.467 to 0.236 s/case on Chiplet. M7 offline totals are 2,531.19 s and 271.42 s. Peak working sets increase to 20,627,955,712 and 10,138,042,368 bytes, versus M6 values of 13,847,257,088 and 6,702,108,672 bytes. Port reduction improves online deployment but is not an offline-memory reduction.

## 13. Break-even analysis

Using M6 full-interface one-step deployment as the matched baseline, the recomputed break-even counts are 74 RRAM26 waveforms and 51 Chiplet waveforms. The calculation includes the additional M7 snapshot, enrichment, full-interface pilot, and projected-Schur setup. The offline pilot alone costs 1,774.57 s and 162.66 s. M7 is therefore justified for repeated deployment, not isolated single queries.

## 14. Scalability and template reuse

Local and port bases are independently owned. Template payload reuse is allowed only under exact fingerprints; several archived configurations correctly report zero reuse because their coupling/input fingerprints differ. M2 is the positive reuse example, with 31.12% serialized-model reduction. Dense port assembly/factorization at dimension 3,200 remains a scalability limit for still larger interface spaces.

## 15. Limitations

- M7 offline pilot and projected-Schur assembly are expensive and raise peak memory.
- Global Block Arnoldi remains faster on the present serial benchmarks.
- Some early milestones did not archive FOM/ROM flux relative-L2, online peak memory, or waveform-count fields; these entries are N/A rather than reconstructed from memory.
- M4 contains a large sampled residual for one pulse diagnostic despite small temperature error; it is retained in the unified accuracy CSV.
- Multi-dt or strongly matrix-varying deployments require new local/port operators or further parametric treatment.

## 16. Recommended final configuration

Use the M7 Local Block Arnoldi + independent physical-interface Port POD + residual/flux enrichment + fixed-dt Reduced Dynamic Schur configuration for repeated RRAM26/Chiplet thermal deployment. Keep M6 full-interface Local-ROM, Global Block Arnoldi, and monolithic PARDISO as explicit baselines. Use corrected mode only when FOM-level residual recovery is required; do not use corrected accuracy to qualify the pure ROM.

## 17. Reproduction and version information

The stable solver tag is baseline-local-port-reduced-dynamic-schur-milestone7 and resolves to d1755bc on codex/feature/local-port-reduced-dynamic-schur. Its archived Release suite passed 51/51 tests. Reproduce only when desired with scripts/run_local_port_reduced_dynamic_schur.ps1; the present summary generator never launches large solvers or reads full-field files. Unified values and field-level provenance are in outputs/milestone1_7_master_summary.csv and outputs/milestone1_7_data_provenance.csv. The independent documentation commit is the commit containing this file; the stable M7 tag intentionally remains on d1755bc.
"@
    Write-Utf8 (Join-Path $docs 'milestone1_7_final_algorithm_results.md') ($finalReport -split "`r?`n")

    # Quality checks: schema width, numeric parseability, units, calculated
    # ratios, source existence, LaTeX structure, and deterministic paths.
    foreach ($entry in $csvOutputs.GetEnumerator()) {
        $path = Join-Path $outputs $entry.Key
        $rows = @(Import-Csv -LiteralPath $path)
        if ($rows.Count -eq 0) { throw "Generated CSV is empty: $path" }
        $width = @($rows[0].PSObject.Properties).Count
        foreach ($row in $rows) {
            if (@($row.PSObject.Properties).Count -ne $width) {
                throw "CSV column count mismatch: $path"
            }
        }
    }
    if ($master.Count -ne 10 -or @($methodEvolution).Count -ne 7 -or @($version).Count -ne 7) {
        throw 'Unexpected M1-M7 row count in generated summaries.'
    }
    if (@($provenance).Count -ne ($master.Count * @($master[0].PSObject.Properties).Count)) {
        throw 'Field-level provenance does not cover every master-summary cell.'
    }
    foreach ($row in $master) {
        foreach ($column in $unitMap.Keys) {
            if ($row.PSObject.Properties.Name -contains $column -and $row.$column -ne $na) {
                [void](D $row.$column)
            }
        }
        if ($row.offline_peak_memory_bytes -ne $na) {
            $expected = (D $row.offline_peak_memory_bytes) / 1MB
            if ([Math]::Abs((D $row.offline_peak_memory_mib)-$expected) -gt 1e-8*$expected) {
                throw "MiB conversion mismatch for $($row.milestone_id)/$($row.case_name)"
            }
        }
    }
    foreach ($row in @($master | Where-Object milestone_id -eq 'M7')) {
        $expectedSpeedup = (D $row.baseline_online_s_per_case)/(D $row.average_online_s_per_case)
        if ([Math]::Abs((D $row.speedup_vs_milestone6)-$expectedSpeedup) -gt 1e-10*$expectedSpeedup) {
            throw "M7 speedup mismatch for $($row.case_name)"
        }
        $pure = ($m7Comparison | Where-Object {
            $_.case_name -eq $row.case_name -and $_.method -eq 'Milestone 7 port-reduced' })[0]
        $baseline = ($m7Comparison | Where-Object {
            $_.case_name -eq $row.case_name -and $_.method -eq 'Milestone 6 full-interface' })[0]
        $expectedBreakEven = Ceil-BreakEven $pure.offline_seconds $baseline.offline_seconds `
            $baseline.online_seconds_per_one_step_rhs $pure.online_seconds_per_one_step_rhs
        if ($row.offline_break_even_waveforms -ne $expectedBreakEven) {
            throw "M7 break-even mismatch for $($row.case_name)"
        }
    }
    $archivedSourceChecks = [System.Collections.Generic.HashSet[string]]::new()
    foreach ($item in $provenance) {
        if ($item.source_file -eq $na) { continue }
        if ($item.source_file.StartsWith('git:')) {
            $parts = $item.source_file.Split(':',3)
            [void]$archivedSourceChecks.Add("$($parts[1])|$($parts[2])")
        } elseif ($item.source_file -ne 'git') {
            $path = Join-Path $project $item.source_file
            if (-not (Test-Path -LiteralPath $path)) {
                throw "Missing provenance source: $($item.source_file)"
            }
            [void]$archivedSourceChecks.Add("$($item.source_commit)|$($item.source_file)")
        }
        if ($item.source_file -match 'results/|build/|temp|debug') {
            throw "Temporary/non-archived source entered provenance: $($item.source_file)"
        }
    }
    foreach ($sourceCheck in $archivedSourceChecks) {
        $parts = $sourceCheck.Split('|',2)
        Assert-GitObject $parts[0] $parts[1]
    }
    foreach ($texName in @('milestone1_7_method_evolution_table.tex',
            'milestone1_7_accuracy_table.tex','milestone1_7_efficiency_table.tex')) {
        $text = Read-SmallText ("outputs/"+$texName)
        foreach ($token in @('\toprule','\midrule','\bottomrule','\begin{tabular}')) {
            if (-not $text.Contains($token)) { throw "Missing $token in $texName" }
        }
        foreach ($line in ($text -split "`r?`n" | Where-Object { $_ -match '\\begin\{tabular\}' })) {
            if ($line.Contains('|')) { throw "Vertical rule found in $texName" }
        }
    }

    Write-Host "Generated $($csvOutputs.Count) CSV files, 3 LaTeX tables, 1 group-meeting summary, and 1 final report."
    Write-Host 'No solver was run and no full-field file was read.'
} finally {
    Pop-Location
}
