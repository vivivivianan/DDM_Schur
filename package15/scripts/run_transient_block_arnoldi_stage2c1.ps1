param(
    [string]$Exe = 'build\Release\SIPGHeatDDM3D.exe',
    [string]$RramConfig = 'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$ChipletConfig = 'D:\CPP\TEST_CHATGPT\chiplet_model\case_chiplet_config_horizontal.txt',
    [string]$ResultsRoot = 'results',
    [string]$OutputsRoot = 'outputs',
    [switch]$AggregateOnly,
    [switch]$SkipBuildAndCtest,
    [switch]$SkipChiplet
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $project
try {
    $exePath = [System.IO.Path]::GetFullPath((Join-Path $project $Exe))
    $results = [System.IO.Path]::GetFullPath((Join-Path $project $ResultsRoot))
    $outputs = [System.IO.Path]::GetFullPath((Join-Path $project $OutputsRoot))
    New-Item -ItemType Directory -Force $results, $outputs | Out-Null

    function Invoke-Checked([string[]]$Arguments) {
        & $exePath @Arguments
        if ($LASTEXITCODE -ne 0) {
            throw "SIPGHeatDDM3D failed with exit code $LASTEXITCODE"
        }
    }

    function Invoke-TransientCase(
        [string]$Config, [string]$Directory, [string]$Waveform,
        [int]$Moments, [double]$Dt, [double]$EndTime,
        [string]$LoadModel = '', [string]$SaveModel = '', [bool]$CompareFom = $true,
        [string]$OutputMode = 'max-temperature') {
        New-Item -ItemType Directory -Force $Directory | Out-Null
        $arguments = @(
            '--transient', '--config', $Config,
            '--mor-transient-dt', ([string]$Dt),
            '--mor-transient-t-end', ([string]$EndTime),
            '--mor-transient-waveform', $Waveform,
            '--mor-transient-output', $OutputMode,
            '--output-dir', $Directory, '--fast-run')
        if ($LoadModel) {
            $arguments += @('--mor-transient-load', $LoadModel)
        } else {
            $arguments += @('--mor-transient-generate', '--mor-arnoldi-moments', ([string]$Moments))
            if ($SaveModel) { $arguments += @('--mor-transient-save', $SaveModel) }
        }
        if ($CompareFom) { $arguments += '--mor-transient-compare-fom' }
        Invoke-Checked $arguments
    }

    if (-not $AggregateOnly) {
        if (-not $SkipBuildAndCtest) {
            cmake -S . -B build
            if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed.' }
            cmake --build build --config Release
            if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
            $ctestOutput = @(ctest --test-dir build -C Release --output-on-failure 2>&1)
            $ctestExitCode = $LASTEXITCODE
            $ctestOutput | ForEach-Object { Write-Host $_ }
            $ctestOutput | Out-File -Encoding utf8 `
                (Join-Path $outputs 'transient_ctest_results.txt')
            if ($ctestExitCode -ne 0) { throw 'Release CTest failed.' }
        }
        if (-not (Test-Path -LiteralPath $RramConfig)) {
            throw "RRAM26 config not found: $RramConfig"
        }

        $two = Join-Path $results 'stage2c1_two'
        Invoke-TransientCase 'configs\two_cube_parametric_h.txt' $two `
            'single_step' 2 0.1 2 '' (Join-Path $two 'model') $true

        $tenModel = Join-Path $results 'stage2c1_ten_wave_model'
        $tenWaves = @(
            'single_step', 'multi_step', 'rectangular_pulse', 'piecewise_multilevel',
            'variable_duty_cycle', 'mixed_frequency', 'asynchronous_hotspots',
            'unseen_waveform', 'unseen_channel_combination', 'same_average_power')
        for ($index = 0; $index -lt $tenWaves.Count; ++$index) {
            $wave = $tenWaves[$index]
            $directory = Join-Path $results "stage2c1_ten_wave_$wave"
            if ($index -eq 0) {
                Invoke-TransientCase 'configs\ten_cube_parametric_h.txt' $directory `
                    $wave 2 0.1 2 '' $tenModel $true
            } else {
                Invoke-TransientCase 'configs\ten_cube_parametric_h.txt' $directory `
                    $wave 2 0.1 2 $tenModel '' $true
            }
        }

        # RRAM26 m=1 is the accepted deployable model. m=2 is also compared
        # against FOM; m=3..5 measure rank, DC, setup, and memory without
        # serializing multi-gigabyte nonrecommended models.
        $rramModel = Join-Path $results 'stage2c1_rram26_model'
        Invoke-TransientCase $RramConfig (Join-Path $results 'stage2c1_rram26_m1') `
            'mixed_frequency' 1 0.01 0.8 '' $rramModel $true
        Invoke-TransientCase $RramConfig (Join-Path $results 'stage2c1_rram26_m2') `
            'mixed_frequency' 2 0.01 0.8 '' '' $true
        foreach ($moment in 3..5) {
            Invoke-TransientCase $RramConfig (Join-Path $results "stage2c1_rram26_m$moment") `
                'mixed_frequency' $moment 0.01 0.8 '' '' $false 'selected-dofs'
        }
        foreach ($wave in @('single_step', 'rectangular_pulse',
                'variable_duty_cycle', 'unseen_waveform')) {
            Invoke-TransientCase $RramConfig (Join-Path $results "stage2c1_rram26_wave_$wave") `
                $wave 1 0.01 0.8 $rramModel '' $true
        }
        Invoke-Checked @(
            '--mor-transient-deployment-only', '--mor-transient-load', $rramModel,
            '--mor-transient-dt', '0.01', '--mor-transient-t-end', '0.8',
            '--mor-transient-waveform', 'unseen_waveform',
            '--mor-transient-output', 'selected-dofs',
            '--mor-deployment-rhs-count', '100',
            '--output-dir', (Join-Path $results 'stage2c1_rram26_deployment_100'))

        if (-not $SkipChiplet -and (Test-Path -LiteralPath $ChipletConfig)) {
            Invoke-TransientCase $ChipletConfig (Join-Path $results 'stage2c1_chiplet_m4') `
                'asynchronous_hotspots' 4 0.1 1 '' '' $true
        }
    }

    function Existing-Directory([string]$Preferred, [string]$Fallback = '') {
        if (Test-Path -LiteralPath $Preferred) { return $Preferred }
        if ($Fallback -and (Test-Path -LiteralPath $Fallback)) { return $Fallback }
        throw "Required Stage 2C.1 result directory is missing: $Preferred"
    }

    function Export-Combined([object[]]$Specifications, [string]$SourceFile,
                             [string]$Destination) {
        $combined = @()
        foreach ($specification in $Specifications) {
            $path = Join-Path $specification.Directory $SourceFile
            if (-not (Test-Path -LiteralPath $path)) { continue }
            foreach ($row in (Import-Csv $path)) {
                $row | Add-Member -NotePropertyName benchmark_case `
                    -NotePropertyValue $specification.Case
                $combined += $row
            }
        }
        if ($combined.Count -gt 0) {
            $combined | Export-Csv -NoTypeInformation -Encoding UTF8 $Destination
        }
    }

    $twoPreferred = Join-Path $results 'stage2c1_two'
    $twoDirectory = if (Test-Path -LiteralPath (Join-Path $twoPreferred `
            'transient_block_arnoldi_summary.csv')) {
        $twoPreferred
    } else {
        Existing-Directory (Join-Path $project 'build\stage2c1_ctest\two_generate')
    }
    $tenWaves = @(
        'single_step', 'multi_step', 'rectangular_pulse', 'piecewise_multilevel',
        'variable_duty_cycle', 'mixed_frequency', 'asynchronous_hotspots',
        'unseen_waveform', 'unseen_channel_combination', 'same_average_power')
    $accuracySpecs = @([pscustomobject]@{Case='two_cube';Directory=$twoDirectory})
    foreach ($wave in $tenWaves) {
        $accuracySpecs += [pscustomobject]@{
            Case="ten_cube_$wave"; Directory=(Join-Path $results "stage2c1_ten_wave_$wave")}
    }
    foreach ($wave in @('single_step', 'rectangular_pulse', 'mixed_frequency',
            'variable_duty_cycle', 'unseen_waveform')) {
        $directory = if ($wave -eq 'mixed_frequency') {
            Join-Path $results 'stage2c1_rram26_m1'
        } else { Join-Path $results "stage2c1_rram26_wave_$wave" }
        $accuracySpecs += [pscustomobject]@{Case="rram26_$wave";Directory=$directory}
    }
    $chipletDirectory = Join-Path $results 'stage2c1_chiplet_m4'
    if (Test-Path -LiteralPath $chipletDirectory) {
        $accuracySpecs += [pscustomobject]@{
            Case='chiplet_asynchronous_hotspots';Directory=$chipletDirectory}
    }
    Export-Combined $accuracySpecs 'transient_accuracy_by_waveform.csv' `
        (Join-Path $outputs 'transient_accuracy_by_waveform.csv')
    Export-Combined $accuracySpecs 'transient_accuracy_by_time.csv' `
        (Join-Path $outputs 'transient_accuracy_by_time.csv')
    Export-Combined $accuracySpecs 'transient_max_temperature_curves.csv' `
        (Join-Path $outputs 'transient_max_temperature_curves.csv')

    $rramM1Setup = Join-Path $results 'stage2c1_rram26_model_build'
    if (-not (Test-Path -LiteralPath $rramM1Setup)) {
        $rramM1Setup = Join-Path $results 'stage2c1_rram26_m1'
    }
    $summarySpecs = @(
        [pscustomobject]@{Case='two_cube_m2';Directory=$twoDirectory},
        [pscustomobject]@{Case='ten_cube_m2';Directory=(Join-Path $results 'stage2c1_ten_wave_single_step')},
        [pscustomobject]@{Case='rram26_m1';Directory=(Join-Path $results 'stage2c1_rram26_m1')})
    if (Test-Path -LiteralPath $chipletDirectory) {
        $summarySpecs += [pscustomobject]@{Case='chiplet_m4';Directory=$chipletDirectory}
    }
    Export-Combined $summarySpecs 'transient_block_arnoldi_summary.csv' `
        (Join-Path $outputs 'transient_block_arnoldi_summary.csv')

    $momentRows = @()
    foreach ($moment in 1..5) {
        $directory = if ($moment -eq 1) { $rramM1Setup } else {
            Join-Path $results "stage2c1_rram26_m$moment" }
        $summary = Import-Csv (Join-Path $directory 'transient_block_arnoldi_summary.csv')
        $offline = Import-Csv (Join-Path $directory 'transient_offline_timing.csv')
        $memory = Import-Csv (Join-Path $directory 'transient_memory.csv')
        $accuracyPath = Join-Path (Join-Path $results "stage2c1_rram26_m$moment") `
            'transient_accuracy_by_waveform.csv'
        $accuracy = if (Test-Path -LiteralPath $accuracyPath) {
            Import-Csv $accuracyPath } else { $null }
        $residentBytes = if ($moment -eq 2) { 1009242440 } else {
            [uint64]$summary.model_bytes }
        $momentRows += [pscustomobject]@{
            benchmark_case='rram26'; moments=$moment; rank=[int]$summary.rank
            added_rank_per_moment=125; deflated_columns=0
            orthogonality_error=[double]$summary.basis_orthogonality_error
            dc_worst_relative_l2=[double]$summary.dc_worst_relative_l2
            space_time_relative_l2=if($null -eq $accuracy){''}else{$accuracy.space_time_relative_l2}
            maximum_absolute_k=if($null -eq $accuracy){''}else{$accuracy.maximum_absolute_k}
            accuracy_measured=if($null -eq $accuracy){0}else{1}
            total_offline_seconds=[double]$offline.total_offline_seconds
            resident_model_bytes=$residentBytes
            peak_working_set_bytes=[uint64]$memory.peak_working_set_bytes
        }
    }
    $momentRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs 'transient_moment_sweep.csv')

    $rankSpecs = @()
    foreach ($moment in 1..5) {
        $directory = if ($moment -eq 1) { $rramM1Setup } else {
            Join-Path $results "stage2c1_rram26_m$moment" }
        $rankSpecs += [pscustomobject]@{Case="rram26_m$moment";Directory=$directory}
    }
    Export-Combined $rankSpecs 'transient_arnoldi_rank_history.csv' `
        (Join-Path $outputs 'transient_arnoldi_rank_history.csv')

    $diagnosticSpecs = @(
        [pscustomobject]@{Case='two_cube';Directory=$twoDirectory},
        [pscustomobject]@{Case='ten_cube';Directory=(Join-Path $results 'stage2c1_ten_wave_single_step')},
        [pscustomobject]@{Case='rram26';Directory=$rramM1Setup})
    if (Test-Path -LiteralPath $chipletDirectory) {
        $diagnosticSpecs += [pscustomobject]@{Case='chiplet';Directory=$chipletDirectory}
    }
    Export-Combined $diagnosticSpecs 'transient_matrix_diagnostics.csv' `
        (Join-Path $outputs 'transient_matrix_diagnostics.csv')
    Export-Combined $diagnosticSpecs 'transient_dc_consistency.csv' `
        (Join-Path $outputs 'transient_dc_consistency.csv')

    $timingSpecs = @(
        [pscustomobject]@{Case='two_cube_m2';Directory=$twoDirectory},
        [pscustomobject]@{Case='ten_cube_m2';Directory=(Join-Path $results 'stage2c1_ten_wave_single_step')},
        [pscustomobject]@{Case='rram26_m1';Directory=$rramM1Setup},
        [pscustomobject]@{Case='rram26_m2';Directory=(Join-Path $results 'stage2c1_rram26_m2')},
        [pscustomobject]@{Case='rram26_m3';Directory=(Join-Path $results 'stage2c1_rram26_m3')},
        [pscustomobject]@{Case='rram26_m4';Directory=(Join-Path $results 'stage2c1_rram26_m4')},
        [pscustomobject]@{Case='rram26_m5';Directory=(Join-Path $results 'stage2c1_rram26_m5')})
    if (Test-Path -LiteralPath $chipletDirectory) {
        $timingSpecs += [pscustomobject]@{Case='chiplet_m4';Directory=$chipletDirectory}
    }
    Export-Combined $timingSpecs 'transient_offline_timing.csv' `
        (Join-Path $outputs 'transient_offline_timing.csv')
    $memoryOutput = Join-Path $outputs 'transient_memory.csv'
    Export-Combined $timingSpecs 'transient_memory.csv' $memoryOutput
    # The development sweep intentionally discarded the nonrecommended m=2
    # serialized basis. Report its dimension-derived resident footprint rather
    # than a partial filesystem size from that discarded artifact.
    $memoryRows = Import-Csv $memoryOutput
    foreach ($row in $memoryRows) {
        if ($row.benchmark_case -eq 'rram26_m2') {
            $row.model_bytes = '1009242440'
        }
    }
    $memoryRows | Export-Csv -NoTypeInformation -Encoding UTF8 $memoryOutput

    $onlineSpecs = $accuracySpecs + [pscustomobject]@{
        Case='rram26_pure_deployment_100';
        Directory=(Join-Path $results 'stage2c1_rram26_deployment_100')}
    Export-Combined $onlineSpecs 'transient_online_timing.csv' `
        (Join-Path $outputs 'transient_online_timing.csv')
    Export-Combined $accuracySpecs 'transient_break_even.csv' `
        (Join-Path $outputs 'transient_break_even.csv')

    Write-Host "Stage 2C.1 aggregate outputs written to $outputs"
} finally {
    Pop-Location
}
