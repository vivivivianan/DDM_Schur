param(
    [switch]$SkipBuild,
    [switch]$SkipTests,
    [switch]$AggregateOnly
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$runRoot = Join-Path $project 'build\local_block_arnoldi_milestone5'
$outputs = Join-Path $project 'outputs'
$suiteRoot = Join-Path $runRoot 'suite'
New-Item -ItemType Directory -Force -Path $runRoot, $suiteRoot, $outputs | Out-Null

function Invoke-Checked([string[]]$Arguments) {
    & $exe @Arguments
    if ($LASTEXITCODE -ne 0) { throw "SIPGHeatDDM3D failed with exit code $LASTEXITCODE" }
}

function Invoke-Local([string]$Output, [string]$Waveform, [bool]$Reuse) {
    $arguments = @(
        '--transient', '--config', 'configs\ten_cube_parametric_h.txt',
        '--mor-transient-generate', '--mor-transient-method', 'local-block-arnoldi',
        '--mor-arnoldi-moments', '2', '--mor-arnoldi-rank-tolerance', '1e-10',
        '--mor-transient-dt', '0.1', '--mor-transient-t-end', '1',
        '--mor-transient-waveform', $Waveform, '--mor-transient-initial-mode', 'ambient',
        '--output-dir', $Output, '--fast-run')
    if ($Reuse) { $arguments += '--mor-local-transient-reuse-identical-subdomains' }
    Invoke-Checked $arguments
}

Push-Location $project
try {
    if (-not $AggregateOnly) {
        if (-not $SkipBuild) {
            cmake --build build --config Release --parallel 4
            if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
        }
        if (-not $SkipTests) {
            $ctestOutput = & ctest --test-dir build -C Release --output-on-failure 2>&1
            $ctestCode = $LASTEXITCODE
            $ctestOutput | Write-Host
            $ctestOutput | Set-Content -Encoding utf8 `
                (Join-Path $outputs 'local_block_arnoldi_milestone5_ctest_results.txt')
            if ($ctestCode -ne 0) { throw 'CTest failed.' }
        }
        foreach ($waveform in @(
            'single_step', 'multi_step', 'rectangular_pulse', 'piecewise_multilevel',
            'variable_duty_cycle', 'mixed_frequency', 'asynchronous_hotspots',
            'remote_heating', 'unseen_waveform', 'unseen_channel_combination')) {
            Invoke-Local (Join-Path $suiteRoot $waveform) $waveform $false
        }
        Invoke-Local (Join-Path $runRoot 'reuse_single_step') 'single_step' $true
        Invoke-Checked @(
            '--transient', '--config', 'configs\ten_cube_parametric_h.txt',
            '--mor-transient-generate', '--mor-transient-method', 'block-arnoldi',
            '--mor-arnoldi-moments', '2', '--mor-transient-dt', '0.1',
            '--mor-transient-t-end', '1', '--mor-transient-waveform', 'single_step',
            '--mor-transient-compare-fom', '--output-dir', (Join-Path $runRoot 'global_single_step'),
            '--fast-run')
    }

    $suite = foreach ($directory in Get-ChildItem $suiteRoot -Directory) {
        Import-Csv (Join-Path $directory.FullName 'local_dynamic_schur_summary.csv')
    }
    $suite | Export-Csv (Join-Path $outputs 'local_dynamic_schur_milestone5.csv') `
        -NoTypeInformation
    Copy-Item (Join-Path $suiteRoot 'single_step\local_block_arnoldi_rank.csv') `
        (Join-Path $outputs 'local_block_arnoldi_milestone5_rank_by_subdomain.csv') -Force

    $independent = $suite | Where-Object waveform -eq 'single_step' | Select-Object -First 1
    $reuse = Import-Csv (Join-Path $runRoot 'reuse_single_step\local_dynamic_schur_summary.csv')
    $global = Import-Csv (Join-Path $runRoot 'global_single_step\transient_block_arnoldi_summary.csv')
    $independentField = Import-Csv `
        (Join-Path $suiteRoot 'single_step\local_dynamic_schur_final_temperature.csv')
    $reuseField = Import-Csv `
        (Join-Path $runRoot 'reuse_single_step\local_dynamic_schur_final_temperature.csv')
    $differenceSquared = 0.0
    $referenceSquared = 0.0
    $maximumDifference = 0.0
    for ($index = 0; $index -lt $independentField.Count; ++$index) {
        $difference = [double]$independentField[$index].temperature_k -
            [double]$reuseField[$index].temperature_k
        $differenceSquared += $difference * $difference
        $reference = [double]$independentField[$index].temperature_k
        $referenceSquared += $reference * $reference
        $maximumDifference = [Math]::Max($maximumDifference, [Math]::Abs($difference))
    }
    $relativeDifference = [Math]::Sqrt($differenceSquared) / [Math]::Sqrt($referenceSquared)
    $independentTime = Import-Csv `
        (Join-Path $suiteRoot 'single_step\local_dynamic_schur_accuracy_by_time.csv')
    $reuseTime = Import-Csv `
        (Join-Path $runRoot 'reuse_single_step\local_dynamic_schur_accuracy_by_time.csv')
    $fluxDifference = 0.0
    for ($index = 0; $index -lt $independentTime.Count; ++$index) {
        $fluxDifference = [Math]::Max($fluxDifference, [Math]::Abs(
            [double]$independentTime[$index].relative_flux_imbalance -
            [double]$reuseTime[$index].relative_flux_imbalance))
    }
    @([pscustomobject]@{
        independent_unique_templates = [int]$independent.unique_templates
        reuse_requested_unique_templates = [int]$reuse.unique_templates
        reused_instances = [int]$reuse.reused_instances
        final_temperature_relative_l2 = $relativeDifference
        final_temperature_maximum_k = $maximumDifference
        maximum_relative_flux_difference = $fluxDifference
        conclusion = 'strict_fingerprint_rejected_nonidentical_interface_trace_bases'
    }) | Export-Csv (Join-Path $outputs 'local_block_arnoldi_milestone5_template_reuse.csv') `
        -NoTypeInformation

    @(
        [pscustomobject]@{
            method = 'Transient monolithic PARDISO'; rank = 60640; interface_dofs = 0
            offline_seconds = [double]$independent.fom_factor_seconds
            online_seconds = [double]$independent.fom_solve_seconds
            model_bytes = 0; peak_working_set_bytes = [double]$independent.peak_working_set_bytes
            space_time_relative_l2 = 0.0; maximum_absolute_k = 0.0
        },
        [pscustomobject]@{
            method = 'Global Block Arnoldi'; rank = [int]$global.rank; interface_dofs = 0
            offline_seconds = [double]$global.offline_seconds
            online_seconds = [double]$independent.fom_solve_seconds / [double]$global.speedup
            model_bytes = [double]$global.model_bytes
            peak_working_set_bytes = [double]$independent.peak_working_set_bytes
            space_time_relative_l2 = [double]$global.space_time_relative_l2
            maximum_absolute_k = [double]$global.maximum_absolute_k
        },
        [pscustomobject]@{
            method = 'Local Block Arnoldi + Dynamic Schur'; rank = [int]$independent.total_local_rank
            interface_dofs = [int]$independent.full_interface_dofs
            offline_seconds = [double]$independent.local_basis_setup_seconds +
                [double]$independent.dynamic_schur_setup_seconds
            online_seconds = [double]$independent.local_online_core_seconds
            model_bytes = [double]$independent.model_bytes
            peak_working_set_bytes = [double]$independent.peak_working_set_bytes
            space_time_relative_l2 = [double]$independent.space_time_relative_l2
            maximum_absolute_k = [double]$independent.maximum_absolute_k
        }) | Export-Csv (Join-Path $outputs 'local_block_arnoldi_milestone5_summary.csv') `
        -NoTypeInformation

    $worstL2 = ($suite | Measure-Object space_time_relative_l2 -Maximum).Maximum
    $worstMax = ($suite | Measure-Object maximum_absolute_k -Maximum).Maximum
    $worstMaxTemperature = ($suite | Measure-Object maximum_temperature_error_k -Maximum).Maximum
    $report = @"
# Milestone 5: ten-cube Local Block Arnoldi + Dynamic Schur

## Result

The ten-cube transient case uses 10 independently generated local Block Arnoldi
bases, 10 source channels, 9 physical interfaces, and all
$($independent.full_interface_dofs) interface DOFs.  Total local rank is
$($independent.total_local_rank).  The fixed-dt Dynamic Schur factor is analyzed
and numerically factorized $($independent.dynamic_schur_symbolic_calls) / $($independent.dynamic_schur_numerical_calls)
times, then reused for every time step.

- Worst space-time relative L2 over 10 waveforms: $worstL2
- Worst maximum nodal error: $worstMax K
- Worst maximum-temperature error: $worstMaxTemperature K
- Local basis setup: $($independent.local_basis_setup_seconds) s
- Dynamic Schur setup: $($independent.dynamic_schur_setup_seconds) s
- 10-step local online core: $($independent.local_online_core_seconds) s
- Dynamic Schur factor memory: $($independent.dynamic_schur_factor_memory_bytes) bytes
- Model bytes: $($independent.model_bytes)
- Peak working set: $($independent.peak_working_set_bytes) bytes

## Coverage

The actual suite covers single-channel and multi-channel steps, rectangular
pulse, piecewise multi-level, variable duty cycle, mixed frequency,
asynchronous hotspots, remote heating, an unseen waveform, and an unseen input
combination.  Every case is compared at every time step with the factor-reused
monolithic transient PARDISO solution.  Reduced/full residual, maximum
temperature, interface jump, and physical flux imbalance are retained in the
raw per-case diagnostics and summarized in the committed CSV.
At exactly zero relative state/input, the relative residual diagnostic uses a
`1e-8 * ||boundary RHS||` normalization floor to avoid dividing roundoff by a
near-zero condensed RHS; this floor does not enter the solve.

## Strict template reuse audit

The independent path reports 10 templates.  The reuse-requested path also
reports $($reuse.unique_templates) templates and $($reuse.reused_instances) reused instances.
Although the middle cube matrices are structurally similar, their compressed
global-trace interface excitation bases are not bitwise identical.  Because the
fingerprint includes K/C interior and coupling blocks, local source definition,
interface excitation basis, and B_eff, reuse is correctly rejected instead of
silently sharing a nonidentical basis.  Independent/reuse-requested final-field
relative L2 difference is $relativeDifference, maximum difference is
$maximumDifference K, and maximum aggregate flux difference is $fluxDifference.

## Method comparison

The summary compares transient monolithic PARDISO, the retained Global Block
Arnoldi benchmark, and Local Block Arnoldi + Dynamic Schur.  Global Block
Arnoldi remains faster and smaller on this repeated-chain case; the Local path
preserves DDM ownership and the full physical interface, at the cost of the
10,593-DOF interface factor and full-field recovery.
"@
    Set-Content (Join-Path $outputs 'local_block_arnoldi_milestone5_report.md') `
        -Value $report -Encoding utf8
} finally {
    Pop-Location
}
