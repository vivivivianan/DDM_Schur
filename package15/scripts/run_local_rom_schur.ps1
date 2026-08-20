param(
    [string]$BuildDirectory = '.\build',
    [string]$ResultsDirectory = '.\results',
    [switch]$SkipBuild,
    [switch]$SkipBaselines
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot

function Invoke-Checked([string]$Label, [scriptblock]$Command) {
    & $Command
    if ($LASTEXITCODE -ne 0) { throw "$Label failed with exit code $LASTEXITCODE." }
}

function Get-DirectoryBytes([string]$Path) {
    return [int64]((Get-ChildItem -LiteralPath $Path -Recurse -File |
        Measure-Object -Property Length -Sum).Sum)
}

function Compare-TemperatureFields([string]$LeftPath, [string]$RightPath) {
    $left = @(Import-Csv $LeftPath)
    $right = @(Import-Csv $RightPath)
    if ($left.Count -ne $right.Count) { throw 'Template comparison field sizes differ.' }
    $differenceSquared = 0.0
    $referenceSquared = 0.0
    $maximum = 0.0
    for ($i = 0; $i -lt $left.Count; ++$i) {
        $reference = [double]$left[$i].temperature_k
        $difference = [double]$right[$i].temperature_k - $reference
        $differenceSquared += $difference * $difference
        $referenceSquared += $reference * $reference
        $maximum = [Math]::Max($maximum, [Math]::Abs($difference))
    }
    return [pscustomobject]@{
        RelativeL2 = [Math]::Sqrt($differenceSquared) / [Math]::Max(1e-300, [Math]::Sqrt($referenceSquared))
        Maximum = $maximum
    }
}

Push-Location $repo
try {
    if (-not $SkipBuild) {
        Invoke-Checked 'Release build' {
            cmake --build $BuildDirectory --config Release -- /m:1 /nodeReuse:false /clp:ErrorsOnly
        }
    }
    $exe = Join-Path $BuildDirectory 'Release\SIPGHeatDDM3D.exe'
    $tenConfig = '.\configs\ten_cube_schur.txt'
    $twoConfig = '.\configs\two_cube_schur.txt'
    $pure = Join-Path $ResultsDirectory 'local_rom_milestone2_ten_cube'
    $reuse = Join-Path $ResultsDirectory 'local_rom_milestone2_ten_cube_reuse'
    $reload = Join-Path $ResultsDirectory 'local_rom_milestone2_ten_cube_reload'
    $corrected = Join-Path $ResultsDirectory 'local_rom_milestone2_ten_cube_corrected'
    $model = Join-Path $ResultsDirectory 'local_rom_milestone2_model'
    $reuseModel = Join-Path $ResultsDirectory 'local_rom_milestone2_reuse_model'
    $twoRegression = Join-Path $ResultsDirectory 'local_rom_milestone1_two_cube_regression'

    Invoke-Checked 'Ten-cube independent Local-ROM' {
        & $exe --steady --config $tenConfig --output-dir $pure `
            --solvers local-rom --local-mor-generate --local-mor-save $model `
            --local-mor-method pod --local-mor-mode pure --local-interface-mode full `
            --local-mor-rank 10 --local-mor-training-cases 20 `
            --local-mor-validation-cases 8 --local-mor-test-cases 8 `
            --local-mor-compare-fom --pcg-tolerance 1e-10 --fast-run
    }

    Invoke-Checked 'Ten-cube identical-subdomain template reuse' {
        & $exe --steady --config $tenConfig --output-dir $reuse `
            --solvers local-rom --local-mor-generate --local-mor-save $reuseModel `
            --local-mor-method pod --local-mor-mode pure --local-interface-mode full `
            --local-mor-rank 10 --local-mor-training-cases 20 `
            --local-mor-validation-cases 0 --local-mor-test-cases 0 `
            --no-local-mor-compare-fom --local-rom-reuse-identical-subdomains --fast-run
    }

    Invoke-Checked 'Ten-cube model reload and multi-RHS reuse' {
        & $exe --steady --config $tenConfig --output-dir $reload `
            --solvers local-rom --local-mor-load $model --local-mor-mode pure `
            --local-interface-mode full --local-mor-validation-cases 0 `
            --local-mor-test-cases 0 --no-local-mor-compare-fom --fast-run
    }

    Invoke-Checked 'Ten-cube corrected Local-ROM' {
        & $exe --steady --config $tenConfig --output-dir $corrected `
            --solvers local-rom --local-mor-load $model --local-mor-mode corrected `
            --local-interface-mode full --local-mor-validation-cases 0 `
            --local-mor-test-cases 0 --local-mor-compare-fom `
            --pcg-tolerance 1e-10 --fast-run
    }

    Invoke-Checked 'Milestone 1 two-cube regression' {
        & $exe --steady --config $twoConfig --output-dir $twoRegression `
            --solvers local-rom --local-mor-generate --local-mor-method pod `
            --local-mor-mode pure --local-interface-mode full `
            --local-mor-rank-per-subdomain '2,2' --local-mor-training-cases 12 `
            --local-mor-validation-cases 2 --local-mor-test-cases 2 `
            --local-mor-compare-fom --pcg-tolerance 1e-10 --fast-run
    }

    $baseline = Join-Path $ResultsDirectory 'local_rom_milestone2_stage1_and_direct'
    $globalRom = Join-Path $ResultsDirectory 'local_rom_milestone2_global_reduced_schur'
    if (-not $SkipBaselines) {
        Invoke-Checked 'Ten-cube Stage 1 and monolithic baselines' {
            & $exe --steady --config $tenConfig --output-dir $baseline `
                --solvers schur,direct --direct-mode spd --pcg-tolerance 1e-10 --fast-run
        }
        Invoke-Checked 'Ten-cube Global Reduced Schur benchmark' {
            & $exe --steady --config $tenConfig --output-dir $globalRom `
                --solvers reduced-schur --mor-generate-model --mor-rank 10 `
                --mor-rank-sweep 10 --mor-training-count 20 `
                --mor-validation-count 4 --mor-test-count 4 `
                --mor-snapshot-solver direct --mor-compare-fom `
                --reduced-schur-mode pure --mor-interior-mode pardiso `
                --pcg-tolerance 1e-10 --fast-run
        }
    }

    $outputs = Join-Path $repo 'outputs'
    New-Item -ItemType Directory -Force -Path $outputs | Out-Null
    foreach ($name in @(
        'local_rom_schur_report.md', 'local_rom_schur_summary.csv',
        'local_rom_rank_by_subdomain.csv', 'local_rom_accuracy_by_case.csv',
        'local_rom_interface_flux.csv', 'local_rom_offline_timing.csv',
        'local_rom_online_timing.csv',
        'local_rom_memory.csv', 'local_rom_subdomain_structure.csv',
        'local_rom_interface_structure.csv', 'local_rom_training_cases.csv',
        'local_rom_validation_cases.csv', 'local_rom_test_cases.csv',
        'local_rom_source_channels.csv')) {
        Copy-Item -LiteralPath (Join-Path $pure $name) -Destination (Join-Path $outputs $name) -Force
    }
    Copy-Item -LiteralPath (Join-Path $reload 'local_rom_multi_rhs_timing.csv') `
        -Destination (Join-Path $outputs 'local_rom_multi_rhs_timing.csv') -Force
    $offlineScenarios = @()
    foreach ($specification in @(
        @{ Scenario='generation'; Directory=$pure },
        @{ Scenario='serialized_model_reload'; Directory=$reload })) {
        $offlineRow = Import-Csv (Join-Path $specification.Directory 'local_rom_offline_timing.csv') |
            Select-Object -First 1
        $offlineScenarios += $offlineRow | Select-Object `
            @{Name='scenario'; Expression={$specification.Scenario}}, *
    }
    $offlineScenarios | Export-Csv -NoTypeInformation `
        (Join-Path $outputs 'local_rom_offline_timing.csv')
    $postprocessing = Import-Csv (Join-Path $pure 'program_timing.csv') |
        Where-Object { $_.stage -eq 'postprocessing' } | Select-Object -First 1
    [pscustomobject]@{
        scope = 'application_fast_run_temperature_and_summary_output_io'
        seconds = $postprocessing.seconds
        included_in_core_online_time = 0
    } | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_output_io_timing.csv')

    $fieldDifference = Compare-TemperatureFields `
        (Join-Path $pure 'temperature_local_pod_schur_rom_pure_nodes.csv') `
        (Join-Path $reuse 'temperature_local_pod_schur_rom_pure_nodes.csv')
    if ($fieldDifference.Maximum -ge 1e-7 -or $fieldDifference.RelativeL2 -ge 1e-10) {
        throw "Template-reuse field gate failed: relative=$($fieldDifference.RelativeL2), max=$($fieldDifference.Maximum) K"
    }
    $template = Import-Csv (Join-Path $reuse 'local_rom_template_reuse.csv') | Select-Object -First 1
    $independentBytes = Get-DirectoryBytes $model
    $reuseBytes = Get-DirectoryBytes $reuseModel
    [pscustomobject]@{
        reuse_enabled = 1
        unique_template_count = $template.unique_template_count
        subdomain_instance_count = $template.subdomain_instance_count
        reused_instance_count = $template.reused_instance_count
        basis_storage_without_reuse_bytes = $template.basis_storage_without_reuse_bytes
        basis_storage_with_reuse_bytes = $template.basis_storage_with_reuse_bytes
        model_size_reduction = 1.0 - $reuseBytes / [Math]::Max(1.0, $independentBytes)
        independent_serialized_model_bytes = $independentBytes
        reused_serialized_model_bytes = $reuseBytes
        online_relative_l2_difference = $fieldDifference.RelativeL2
        maximum_temperature_difference_k = $fieldDifference.Maximum
    } | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_template_reuse.csv')

    $comparison = @()
    $local = Import-Csv (Join-Path $pure 'local_rom_schur_summary.csv') | Select-Object -First 1
    $correctedSummary = Import-Csv (Join-Path $corrected 'local_rom_schur_summary.csv') | Select-Object -First 1
    foreach ($specification in @(
        @{ Method='Local-POD-Schur-ROM'; Row=$local; Scope='10 independent local bases'; Interface='full Schur'; Iterations=0 },
        @{ Method='Local-POD-Schur-ROM-Corrected'; Row=$correctedSummary; Scope='10 independent local bases'; Interface='full exact Schur correction'; Iterations=[int]$correctedSummary.local_rom_guess_iterations })) {
        $comparison += [pscustomobject]@{
            method=$specification.Method; basis_scope=$specification.Scope
            interface_scope=$specification.Interface; rank_or_coarse_dim=$specification.Row.total_local_rank
            setup_seconds=$specification.Row.offline_seconds
            solve_seconds=$specification.Row.nominal_online_seconds
            iterations=$specification.Iterations; relative_l2=$specification.Row.relative_l2
            max_node_error_k=$specification.Row.max_node_error_k; status=$specification.Row.status
        }
    }
    if (-not $SkipBaselines) {
        foreach ($specification in @(
            @{ Directory=$globalRom; Pattern='Reduced-Schur'; Scope='global interface basis'; Interface='reduced' },
            @{ Directory=$baseline; Pattern='DDM-Schur'; Scope='none/exact FOM'; Interface='full Schur' },
            @{ Directory=$baseline; Pattern='Global-PARDISO-SPD'; Scope='none/monolithic FOM'; Interface='monolithic' })) {
            $row = Import-Csv (Join-Path $specification.Directory 'solver_comparison.csv') |
                Where-Object { $_.solver -like ('*' + $specification.Pattern + '*') } | Select-Object -First 1
            if ($null -ne $row) {
                $relativeL2 = $row.relative_l2_diff_vs_global
                $maximumError = $row.max_abs_diff_vs_global
                if ($specification.Pattern -eq 'Reduced-Schur') {
                    $heldOut = Import-Csv (Join-Path $globalRom 'mor_case_results.csv') |
                        Where-Object { $_.split -in @('validation', 'test') }
                    $relativeL2 = ($heldOut | ForEach-Object {
                        [double]$_.full_field_relative_l2
                    } | Measure-Object -Maximum).Maximum
                    $maximumError = ($heldOut | ForEach-Object {
                        [double]$_.maximum_absolute_error
                    } | Measure-Object -Maximum).Maximum
                }
                $comparison += [pscustomobject]@{
                    method=$row.solver; basis_scope=$specification.Scope; interface_scope=$specification.Interface
                    rank_or_coarse_dim=$row.coarse_dim; setup_seconds=$row.setup_seconds
                    solve_seconds=$row.solve_seconds; iterations=$row.iterations
                    relative_l2=$relativeL2
                    max_node_error_k=$maximumError; status=$row.status
                }
            }
        }
    }
    $comparison | Export-Csv -NoTypeInformation (Join-Path $outputs 'local_rom_vs_global_rom.csv')

    $reportPath = Join-Path $outputs 'local_rom_schur_report.md'
    Add-Content -Encoding UTF8 -Path $reportPath -Value @(
        '', '## Corrected and template-reuse validation', '',
        '| zero-guess iterations | Local-ROM-guess iterations | initial true residual | final true residual | corrected solve (s) |',
        '|---:|---:|---:|---:|---:|',
        ('| {0} | {1} | {2} | {3} | {4} |' -f $correctedSummary.zero_guess_iterations,
            $correctedSummary.local_rom_guess_iterations, $correctedSummary.initial_true_residual,
            $correctedSummary.final_true_residual, $correctedSummary.nominal_online_seconds), '',
        ('Strict fingerprints produced {0} unique templates for 10 instances; {1} instances reused payloads.' -f
            $template.unique_template_count, $template.reused_instance_count),
        ('Independent/reused field difference: relative L2 {0}, maximum {1} K. Serialized model reduction: {2:P2}.' -f
            $fieldDifference.RelativeL2, $fieldDifference.Maximum,
            (1.0 - $reuseBytes / [Math]::Max(1.0, $independentBytes)))
    )
    Add-Content -Encoding UTF8 -Path $reportPath -Value @(
        '', '## Serialized-model multi-RHS deployment', '',
        '| RHS count | source projection total (s) | core online total (s) | average online (s) |',
        '|---:|---:|---:|---:|'
    )
    foreach ($timing in @(Import-Csv (Join-Path $reload 'local_rom_multi_rhs_timing.csv'))) {
        Add-Content -Encoding UTF8 -Path $reportPath -Value `
            ('| {0} | {1} | {2} | {3} |' -f $timing.rhs_count,
                $timing.source_projection_total_seconds, $timing.online_total_seconds,
                $timing.average_online_seconds)
    }
    Add-Content -Encoding UTF8 -Path $reportPath -Value @(
        '', ('Serialized model load: {0} s; one full interface setup/factorization is reused by all RHS.' -f
            ((Import-Csv (Join-Path $reload 'local_rom_offline_timing.csv') | Select-Object -First 1).model_load_seconds)), '',
        '## Method comparison', '',
        '| method | basis scope | interface scope | rank/coarse dim | setup (s) | solve (s) | iterations | relative L2 | max-node error (K) |',
        '|---|---|---|---:|---:|---:|---:|---:|---:|'
    )
    foreach ($item in $comparison) {
        Add-Content -Encoding UTF8 -Path $reportPath -Value `
            ('| {0} | {1} | {2} | {3} | {4} | {5} | {6} | {7} | {8} |' -f
                $item.method, $item.basis_scope, $item.interface_scope, $item.rank_or_coarse_dim,
                $item.setup_seconds, $item.solve_seconds, $item.iterations,
                $item.relative_l2, $item.max_node_error_k)
    }

    $ctestOutput = & ctest --test-dir $BuildDirectory -C Release --output-on-failure 2>&1
    $ctestOutput | Set-Content -Encoding UTF8 (Join-Path $outputs 'local_rom_ctest_results.txt')
    if ($LASTEXITCODE -ne 0) { throw 'Full CTest regression failed.' }
}
finally {
    Pop-Location
}
