param(
    [Parameter(Mandatory=$true)][string]$Mode,
    [Parameter(Mandatory=$true)][string]$Exe,
    [Parameter(Mandatory=$true)][string]$Root
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$config = Join-Path $repo 'configs\two_cube_schur.txt'
$fingerprintMismatchConfig = Join-Path $repo `
    'configs\two_cube_parametric_h_fingerprint_mismatch.txt'
$generate = Join-Path $Root 'generate'
$reload = Join-Path $Root 'reload'
$corrected = Join-Path $Root 'corrected'
$model = Join-Path $Root 'model'

function Invoke-Checked([string[]]$Arguments) {
    & $Exe @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Local-ROM command failed with exit code $LASTEXITCODE" }
}

function Require([bool]$Condition, [string]$Message) {
    if (-not $Condition) { throw $Message }
}

function Summary([string]$Directory) {
    return Import-Csv (Join-Path $Directory 'local_rom_schur_summary.csv') | Select-Object -First 1
}

New-Item -ItemType Directory -Force -Path $Root | Out-Null

switch ($Mode.ToLowerInvariant()) {
    'generate' {
        Invoke-Checked @(
            '--steady', '--config', $config, '--output-dir', $generate,
            '--solvers', 'local-rom', '--local-mor-generate',
            '--local-mor-save', $model, '--local-mor-method', 'pod',
            '--local-mor-mode', 'pure', '--local-interface-mode', 'full',
            '--local-mor-rank-per-subdomain', '2,2',
            '--local-mor-training-cases', '12',
            '--local-mor-validation-cases', '2', '--local-mor-test-cases', '2',
            '--local-mor-compare-fom', '--pcg-tolerance', '1e-10',
            '--max-pcg-iterations', '500', '--gmres-restart', '100', '--fast-run')
        $summary = Summary $generate
        Require ($summary.status -eq 'success') 'Local ROM accuracy gate failed.'
        Require ([int]$summary.subdomains -eq 2) 'Expected two independently reduced subdomains.'
        Require ([int]$summary.interface_dofs -gt 0) 'Full physical interface was not retained.'
        Require ([int]$summary.total_local_rank -eq 4) 'Expected two independent rank-2 bases.'
        Require ([double]$summary.relative_l2 -lt 1e-10) 'Two-cube local ROM is not exact enough.'
        Require ([double]$summary.max_node_error_k -lt 1e-6) 'Two-cube max-node error gate failed.'
        Require (Test-Path (Join-Path $model 'local_rom_model.bin')) 'Serialized model is missing.'
    }
    'assembly' {
        $summary = Summary $generate
        $ranks = @(Import-Csv (Join-Path $generate 'local_rom_rank_by_subdomain.csv'))
        Require ($ranks.Count -eq 2) 'Reduced block report must contain two subdomains.'
        Require ((($ranks | Measure-Object -Property selected_rank -Sum).Sum) -eq 4) `
            'Local reduced ranks were not assembled independently.'
        Require ([int]$summary.total_local_rank -eq 4) 'Reduced Schur total rank is inconsistent.'
    }
    'mapping' {
        $summary = Summary $generate
        $ranks = @(Import-Csv (Join-Path $generate 'local_rom_rank_by_subdomain.csv'))
        $mapped = ($ranks | Measure-Object -Property interface_dofs -Sum).Sum
        Require ([int]$mapped -eq [int]$summary.interface_dofs) `
            'Local/full interface mapping does not cover every interface DOF exactly once.'
    }
    'flux' {
        $rows = @(Import-Csv (Join-Path $generate 'local_rom_interface_flux.csv'))
        Require ($rows.Count -gt 0) 'Nonmatching interface flux report is empty.'
        $maxJump = ($rows | Measure-Object -Property temperature_jump_rms_k -Maximum).Maximum
        $maxImbalance = ($rows | Measure-Object -Property maximum_flux_imbalance_w_m2 -Maximum).Maximum
        Require ([double]$maxJump -lt 1e-8) 'Interface temperature-jump regression failed.'
        Require ([double]$maxImbalance -lt 1e-2) 'Physical flux-imbalance regression failed.'
    }
    'serialization' {
        Invoke-Checked @(
            '--steady', '--config', $config, '--output-dir', $reload,
            '--solvers', 'local-rom', '--local-mor-load', $model,
            '--local-mor-mode', 'pure', '--local-interface-mode', 'full',
            '--local-mor-validation-cases', '0', '--local-mor-test-cases', '0',
            '--no-local-mor-compare-fom', '--fast-run')
        $left = @(Import-Csv (Join-Path $generate 'temperature_local_pod_schur_rom_pure_nodes.csv'))
        $right = @(Import-Csv (Join-Path $reload 'temperature_local_pod_schur_rom_pure_nodes.csv'))
        Require ($left.Count -eq $right.Count) 'Reloaded full field has the wrong size.'
        $maximum = 0.0
        for ($i = 0; $i -lt $left.Count; ++$i) {
            $difference = [Math]::Abs([double]$left[$i].temperature_k - [double]$right[$i].temperature_k)
            if ($difference -gt $maximum) { $maximum = $difference }
        }
        Require ($maximum -le 1e-12) "Serialization changed the local ROM field: $maximum"
    }
    'fingerprint' {
        $output = Join-Path $Root 'fingerprint'
        Invoke-Checked @(
            '--steady', '--config', $fingerprintMismatchConfig,
            '--output-dir', $output, '--solvers', 'local-rom',
            '--local-mor-load', $model, '--local-mor-mode', 'pure',
            '--local-interface-mode', 'full', '--no-local-mor-compare-fom',
            '--local-mor-validation-cases', '0', '--local-mor-test-cases', '0',
            '--fast-run')
        $row = Import-Csv (Join-Path $output 'solver_comparison.csv') |
            Where-Object { $_.solver -like 'Local-POD-Schur-ROM*' } |
            Select-Object -First 1
        Require ($row.status -eq 'failed') 'Changed Local-ROM fingerprint was not rejected.'
        Require ($row.failure_reason -match 'fingerprint mismatch') `
            'Fingerprint rejection did not report the mismatch.'
    }
    'corrected' {
        Invoke-Checked @(
            '--steady', '--config', $config, '--output-dir', $corrected,
            '--solvers', 'local-rom', '--local-mor-load', $model,
            '--local-mor-mode', 'corrected', '--local-interface-mode', 'full',
            '--local-mor-validation-cases', '0', '--local-mor-test-cases', '0',
            '--local-mor-compare-fom', '--pcg-tolerance', '1e-10',
            '--max-pcg-iterations', '500', '--gmres-restart', '100', '--fast-run')
        $summary = Summary $corrected
        Require ([int]$summary.zero_guess_iterations -gt 0) 'Zero-guess Schur baseline was not measured.'
        Require ([int]$summary.local_rom_guess_iterations -eq 0) 'Exact local ROM initial guess should need zero correction iterations.'
        Require ([double]$summary.final_true_residual -lt 1e-10) 'Corrected true-residual gate failed.'
    }
    'matrixfree' {
        $output = Join-Path $Root 'matrix_free'
        Invoke-Checked @(
            '--steady', '--config', $config, '--output-dir', $output,
            '--solvers', 'local-rom', '--local-mor-load', $model,
            '--local-mor-mode', 'pure', '--local-interface-mode', 'full',
            '--local-mor-matrix-free-threshold', '0',
            '--schur-proxy-block-size', '8',
            '--local-mor-validation-cases', '0', '--local-mor-test-cases', '0',
            '--no-local-mor-compare-fom', '--pcg-tolerance', '1e-10',
            '--max-pcg-iterations', '500', '--gmres-restart', '100', '--fast-run')
        $summary = Summary $output
        Require ($summary.status -eq 'success') 'Matrix-free Local-ROM solve failed.'
        Require ($summary.interface_solver -eq 'matrix-free-fgmres') `
            'Large-interface matrix-free path was not selected.'
        Require ([int]$summary.interface_fgmres_iterations -gt 0) `
            'Matrix-free Local-ROM did not execute FGMRES.'
        Require ([double]$summary.interface_fgmres_true_residual -lt 1e-10) `
            'Matrix-free Local-ROM true residual gate failed.'
        Require ([int]$summary.proxy_colors -gt 0) `
            'Local-ROM proxy was not rebuilt from the matrix-free operator.'
        $left = @(Import-Csv (Join-Path $generate 'temperature_local_pod_schur_rom_pure_nodes.csv'))
        $right = @(Import-Csv (Join-Path $output 'temperature_local_pod_schur_rom_pure_nodes.csv'))
        Require ($left.Count -eq $right.Count) 'Matrix-free full field has the wrong size.'
        $maximum = 0.0
        for ($i = 0; $i -lt $left.Count; ++$i) {
            $difference = [Math]::Abs([double]$left[$i].temperature_k - [double]$right[$i].temperature_k)
            if ($difference -gt $maximum) { $maximum = $difference }
        }
        Require ($maximum -lt 1e-6) `
            "Matrix-free and explicit Local-ROM fields differ: $maximum K"
    }
    default { throw "Unknown Local-ROM test mode: $Mode" }
}
