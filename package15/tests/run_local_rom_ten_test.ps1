param(
    [Parameter(Mandatory=$true)][string]$Mode,
    [Parameter(Mandatory=$true)][string]$Exe,
    [Parameter(Mandatory=$true)][string]$Root
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$config = Join-Path $repo 'configs\ten_cube_schur.txt'
$mismatchConfig = Join-Path $repo 'configs\ten_cube_schur_fingerprint_mismatch.txt'
$generate = Join-Path $Root 'generate'
$reload = Join-Path $Root 'reload'
$corrected = Join-Path $Root 'corrected'
$reuse = Join-Path $Root 'reuse'
$model = Join-Path $Root 'model'
$reuseModel = Join-Path $Root 'reuse_model'

function Invoke-Checked([string[]]$Arguments) {
    & $Exe @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Ten-cube Local-ROM command failed: $LASTEXITCODE" }
}

function Require([bool]$Condition, [string]$Message) {
    if (-not $Condition) { throw $Message }
}

function Summary([string]$Directory) {
    return Import-Csv (Join-Path $Directory 'local_rom_schur_summary.csv') | Select-Object -First 1
}

function Maximum-Temperature-Difference([string]$LeftDirectory, [string]$RightDirectory) {
    $left = @(Import-Csv (Join-Path $LeftDirectory 'temperature_local_pod_schur_rom_pure_nodes.csv'))
    $right = @(Import-Csv (Join-Path $RightDirectory 'temperature_local_pod_schur_rom_pure_nodes.csv'))
    Require ($left.Count -eq $right.Count) 'Compared ten-cube fields have different sizes.'
    $maximum = 0.0
    for ($i = 0; $i -lt $left.Count; ++$i) {
        $difference = [Math]::Abs([double]$left[$i].temperature_k - [double]$right[$i].temperature_k)
        if ($difference -gt $maximum) { $maximum = $difference }
    }
    return $maximum
}

New-Item -ItemType Directory -Force -Path $Root | Out-Null

switch ($Mode.ToLowerInvariant()) {
    'generate' {
        Invoke-Checked @(
            '--steady', '--config', $config, '--output-dir', $generate,
            '--solvers', 'local-rom', '--local-mor-generate', '--local-mor-save', $model,
            '--local-mor-method', 'pod', '--local-mor-mode', 'pure',
            '--local-interface-mode', 'full', '--local-mor-rank', '10',
            '--local-mor-training-cases', '20', '--local-mor-validation-cases', '2',
            '--local-mor-test-cases', '2', '--local-mor-compare-fom',
            '--pcg-tolerance', '1e-10', '--fast-run')
        $summary = Summary $generate
        Require ($summary.status -eq 'success') 'Ten-cube pure Local-ROM accuracy gate failed.'
        Require ([int]$summary.subdomains -eq 10) 'Expected ten physical subdomains.'
        Require ([int]$summary.interface_dofs -eq 10593) 'Unexpected full interface dimension.'
        Require ([int]$summary.total_local_rank -eq 100) 'Expected ten independent rank-10 bases.'
        Require (Test-Path (Join-Path $model 'local_rom_model.bin')) 'Ten-cube model is missing.'
    }
    'structure' {
        $domains = @(Import-Csv (Join-Path $generate 'local_rom_subdomain_structure.csv'))
        $interfaces = @(Import-Csv (Join-Path $generate 'local_rom_interface_structure.csv'))
        Require ($domains.Count -eq 10) 'Subdomain structure table must contain ten rows.'
        Require ($interfaces.Count -eq 9) 'Interface structure table must contain nine physical pairs.'
        Require ((($interfaces | Measure-Object face_pair_count -Sum).Sum) -eq 4554) 'Face-pair count changed.'
        Require ((($interfaces | Measure-Object triangle_overlap_count -Sum).Sum) -eq 6840) 'Triangle-overlap count changed.'
        foreach ($row in $domains) {
            Require ([int]$row.global_unique_interface_dofs -eq 10593) 'Global unique interface count changed.'
            Require ([int]$row.unique_trace_union -eq 10593) 'Explicit interface union changed.'
            Require ([int]$row.all_local_trace_references -eq 10593) 'Trace-reference diagnostic changed.'
        }
    }
    'localbases' {
        $ranks = @(Import-Csv (Join-Path $generate 'local_rom_rank_by_subdomain.csv'))
        Require ($ranks.Count -eq 10) 'Rank table must contain ten independent local bases.'
        Require ((($ranks | Measure-Object selected_rank -Sum).Sum) -eq 100) 'Local rank sum changed.'
        Require (@($ranks | Where-Object { [int]$_.template_reused -ne 0 }).Count -eq 0) `
            'Independent run silently reused a local basis.'
        $families = @((Import-Csv (Join-Path $generate 'local_rom_training_cases.csv')).family | Select-Object -Unique)
        foreach ($required in @('unit_channel','adjacent_pair','first_last_pair','middle_source',
                'multiple_hotspots','all_uniform','same_total_endpoint_skew','remote_interface_driven')) {
            Require ($families -contains $required) "Missing training family: $required"
        }
        Invoke-Checked @(
            '--steady', '--config', $config, '--output-dir', $reuse,
            '--solvers', 'local-rom', '--local-mor-generate', '--local-mor-save', $reuseModel,
            '--local-mor-method', 'pod', '--local-mor-mode', 'pure',
            '--local-interface-mode', 'full', '--local-mor-rank', '10',
            '--local-mor-training-cases', '20', '--local-mor-validation-cases', '0',
            '--local-mor-test-cases', '0', '--no-local-mor-compare-fom',
            '--local-rom-reuse-identical-subdomains', '--fast-run')
        $template = Import-Csv (Join-Path $reuse 'local_rom_template_reuse.csv') | Select-Object -First 1
        Require ([int]$template.unique_template_count -eq 3) 'Expected exactly three fingerprint templates.'
        Require ([int]$template.reused_instance_count -eq 7) 'Expected seven reused middle instances.'
        $reuseRanks = @(Import-Csv (Join-Path $reuse 'local_rom_rank_by_subdomain.csv'))
        $worstConsistency = ($reuseRanks | Measure-Object template_consistency_difference -Maximum).Maximum
        Require ([double]$worstConsistency -lt 1e-10) 'Template reduced blocks are inconsistent.'
        $fieldDifference = Maximum-Temperature-Difference $generate $reuse
        Require ($fieldDifference -lt 1e-7) "Template reuse changed the ten-cube field: $fieldDifference"
    }
    'reducedblocks' {
        $ranks = @(Import-Csv (Join-Path $generate 'local_rom_rank_by_subdomain.csv'))
        foreach ($row in $ranks) {
            Require ([double]$row.reduced_interior_symmetry_error -lt 1e-10) 'Reduced A_II lost symmetry.'
            Require ([double]$row.coupling_symmetry_error -lt 1e-10) 'Reduced coupling blocks lost transpose symmetry.'
            Require ([double]$row.local_schur_symmetry_error -lt 1e-10) 'Local reduced Schur correction lost symmetry.'
            Require ([double]$row.reduced_interior_min_eigenvalue -gt 0) 'Reduced A_II is not SPD.'
            Require ([double]$row.reduced_interior_condition_estimate -lt 1e8) 'Reduced A_II conditioning gate failed.'
        }
    }
    'interfaceassembly' {
        $summary = Summary $generate
        Require ($summary.interface_solver -eq 'sparse-pardiso') 'Large full interface did not use exact sparse PARDISO.'
        Require ([int64]$summary.interface_nnz -gt 1000000) 'Reduced Schur interface assembly is unexpectedly empty.'
        $offline = Import-Csv (Join-Path $generate 'local_rom_offline_timing.csv') | Select-Object -First 1
        Require ([double]$offline.interface_numerical_factorization_seconds -gt 0) 'Interface numerical factorization was not timed.'
    }
    'mapping' {
        $domains = @(Import-Csv (Join-Path $generate 'local_rom_subdomain_structure.csv'))
        $localTraceSum = ($domains | Measure-Object local_interface_dofs -Sum).Sum
        Require ([int]$localTraceSum -eq 10593) 'Local trace map does not cover the full DG interface.'
        Require ((@($domains.mesh_fingerprint | Select-Object -Unique)).Count -eq 1) 'Translation-invariant mesh fingerprint changed.'
        Require ((@($domains.template_id | Select-Object -Unique)).Count -eq 3) 'Boundary/interface fingerprints did not distinguish endpoints.'
    }
    'pure' {
        $rows = @(Import-Csv (Join-Path $generate 'local_rom_accuracy_by_case.csv'))
        Require ($rows.Count -eq 5) 'Expected nominal plus four held-out cases.'
        $worstL2 = ($rows | Measure-Object relative_l2 -Maximum).Maximum
        $worstMax = ($rows | Measure-Object max_node_error_k -Maximum).Maximum
        Require ([double]$worstL2 -lt 1e-5) 'Ten-cube held-out relative-L2 gate failed.'
        Require ([double]$worstMax -lt 0.1) 'Ten-cube held-out maximum-error gate failed.'
        $timing = @(Import-Csv (Join-Path $generate 'local_rom_multi_rhs_timing.csv'))
        Require (($timing.rhs_count -join ',') -eq '1,10,100') 'Missing 1/10/100 RHS timing sweep.'
        $flux = @(Import-Csv (Join-Path $generate 'local_rom_interface_flux.csv'))
        Require ($flux.Count -gt 0) 'Ten-cube flux diagnostics are empty.'
        Require ($null -ne $flux[0].worst_integration_triangle) 'Worst integration triangle was not reported.'
    }
    'serialization' {
        Invoke-Checked @(
            '--steady', '--config', $config, '--output-dir', $reload,
            '--solvers', 'local-rom', '--local-mor-load', $model,
            '--local-mor-mode', 'pure', '--local-interface-mode', 'full',
            '--local-mor-validation-cases', '0', '--local-mor-test-cases', '0',
            '--no-local-mor-compare-fom', '--fast-run')
        $difference = Maximum-Temperature-Difference $generate $reload
        Require ($difference -le 1e-8) "Ten-cube serialization changed the field: $difference"
    }
    'fingerprint' {
        $output = Join-Path $Root 'fingerprint'
        Invoke-Checked @(
            '--steady', '--config', $mismatchConfig, '--output-dir', $output,
            '--solvers', 'local-rom', '--local-mor-load', $model,
            '--local-mor-mode', 'pure', '--local-interface-mode', 'full',
            '--local-mor-validation-cases', '0', '--local-mor-test-cases', '0',
            '--no-local-mor-compare-fom', '--fast-run')
        $row = Import-Csv (Join-Path $output 'solver_comparison.csv') |
            Where-Object { $_.solver -like 'Local-POD-Schur-ROM*' } | Select-Object -First 1
        Require ($row.status -eq 'failed') 'Changed ten-cube physics fingerprint was not rejected.'
        Require ($row.failure_reason -match 'fingerprint mismatch') 'Fingerprint rejection reason is missing.'
    }
    'corrected' {
        Invoke-Checked @(
            '--steady', '--config', $config, '--output-dir', $corrected,
            '--solvers', 'local-rom', '--local-mor-load', $model,
            '--local-mor-mode', 'corrected', '--local-interface-mode', 'full',
            '--local-mor-validation-cases', '0', '--local-mor-test-cases', '0',
            '--local-mor-compare-fom', '--pcg-tolerance', '1e-10', '--fast-run')
        $summary = Summary $corrected
        Require ([int]$summary.zero_guess_iterations -gt 0) 'Ten-cube zero-guess baseline was not measured.'
        Require ([int]$summary.local_rom_guess_iterations -le 2) 'Local-ROM guess did not eliminate correction work.'
        Require ([double]$summary.final_true_residual -lt 1e-10) 'Corrected true-residual gate failed.'
    }
    default { throw "Unknown ten-cube Local-ROM test mode: $Mode" }
}
