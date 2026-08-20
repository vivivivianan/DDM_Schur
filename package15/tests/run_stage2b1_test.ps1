param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('Affine', 'Generate', 'Reload', 'OutOfRange', 'Fingerprint', 'Corrected', 'Stage2A1')]
    [string]$Mode,
    [Parameter(Mandatory = $true)]
    [string]$Exe,
    [Parameter(Mandatory = $true)]
    [string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [System.IO.Path]::GetFullPath($Root)
if (-not $rootFull.StartsWith([System.IO.Path]::GetFullPath((Join-Path $project 'build')),
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under the project build directory: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null
$model = Join-Path $rootFull 'model'
$generateOutput = Join-Path $rootFull 'generate'

function Invoke-Checked([string[]]$Arguments) {
    & $Exe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "SIPGHeatDDM3D failed with exit code $LASTEXITCODE"
    }
}

function Common-Parametric-Arguments([string]$Config, [string]$Output) {
    return @(
        '--steady', '--config', $Config, '--solvers', 'reduced-schur',
        '--mor-matrix-parameter', 'convection-h',
        '--mor-parameter-subdomain', '1', '--mor-parameter-region-id', '3',
        '--mor-parameter-min', '50', '--mor-parameter-max', '200',
        '--mor-parameter-reference', '100', '--mor-parameter-training-count', '3',
        '--mor-parameter-validation-count', '1', '--mor-parameter-test-count', '1',
        '--mor-interface-rank', '8', '--mor-local-rank', '8',
        '--output-dir', $Output, '--fast-run'
    )
}

Push-Location $project
try {
    switch ($Mode) {
        'Affine' {
            $output = Join-Path $rootFull 'affine'
            $arguments = Common-Parametric-Arguments 'configs\two_cube_parametric_h.txt' $output
            Invoke-Checked ($arguments + @('--mor-parametric-generate', '--mor-parametric-affine-validation-only'))
            $rows = Import-Csv (Join-Path $output 'parametric_affine_validation.csv')
            if ($rows.Count -ne 3) { throw 'Expected three affine validation points.' }
            foreach ($row in $rows) {
                if ([double]$row.matrix_relative_difference -ge 1.0e-12 -or
                    [double]$row.rhs_relative_difference -ge 1.0e-12) {
                    throw "Affine decomposition is not machine-precision at mu=$($row.parameter_value)."
                }
            }
        }
        'Generate' {
            $arguments = Common-Parametric-Arguments 'configs\two_cube_parametric_h.txt' $generateOutput
            Invoke-Checked ($arguments + @('--mor-parametric-generate', '--mor-parametric-save', $model,
                    '--mor-parametric-mode', 'pure'))
            $summary = Import-Csv (Join-Path $generateOutput 'parametric_stage2b1_summary.csv')
            if ($summary.status -ne 'success') { throw 'Parametric model generation did not succeed.' }
            if ([double]$summary.worst_relative_l2 -ge 1.0e-10) {
                throw "Full-rank two-cube training/validation exactness failed: $($summary.worst_relative_l2)"
            }
            if (-not (Test-Path (Join-Path $model 'model.bin'))) { throw 'Serialized model is missing.' }
        }
        'Reload' {
            $output = Join-Path $rootFull 'reload'
            Invoke-Checked (@('--mor-parametric-deployment-only', '--mor-parametric-load', $model,
                '--mor-parameter-value', '100', '--output-dir', $output))
            $expected = Import-Csv (Join-Path $generateOutput 'parametric_temperature_nodes.csv')
            $actual = Import-Csv (Join-Path $output 'parametric_temperature_nodes.csv')
            if ($expected.Count -ne $actual.Count) { throw 'Reloaded temperature size mismatch.' }
            $maximum = 0.0
            for ($index = 0; $index -lt $expected.Count; ++$index) {
                $maximum = [Math]::Max($maximum,
                    [Math]::Abs([double]$expected[$index].temperature_k - [double]$actual[$index].temperature_k))
            }
            if ($maximum -ge 1.0e-12) { throw "Serialization/reload changed the solution by $maximum K." }
        }
        'OutOfRange' {
            $output = Join-Path $rootFull 'out_of_range'
            & $Exe --mor-parametric-deployment-only --mor-parametric-load $model `
                --mor-parameter-value 250 --output-dir $output
            if ($LASTEXITCODE -eq 0) { throw 'Pure ROM silently accepted an out-of-range parameter.' }
        }
        'Fingerprint' {
            $output = Join-Path $rootFull 'fingerprint'
            $arguments = Common-Parametric-Arguments `
                'configs\two_cube_parametric_h_fingerprint_mismatch.txt' $output
            Invoke-Checked ($arguments + @('--mor-parametric-load', $model, '--mor-parametric-mode', 'pure'))
            $row = Import-Csv (Join-Path $output 'solver_comparison.csv') | Where-Object solver -Like 'Reduced-Schur*'
            if ($row.status -ne 'failed' -or $row.failure_reason -notmatch 'fingerprint rejection') {
                throw 'Changed boundary/operator fingerprint was not rejected.'
            }
        }
        'Corrected' {
            $output = Join-Path $rootFull 'corrected'
            # Keep the shared model, but suppress validation/test FOM work in this residual-gate test.
            Invoke-Checked (@('--steady', '--config', 'configs\two_cube_parametric_h.txt',
                '--solvers', 'reduced-schur', '--mor-parametric-load', $model,
                '--mor-matrix-parameter', 'convection-h', '--mor-parameter-subdomain', '1',
                '--mor-parameter-region-id', '3', '--mor-parameter-min', '50',
                '--mor-parameter-max', '200', '--mor-parameter-reference', '100',
                '--mor-parameter-training-count', '3', '--mor-parameter-validation-count', '0',
                '--mor-parameter-test-count', '0', '--mor-parametric-mode', 'corrected',
                '--pcg-tolerance', '1e-10', '--max-pcg-iterations', '200',
                '--output-dir', $output, '--fast-run'))
            $comparison = Import-Csv (Join-Path $output 'parametric_corrected_comparison.csv')
            if ([double]$comparison.true_interface_residual -ge 1.0e-8) {
                throw "Corrected residual gate failed: $($comparison.true_interface_residual)"
            }
        }
        'Stage2A1' {
            $output = Join-Path $rootFull 'stage2a1_regression'
            Invoke-Checked (@('--steady', '--config', 'configs\two_cube_schur.txt',
                '--solvers', 'reduced-schur', '--mor-generate-model',
                '--mor-save-model', (Join-Path $output 'model'), '--mor-rank', '2',
                '--mor-rank-sweep', '2', '--mor-training-count', '4',
                '--mor-validation-count', '2', '--mor-test-count', '2',
                '--mor-snapshot-solver', 'direct', '--mor-compare-fom',
                '--reduced-schur-mode', 'pure', '--mor-interior-mode', 'exact-response',
                '--mor-precompute-power-response', '--output-dir', $output, '--fast-run'))
            $summary = Get-Content (Join-Path $output 'mor_summary.json') -Raw | ConvertFrom-Json
            if ([double]$summary.nominal_global_residual -ge 1.0e-10) {
                throw "Stage 2A.1 fixed-matrix residual regression: $($summary.nominal_global_residual)"
            }
        }
    }
} finally {
    Pop-Location
}
