param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('Matrix', 'Dc', 'Step', 'Serialization', 'InitialCondition',
        'Fingerprint', 'ReducedSpd', 'TenCube')]
    [string]$Mode,
    [Parameter(Mandatory = $true)]
    [string]$Exe,
    [Parameter(Mandatory = $true)]
    [string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [System.IO.Path]::GetFullPath($Root)
$buildRoot = [System.IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under the project build directory: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null
$twoOutput = Join-Path $rootFull 'two_generate'
$twoModel = Join-Path $rootFull 'two_model'

function Invoke-Checked([string[]]$Arguments) {
    & $Exe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "SIPGHeatDDM3D failed with exit code $LASTEXITCODE"
    }
}

function Assert-Less([double]$Value, [double]$Limit, [string]$Message) {
    if (-not ($Value -lt $Limit)) {
        throw "$Message (value=$Value, limit=$Limit)"
    }
}

Push-Location $project
try {
    switch ($Mode) {
        'Matrix' {
            Invoke-Checked @(
                '--transient', '--config', 'configs\two_cube_parametric_h.txt',
                '--mor-transient-generate', '--mor-transient-save', $twoModel,
                '--mor-arnoldi-moments', '2', '--mor-transient-dt', '0.1',
                '--mor-transient-t-end', '2', '--mor-transient-waveform', 'single_step',
                '--mor-transient-compare-fom', '--mor-transient-output', 'max-temperature',
                '--output-dir', $twoOutput, '--fast-run')
            $matrix = Import-Csv (Join-Path $twoOutput 'transient_matrix_diagnostics.csv')
            if ($matrix.Count -ne 2) { throw 'Expected separate capacity and conductivity diagnostics.' }
            foreach ($row in $matrix) {
                Assert-Less ([double]$row.symmetry_error) 1.0e-12 "$($row.matrix) lost symmetry"
                if ($row.positive_on_free_dofs -ne '1') {
                    throw "$($row.matrix) has a non-positive free diagonal."
                }
            }
            if (($matrix | Where-Object matrix -eq 'capacity').units -ne 'J/K') {
                throw 'Capacity units were not reported as J/K.'
            }
            if (($matrix | Where-Object matrix -eq 'conductivity').units -ne 'W/K') {
                throw 'Conductivity units were not reported as W/K.'
            }
        }
        'Dc' {
            $rows = Import-Csv (Join-Path $twoOutput 'transient_dc_consistency.csv')
            if ($rows.Count -lt 2) { throw 'DC consistency did not cover nominal and channel cases.' }
            foreach ($row in $rows) {
                Assert-Less ([double]$row.relative_l2) 1.0e-11 "DC response was not moment-matched"
                Assert-Less ([double]$row.maximum_absolute_k) 1.0e-8 "DC maximum error was not exact"
            }
        }
        'Step' {
            $row = Import-Csv (Join-Path $twoOutput 'transient_accuracy_by_waveform.csv')
            Assert-Less ([double]$row.space_time_relative_l2) 1.0e-4 'Two-cube step relative L2 failed'
            Assert-Less ([double]$row.maximum_absolute_k) 0.1 'Two-cube step maximum error failed'
            if ([double]$row.speedup -le 10.0) { throw "Two-cube ROM speedup is only $($row.speedup)." }
        }
        'Serialization' {
            $output = Join-Path $rootFull 'pure_reload'
            Invoke-Checked @(
                '--mor-transient-deployment-only', '--mor-transient-load', $twoModel,
                '--mor-transient-dt', '0.1', '--mor-transient-t-end', '2',
                '--mor-transient-waveform', 'single_step',
                '--mor-transient-output', 'max-temperature', '--output-dir', $output)
            $expected = Import-Csv (Join-Path $twoOutput 'transient_max_temperature_curves.csv')
            $actual = Import-Csv (Join-Path $output 'transient_max_temperature_curves.csv')
            if ($expected.Count -ne $actual.Count) { throw 'Reload trajectory length changed.' }
            for ($index = 0; $index -lt $expected.Count; ++$index) {
                Assert-Less ([Math]::Abs([double]$expected[$index].rom_maximum_k -
                    [double]$actual[$index].rom_maximum_k)) 1.0e-12 'Serialization changed the trajectory'
            }
        }
        'InitialCondition' {
            $output = Join-Path $rootFull 'nonzero_initial'
            Invoke-Checked @(
                '--transient', '--config', 'configs\two_cube_parametric_h.txt',
                '--mor-transient-load', $twoModel, '--mor-transient-dt', '0.1',
                '--mor-transient-t-end', '0.5', '--mor-transient-waveform', 'single_step',
                '--mor-transient-initial-temperature', '300', '--mor-transient-compare-fom',
                '--mor-transient-output', 'max-temperature', '--output-dir', $output, '--fast-run')
            $projection = Import-Csv (Join-Path $output 'transient_initial_projection.csv')
            if ($projection.projection -ne 'C_weighted') { throw 'Initial state did not use C-weighting.' }
            Assert-Less ([double]$projection.c_weighted_orthogonality_error) 1.0e-10 `
                'Nonzero initial projection violates C-orthogonality'
            $accuracy = Import-Csv (Join-Path $output 'transient_accuracy_by_waveform.csv')
            # The input-Krylov space is not enriched with this deliberately
            # unseen initial field; this regression checks the required
            # C-weighted projection and a bounded, stable transient response.
            Assert-Less ([double]$accuracy.space_time_relative_l2) 5.0e-3 `
                'Nonzero initial transient regression failed'
        }
        'Fingerprint' {
            $output = Join-Path $rootFull 'fingerprint_reject'
            & $Exe --transient --config configs\two_cube_parametric_h_fingerprint_mismatch.txt `
                --mor-transient-load $twoModel --mor-transient-dt 0.1 `
                --mor-transient-t-end 0.2 --output-dir $output --fast-run
            if ($LASTEXITCODE -eq 0) { throw 'Transient load silently accepted changed fingerprints.' }
        }
        'ReducedSpd' {
            $summary = Import-Csv (Join-Path $twoOutput 'transient_block_arnoldi_summary.csv')
            Assert-Less ([double]$summary.basis_orthogonality_error) 1.0e-10 `
                'Block Arnoldi basis lost orthogonality'
            if ([double]$summary.capacity_min_eigenvalue -le 0.0 -or
                [double]$summary.conductivity_min_eigenvalue -le 0.0) {
                throw 'Reduced C or K is not positive definite.'
            }
            $metadata = Get-Content (Join-Path $twoModel 'metadata.json') -Raw | ConvertFrom-Json
            if ($metadata.capacity_cholesky -ne 1 -or $metadata.conductivity_cholesky -ne 1) {
                throw 'Reduced C/K Cholesky diagnostic failed.'
            }
        }
        'TenCube' {
            $output = Join-Path $rootFull 'ten_cube'
            Invoke-Checked @(
                '--transient', '--config', 'configs\ten_cube_parametric_h.txt',
                '--mor-transient-generate', '--mor-transient-save', (Join-Path $rootFull 'ten_model'),
                '--mor-arnoldi-moments', '2', '--mor-transient-dt', '0.1',
                '--mor-transient-t-end', '1', '--mor-transient-waveform', 'rectangular_pulse',
                '--mor-transient-compare-fom', '--mor-transient-output', 'max-temperature',
                '--output-dir', $output, '--fast-run')
            $summary = Import-Csv (Join-Path $output 'transient_block_arnoldi_summary.csv')
            if ($summary.status -ne 'success' -or [int]$summary.source_channels -ne 10) {
                throw 'Ten-cube transient Block Arnoldi regression failed.'
            }
            Assert-Less ([double]$summary.space_time_relative_l2) 1.0e-4 `
                'Ten-cube transient relative L2 failed'
            Assert-Less ([double]$summary.maximum_absolute_k) 0.1 `
                'Ten-cube transient maximum error failed'
        }
    }
} finally {
    Pop-Location
}
