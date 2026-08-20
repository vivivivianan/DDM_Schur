param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('HeterogeneousSipg', 'PenaltyScaling', 'NipgReject',
        'FullResidualGate')]
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
if (-not $rootFull.StartsWith(
        $buildRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

function Assert-Less([double]$Value, [double]$Limit, [string]$Message) {
    if (-not ($Value -lt $Limit)) {
        throw "$Message (value=$Value, limit=$Limit)"
    }
}

function Invoke-Checked([string[]]$Arguments) {
    $log = & $Exe @Arguments 2>&1
    if ($LASTEXITCODE -ne 0) {
        $log | ForEach-Object { Write-Host $_ }
        throw "SIPGHeatDDM3D failed with exit code $LASTEXITCODE"
    }
}

function New-TwoCubeConfig(
    [string]$Path,
    [double]$LeftConductivity,
    [double]$RightConductivity,
    [double]$PenaltyFactor,
    [string]$PenaltyMode,
    [string]$PenaltyScaling,
    [string]$InterfaceScheme = 'sipg') {
    $source = Join-Path $project 'configs\two_cube_schur.txt'
    [string[]]$lines = Get-Content -LiteralPath $source
    $domain = 0
    $fineMesh = (Join-Path $project 'data\mesh\cube_fine.mphtxt')
    $coarseMesh = (Join-Path $project 'data\mesh\cube_coarse_top.mphtxt')
    for ($index = 0; $index -lt $lines.Count; ++$index) {
        if ($lines[$index] -match '^\s*name\s*=') {
            $lines[$index] = "name = numerical_correctness_$([IO.Path]::GetFileNameWithoutExtension($Path))"
        } elseif ($lines[$index] -match '^\s*output_dir\s*=') {
            $lines[$index] = "output_dir = $rootFull"
        } elseif ($lines[$index] -match '^\s*penalty_factor\s*=') {
            $lines[$index] = "penalty_factor = $($PenaltyFactor.ToString('R', [Globalization.CultureInfo]::InvariantCulture))"
        } elseif ($lines[$index] -match '^\s*penalty_mode\s*=') {
            $lines[$index] = "penalty_mode = $PenaltyMode"
        } elseif ($lines[$index] -match '^\s*penalty_scaling\s*=') {
            $lines[$index] = "penalty_scaling = $PenaltyScaling"
        } elseif ($lines[$index] -match '^\s*interface_scheme\s*=') {
            $lines[$index] = "interface_scheme = $InterfaceScheme"
        } elseif ($lines[$index] -match '^\s*domain\s*=') {
            $parts = @($lines[$index].Split(','))
            if ($domain -eq 0) {
                $parts[0] = "domain = $fineMesh"
                $parts[2] = " $($LeftConductivity.ToString('R', [Globalization.CultureInfo]::InvariantCulture))"
            } else {
                $parts[0] = "domain = $coarseMesh"
                $parts[2] = " $($RightConductivity.ToString('R', [Globalization.CultureInfo]::InvariantCulture))"
            }
            $lines[$index] = $parts -join ','
            ++$domain
        }
    }
    [IO.File]::WriteAllLines($Path, $lines)
}

function Read-Metric([string]$Path, [string]$Name) {
    $row = Import-Csv -LiteralPath $Path | Where-Object metric -eq $Name
    if ($null -eq $row) { throw "Missing metric '$Name' in $Path" }
    return $row.value
}

Push-Location $project
try {
    switch ($Mode) {
        'HeterogeneousSipg' {
            $results = @()
            foreach ($ratio in @(1.0, 130.0, 1000.0)) {
                $config = Join-Path $rootFull "heterogeneous_$([int]$ratio).txt"
                $output = Join-Path $rootFull "heterogeneous_$([int]$ratio)"
                New-TwoCubeConfig $config $ratio 1.0 20.0 'harmonic' 'p1_squared'
                Invoke-Checked @(
                    '--steady', '--config', $config, '--output-dir', $output,
                    '--solvers', 'schur,direct', '--direct-mode', 'spd',
                    '--pcg-tolerance', '1e-10', '--max-pcg-iterations', '500',
                    '--gmres-restart', '100')

                $symmetry = [double](Read-Metric `
                    (Join-Path $output 'final_matrix_spd_diagnostics.csv') `
                    'symmetry_error')
                Assert-Less $symmetry 1.0e-12 `
                    "SIPG matrix lost symmetry for kLeft/kRight=$ratio"
                if ((Read-Metric (Join-Path $output `
                        'final_matrix_spd_diagnostics.csv') `
                        'pardiso_spd_factorization') -ne 'success') {
                    throw "PARDISO SPD failed for kLeft/kRight=$ratio"
                }
                if ([int](Read-Metric (Join-Path $output `
                        'final_matrix_spd_diagnostics.csv') `
                        'ldlt_negative_count') -ne 0) {
                    throw "Negative LDLT curvature detected for kLeft/kRight=$ratio"
                }
                $comparison = Import-Csv (Join-Path $output 'solver_comparison.csv') |
                    Where-Object solver -eq 'DDM-Schur-FGMRES'
                if ($comparison.status -ne 'success') {
                    throw "Schur solve failed for kLeft/kRight=$ratio"
                }
                Assert-Less ([double]$comparison.relative_l2_diff_vs_global) `
                    1.0e-10 "Schur/monolithic mismatch for kLeft/kRight=$ratio"
                $schur = Import-Csv (Join-Path $output 'schur_solver_summary.csv')
                $results += [pscustomobject]@{
                    conductivity_ratio = $ratio
                    symmetry_error = $symmetry
                    pardiso_spd = 'success'
                    ldlt_negative_count = 0
                    schur_iterations = [int]$schur.iterations
                    schur_monolithic_relative_l2 =
                        [double]$comparison.relative_l2_diff_vs_global
                }
            }
            $results | Export-Csv (Join-Path $rootFull `
                'heterogeneous_sipg_results.csv') -NoTypeInformation
        }

        'PenaltyScaling' {
            $rows = @()
            foreach ($penaltyMode in @('max', 'harmonic')) {
                foreach ($scaling in @('p_p1', 'p1_squared')) {
                    $caseName = "penalty_$($penaltyMode)_$scaling"
                    $config = Join-Path $rootFull "$caseName.txt"
                    $output = Join-Path $rootFull $caseName
                    New-TwoCubeConfig $config 40.0 130.0 20.0 `
                        $penaltyMode $scaling
                    Invoke-Checked @(
                        '--steady', '--config', $config, '--output-dir', $output,
                        '--diagnostics-only', '--solvers', 'direct',
                        '--direct-mode', 'spd')
                    $metadata = Import-Csv (Join-Path $output 'run_metadata.csv')
                    if ([double]$metadata.penalty_factor -ne 20.0) {
                        throw 'Config penalty_factor=20 was overwritten.'
                    }
                    if ($metadata.penalty_mode -ne $penaltyMode) {
                        throw "Config penalty_mode=$penaltyMode was overwritten."
                    }
                    if ($metadata.penalty_scaling -ne $scaling) {
                        throw "Config penalty_scaling=$scaling was overwritten."
                    }
                    $face = Import-Csv (Join-Path $output `
                        'rram_interface_face_pair_diagnostics.csv') |
                        Select-Object -First 1
                    $pScale = if ($scaling -eq 'p1_squared') { 9.0 } else { 6.0 }
                    $alphaLeft = [double]$face.k_left / [double]$face.h_left
                    $alphaRight = [double]$face.k_right / [double]$face.h_right
                    $alphaScale = if ($penaltyMode -eq 'max') {
                        [Math]::Max($alphaLeft, $alphaRight)
                    } else {
                        2.0 * $alphaLeft * $alphaRight / ($alphaLeft + $alphaRight)
                    }
                    $expected = 20.0 * $pScale * $alphaScale
                    Assert-Less ([Math]::Abs([double]$face.sigma - $expected) /
                        [Math]::Max(1.0, [Math]::Abs($expected))) 1.0e-12 `
                        "$penaltyMode penalty did not use pScale=$pScale and factor 20"
                    $stats = Import-Csv (Join-Path $output `
                        'interface_penalty_stats.csv')
                    $rows += [pscustomobject]@{
                        penalty_mode = $penaltyMode
                        penalty_scaling = $scaling
                        p_scale = $pScale
                        penalty_factor = [double]$metadata.penalty_factor
                        average_penalty = [double]$stats.avg_eta_F
                    }
                }
            }
            foreach ($penaltyMode in @('max', 'harmonic')) {
                $modeRows = @($rows | Where-Object penalty_mode -eq $penaltyMode)
                $ratio = [double]$modeRows[1].average_penalty /
                    [double]$modeRows[0].average_penalty
                Assert-Less ([Math]::Abs($ratio - 1.5)) 1.0e-12 `
                    "$penaltyMode P2 p1_squared/p_p1 ratio is not 9/6."
            }
            $rows | Export-Csv (Join-Path $rootFull `
                'penalty_scaling_results.csv') -NoTypeInformation
        }

        'NipgReject' {
            $config = Join-Path $rootFull 'nipg_reject.txt'
            $output = Join-Path $rootFull 'nipg_output_must_not_exist'
            New-TwoCubeConfig $config 40.0 130.0 20.0 'harmonic' 'p_p1' 'nipg'
            $savedErrorAction = $ErrorActionPreference
            $ErrorActionPreference = 'Continue'
            $log = & $Exe --steady --config $config --output-dir $output `
                --solvers schur --direct-mode spd 2>&1
            $nipgExitCode = $LASTEXITCODE
            $ErrorActionPreference = $savedErrorAction
            if ($nipgExitCode -eq 0) {
                throw 'NIPG Schur path silently reported success.'
            }
            $text = $log -join "`n"
            if ($text -notmatch 'NIPG produces a nonsymmetric system' -or
                $text -notmatch 'Use interface_scheme=sipg') {
                throw "NIPG rejection message was unclear: $text"
            }
            if (Test-Path $output) {
                $temperatures = Get-ChildItem -LiteralPath $output -Filter `
                    'temperature*.csv' -ErrorAction SilentlyContinue
                if (@($temperatures).Count -ne 0) {
                    throw 'NIPG rejection wrote a successful temperature result.'
                }
            }
        }

        'FullResidualGate' {
            $low = Join-Path $rootFull 'low_rank_fallback'
            Invoke-Checked @(
                '--transient', '--config', 'configs\two_cube_parametric_h.txt',
                '--mor-transient-generate', '--mor-transient-method',
                'local-block-arnoldi', '--mor-arnoldi-moments', '1',
                '--mor-transient-dt', '0.1', '--mor-transient-t-end', '0.2',
                '--mor-transient-waveform', 'single_step',
                '--mor-full-residual-tolerance', '1e-8',
                '--mor-full-residual-fallback',
                '--output-dir', $low, '--fast-run')
            $lowSummary = Import-Csv (Join-Path $low `
                'local_dynamic_schur_summary.csv')
            if ($lowSummary.residual_gate_passed -ne '1' -or
                [int]$lowSummary.residual_fallback_steps -lt 1) {
                throw 'Low-rank ROM did not trigger and pass FOM residual fallback.'
            }
            if (-not ([double]$lowSummary.maximum_full_residual_before_gate -gt
                    [double]$lowSummary.full_residual_tolerance)) {
                throw 'Low-rank pre-gate residual did not exceed the tolerance.'
            }
            Assert-Less ([double]$lowSummary.maximum_full_residual) `
                ([double]$lowSummary.full_residual_tolerance) `
                'FOM fallback did not satisfy the residual gate'
            $lowByTime = Import-Csv (Join-Path $low `
                'local_dynamic_schur_accuracy_by_time.csv')
            foreach ($field in @('full_residual_before_gate',
                    'full_residual_tolerance', 'residual_gate_passed',
                    'residual_fallback_used', 'full_residual_after_fallback')) {
                if (-not $lowByTime[0].PSObject.Properties[$field]) {
                    throw "Per-step residual gate field '$field' is missing."
                }
            }

            $high = Join-Path $rootFull 'sufficient_rank_pass'
            Invoke-Checked @(
                '--transient', '--config', 'configs\two_cube_parametric_h.txt',
                '--mor-transient-generate', '--mor-transient-method',
                'local-block-arnoldi', '--mor-arnoldi-moments', '3',
                '--mor-transient-dt', '0.1', '--mor-transient-t-end', '0.2',
                '--mor-transient-waveform', 'single_step',
                '--mor-full-residual-tolerance', '1e-5',
                '--output-dir', $high, '--fast-run')
            $highSummary = Import-Csv (Join-Path $high `
                'local_dynamic_schur_summary.csv')
            if ($highSummary.residual_gate_passed -ne '1' -or
                [int]$highSummary.residual_fallback_steps -ne 0) {
                throw 'Sufficient-rank ROM did not pass the residual gate directly.'
            }
        }
    }
} finally {
    Pop-Location
}
