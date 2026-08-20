param(
    [ValidateSet('tolerance', 'block', 'multi-rhs', 'all')]
    [string]$Mode = 'all',
    [string]$Configuration = 'Release',
    [string]$RramConfig = 'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [switch]$ValidateBlockEquivalence,
    [switch]$SummarizeOnly
)

$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$executable = Join-Path $projectRoot "build\$Configuration\SIPGHeatDDM3D.exe"
$resultsRoot = Join-Path $projectRoot 'results'
$outputsRoot = Join-Path $projectRoot 'outputs'

if (-not (Test-Path -LiteralPath $executable)) {
    throw "Executable not found: $executable"
}
if (-not (Test-Path -LiteralPath $RramConfig)) {
    throw "RRAM26 config not found: $RramConfig"
}
New-Item -ItemType Directory -Force -Path $resultsRoot, $outputsRoot | Out-Null

function Invoke-RramSchur(
    [string]$Name,
    [double]$Tolerance,
    [int]$BlockSize,
    [bool]$WithDirect,
    [int]$MultiRhsCount = 0,
    [bool]$ValidateBlock = $false
) {
    $output = Join-Path $resultsRoot $Name
    $solvers = if ($WithDirect) { 'schur,direct' } else { 'schur' }
    $arguments = @(
        '--steady', '--config', $RramConfig,
        '--output-dir', $output,
        '--solvers', $solvers,
        '--direct-mode', 'spd',
        '--schur-proxy-ring', '1',
        '--schur-proxy-use-material-connectivity',
        '--schur-proxy-high-k-threshold', '20',
        '--schur-proxy-block-size', $BlockSize,
        '--pcg-tolerance', $Tolerance.ToString('R', [Globalization.CultureInfo]::InvariantCulture),
        '--max-pcg-iterations', '500',
        '--gmres-restart', '100',
        '--fast-run'
    )
    if ($MultiRhsCount -gt 0) {
        $arguments += @('--schur-multi-rhs-count', $MultiRhsCount)
    }
    if ($ValidateBlock) {
        $arguments += '--schur-proxy-validate-block-equivalence'
    }
    & $executable @arguments | ForEach-Object { Write-Host $_ }
    if ($LASTEXITCODE -ne 0) {
        throw "$Name failed with exit code $LASTEXITCODE"
    }
    return $output
}

if ($Mode -eq 'tolerance' -or $Mode -eq 'all') {
    $rows = @()
    foreach ($tolerance in @(1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8)) {
        $label = $tolerance.ToString('0e+0', [Globalization.CultureInfo]::InvariantCulture).Replace('+', '')
        $directory = Join-Path $resultsRoot "rram26_stage2_tol_$label"
        if (-not $SummarizeOnly) {
            $directory = Invoke-RramSchur "rram26_stage2_tol_$label" $tolerance 64 $true
        }
        $schur = Import-Csv (Join-Path $directory 'schur_solver_summary.csv')
        $comparison = Import-Csv (Join-Path $directory 'solver_comparison.csv') |
            Where-Object solver -eq 'DDM-Schur-FGMRES'
        $residual = Import-Csv (Join-Path $directory 'rram_true_residual_check.csv') |
            Where-Object solver -eq 'DDM-Schur-FGMRES'
        $rows += [pscustomobject]@{
            tolerance = $tolerance
            coarse_dimension = $schur.coarse_dimension
            iterations = $schur.iterations
            interface_residual = $schur.relative_residual
            global_residual = $residual.true_relative_residual
            maximum_temperature_difference = $comparison.max_abs_diff_vs_global
            relative_l2 = $comparison.relative_l2_diff_vs_global
            setup_seconds = $schur.setup_seconds
            solve_seconds = $schur.total_solve_seconds
            total_seconds = $schur.total_seconds
            proxy_colors = $schur.proxy_colors
            proxy_block_size = $schur.proxy_probing_block_size
            status = $schur.status
        }
    }
    $rows | Export-Csv (Join-Path $outputsRoot 'schur_stage2_tolerance_sweep.csv') -NoTypeInformation
}

if ($Mode -eq 'block' -or $Mode -eq 'all') {
    $rows = @()
    foreach ($blockSize in @(8, 16, 32, 64)) {
        $directory = Join-Path $resultsRoot "rram26_stage2_block_$blockSize"
        if (-not $SummarizeOnly) {
            $directory = Invoke-RramSchur "rram26_stage2_block_$blockSize" 1.0e-5 $blockSize $false 0 $false
        }
        $schur = Import-Csv (Join-Path $directory 'schur_solver_summary.csv')
        $rows += [pscustomobject]@{
            block_size = $blockSize
            proxy_colors = $schur.proxy_colors
            exact_probing_vectors = $schur.proxy_probing_schur_applies
            phase33_multi_rhs_calls = $schur.proxy_probing_block_calls
            proxy_nnz = $schur.proxy_nnz
            proxy_value_hash = $schur.proxy_value_hash
            iterations = $schur.iterations
            setup_seconds = $schur.setup_seconds
            solve_seconds = $schur.total_solve_seconds
            total_seconds = $schur.total_seconds
            status = $schur.status
        }
    }
    $rows | Export-Csv (Join-Path $outputsRoot 'schur_stage2_block_probing_sweep.csv') -NoTypeInformation

    if ($ValidateBlockEquivalence) {
        $validationRows = @()
        foreach ($blockSize in @(8, 16, 32, 64)) {
            $directory = Join-Path $resultsRoot "rram26_stage2_block_${blockSize}_validation"
            if (-not $SummarizeOnly) {
                $directory = Invoke-RramSchur "rram26_stage2_block_${blockSize}_validation" 1.0e-5 $blockSize $false 0 $true
            }
            $schur = Import-Csv (Join-Path $directory 'schur_solver_summary.csv')
            $validationRows += [pscustomobject]@{
                block_size = $blockSize
                proxy_colors = $schur.proxy_colors
                block_exact_vectors = $schur.proxy_probing_schur_applies
                scalar_validation_vectors = $schur.proxy_validation_schur_applies
                maximum_matrix_difference = $schur.proxy_block_maximum_difference
                relative_matrix_difference = $schur.proxy_block_relative_difference
                proxy_nnz = $schur.proxy_nnz
                status = $schur.status
            }
        }
        $validationRows | Export-Csv `
            (Join-Path $outputsRoot 'schur_stage2_block_equivalence.csv') -NoTypeInformation
    }
}

if ($Mode -eq 'multi-rhs' -or $Mode -eq 'all') {
    $directory = Join-Path $resultsRoot 'rram26_stage2_multi_rhs_100'
    if (-not $SummarizeOnly) {
        $directory = Invoke-RramSchur 'rram26_stage2_multi_rhs_100' 1.0e-5 64 $false 100 $false
    }
    Copy-Item -LiteralPath (Join-Path $directory 'schur_multi_rhs_summary.csv') `
        -Destination (Join-Path $outputsRoot 'schur_stage2_multi_rhs_summary.csv') -Force
    Copy-Item -LiteralPath (Join-Path $directory 'schur_multi_rhs_runs.csv') `
        -Destination (Join-Path $outputsRoot 'schur_stage2_multi_rhs_runs.csv') -Force
}
