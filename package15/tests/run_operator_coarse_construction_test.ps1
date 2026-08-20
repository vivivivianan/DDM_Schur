param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [IO.Path]::GetFullPath($Root)
$buildRoot = [IO.Path]::GetFullPath((Join-Path $project 'build'))
$prefix = $buildRoot.TrimEnd('\', '/') + [IO.Path]::DirectorySeparatorChar
if (-not $rootFull.StartsWith(
        $prefix, [StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
if (Test-Path -LiteralPath $rootFull) {
    Remove-Item -LiteralPath $rootFull -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

$model = Join-Path $rootFull 'model'
$coldOutput = Join-Path $rootFull 'cold'
$warmOutput = Join-Path $rootFull 'warm'
$common = @(
    '--transient', '--config', 'configs\two_cube_parametric_h.txt',
    '--mor-transient-generate', '--mor-transient-method',
    'local-block-arnoldi',
    '--mor-construction-traces', 'operator-coarse',
    '--mor-interface-rank', '19',
    '--mor-arnoldi-moments', '3',
    '--mor-transient-dt', '0.1', '--mor-transient-t-end', '0.2',
    '--mor-transient-waveform', 'single_step',
    '--mor-transient-initial-mode', 'ambient',
    '--mor-interface-krylov', 'augmented-direct',
    '--local-mor-matrix-free-threshold', '0',
    '--schur-local-solve-threads', '2',
    '--schur-local-pardiso-threads', '1',
    '--max-pcg-iterations', '500', '--pcg-tolerance', '1e-10',
    '--no-mor-full-residual-fallback', '--fast-run'
)

Push-Location $project
try {
    & $Exe @common --mor-transient-save $model `
        --mor-transient-production --mor-transient-compare-fom-summary-only `
        --output-dir $coldOutput | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw 'Operator-coarse cold construction failed.'
    }
    $cold = Import-Csv (
        Join-Path $coldOutput 'local_dynamic_schur_summary.csv')
    if ($cold.status -ne 'success' -or
        $cold.construction_trace_mode -ne 'operator-coarse' -or
        [int]$cold.analytic_reference_used -ne 1 -or
        [int]$cold.global_construction_factor_used -ne 0 -or
        [int]$cold.operator_coarse_trace_dimension -le 0 -or
        [int]$cold.operator_trace_krylov_iterations -le 0 -or
        [double]$cold.operator_trace_krylov_maximum_relative_residual `
            -ge 1e-9 -or
        [double]$cold.space_time_relative_l2 -ge 1e-8 -or
        [double]$cold.maximum_temperature_error_k -ge 1e-5 -or
        [double]$cold.maximum_full_residual -ge 1e-5 -or
        [int]$cold.residual_fallback_steps -ne 0) {
        throw 'Operator-coarse cold construction accuracy/provenance gate failed.'
    }

    & $Exe @common --mor-transient-load $model `
        --mor-transient-production --output-dir $warmOutput | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw 'Operator-coarse warm cache run failed.'
    }
    $warm = Import-Csv (
        Join-Path $warmOutput 'local_dynamic_schur_summary.csv')
    if ($warm.status -ne 'success' -or
        [int]$warm.descriptor_cache_hit -ne 1 -or
        [int]$warm.reference_cache_hit -ne 1 -or
        [int]$warm.local_model_cache_hit -ne 1 -or
        [double]$warm.construction_trace_setup_seconds -ne 0.0 -or
        [int]$warm.global_construction_factor_used -ne 0) {
        throw 'Operator-coarse warm cache bypass gate failed.'
    }
} finally {
    Pop-Location
}
