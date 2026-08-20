param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [IO.Path]::GetFullPath($Root)
$trimSeparators = [char[]]@(
    [IO.Path]::DirectorySeparatorChar, [IO.Path]::AltDirectorySeparatorChar)
$buildRoot = [IO.Path]::GetFullPath((Join-Path $project 'build')).TrimEnd(
    $trimSeparators)
$prefix = $buildRoot + [IO.Path]::DirectorySeparatorChar
if (-not $rootFull.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
if (Test-Path -LiteralPath $rootFull) {
    Remove-Item -LiteralPath $rootFull -Recurse -Force
}

Push-Location $project
try {
    & (Join-Path $project 'scripts\run_package15_fastest_dynamic_schur.ps1') `
        -Profile smoke -Stage Cold -Exe $Exe -Root $rootFull `
        -ValidationSteps 1 -Force | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw 'Package15 fastest production preset failed.'
    }

    $output = Join-Path $rootFull 'dynamic\cold_validation'
    $summary = Import-Csv (
        Join-Path $output 'local_dynamic_schur_summary.csv')
    if ($summary.status -ne 'success' -or
        $summary.construction_trace_mode -ne 'global-fom' -or
        $summary.interface_solver -ne 'augmented-pardiso' -or
        $summary.interface_krylov_requested -ne 'augmented-direct' -or
        $summary.interface_krylov_actual -ne 'augmented-direct' -or
        [int]$summary.local_solve_threads -ne 16 -or
        [int]$summary.local_pardiso_threads -ne 1 -or
        [int]$summary.interface_factor_threads -ne 16 -or
        [int]$summary.analytic_reference_used -ne 1 -or
        [int]$summary.global_construction_factor_used -ne 1 -or
        [int]$summary.residual_fallback_enabled -ne 0 -or
        [int]$summary.residual_fallback_steps -ne 0 -or
        [int]$summary.port_reduction -ne 0 -or
        [int]$summary.fom_comparison_enabled -ne 0 -or
        [int]$summary.steps -ne 1) {
        throw 'Package15 fastest preset changed its validated algorithm contract.'
    }
    if (Test-Path -LiteralPath (
            Join-Path $output 'local_dynamic_schur_final_temperature.csv')) {
        throw 'Package15 fastest preset unexpectedly wrote a full temperature field.'
    }

    & (Join-Path $project 'scripts\run_package15_fastest_dynamic_schur.ps1') `
        -Profile smoke -Stage Warm -Exe $Exe -Root $rootFull `
        -WarmSteps 2 | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw 'Package15 fastest warm-cache preset failed.'
    }
    $warmOutput = Join-Path $rootFull 'dynamic\warm_2_steps'
    $warm = Import-Csv (
        Join-Path $warmOutput 'local_dynamic_schur_summary.csv')
    if ($warm.status -ne 'success' -or
        [int]$warm.descriptor_cache_hit -ne 1 -or
        [int]$warm.reference_cache_hit -ne 1 -or
        [int]$warm.local_model_cache_hit -ne 1 -or
        [double]$warm.construction_trace_setup_seconds -ne 0.0 -or
        [int]$warm.global_construction_factor_used -ne 0 -or
        $warm.interface_krylov_actual -ne 'augmented-direct' -or
        [int]$warm.interface_factor_threads -ne 16 -or
        [int]$warm.residual_fallback_steps -ne 0) {
        throw 'Package15 fastest warm-cache contract failed.'
    }
    if (Test-Path -LiteralPath (
            Join-Path $warmOutput 'local_dynamic_schur_final_temperature.csv')) {
        throw 'Package15 fastest warm preset unexpectedly wrote a full field.'
    }
} finally {
    Pop-Location
}
