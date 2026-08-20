param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root,
    [Parameter(Mandatory = $true)][string]$BinaryRoot
)

# End-to-end production contract test. It deliberately exercises one cold step
# and two warm steps: cold checks numerical residual and cache creation; warm
# checks strict cache reuse and multi-step control flow. Full multi-step
# accuracy is validated separately against summary-only FOM.
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [IO.Path]::GetFullPath($Root)
$trimSeparators = [char[]]@(
    [IO.Path]::DirectorySeparatorChar,
    [IO.Path]::AltDirectorySeparatorChar)
$buildRoot = [IO.Path]::GetFullPath($BinaryRoot).TrimEnd($trimSeparators)
$prefix = $buildRoot + [IO.Path]::DirectorySeparatorChar
if (-not $rootFull.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest root must remain under build: $rootFull"
}
if (Test-Path -LiteralPath $rootFull) {
    Remove-Item -LiteralPath $rootFull -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

$cache = Join-Path $rootFull 'cache'
$cold = Join-Path $rootFull 'cold'
$warm = Join-Path $rootFull 'warm'
$config = Join-Path $project 'configs\package15_smoke.txt'

Push-Location $project
try {
    & $Exe --config $config --output $cold --cache $cache --mode cold `
        --steps 1 --dt 0.05 --threads 16 --local-mkl 1 | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'Cold production smoke failed.' }

    $coldSummary = Import-Csv (
        Join-Path $cold 'local_dynamic_schur_summary.csv')
    if ($coldSummary.status -ne 'success' -or
        $coldSummary.construction_trace_mode -ne 'global-fom' -or
        $coldSummary.interface_krylov_actual -ne 'augmented-direct' -or
        [int]$coldSummary.interface_factor_threads -ne 16 -or
        [int]$coldSummary.local_pardiso_threads -ne 1 -or
        [int]$coldSummary.residual_fallback_enabled -ne 0 -or
        [int]$coldSummary.fom_comparison_enabled -ne 0 -or
        [int]$coldSummary.port_reduction -ne 0 -or
        [double]$coldSummary.maximum_full_residual -gt 1.0e-4) {
        throw 'Cold run violated the production algorithm contract.'
    }

    & $Exe --config $config --output $warm --cache $cache --mode warm `
        --steps 2 --dt 0.05 --threads 16 --local-mkl 1 | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'Warm production smoke failed.' }

    $warmSummary = Import-Csv (
        Join-Path $warm 'local_dynamic_schur_summary.csv')
    if ($warmSummary.status -ne 'success' -or
        [int]$warmSummary.descriptor_cache_hit -ne 1 -or
        [int]$warmSummary.reference_cache_hit -ne 1 -or
        [int]$warmSummary.local_model_cache_hit -ne 1 -or
        [double]$warmSummary.construction_trace_setup_seconds -ne 0.0 -or
        [int]$warmSummary.global_construction_factor_used -ne 0 -or
        $warmSummary.interface_krylov_actual -ne 'augmented-direct') {
        throw 'Warm run violated the cache/solver contract.'
    }

    foreach ($output in @($cold, $warm)) {
        if (Test-Path -LiteralPath (
                Join-Path $output 'local_dynamic_schur_final_temperature.csv')) {
            throw 'Production run unexpectedly wrote a full field.'
        }
    }
} finally {
    Pop-Location
}
