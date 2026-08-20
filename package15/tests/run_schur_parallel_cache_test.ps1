param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [IO.Path]::GetFullPath($Root)
$buildRoot = [IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null
$cache = Join-Path $rootFull 'two_cube.proxycache'
if (Test-Path -LiteralPath $cache) { Remove-Item -LiteralPath $cache -Force }

function Invoke-Schur([string]$Name, [int]$OuterThreads, [double]$HighKThreshold = 100.0) {
    $output = Join-Path $rootFull $Name
    if (Test-Path -LiteralPath $output) { Remove-Item -LiteralPath $output -Recurse -Force }
    & $Exe --steady --config configs\two_cube_schur.txt --solvers schur,direct `
        --schur-local-solve-threads $OuterThreads --schur-local-pardiso-threads 1 `
        --schur-proxy-ring 1 --schur-proxy-high-k-threshold $HighKThreshold `
        --schur-proxy-cache $cache --output-dir $output
    if ($LASTEXITCODE -ne 0) { throw "$Name Schur run failed." }
    return Import-Csv (Join-Path $output 'schur_solver_summary.csv')
}

Push-Location $project
try {
    $miss = Invoke-Schur 'miss_outer1' 1
    $hit = Invoke-Schur 'hit_outer2' 2
    if ([int]$miss.proxy_matrix_cache_hit -ne 0 -or [int]$hit.proxy_matrix_cache_hit -ne 1) {
        throw 'Persistent proxy cache miss/hit sequence was not reported.'
    }
    if ([int]$hit.proxy_probing_schur_applies -ne 0) {
        throw 'A proxy matrix cache hit unexpectedly repeated exact color probing.'
    }
    if ([int]$miss.local_symbolic_analysis_calls -ne 2 -or
        [int]$miss.local_numerical_factorization_calls -ne 2) {
        throw 'Unused full blockSolvers were constructed on the proxy/non-energy path.'
    }
    $domains = Import-Csv (Join-Path $rootFull 'hit_outer2\schur_subdomain_performance.csv')
    $domainSummary = Import-Csv (Join-Path $rootFull 'hit_outer2\schur_subdomain_performance_summary.csv')
    if (@($domains).Count -ne 2 -or [int]$domainSummary.phase33_calls -le 0) {
        throw 'Per-subdomain phase-33 performance output is incomplete.'
    }
    [double]$sum = 0.0
    foreach ($domain in $domains) { $sum += [double]$domain.phase33_total_time }
    if ([Math]::Abs($sum - [double]$domainSummary.phase33_total_time) -gt 1e-10) {
        throw 'Per-domain and aggregate phase-33 times are inconsistent.'
    }
    $comparison = Import-Csv (Join-Path $rootFull 'hit_outer2\solver_comparison.csv')
    $schur = $comparison | Where-Object { $_.solver -eq 'ddm_schur_fgmres' }
    if ([double]$schur.relative_l2_vs_direct -ge 1e-10) {
        throw "Parallel local solves changed the Schur solution: $($schur.relative_l2_vs_direct)"
    }

    $mismatch = Invoke-Schur 'fingerprint_mismatch' 2 9999.0
    if ([int]$mismatch.proxy_matrix_cache_hit -ne 0) {
        throw 'A proxy configuration fingerprint mismatch was incorrectly reused.'
    }
} finally {
    Pop-Location
}
