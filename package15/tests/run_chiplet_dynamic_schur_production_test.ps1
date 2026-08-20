param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root,
    [Parameter(Mandatory = $true)][string]$Config
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [IO.Path]::GetFullPath($Root)
$buildRoot = [IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
if (Test-Path -LiteralPath $rootFull) {
    Remove-Item -LiteralPath $rootFull -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

$model = Join-Path $rootFull 'local_dynamic_model'
$proxy = Join-Path $rootFull 'dynamic_schur.proxycache'
$coarse = Join-Path $rootFull 'operator_coarse_r4.bin'

function Read-Summary([string]$Output) {
    return Import-Csv (Join-Path $Output 'local_dynamic_schur_summary.csv')
}

Push-Location $project
try {
    $coldOutput = Join-Path $rootFull 'cold_enriched'
    & $Exe --transient --config $Config `
        --mor-transient-generate --mor-transient-method local-block-arnoldi `
        --mor-transient-save $model --mor-arnoldi-moments 2 `
        --mor-arnoldi-second-moment-energy 1 `
        --mor-arnoldi-second-moment-max-columns 4 `
        --local-port-enrichment-rounds 11 --mor-transient-dt 0.1 `
        --mor-transient-t-end 1.0 --mor-transient-waveform asynchronous_hotspots `
        --mor-transient-initial-mode ambient --mor-transient-production `
        --mor-transient-compare-fom-summary-only `
        --mor-interface-initial-guess previous --mor-interface-krylov fgmres `
        --local-mor-matrix-free-threshold 0 --max-pcg-iterations 500 `
        --gmres-restart 200 --pcg-tolerance 1e-10 `
        --mor-adaptive-interface-tolerance 1e-10 --schur-proxy-ring 1 `
        --schur-proxy-block-size 64 --schur-local-solve-threads 2 `
        --schur-local-pardiso-threads 1 --schur-proxy-cache $proxy `
        --schur-interface-operator-coarse-rank 4 `
        --schur-interface-operator-coarse-sweeps 2 `
        --schur-interface-operator-coarse-cache $coarse `
        --output-dir $coldOutput --fast-run | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'Chiplet cold enriched run failed.' }
    $cold = Read-Summary $coldOutput

    $warmOutput = Join-Path $rootFull 'warm_reload'
    & (Join-Path $project 'scripts\run_cached_local_dynamic_schur.ps1') `
        -Exe $Exe -Config $Config -ModelCache $model -ProxyCache $proxy `
        -OutputDirectory $warmOutput -Steps 10 -Dt 0.1 -OuterThreads 2 `
        -MklThreads 1 -Restart 200 -MaxIterations 500 `
        -AdaptiveTolerance 1e-10 -OperatorCoarseRank 4 `
        -OperatorCoarseSweeps 2 -OperatorCoarseCache $coarse `
        -ArnoldiMoments 2 -SecondMomentEnergy 1 `
        -SecondMomentMaxColumns 4 -Waveform asynchronous_hotspots `
        -CompareFomSummaryOnly | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'Chiplet warm cache reload failed.' }
    $warm = Read-Summary $warmOutput

    foreach ($summary in @($cold, $warm)) {
        if ($summary.status -ne 'success' -or
            [int]$summary.residual_fallback_steps -ne 0 -or
            [double]$summary.maximum_full_residual -ge 2e-6 -or
            [double]$summary.maximum_interface_relative_residual -ge 1e-9) {
            throw 'Chiplet production residual/fallback gate failed.'
        }
    }
    if ([int]$cold.total_local_rank -ne 137 -or
        [int]$warm.total_local_rank -ne [int]$cold.total_local_rank -or
        [int]$warm.local_model_cache_hit -ne 1 -or
        [int]$warm.descriptor_cache_hit -ne 1 -or
        [int]$warm.reference_cache_hit -ne 1) {
        throw 'Enriched Chiplet local-model cache did not reload intact.'
    }
} finally {
    Pop-Location
}
