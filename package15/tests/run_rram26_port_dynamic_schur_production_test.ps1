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

function Read-Summary([string]$Output) {
    return Import-Csv (Join-Path $Output 'local_dynamic_schur_summary.csv')
}

Push-Location $project
try {
    $coldOutput = Join-Path $rootFull 'cold_port_model'
    & $Exe --transient --config $Config `
        --mor-transient-generate --mor-transient-method local-port-block-arnoldi `
        --mor-transient-save $model --mor-arnoldi-moments 1 `
        --local-port-rank 16 --local-port-energy-tolerance 1e-20 `
        --local-port-enrichment-rounds 1 --mor-transient-dt 0.01 `
        --mor-transient-t-end 0.01 --mor-transient-waveform mixed_frequency `
        --mor-transient-initial-mode ambient --mor-transient-production `
        --mor-transient-compare-fom-summary-only `
        --schur-proxy-cache $proxy --output-dir $coldOutput --fast-run | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'RRAM26 cold port-model run failed.' }
    $cold = Read-Summary $coldOutput

    $warmOutput = Join-Path $rootFull 'warm_port_reload'
    & (Join-Path $project 'scripts\run_cached_local_dynamic_schur.ps1') `
        -Exe $Exe -Config $Config -ModelCache $model -ProxyCache $proxy `
        -OutputDirectory $warmOutput -Method local-port-block-arnoldi `
        -PortRank 16 -PortEnergyTolerance 1e-20 -PortEnrichmentRounds 1 `
        -ArnoldiMoments 1 -Steps 1 -Dt 0.01 -Waveform mixed_frequency `
        -CompareFomSummaryOnly | Out-Null
    if ($LASTEXITCODE -ne 0) { throw 'RRAM26 warm port-model reload failed.' }
    $warm = Read-Summary $warmOutput

    foreach ($summary in @($cold, $warm)) {
        if ($summary.status -ne 'success' -or [int]$summary.port_reduction -ne 1 -or
            [int]$summary.port_dimension -le 0 -or
            [int]$summary.port_dimension -ge [int]$summary.full_interface_dofs -or
            [double]$summary.maximum_full_residual -ge 1e-6) {
            throw 'RRAM26 port production residual/compression gate failed.'
        }
    }
    if (-not (Test-Path -LiteralPath (Join-Path $model 'local_port_basis.bin') `
        -PathType Leaf)) {
        throw 'RRAM26 production cache did not reload the serialized port model.'
    }
} finally {
    Pop-Location
}
