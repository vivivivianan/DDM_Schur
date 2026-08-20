param(
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [string]$ResultsDirectory =
        'results\milestone8_adaptive_port_basis_rram26',
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$configPath = (Resolve-Path -LiteralPath $Config).Path
$result = if ([IO.Path]::IsPathRooted($ResultsDirectory)) {
    [IO.Path]::GetFullPath($ResultsDirectory)
} else {
    [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite M8.9 basis directory: $result"
}
if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
    throw "Release executable not found: $exe"
}

$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
New-Item -ItemType Directory -Path $result | Out-Null

Push-Location $project
try {
    & $exe --transient --config $configPath `
        --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --port-basis-method hybrid-randomized `
        --mor-arnoldi-moments 1 `
        --randomized-port-rank 16 `
        --randomized-oversampling 5 `
        --randomized-power-iterations 1 `
        --randomized-seed 12345 `
        --residual-krylov-max-rank 4 `
        --residual-krylov-max-sweeps 2 `
        --residual-krylov-tol 1e-4 `
        --residual-krylov-block-size 4 `
        --residual-krylov-probe-mode operator-geometry `
        --residual-krylov-inner-solver woodbury-exact `
        --optimal-port-source-mode trace-only `
        --optimal-port-inner-solver woodbury-exact `
        --optimal-port-inner-tol 1e-10 `
        --optimal-port-inner-refinement-max-iters 3 `
        --optimal-port-inner-refinement-tol 1e-10 `
        --history-compression-method deterministic-rrqr `
        --history-compression-rank 64 `
        --history-compression-tolerance 1e-12 `
        --milestone8-adaptive-production `
        --milestone8-production-basis-only `
        --mor-transient-dt 0.01 `
        --mor-transient-t-end 0.01 `
        --mor-transient-waveform single_step `
        --mor-transient-initial-mode ambient `
        --mor-transient-save $result `
        --output-dir $result --fast-run
    if ($LASTEXITCODE -ne 0) {
        throw 'M8.9 adaptive all-interface basis build failed.'
    }

    $summaryPath = Join-Path $result `
        'milestone8_adaptive_port_basis_summary.csv'
    $diagnosticsPath = Join-Path $result `
        'milestone8_adaptive_port_interface_diagnostics.csv'
    $summary = Import-Csv -LiteralPath $summaryPath |
        Select-Object -First 1
    $diagnostics = @(Import-Csv -LiteralPath $diagnosticsPath)
    if ($summary.status -ne 'passed' -or
        [int]$summary.physical_interfaces -ne 25 -or
        [int]$summary.adaptive_interfaces -ne 5 -or
        [int]$summary.default_interfaces -ne 20 -or
        [int]$summary.snapshot_used -ne 0 -or
        [int]$summary.fom_used_for_basis -ne 0 -or
        [int]$summary.pod_used -ne 0 -or
        [int]$summary.svd_used -ne 0 -or
        [int]$summary.full_field_read -ne 0 -or
        [int]$summary.transient_advanced -ne 0 -or
        $diagnostics.Count -ne 25 -or
        @($diagnostics | Where-Object status -ne 'passed').Count -ne 0) {
        throw 'M8.9 adaptive basis provenance/resource gate failed.'
    }
    if (-not (Test-Path -LiteralPath (
            Join-Path $result 'local_port_basis.bin') -PathType Leaf)) {
        throw 'Serialized M8.9 adaptive basis is missing.'
    }
    Copy-Item -LiteralPath $summaryPath -Destination (
        Join-Path $project `
            'outputs\milestone8_adaptive_port_basis_summary.csv')
    Copy-Item -LiteralPath $diagnosticsPath -Destination (
        Join-Path $project `
            'outputs\milestone8_adaptive_port_interface_diagnostics.csv')
} finally {
    Pop-Location
}
