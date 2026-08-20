param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [System.IO.Path]::GetFullPath($Root)
$buildRoot = [System.IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

function Invoke-Case([string]$Name, [string[]]$Extra) {
    $output = Join-Path $rootFull $Name
    $arguments = @(
        '--transient', '--config', 'configs\two_cube_parametric_h.txt',
        '--mor-transient-generate', '--mor-transient-method', 'local-block-arnoldi',
        '--mor-arnoldi-moments', '3', '--mor-transient-dt', '0.1',
        '--mor-transient-t-end', '0.2', '--mor-transient-waveform', 'single_step',
        '--mor-transient-initial-mode', 'ambient', '--mor-interface-initial-guess',
        'previous', '--local-mor-matrix-free-threshold', '0',
        '--max-pcg-iterations', '200', '--gmres-restart', '30',
        '--pcg-tolerance', '1e-10', '--schur-proxy-ring', '1',
        '--schur-proxy-block-size', '8', '--output-dir', $output, '--fast-run')
    $arguments += $Extra
    & $Exe @arguments | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "$Name failed." }
    return Import-Csv (Join-Path $output 'local_dynamic_schur_summary.csv')
}

Push-Location $project
try {
    $validation = Invoke-Case 'validation' @()
    $production = Invoke-Case 'production' @('--mor-transient-production')
    $fallback = Invoke-Case 'fallback' @(
        '--mor-transient-production', '--mor-full-residual-tolerance', '1e-8',
        '--mor-full-residual-fallback')
    $monitor = Invoke-Case 'monitor_only' @(
        '--mor-transient-production', '--mor-full-residual-tolerance', '1e-8',
        '--no-mor-full-residual-fallback')

    if ($validation.status -ne 'success' -or
        [int]$validation.fom_comparison_enabled -ne 1 -or
        [uint64]$validation.fom_factor_memory_bytes -eq 0) {
        throw 'Validation mode did not retain the FOM comparison.'
    }
    if ($production.status -ne 'success' -or
        [int]$production.fom_comparison_enabled -ne 0 -or
        [uint64]$production.fom_factor_memory_bytes -ne 0 -or
        [int]$production.residual_fallback_steps -ne 0) {
        throw 'Production mode eagerly created or used the FOM factor.'
    }
    if (Test-Path (Join-Path $rootFull 'production\local_dynamic_schur_interface_flux.csv')) {
        throw 'Production mode wrote validation-only detailed flux output.'
    }
    if (Test-Path (Join-Path $rootFull 'production\local_dynamic_schur_final_temperature.csv')) {
        throw 'Production max-temperature mode wrote the full field.'
    }
    if ($fallback.status -ne 'success' -or
        [int]$fallback.fom_comparison_enabled -ne 0 -or
        [uint64]$fallback.fom_factor_memory_bytes -eq 0 -or
        [int]$fallback.residual_fallback_steps -ne 2 -or
        [double]$fallback.maximum_full_residual -gt 1e-8) {
        throw 'Lazy FOM residual fallback did not satisfy the gate.'
    }
    if ($monitor.status -ne 'success' -or
        [int]$monitor.fom_comparison_enabled -ne 0 -or
        [int]$monitor.residual_fallback_enabled -ne 0 -or
        [uint64]$monitor.fom_factor_memory_bytes -ne 0 -or
        [int]$monitor.residual_fallback_steps -ne 0 -or
        [int]$monitor.residual_tolerance_violation_steps -ne 2 -or
        [int]$monitor.residual_gate_passed -ne 0 -or
        [double]$monitor.maximum_full_residual -le 1e-8) {
        throw 'Monitor-only residual policy invoked FOM or hid a violation.'
    }
    $timing = Import-Csv (Join-Path $rootFull `
        'production\local_dynamic_schur_reduced_timing.csv')
    if (@($timing).Count -ne 2 -or
        $timing[0].interface_initial_guess -ne 'previous') {
        throw 'Production per-step initial-guess timing is incomplete.'
    }
} finally {
    Pop-Location
}
