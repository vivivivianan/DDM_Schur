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

function Invoke-Path(
    [string]$Name,
    [int]$Threshold,
    [string]$InterfaceSolver = 'fgmres',
    [string]$PortCoreCache = '') {
    $output = Join-Path $rootFull $Name
    $cacheArguments = if ([string]::IsNullOrWhiteSpace($PortCoreCache)) {
        @()
    } else {
        @('--mor-port-core-cache', $PortCoreCache)
    }
    & $Exe --transient --config configs\two_cube_parametric_h.txt `
        --mor-transient-generate --mor-transient-method local-block-arnoldi `
        --mor-arnoldi-moments 3 --mor-transient-dt 0.1 --mor-transient-t-end 0.2 `
        --mor-transient-waveform single_step --mor-transient-initial-mode ambient `
        --mor-interface-krylov $InterfaceSolver @cacheArguments `
        --local-mor-matrix-free-threshold $Threshold --max-pcg-iterations 200 `
        --gmres-restart 30 --pcg-tolerance 1e-10 --schur-proxy-ring 1 `
        --schur-proxy-block-size 8 --output-dir $output --fast-run
    if ($LASTEXITCODE -ne 0) { throw "$Name Dynamic Schur path failed." }
    return Import-Csv (Join-Path $output 'local_dynamic_schur_summary.csv')
}

Push-Location $project
try {
    $matrixFree = Invoke-Path 'matrix_free' 0
    $explicit = Invoke-Path 'explicit' 20000
    $augmented = Invoke-Path 'augmented_direct' 0 'augmented-direct'
    $portCoreCache = Join-Path $rootFull 'port_core_artifacts.bin'
    if (Test-Path -LiteralPath $portCoreCache -PathType Leaf) {
        Remove-Item -LiteralPath $portCoreCache -Force
    }
    $portCore = Invoke-Path 'port_core' 0 'port-core' $portCoreCache
    $portCoreCached = Invoke-Path `
        'port_core_cached' 0 'port-core' $portCoreCache
    if ($matrixFree.status -ne 'success' -or
        $matrixFree.interface_solver -ne 'matrix-free-fgmres') {
        throw 'Matrix-free Dynamic Schur did not report success.'
    }
    if ([int]$matrixFree.interface_iterations_maximum -le 0 -or
        [int]$matrixFree.interface_iterations_maximum -ge 100) {
        throw 'Matrix-free Dynamic Schur iteration regression.'
    }
    if ([double]$matrixFree.maximum_interface_relative_residual -ge 1e-9) {
        throw 'Matrix-free Dynamic Schur true residual gate failed.'
    }
    if ([int]$matrixFree.coarse_dimension -ne 8 -or
        [int]$matrixFree.proxy_colors -le 0 -or
        [int]$matrixFree.proxy_nonzeros -le 0) {
        throw 'Matrix-free Dynamic Schur proxy/coarse setup is incomplete.'
    }
    if ($explicit.interface_solver -ne 'sparse-pardiso') {
        throw 'Explicit Dynamic Schur regression path was not retained.'
    }
    if ($augmented.status -ne 'success' -or
        $augmented.interface_solver -ne 'augmented-pardiso' -or
        $augmented.interface_krylov_actual -ne 'augmented-direct' -or
        [int]$augmented.interface_iterations_maximum -ne 1 -or
        [double]$augmented.maximum_interface_relative_residual -ge 1e-9) {
        throw 'Augmented direct reduced solve failed its residual gate.'
    }
    $augmentedSetup = Import-Csv (
        Join-Path $rootFull 'augmented_direct\augmented_direct_summary.csv')
    if ([int]$augmentedSetup.dimension -le [int]$augmented.full_interface_dofs -or
        [int]$augmentedSetup.reduced_dofs -ne [int]$augmented.total_local_rank -or
        [int64]$augmentedSetup.nonzeros -le 0) {
        throw 'Augmented direct setup diagnostics are incomplete.'
    }
    if ($portCore.status -ne 'success' -or
        $portCore.interface_solver -ne 'port-core-direct' -or
        $portCore.interface_krylov_actual -ne 'port-core-direct' -or
        [int]$portCore.interface_iterations_maximum -ne 1 -or
        [double]$portCore.maximum_interface_relative_residual -ge 1e-9) {
        throw 'Parallel port/core direct solve failed its reduced residual gate.'
    }
    if ([double]$portCore.port_forward_solve_seconds -le 0.0 -or
        [double]$portCore.port_core_solve_seconds -le 0.0 -or
        [double]$portCore.port_back_substitution_seconds -le 0.0 -or
        [double]$portCore.proxy_solve_seconds -ne 0.0 -or
        [double]$portCore.coarse_solve_seconds -ne 0.0) {
        throw 'Parallel port/core timing was not reported in dedicated fields.'
    }
    if ([int]$portCore.port_core_cache_hit -ne 0 -or
        [double]$portCore.port_core_cache_save_seconds -le 0.0 -or
        [int]$portCoreCached.port_core_cache_hit -ne 1 -or
        [double]$portCoreCached.port_core_cache_load_seconds -le 0.0 -or
        [double]$portCoreCached.maximum_interface_relative_residual -ge 1e-9) {
        throw 'Exact port/core artifact cache cold/warm gate failed.'
    }
    $portCoreSetup = Import-Csv (
        Join-Path $rootFull 'port_core\port_core_summary.csv')
    if ([int]$portCoreSetup.physical_ports -ne 1 -or
        [int]$portCoreSetup.leaf_dofs -le 0 -or
        [int]$portCoreSetup.core_dimension -le 0) {
        throw 'Parallel port/core partition diagnostics are incomplete.'
    }

    $left = Import-Csv (Join-Path $rootFull 'matrix_free\local_dynamic_schur_final_temperature.csv')
    $right = Import-Csv (Join-Path $rootFull 'explicit\local_dynamic_schur_final_temperature.csv')
    if ($left.Count -ne $right.Count) { throw 'Temperature vector sizes differ.' }
    [double]$error2 = 0.0
    [double]$reference2 = 0.0
    [double]$maximum = 0.0
    for ($row = 0; $row -lt $left.Count; ++$row) {
        $a = [double]$left[$row].temperature_k
        $b = [double]$right[$row].temperature_k
        $difference = $a - $b
        $error2 += $difference * $difference
        $reference2 += $b * $b
        $maximum = [Math]::Max($maximum, [Math]::Abs($difference))
    }
    if ([Math]::Sqrt($error2 / $reference2) -ge 1e-10 -or $maximum -ge 1e-7) {
        throw "Matrix-free/explicit mismatch: relative=$([Math]::Sqrt($error2 / $reference2)), max=$maximum"
    }

    $left = Import-Csv (Join-Path $rootFull 'port_core\local_dynamic_schur_final_temperature.csv')
    if ($left.Count -ne $right.Count) { throw 'Port-core temperature vector size differs.' }
    [double]$error2 = 0.0
    [double]$reference2 = 0.0
    [double]$maximum = 0.0
    for ($row = 0; $row -lt $left.Count; ++$row) {
        $a = [double]$left[$row].temperature_k
        $b = [double]$right[$row].temperature_k
        $difference = $a - $b
        $error2 += $difference * $difference
        $reference2 += $b * $b
        $maximum = [Math]::Max($maximum, [Math]::Abs($difference))
    }
    if ([Math]::Sqrt($error2 / $reference2) -ge 1e-10 -or $maximum -ge 1e-7) {
        throw "Port-core/explicit mismatch: relative=$([Math]::Sqrt($error2 / $reference2)), max=$maximum"
    }
    $portCoreTemperature = $left
    $left = Import-Csv (Join-Path $rootFull 'augmented_direct\local_dynamic_schur_final_temperature.csv')
    if ($left.Count -ne $right.Count) {
        throw 'Augmented-direct temperature vector size differs.'
    }
    [double]$error2 = 0.0
    [double]$reference2 = 0.0
    [double]$maximum = 0.0
    for ($row = 0; $row -lt $left.Count; ++$row) {
        $a = [double]$left[$row].temperature_k
        $b = [double]$right[$row].temperature_k
        $difference = $a - $b
        $error2 += $difference * $difference
        $reference2 += $b * $b
        $maximum = [Math]::Max($maximum, [Math]::Abs($difference))
    }
    if ([Math]::Sqrt($error2 / $reference2) -ge 1e-10 -or $maximum -ge 1e-7) {
        throw "Augmented-direct/explicit mismatch: relative=$([Math]::Sqrt($error2 / $reference2)), max=$maximum"
    }
    $cached = Import-Csv (
        Join-Path $rootFull 'port_core_cached\local_dynamic_schur_final_temperature.csv')
    if ($cached.Count -ne $portCoreTemperature.Count) {
        throw 'Cached port-core temperature vector size differs.'
    }
    [double]$error2 = 0.0
    [double]$reference2 = 0.0
    for ($row = 0; $row -lt $portCoreTemperature.Count; ++$row) {
        $a = [double]$cached[$row].temperature_k
        $b = [double]$portCoreTemperature[$row].temperature_k
        $difference = $a - $b
        $error2 += $difference * $difference
        $reference2 += $b * $b
    }
    if ([Math]::Sqrt($error2 / $reference2) -ge 1e-13) {
        throw 'Cached and cold port-core trajectories differ.'
    }
    $flux = Import-Csv (Join-Path $rootFull 'matrix_free\local_dynamic_schur_interface_flux.csv')
    if (@($flux).Count -le 0 -or
        -not $flux[0].PSObject.Properties['fom_sipg_numerical_flux_w_m2']) {
        throw 'Detailed transient FOM/ROM interface flux output is missing.'
    }
} finally {
    Pop-Location
}
