[CmdletBinding()]
param(
    [string]$Config = 'E:\tsv_pdn4_ddm32\ddm_input.txt',
    [string]$Root = 'E:\tsv_pdn4_ddm32_augmented_direct_cn',
    [ValidateSet('Cold', 'Warm', 'Both')]
    [string]$Stage = 'Cold',
    [ValidateRange(1, 1000)][int]$ValidationSteps = 1,
    [ValidateRange(1, 1000000)][int]$WarmSteps = 200,
    [ValidateRange(1, 20)][int]$ArnoldiMoments = 3,
    [ValidateRange(1.0e-15, 1.0)][double]$Dt = 1.0e-9,
    [ValidateRange(1.0e-12, 1.0e6)][double]$FullResidualTolerance = 1.0,
    [ValidateSet('standard', 'shift-invert')][string]$Basis = 'standard',
    [double]$BasisShift = 0.0,
    [switch]$BoundaryAwareBasis,
    [ValidateSet('max-temperature', 'full-field')][string]$OutputMode = 'max-temperature',
    [switch]$ValidateWithFom,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$configFull = [IO.Path]::GetFullPath($Config)
$rootFull = [IO.Path]::GetFullPath($Root)
if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) { throw "Missing executable: $exe" }
if (-not (Test-Path -LiteralPath $configFull -PathType Leaf)) { throw "Missing config: $configFull" }
if ($Basis -eq 'shift-invert' -and -not ($BasisShift -gt 0.0)) {
    throw "BasisShift must be positive for shift-invert."
}

function Reset-Child([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    $prefix = $rootFull.TrimEnd('\') + '\'
    if (-not $full.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to reset outside case root: $full"
    }
    if (Test-Path -LiteralPath $full) { Remove-Item -LiteralPath $full -Recurse -Force }
    New-Item -ItemType Directory -Force -Path $full | Out-Null
}

New-Item -ItemType Directory -Force -Path $rootFull | Out-Null
$waveform = Join-Path $rootFull 'tsv_constant_heat_sources.csv'
& (Join-Path $PSScriptRoot 'generate_tsv_constant_heat_waveform.ps1') `
    -Config $configFull -Output $waveform -EndTime ($Dt * [Math]::Max($ValidationSteps, $WarmSteps))

$cache = Join-Path $rootFull 'local_dynamic_model'
function Invoke-Case([string]$Name, [int]$Steps, [bool]$LoadCache) {
    $output = Join-Path $rootFull $Name
    if ($Force) { Reset-Child $output } else { New-Item -ItemType Directory -Force -Path $output | Out-Null }
    $args = [Collections.Generic.List[string]]@(
        '--transient', '--config', $configFull,
        '--mor-transient-generate', '--mor-transient-method', 'local-block-arnoldi',
        '--mor-arnoldi-moments', $ArnoldiMoments.ToString(),
        '--mor-construction-traces', 'global-fom', '--mor-interface-rank', '19',
        '--mor-arnoldi-rank-tolerance', '1e-6',
        '--mor-basis', $Basis,
        '--mor-transient-dt', $Dt.ToString('R', [Globalization.CultureInfo]::InvariantCulture),
        '--mor-transient-t-end', ($Dt * $Steps).ToString('R', [Globalization.CultureInfo]::InvariantCulture),
        '--mor-transient-integrator', 'crank-nicolson',
        '--mor-transient-input', $waveform,
        '--mor-transient-output', $OutputMode,
        '--mor-transient-initial-mode', 'uniform', '--mor-transient-initial-temperature', '293.15',
        '--mor-transient-production',
        '--mor-interface-krylov', 'augmented-direct',
        '--local-mor-matrix-free-threshold', '0',
        '--schur-local-solve-threads', '8', '--schur-local-pardiso-threads', '16',
        '--no-mor-native-reduced-history', '--no-mor-full-residual-fallback',
        '--mor-full-residual-tolerance', $FullResidualTolerance.ToString('R', [Globalization.CultureInfo]::InvariantCulture),
        '--output-dir', $output, '--fast-run'
    )
    if ($Basis -eq 'shift-invert') {
        [void]$args.Add('--mor-basis-shift')
        [void]$args.Add($BasisShift.ToString('R', [Globalization.CultureInfo]::InvariantCulture))
    }
    if ($BoundaryAwareBasis) { [void]$args.Add('--mor-boundary-aware-basis') }
    if ($LoadCache) { [void]$args.Add('--mor-transient-load'); [void]$args.Add($cache) }
    else { [void]$args.Add('--mor-transient-save'); [void]$args.Add($cache) }
    if ($ValidateWithFom) { [void]$args.Add('--mor-transient-compare-fom-summary-only') }
    $env:OMP_NUM_THREADS = '8'; $env:MKL_NUM_THREADS = '16'; $env:MKL_DYNAMIC = 'FALSE'; $env:SIPG_SOLVER_WORKERS = '8'
    & $exe @args 2>&1 | Tee-Object -FilePath (Join-Path $output 'run.log')
    if ($LASTEXITCODE -ne 0) { throw "Augmented-direct CN run failed: exit code $LASTEXITCODE" }
}

if ($Stage -in @('Cold', 'Both')) { Invoke-Case 'cold_validation' $ValidationSteps $false }
if ($Stage -in @('Warm', 'Both')) {
    if (-not (Test-Path -LiteralPath (Join-Path $cache 'local_dynamic_interior_model.bin') -PathType Leaf)) {
        throw "Warm stage requires cold cache: $cache"
    }
    Invoke-Case "warm_$($WarmSteps)_steps" $WarmSteps $true
}
