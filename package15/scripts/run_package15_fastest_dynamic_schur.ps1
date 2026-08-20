[CmdletBinding()]
param(
    [ValidateSet('smoke', 'medium', 'large')]
    [string]$Profile = 'large',
    [ValidateSet('Cold', 'Warm', 'Both')]
    [string]$Stage = 'Cold',
    [string]$Exe = '',
    [string]$Root = '',
    [ValidateRange(1, 1000)][int]$ValidationSteps = 1,
    [ValidateRange(1, 1000000)][int]$WarmSteps = 100,
    [ValidateRange(1.0e-12, 1.0e6)][double]$Dt = 0.05,
    [switch]$ValidateWithFom,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$driver = Join-Path $PSScriptRoot 'run_package15_dynamic_schur_benchmark.ps1'

if ([string]::IsNullOrWhiteSpace($Root)) {
    $Root = Join-Path $project "build\package15_fastest_$Profile"
}

# This is the measured single-node optimum for the current workstation:
# - the driver bounds independent construction/local work to eight workers;
# - each concurrent local Arnoldi factor receives one MKL thread;
# - augmented-direct receives all 16 outer threads for its one global factor.
$fixed = @{
    Profile = $Profile
    ArnoldiMoments = 1
    ConstructionTraces = 'global-fom'
    ConstructionTraceRank = 19
    SecondMomentMaxColumns = 0
    RankTolerance = 1.0e-6
    EnrichmentRounds = 0
    ValidationSteps = $ValidationSteps
    WarmSteps = $WarmSteps
    Dt = $Dt
    OuterThreads = 16
    MklThreads = 1
    InterfaceSolver = 'augmented-direct'
    Restart = 100
    MaxIterations = 500
    Root = $Root
}
if (-not [string]::IsNullOrWhiteSpace($Exe)) {
    $fixed['Exe'] = $Exe
}
if (-not $ValidateWithFom) {
    # Production timing excludes the optional summary-only monolithic validator.
    $fixed['NoFomComparison'] = $true
}
if ($Force) {
    $fixed['Force'] = $true
}

$stages = if ($Stage -eq 'Both') { @('Cold', 'Warm') } else { @($Stage) }
foreach ($currentStage in $stages) {
    Write-Host (
        "Package15 fastest preset: profile=$Profile stage=$currentStage " +
        'construction=global-fom local=M1 interface=augmented-direct ' +
        'workers=8 factor_threads=16 local_mkl_threads=1 fallback=off')
    $fixed['Stage'] = $currentStage
    & $driver @fixed
}
