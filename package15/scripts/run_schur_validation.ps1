param(
    [ValidateSet('two', 'ten', 'all')]
    [string]$Case = 'all',
    [string]$Configuration = 'Release'
)

$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$executable = Join-Path $projectRoot "build\$Configuration\SIPGHeatDDM3D.exe"
if (-not (Test-Path -LiteralPath $executable)) {
    throw "Executable not found: $executable"
}

function Invoke-SchurCase([string]$name, [string]$config, [int]$restart) {
    $output = Join-Path $projectRoot "results\${name}_verified"
    & $executable `
        --steady `
        --config (Join-Path $projectRoot $config) `
        --output-dir $output `
        --solvers schur,direct `
        --direct-mode spd `
        --schur-proxy-ring 1 `
        --schur-proxy-use-material-connectivity `
        --schur-proxy-high-k-threshold 20 `
        --schur-proxy-block-size 64 `
        --pcg-tolerance 1e-10 `
        --max-pcg-iterations 500 `
        --gmres-restart $restart
    if ($LASTEXITCODE -ne 0) {
        throw "$name validation failed with exit code $LASTEXITCODE"
    }
}

if ($Case -eq 'two' -or $Case -eq 'all') {
    Invoke-SchurCase 'two_cube' 'configs\two_cube_schur.txt' 766
}
if ($Case -eq 'ten' -or $Case -eq 'all') {
    Invoke-SchurCase 'ten_cube' 'configs\ten_cube_schur.txt' 400
}
