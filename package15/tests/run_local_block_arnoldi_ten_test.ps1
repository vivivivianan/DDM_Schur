param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('Generate', 'Reuse')]
    [string]$Mode,
    [Parameter(Mandatory = $true)]
    [string]$Exe,
    [Parameter(Mandatory = $true)]
    [string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [System.IO.Path]::GetFullPath($Root)
$buildRoot = [System.IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null
$independent = Join-Path $rootFull 'independent'
$reuse = Join-Path $rootFull 'reuse'

function Invoke-Ten([string]$Output, [bool]$Reuse) {
    $arguments = @(
        '--transient', '--config', 'configs\ten_cube_parametric_h.txt',
        '--mor-transient-generate', '--mor-transient-method', 'local-block-arnoldi',
        '--mor-arnoldi-moments', '2', '--mor-transient-dt', '0.1',
        '--mor-transient-t-end', '0.2', '--mor-transient-waveform', 'single_step',
        '--mor-transient-initial-mode', 'ambient', '--output-dir', $Output, '--fast-run')
    if ($Reuse) { $arguments += '--mor-local-transient-reuse-identical-subdomains' }
    & $Exe @arguments
    if ($LASTEXITCODE -ne 0) { throw 'Ten-cube Local Dynamic Schur execution failed.' }
}

function Assert-Less([double]$Value, [double]$Limit, [string]$Message) {
    if (-not ($Value -lt $Limit)) { throw "$Message (value=$Value, limit=$Limit)" }
}

Push-Location $project
try {
    if ($Mode -eq 'Generate') {
        Invoke-Ten $independent $false
        $summary = Import-Csv (Join-Path $independent 'local_dynamic_schur_summary.csv')
        if ($summary.status -ne 'success' -or [int]$summary.subdomains -ne 10 -or
            [int]$summary.full_interface_dofs -ne 10593) {
            throw 'Ten-cube did not retain ten local bases and the 10,593-DOF interface.'
        }
        Assert-Less ([double]$summary.space_time_relative_l2) 1.0e-5 `
            'Ten-cube space-time accuracy failed'
        Assert-Less ([double]$summary.maximum_absolute_k) 0.1 `
            'Ten-cube maximum error failed'
        if ([int]$summary.dynamic_schur_symbolic_calls -ne 1 -or
            [int]$summary.dynamic_schur_numerical_calls -ne 1) {
            throw 'Ten-cube Dynamic Schur factor was rebuilt across time steps.'
        }
        $rank = Import-Csv (Join-Path $independent 'local_block_arnoldi_rank.csv')
        if (@($rank | Select-Object -ExpandProperty subdomain -Unique).Count -ne 10) {
            throw 'Ten-cube local rank output does not contain ten subdomains.'
        }
    } else {
        if (-not (Test-Path (Join-Path $independent 'local_dynamic_schur_summary.csv'))) {
            throw 'Independent fixture is missing.'
        }
        Invoke-Ten $reuse $true
        $summary = Import-Csv (Join-Path $reuse 'local_dynamic_schur_summary.csv')
        if ([int]$summary.reuse_requested -ne 1 -or [int]$summary.unique_templates -ne 10 -or
            [int]$summary.reused_instances -ne 0) {
            throw 'Strict fingerprint reuse audit did not reject the distinct trace bases.'
        }
        $a = Import-Csv (Join-Path $independent 'local_dynamic_schur_final_temperature.csv')
        $b = Import-Csv (Join-Path $reuse 'local_dynamic_schur_final_temperature.csv')
        if ($a.Count -ne $b.Count) { throw 'Reuse changed the final field size.' }
        $maximum = 0.0
        for ($index = 0; $index -lt $a.Count; ++$index) {
            $maximum = [Math]::Max($maximum,
                [Math]::Abs([double]$a[$index].temperature_k - [double]$b[$index].temperature_k))
        }
        Assert-Less $maximum 1.0e-8 'Reuse audit changed the final temperature'
    }
} finally {
    Pop-Location
}
