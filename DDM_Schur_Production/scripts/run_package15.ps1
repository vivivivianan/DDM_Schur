[CmdletBinding()]
param(
    [ValidateSet('smoke', 'medium', 'large')]
    [string]$Profile = 'large',
    [ValidateSet('Cold', 'Warm', 'Both')]
    [string]$Stage = 'Cold',
    [ValidateRange(1, 1000000)][int]$Steps = 1,
    [ValidateRange(1.0e-12, 1.0e6)][double]$Dt = 0.05,
    [ValidateRange(1, 1024)][int]$Threads = 16,
    [ValidateRange(1, 1024)][int]$LocalMklThreads = 1,
    [switch]$ValidateWithFom,
    [switch]$Force
)

# Supported user-facing launcher. It maps a mesh profile and cold/warm stage to
# the deliberately small executable CLI, protects output deletion to runs/,
# and prints selected timing/cache fields from the authoritative summary CSV.
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\DynamicSchurProduction.exe'
$config = Join-Path $project "configs\package15_$Profile.txt"
$cache = Join-Path $project "cache\package15_$Profile"

if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
    throw "Release executable not found. Build the project first: $exe"
}

function Reset-RunDirectory([string]$Path) {
    $trimSeparators = [char[]]@(
        [IO.Path]::DirectorySeparatorChar,
        [IO.Path]::AltDirectorySeparatorChar)
    $runsRoot = [IO.Path]::GetFullPath(
        (Join-Path $project 'runs')).TrimEnd($trimSeparators)
    $target = [IO.Path]::GetFullPath($Path).TrimEnd($trimSeparators)
    $prefix = $runsRoot + [IO.Path]::DirectorySeparatorChar
    if (-not $target.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to reset a path outside the project runs directory: $target"
    }
    if (Test-Path -LiteralPath $target) {
        Remove-Item -LiteralPath $target -Recurse -Force
    }
    New-Item -ItemType Directory -Force -Path $target | Out-Null
}

function Invoke-Stage([string]$Mode) {
    $output = Join-Path $project (
        "runs\package15_${Profile}_${Mode}_${Steps}steps")
    if ($Force) {
        Reset-RunDirectory $output
    } elseif (Test-Path -LiteralPath $output) {
        if (-not (Get-ChildItem -LiteralPath $output -Force | Select-Object -First 1)) {
            # The executable accepts an existing empty output directory.
        } else {
            throw "Output already exists and is nonempty: $output. Use -Force."
        }
    }

    $arguments = @(
        '--config', $config,
        '--output', $output,
        '--cache', $cache,
        '--mode', $Mode.ToLowerInvariant(),
        '--steps', $Steps,
        '--dt', $Dt.ToString('R', [Globalization.CultureInfo]::InvariantCulture),
        '--threads', $Threads,
        '--local-mkl', $LocalMklThreads)
    if ($ValidateWithFom) {
        $arguments += '--validate-fom'
    }

    Push-Location $project
    try {
        & $exe @arguments
        if ($LASTEXITCODE -ne 0) {
            throw "DynamicSchurProduction failed with exit code $LASTEXITCODE."
        }
    } finally {
        Pop-Location
    }

    $summary = Import-Csv (
        Join-Path $output 'local_dynamic_schur_summary.csv')
    $summary | Select-Object status,global_dofs,full_interface_dofs,
        total_local_rank,steps,descriptor_cache_hit,reference_cache_hit,
        local_model_cache_hit,construction_trace_setup_seconds,
        dynamic_schur_setup_seconds,time_stepping_seconds,
        local_online_core_seconds,total_seconds,maximum_full_residual |
        Format-List
}

$requestedStages = if ($Stage -eq 'Both') {
    @('cold', 'warm')
} else {
    @($Stage.ToLowerInvariant())
}
foreach ($requested in $requestedStages) {
    Invoke-Stage $requested
}
