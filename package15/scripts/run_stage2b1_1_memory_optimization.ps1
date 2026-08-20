param(
    [string]$ModelDirectory = 'results\stage2b1_rram26_r130_model\model',
    [string]$OutputDirectory = 'results\stage2b1_1_memory_optimization',
    [string]$ReportDirectory = 'outputs',
    [int]$Seed = 20260721,
    [switch]$SkipBuild
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$model = [System.IO.Path]::GetFullPath((Join-Path $project $ModelDirectory))
$output = [System.IO.Path]::GetFullPath((Join-Path $project $OutputDirectory))
$report = [System.IO.Path]::GetFullPath((Join-Path $project $ReportDirectory))
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'

function Invoke-Checked([string[]]$Arguments) {
    & $exe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "SIPGHeatDDM3D failed with exit code $LASTEXITCODE"
    }
}

Push-Location $project
try {
    if (-not $SkipBuild) {
        cmake -S . -B build
        if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed.' }
        cmake --build build --config Release --parallel
        if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
        ctest --test-dir build -C Release --output-on-failure
        if ($LASTEXITCODE -ne 0) { throw 'Stage 2B.1 regression tests failed.' }
    }
    if (-not (Test-Path -LiteralPath (Join-Path $model 'model.bin'))) {
        throw "Serialized rank-130 model not found: $model"
    }

    New-Item -ItemType Directory -Path $output -Force | Out-Null
    New-Item -ItemType Directory -Path $report -Force | Out-Null
    $rows = @()
    foreach ($count in @(10, 50, 100)) {
        $caseOutput = Join-Path $output ("rram26_batch_{0}" -f $count)
        Invoke-Checked (@('--mor-parametric-deployment-only',
            '--mor-parametric-load', $model,
            '--mor-deployment-rhs-count', [string]$count,
            '--mor-seed', [string]$Seed,
            '--output-dir', $caseOutput))
        $summary = Import-Csv (Join-Path $caseOutput 'stage2b1_1_batch_summary.csv')
        $rows += [pscustomobject]@{
            deployment_cases = $summary.deployment_cases
            load_seconds = $summary.load_seconds
            average_online_seconds = $summary.optimized_average_online_seconds
            batch_wall_seconds = $summary.batch_wall_seconds
            peak_working_set_mib = $summary.optimized_peak_working_set_mib
            final_private_mib = $summary.final_private_mib
            model_resident_mib = $summary.model_resident_mib
            workspace_mib = $summary.workspace_mib
            reconstruction_mib = $summary.reconstruction_mib
            result_cache_mib = $summary.result_cache_mib
            status = $summary.status
        }
    }
    $rows | Export-Csv (Join-Path $report 'stage2b1_1_batch_summary.csv') -NoTypeInformation
    Copy-Item -LiteralPath (Join-Path $output `
        'rram26_batch_100\stage2b1_1_memory_timeline.csv') `
        -Destination (Join-Path $report 'stage2b1_1_memory_timeline.csv') -Force

    $hundred = $rows | Where-Object { [int]$_.deployment_cases -eq 100 }
    if ([double]$hundred.peak_working_set_mib -ge 1024.0) {
        throw "100-case peak memory target failed: $($hundred.peak_working_set_mib) MiB"
    }
    if ([double]$hundred.average_online_seconds -gt 0.09) {
        throw "100-case throughput target failed: $($hundred.average_online_seconds) s/case"
    }
    Write-Host ("Stage 2B.1.1 PASS: {0} MiB peak, {1} s/case" -f `
        $hundred.peak_working_set_mib, $hundred.average_online_seconds)
} finally {
    Pop-Location
}
