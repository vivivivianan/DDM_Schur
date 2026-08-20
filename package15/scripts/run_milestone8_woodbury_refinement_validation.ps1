param(
    [string]$ResultsDirectory =
        'results\milestone8_woodbury_refinement_validation',
    [string]$OutputsDirectory = 'outputs',
    [int]$CpuThreads = 8,
    [switch]$SkipBuild
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$result = [System.IO.Path]::GetFullPath(
    (Join-Path $project $ResultsDirectory))
$outputs = [System.IO.Path]::GetFullPath(
    (Join-Path $project $OutputsDirectory))
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'

if ($CpuThreads -le 0) {
    throw 'CPU thread count must be positive.'
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite refinement validation: $result"
}

Push-Location $project
try {
    if (-not $SkipBuild) {
        cmake -S . -B build
        if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed.' }
        cmake --build build --config Release -j $CpuThreads
        if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
    }
    if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
        throw "Release executable not found: $exe"
    }

    New-Item -ItemType Directory -Path $result | Out-Null
    New-Item -ItemType Directory -Force -Path $outputs | Out-Null
    $env:OMP_NUM_THREADS = "$CpuThreads"
    $env:MKL_NUM_THREADS = "$CpuThreads"
    $env:MKL_DYNAMIC = 'FALSE'

    $allRows = @()
    $cases = @(
        [pscustomobject]@{
            Label = 'two-cube'
            Config = 'configs\two_cube_parametric_h.txt'
        },
        [pscustomobject]@{
            Label = 'ten-cube'
            Config = 'configs\ten_cube_parametric_h.txt'
        })
    foreach ($case in $cases) {
        $caseOutput = Join-Path $result $case.Label
        $arguments = @(
            '--transient','--config',$case.Config,
            '--mor-transient-generate',
            '--mor-transient-method','local-port-block-arnoldi',
            '--port-basis-method','optimal-transfer',
            '--mor-arnoldi-moments','1',
            '--optimal-port-rank','6',
            '--optimal-port-inner-solver','woodbury-exact',
            '--optimal-port-inner-tol','1e-10',
            '--optimal-port-inner-max-iters','1000',
            '--optimal-port-inner-refinement-max-iters','3',
            '--optimal-port-inner-refinement-tol','1e-10',
            '--optimal-port-eigen-max-iters','200',
            '--optimal-port-eigen-tol','1e-8',
            '--optimal-port-source-mode','generalized-dynamic',
            '--optimal-port-ablation','constant-geometry-generalized',
            '--optimal-port-refinement-validation',
            '--mor-transient-dt','0.1',
            '--mor-transient-t-end','0.1',
            '--mor-transient-waveform','single_step',
            '--mor-transient-initial-mode','ambient',
            '--output-dir',$caseOutput,'--fast-run')
        & $exe @arguments
        if ($LASTEXITCODE -ne 0) {
            throw "$($case.Label) refinement validation failed."
        }
        if (Test-Path -LiteralPath (
                Join-Path $caseOutput 'local_dynamic_schur_summary.csv')) {
            throw "$($case.Label) advanced beyond the no-transient stop."
        }
        $rows = @(Import-Csv -LiteralPath (
            Join-Path $caseOutput (
                'milestone8_woodbury_refinement_validation.csv')))
        if ($rows.Count -ne 2 -or
            @($rows | Where-Object { $_.status -ne 'pass' }).Count -ne 0) {
            throw "$($case.Label) refinement validation gate failed."
        }
        $allRows += $rows
    }
    if ($allRows.Count -ne 4) {
        throw 'Expected four two-/ten-cube refinement rows.'
    }
    $allRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        -LiteralPath (Join-Path $outputs (
            'milestone8_woodbury_refinement_validation.csv'))
    [pscustomobject]@{
        build_type = 'Release'
        cpu_threads = $CpuThreads
        omp_num_threads = $env:OMP_NUM_THREADS
        mkl_num_threads = $env:MKL_NUM_THREADS
        mkl_dynamic = $env:MKL_DYNAMIC
        transient_advanced = 0
        full_field_read = 0
        snapshot_used = 0
        fom_used_for_basis = 0
    } | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (
        Join-Path $result 'environment.csv')
} finally {
    Pop-Location
}

Write-Host 'Two-/ten-cube Woodbury refinement validation passed.'
