param(
    [string]$ResultsDirectory =
        'results\milestone8_target_solver_comparison',
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
    throw 'CpuThreads must be positive.'
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite an existing comparison: $result"
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

    $cases = @(
        [pscustomobject]@{
            Name = 'two_cube'
            Config = 'configs\two_cube_parametric_h.txt'
            Rank = 7
            SourceMode = 'generalized-dynamic'
            Ablation = 'constant-geometry-generalized'
        },
        [pscustomobject]@{
            Name = 'ten_cube'
            Config = 'configs\ten_cube_parametric_h.txt'
            Rank = 8
            SourceMode = 'trace-only'
            Ablation = 'constant-geometry-trace'
        })

    $allRows = @()
    foreach ($case in $cases) {
        $caseResult = Join-Path $result $case.Name
        New-Item -ItemType Directory -Path $caseResult | Out-Null
        $log = Join-Path $caseResult 'console.log'
        $console = @(& $exe --transient --config $case.Config `
            --mor-transient-generate `
            --mor-transient-method local-port-block-arnoldi `
            --port-basis-method optimal-transfer `
            --mor-arnoldi-moments 2 `
            --optimal-port-rank $case.Rank `
            --optimal-port-source-mode $case.SourceMode `
            --optimal-port-ablation $case.Ablation `
            --optimal-port-target-solver-comparison `
            --optimal-port-inner-tol 1e-12 `
            --optimal-port-inner-max-iters 1000 `
            --optimal-port-eigen-tol 1e-8 `
            --optimal-port-eigen-max-iters 1000 `
            --mor-transient-dt 0.1 `
            --mor-transient-t-end 0.1 `
            --mor-transient-waveform single_step `
            --mor-transient-initial-mode ambient `
            --output-dir $caseResult --fast-run 2>&1)
        $exitCode = $LASTEXITCODE
        $console | Set-Content -Encoding UTF8 -LiteralPath $log
        if ($exitCode -ne 0) {
            throw "$($case.Name) solver comparison failed. See $log"
        }
        $csv = Join-Path $caseResult (
            'milestone8_target_solver_comparison.csv')
        $rows = @(Import-Csv -LiteralPath $csv)
        if ($rows.Count -ne 3 -or
            @($rows | Where-Object { $_.status -ne 'pass' }).Count -ne 0) {
            throw "$($case.Name) target-solver gates did not all pass."
        }
        $allRows += $rows
        foreach ($forbidden in @(
                'local_dynamic_schur_summary.csv',
                'local_dynamic_schur_accuracy_by_time.csv',
                'temperature_local_port_rom_nodes.csv')) {
            if (Test-Path -LiteralPath (Join-Path $caseResult $forbidden)) {
                throw "$($case.Name) advanced beyond the comparison stop gate."
            }
        }
    }

    $allRows | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (
        Join-Path $outputs 'milestone8_target_solver_comparison.csv')
    [pscustomobject]@{
        build_type = 'Release'
        cpu_threads = $CpuThreads
        omp_num_threads = $env:OMP_NUM_THREADS
        mkl_num_threads = $env:MKL_NUM_THREADS
        mkl_dynamic = $env:MKL_DYNAMIC
        inner_solver_tolerance = 1e-12
        eigensolver_tolerance = 1e-8
        eigensolver_maximum_iterations = 1000
        transient_advanced = 0
        full_field_read = 0
        snapshot_basis_used = 0
    } | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (
        Join-Path $result 'environment.csv')
} finally {
    Pop-Location
}

Write-Host 'Two/ten-cube exact target-solver comparison passed.'
