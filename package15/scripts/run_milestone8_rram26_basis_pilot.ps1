param(
    [string]$RramConfig =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$ResultsDirectory =
        'results\milestone8_rram26_basis_pilot',
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

if (-not (Test-Path -LiteralPath $RramConfig -PathType Leaf)) {
    throw "RRAM26 config not found: $RramConfig"
}
if ($CpuThreads -le 0) {
    throw 'CpuThreads must be positive.'
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite an existing pilot: $result"
}

Push-Location $project
try {
    if (-not $SkipBuild) {
        cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
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

    $log = Join-Path $result 'console.log'
    $console = @(& $exe --transient --config $RramConfig `
        --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --port-basis-method optimal-transfer `
        --mor-arnoldi-moments 1 `
        --optimal-port-source-mode trace-only `
        --optimal-port-ablation constant-geometry-trace `
        --optimal-port-basis-pilot `
        --optimal-port-inner-solver auto `
        --optimal-port-inner-tol 1e-10 `
        --optimal-port-inner-max-iters 1000 `
        --optimal-port-eigen-tol 1e-8 `
        --optimal-port-eigen-max-iters 1000 `
        --mor-transient-dt 0.1 `
        --mor-transient-t-end 0.1 `
        --mor-transient-waveform single_step `
        --mor-transient-initial-mode ambient `
        --output-dir $result --fast-run 2>&1)
    $exitCode = $LASTEXITCODE
    $console | Set-Content -Encoding UTF8 -LiteralPath $log
    if ($exitCode -ne 0) {
        throw "RRAM26 basis pilot failed with exit code $exitCode. See $log"
    }

    $pilot = Join-Path $result 'milestone8_rram26_basis_pilot.csv'
    $stop = Join-Path $result 'milestone8_rram26_basis_pilot_stop.csv'
    if (-not (Test-Path -LiteralPath $pilot -PathType Leaf) -or
        -not (Test-Path -LiteralPath $stop -PathType Leaf)) {
        throw 'Pilot stop artifacts are missing.'
    }
    $rows = @(Import-Csv -LiteralPath $pilot)
    if ($rows.Count -ne 9) {
        throw "Expected 9 pilot rows, found $($rows.Count)."
    }
    if (@($rows | Where-Object {
            [int]$_.total_rank_budget -eq 24 }).Count -ne 0) {
        throw 'Forbidden rank 24 was executed.'
    }
    foreach ($forbidden in @(
            'local_dynamic_schur_summary.csv',
            'local_dynamic_schur_accuracy_by_time.csv',
            'temperature_local_port_rom_nodes.csv')) {
        if (Test-Path -LiteralPath (Join-Path $result $forbidden)) {
            throw "Pilot advanced beyond its stop gate: $forbidden"
        }
    }

    Copy-Item -LiteralPath $pilot -Destination (
        Join-Path $outputs 'milestone8_rram26_basis_pilot.csv')
    [pscustomobject]@{
        build_type = 'Release'
        cpu_threads = $CpuThreads
        omp_num_threads = $env:OMP_NUM_THREADS
        mkl_num_threads = $env:MKL_NUM_THREADS
        mkl_dynamic = $env:MKL_DYNAMIC
        eigensolver_tolerance = 1e-8
        eigensolver_maximum_iterations = 1000
        inner_solver_tolerance = 1e-10
        preconditioner_diagonal = 'assembled-interface'
        rank_budgets = '8;12;16'
        transient_advanced = 0
        full_field_read = 0
        snapshot_basis_used = 0
    } | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (
        Join-Path $result 'milestone8_rram26_basis_pilot_environment.csv')
} finally {
    Pop-Location
}

Write-Host 'RRAM26 three-interface basis pilot finished at the required stop gate.'
