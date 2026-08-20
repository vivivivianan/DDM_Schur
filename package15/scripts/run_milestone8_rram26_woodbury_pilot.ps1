param(
    [string]$RramConfig =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$ResultsDirectory =
        'results\milestone8_rram26_woodbury_pilot',
    [string]$OutputsDirectory = 'outputs',
    [int]$CpuThreads = 8,
    [int]$ExternalTimeoutSeconds = 630,
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
if ($CpuThreads -le 0 -or $ExternalTimeoutSeconds -lt 600) {
    throw 'Invalid CPU thread count or external timeout.'
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite an existing pilot: $result"
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

    $stdout = Join-Path $result 'console.log'
    $stderr = Join-Path $result 'console.err.log'
    $arguments = @(
        '--transient','--config',$RramConfig,
        '--mor-transient-generate',
        '--mor-transient-method','local-port-block-arnoldi',
        '--port-basis-method','optimal-transfer',
        '--mor-arnoldi-moments','1',
        '--optimal-port-source-mode','trace-only',
        '--optimal-port-ablation','constant-geometry-trace',
        '--optimal-port-woodbury-pilot',
        '--optimal-port-inner-solver','woodbury-exact',
        '--optimal-port-inner-tol','1e-10',
        '--optimal-port-inner-max-iters','1000',
        '--optimal-port-eigen-tol','1e-8',
        '--optimal-port-eigen-max-iters','1000',
        '--mor-transient-dt','0.1',
        '--mor-transient-t-end','0.1',
        '--mor-transient-waveform','single_step',
        '--mor-transient-initial-mode','ambient',
        '--output-dir',$result,'--fast-run')
    $argumentLine = (($arguments | ForEach-Object {
        $value = [string]$_
        if ($value -match '\s|"') {
            '"' + $value.Replace('"', '\"') + '"'
        } else {
            $value
        }
    }) -join ' ')
    $process = Start-Process -FilePath $exe -ArgumentList $argumentLine `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr `
        -WindowStyle Hidden -PassThru
    if (-not $process.WaitForExit($ExternalTimeoutSeconds * 1000)) {
        Stop-Process -Id $process.Id -Force
        throw "RRAM26 pilot exceeded external timeout: $ExternalTimeoutSeconds s"
    }
    # Complete redirected stream draining before querying ExitCode. Without
    # this second wait Windows PowerShell can expose a null ExitCode even
    # though the child and both redirected streams have finished.
    $process.WaitForExit()
    $process.Refresh()
    $exitCode = $process.ExitCode
    if ($exitCode -ne 0) {
        throw "RRAM26 pilot failed with exit code $exitCode."
    }

    $pilot = Join-Path $result 'milestone8_rram26_woodbury_pilot.csv'
    $memory = Join-Path $result 'milestone8_target_solver_memory.csv'
    if (-not (Test-Path -LiteralPath $pilot -PathType Leaf) -or
        -not (Test-Path -LiteralPath $memory -PathType Leaf)) {
        throw 'Woodbury pilot artifacts are missing.'
    }
    $row = Import-Csv -LiteralPath $pilot
    if ([int]$row.interface_id -ne 16 -or
        [int]$row.target_dofs -ne 3299 -or
        [int]$row.source_dofs -ne 4169) {
        throw 'RRAM26 interface-16 topology changed.'
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
        Join-Path $outputs 'milestone8_rram26_woodbury_pilot.csv')
    Copy-Item -LiteralPath $memory -Destination (
        Join-Path $outputs 'milestone8_target_solver_memory.csv')
    [pscustomobject]@{
        build_type = 'Release'
        cpu_threads = $CpuThreads
        omp_num_threads = $env:OMP_NUM_THREADS
        mkl_num_threads = $env:MKL_NUM_THREADS
        mkl_dynamic = $env:MKL_DYNAMIC
        interface_id = 16
        transfer_rank = 4
        transient_advanced = 0
        full_field_read = 0
        snapshot_basis_used = 0
    } | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (
        Join-Path $result 'environment.csv')
} finally {
    Pop-Location
}

Write-Host 'RRAM26 interface-16 Woodbury pilot reached the required stop gate.'
