param(
    [string]$RramConfig =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$CaseName = 'rram26_from_sim_parameter_1',
    [string]$TopologyAuditCsv =
        'outputs\milestone8_large_case_topology_audit.csv',
    [string]$ResultsDirectory =
        'results\milestone8_rram26_all_interface_rank4',
    [string]$OutputsDirectory = 'outputs',
    [int]$CpuThreads = 8,
    [int]$ExternalTimeoutSeconds = 16200,
    [switch]$RunAllInterfaces,
    [switch]$SkipBuild
)

$ErrorActionPreference = 'Stop'
if (-not $RunAllInterfaces) {
    throw 'Stage B is disabled by default. Pass -RunAllInterfaces only after explicit confirmation.'
}
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$result = [System.IO.Path]::GetFullPath(
    (Join-Path $project $ResultsDirectory))
$outputs = [System.IO.Path]::GetFullPath(
    (Join-Path $project $OutputsDirectory))
$auditPath = [System.IO.Path]::GetFullPath(
    (Join-Path $project $TopologyAuditCsv))
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'

if (-not (Test-Path -LiteralPath $RramConfig -PathType Leaf)) {
    throw "RRAM26 config not found: $RramConfig"
}
if (-not (Test-Path -LiteralPath $auditPath -PathType Leaf)) {
    throw "Topology audit CSV not found: $auditPath"
}
if ($CpuThreads -le 0 -or $ExternalTimeoutSeconds -lt 15000) {
    throw 'Invalid CPU thread count or Stage-B external timeout.'
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite an existing Stage-B run: $result"
}
$auditRows = @(Import-Csv -LiteralPath $auditPath |
    Where-Object { $_.case -eq $CaseName })
if ($auditRows.Count -ne 25) {
    throw "Expected 25 audited RRAM26 interfaces; found $($auditRows.Count)."
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
        '--optimal-port-rank','8',
        '--optimal-port-source-mode','trace-only',
        '--optimal-port-ablation','constant-geometry-trace',
        '--optimal-port-all-interface-basis',
        '--optimal-port-topology-audit-csv',$auditPath,
        '--optimal-port-inner-solver','woodbury-exact',
        '--optimal-port-inner-tol','1e-10',
        '--optimal-port-inner-max-iters','1000',
        '--optimal-port-eigen-tol','1e-8',
        '--optimal-port-eigen-max-iters','500',
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
    $wall = [System.Diagnostics.Stopwatch]::StartNew()
    $process = Start-Process -FilePath $exe -ArgumentList $argumentLine `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr `
        -WindowStyle Hidden -PassThru
    if (-not $process.WaitForExit($ExternalTimeoutSeconds * 1000)) {
        Stop-Process -Id $process.Id -Force
        throw "Stage B exceeded external timeout: $ExternalTimeoutSeconds s"
    }
    $process.WaitForExit()
    $process.Refresh()
    $wall.Stop()
    $exitCode = $process.ExitCode
    if ($null -eq $exitCode) {
        $exitCode = if ((Get-Item -LiteralPath $stderr).Length -gt 0) {
            1
        } else {
            0
        }
    }
    if ($exitCode -ne 0) {
        throw "Stage B failed with exit code $exitCode."
    }

    $basisPath = Join-Path $result (
        'milestone8_rram26_all_interface_basis.csv')
    $summaryPath = Join-Path $result (
        'milestone8_rram26_all_interface_summary.csv')
    $timingPath = Join-Path $result (
        'milestone8_rram26_basis_timing.csv')
    $basisRows = @(Import-Csv -LiteralPath $basisPath)
    $summaryRows = @(Import-Csv -LiteralPath $summaryPath)
    if ($basisRows.Count -ne 25 -or
        @($basisRows | Where-Object { $_.status -ne 'passed' }).Count -ne 0) {
        throw 'Stage B did not complete all 25 interfaces successfully.'
    }
    if ($summaryRows.Count -ne 1 -or
        $summaryRows[0].status -ne 'passed' -or
        [int]$summaryRows[0].successful_interfaces -ne 25) {
        throw 'Stage-B summary did not pass.'
    }
    Add-Content -LiteralPath $timingPath -Encoding UTF8 -Value (
        "all-interface,all,-1,4,complete_process_wall,$($wall.Elapsed.TotalSeconds.ToString('R',[Globalization.CultureInfo]::InvariantCulture))")
    foreach ($forbidden in @(
            'local_dynamic_schur_summary.csv',
            'local_dynamic_schur_accuracy_by_time.csv',
            'temperature_local_port_rom_nodes.csv')) {
        if (Test-Path -LiteralPath (Join-Path $result $forbidden)) {
            throw "Stage B advanced beyond its stop gate: $forbidden"
        }
    }
    [pscustomobject]@{
        build_type = 'Release'
        cpu_threads = $CpuThreads
        representative_stage_required = 1
        run_all_interfaces_explicit = 1
        physical_interfaces = 25
        requested_transfer_rank = 4
        transient_advanced = 0
        full_field_read = 0
        snapshot_used = 0
        fom_used_for_basis = 0
        complete_process_wall_time_s = $wall.Elapsed.TotalSeconds
    } | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (
        Join-Path $result 'environment.csv')

    $artifacts = @(
        'milestone8_rram26_all_interface_basis.csv',
        'milestone8_rram26_all_interface_summary.csv',
        'milestone8_rram26_interface_eigenvalues.csv',
        'milestone8_rram26_basis_timing.csv',
        'milestone8_rram26_basis_memory.csv',
        'milestone8_rram26_basis_failures.csv')
    foreach ($artifact in $artifacts) {
        Copy-Item -LiteralPath (Join-Path $result $artifact) `
            -Destination (Join-Path $outputs $artifact)
    }
} finally {
    Pop-Location
}

Write-Host 'Explicit RRAM26 all-interface rank-4 basis build passed; transient skipped.'
