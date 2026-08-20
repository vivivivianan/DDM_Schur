param(
    [string]$RramConfig =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$CaseName = 'rram26_from_sim_parameter_1',
    [string]$TopologyAuditCsv =
        'outputs\milestone8_large_case_topology_audit.csv',
    [string]$ResultsDirectory =
        'results\milestone8_rram26_representative_interfaces',
    [string]$OutputsDirectory = 'outputs',
    [int]$CpuThreads = 8,
    [int]$ExternalTimeoutSeconds = 4200,
    [switch]$RunAllInterfaces,
    [switch]$SkipBuild
)

$ErrorActionPreference = 'Stop'
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
if ($CpuThreads -le 0 -or $ExternalTimeoutSeconds -lt 3600) {
    throw 'Invalid CPU thread count or Stage-A external timeout.'
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite an existing Stage-A run: $result"
}

$auditRows = @(Import-Csv -LiteralPath $auditPath |
    Where-Object { $_.case -eq $CaseName } |
    Sort-Object @{ Expression = { [int]$_.target_dofs } },
                @{ Expression = { [int]$_.interface_id } })
if ($auditRows.Count -ne 25) {
    throw "Expected 25 audited RRAM26 interfaces; found $($auditRows.Count)."
}
$selected = @(
    [pscustomobject]@{
        selection = 'minimum'
        interface_id = [int]$auditRows[0].interface_id
        target_dofs = [int]$auditRows[0].target_dofs
        source_dofs = [int]$auditRows[0].source_dofs
    },
    [pscustomobject]@{
        selection = 'median'
        interface_id = [int]$auditRows[[math]::Floor(
            $auditRows.Count / 2)].interface_id
        target_dofs = [int]$auditRows[[math]::Floor(
            $auditRows.Count / 2)].target_dofs
        source_dofs = [int]$auditRows[[math]::Floor(
            $auditRows.Count / 2)].source_dofs
    },
    [pscustomobject]@{
        selection = 'maximum'
        interface_id = [int]$auditRows[-1].interface_id
        target_dofs = [int]$auditRows[-1].target_dofs
        source_dofs = [int]$auditRows[-1].source_dofs
    })

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
    $selected | Export-Csv -NoTypeInformation -Encoding UTF8 `
        -LiteralPath (Join-Path $result 'selected_interfaces.csv')
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
        '--optimal-port-rank','12',
        '--optimal-port-source-mode','trace-only',
        '--optimal-port-ablation','constant-geometry-trace',
        '--optimal-port-representative-interface-pilot',
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
        throw "Stage A exceeded external timeout: $ExternalTimeoutSeconds s"
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
        throw "Stage A failed with exit code $exitCode."
    }

    $pilotPath = Join-Path $result (
        'milestone8_rram26_representative_interface_pilot.csv')
    $timingPath = Join-Path $result (
        'milestone8_rram26_basis_timing.csv')
    $summaryPath = Join-Path $result (
        'milestone8_rram26_all_interface_summary.csv')
    Add-Content -LiteralPath $timingPath -Encoding UTF8 -Value (
        "representative,all,-1,0,complete_process_wall,$($wall.Elapsed.TotalSeconds.ToString('R',[Globalization.CultureInfo]::InvariantCulture))")
    $pilotRows = @(Import-Csv -LiteralPath $pilotPath)
    if ($pilotRows.Count -lt 1 -or $pilotRows.Count -gt 6) {
        throw "Stage A returned an invalid row count: $($pilotRows.Count)."
    }
    $expectedIds = @($selected | ForEach-Object { $_.interface_id })
    $actualIds = @($pilotRows | Select-Object -ExpandProperty interface_id -Unique |
        ForEach-Object { [int]$_ })
    $unexpectedIds = @($actualIds | Where-Object { $_ -notin $expectedIds })
    if ($unexpectedIds.Count -ne 0) {
        throw 'C++ representative selection is not a topology-audit subset.'
    }
    $actualRanks = @($pilotRows |
        Select-Object -ExpandProperty requested_transfer_rank -Unique |
        ForEach-Object { [int]$_ } | Sort-Object)
    $unexpectedRanks = @($actualRanks | Where-Object { $_ -notin @(4, 8) })
    if ($unexpectedRanks.Count -ne 0) {
        throw 'Stage A used a transfer rank other than 4 or 8.'
    }
    $failed = @($pilotRows | Where-Object { $_.status -ne 'passed' })
    $allSummary = @(Import-Csv -LiteralPath $summaryPath)
    if ($allSummary.Count -ne 1 -or
        $allSummary[0].status -ne 'not_run_stage_a_stop') {
        throw 'Stage A did not preserve the all-interface stop gate.'
    }

    foreach ($forbidden in @(
            'local_dynamic_schur_summary.csv',
            'local_dynamic_schur_accuracy_by_time.csv',
            'temperature_local_port_rom_nodes.csv')) {
        if (Test-Path -LiteralPath (Join-Path $result $forbidden)) {
            throw "Stage A advanced beyond its stop gate: $forbidden"
        }
    }
    [pscustomobject]@{
        build_type = 'Release'
        cpu_threads = $CpuThreads
        omp_num_threads = $env:OMP_NUM_THREADS
        mkl_num_threads = $env:MKL_NUM_THREADS
        mkl_dynamic = $env:MKL_DYNAMIC
        representative_selection_source = $auditPath
        representative_interfaces = ($selected.interface_id -join ';')
        transfer_ranks = '4;8'
        transient_advanced = 0
        full_field_read = 0
        snapshot_used = 0
        fom_used_for_basis = 0
        complete_process_wall_time_s = $wall.Elapsed.TotalSeconds
        all_interfaces_requested = [int]$RunAllInterfaces.IsPresent
    } | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (
        Join-Path $result 'environment.csv')

    $artifacts = @(
        'milestone8_rram26_representative_interface_pilot.csv',
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
    if ($failed.Count -ne 0) {
        throw "Stage A stopped at failed interface $($failed[0].interface_id), rank $($failed[0].requested_transfer_rank): $($failed[0].failure_reason)"
    }
    if ($pilotRows.Count -ne 6) {
        throw "Stage A stopped without a failed row after $($pilotRows.Count) rows."
    }
} finally {
    Pop-Location
}

Write-Host 'RRAM26 representative-interface Stage A passed and stopped.'

if ($RunAllInterfaces) {
    & (Join-Path $PSScriptRoot 'run_milestone8_rram26_all_interface_basis.ps1') `
        -RramConfig $RramConfig -CaseName $CaseName `
        -TopologyAuditCsv $TopologyAuditCsv `
        -OutputsDirectory $OutputsDirectory -CpuThreads $CpuThreads `
        -RunAllInterfaces -SkipBuild
    if ($LASTEXITCODE -ne 0) {
        throw 'Explicit Stage B invocation failed.'
    }
}
