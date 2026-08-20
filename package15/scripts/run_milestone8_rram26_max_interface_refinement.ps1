param(
    [string]$RramConfig =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$CaseName = 'rram26_from_sim_parameter_1',
    [string]$TopologyAuditCsv =
        'outputs\milestone8_large_case_topology_audit.csv',
    [string]$ResultsDirectory =
        'results\milestone8_rram26_max_interface_refinement',
    [string]$OutputsDirectory = 'outputs',
    [int]$CpuThreads = 8,
    [int]$ExternalTimeoutSeconds = 1500,
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
if ($CpuThreads -le 0 -or $ExternalTimeoutSeconds -lt 1200) {
    throw 'Invalid CPU thread count or maximum-interface timeout.'
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite maximum-interface run: $result"
}

$auditRows = @(Import-Csv -LiteralPath $auditPath |
    Where-Object { $_.case -eq $CaseName } |
    Sort-Object @{ Expression = { [int]$_.target_dofs } },
                @{ Expression = { [int]$_.interface_id } })
if ($auditRows.Count -ne 25) {
    throw "Expected 25 audited RRAM26 interfaces; found $($auditRows.Count)."
}
$maximum = $auditRows[-1]
if ([int]$maximum.interface_id -ne 24 -or
    [int]$maximum.target_dofs -ne 38921 -or
    [int]$maximum.source_dofs -ne 18833) {
    throw 'Audited RRAM26 maximum interface is no longer interface 24.'
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
        '--optimal-port-rank','12',
        '--optimal-port-source-mode','trace-only',
        '--optimal-port-ablation','constant-geometry-trace',
        '--optimal-port-maximum-interface-refinement-pilot',
        '--optimal-port-topology-audit-csv',$auditPath,
        '--optimal-port-inner-solver','woodbury-exact',
        '--optimal-port-inner-tol','1e-10',
        '--optimal-port-inner-max-iters','1000',
        '--optimal-port-inner-refinement-max-iters','3',
        '--optimal-port-inner-refinement-tol','1e-10',
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
        throw "Maximum-interface run exceeded $ExternalTimeoutSeconds s."
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
        throw "Maximum-interface process failed: $exitCode."
    }

    $pilotPath = Join-Path $result (
        'milestone8_rram26_representative_interface_pilot.csv')
    $breakdownPath = Join-Path $result (
        'milestone8_woodbury_residual_breakdown.csv')
    $failurePath = Join-Path $result (
        'milestone8_rram26_basis_failures.csv')
    $pilotRows = @(Import-Csv -LiteralPath $pilotPath)
    $breakdownRows = @(Import-Csv -LiteralPath $breakdownPath)
    if ($pilotRows.Count -lt 1 -or $pilotRows.Count -gt 2 -or
        @($pilotRows | Where-Object {
            [int]$_.interface_id -ne [int]$maximum.interface_id
        }).Count -ne 0) {
        throw 'Maximum-interface mode returned an invalid interface set.'
    }
    $rank4 = $pilotRows | Where-Object {
        [int]$_.requested_transfer_rank -eq 4 }
    if ($null -eq $rank4) {
        throw 'Maximum-interface rank 4 row is missing.'
    }
    $rank4Passed = $rank4.status -eq 'passed'
    $rank8 = $pilotRows | Where-Object {
        [int]$_.requested_transfer_rank -eq 8 }
    if (-not $rank4Passed -and $null -ne $rank8) {
        throw 'Rank 8 ran after the required rank-4 stop gate failed.'
    }

    $joined = foreach ($pilot in $pilotRows) {
        $residual = $breakdownRows | Where-Object {
            [int]$_.interface_id -eq [int]$pilot.interface_id -and
            [int]$_.requested_transfer_rank -eq
                [int]$pilot.requested_transfer_rank
        } | Select-Object -First 1
        if ($null -eq $residual) {
            throw "Residual breakdown missing for rank $($pilot.requested_transfer_rank)."
        }
        [pscustomobject]@{
            selection = $pilot.selection
            interface_id = $pilot.interface_id
            adjacent_subdomains = $pilot.adjacent_subdomains
            target_dofs = $pilot.target_dofs
            source_dofs = $pilot.source_dofs
            mandatory_rank = $pilot.mandatory_rank
            requested_transfer_rank = $pilot.requested_transfer_rank
            converged_transfer_rank = $pilot.converged_transfer_rank
            total_port_rank = $pilot.total_port_rank
            A_tt_nnz = $pilot.A_tt_nnz
            woodbury_setup_time_s = $pilot.woodbury_setup_time_s
            mean_target_solve_time_s = $pilot.mean_target_solve_time_s
            A_tt_solve_relative_residual =
                $residual.A_tt_solve_relative_residual
            Q_solve_pre_refinement_residual =
                $residual.Q_solve_pre_refinement_residual
            Q_solve_relative_residual =
                $residual.Q_solve_relative_residual
            Q_refinement_iterations =
                $residual.Q_refinement_iterations
            woodbury_pre_refinement_residual =
                $residual.woodbury_pre_refinement_residual
            woodbury_post_refinement_residual =
                $residual.woodbury_post_refinement_residual
            refinement_iterations = $residual.refinement_iterations
            refinement_residual_0 = $residual.refinement_residual_0
            refinement_residual_1 = $residual.refinement_residual_1
            refinement_residual_2 = $residual.refinement_residual_2
            refinement_residual_3 = $residual.refinement_residual_3
            refinement_correction_relative_norm =
                $residual.refinement_correction_relative_norm
            refinement_reduction_factor =
                $residual.refinement_reduction_factor
            refinement_converged = $residual.refinement_converged
            refinement_triggered_solve_calls =
                $residual.refinement_triggered_solve_calls
            pardiso_internal_refinement_steps =
                $residual.pardiso_internal_refinement_steps
            Q_min_abs_factor_diagonal =
                $residual.Q_min_abs_factor_diagonal
            Q_max_abs_factor_diagonal =
                $residual.Q_max_abs_factor_diagonal
            Q_factor_diagonal_ratio =
                $residual.Q_factor_diagonal_ratio
            woodbury_cancellation_factor =
                $residual.woodbury_cancellation_factor
            target_solve_residual = $pilot.target_solve_residual
            weighted_adjoint_error = $pilot.weighted_adjoint_error
            eigensolver_iterations = $pilot.eigensolver_iterations
            eigensolver_residual = $pilot.eigensolver_residual
            transfer_eigenvalues = $pilot.transfer_eigenvalues
            total_interface_basis_time_s =
                $pilot.total_interface_basis_time_s
            incremental_workspace_bytes =
                $pilot.incremental_workspace_bytes
            process_peak_memory_bytes =
                $pilot.process_peak_memory_bytes
            status = $pilot.status
            failure_reason = $pilot.failure_reason
        }
    }
    $joined | Export-Csv -NoTypeInformation -Encoding UTF8 `
        -LiteralPath (Join-Path $outputs (
            'milestone8_rram26_max_interface_refinement.csv'))
    Copy-Item -LiteralPath $breakdownPath -Destination (
        Join-Path $outputs 'milestone8_woodbury_residual_breakdown.csv')

    $representativeOutput = Join-Path $outputs (
        'milestone8_rram26_representative_interface_pilot.csv')
    $previousRepresentative = @()
    if (Test-Path -LiteralPath $representativeOutput) {
        $previousRepresentative = @(Import-Csv -LiteralPath (
            $representativeOutput) | Where-Object {
                [int]$_.interface_id -ne [int]$maximum.interface_id
            })
    }
    @($previousRepresentative + $pilotRows) |
        Sort-Object @{ Expression = {
            switch ($_.selection) {
                'minimum' { 0 }
                'median' { 1 }
                default { 2 }
            }}},
            @{ Expression = { [int]$_.requested_transfer_rank }} |
        Export-Csv -NoTypeInformation -Encoding UTF8 `
            -LiteralPath $representativeOutput

    $failureOutput = Join-Path $outputs (
        'milestone8_rram26_basis_failures.csv')
    $previousFailures = @()
    if (Test-Path -LiteralPath $failureOutput) {
        $previousFailures = @(Import-Csv -LiteralPath $failureOutput |
            Where-Object {
                [int]$_.interface_id -ne [int]$maximum.interface_id
            })
    }
    $newFailures = @(Import-Csv -LiteralPath $failurePath)
    $combinedFailures = @($previousFailures + $newFailures)
    if ($combinedFailures.Count -gt 0) {
        $combinedFailures | Export-Csv -NoTypeInformation -Encoding UTF8 `
            -LiteralPath $failureOutput
    } else {
        Set-Content -LiteralPath $failureOutput -Encoding UTF8 -Value (
            'scope,selection,interface_id,requested_transfer_rank,stage,status,failure_reason')
    }

    foreach ($forbidden in @(
            'local_dynamic_schur_summary.csv',
            'local_dynamic_schur_accuracy_by_time.csv',
            'temperature_local_port_rom_nodes.csv',
            'milestone8_rram26_rank4_port_basis.bin')) {
        if (Test-Path -LiteralPath (Join-Path $result $forbidden)) {
            throw "Maximum-interface run advanced beyond its stop: $forbidden"
        }
    }
    [pscustomobject]@{
        build_type = 'Release'
        cpu_threads = $CpuThreads
        omp_num_threads = $env:OMP_NUM_THREADS
        mkl_num_threads = $env:MKL_NUM_THREADS
        mkl_dynamic = $env:MKL_DYNAMIC
        interface_id = [int]$maximum.interface_id
        transfer_ranks = if ($null -ne $rank8) { '4;8' } else { '4' }
        rank4_passed = [int]$rank4Passed
        all_interface_rank4_qualified = [int]$rank4Passed
        transient_advanced = 0
        full_field_read = 0
        snapshot_used = 0
        fom_used_for_basis = 0
        complete_process_wall_time_s = $wall.Elapsed.TotalSeconds
    } | Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath (
        Join-Path $result 'environment.csv')

    if (-not $rank4Passed) {
        throw "Maximum interface rank 4 failed: $($rank4.failure_reason)"
    }
    if ($null -eq $rank8) {
        throw 'Rank 4 passed but the required separate rank-8 run is absent.'
    }
    if ($rank8.status -ne 'passed') {
        Write-Warning (
            "Maximum-interface rank 8 extension failed: " +
            $rank8.failure_reason)
    }
} finally {
    Pop-Location
}

Write-Host (
    'RRAM26 maximum-interface rank 4 passed; rank 8 extension recorded. ' +
    'No all-interface basis or transient was run.')
