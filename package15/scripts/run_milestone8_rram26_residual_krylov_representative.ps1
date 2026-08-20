param(
    [string]$RramConfig =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$TopologyAuditCsv =
        'outputs\milestone8_large_case_topology_audit.csv',
    [string]$WorkDirectory =
        'results\milestone8_rram26_residual_krylov_representative',
    [string]$OutputsDirectory = 'outputs',
    [int]$CpuThreads = 8,
    [switch]$FinalizeExistingRankZero,
    [switch]$SkipBuild
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$work = [IO.Path]::GetFullPath((Join-Path $project $WorkDirectory))
$outputs = [IO.Path]::GetFullPath((Join-Path $project $OutputsDirectory))
$audit = [IO.Path]::GetFullPath((Join-Path $project $TopologyAuditCsv))
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
if ((Test-Path -LiteralPath $work) -and
    -not $FinalizeExistingRankZero) {
    throw "Refusing to overwrite existing work directory: $work"
}
if (-not (Test-Path -LiteralPath $RramConfig -PathType Leaf)) {
    throw "RRAM26 config not found: $RramConfig"
}
if (-not (Test-Path -LiteralPath $audit -PathType Leaf)) {
    throw "Topology audit not found: $audit"
}

Push-Location $project
try {
    if ($FinalizeExistingRankZero) {
        $existing = Join-Path $work `
            'rank_0\milestone8_rram26_residual_krylov_representative.csv'
        if (-not (Test-Path -LiteralPath $existing -PathType Leaf)) {
            throw "Existing rank-zero pilot is unavailable: $existing"
        }
        New-Item -ItemType Directory -Force -Path $outputs | Out-Null
        $all = @(Import-Csv -LiteralPath $existing)
        foreach ($row in $all) {
            if ([double]$row.target_residual -gt 1e-9) {
                $row.status = 'target_residual_gate_failed'
            } elseif ([double]$row.weighted_adjoint_error -gt 1e-8) {
                $row.status = 'weighted_adjoint_gate_failed'
            } elseif ([double]$row.setup_time_s -gt 120.0) {
                $row.status = 'single_interface_time_gate_failed'
            } elseif ([uint64]$row.peak_incremental_memory_bytes -gt 1GB) {
                $row.status = 'incremental_memory_gate_failed'
            } else {
                $row.status = 'success'
            }
        }
        $all | Export-Csv -NoTypeInformation -Encoding UTF8 `
            (Join-Path $outputs `
                'milestone8_rram26_residual_krylov_representative.csv')
        @($all | Where-Object { $_.status -ne 'success' }) |
            Select-Object case,interface_id,
                requested_enrichment_rank,status |
            Export-Csv -NoTypeInformation -Encoding UTF8 `
                (Join-Path $outputs `
                    'milestone8_residual_krylov_failures.csv')
        $timingPath = Join-Path $outputs `
            'milestone8_residual_krylov_timing.csv'
        $timing = if (Test-Path -LiteralPath $timingPath) {
            @(Import-Csv -LiteralPath $timingPath |
                Where-Object { $_.case -ne
                    'rram26_from_sim_parameter_1' })
        } else { @() }
        $timing += @($all | ForEach-Object {
            [pscustomobject]@{
                case = $_.case
                method = 'mandatory-only'
                requested_enrichment_rank =
                    $_.requested_enrichment_rank
                target_solve_count = $_.target_solve_count
                basis_build_time_s = $_.setup_time_s
                total_time_s = $_.setup_time_s
            }
        })
        $timing | Export-Csv -NoTypeInformation -Encoding UTF8 `
            $timingPath
        $memoryPath = Join-Path $outputs `
            'milestone8_residual_krylov_memory.csv'
        $memory = if (Test-Path -LiteralPath $memoryPath) {
            @(Import-Csv -LiteralPath $memoryPath |
                Where-Object { $_.case -ne
                    'rram26_from_sim_parameter_1' })
        } else { @() }
        $memory += @($all | ForEach-Object {
            [pscustomobject]@{
                case = $_.case
                method = 'mandatory-only'
                requested_enrichment_rank =
                    $_.requested_enrichment_rank
                peak_incremental_memory_bytes =
                    $_.peak_incremental_memory_bytes
            }
        })
        $memory | Export-Csv -NoTypeInformation -Encoding UTF8 `
            $memoryPath
        [pscustomobject]@{
            case = 'rram26_from_sim_parameter_1'
            selection = 'all'
            interface_id = -1
            requested_enrichment_rank = 4
            accepted_enrichment_rank = 0
            status = 'not_run_stage_b_gate_failed'
            transient_advanced = 0
            full_field_read = 0
            snapshot_basis_used = 0
        } | Export-Csv -NoTypeInformation -Encoding UTF8 `
            (Join-Path $outputs `
                'milestone8_rram26_residual_krylov_all_interfaces.csv')
        $comparisonPath = Join-Path $outputs `
            'milestone8_residual_krylov_vs_optimal_transfer.csv'
        $comparison = if (Test-Path -LiteralPath $comparisonPath) {
            @(Import-Csv -LiteralPath $comparisonPath |
                Where-Object { $_.case -ne
                    'rram26_from_sim_parameter_1' })
        } else { @() }
        $comparison += [pscustomobject]@{
            case = 'rram26_from_sim_parameter_1'
            method = 'optimal-transfer-interface-24'
            requested_enrichment_rank = 4
            accepted_enrichment_rank = 2
            target_solve_count = 2670
            basis_build_time_s = 600.35
            port_rank = ''
            interface_residual = ''
            flux_relative_l2 = ''
            full_residual = ''
            peak_incremental_memory_bytes = 608583680
        }
        $comparison += @($all | ForEach-Object {
            [pscustomobject]@{
                case = $_.case
                method =
                    "mandatory-only-interface-$($_.interface_id)"
                requested_enrichment_rank = 0
                accepted_enrichment_rank = 0
                target_solve_count = $_.target_solve_count
                basis_build_time_s = $_.setup_time_s
                port_rank = $_.mandatory_rank_total
                interface_residual =
                    $_.initial_max_probe_residual
                flux_relative_l2 = ''
                full_residual = ''
                peak_incremental_memory_bytes =
                    $_.peak_incremental_memory_bytes
            }
        })
        $comparison | Export-Csv -NoTypeInformation -Encoding UTF8 `
            $comparisonPath
        return
    }
    if (-not $SkipBuild) {
        cmake -S . -B build
        if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed.' }
        cmake --build build --config Release -j $CpuThreads
        if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
    }
    New-Item -ItemType Directory -Force -Path $work | Out-Null
    New-Item -ItemType Directory -Force -Path $outputs | Out-Null
    $env:OMP_NUM_THREADS = "$CpuThreads"
    $env:MKL_NUM_THREADS = "$CpuThreads"
    $env:MKL_DYNAMIC = 'FALSE'
    $all = @()
    foreach ($rank in @(0,1,2,4)) {
        $method = if ($rank -eq 0) {
            'mandatory-only'
        } else {
            'residual-krylov'
        }
        $directory = Join-Path $work "rank_$rank"
        $arguments = @(
            '--transient','--config',$RramConfig,
            '--mor-transient-generate',
            '--mor-transient-method','local-port-block-arnoldi',
            '--port-basis-method',$method,
            '--mor-arnoldi-moments','1',
            '--residual-krylov-max-rank',"$rank",
            '--residual-krylov-max-sweeps','2',
            '--residual-krylov-tol','1e-4',
            '--residual-krylov-block-size','2',
            '--residual-krylov-probe-mode','operator-geometry',
            '--residual-krylov-inner-solver','woodbury-exact',
            '--residual-krylov-representative-pilot',
            '--optimal-port-topology-audit-csv',$audit,
            '--optimal-port-inner-solver','woodbury-exact',
            '--optimal-port-inner-tol','1e-10',
            '--optimal-port-inner-refinement-max-iters','3',
            '--optimal-port-inner-refinement-tol','1e-10',
            '--mor-transient-dt','0.1',
            '--mor-transient-t-end','0.1',
            '--mor-transient-waveform','single_step',
            '--mor-transient-initial-mode','ambient',
            '--output-dir',$directory,'--fast-run')
        & $exe @arguments
        if ($LASTEXITCODE -ne 0) {
            throw "RRAM26 representative rank $rank failed."
        }
        $rows = @(Import-Csv (
            Join-Path $directory `
                'milestone8_rram26_residual_krylov_representative.csv'))
        if ($rows.Count -ne 3) {
            throw "Expected three representative rows for rank $rank."
        }
        $all += $rows
        $gateFailed = @($rows | Where-Object {
            $_.status -ne 'success' -or
            [double]$_.target_residual -gt 1e-9 -or
            [double]$_.weighted_adjoint_error -gt 1e-8 -or
            [double]$_.setup_time_s -gt 120.0 -or
            [uint64]$_.peak_incremental_memory_bytes -gt 1GB
        })
        if ($gateFailed.Count) {
            break
        }
    }
    $destination = Join-Path $outputs `
        'milestone8_rram26_residual_krylov_representative.csv'
    $all | Export-Csv -NoTypeInformation -Encoding UTF8 $destination
    $failed = @($all | Where-Object { $_.status -ne 'success' })
    if ($all.Count -eq 12 -and $failed.Count -eq 0) {
        [pscustomobject]@{
            status = 'stage_b_passed'
            configurations = 12
            transient_advanced = 0
            full_field_read = 0
            snapshot_basis_used = 0
        } | Export-Csv -NoTypeInformation -Encoding UTF8 `
            (Join-Path $outputs `
                'milestone8_rram26_residual_krylov_stage_b_gate.csv')
    } else {
        $failed | Select-Object case,interface_id,
            requested_enrichment_rank,status |
            Export-Csv -NoTypeInformation -Encoding UTF8 `
                (Join-Path $outputs `
                    'milestone8_residual_krylov_failures.csv')
    }
} finally {
    Pop-Location
}
