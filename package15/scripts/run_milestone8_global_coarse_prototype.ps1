param(
    [string]$Config =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$LocalPortModel =
        'D:\AI agent\kimi\DDM_Schur_m8_runs\milestone8_adaptive_port_basis_rram26\local_port_basis.bin',
    [string]$M89Accuracy =
        'D:\AI agent\kimi\DDM_Schur_m8_runs\milestone8_adaptive_port_accuracy_final_diagnostics_rram26\milestone8_adaptive_port_accuracy.csv',
    [string]$ResultsDirectory =
        'results\milestone8_global_coarse_prototype',
    [int[]]$Ranks = @(4, 8, 16, 32),
    [int[]]$InterfaceIds = @(13, 23, 18, 4, 10),
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$configPath = (Resolve-Path -LiteralPath $Config).Path
$modelPath = (Resolve-Path -LiteralPath $LocalPortModel).Path
$accuracyPath = (Resolve-Path -LiteralPath $M89Accuracy).Path
$result = if ([IO.Path]::IsPathRooted($ResultsDirectory)) {
    [IO.Path]::GetFullPath($ResultsDirectory)
} else {
    [IO.Path]::GetFullPath(
        (Join-Path $project $ResultsDirectory))
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite M8.10.1 directory: $result"
}
if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
    throw "Release executable not found: $exe"
}
$supportedRanks = @(4, 8, 16, 32)
foreach ($rank in $Ranks) {
    if ($rank -notin $supportedRanks) {
        throw "Unsupported global coarse rank: $rank"
    }
}
if (($InterfaceIds | Sort-Object -Unique).Count -ne
    $InterfaceIds.Count) {
    throw 'Global coarse interface ids must be unique.'
}

$baselineRows = @(Import-Csv -LiteralPath $accuracyPath)
$baseline = $baselineRows |
    Where-Object {
        $_.method -like '*M8.9*' -or
        $_.method -like '*Adaptive*'
    } |
    Select-Object -Last 1
if (-not $baseline) {
    throw 'M8.9 adaptive accuracy baseline row is missing.'
}
$baselineL2 = [double]$baseline.temperature_relative_l2
$interfaceCsv = $InterfaceIds -join ','
$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
New-Item -ItemType Directory -Path $result | Out-Null

$summaryRows = @()
$interfaceRows = @()
$stoppedAtGate = $false
Push-Location $project
try {
    foreach ($rank in $Ranks) {
        $run = Join-Path $result "rank_$rank"
        New-Item -ItemType Directory -Path $run | Out-Null
        $candidateDimension = [Math]::Max(64, 3 * $rank)
        & $exe --transient --config $configPath `
            --mor-transient-generate `
            --mor-transient-method local-port-block-arnoldi `
            --port-basis-method hybrid-randomized `
            --mor-arnoldi-moments 1 `
            --mor-interface-rank 0 `
            --randomized-port-rank 16 `
            --randomized-oversampling 5 `
            --randomized-power-iterations 1 `
            --randomized-seed 12345 `
            --residual-krylov-max-rank 4 `
            --residual-krylov-max-sweeps 2 `
            --residual-krylov-tol 1e-4 `
            --residual-krylov-block-size 4 `
            --residual-krylov-probe-mode operator-geometry `
            --residual-krylov-inner-solver woodbury-exact `
            --optimal-port-source-mode trace-only `
            --optimal-port-inner-solver woodbury-exact `
            --optimal-port-inner-tol 1e-10 `
            --optimal-port-inner-refinement-max-iters 3 `
            --optimal-port-inner-refinement-tol 1e-10 `
            --history-compression-method deterministic-rrqr `
            --history-compression-rank 64 `
            --history-compression-tolerance 1e-12 `
            --milestone8-adaptive-production `
            --global-interface-coarse-prototype `
            --global-interface-coarse-rank $rank `
            --global-interface-coarse-candidate-dimension $candidateDimension `
            --global-interface-coarse-max-iters 12 `
            --global-interface-coarse-krylov-sweeps 2 `
            --global-interface-coarse-tol 1e-8 `
            --global-interface-coarse-interface-ids $interfaceCsv `
            --mor-transient-load $modelPath `
            --mor-transient-dt 0.01 `
            --mor-transient-t-end 0.01 `
            --mor-transient-waveform mixed_frequency `
            --mor-transient-initial-mode ambient `
            --output-dir $run --fast-run
        if ($LASTEXITCODE -ne 0) {
            throw "M8.10.1 rank $rank run failed."
        }

        $coarse = Import-Csv -LiteralPath (
            Join-Path $run 'milestone8_global_coarse_prototype.csv') |
            Select-Object -First 1
        $accuracy = Import-Csv -LiteralPath (
            Join-Path $run 'milestone8_global_coarse_single_step.csv') |
            Select-Object -First 1
        $physical = @(
            Import-Csv -LiteralPath (
                Join-Path $run `
                    'milestone8_adaptive_port_physical_diagnostics.csv') |
            Where-Object {
                [int]$_.step -eq 1 -and
                [int]$_.interface_id -in $InterfaceIds
            })
        $algebraPassed =
            [double]$coarse.ritz_residual -lt 1e-8 -and
            [double]$coarse.local_coarse_orthogonality -lt 1e-10 -and
            [double]$coarse.basis_time_s -lt 300.0 -and
            [double]$coarse.peak_incremental_memory_bytes -lt 1GB
        $improvement =
            1.0 - [double]$accuracy.temperature_relative_l2 /
                $baselineL2
        $physicalPassed = $improvement -ge 0.05
        $summaryRows += [pscustomobject]@{
            case = 'rram26'
            selected_interfaces = $interfaceCsv
            local_port_rank = $accuracy.local_port_rank
            coarse_rank = $accuracy.coarse_rank
            augmented_rank = $accuracy.augmented_rank
            ritz_residual = $coarse.ritz_residual
            schur_residual = $coarse.schur_residual
            local_coarse_orthogonality =
                $coarse.local_coarse_orthogonality
            symmetry_error = $coarse.symmetry_error
            basis_time_s = $coarse.basis_time_s
            peak_incremental_memory_bytes =
                $coarse.peak_incremental_memory_bytes
            temperature_relative_l2 =
                $accuracy.temperature_relative_l2
            baseline_temperature_relative_l2 = $baselineL2
            temperature_improvement_fraction = $improvement
            max_error_k = $accuracy.max_error_k
            flux_relative_l2 = $accuracy.flux_relative_l2
            interface_residual = $accuracy.interface_residual
            full_residual = $accuracy.full_residual
            algebra_gate = if ($algebraPassed) {'passed'} else {'failed'}
            physical_gate = if ($physicalPassed) {'passed'} else {'failed'}
            corrected = 0
            snapshot_used = 0
            fom_used_for_basis = 0
            pod_used = 0
            svd_used = 0
        }
        foreach ($row in $physical) {
            $interfaceRows += [pscustomobject]@{
                coarse_rank = $rank
                interface_id = $row.interface_id
                target_dofs = $row.target_dofs
                temperature_relative_l2 =
                    $row.temperature_relative_l2
                max_nodal_error_k = $row.max_nodal_error_k
                local_full_residual = $row.local_full_residual
                flux_relative_l2 = $row.flux_relative_l2
                temperature_jump_rms_k =
                    $row.temperature_jump_rms_k
                max_relative_flux_imbalance =
                    $row.max_relative_flux_imbalance
                corrected = 0
            }
        }
        $summaryRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
            Join-Path $result 'milestone8_global_coarse_rank_sweep.csv')
        $interfaceRows | Export-Csv -NoTypeInformation -Encoding UTF8 (
            Join-Path $result `
                'milestone8_global_coarse_interface_validation.csv')
        if (-not $algebraPassed) {
            Write-Warning (
                "Rank $rank failed an algebra/resource gate; " +
                "stopping the rank sweep without discarding diagnostics.")
            $stoppedAtGate = $true
            break
        }
    }
} finally {
    Pop-Location
}

Copy-Item -LiteralPath (
    Join-Path $result 'milestone8_global_coarse_rank_sweep.csv'
) -Destination (
    Join-Path $project `
        'outputs\milestone8_global_coarse_rank_sweep.csv')
Copy-Item -LiteralPath (
    Join-Path $result `
        'milestone8_global_coarse_interface_validation.csv'
) -Destination (
    Join-Path $project `
        'outputs\milestone8_global_coarse_interface_validation.csv')
if ($stoppedAtGate) {
    Write-Warning (
        'M8.10.1 stopped at the first failed algebra/resource gate. ' +
        'No larger rank was run.')
}
