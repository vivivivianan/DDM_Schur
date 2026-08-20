param(
    [string]$RramConfig =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$ChipletConfig =
        'D:\CPP\TEST_CHATGPT\chiplet_model\case_chiplet_config_horizontal.txt',
    [string]$ResultsDirectory = 'results\milestone8_readiness',
    [string]$OutputsDirectory = 'outputs',
    [switch]$SkipBuild,
    [switch]$SkipTenCube,
    [switch]$SkipLargeCaseAudit
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$results = [System.IO.Path]::GetFullPath(
    (Join-Path $project $ResultsDirectory))
$outputs = [System.IO.Path]::GetFullPath(
    (Join-Path $project $OutputsDirectory))
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$tenConfig = Join-Path $project 'configs\ten_cube_parametric_h.txt'
$rankOutput = Join-Path $outputs 'milestone8_ten_cube_rank_sweep.csv'
$ablationOutput = Join-Path $outputs 'milestone8_ten_cube_ablation.csv'
$topologyOutput = Join-Path $outputs 'milestone8_large_case_topology_audit.csv'
$memoryOutput = Join-Path $outputs 'milestone8_large_case_memory_estimate.csv'

$productionEigenTolerance = 1.0e-8
$productionInnerTolerance = 1.0e-10
$productionEigenMaximumIterations = 1000
$ablationRankPerPort = 12

function Assert-File([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        throw "Required file is missing: $Path"
    }
}

function Read-One([string]$Path) {
    Assert-File $Path
    return Import-Csv -LiteralPath $Path | Select-Object -First 1
}

function Sum-Column($Rows, [string]$Column) {
    if (@($Rows).Count -eq 0) { return 0 }
    return [double](($Rows | Measure-Object -Property $Column -Sum).Sum)
}

function Max-Column($Rows, [string]$Column) {
    if (@($Rows).Count -eq 0) { return 0 }
    return [double](($Rows | Measure-Object -Property $Column -Maximum).Maximum)
}

function Join-Unique($Rows, [string]$Column) {
    if (@($Rows).Count -eq 0) { return 'not_applicable' }
    return (@($Rows | ForEach-Object { $_.$Column } |
        Where-Object { $_ -and $_ -ne '' } |
        Sort-Object -Unique) -join ';')
}

function Invoke-ReadinessCase(
    [string]$Name,
    [string]$Config,
    [string]$Output,
    [string[]]$Arguments) {
    if (Test-Path -LiteralPath $Output) {
        throw "Refusing to overwrite an existing readiness run: $Output"
    }
    New-Item -ItemType Directory -Path $Output | Out-Null
    $log = Join-Path $Output 'console.log'
    Write-Host "Running $Name"
    $console = @(& $exe --transient --config $Config `
        --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --mor-arnoldi-moments 1 `
        --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
        --mor-transient-waveform single_step `
        --mor-transient-initial-mode ambient `
        --output-dir $Output --fast-run @Arguments 2>&1)
    $exitCode = $LASTEXITCODE
    $console | Set-Content -Encoding UTF8 -LiteralPath $log
    if ($exitCode -ne 0) {
        throw "$Name failed with exit code $exitCode. See $log"
    }
}

function Get-OptimalRunRow(
    [string]$RunName,
    [string]$Contribution,
    [int]$RequestedRank,
    [int]$EigenMaximumIterations,
    [string]$Directory) {
    $summary = Read-One (
        Join-Path $Directory 'local_dynamic_schur_summary.csv')
    $ports = @(Import-Csv -LiteralPath (
        Join-Path $Directory 'optimal_port_rank_by_interface.csv'))
    $operators = @(Import-Csv -LiteralPath (
        Join-Path $Directory 'optimal_port_operator_diagnostics.csv'))
    $inner = @(Import-Csv -LiteralPath (
        Join-Path $Directory 'optimal_port_inner_solver.csv'))
    $timing = Read-One (
        Join-Path $Directory 'optimal_port_timing.csv')
    $residualPath = Join-Path $Directory 'optimal_port_eigenpair_residual.csv'
    $eigenResiduals = if ((Get-Item -LiteralPath $residualPath).Length -gt 64) {
        @(Import-Csv -LiteralPath $residualPath)
    } else {
        @()
    }
    $transferColumn = if (
        $ports[0].PSObject.Properties.Name -contains
            'converged_transfer_rank') {
        'converged_transfer_rank'
    } else {
        'transfer_modes'
    }
    $mandatoryColumn = if (
        $ports[0].PSObject.Properties.Name -contains 'mandatory_rank') {
        'mandatory_rank'
    } else {
        'mandatory_modes'
    }
    $activeTransferInterfaces = @($ports |
        Where-Object { [int]($_.$transferColumn) -gt 0 })
    $activeIds = @($activeTransferInterfaces |
        ForEach-Object { [int]$_.interface_id })
    $failedEigen = @($operators | Where-Object {
        $activeIds -contains [int]$_.interface_id -and
        [int]$_.eigen_converged -ne 1 })
    $fallback = @($inner |
        Where-Object { [int]$_.fallback_triggered -ne 0 })
    $fallbackReasons = @($fallback |
        ForEach-Object { $_.fallback_reason } |
        Where-Object { $_ -and $_ -ne '' } |
        Sort-Object -Unique) -join ';'
    $actualSolver = Join-Unique $inner 'actual_solver'
    $iterativeMethod = if ($actualSolver -match 'fgmres') {
        'FGMRES'
    } elseif ($actualSolver -match 'pcg') {
        'PCG'
    } else {
        'not_used'
    }
    [pscustomobject]@{
        run_name = $RunName
        method = 'optimal-transfer'
        contribution = $Contribution
        requested_rank_per_port = $RequestedRank
        physical_interfaces = $ports.Count
        coarse_dimension = [int]$summary.port_dimension
        mandatory_modes = [int](Sum-Column $ports $mandatoryColumn)
        transfer_modes = [int](Sum-Column $ports $transferColumn)
        eigensolver_tolerance = $productionEigenTolerance
        eigensolver_max_iterations = $EigenMaximumIterations
        eigen_iterations_maximum = [int](Max-Column $operators 'eigen_iterations')
        eigen_operator_applies = [int](Sum-Column $operators 'eigen_operator_applies')
        eigen_converged = if ($activeIds.Count -eq 0) {
            'not_applicable'
        } elseif ($failedEigen.Count -eq 0) {
            'yes'
        } else {
            'no'
        }
        maximum_eigenpair_relative_residual = Max-Column `
            $eigenResiduals 'relative_residual'
        inner_tolerance = $productionInnerTolerance
        requested_solver = Join-Unique $inner 'requested_solver'
        actual_solver = $actualSolver
        pcg_or_fgmres = $iterativeMethod
        solver_path = Join-Unique $inner 'solver_path'
        inner_solver_status = Join-Unique $inner 'status'
        fallback_triggered = if ($fallback.Count -gt 0) { 1 } else { 0 }
        fallback_reason = $fallbackReasons
        inner_iterations_total = [int](Sum-Column $inner 'total_iterations')
        inner_iterations_maximum = [int](Max-Column $inner 'max_iterations')
        inner_setup_seconds = Sum-Column $inner 'setup_seconds'
        repeated_inner_solve_seconds = Sum-Column $inner 'total_solve_seconds'
        maximum_inner_relative_residual = Max-Column `
            $inner 'final_relative_residual'
        setup_seconds =
            [double]$summary.port_snapshot_seconds +
            [double]$summary.port_basis_seconds +
            [double]$summary.port_schur_assembly_seconds +
            [double]$summary.port_schur_factor_seconds
        port_setup_seconds = [double]$timing.total_port_offline_seconds
        eigensolve_seconds = [double]$timing.eigen_solve_seconds
        interface_iterations_maximum =
            [int]$summary.interface_iterations_maximum
        interface_solve_seconds = [double]$summary.interface_solve_seconds
        interface_residual = [double]$summary.maximum_interface_relative_residual
        global_residual = [double]$summary.maximum_full_residual
        maximum_temperature_difference_k = [double]$summary.maximum_absolute_k
        relative_l2 = [double]$summary.space_time_relative_l2
        total_seconds = [double]$summary.total_seconds
        status = $summary.status
    }
}

function Get-GeneralRunRow(
    [string]$RunName,
    [string]$Method,
    [string]$Contribution,
    [int]$RequestedRank,
    [int]$PhysicalInterfaces,
    [string]$Directory) {
    $summary = Read-One (
        Join-Path $Directory 'local_dynamic_schur_summary.csv')
    $rankPath = Join-Path $Directory 'local_port_rank_by_interface.csv'
    $ports = if (Test-Path -LiteralPath $rankPath) {
        @(Import-Csv -LiteralPath $rankPath)
    } else {
        @()
    }
    $mandatory = if ($ports.Count -gt 0 -and
        $ports[0].PSObject.Properties.Name -contains 'mandatory_rank') {
        [int](Sum-Column $ports 'mandatory_rank')
    } elseif ($ports.Count -gt 0 -and
        $ports[0].PSObject.Properties.Name -contains 'mandatory_modes') {
        [int](Sum-Column $ports 'mandatory_modes')
    } else {
        0
    }
    $spectral = if ($ports.Count -gt 0 -and
        $ports[0].PSObject.Properties.Name -contains
            'converged_transfer_rank') {
        [int](Sum-Column $ports 'converged_transfer_rank')
    } elseif ($ports.Count -gt 0 -and
        $ports[0].PSObject.Properties.Name -contains 'spectral_modes') {
        [int](Sum-Column $ports 'spectral_modes')
    } else {
        0
    }
    [pscustomobject]@{
        run_name = $RunName
        method = $Method
        contribution = $Contribution
        requested_rank_per_port = $RequestedRank
        physical_interfaces = $PhysicalInterfaces
        coarse_dimension = if ([int]$summary.port_reduction -ne 0) {
            [int]$summary.port_dimension
        } else {
            [int]$summary.full_interface_dofs
        }
        mandatory_modes = $mandatory
        transfer_modes = $spectral
        eigensolver_tolerance = 'not_applicable'
        eigensolver_max_iterations = 'not_applicable'
        eigen_iterations_maximum = 0
        eigen_operator_applies = 0
        eigen_converged = 'not_applicable'
        maximum_eigenpair_relative_residual = 0
        inner_tolerance = 'not_applicable'
        requested_solver = 'not_applicable'
        actual_solver = 'not_applicable'
        pcg_or_fgmres = 'not_applicable'
        solver_path = 'not_applicable'
        inner_solver_status = 'not_applicable'
        fallback_triggered = 0
        fallback_reason = ''
        inner_iterations_total = 0
        inner_iterations_maximum = 0
        inner_setup_seconds = 0
        repeated_inner_solve_seconds = 0
        maximum_inner_relative_residual = 0
        setup_seconds = if ([int]$summary.port_reduction -ne 0) {
            [double]$summary.port_snapshot_seconds +
            [double]$summary.port_basis_seconds +
            [double]$summary.port_schur_assembly_seconds +
            [double]$summary.port_schur_factor_seconds
        } else {
            [double]$summary.dynamic_schur_setup_seconds +
            [double]$summary.dynamic_schur_factor_seconds
        }
        port_setup_seconds = [double]$summary.port_basis_seconds
        eigensolve_seconds = 0
        interface_iterations_maximum =
            [int]$summary.interface_iterations_maximum
        interface_solve_seconds = [double]$summary.interface_solve_seconds
        interface_residual = [double]$summary.maximum_interface_relative_residual
        global_residual = [double]$summary.maximum_full_residual
        maximum_temperature_difference_k = [double]$summary.maximum_absolute_k
        relative_l2 = [double]$summary.space_time_relative_l2
        total_seconds = [double]$summary.total_seconds
        status = $summary.status
    }
}

function Assert-TopologyOnly([string]$Directory) {
    Assert-File (Join-Path $Directory 'optimal_port_topology_audit.csv')
    Assert-File (
        Join-Path $Directory 'optimal_port_topology_memory_estimate.csv')
    $forbidden = @(
        'local_dynamic_schur_summary.csv',
        'optimal_port_rank_by_interface.csv',
        'optimal_port_eigenpair_residual.csv',
        'transient_temperature_comparison.csv')
    foreach ($name in $forbidden) {
        if (Test-Path -LiteralPath (Join-Path $Directory $name)) {
            throw "Topology-only guard failed; unexpected solve artifact: $name"
        }
    }
}

Push-Location $project
try {
    New-Item -ItemType Directory -Force -Path $results,$outputs | Out-Null
    Assert-File $tenConfig
    if (-not $SkipBuild) {
        cmake -S . -B build -DBUILD_TESTING=ON
        if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed.' }
        cmake --build build --config Release --target SIPGHeatDDM3D --parallel 2
        if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
    }
    Assert-File $exe

    if (-not $SkipTenCube) {
        $rankRows = @()
        foreach ($rank in @(6,8,10,12,16,24)) {
            $directory = Join-Path $results "ten_rank_$rank"
            Invoke-ReadinessCase "ten-cube production rank $rank" `
                $tenConfig $directory @(
                    '--port-basis-method','optimal-transfer',
                    '--optimal-port-ablation','mandatory-transfer',
                    '--optimal-port-rank',"$rank",
                    '--optimal-port-inner-solver','auto',
                    '--optimal-port-inner-tol',"$productionInnerTolerance",
                    '--optimal-port-inner-max-iters','4000',
                    '--optimal-port-eigen-tol',"$productionEigenTolerance",
                    '--optimal-port-eigen-max-iters',
                        "$productionEigenMaximumIterations")
            $rankRows += Get-OptimalRunRow "rank-$rank" `
                'mandatory + transfer' $rank `
                $productionEigenMaximumIterations $directory
        }
        $rankRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
            -LiteralPath $rankOutput

        $ablationRows = @()
        $optimalAblations = @(
            [pscustomobject]@{
                Name='mandatory-only'; Cli='mandatory-only'
                Contribution='constant + geometry + particular'},
            [pscustomobject]@{
                Name='transfer-only'; Cli='transfer-only'
                Contribution='transfer only'},
            [pscustomobject]@{
                Name='constant-only'; Cli='constant-only'
                Contribution='constant only'},
            [pscustomobject]@{
                Name='geometry + particular'; Cli='geometry-particular'
                Contribution='geometry + particular'},
            [pscustomobject]@{
                Name='mandatory + transfer'; Cli='mandatory-transfer'
                Contribution='mandatory + transfer'})
        foreach ($ablation in $optimalAblations) {
            $slug = $ablation.Cli.Replace('-', '_')
            $directory = Join-Path $results "ten_ablation_$slug"
            Invoke-ReadinessCase "ten-cube ablation $($ablation.Name)" `
                $tenConfig $directory @(
                    '--port-basis-method','optimal-transfer',
                    '--optimal-port-ablation',$ablation.Cli,
                    '--optimal-port-rank',"$ablationRankPerPort",
                    '--optimal-port-inner-solver','auto',
                    '--optimal-port-inner-tol',"$productionInnerTolerance",
                    '--optimal-port-inner-max-iters','4000',
                    '--optimal-port-eigen-tol',"$productionEigenTolerance",
                    '--optimal-port-eigen-max-iters',
                        "$productionEigenMaximumIterations")
            $row = Get-OptimalRunRow $ablation.Name `
                $ablation.Contribution $ablationRankPerPort `
                $productionEigenMaximumIterations $directory
            $ablationRows += $row
        }

        $steklov = Join-Path $results 'ten_ablation_steklov_schur'
        Invoke-ReadinessCase 'ten-cube ablation steklov-schur' `
            $tenConfig $steklov @(
                '--port-basis-method','steklov-schur',
                '--local-port-rank',"$ablationRankPerPort")
        $ablationRows += Get-GeneralRunRow 'steklov-schur' `
            'steklov-schur' 'local Schur eigenmodes' `
            $ablationRankPerPort 9 $steklov

        $portPod = Join-Path $results 'ten_ablation_m7_port_pod'
        Invoke-ReadinessCase 'ten-cube ablation M7 port-pod' `
            $tenConfig $portPod @(
                '--port-basis-method','port-pod',
                '--local-port-rank',"$ablationRankPerPort",
                '--local-port-enrichment-rounds','0')
        $ablationRows += Get-GeneralRunRow 'M7 port-pod' `
            'port-pod' 'snapshot/POD baseline' `
            $ablationRankPerPort 9 $portPod

        $full = Join-Path $results 'ten_ablation_full_interface'
        Invoke-ReadinessCase 'ten-cube ablation full-interface' `
            $tenConfig $full @(
                '--port-basis-method','full-interface')
        $ablationRows += Get-GeneralRunRow 'full-interface' `
            'full-interface' 'unreduced interface reference' 0 9 $full

        $ablationRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
            -LiteralPath $ablationOutput
    }

    if (-not $SkipLargeCaseAudit) {
        Assert-File $RramConfig
        Assert-File $ChipletConfig
        $twoCubeAudit = Join-Path $results 'two_cube_topology_only'
        Invoke-ReadinessCase 'two-cube topology audit control' `
            (Join-Path $project 'configs\two_cube_parametric_h.txt') `
            $twoCubeAudit @(
                '--port-basis-method','optimal-transfer',
                '--optimal-port-topology-audit',
                '--optimal-port-rank',"$ablationRankPerPort",
                '--optimal-port-inner-solver','auto')
        Assert-TopologyOnly $twoCubeAudit
        $twoCubeTopology = @(Import-Csv -LiteralPath (
            Join-Path $twoCubeAudit 'optimal_port_topology_audit.csv'))
        if ($twoCubeTopology.Count -ne 1 -or
            [int]$twoCubeTopology[0].source_empty -ne 1) {
            throw 'Two-cube topology control no longer has one empty transfer source.'
        }

        $auditCases = @(
            [pscustomobject]@{
                Name='RRAM26'; Config=$RramConfig
                Directory=(Join-Path $results 'rram26_topology_only')},
            [pscustomobject]@{
                Name='Chiplet'; Config=$ChipletConfig
                Directory=(Join-Path $results 'chiplet_topology_only')})
        foreach ($case in $auditCases) {
            Invoke-ReadinessCase "$($case.Name) topology audit only" `
                $case.Config $case.Directory @(
                    '--port-basis-method','optimal-transfer',
                    '--optimal-port-topology-audit',
                    '--optimal-port-rank',"$ablationRankPerPort",
                    '--optimal-port-inner-solver','auto')
            Assert-TopologyOnly $case.Directory
        }
        $topologyRows = foreach ($case in $auditCases) {
            Import-Csv -LiteralPath (
                Join-Path $case.Directory 'optimal_port_topology_audit.csv')
        }
        $memoryRows = foreach ($case in $auditCases) {
            Import-Csv -LiteralPath (
                Join-Path $case.Directory `
                    'optimal_port_topology_memory_estimate.csv')
        }
        $topologyRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
            -LiteralPath $topologyOutput
        $memoryRows | Export-Csv -NoTypeInformation -Encoding UTF8 `
            -LiteralPath $memoryOutput
    }
} finally {
    Pop-Location
}

Write-Host 'M8.5 readiness runs completed without large-case eigensolve/transient.'
