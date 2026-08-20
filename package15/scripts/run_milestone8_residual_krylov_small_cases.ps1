param(
    [string]$Exe = 'build\Release\SIPGHeatDDM3D.exe',
    [string]$WorkDirectory =
        'results\milestone8_residual_krylov_small_cases',
    [string]$OutputsDirectory = 'outputs',
    [int]$CpuThreads = 8,
    [switch]$SkipBuild
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exePath = [IO.Path]::GetFullPath((Join-Path $project $Exe))
$work = [IO.Path]::GetFullPath((Join-Path $project $WorkDirectory))
$outputs = [IO.Path]::GetFullPath((Join-Path $project $OutputsDirectory))
if (Test-Path -LiteralPath $work) {
    throw "Refusing to overwrite existing work directory: $work"
}

function Invoke-Case(
    [string]$CaseName, [string]$Config,
    [string]$Method, [int]$Rank) {
    $name = "${CaseName}_${Method}_r${Rank}"
    $directory = Join-Path $work $name
    New-Item -ItemType Directory -Force -Path $directory | Out-Null
    $arguments = @(
        '--transient','--config',$Config,
        '--mor-transient-generate',
        '--mor-transient-method','local-port-block-arnoldi',
        '--port-basis-method',$Method,
        '--mor-arnoldi-moments','2',
        '--mor-transient-dt','0.1',
        '--mor-transient-t-end','0.1',
        '--mor-transient-waveform','single_step',
        '--mor-transient-initial-mode','ambient',
        '--mor-transient-output','max-temperature',
        '--output-dir',$directory,'--fast-run')
    if ($Method -in @('mandatory-only','residual-krylov')) {
        $arguments += @(
            '--residual-krylov-max-rank',"$Rank",
            '--residual-krylov-max-sweeps','2',
            '--residual-krylov-tol','1e-12',
            '--residual-krylov-block-size','2',
            '--residual-krylov-probe-mode','operator-geometry',
            '--residual-krylov-inner-solver','woodbury-exact',
            '--optimal-port-inner-solver','woodbury-exact',
            '--optimal-port-inner-tol','1e-10',
            '--optimal-port-inner-refinement-max-iters','3',
            '--optimal-port-inner-refinement-tol','1e-10')
    } elseif ($Method -eq 'optimal-transfer') {
        $arguments += @(
            '--optimal-port-rank','8',
            '--optimal-port-inner-solver','direct',
            '--optimal-port-ablation','original-mandatory-trace',
            '--optimal-port-eigen-max-iters','80',
            '--optimal-port-eigen-tol','1e-8')
    } elseif ($Method -eq 'steklov-schur') {
        $arguments += @('--local-port-rank','8')
    } elseif ($Method -eq 'port-pod') {
        $arguments += @('--local-port-rank','8')
    }
    & $exePath @arguments | ForEach-Object { Write-Host $_ }
    if ($LASTEXITCODE -ne 0) {
        throw "Small-case validation failed: $name"
    }
    $summary = Import-Csv (
        Join-Path $directory 'local_dynamic_schur_summary.csv')
    $diagnostics = $null
    if ($Method -in @('mandatory-only','residual-krylov')) {
        $diagnostics = @(Import-Csv (
            Join-Path $directory `
                'residual_krylov_interface_diagnostics.csv'))
    }
    $targetSolves = 0
    $accepted = 0
    $initialResidual = 0.0
    $finalResidual = 0.0
    $incrementalMemory = 0
    if ($null -ne $diagnostics) {
        $targetSolves = ($diagnostics |
            Measure-Object target_solve_count -Sum).Sum
        $accepted = ($diagnostics |
            Measure-Object accepted_enrichment_rank -Sum).Sum
        $initialResidual = ($diagnostics |
            Measure-Object initial_max_probe_residual -Maximum).Maximum
        $finalResidual = ($diagnostics |
            Measure-Object final_max_probe_residual -Maximum).Maximum
        $incrementalMemory = ($diagnostics |
            Measure-Object peak_incremental_memory_bytes -Maximum).Maximum
    } elseif ($Method -eq 'optimal-transfer') {
        $inner = @(Import-Csv (
            Join-Path $directory 'optimal_port_inner_solver.csv'))
        $targetSolves = ($inner |
            Measure-Object solve_calls -Sum).Sum
        $incrementalMemory = ($inner |
            Measure-Object peak_incremental_memory_bytes -Maximum).Maximum
    }
    return [pscustomobject]@{
        case = $CaseName
        method = $Method
        requested_enrichment_rank = $Rank
        accepted_enrichment_rank = $accepted
        port_rank = [int]$summary.port_dimension
        target_solve_count = [int]$targetSolves
        initial_max_probe_residual = [double]$initialResidual
        final_max_probe_residual = [double]$finalResidual
        relative_l2 = [double]$summary.space_time_relative_l2
        max_node_error_k = [double]$summary.maximum_absolute_k
        flux_relative_l2 =
            [double]$summary.maximum_fom_rom_flux_relative_l2
        interface_residual =
            [double]$summary.maximum_interface_relative_residual
        full_residual = [double]$summary.maximum_full_residual
        basis_build_time_s = [double]$summary.port_basis_seconds
        total_time_s = [double]$summary.total_seconds
        peak_incremental_memory_bytes = [uint64]$incrementalMemory
        corrected = $false
        snapshot_used = [int]$summary.port_snapshot_used
        pod_used = [int]$summary.port_fom_used_for_basis
        svd_used = [int]$summary.port_fom_used_for_basis
    }
}

Push-Location $project
try {
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
    $rows = @()
    foreach ($case in @(
        @{ Name='two-cube'; Config='configs\two_cube_parametric_h.txt' },
        @{ Name='ten-cube'; Config='configs\ten_cube_parametric_h.txt' })) {
        $rows += Invoke-Case $case.Name $case.Config 'mandatory-only' 0
        foreach ($rank in @(1,2,4)) {
            $rows += Invoke-Case $case.Name $case.Config `
                'residual-krylov' $rank
        }
        $rows += Invoke-Case $case.Name $case.Config `
            'optimal-transfer' 0
        $rows += Invoke-Case $case.Name $case.Config 'steklov-schur' 0
        $rows += Invoke-Case $case.Name $case.Config 'port-pod' 0
        $rows += Invoke-Case $case.Name $case.Config 'full-interface' 0
    }
    $rows | Export-Csv -NoTypeInformation -Encoding UTF8 `
        (Join-Path $outputs `
            'milestone8_residual_krylov_two_ten_cube.csv')
    $rows | Select-Object case,method,requested_enrichment_rank,
        accepted_enrichment_rank,port_rank,relative_l2,
        max_node_error_k,flux_relative_l2,interface_residual,
        full_residual,corrected,snapshot_used,pod_used,svd_used |
        Export-Csv -NoTypeInformation -Encoding UTF8 `
            (Join-Path $outputs `
                'milestone8_residual_krylov_accuracy.csv')
    $rows | Select-Object case,method,requested_enrichment_rank,
        target_solve_count,basis_build_time_s,total_time_s |
        Export-Csv -NoTypeInformation -Encoding UTF8 `
            (Join-Path $outputs `
                'milestone8_residual_krylov_timing.csv')
    $rows | Select-Object case,method,requested_enrichment_rank,
        peak_incremental_memory_bytes |
        Export-Csv -NoTypeInformation -Encoding UTF8 `
            (Join-Path $outputs `
                'milestone8_residual_krylov_memory.csv')
    $rows | Where-Object {
        $_.method -in @('residual-krylov','optimal-transfer') } |
        Select-Object case,method,requested_enrichment_rank,
            accepted_enrichment_rank,target_solve_count,
            basis_build_time_s,port_rank,interface_residual,
            flux_relative_l2,full_residual,
            peak_incremental_memory_bytes |
        Export-Csv -NoTypeInformation -Encoding UTF8 `
            (Join-Path $outputs `
                'milestone8_residual_krylov_vs_optimal_transfer.csv')
} finally {
    Pop-Location
}
