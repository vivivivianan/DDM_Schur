param(
    [string]$RramConfig =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$TopologyAuditCsv =
        'outputs\milestone8_large_case_topology_audit.csv',
    [string]$ResultsDirectory =
        'results\milestone8_rram26_residual_krylov_all_interfaces',
    [int]$CpuThreads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$gate = Join-Path $project `
    'outputs\milestone8_rram26_residual_krylov_stage_b_gate.csv'
if (-not (Test-Path -LiteralPath $gate -PathType Leaf) -or
    -not ((Get-Content -LiteralPath $gate -Raw) -match
        'stage_b_passed')) {
    throw 'Stage C is locked until all 12 Stage-B rows pass.'
}
$result = [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite existing Stage-C directory: $result"
}
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$env:OMP_NUM_THREADS = "$CpuThreads"
$env:MKL_NUM_THREADS = "$CpuThreads"
$env:MKL_DYNAMIC = 'FALSE'
Push-Location $project
try {
    & $exe --transient --config $RramConfig `
        --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --port-basis-method residual-krylov `
        --mor-arnoldi-moments 1 `
        --residual-krylov-max-rank 4 `
        --residual-krylov-max-sweeps 2 `
        --residual-krylov-tol 1e-4 `
        --residual-krylov-block-size 2 `
        --residual-krylov-probe-mode operator-geometry `
        --residual-krylov-inner-solver woodbury-exact `
        --residual-krylov-all-interface-basis `
        --optimal-port-topology-audit-csv $TopologyAuditCsv `
        --optimal-port-inner-solver woodbury-exact `
        --optimal-port-inner-tol 1e-10 `
        --optimal-port-inner-refinement-max-iters 3 `
        --optimal-port-inner-refinement-tol 1e-10 `
        --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
        --mor-transient-waveform single_step `
        --mor-transient-initial-mode ambient `
        --output-dir $result --fast-run
    if ($LASTEXITCODE -ne 0) {
        throw 'RRAM26 all-interface residual-Krylov build failed.'
    }
    Copy-Item -LiteralPath (
        Join-Path $result `
            'milestone8_rram26_residual_krylov_all_interfaces.csv') `
        -Destination (
            Join-Path $project `
                'outputs\milestone8_rram26_residual_krylov_all_interfaces.csv')
} finally {
    Pop-Location
}
