param(
    [string]$Config =
        'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$LocalPortModel =
        'D:\AI agent\kimi\DDM_Schur_m8_runs\milestone8_adaptive_port_basis_rram26\local_port_basis.bin',
    [string]$SmallCaseGate =
        'outputs\milestone8_explicit_vs_matrixfree_coarse.csv',
    [string]$ResultsDirectory =
        'results\milestone8_rram26_exact_coarse_corridor',
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$gatePath = (Resolve-Path -LiteralPath $SmallCaseGate).Path
$gateRows = @(Import-Csv -LiteralPath $gatePath)
$failedSmallCase = @($gateRows | Where-Object {
    $_.case -in @('two-cube', 'ten-cube') -and
    $_.algebra_gate -ne 'passed'
})
if ($failedSmallCase.Count -gt 0) {
    throw (
        'RRAM26 exact corridor is not authorized: ' +
        'a required small-case algebra gate failed.')
}
if ($gateRows.Count -lt 2) {
    throw 'RRAM26 exact corridor is not authorized: gate rows are missing.'
}

$configPath = (Resolve-Path -LiteralPath $Config).Path
$modelPath = (Resolve-Path -LiteralPath $LocalPortModel).Path
$result = if ([IO.Path]::IsPathRooted($ResultsDirectory)) {
    [IO.Path]::GetFullPath($ResultsDirectory)
} else {
    [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite exact-corridor results: $result"
}
New-Item -ItemType Directory -Path $result | Out-Null
$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'

Push-Location $project
try {
    & $exe --transient --config $configPath `
        --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --port-basis-method hybrid-randomized `
        --mor-arnoldi-moments 1 --mor-interface-rank 0 `
        --randomized-port-rank 16 --randomized-oversampling 5 `
        --randomized-power-iterations 1 --randomized-seed 12345 `
        --residual-krylov-max-rank 4 `
        --residual-krylov-max-sweeps 2 `
        --residual-krylov-tol 1e-4 `
        --residual-krylov-block-size 4 `
        --residual-krylov-probe-mode operator-geometry `
        --residual-krylov-inner-solver woodbury-exact `
        --optimal-port-source-mode trace-only `
        --optimal-port-inner-solver woodbury-exact `
        --optimal-port-inner-tol 1e-10 `
        --history-compression-method deterministic-rrqr `
        --history-compression-rank 64 `
        --history-compression-tolerance 1e-12 `
        --milestone8-adaptive-production `
        --global-interface-coarse-prototype `
        --global-interface-coarse-inverse-mode exact-pcg `
        --global-interface-coarse-rank 4 `
        --global-interface-coarse-candidate-dimension 16 `
        --global-interface-coarse-max-iters 3 `
        --global-interface-coarse-krylov-sweeps 2 `
        --global-interface-coarse-tol 1e-8 `
        --global-interface-coarse-inner-max-iters 1000 `
        --global-interface-coarse-inner-tol 1e-10 `
        --global-interface-coarse-interface-ids '13,23,18,4,10' `
        --mor-transient-load $modelPath `
        --mor-transient-dt 0.01 `
        --mor-transient-t-end 0.01 `
        --mor-transient-waveform mixed_frequency `
        --mor-transient-initial-mode ambient `
        --output-dir $result --fast-run
    if ($LASTEXITCODE -ne 0) {
        throw 'RRAM26 exact rank-4 corridor run failed.'
    }
} finally {
    Pop-Location
}
