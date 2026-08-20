param(
    [string]$RramConfig = 'D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt',
    [string]$ResultsDirectory = 'results',
    [switch]$SkipBuild,
    [switch]$SkipRram,
    [switch]$RunCorrected
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$results = [System.IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'

function Invoke-Solver([string[]]$Arguments) {
    & $exe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "SIPGHeatDDM3D failed with exit code $LASTEXITCODE"
    }
}

Push-Location $project
try {
    if (-not $SkipBuild) {
        cmake -S . -B build
        if ($LASTEXITCODE -ne 0) { throw 'CMake configure failed.' }
        cmake --build build --config Release --parallel
        if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
        ctest --test-dir build -C Release --output-on-failure
        if ($LASTEXITCODE -ne 0) { throw 'CTest failed.' }
    }

    $two = Join-Path $results 'stage2b1_two_cube_repro'
    Invoke-Solver (@('--steady', '--config', 'configs\two_cube_parametric_h.txt',
        '--solvers', 'reduced-schur', '--mor-parametric-generate',
        '--mor-parametric-save', (Join-Path $two 'model'),
        '--mor-matrix-parameter', 'convection-h', '--mor-parameter-subdomain', '1',
        '--mor-parameter-region-id', '3', '--mor-parameter-min', '50',
        '--mor-parameter-max', '200', '--mor-parameter-reference', '100',
        '--mor-parameter-training-count', '3', '--mor-parameter-validation-count', '1',
        '--mor-parameter-test-count', '1', '--mor-interface-rank', '8',
        '--mor-local-rank', '8', '--mor-parametric-mode', 'pure',
        '--output-dir', $two, '--fast-run'))

    $ten = Join-Path $results 'stage2b1_ten_cube_repro'
    Invoke-Solver (@('--steady', '--config', 'configs\ten_cube_parametric_h.txt',
        '--solvers', 'reduced-schur', '--mor-parametric-generate',
        '--mor-parametric-save', (Join-Path $ten 'model'),
        '--mor-matrix-parameter', 'convection-h', '--mor-parameter-subdomain', '9',
        '--mor-parameter-region-id', '5', '--mor-parameter-min', '50',
        '--mor-parameter-max', '200', '--mor-parameter-reference', '100',
        '--mor-parameter-training-count', '5', '--mor-parameter-validation-count', '2',
        '--mor-parameter-test-count', '2', '--mor-interface-rank', '32',
        '--mor-local-rank', '32', '--mor-parametric-mode', 'pure',
        '--output-dir', $ten, '--fast-run'))
    Invoke-Solver (@('--mor-parametric-deployment-only', '--mor-parametric-load',
        (Join-Path $ten 'model'), '--mor-parameter-value', '125',
        '--output-dir', (Join-Path $results 'stage2b1_ten_cube_deployment_repro')))

    if (-not $SkipRram) {
        if (-not (Test-Path -LiteralPath $RramConfig)) {
            throw "RRAM26 config not found: $RramConfig"
        }
        $rramOffline = Join-Path $results 'stage2b1_rram26_r160_repro'
        Invoke-Solver (@('--steady', '--config', $RramConfig,
            '--solvers', 'reduced-schur', '--mor-parametric-generate',
            '--mor-parametric-save', (Join-Path $rramOffline 'model'),
            '--mor-matrix-parameter', 'material-k', '--mor-parameter-subdomain', '4',
            '--mor-parameter-region-id', '160', '--mor-parameter-min', '0.5',
            '--mor-parameter-max', '2.0', '--mor-parameter-reference', '1.0',
            '--mor-parameter-training-count', '5', '--mor-parameter-validation-count', '1',
            '--mor-parameter-test-count', '1', '--mor-interface-rank', '160',
            '--mor-local-rank', '160', '--mor-parametric-mode', 'pure',
            '--output-dir', $rramOffline, '--fast-run'))

        $rramAccepted = Join-Path $results 'stage2b1_rram26_r130_validation_repro'
        Invoke-Solver (@('--steady', '--config', $RramConfig,
            '--solvers', 'reduced-schur', '--mor-parametric-load',
            (Join-Path $rramOffline 'model'), '--mor-matrix-parameter', 'material-k',
            '--mor-parameter-subdomain', '4', '--mor-parameter-region-id', '160',
            '--mor-parameter-min', '0.5', '--mor-parameter-max', '2.0',
            '--mor-parameter-reference', '1.0', '--mor-parameter-training-count', '5',
            '--mor-parameter-validation-count', '1', '--mor-parameter-test-count', '1',
            '--mor-interface-rank', '130', '--mor-local-rank', '130',
            '--mor-parametric-mode', 'pure', '--output-dir', $rramAccepted, '--fast-run'))
        $rramDeploymentModel = Join-Path $results 'stage2b1_rram26_r130_model\model'
        Invoke-Solver (@('--mor-parametric-deployment-only', '--mor-parametric-load',
            (Join-Path $rramOffline 'model'), '--mor-interface-rank', '130',
            '--mor-local-rank', '130', '--mor-parametric-save', $rramDeploymentModel,
            '--mor-parameter-value', '1.0',
            '--output-dir', (Join-Path $results 'stage2b1_rram26_r130_model_packaging')))
        Invoke-Solver (@('--mor-parametric-deployment-only', '--mor-parametric-load',
            $rramDeploymentModel, '--mor-deployment-rhs-count', '100', '--mor-seed', '20260721',
            '--output-dir', (Join-Path $results 'stage2b1_rram26_r130_deployment_repro')))

        if ($RunCorrected) {
            Invoke-Solver (@('--steady', '--config', $RramConfig,
                '--solvers', 'reduced-schur', '--mor-parametric-load',
                (Join-Path $rramOffline 'model'), '--mor-matrix-parameter', 'material-k',
                '--mor-parameter-subdomain', '4', '--mor-parameter-region-id', '160',
                '--mor-parameter-min', '0.5', '--mor-parameter-max', '2.0',
                '--mor-parameter-reference', '1.0', '--mor-parameter-value', '0.625',
                '--mor-parameter-training-count', '5', '--mor-parameter-validation-count', '0',
                '--mor-parameter-test-count', '0', '--mor-interface-rank', '130',
                '--mor-local-rank', '130', '--mor-parametric-mode', 'corrected',
                '--pcg-tolerance', '1e-8', '--max-pcg-iterations', '500',
                '--gmres-restart', '100',
                '--output-dir', (Join-Path $results 'stage2b1_rram26_r130_corrected_repro'),
                '--fast-run'))
        }
    }
} finally {
    Pop-Location
}
