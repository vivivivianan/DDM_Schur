param(
    [string]$BuildDirectory = ".\build",
    [string]$ResultsDirectory = ".\results",
    [string]$RramConfig = "D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt",
    [string]$ChipletConfig = "D:\CPP\TEST_CHATGPT\chiplet_model\case_chiplet_config_horizontal.txt",
    [switch]$SkipBuild,
    [switch]$SkipLargeCases
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
Push-Location $repo
try {
    if (-not $SkipBuild) {
        cmake --build $BuildDirectory --config Release -- /m:1 /clp:ErrorsOnly
        if ($LASTEXITCODE -ne 0) { throw "Release build failed." }
    }
    $exe = Join-Path $BuildDirectory "Release\SIPGHeatDDM3D.exe"
    if (-not (Test-Path -LiteralPath $exe)) { throw "Executable not found: $exe" }

    function Invoke-ExactResponseCase {
        param(
            [string]$Name,
            [string]$Config,
            [int]$InterfaceRank,
            [int]$Training,
            [int]$Validation,
            [int]$Test,
            [double]$Tolerance
        )
        $output = Join-Path $ResultsDirectory $Name
        $model = Join-Path $output "mor_model"
        & $exe --steady --config $Config --output-dir $output `
            --solvers reduced-schur --mor-generate-model --mor-save-model $model `
            --mor-rank $InterfaceRank --mor-rank-sweep $InterfaceRank `
            --mor-training-count $Training --mor-validation-count $Validation `
            --mor-test-count $Test --mor-random-seed 20260721 `
            --mor-snapshot-solver direct --mor-compare-fom `
            --reduced-schur-mode pure --mor-interior-mode exact-response `
            --mor-precompute-power-response --mor-storage-precision float64 `
            --mor-report-io-time --pcg-tolerance $Tolerance `
            --max-pcg-iterations 500 --gmres-restart 100 --fast-run
        if ($LASTEXITCODE -ne 0) { throw "Exact-response case failed: $Name" }

        $deployment = Join-Path $ResultsDirectory ($Name + "_deployment")
        & $exe --mor-deployment-only --mor-load-model $model `
            --mor-deployment-rhs-count 100 --mor-report-io-time `
            --output-dir $deployment
        if ($LASTEXITCODE -ne 0) { throw "Deployment case failed: $Name" }
    }

    Invoke-ExactResponseCase -Name "stage2a1_two_exact_repro" `
        -Config ".\configs\two_cube_schur.txt" -InterfaceRank 2 `
        -Training 12 -Validation 8 -Test 8 -Tolerance 1e-10
    Invoke-ExactResponseCase -Name "stage2a1_ten_exact_repro" `
        -Config ".\configs\ten_cube_schur.txt" -InterfaceRank 10 `
        -Training 20 -Validation 20 -Test 20 -Tolerance 1e-10

    # Stage 1, Stage 2A PARDISO-interior, and corrected-mode regressions.
    & $exe --steady --config ".\configs\two_cube_schur.txt" `
        --output-dir (Join-Path $ResultsDirectory "stage2a1_stage1_regression_repro") `
        --solvers schur,direct --pcg-tolerance 1e-10 `
        --max-pcg-iterations 500 --gmres-restart 100 --fast-run
    if ($LASTEXITCODE -ne 0) { throw "Stage 1 regression failed." }

    $twoModel = Join-Path $ResultsDirectory "stage2a1_two_exact_repro\mor_model"
    foreach ($mode in @("pure", "corrected")) {
        & $exe --steady --config ".\configs\two_cube_schur.txt" `
            --output-dir (Join-Path $ResultsDirectory ("stage2a1_" + $mode + "_regression_repro")) `
            --solvers reduced-schur --mor-load-model $twoModel `
            --mor-rank 2 --mor-rank-sweep 2 --mor-training-count 0 `
            --mor-validation-count 1 --mor-test-count 1 `
            --reduced-schur-mode $mode `
            --mor-interior-mode $(if ($mode -eq "pure") { "pardiso" } else { "exact-response" }) `
            --pcg-tolerance 1e-10 --max-pcg-iterations 500 --gmres-restart 100 --fast-run
        if ($LASTEXITCODE -ne 0) { throw "$mode regression failed." }
    }

    if (-not $SkipLargeCases) {
        if (-not (Test-Path -LiteralPath $RramConfig)) {
            throw "RRAM26 config not found: $RramConfig"
        }
        if (-not (Test-Path -LiteralPath $ChipletConfig)) {
            throw "Chiplet config not found: $ChipletConfig"
        }
        Invoke-ExactResponseCase -Name "stage2a1_rram26_exact_repro" `
            -Config $RramConfig -InterfaceRank 125 `
            -Training 133 -Validation 20 -Test 20 -Tolerance 1e-5
        Invoke-ExactResponseCase -Name "stage2a1_chiplet_exact_repro" `
            -Config $ChipletConfig -InterfaceRank 4 `
            -Training 12 -Validation 20 -Test 20 -Tolerance 1e-5

        $exactModel = Join-Path $ResultsDirectory "stage2a1_rram26_exact_repro\mor_model"
        $compressedModel = Join-Path $ResultsDirectory "stage2a1_rram26_compressed_r125_repro\mor_model"
        & $exe --mor-deployment-only --mor-load-model $exactModel `
            --mor-save-model $compressedModel --mor-interior-mode compressed-rb `
            --mor-interior-rank 125 `
            --mor-interior-rank-sweep 5,10,20,40,60,80,100,105,110,115,120,121,122,123,124,125 `
            --mor-deployment-rhs-count 100 --mor-compare-interior-modes `
            --mor-report-io-time `
            --output-dir (Join-Path $ResultsDirectory "stage2a1_rram26_compression_sweep_repro")
        if ($LASTEXITCODE -ne 0) { throw "RRAM26 compressed-rank sweep failed." }

        & $exe --mor-deployment-only --mor-load-model $compressedModel `
            --mor-deployment-rhs-count 100 --mor-report-io-time `
            --output-dir (Join-Path $ResultsDirectory "stage2a1_rram26_compressed_deployment_repro")
        if ($LASTEXITCODE -ne 0) { throw "RRAM26 compressed deployment failed." }

        foreach ($energy in @(
            @{ Name = "e8"; Tolerance = "1e-8" },
            @{ Name = "e6"; Tolerance = "1e-6" },
            @{ Name = "e4"; Tolerance = "1e-4" }
        )) {
            & $exe --mor-deployment-only --mor-load-model $exactModel `
                --mor-interior-mode compressed-rb `
                --mor-interior-energy-tolerance $energy.Tolerance `
                --mor-interior-singular-value-tolerance 0 `
                --mor-interior-rank-sweep 0 --mor-deployment-rhs-count 100 `
                --mor-compare-interior-modes --mor-report-io-time `
                --output-dir (Join-Path $ResultsDirectory ("stage2a1_rram26_energy_" + $energy.Name + "_repro"))
            if ($LASTEXITCODE -ne 0) { throw "RRAM26 energy audit failed: $($energy.Name)" }
        }

        $float32Model = Join-Path $ResultsDirectory "stage2a1_rram26_exact_float32_repro\mor_model"
        & $exe --mor-deployment-only --mor-load-model $exactModel `
            --mor-save-model $float32Model --mor-storage-precision float32 `
            --mor-deployment-rhs-count 1 --mor-report-io-time `
            --output-dir (Join-Path $ResultsDirectory "stage2a1_rram26_float32_conversion_repro")
        if ($LASTEXITCODE -ne 0) { throw "RRAM26 float32 serialization audit failed." }
        & $exe --mor-deployment-only --mor-load-model $float32Model `
            --mor-deployment-rhs-count 100 --mor-report-io-time `
            --output-dir (Join-Path $ResultsDirectory "stage2a1_rram26_float32_deployment_repro")
        if ($LASTEXITCODE -ne 0) { throw "RRAM26 float32 deployment audit failed." }
    }
}
finally {
    Pop-Location
}
