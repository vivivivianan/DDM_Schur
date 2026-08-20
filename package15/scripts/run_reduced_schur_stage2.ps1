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

    function Invoke-ReducedSchurCase {
        param(
            [string]$Name,
            [string]$Config,
            [int]$Rank,
            [int]$Training,
            [int]$Validation,
            [int]$Test,
            [double]$Tolerance
        )
        $output = Join-Path $ResultsDirectory $Name
        & $exe --steady --config $Config --output-dir $output `
            --solvers reduced-schur --mor-generate-model `
            --mor-save-model (Join-Path $output "mor_model") `
            --mor-rank $Rank --mor-rank-sweep 5,10,20,40,60,80,100,125 `
            --mor-training-count $Training --mor-validation-count $Validation `
            --mor-test-count $Test --mor-random-seed 20260721 `
            --mor-snapshot-solver auto `
            --mor-compare-fom --reduced-schur-mode pure `
            --mor-exact-interior-recovery --pcg-tolerance $Tolerance `
            --max-pcg-iterations 500 --gmres-restart 100 --fast-run
        if ($LASTEXITCODE -ne 0) { throw "Reduced Schur case failed: $Name" }
    }

    Invoke-ReducedSchurCase -Name "stage2_two_cube_repro" `
        -Config ".\configs\two_cube_schur.txt" -Rank 2 `
        -Training 12 -Validation 8 -Test 8 -Tolerance 1e-10
    Invoke-ReducedSchurCase -Name "stage2_ten_cube_repro" `
        -Config ".\configs\ten_cube_schur.txt" -Rank 10 `
        -Training 20 -Validation 20 -Test 20 -Tolerance 1e-10

    if (-not $SkipLargeCases) {
        if (-not (Test-Path -LiteralPath $RramConfig)) {
            throw "RRAM26 config not found: $RramConfig"
        }
        if (-not (Test-Path -LiteralPath $ChipletConfig)) {
            throw "Chiplet config not found: $ChipletConfig"
        }
        # Unit-channel coverage takes precedence over the requested training
        # count; the current RRAM26 input contains 125 independent channels.
        # Eight additional deterministic entries exercise representative
        # multi-source families in the training split.
        Invoke-ReducedSchurCase -Name "stage2_rram26_repro" `
            -Config $RramConfig -Rank 125 `
            -Training 133 -Validation 20 -Test 20 -Tolerance 1e-5
        Invoke-ReducedSchurCase -Name "stage2_chiplet_repro" `
            -Config $ChipletConfig -Rank 4 `
            -Training 12 -Validation 20 -Test 20 -Tolerance 1e-5
    }
}
finally {
    Pop-Location
}
