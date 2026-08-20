param(
    [string]$ResultsDirectory =
        'results\milestone8_exact_projected_coarse_small',
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$result = if ([IO.Path]::IsPathRooted($ResultsDirectory)) {
    [IO.Path]::GetFullPath($ResultsDirectory)
} else {
    [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite exact-coarse results: $result"
}
if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
    throw "Release executable not found: $exe"
}
$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
New-Item -ItemType Directory -Path $result | Out-Null

function Common-HybridArguments(
    [string]$Config,
    [string]$Output) {
    return @(
        '--transient', '--config', $Config,
        '--mor-transient-generate',
        '--mor-transient-method', 'local-port-block-arnoldi',
        '--port-basis-method', 'hybrid-randomized',
        '--mor-arnoldi-moments', '1',
        '--mor-interface-rank', '0',
        '--randomized-port-rank', '16',
        '--randomized-oversampling', '5',
        '--randomized-power-iterations', '1',
        '--randomized-seed', '12345',
        '--residual-krylov-max-rank', '4',
        '--residual-krylov-max-sweeps', '2',
        '--residual-krylov-tol', '1e-4',
        '--residual-krylov-block-size', '4',
        '--residual-krylov-probe-mode', 'operator-geometry',
        '--residual-krylov-inner-solver', 'woodbury-exact',
        '--optimal-port-inner-solver', 'woodbury-exact',
        '--optimal-port-inner-tol', '1e-10',
        '--history-compression-method', 'deterministic-rrqr',
        '--history-compression-rank', '64',
        '--history-compression-tolerance', '1e-12',
        '--mor-transient-dt', '0.1',
        '--mor-transient-t-end', '0.1',
        '--mor-transient-waveform', 'single_step',
        '--mor-transient-initial-mode', 'ambient',
        '--output-dir', $Output, '--fast-run'
    )
}

$cases = @(
    [pscustomobject]@{
        Name = 'two-cube'
        Config = 'configs\two_cube_parametric_h.txt'
        SourceMode = 'generalized-dynamic'
        Interfaces = '0'
        CandidateDimension = 32
        MaximumIterations = 8
    },
    [pscustomobject]@{
        Name = 'ten-cube'
        Config = 'configs\ten_cube_parametric_h.txt'
        SourceMode = 'trace-only'
        Interfaces = '0,1,2,3,4,5,6,7,8'
        CandidateDimension = 16
        MaximumIterations = 3
    }
)

$summary = @()
$eigenpairs = @()
Push-Location $project
try {
    foreach ($case in $cases) {
        $modelRun = Join-Path $result "$($case.Name)_model"
        $exactRun = Join-Path $result "$($case.Name)_exact"
        $modelPath = Join-Path $modelRun 'model'
        New-Item -ItemType Directory -Path $modelRun | Out-Null
        $generate = Common-HybridArguments $case.Config $modelRun
        $generate += @(
            '--optimal-port-source-mode', $case.SourceMode,
            '--mor-transient-save', $modelPath)
        & $exe @generate
        if ($LASTEXITCODE -ne 0) {
            throw "$($case.Name) hybrid model generation failed."
        }

        New-Item -ItemType Directory -Path $exactRun | Out-Null
        $exact = Common-HybridArguments $case.Config $exactRun
        $exact += @(
            '--optimal-port-source-mode', $case.SourceMode,
            '--global-interface-coarse-prototype',
            '--global-interface-coarse-inverse-mode', 'exact-pcg',
            '--global-interface-coarse-explicit-reference',
            '--global-interface-coarse-rank', '4',
            '--global-interface-coarse-candidate-dimension',
                "$($case.CandidateDimension)",
            '--global-interface-coarse-max-iters',
                "$($case.MaximumIterations)",
            '--global-interface-coarse-krylov-sweeps', '2',
            '--global-interface-coarse-tol', '1e-8',
            '--global-interface-coarse-inner-max-iters', '1000',
            '--global-interface-coarse-inner-tol', '1e-10',
            '--global-interface-coarse-interface-ids',
                $case.Interfaces,
            '--mor-transient-load', $modelPath)
        & $exe @exact
        if ($LASTEXITCODE -ne 0) {
            throw "$($case.Name) exact projected coarse run failed."
        }
        $row = Import-Csv -LiteralPath (Join-Path $exactRun `
            'milestone8_global_coarse_prototype.csv')
        $row | Add-Member -NotePropertyName validation_case `
            -NotePropertyValue $case.Name
        $summary += $row
        foreach ($mode in @(
            Import-Csv -LiteralPath (Join-Path $exactRun `
                'milestone8_global_coarse_modes.csv'))) {
            $mode | Add-Member -NotePropertyName validation_case `
                -NotePropertyValue $case.Name
            $eigenpairs += $mode
        }
    }
} finally {
    Pop-Location
}

$summary | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $result 'milestone8_exact_coarse_small_cases.csv')
$eigenpairs | Export-Csv -NoTypeInformation -Encoding UTF8 (
    Join-Path $result 'milestone8_exact_coarse_eigenpairs.csv')
