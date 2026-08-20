param(
    [switch]$SkipBuild,
    [switch]$SkipTests,
    [switch]$AggregateOnly
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$runRoot = Join-Path $project 'build\local_block_arnoldi_milestone4'
$outputs = Join-Path $project 'outputs'
New-Item -ItemType Directory -Force -Path $runRoot, $outputs | Out-Null

function Invoke-Checked([string[]]$Arguments) {
    & $exe @Arguments
    if ($LASTEXITCODE -ne 0) { throw "SIPGHeatDDM3D failed with exit code $LASTEXITCODE" }
}

function Invoke-LocalCase([string]$Name, [string]$Waveform, [string]$InitialMode) {
    Invoke-Checked @(
        '--transient', '--config', 'configs\two_cube_parametric_h.txt',
        '--mor-transient-generate', '--mor-transient-method', 'local-block-arnoldi',
        '--mor-arnoldi-moments', '3', '--mor-arnoldi-rank-tolerance', '1e-10',
        '--mor-transient-dt', '0.1', '--mor-transient-t-end', '1',
        '--mor-transient-waveform', $Waveform,
        '--mor-transient-initial-mode', $InitialMode,
        '--output-dir', (Join-Path $runRoot $Name), '--fast-run')
}

Push-Location $project
try {
    if (-not $AggregateOnly) {
        if (-not $SkipBuild) {
            cmake --build build --config Release --parallel 4
            if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }
        }
        if (-not $SkipTests) {
            $ctestOutput = & ctest --test-dir build -C Release --output-on-failure 2>&1
            $ctestCode = $LASTEXITCODE
            $ctestOutput | Write-Host
            $ctestOutput | Set-Content -Encoding utf8 `
                (Join-Path $outputs 'local_block_arnoldi_milestone4_ctest_results.txt')
            if ($ctestCode -ne 0) { throw 'CTest failed.' }
        }
        foreach ($case in @(
            @('ambient_step', 'single_step', 'ambient'),
            @('uniform_step', 'single_step', 'uniform'),
            @('steady_dc', 'multi_step', 'steady'),
            @('ambient_pulse', 'rectangular_pulse', 'ambient'),
            @('ambient_multilevel', 'piecewise_multilevel', 'ambient'),
            @('ambient_duty', 'variable_duty_cycle', 'ambient'),
            @('ambient_async', 'asynchronous_hotspots', 'ambient'))) {
            Invoke-LocalCase $case[0] $case[1] $case[2]
        }
        Invoke-Checked @(
            '--transient', '--config', 'configs\two_cube_parametric_h.txt',
            '--mor-transient-generate', '--mor-transient-method', 'block-arnoldi',
            '--mor-arnoldi-moments', '2', '--mor-transient-dt', '0.1',
            '--mor-transient-t-end', '1', '--mor-transient-waveform', 'single_step',
            '--mor-transient-compare-fom', '--output-dir', (Join-Path $runRoot 'global_step'),
            '--fast-run')
    }

    $localRows = foreach ($directory in Get-ChildItem $runRoot -Directory |
        Where-Object Name -ne 'global_step') {
        $path = Join-Path $directory.FullName 'local_dynamic_schur_summary.csv'
        if (Test-Path $path) { Import-Csv $path }
    }
    $localRows | Export-Csv (Join-Path $outputs 'local_dynamic_schur_milestone4.csv') `
        -NoTypeInformation
    Copy-Item (Join-Path $runRoot 'ambient_step\local_block_arnoldi_rank.csv') `
        (Join-Path $outputs 'local_block_arnoldi_milestone4_rank.csv') -Force

    $canonical = $localRows | Where-Object {
        $_.waveform -eq 'single_step' -and $_.initial_mode -eq 'ambient' } |
        Select-Object -First 1
    $global = Import-Csv (Join-Path $runRoot 'global_step\transient_block_arnoldi_summary.csv')
    $comparison = @(
        [pscustomobject]@{
            method = 'Transient monolithic PARDISO'; rank = 7214; interface_dofs = 0
            offline_seconds = [double]$canonical.fom_factor_seconds
            online_seconds = [double]$canonical.fom_solve_seconds
            space_time_relative_l2 = 0.0; maximum_absolute_k = 0.0
        },
        [pscustomobject]@{
            method = 'Global Block Arnoldi'; rank = [int]$global.rank; interface_dofs = 0
            offline_seconds = [double]$global.offline_seconds
            online_seconds = [double]$canonical.fom_solve_seconds / [double]$global.speedup
            space_time_relative_l2 = [double]$global.space_time_relative_l2
            maximum_absolute_k = [double]$global.maximum_absolute_k
        },
        [pscustomobject]@{
            method = 'Local Block Arnoldi + Dynamic Schur'
            rank = [int]$canonical.total_local_rank
            interface_dofs = [int]$canonical.full_interface_dofs
            offline_seconds = [double]$canonical.local_basis_setup_seconds +
                [double]$canonical.dynamic_schur_setup_seconds
            online_seconds = [double]$canonical.local_online_core_seconds
            space_time_relative_l2 = [double]$canonical.space_time_relative_l2
            maximum_absolute_k = [double]$canonical.maximum_absolute_k
        })
    $comparison | Export-Csv (Join-Path $outputs 'local_block_arnoldi_milestone4_summary.csv') `
        -NoTypeInformation

    $worstL2 = ($localRows | Measure-Object -Property space_time_relative_l2 -Maximum).Maximum
    $worstMax = ($localRows | Measure-Object -Property maximum_absolute_k -Maximum).Maximum
    $steady = $localRows | Where-Object initial_mode -eq 'steady' | Select-Object -First 1
    $report = @"
# Milestone 4: two-cube Local Block Arnoldi + Dynamic Schur

## Result

The transient mainline now uses two independently generated local Block Arnoldi
bases and retains all $($canonical.full_interface_dofs) physical interface DOFs.
No global basis is sliced or relabelled as local.  The existing Global Block
Arnoldi implementation is retained only as a benchmark.

- Subdomains: $($canonical.subdomains)
- Total local rank: $($canonical.total_local_rank)
- Worst space-time relative L2 across the required initial/waveform suite: $worstL2
- Worst maximum nodal error: $worstMax K
- Steady-initial DC consistency: $($steady.space_time_relative_l2)
- Fixed-dt Dynamic Schur symbolic/numerical calls: $($canonical.dynamic_schur_symbolic_calls) / $($canonical.dynamic_schur_numerical_calls)

## Algorithm audit

Each subdomain uses its own interior matrices and PARDISO factor.  Its initial
Block Arnoldi block contains local physical power, `-K_IGamma E_Gamma`, and
`-C_IGamma E_Gamma`.  `E_Gamma` is compressed by two-pass rank-revealing QR from constant, steady FOM
interface traces, and first dynamic-moment traces.  It is used only offline;
the online interface temperature remains full order.

For fixed `dt=0.1 s`, each reduced local `C_II/dt+K_II` factor and the explicit
two-cube Dynamic Schur PARDISO factor are constructed once.  Every time step
performs only condensed-RHS assembly, interface solve, and local recovery.

## Validation coverage

Ambient, uniform nonzero, and steady-state initial conditions were tested with
single step, pulse, piecewise multi-level, variable duty-cycle, and asynchronous
source waveforms.  The reduced residual is reported separately from the sampled
full residual.  Interface temperature jump and physical flux imbalance are
reported at every step in the raw per-case CSV; physical-gradient flux imbalance
is a diagnostic of the underlying SIPG field and is not hidden as a ROM error.

## Performance interpretation

The monolithic and Global Block Arnoldi rows are included in the summary CSV.
For this small two-domain problem Global Block Arnoldi has the smallest online
state.  The Local method pays for the unreduced 766-DOF interface; its purpose is
the partition-preserving path that extends to independently reusable subdomains.
"@
    Set-Content -Path (Join-Path $outputs 'local_block_arnoldi_milestone4_report.md') `
        -Value $report -Encoding utf8
} finally {
    Pop-Location
}
