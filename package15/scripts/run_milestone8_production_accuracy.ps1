param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('rram26','chiplet')]
    [string]$Case,
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [Parameter(Mandatory = $true)]
    [string]$M7Model,
    [Parameter(Mandatory = $true)]
    [string]$M8Model,
    [string]$ResultsDirectory = '',
    [int]$Threads = 8,
    [switch]$AggregateOnly
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$configPath = (Resolve-Path -LiteralPath $Config).Path
$m7ModelPath = (Resolve-Path -LiteralPath $M7Model).Path
$m8ModelPath = (Resolve-Path -LiteralPath $M8Model).Path
if ([string]::IsNullOrWhiteSpace($ResultsDirectory)) {
    $ResultsDirectory = "results\milestone8_accuracy_$Case"
}
$root = [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
if ((Test-Path -LiteralPath $root) -and -not $AggregateOnly) {
    throw "Refusing to overwrite production accuracy directory: $root"
}
if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
    throw "Release executable not found: $exe"
}

$dt = if ($Case -eq 'rram26') { '0.01' } else { '0.1' }
$moments = if ($Case -eq 'rram26') { '1' } else { '4' }
$waveform = if ($Case -eq 'rram26') {
    'mixed_frequency'
} else {
    'asynchronous_hotspots'
}
$sourceMode = if ($Case -eq 'rram26') {
    'trace-only'
} else {
    'generalized-dynamic'
}
$m7PortRank = if ($Case -eq 'rram26') { '128' } else { '200' }
$m7TemperatureWeight = if ($Case -eq 'rram26') { '100' } else { '1' }

function Invoke-Checked([string]$Label, [string[]]$Arguments) {
    $directory = Join-Path $root $Label
    New-Item -ItemType Directory -Path $directory | Out-Null
    & $exe @Arguments --output-dir $directory --fast-run | Out-Host
    if ($LASTEXITCODE -ne 0) {
        throw "M8 production accuracy run failed: $Case/$Label"
    }
    return $directory
}

function Read-Summary([string]$Directory) {
    $path = Join-Path $Directory 'local_dynamic_schur_summary.csv'
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Missing accuracy summary: $path"
    }
    return Import-Csv -LiteralPath $path | Select-Object -First 1
}

$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
if (-not (Test-Path -LiteralPath $root)) {
    New-Item -ItemType Directory -Path $root | Out-Null
}

$common = @(
    '--transient', '--config', $configPath,
    '--mor-transient-generate',
    '--mor-transient-method', 'local-port-block-arnoldi',
    '--mor-arnoldi-moments', $moments,
    '--mor-interface-rank', '0',
    '--mor-transient-dt', $dt,
    '--mor-transient-t-end', $dt,
    '--mor-transient-waveform', $waveform,
    '--mor-transient-initial-mode', 'ambient')

Push-Location $project
try {
    if ($AggregateOnly) {
        $fullDirectory = Join-Path $root 'full_interface'
        $m7Directory = Join-Path $root 'm7_port_pod'
        $m8Directory = Join-Path $root 'm8_hybrid_operator_port'
    } else {
        $fullDirectory = Invoke-Checked 'full_interface' (
            $common + @(
                '--port-basis-method', 'full-interface',
                '--local-mor-matrix-free-threshold', '20000',
                '--max-pcg-iterations', '500',
                '--gmres-restart', '100',
                '--pcg-tolerance', '1e-5',
                '--schur-proxy-ring', '1',
                '--schur-proxy-block-size', '64'))

        $m7Directory = Invoke-Checked 'm7_port_pod' (
            $common + @(
                '--port-basis-method', 'port-pod',
                '--local-port-rank', $m7PortRank,
                '--local-port-energy-tolerance', '1e-20',
                '--local-port-temperature-weight', $m7TemperatureWeight,
                '--local-port-flux-weight', '1',
                '--local-port-residual-weight', '1',
                '--local-port-enrichment-rounds', '5',
                '--mor-transient-load', $m7ModelPath))

        $m8Directory = Invoke-Checked 'm8_hybrid_operator_port' (
            $common + @(
                '--port-basis-method', 'hybrid-randomized',
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
                '--optimal-port-source-mode', $sourceMode,
                '--optimal-port-inner-solver', 'woodbury-exact',
                '--optimal-port-inner-tol', '1e-10',
                '--optimal-port-inner-refinement-max-iters', '3',
                '--optimal-port-inner-refinement-tol', '1e-10',
                '--history-compression-method', 'deterministic-rrqr',
                '--history-compression-rank', '64',
                '--history-compression-tolerance', '1e-12',
                '--mor-transient-load', $m8ModelPath))
    }

    $runs = @(
        [pscustomobject]@{
            method = 'Full Interface'
            directory = $fullDirectory
            snapshot_used = 0
            fom_used_for_basis = 0
        },
        [pscustomobject]@{
            method = 'M7 Port POD'
            directory = $m7Directory
            snapshot_used = 1
            fom_used_for_basis = 1
        },
        [pscustomobject]@{
            method = 'M8 Hybrid Operator Port'
            directory = $m8Directory
            snapshot_used = 0
            fom_used_for_basis = 0
        })

    $rows = foreach ($run in $runs) {
        $summary = Read-Summary $run.directory
        if ($run.method -eq 'M8 Hybrid Operator Port' -and (
            [int]$summary.port_snapshot_used -ne 0 -or
            [int]$summary.port_fom_used_for_basis -ne 0)) {
            throw 'M8 basis provenance gate failed during production accuracy.'
        }
        [pscustomobject]@{
            case = $Case
            method = $run.method
            scope = 'one_step'
            port_dimension = $summary.port_dimension
            total_local_rank = $summary.total_local_rank
            temperature_relative_l2 = $summary.space_time_relative_l2
            max_error_k = $summary.maximum_absolute_k
            max_temperature_error_k =
                $summary.maximum_temperature_error_k
            flux_relative_l2 =
                $summary.maximum_fom_rom_flux_relative_l2
            interface_residual =
                $summary.maximum_interface_relative_residual
            full_residual = $summary.maximum_full_residual
            online_time_s = $summary.local_online_core_seconds
            total_time_s = $summary.total_seconds
            peak_memory_bytes = $summary.peak_working_set_bytes
            snapshot_used = $run.snapshot_used
            fom_used_for_basis = $run.fom_used_for_basis
            corrected = 0
            status = $summary.status
        }
    }
    $rows | Export-Csv -NoTypeInformation -Encoding UTF8 (
        Join-Path $root 'milestone8_accuracy.csv')
} finally {
    Pop-Location
}
