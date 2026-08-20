param(
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [Parameter(Mandatory = $true)]
    [string]$AdaptiveModel,
    [Parameter(Mandatory = $true)]
    [string]$RomDiagnostics,
    [Parameter(Mandatory = $true)]
    [string]$RomAccuracy,
    [string]$ResultsDirectory =
        'results\milestone8_projection_diagnosis_rram26',
    [string]$ProjectionInterfaceIds = '13,23,18,4,10',
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$configPath = (Resolve-Path -LiteralPath $Config).Path
$modelPath = (Resolve-Path -LiteralPath $AdaptiveModel).Path
$romDiagnosticsPath = (Resolve-Path -LiteralPath $RomDiagnostics).Path
$romAccuracyPath = (Resolve-Path -LiteralPath $RomAccuracy).Path
$result = if ([IO.Path]::IsPathRooted($ResultsDirectory)) {
    [IO.Path]::GetFullPath($ResultsDirectory)
} else {
    [IO.Path]::GetFullPath((Join-Path $project $ResultsDirectory))
}
if (Test-Path -LiteralPath $result) {
    throw "Refusing to overwrite M8.12 results: $result"
}
if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) {
    throw "Release executable not found: $exe"
}
$selected = @(
    $ProjectionInterfaceIds.Split(',') |
        ForEach-Object { [int]$_.Trim() })
if ($selected.Count -ne 5 -or
    @($selected | Select-Object -Unique).Count -ne 5) {
    throw 'M8.12 requires exactly five unique difficult interfaces.'
}

$env:OMP_NUM_THREADS = "$Threads"
$env:MKL_NUM_THREADS = "$Threads"
$env:MKL_DYNAMIC = 'FALSE'
New-Item -ItemType Directory -Path $result | Out-Null

Push-Location $project
try {
    & $exe --transient --config $configPath `
        --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --port-basis-method hybrid-randomized `
        --mor-arnoldi-moments 1 `
        --mor-interface-rank 0 `
        --randomized-port-rank 16 `
        --randomized-oversampling 5 `
        --randomized-power-iterations 1 `
        --randomized-seed 12345 `
        --residual-krylov-max-rank 4 `
        --residual-krylov-max-sweeps 2 `
        --residual-krylov-tol 1e-4 `
        --residual-krylov-block-size 4 `
        --residual-krylov-probe-mode operator-geometry `
        --residual-krylov-inner-solver woodbury-exact `
        --optimal-port-source-mode trace-only `
        --optimal-port-inner-solver woodbury-exact `
        --optimal-port-inner-tol 1e-10 `
        --optimal-port-inner-refinement-max-iters 3 `
        --optimal-port-inner-refinement-tol 1e-10 `
        --history-compression-method deterministic-rrqr `
        --history-compression-rank 64 `
        --history-compression-tolerance 1e-12 `
        --milestone8-adaptive-production `
        --milestone8-projection-diagnosis `
        --projection-interface-ids $ProjectionInterfaceIds `
        --mor-transient-load $modelPath `
        --mor-transient-dt 0.01 `
        --mor-transient-t-end 0.01 `
        --mor-transient-waveform mixed_frequency `
        --mor-transient-initial-mode ambient `
        --output-dir $result --fast-run
    if ($LASTEXITCODE -ne 0) {
        throw 'M8.12 five-interface projection diagnosis failed.'
    }

    $temperaturePath = Join-Path $result `
        'milestone8_projection_temperature.csv'
    $fluxPath = Join-Path $result `
        'milestone8_projection_flux.csv'
    $forcingPath = Join-Path $result `
        'milestone8_projection_forcing.csv'
    $summaryPath = Join-Path $result `
        'milestone8_projection_summary.csv'
    $temperature = @(Import-Csv -LiteralPath $temperaturePath)
    $flux = @(Import-Csv -LiteralPath $fluxPath)
    $forcing = @(Import-Csv -LiteralPath $forcingPath)
    $summary = Import-Csv -LiteralPath $summaryPath |
        Select-Object -First 1
    if ($summary.status -ne 'passed' -or
        [int]$summary.requested_interfaces -ne 5 -or
        [int]$summary.completed_interfaces -ne 5 -or
        [int]$summary.transient_steps -ne 1 -or
        [int]$summary.full_field_read -ne 0 -or
        [int]$summary.snapshot_used -ne 0 -or
        [int]$summary.fom_used_for_basis -ne 0 -or
        [int]$summary.pod_used -ne 0 -or
        [int]$summary.svd_used -ne 0 -or
        $temperature.Count -ne 5 -or
        $flux.Count -ne 5 -or
        $forcing.Count -ne 5) {
        throw 'M8.12 diagnostic/provenance gate failed.'
    }
    $romRows = @(Import-Csv -LiteralPath $romDiagnosticsPath)
    $romOverall = Import-Csv -LiteralPath $romAccuracyPath |
        Where-Object method -eq 'M8.9 Adaptive Operator Port' |
        Select-Object -First 1
    if (-not $romOverall) {
        throw 'M8.9 one-step accuracy reference is missing.'
    }
    $fluxById = @{}
    $forcingById = @{}
    $romById = @{}
    foreach ($row in $flux) {
        $fluxById[[int]$row.interface_id] = $row
    }
    foreach ($row in $forcing) {
        $forcingById[[int]$row.interface_id] = $row
    }
    foreach ($row in $romRows) {
        if ([int]$row.step -eq 1 -and
            [int]$row.corrected -eq 0) {
            $romById[[int]$row.interface_id] = $row
        }
    }

    $comparison = foreach ($temperatureRow in $temperature) {
        $id = [int]$temperatureRow.interface_id
        if (-not $fluxById.ContainsKey($id) -or
            -not $forcingById.ContainsKey($id) -or
            -not $romById.ContainsKey($id)) {
            throw "M8.12 comparison row is incomplete for interface $id."
        }
        $fluxRow = $fluxById[$id]
        $forcingRow = $forcingById[$id]
        $rom = $romById[$id]
        $temperatureProjection =
            [double]$temperatureRow.temperature_projection_error
        $temperatureRom =
            [double]$rom.temperature_relative_l2
        $fluxProjection =
            [double]$fluxRow.flux_projection_error
        $fluxRom = [double]$rom.flux_relative_l2
        $forcingProjection =
            [double]$forcingRow.forcing_projection_error
        $classification = @()
        if ($temperatureProjection -ge 1e-2 -or
            $temperatureProjection -ge 0.25 * $temperatureRom) {
            $classification += 'A_port_space'
        }
        if ($fluxProjection -ge 1e-2) {
            $classification += 'C_flux_representation'
        }
        if ($temperatureProjection -lt 0.1 * $temperatureRom -and
            $forcingProjection -ge 1e-2) {
            $classification += 'B_reduced_forcing'
        } elseif ($forcingProjection -ge 1e-2) {
            $classification += 'forcing_secondary'
        }
        if ($classification.Count -eq 0) {
            $classification += 'projection_passed'
        }
        [pscustomobject]@{
            case = 'rram26'
            interface_id = $id
            interface_dofs = [int]$temperatureRow.interface_dofs
            local_port_rank = [int]$temperatureRow.local_port_rank
            energy_projection_rank =
                [int]$temperatureRow.energy_projection_rank
            temperature_projection_error = $temperatureProjection
            temperature_ROM_error = $temperatureRom
            temperature_projection_to_ROM_ratio =
                $temperatureProjection /
                [Math]::Max(1e-300, $temperatureRom)
            flux_projection_error = $fluxProjection
            flux_ROM_error = $fluxRom
            flux_projection_to_ROM_ratio =
                $fluxProjection / [Math]::Max(1e-300, $fluxRom)
            flux_triangle_count = [int]$fluxRow.flux_triangle_count
            flux_area = [double]$fluxRow.flux_area
            forcing_projection_error = $forcingProjection
            input_projection_error =
                [double]$forcingRow.input_projection_error
            boundary_projection_error =
                [double]$forcingRow.boundary_projection_error
            history_projection_error =
                [double]$forcingRow.history_projection_error
            interface_residual_ROM =
                [double]$romOverall.interface_residual
            interface_residual_ROM_scope = 'global_one_step'
            full_residual_ROM = [double]$rom.local_full_residual
            global_full_residual_ROM =
                [double]$romOverall.full_residual
            temperature_jump_rms_ROM =
                [double]$rom.temperature_jump_rms_k
            max_flux_imbalance_ROM =
                [double]$rom.max_relative_flux_imbalance
            comparison_norm_note =
                'projection=S-energy; ROM=local-volume-relative-L2'
            classification = $classification -join '+'
            corrected = 0
            full_field_read = 0
            snapshot_used = 0
            fom_used_for_basis = 0
            status = $temperatureRow.status
        }
    }
    $comparison = @($comparison | Sort-Object interface_id)
    $comparison | Export-Csv -NoTypeInformation -Encoding UTF8 (
        Join-Path $result 'milestone8_projection_vs_rom_error.csv')
    $comparison | Export-Csv -NoTypeInformation -Encoding UTF8 (
        Join-Path $result `
            'milestone8_projection_difficult_interfaces.csv')

    foreach ($name in @(
        'milestone8_projection_temperature.csv',
        'milestone8_projection_flux.csv',
        'milestone8_projection_forcing.csv',
        'milestone8_projection_vs_rom_error.csv',
        'milestone8_projection_difficult_interfaces.csv')) {
        Copy-Item -LiteralPath (Join-Path $result $name) `
            -Destination (Join-Path $project "outputs\$name")
    }
} finally {
    Pop-Location
}
