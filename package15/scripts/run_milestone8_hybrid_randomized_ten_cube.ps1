param(
    [string]$Exe = ".\build\Release\SIPGHeatDDM3D.exe",
    [string]$RawRoot = ".\outputs\milestone8_hybrid_runs",
    [int]$Threads = 8
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $project
try {
    $exePath = (Resolve-Path $Exe).Path
    $rawPath = [System.IO.Path]::GetFullPath($RawRoot)
    New-Item -ItemType Directory -Force -Path $rawPath | Out-Null
    $env:OMP_NUM_THREADS = "$Threads"
    $env:MKL_NUM_THREADS = "$Threads"

    $cases = @(
        [pscustomobject]@{
            label = 'mandatory-only'; method = 'mandatory-only'
            randomizedRank = 8; enrichmentRank = 0
        },
        [pscustomobject]@{
            label = 'randomized-transfer-only'
            method = 'randomized-transfer'
            randomizedRank = 8; enrichmentRank = 0
        },
        [pscustomobject]@{
            label = 'mandatory+randomized'
            method = 'hybrid-randomized'
            randomizedRank = 8; enrichmentRank = 0
        },
        [pscustomobject]@{
            label = 'mandatory+randomized+residual'
            method = 'hybrid-randomized'
            randomizedRank = 8; enrichmentRank = 4
        }
    )
    $comparison = @()
    $timing = @()
    $memory = @()
    foreach ($case in $cases) {
        $safe = $case.label.Replace('+', '_').Replace('-', '_')
        $output = Join-Path $rawPath "ten_cube_$safe"
        & $exePath --transient `
            --config configs\ten_cube_parametric_h.txt `
            --mor-transient-generate `
            --mor-transient-method local-port-block-arnoldi `
            --port-basis-method $case.method `
            --mor-arnoldi-moments 1 `
            --randomized-port-rank $case.randomizedRank `
            --randomized-oversampling 5 `
            --randomized-power-iterations 1 `
            --randomized-seed 12345 `
            --residual-krylov-max-rank $case.enrichmentRank `
            --residual-krylov-max-sweeps 2 `
            --residual-krylov-tol 1e-4 `
            --residual-krylov-block-size 4 `
            --residual-krylov-probe-mode operator-geometry `
            --optimal-port-source-mode trace-only `
            --optimal-port-inner-solver woodbury-exact `
            --optimal-port-inner-tol 1e-10 `
            --optimal-port-inner-refinement-max-iters 3 `
            --optimal-port-inner-refinement-tol 1e-10 `
            --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
            --mor-transient-waveform single_step `
            --mor-transient-initial-mode ambient `
            --output-dir $output --fast-run
        if ($LASTEXITCODE -ne 0) {
            throw "ten-cube hybrid case failed: $($case.label)"
        }

        $summary = Import-Csv (
            Join-Path $output 'local_dynamic_schur_summary.csv')
        $rankRows = @(Import-Csv (
            Join-Path $output 'local_port_rank_by_interface.csv'))
        $residual = @()
        $randomized = @()
        if (Test-Path (
                Join-Path $output `
                    'residual_krylov_interface_diagnostics.csv')) {
            $residual = @(Import-Csv (
                Join-Path $output `
                    'residual_krylov_interface_diagnostics.csv'))
        }
        if (Test-Path (
                Join-Path $output `
                    'randomized_transfer_interface_diagnostics.csv')) {
            $randomized = @(Import-Csv (
                Join-Path $output `
                    'randomized_transfer_interface_diagnostics.csv'))
        }
        $portRank = ($rankRows |
            Measure-Object total_port_rank -Sum).Sum
        $mandatoryRank = ($residual |
            Measure-Object mandatory_rank_total -Sum).Sum
        $acceptedRandomized = if (
            $case.method -eq 'randomized-transfer') {
            ($randomized | Measure-Object accepted_rank -Sum).Sum
        } else {
            ($residual |
                Measure-Object accepted_randomized_rank -Sum).Sum
        }
        $acceptedResidual = ($residual |
            Measure-Object accepted_enrichment_rank -Sum).Sum
        $randomizedSolves = ($randomized |
            Measure-Object target_solve_count -Sum).Sum
        $residualSolves = ($residual |
            Measure-Object target_solve_count -Sum).Sum
        $peakIncremental = [uint64]0
        foreach ($row in $randomized) {
            $peakIncremental = [Math]::Max(
                $peakIncremental,
                [uint64]$row.peak_incremental_memory_bytes)
            $memory += [pscustomobject]@{
                case = 'ten-cube'; configuration = $case.label
                interface_id = [int]$row.interface_id
                component = 'randomized-transfer'
                bytes = [uint64]$row.peak_incremental_memory_bytes
            }
        }
        foreach ($row in $residual) {
            $peakIncremental = [Math]::Max(
                $peakIncremental,
                [uint64]$row.peak_incremental_memory_bytes)
            $memory += [pscustomobject]@{
                case = 'ten-cube'; configuration = $case.label
                interface_id = [int]$row.interface_id
                component = 'mandatory-residual'
                bytes = [uint64]$row.peak_incremental_memory_bytes
            }
        }
        $comparison += [pscustomobject]@{
            case = 'ten-cube'; configuration = $case.label
            method = $case.method
            requested_randomized_rank = $case.randomizedRank
            requested_residual_rank = $case.enrichmentRank
            mandatory_rank = [int]$mandatoryRank
            accepted_randomized_rank = [int]$acceptedRandomized
            accepted_residual_rank = [int]$acceptedResidual
            port_rank = [int]$portRank
            temperature_relative_L2 =
                [double]$summary.space_time_relative_l2
            max_node_error_k =
                [double]$summary.maximum_temperature_error_k
            flux_relative_L2 =
                [double]$summary.maximum_fom_rom_flux_relative_l2
            interface_residual =
                [double]$summary.maximum_interface_relative_residual
            full_residual = [double]$summary.maximum_full_residual
            basis_build_time_s = [double]$summary.port_basis_seconds
            target_solve_count =
                [int]$randomizedSolves + [int]$residualSolves
            peak_incremental_memory_bytes = $peakIncremental
            process_peak_memory_bytes =
                [uint64]$summary.peak_working_set_bytes
            corrected = $false
            snapshot_used = $false
            fom_used_for_basis = $false
            status = $summary.status
        }
        $timing += [pscustomobject]@{
            case = 'ten-cube'; configuration = $case.label
            randomized_build_time_s = if ($randomized.Count -gt 0) {
                [double](Import-Csv (
                    Join-Path $output `
                        'randomized_transfer_timing.csv')).total_basis_time_s
            } else { 0.0 }
            mandatory_residual_build_time_s =
                if ($residual.Count -gt 0) {
                    [double](Import-Csv (
                        Join-Path $output `
                            'residual_krylov_build_timing.csv')).total_basis_seconds
                } else { 0.0 }
            total_basis_build_time_s =
                [double]$summary.port_basis_seconds
            online_interface_s =
                [double]$summary.interface_solve_seconds
            recovery_s = [double]$summary.local_recovery_seconds
        }
    }
    $comparison | Export-Csv -NoTypeInformation -Encoding utf8 `
        outputs\milestone8_hybrid_port_comparison.csv
    $comparison | Export-Csv -NoTypeInformation -Encoding utf8 `
        outputs\milestone8_hybrid_rank_accuracy.csv
    $timing | Export-Csv -NoTypeInformation -Encoding utf8 `
        outputs\milestone8_hybrid_timing.csv
    $memory | Export-Csv -NoTypeInformation -Encoding utf8 `
        outputs\milestone8_hybrid_memory.csv
} finally {
    Pop-Location
}
