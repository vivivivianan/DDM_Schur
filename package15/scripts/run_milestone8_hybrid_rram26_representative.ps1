param(
    [string]$Exe = ".\build\Release\SIPGHeatDDM3D.exe",
    [string]$Config =
        "D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt",
    [string]$Audit =
        ".\outputs\milestone8_large_case_topology_audit.csv",
    [string]$RawRoot =
        ".\outputs\milestone8_hybrid_rram26_runs",
    [int]$Threads = 8,
    [switch]$AggregateOnly
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $project
try {
    $exePath = (Resolve-Path $Exe).Path
    $configPath = (Resolve-Path $Config).Path
    $auditPath = (Resolve-Path $Audit).Path
    $rawPath = [System.IO.Path]::GetFullPath($RawRoot)
    New-Item -ItemType Directory -Force -Path $rawPath | Out-Null
    $env:OMP_NUM_THREADS = "$Threads"
    $env:MKL_NUM_THREADS = "$Threads"
    $all = @()
    foreach ($rank in 8,16) {
        $output = Join-Path $rawPath "rank_$rank"
        if (-not $AggregateOnly) {
            & $exePath --transient --config $configPath `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method hybrid-randomized `
                --mor-arnoldi-moments 1 `
                --randomized-port-rank $rank `
                --randomized-oversampling 5 `
                --randomized-power-iterations 1 `
                --randomized-seed 12345 `
                --randomized-port-representative-pilot `
                --residual-krylov-max-rank 4 `
                --residual-krylov-max-sweeps 2 `
                --residual-krylov-tol 1e-4 `
                --residual-krylov-block-size 4 `
                --residual-krylov-probe-mode operator-geometry `
                --optimal-port-source-mode trace-only `
                --optimal-port-inner-solver woodbury-exact `
                --optimal-port-inner-tol 1e-10 `
                --optimal-port-inner-refinement-max-iters 3 `
                --optimal-port-inner-refinement-tol 1e-10 `
                --optimal-port-topology-audit-csv $auditPath `
                --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
                --mor-transient-waveform single_step `
                --mor-transient-initial-mode ambient `
                --output-dir $output --fast-run
            if ($LASTEXITCODE -ne 0) {
                throw "RRAM26 hybrid rank $rank process failed."
            }
        }
        $rows = @(Import-Csv (
            Join-Path $output `
                'milestone8_hybrid_rram26_representative.csv'))
        foreach ($row in $rows) {
            $row | Add-Member -NotePropertyName rank_budget `
                -NotePropertyValue $rank
            $all += $row
        }
        if ($rows.Count -ne 3 -or
            @($rows | Where-Object {
                $_.status -notlike 'success*'
            }).Count -gt 0) {
            break
        }
    }
    $all | Export-Csv -NoTypeInformation -Encoding utf8 `
        outputs\milestone8_hybrid_rram26_representative.csv
    $timing = @()
    if (Test-Path outputs\milestone8_hybrid_timing.csv) {
        $timing = @(Import-Csv outputs\milestone8_hybrid_timing.csv |
            Where-Object { $_.case -ne 'RRAM26' })
    }
    $memory = @()
    if (Test-Path outputs\milestone8_hybrid_memory.csv) {
        $memory = @(Import-Csv outputs\milestone8_hybrid_memory.csv |
            Where-Object { $_.case -ne 'RRAM26' })
    }
    foreach ($row in $all) {
        $timing += [pscustomobject]@{
            case = 'RRAM26'
            configuration =
                "rank$($row.rank_budget)-$($row.selection)"
            randomized_build_time_s =
                [double]$row.randomized_build_time_s
            mandatory_residual_build_time_s =
                [double]$row.residual_build_time_s
            total_basis_build_time_s =
                [double]$row.total_basis_build_time_s
            online_interface_s = 0.0
            recovery_s = 0.0
        }
        $memory += [pscustomobject]@{
            case = 'RRAM26'
            configuration =
                "rank$($row.rank_budget)-$($row.selection)"
            interface_id = [int]$row.interface_id
            component = 'hybrid-total'
            bytes = [uint64]$row.incremental_memory_bytes
        }
    }
    $timing | Export-Csv -NoTypeInformation -Encoding utf8 `
        outputs\milestone8_hybrid_timing.csv
    $memory | Export-Csv -NoTypeInformation -Encoding utf8 `
        outputs\milestone8_hybrid_memory.csv
} finally {
    Pop-Location
}
