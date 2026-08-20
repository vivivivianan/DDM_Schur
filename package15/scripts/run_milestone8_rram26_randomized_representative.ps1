param(
    [string]$Exe = ".\build\Release\SIPGHeatDDM3D.exe",
    [string]$Config =
        "D:\CPP\TEST_CHATGPT\configs\rram26_from_sim_parameter_1_bshift_minus1.txt",
    [string]$Audit =
        ".\outputs\milestone8_large_case_topology_audit.csv",
    [string]$RawRoot =
        ".\outputs\milestone8_rram26_randomized_runs",
    [int]$Threads = 8
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
    $allRows = @()
    $failureRows = @()
    foreach ($rank in 8,16) {
        $output = Join-Path $rawPath "rank_$rank"
        & $exePath --transient --config $configPath `
            --mor-transient-generate `
            --mor-transient-method local-port-block-arnoldi `
            --port-basis-method randomized-transfer `
            --mor-arnoldi-moments 1 `
            --randomized-port-rank $rank `
            --randomized-oversampling 5 `
            --randomized-power-iterations 1 `
            --randomized-seed 12345 `
            --randomized-port-representative-pilot `
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
            $failureRows += [pscustomobject]@{
                rank = $rank; selection = ''
                interface_id = ''; reason = "process_exit_$LASTEXITCODE"
            }
            break
        }
        $rows = @(Import-Csv (
            Join-Path $output `
                'milestone8_rram26_randomized_representative.csv'))
        $allRows += $rows
        $failed = @($rows | Where-Object {
            $_.status -notlike 'success*' })
        if ($failed.Count -gt 0 -or $rows.Count -ne 3) {
            foreach ($item in $failed) {
                $failureRows += [pscustomobject]@{
                    rank = $rank; selection = $item.selection
                    interface_id = $item.interface_id
                    reason = $item.status
                }
            }
            if ($rows.Count -ne 3) {
                $failureRows += [pscustomobject]@{
                    rank = $rank; selection = ''
                    interface_id = ''
                    reason = 'representative_run_stopped_early'
                }
            }
            break
        }
    }
    $allRows | Export-Csv -NoTypeInformation -Encoding utf8 `
        'outputs\milestone8_rram26_randomized_representative.csv'
    if ($failureRows.Count -gt 0) {
        $failureRows | Export-Csv -NoTypeInformation -Encoding utf8 `
            'outputs\milestone8_randomized_failures_rram26.csv'
        throw 'RRAM26 randomized representative gate failed; stopped.'
    }
    @([pscustomobject]@{
        rank = ''; selection = ''; interface_id = ''; reason = 'none'
    }) | Export-Csv -NoTypeInformation -Encoding utf8 `
        'outputs\milestone8_randomized_failures_rram26.csv'
} finally {
    Pop-Location
}
