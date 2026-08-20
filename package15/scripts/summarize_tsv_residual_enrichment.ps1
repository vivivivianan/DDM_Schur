[CmdletBinding()]
param(
    [string]$Root = 'E:\tsv_pdn4_ddm32_basis_enrichment',
    [string]$Output = ''
)

$ErrorActionPreference = 'Stop'
if ([string]::IsNullOrWhiteSpace($Output)) { $Output = Join-Path $Root 'summary.csv' }

function Get-TemperatureStatistics([string]$Path) {
    $minimum = [double]::PositiveInfinity
    $maximum = [double]::NegativeInfinity
    $belowInitial = 0
    $below292 = 0
    Get-Content -LiteralPath $Path -ReadCount 10000 | ForEach-Object {
        foreach ($line in $_) {
            if ($line.StartsWith('global_dof,')) { continue }
            $temperature = [double]($line.Split(',')[5])
            if ($temperature -lt $minimum) { $minimum = $temperature }
            if ($temperature -gt $maximum) { $maximum = $temperature }
            if ($temperature -lt 293.15) { ++$belowInitial }
            if ($temperature -lt 292.0) { ++$below292 }
        }
    }
    [pscustomobject]@{ Tmin = $minimum; Tmax = $maximum; BelowInitial = $belowInitial; Below292 = $below292 }
}

function Get-LastCsvRow([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    @(Import-Csv -LiteralPath $Path | Select-Object -Last 1)[0]
}

$cases = @(
    [pscustomobject]@{ Method = 'M3 baseline'; TopN = 0; Iteration = 0; Directory = 'E:\\tsv_pdn4_ddm32_augmented_direct_cn\\m3_operational\\warm_200_full_field' },
    [pscustomobject]@{ Method = 'Residual Top2'; TopN = 2; Iteration = 1; Directory = 'residual_top2' },
    [pscustomobject]@{ Method = 'Residual Top4'; TopN = 4; Iteration = 1; Directory = 'residual_top4' },
    [pscustomobject]@{ Method = 'Residual Top8'; TopN = 8; Iteration = 1; Directory = 'residual_top8' },
    [pscustomobject]@{ Method = 'Residual Top4 two rounds'; TopN = 4; Iteration = 2; Directory = 'residual_top4_round2' },
    [pscustomobject]@{ Method = 'Residual Top4 FOM compare'; TopN = 4; Iteration = 1; Directory = 'residual_top4_fom_compare' }
)

$rows = foreach ($case in $cases) {
    $directory = if ([IO.Path]::IsPathRooted($case.Directory)) { $case.Directory } else { Join-Path $Root $case.Directory }
    $summary = Get-LastCsvRow (Join-Path $directory 'local_dynamic_schur_summary.csv')
    $accuracy = Get-LastCsvRow (Join-Path $directory 'local_dynamic_schur_accuracy_by_time.csv')
    $residual = Join-Path $directory 'rom_residual_enrichment.csv'
    $selected = if (Test-Path -LiteralPath $residual) {
        ((Import-Csv -LiteralPath $residual | Where-Object { $_.selected -eq '1' } |
            ForEach-Object { "SD$([int]$_.subdomain_id + 1)" }) -join ';')
    } else { '' }
    $field = Join-Path $directory 'local_dynamic_schur_final_temperature.csv'
    $stats = if (Test-Path -LiteralPath $field) { Get-TemperatureStatistics $field } else { $null }
    [pscustomobject]@{
        method = $case.Method
        top_n = $case.TopN
        iteration = $case.Iteration
        selected_subdomains = $selected
        rank_before = if ($summary) { [int]$summary.total_local_rank - [int]$summary.enrichment_added_rank } else { $null }
        rank_after = if ($summary) { [int]$summary.total_local_rank } else { $null }
        added_rank = if ($summary) { [int]$summary.enrichment_added_rank } else { $null }
        basis_update_time = if ($summary) { [double]$summary.enrichment_orthogonalization_seconds } else { $null }
        # The core report records the whole ROM-only enrichment phase.  The
        # remainder after local factors, correction solves, and mass
        # orthogonalization is pilot-ROM/residual evaluation plus projection.
        residual_time = if ($summary -and [double]$summary.enrichment_total_seconds -gt 0.0) {
            [double]$summary.enrichment_total_seconds - [double]$summary.enrichment_factorization_seconds -
                [double]$summary.enrichment_solve_seconds - [double]$summary.enrichment_orthogonalization_seconds
        } else { 0.0 }
        factorization_time = if ($summary) { [double]$summary.enrichment_factorization_seconds } else { $null }
        solve_time = if ($summary) { [double]$summary.time_stepping_seconds } else { $null }
        total_runtime = if ($summary) { [double]$summary.total_seconds } else { $null }
        Tmin_200ns = if ($stats) { $stats.Tmin } else { $null }
        Tmax_200ns = if ($stats) { $stats.Tmax } else { $null }
        L2_error = if ($accuracy -and $accuracy.relative_l2 -ne 'nan') { [double]$accuracy.relative_l2 } else { $null }
        max_point_error = if ($accuracy -and $accuracy.maximum_absolute_k -ne 'nan') { [double]$accuracy.maximum_absolute_k } else { $null }
        undershoot_nodes = if ($stats) { $stats.Below292 } else { $null }
        nodes_below_293K = if ($stats) { $stats.BelowInitial } else { $null }
        status = if ($summary) { $summary.status } else { 'not_run' }
    }
}
$rows | Export-Csv -LiteralPath $Output -NoTypeInformation
Write-Host "Wrote $Output"
