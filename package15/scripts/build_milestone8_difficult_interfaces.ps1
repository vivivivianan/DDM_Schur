param(
    [string]$FullRunDirectory =
        'results\milestone8_efficiency_rram26',
    [string]$BasisDirectory =
        'results\milestone8_production_basis_rram26',
    [string]$Output =
        'outputs\milestone8_difficult_interfaces.csv'
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$fullRun = [IO.Path]::GetFullPath((Join-Path $project $FullRunDirectory))
$basis = [IO.Path]::GetFullPath((Join-Path $project $BasisDirectory))
$outputPath = [IO.Path]::GetFullPath((Join-Path $project $Output))
$fluxPath = Join-Path $fullRun 'local_dynamic_schur_interface_flux.csv'
if (-not (Test-Path -LiteralPath $fluxPath -PathType Leaf)) {
    throw "Missing interface-only flux diagnostics: $fluxPath"
}

function Read-ByInterface([string]$Path) {
    $rows = @{}
    foreach ($row in Import-Csv -LiteralPath $Path) {
        $rows[[int]$row.interface_id] = $row
    }
    return $rows
}

$rankByInterface = Read-ByInterface (
    Join-Path $fullRun 'local_port_rank_by_interface.csv')
$residualByInterface = Read-ByInterface (
    Join-Path $basis 'residual_krylov_interface_diagnostics.csv')
$randomizedByInterface = Read-ByInterface (
    Join-Path $basis 'randomized_transfer_interface_diagnostics.csv')
$interfaceByPair = @{}
foreach ($rank in $rankByInterface.Values) {
    $left = [int]$rank.left_subdomain
    $right = [int]$rank.right_subdomain
    $interfaceByPair["$left`:$right"] = [int]$rank.interface_id
    $interfaceByPair["$right`:$left"] = [int]$rank.interface_id
}
$summary = Import-Csv -LiteralPath (
    Join-Path $fullRun 'local_dynamic_schur_summary.csv') |
    Select-Object -First 1

$enumerator = [IO.File]::ReadLines($fluxPath).GetEnumerator()
if (-not $enumerator.MoveNext()) {
    throw 'Interface flux diagnostics are empty.'
}
$header = $enumerator.Current.Split(',')
$column = @{}
for ($index = 0; $index -lt $header.Count; ++$index) {
    $column[$header[$index]] = $index
}
$required = @(
    'step', 'left_subdomain', 'right_subdomain', 'area_m2',
    'rom_temperature_jump_rms_k', 'rom_flux_imbalance_l2_w_m2',
    'rom_relative_flux_imbalance', 'rom_sipg_numerical_flux_w_m2')
foreach ($name in $required) {
    if (-not $column.ContainsKey($name)) {
        throw "Missing interface diagnostic column: $name"
    }
}

$culture = [Globalization.CultureInfo]::InvariantCulture
$styles = [Globalization.NumberStyles]::Float
$aggregates = @{}
$diagnosticStep = -1
while ($enumerator.MoveNext()) {
    $values = $enumerator.Current.Split(',')
    $leftSubdomain = [int]::Parse(
        $values[$column['left_subdomain']], $culture)
    $rightSubdomain = [int]::Parse(
        $values[$column['right_subdomain']], $culture)
    $pair = "$leftSubdomain`:$rightSubdomain"
    if (-not $interfaceByPair.ContainsKey($pair)) {
        throw "Interface flux row has no physical-port mapping: $pair"
    }
    $interfaceId = $interfaceByPair[$pair]
    $step = [int]::Parse($values[$column['step']], $culture)
    $diagnosticStep = [Math]::Max($diagnosticStep, $step)
    if (-not $aggregates.ContainsKey($interfaceId)) {
        $aggregates[$interfaceId] = [pscustomobject]@{
            area = 0.0
            temperatureJumpEnergy = 0.0
            fluxImbalanceEnergy = 0.0
            fluxScaleEnergy = 0.0
            relativeFluxEnergy = 0.0
            maximumRelativeFluxImbalance = 0.0
            triangles = 0
        }
    }
    $aggregate = $aggregates[$interfaceId]
    $area = [double]::Parse(
        $values[$column['area_m2']], $styles, $culture)
    $jump = [double]::Parse(
        $values[$column['rom_temperature_jump_rms_k']],
        $styles, $culture)
    $imbalance = [double]::Parse(
        $values[$column['rom_flux_imbalance_l2_w_m2']],
        $styles, $culture)
    $relativeImbalance = [double]::Parse(
        $values[$column['rom_relative_flux_imbalance']],
        $styles, $culture)
    $flux = [double]::Parse(
        $values[$column['rom_sipg_numerical_flux_w_m2']],
        $styles, $culture)
    $aggregate.area += $area
    $aggregate.temperatureJumpEnergy += $area * $jump * $jump
    $aggregate.fluxImbalanceEnergy +=
        $area * $imbalance * $imbalance
    $aggregate.fluxScaleEnergy += $area * $flux * $flux
    $aggregate.relativeFluxEnergy +=
        $area * $relativeImbalance * $relativeImbalance
    $aggregate.maximumRelativeFluxImbalance = [Math]::Max(
        $aggregate.maximumRelativeFluxImbalance,
        [Math]::Abs($relativeImbalance))
    ++$aggregate.triangles
}
$enumerator.Dispose()

$totalTemperatureEnergy = (
    $aggregates.Values |
    Measure-Object temperatureJumpEnergy -Sum).Sum
$totalResidualEnergy = (
    $aggregates.Values |
    Measure-Object fluxImbalanceEnergy -Sum).Sum
$maximumFluxImbalance = (
    $aggregates.Values |
    Measure-Object maximumRelativeFluxImbalance -Maximum).Maximum

$allRows = foreach ($interfaceId in ($aggregates.Keys | Sort-Object)) {
    $aggregate = $aggregates[$interfaceId]
    $rank = $rankByInterface[$interfaceId]
    $residual = $residualByInterface[$interfaceId]
    $randomized = $randomizedByInterface[$interfaceId]
    if ($null -eq $rank -or $null -eq $residual -or
        $null -eq $randomized) {
        throw "Incomplete basis diagnostics for interface $interfaceId."
    }
    $temperatureContribution =
        $aggregate.temperatureJumpEnergy /
        [Math]::Max($totalTemperatureEnergy, 1.0e-300)
    $residualContribution =
        $aggregate.fluxImbalanceEnergy /
        [Math]::Max($totalResidualEnergy, 1.0e-300)
    $fluxScore = $aggregate.maximumRelativeFluxImbalance /
        [Math]::Max($maximumFluxImbalance, 1.0e-300)
    [pscustomobject]@{
        case = 'RRAM26'
        interface_id = $interfaceId
        left_subdomain = [int]$rank.left_subdomain
        right_subdomain = [int]$rank.right_subdomain
        target_dofs = [int]$residual.target_dofs
        source_dofs = [int]$residual.source_dofs
        diagnostic_step = $diagnosticStep
        temperature_contribution =
            $temperatureContribution
        temperature_jump_rms_k = [Math]::Sqrt(
            $aggregate.temperatureJumpEnergy /
            [Math]::Max($aggregate.area, 1.0e-300))
        interface_residual_proxy = [Math]::Sqrt(
            $aggregate.fluxImbalanceEnergy /
            [Math]::Max($aggregate.fluxScaleEnergy, 1.0e-300))
        interface_residual_contribution =
            $residualContribution
        flux_imbalance_rms = [Math]::Sqrt(
            $aggregate.relativeFluxEnergy /
            [Math]::Max($aggregate.area, 1.0e-300))
        flux_imbalance_max =
            $aggregate.maximumRelativeFluxImbalance
        total_port_rank = [int]$rank.total_port_rank
        history_rank = [int]$residual.compressed_history_rank
        requested_transfer_rank =
            [int]$randomized.requested_rank
        transfer_rank = [int]$randomized.accepted_rank
        requested_residual_rank =
            [int]$residual.requested_enrichment_rank
        residual_rank =
            [int]$residual.accepted_enrichment_rank
        history_residual =
            [double]$residual.history_compression_relative_error
        transfer_error =
            [double]$randomized.basis_error_indicator
        operator_probe_residual =
            [double]$residual.final_max_probe_residual
        difficulty_score =
            0.40 * $temperatureContribution +
            0.35 * $residualContribution +
            0.25 * $fluxScore
        diagnostic_source =
            'ROM interface trace/flux only; no FOM field or FOM error used'
        snapshot_used = 0
        fom_used_for_basis = 0
    }
}

$top = @($allRows |
    Sort-Object difficulty_score -Descending |
    Select-Object -First 5)
for ($index = 0; $index -lt $top.Count; ++$index) {
    $top[$index] | Add-Member -NotePropertyName difficulty_rank `
        -NotePropertyValue ($index + 1)
}
$top |
    Export-Csv -NoTypeInformation -Encoding UTF8 -LiteralPath $outputPath
