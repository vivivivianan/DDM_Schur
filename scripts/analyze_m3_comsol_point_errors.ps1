param(
    [Parameter(Mandatory = $true)][string]$M3RawCsv,
    [Parameter(Mandatory = $true)][string]$PointComparisonCsv,
    [Parameter(Mandatory = $true)][string]$OutputDirectory
)

$ErrorActionPreference = 'Stop'
$inv = [Globalization.CultureInfo]::InvariantCulture
function ToDouble([string]$value) {
    return [double]::Parse($value, [Globalization.NumberStyles]::Float, $inv)
}

New-Item -ItemType Directory -Force -Path $OutputDirectory | Out-Null
$hotPath = Join-Path $OutputDirectory 'large_error_points_abs_ge_10K.csv'
$statsPath = Join-Path $OutputDirectory 'error_distribution_by_subdomain.csv'
$summaryPath = Join-Path $OutputDirectory 'error_hotspot_summary.txt'

$stats = @(0..31 | ForEach-Object {
    [pscustomobject]@{
        Nodes = 0; Sum = 0.0; SumSq = 0.0
        Min = [double]::PositiveInfinity; Max = [double]::NegativeInfinity; MaxAbs = 0.0
        Ge5 = 0; Ge10 = 0; Ge20 = 0; Ge30 = 0
        XMin = [double]::PositiveInfinity; XMax = [double]::NegativeInfinity
        YMin = [double]::PositiveInfinity; YMax = [double]::NegativeInfinity
        ZMin = [double]::PositiveInfinity; ZMax = [double]::NegativeInfinity
    }
})

$rawReader = [IO.StreamReader]::new($M3RawCsv)
$comparisonReader = [IO.StreamReader]::new($PointComparisonCsv)
$hotWriter = [IO.StreamWriter]::new($hotPath, $false, [Text.UTF8Encoding]::new($false))
try {
    if ($rawReader.ReadLine() -ne 'global_dof,x_m,y_m,z_m,subdomain,temperature_k') {
        throw 'Unexpected M3 raw CSV header.'
    }
    if ($comparisonReader.ReadLine() -ne 'x_m,y_m,z_m,ddm_m3_temperature_k,comsol_temperature_k,error_ddm_minus_comsol_k') {
        throw 'Unexpected point-comparison CSV header.'
    }
    $hotWriter.WriteLine('x_m,y_m,z_m,subdomain_id_1based,ddm_m3_temperature_k,comsol_temperature_k,error_ddm_minus_comsol_k')
    $count = 0; $sum = 0.0; $sumSq = 0.0; $globalMaxAbs = -1.0
    $globalMax = ''; $ddmMin = [double]::PositiveInfinity; $ddmMinRow = ''
    $comsolMin = [double]::PositiveInfinity; $comsolMinRow = ''
    while ($true) {
        $rawLine = $rawReader.ReadLine(); $comparisonLine = $comparisonReader.ReadLine()
        if ($null -eq $rawLine -and $null -eq $comparisonLine) { break }
        if ($null -eq $rawLine -or $null -eq $comparisonLine) { throw "Input row count mismatch after $count rows." }
        $raw = $rawLine.Split(','); $point = $comparisonLine.Split(',')
        if ($raw.Count -ne 6 -or $point.Count -ne 6) { throw "Malformed row at index $count." }
        $subdomain = [int]$raw[4]
        if ($subdomain -lt 0 -or $subdomain -ge 32) { throw "Invalid subdomain at index $count." }
        $x = ToDouble $point[0]; $y = ToDouble $point[1]; $z = ToDouble $point[2]
        $ddm = ToDouble $point[3]; $comsol = ToDouble $point[4]; $deltaK = ToDouble $point[5]
        $absolute = [Math]::Abs($deltaK); $s = $stats[$subdomain]
        $s.Nodes++; $s.Sum += $deltaK; $s.SumSq += $deltaK * $deltaK
        $s.Min = [Math]::Min($s.Min, $deltaK); $s.Max = [Math]::Max($s.Max, $deltaK)
        $s.MaxAbs = [Math]::Max($s.MaxAbs, $absolute)
        if ($absolute -ge 5.0) { $s.Ge5++ }
        if ($absolute -ge 10.0) {
            $s.Ge10++; $s.XMin = [Math]::Min($s.XMin, $x); $s.XMax = [Math]::Max($s.XMax, $x)
            $s.YMin = [Math]::Min($s.YMin, $y); $s.YMax = [Math]::Max($s.YMax, $y)
            $s.ZMin = [Math]::Min($s.ZMin, $z); $s.ZMax = [Math]::Max($s.ZMax, $z)
            $hotWriter.WriteLine("$($point[0]),$($point[1]),$($point[2]),$($subdomain + 1),$($point[3]),$($point[4]),$($point[5])")
        }
        if ($absolute -ge 20.0) { $s.Ge20++ }
        if ($absolute -ge 30.0) { $s.Ge30++ }
        $sum += $deltaK; $sumSq += $deltaK * $deltaK
        if ($absolute -gt $globalMaxAbs) { $globalMaxAbs = $absolute; $globalMax = "$x,$y,$z,$($subdomain + 1),$ddm,$comsol,$deltaK" }
        if ($ddm -lt $ddmMin) { $ddmMin = $ddm; $ddmMinRow = "$x,$y,$z,$($subdomain + 1),$ddm,$comsol,$deltaK" }
        if ($comsol -lt $comsolMin) { $comsolMin = $comsol; $comsolMinRow = "$x,$y,$z,$($subdomain + 1),$ddm,$comsol,$deltaK" }
        $count++
    }
} finally {
    $rawReader.Dispose(); $comparisonReader.Dispose(); $hotWriter.Dispose()
}
if ($count -ne 802737) { throw "Expected 802737 rows, found $count." }

$statsWriter = [IO.StreamWriter]::new($statsPath, $false, [Text.UTF8Encoding]::new($false))
try {
    $statsWriter.WriteLine('subdomain_id_1based,node_count,mean_error_k,rms_error_k,min_error_k,max_error_k,max_abs_error_k,n_abs_ge_5k,n_abs_ge_10k,n_abs_ge_20k,n_abs_ge_30k,hotspot_xmin_m,hotspot_xmax_m,hotspot_ymin_m,hotspot_ymax_m,hotspot_zmin_m,hotspot_zmax_m')
    for ($index = 0; $index -lt 32; $index++) {
        $s = $stats[$index]; $rms = [Math]::Sqrt($s.SumSq / $s.Nodes)
        if ($s.Ge10 -eq 0) { $bbox = ',,,,,' } else { $bbox = "$($s.XMin),$($s.XMax),$($s.YMin),$($s.YMax),$($s.ZMin),$($s.ZMax)" }
        $statsWriter.WriteLine("$($index + 1),$($s.Nodes),$($s.Sum / $s.Nodes),$rms,$($s.Min),$($s.Max),$($s.MaxAbs),$($s.Ge5),$($s.Ge10),$($s.Ge20),$($s.Ge30),$bbox")
    }
} finally { $statsWriter.Dispose() }

$summaryWriter = [IO.StreamWriter]::new($summaryPath, $false, [Text.UTF8Encoding]::new($false))
try {
    $summaryWriter.WriteLine("POINT_COUNT=$count")
    $summaryWriter.WriteLine("ERROR_MEAN_K=$($sum / $count)")
    $summaryWriter.WriteLine("ERROR_RMS_K=$([Math]::Sqrt($sumSq / $count))")
    $summaryWriter.WriteLine("MAX_ABSOLUTE_ERROR_XYZ_M_SUBDOMAIN_DDM_COMSOL_ERROR=$globalMax")
    $summaryWriter.WriteLine("DDM_TMIN_LOCATION_XYZ_M_SUBDOMAIN_DDM_COMSOL_ERROR=$ddmMinRow")
    $summaryWriter.WriteLine("COMSOL_SAMPLED_TMIN_LOCATION_XYZ_M_SUBDOMAIN_DDM_COMSOL_ERROR=$comsolMinRow")
} finally { $summaryWriter.Dispose() }

Write-Output "Completed $count point comparisons in $OutputDirectory"
