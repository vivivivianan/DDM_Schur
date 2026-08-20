param(
    [Parameter(Mandatory = $true)][string]$M3RawCsv,
    [Parameter(Mandatory = $true)][double]$ComsolNativeTminK,
    [Parameter(Mandatory = $true)][string]$OutputDirectory
)

$ErrorActionPreference = 'Stop'
$inv = [Globalization.CultureInfo]::InvariantCulture
function ToDouble([string]$value) { [double]::Parse($value, [Globalization.NumberStyles]::Float, $inv) }
New-Item -ItemType Directory -Force -Path $OutputDirectory | Out-Null
$hotPath = Join-Path $OutputDirectory 'm3_points_below_comsol_native_tmin.csv'
$statsPath = Join-Path $OutputDirectory 'm3_undershoot_by_subdomain.csv'
$summaryPath = Join-Path $OutputDirectory 'm3_native_range_comparison.txt'
$stats = @(0..31 | ForEach-Object {
    [pscustomobject]@{ Nodes = 0; Below0 = 0; Below5 = 0; Below10 = 0; Below15 = 0; Tmin = [double]::PositiveInfinity; XMin = [double]::PositiveInfinity; XMax = [double]::NegativeInfinity; YMin = [double]::PositiveInfinity; YMax = [double]::NegativeInfinity; ZMin = [double]::PositiveInfinity; ZMax = [double]::NegativeInfinity }
})
$reader = [IO.StreamReader]::new($M3RawCsv); $writer = [IO.StreamWriter]::new($hotPath, $false, [Text.UTF8Encoding]::new($false))
try {
    if ($reader.ReadLine() -ne 'global_dof,x_m,y_m,z_m,subdomain,temperature_k') { throw 'Unexpected M3 raw CSV header.' }
    $writer.WriteLine('x_m,y_m,z_m,subdomain_id_1based,m3_temperature_k,minimum_possible_error_vs_comsol_k')
    $count = 0; $globalTmin = [double]::PositiveInfinity; $globalMinRow = ''
    while (($line = $reader.ReadLine()) -ne $null) {
        $p = $line.Split(','); if ($p.Count -ne 6) { throw "Malformed row $count." }
        $x = ToDouble $p[1]; $y = ToDouble $p[2]; $z = ToDouble $p[3]; $sid = [int]$p[4]; $temperature = ToDouble $p[5]
        if ($sid -lt 0 -or $sid -ge 32) { throw "Invalid subdomain at row $count." }
        $s = $stats[$sid]; $s.Nodes++; $s.Tmin = [Math]::Min($s.Tmin, $temperature)
        $undershoot = $ComsolNativeTminK - $temperature
        if ($undershoot -gt 0.0) {
            $s.Below0++; $s.XMin = [Math]::Min($s.XMin, $x); $s.XMax = [Math]::Max($s.XMax, $x); $s.YMin = [Math]::Min($s.YMin, $y); $s.YMax = [Math]::Max($s.YMax, $y); $s.ZMin = [Math]::Min($s.ZMin, $z); $s.ZMax = [Math]::Max($s.ZMax, $z)
            if ($undershoot -ge 5.0) { $s.Below5++ }; if ($undershoot -ge 10.0) { $s.Below10++ }; if ($undershoot -ge 15.0) { $s.Below15++ }
            $writer.WriteLine("$($p[1]),$($p[2]),$($p[3]),$($sid + 1),$($p[5]),$undershoot")
        }
        if ($temperature -lt $globalTmin) { $globalTmin = $temperature; $globalMinRow = "$x,$y,$z,$($sid + 1),$temperature" }
        $count++
    }
} finally { $reader.Dispose(); $writer.Dispose() }
if ($count -ne 802737) { throw "Expected 802737 rows, found $count." }
$sw = [IO.StreamWriter]::new($statsPath, $false, [Text.UTF8Encoding]::new($false))
try {
    $sw.WriteLine('subdomain_id_1based,node_count,tmin_m3_k,n_below_comsol_tmin,n_at_least_5k_below,n_at_least_10k_below,n_at_least_15k_below,xmin_m,xmax_m,ymin_m,ymax_m,zmin_m,zmax_m')
    for ($i = 0; $i -lt 32; $i++) { $s = $stats[$i]; $bbox = if ($s.Below0 -eq 0) { ',,,,,' } else { "$($s.XMin),$($s.XMax),$($s.YMin),$($s.YMax),$($s.ZMin),$($s.ZMax)" }; $sw.WriteLine("$($i + 1),$($s.Nodes),$($s.Tmin),$($s.Below0),$($s.Below5),$($s.Below10),$($s.Below15),$bbox") }
} finally { $sw.Dispose() }
$allBelow = ($stats | Measure-Object -Property Below0 -Sum).Sum; $below5 = ($stats | Measure-Object -Property Below5 -Sum).Sum; $below10 = ($stats | Measure-Object -Property Below10 -Sum).Sum; $below15 = ($stats | Measure-Object -Property Below15 -Sum).Sum
[IO.File]::WriteAllText($summaryPath, "COMSOL_NATIVE_TMIN_K=$ComsolNativeTminK`nM3_TMIN_K=$globalTmin`nM3_MINUS_COMSOL_TMIN_K=$($globalTmin-$ComsolNativeTminK)`nM3_TMIN_XYZ_M_SUBDOMAIN=$globalMinRow`nPOINT_COUNT=$count`nPOINTS_BELOW_COMSOL_NATIVE_TMIN=$allBelow`nPOINTS_AT_LEAST_5K_BELOW=$below5`nPOINTS_AT_LEAST_10K_BELOW=$below10`nPOINTS_AT_LEAST_15K_BELOW=$below15`n", [Text.UTF8Encoding]::new($false))
Write-Output "Completed native-range undershoot analysis for $count M3 points."
