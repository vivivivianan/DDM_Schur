[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$Config,
    [Parameter(Mandatory = $true)]
    [string]$Output,
    [Parameter(Mandatory = $true)]
    [double]$EndTime
)

$ErrorActionPreference = 'Stop'
$invariant = [Globalization.CultureInfo]::InvariantCulture
$configFull = [IO.Path]::GetFullPath($Config)
$outputFull = [IO.Path]::GetFullPath($Output)
if (-not (Test-Path -LiteralPath $configFull -PathType Leaf)) {
    throw "DDM input file does not exist: $configFull"
}
if (-not ($EndTime -gt 0.0)) {
    throw 'EndTime must be positive.'
}

# Channel order is exactly the input order parsed by config_io.hpp and used by
# source_parameterization.cpp.  Values are absolute W, not multipliers.
$powers = [Collections.Generic.List[double]]::new()
foreach ($line in Get-Content -LiteralPath $configFull) {
    if ($line -match '^\s*heat_source\s*=\s*(.+)$') {
        $fields = $Matches[1].Split(',') | ForEach-Object { $_.Trim() }
        if ($fields.Count -lt 3) {
            throw "Malformed heat_source entry: $line"
        }
        [void]$powers.Add([double]::Parse(
            $fields[2], [Globalization.CultureInfo]::InvariantCulture))
    }
}
if ($powers.Count -eq 0) {
    throw 'No heat_source entries were found; refusing to create an empty waveform.'
}

$parent = Split-Path -Parent $outputFull
New-Item -ItemType Directory -Force -Path $parent | Out-Null
$header = @('time_s') + (0..($powers.Count - 1) | ForEach-Object { "heat_source_$($_ + 1)_w" })
$powerStrings = $powers | ForEach-Object { $_.ToString('R', $invariant) }
$rows = @(
    ($header -join ','),
    ((@('0') + @($powerStrings)) -join ','),
    ((@($EndTime.ToString('R', $invariant)) + @($powerStrings)) -join ',')
)
[int]$expectedColumns = $powers.Count + 1
foreach ($row in $rows) {
    if ($row.Split(',').Count -ne $expectedColumns) {
        throw "Generated waveform row has an unexpected column count: $($row.Split(',').Count), expected $expectedColumns."
    }
}
[IO.File]::WriteAllLines($outputFull, $rows, [Text.Encoding]::ASCII)
Write-Host "Generated constant physical heat waveform: $outputFull"
Write-Host "Source channels: $($powers.Count); total power W: $(($powers | Measure-Object -Sum).Sum.ToString('R', $invariant))"
