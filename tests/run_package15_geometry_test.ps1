param(
    [string]$Exe = '',
    [string]$Root = '',
    [ValidateSet('smoke', 'medium')][string]$Profile = 'smoke'
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
if ([string]::IsNullOrWhiteSpace($Exe)) {
    $Exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
}
if ([string]::IsNullOrWhiteSpace($Root)) {
    $Root = Join-Path $project "build\package15_geometry_$Profile"
}
$exeFull = [IO.Path]::GetFullPath($Exe)
$rootFull = [IO.Path]::GetFullPath($Root)
$buildRoot = [IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Geometry-test output must remain under build: $rootFull"
}
if (Test-Path -LiteralPath $rootFull) {
    $resolved = [IO.Path]::GetFullPath((Resolve-Path -LiteralPath $rootFull).Path)
    if (-not $resolved.StartsWith($buildRoot, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to remove output outside build: $resolved"
    }
    Remove-Item -LiteralPath $resolved -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

if (-not (Test-Path -LiteralPath $exeFull -PathType Leaf)) {
    throw "Release executable not found: $exeFull"
}

$generator = Join-Path $project 'tools\generate_package15_mesh.py'
& python $generator --profile $Profile
if ($LASTEXITCODE -ne 0) {
    throw "package15 $Profile mesh generation failed."
}

$generatedManifestPath = Join-Path $project (
    "data\generated\package15\$Profile\mesh_manifest.json")
$geometryManifestPath = Join-Path $project 'data\manifests\package15_geometry.json'
$generated = Get-Content -LiteralPath $generatedManifestPath -Raw | ConvertFrom-Json
$geometry = Get-Content -LiteralPath $geometryManifestPath -Raw | ConvertFrom-Json
$config = Join-Path $project "configs\package15_$Profile.txt"
$solverOutput = Join-Path $rootFull 'loader'

$arguments = @(
    '--transient', '--config', $config,
    '--solvers', 'schur-direct-exact', '--schur-direct-pattern-only',
    '--schur-direct-verify-operator', 'false',
    '--output-dir', $solverOutput, '--fast-run'
)
$lines = @(& $exeFull @arguments 2>&1 | ForEach-Object { [string]$_ })
$lines | Set-Content -LiteralPath (Join-Path $rootFull 'loader.log') -Encoding ascii
if ($LASTEXITCODE -ne 0) {
    throw "package15 loader validation failed with exit code $LASTEXITCODE."
}

$countLine = $lines | Where-Object { $_ -match '^Nodes\(P2\):' } | Select-Object -First 1
if ($null -eq $countLine -or
    $countLine -notmatch '^Nodes\(P2\):\s*(\d+),\s*tets:\s*(\d+),\s*boundary triangles:\s*(\d+)') {
    throw 'Could not parse package15 mesh counts from the executable output.'
}
$actualDofs = [int64]$Matches[1]
$actualTets = [int64]$Matches[2]
$actualBoundaryTriangles = [int64]$Matches[3]
if ($actualDofs -ne [int64]$generated.total_expected_p2_dofs -or
    $actualTets -ne [int64]$generated.total_tetrahedra) {
    throw "Loaded mesh count mismatch: DOFs=$actualDofs, tets=$actualTets."
}

$interfaceLines = @($lines | Where-Object { $_ -match '^\s+interface\s+\d+:' })
if ($interfaceLines.Count -ne 18) {
    throw "Expected 18 interface summaries, found $($interfaceLines.Count)."
}
$areaTolerance = 1.0e-10
$rows = @()
for ($index = 0; $index -lt $interfaceLines.Count; ++$index) {
    $line = $interfaceLines[$index]
    if ($line -notmatch 'interface\s+(\d+):.*matched_overlap=([0-9eE+.-]+).*overlap_ratio=([0-9eE+.-]+).*normal_dot\[min/avg/max\]=([0-9eE+.-]+)/([0-9eE+.-]+)/([0-9eE+.-]+)') {
        throw "Could not parse interface summary: $line"
    }
    $interfaceId = [int]$Matches[1]
    $matchedAreaM2 = [double]$Matches[2]
    $overlapRatio = [double]$Matches[3]
    $normalMinimum = [double]$Matches[4]
    $normalAverage = [double]$Matches[5]
    $normalMaximum = [double]$Matches[6]
    $expected = $geometry.interfaces[$interfaceId]
    $expectedAreaM2 = [double]$expected.area_mm2 * 1.0e-6
    if ([Math]::Abs($matchedAreaM2 - $expectedAreaM2) -gt $areaTolerance -or
        $overlapRatio -lt 0.999999 -or
        [Math]::Abs($normalMinimum + 1.0) -gt 1.0e-10 -or
        [Math]::Abs($normalAverage + 1.0) -gt 1.0e-10 -or
        [Math]::Abs($normalMaximum + 1.0) -gt 1.0e-10) {
        throw "Interface $interfaceId geometry/normal gate failed: $line"
    }
    $rows += [pscustomobject]@{
        interface_id = $interfaceId
        name = [string]$expected.name
        matched_area_m2 = $matchedAreaM2
        expected_area_m2 = $expectedAreaM2
        overlap_ratio = $overlapRatio
        normal_dot_average = $normalAverage
    }
}
$rows | Export-Csv -LiteralPath (Join-Path $rootFull 'interfaces.csv') -NoTypeInformation

$summary = [pscustomobject]@{
    status = 'success'
    profile = $Profile
    domains = 15
    interfaces = $interfaceLines.Count
    p2_dofs = $actualDofs
    tetrahedra = $actualTets
    boundary_triangles = $actualBoundaryTriangles
    total_interface_area_mm2 = [double]$generated.total_interface_area_mm2
}
$summary | Export-Csv -LiteralPath (Join-Path $rootFull 'summary.csv') -NoTypeInformation
$summary | Format-List
