param(
    [string]$Exe = ".\build\Release\SIPGHeatDDM3D.exe",
    [string]$Root = "results\schur_direct_exact_test"
)

$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force -Path $Root | Out-Null
$output = Join-Path $Root 'two_cube'
if (Test-Path -LiteralPath $output) { Remove-Item -LiteralPath $output -Recurse -Force }
& $Exe --steady --config configs\two_cube_schur.txt --solvers schur-direct-exact `
    --output-dir $output --schur-direct-verify-operator true --schur-direct-random-checks 3
if ($LASTEXITCODE -ne 0) { throw "schur-direct-exact process failed: $LASTEXITCODE" }
$summary = Import-Csv (Join-Path $output 'schur_direct_exact_summary.csv')
$status = ($summary | Where-Object field -eq 'schur_direct_status').value
$operator = [double](($summary | Where-Object field -eq 'operator_relative_error').value)
if ($status -ne 'ready' -and $status -ne 'success') { throw "exact Schur status was $status" }
if ($operator -gt 1.0e-10) { throw "explicit/matrix-free operator check failed: $operator" }
Write-Host "schur-direct-exact two-cube test passed; operator relative error=$operator"
