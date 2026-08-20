[CmdletBinding()]
param([ValidateRange(1, 64)][int]$Jobs = 8)

# Configure and build only the Release production target. CMake locates oneAPI
# in MKLROOT first and then in the standard Windows installation directory.
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path

cmake -S $project -B (Join-Path $project 'build')
if ($LASTEXITCODE -ne 0) { throw 'CMake configuration failed.' }

cmake --build (Join-Path $project 'build') --config Release --parallel $Jobs
if ($LASTEXITCODE -ne 0) { throw 'Release build failed.' }

Write-Host "Built: $project\build\Release\DynamicSchurProduction.exe"
