param(
  [string]$ComsolRoot = 'D:\Program Files\COMSOL\COMSOL61\Multiphysics',
  [string]$OutputDir = (Join-Path $PSScriptRoot '..\data\generated\two_by_three_comsol61\smoke'),
  [switch]$UseDefaultPreferences
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path $ComsolRoot).Path
$output = [IO.Path]::GetFullPath($OutputDir)
$runtime = Join-Path (Split-Path $output -Parent) 'runtime'
foreach ($dir in @($output, (Join-Path $runtime 'prefs'), (Join-Path $runtime 'recovery'), (Join-Path $runtime 'tmp'), (Join-Path $output 'logs'))) {
  New-Item -ItemType Directory -Force -Path $dir | Out-Null
}
& (Join-Path $root 'bin\win64\comsolcompile.exe') (Join-Path $PSScriptRoot 'ComsolBatchSmoke.java')
if ($LASTEXITCODE -ne 0) { throw "COMSOL smoke compilation failed with exit code $LASTEXITCODE" }
$batchExe = Join-Path $root 'bin\win64\comsolbatch.exe'
$batchArgs = @('-recoverydir', (Join-Path $runtime 'recovery'), '-tmpdir', (Join-Path $runtime 'tmp'), '-inputfile', (Join-Path $PSScriptRoot 'ComsolBatchSmoke.class'), '-outputfile', (Join-Path $output 'smoke_batch.mph'))
if (-not $UseDefaultPreferences) { $batchArgs = @('-prefsdir', (Join-Path $runtime 'prefs')) + $batchArgs }
$commandText = '"' + $batchExe + '" ' + (($batchArgs | ForEach-Object { '"' + $_ + '"' }) -join ' ')
$commandText | Tee-Object -FilePath (Join-Path $output 'logs\comsolbatch_command.txt')
& $batchExe @batchArgs 1> (Join-Path $output 'logs\stdout.txt') 2> (Join-Path $output 'logs\stderr.txt')
$exitCode = $LASTEXITCODE
"exit_code=$exitCode" | Tee-Object -FilePath (Join-Path $output 'logs\exit_code.txt')
Get-ChildItem (Join-Path $output 'smoke_cube.mph'), (Join-Path $output 'smoke_cube.mphtxt') -ErrorAction SilentlyContinue | Select-Object FullName,Length,LastWriteTime | Tee-Object -FilePath (Join-Path $output 'logs\output_files.txt')
if ($exitCode -ne 0) { exit $exitCode }
if (-not (Test-Path (Join-Path $output 'smoke_cube.mphtxt'))) { throw 'smoke_cube.mphtxt was not created.' }
if ((Get-Item (Join-Path $output 'smoke_cube.mphtxt')).Length -le 0) { throw 'smoke_cube.mphtxt is empty.' }
Get-FileHash (Join-Path $output 'smoke_cube.mphtxt') -Algorithm SHA256
