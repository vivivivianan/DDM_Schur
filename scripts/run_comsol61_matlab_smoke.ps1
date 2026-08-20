param(
  [string]$ComsolRoot = 'D:\Program Files\COMSOL\COMSOL61\Multiphysics',
  [string]$MatlabRoot = 'F:\Program Files\MATLAB\R2022b',
  [string]$OutputDir = (Join-Path $PSScriptRoot '..\data\generated\two_by_three_comsol61\smoke'),
  [int]$Port = 20361
)
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path -LiteralPath $ComsolRoot).Path
$matlab = Join-Path ((Resolve-Path -LiteralPath $MatlabRoot).Path) 'bin\matlab.exe'
$output = [IO.Path]::GetFullPath($OutputDir)
$runtime = Join-Path (Split-Path $output -Parent) 'matlab_runtime'
foreach ($dir in @($output, $runtime, (Join-Path $output 'logs'))) { New-Item -ItemType Directory -Force -Path $dir | Out-Null }
$server = Join-Path $root 'bin\win64\comsolmphserver.exe'
$serverOut = Join-Path $output 'logs\matlab_comsol_server_stdout.txt'; $serverErr = Join-Path $output 'logs\matlab_comsol_server_stderr.txt'
$recovery = Join-Path $runtime 'recovery'; $tmp = Join-Path $runtime 'tmp'
New-Item -ItemType Directory -Force -Path $recovery, $tmp | Out-Null
# Start-Process joins ArgumentList; quote paths because this workspace contains spaces.
$serverArgs = "-port $Port -silent -recoverydir `"$recovery`" -tmpdir `"$tmp`""
$proc = Start-Process -FilePath $server -ArgumentList $serverArgs -RedirectStandardOutput $serverOut -RedirectStandardError $serverErr -PassThru -WindowStyle Hidden
try {
  $deadline = (Get-Date).AddSeconds(45); $ready = $false
  while ((Get-Date) -lt $deadline) { if ((Test-NetConnection -ComputerName localhost -Port $Port -InformationLevel Quiet -WarningAction SilentlyContinue)) { $ready = $true; break }; Start-Sleep -Milliseconds 500 }
  if (-not $ready) { throw "COMSOL mphserver did not listen on port $Port." }
  $script = Join-Path $PSScriptRoot 'comsol61_matlab_smoke.m'
  $expr = "addpath('$($PSScriptRoot.Replace('\','/'))'); comsol61_matlab_smoke('$($output.Replace('\','/'))','$((Join-Path $root 'mli').Replace('\','/'))',$Port);"
  $cmd = '"' + $matlab + '" -batch "' + $expr + '"'
  $cmd | Tee-Object -FilePath (Join-Path $output 'logs\matlab_smoke_command.txt')
  & $matlab -batch $expr 1> (Join-Path $output 'logs\matlab_smoke_stdout.txt') 2> (Join-Path $output 'logs\matlab_smoke_stderr.txt')
  $code = $LASTEXITCODE; "exit_code=$code" | Tee-Object -FilePath (Join-Path $output 'logs\matlab_smoke_exit_code.txt')
  if ($code -ne 0) { exit $code }
  Get-Item -LiteralPath (Join-Path $output 'matlab_smoke.mphtxt'), (Join-Path $output 'matlab_smoke.mph') | Select-Object FullName,Length,LastWriteTime | Tee-Object -FilePath (Join-Path $output 'logs\matlab_smoke_files.txt')
} finally { if (-not $proc.HasExited) { Stop-Process -Id $proc.Id -Force } }
