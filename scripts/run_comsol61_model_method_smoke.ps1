param(
  [Parameter(Mandatory = $true)] [string]$TemplateMph,
  [Parameter(Mandatory = $true)] [string]$MethodCallTag,
  [string]$ComsolRoot = 'D:\Program Files\COMSOL\COMSOL61\Multiphysics',
  [string]$OutputDir = (Join-Path $PSScriptRoot '..\data\generated\two_by_three_comsol61\smoke')
)

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path -LiteralPath $ComsolRoot).Path
$template = (Resolve-Path -LiteralPath $TemplateMph).Path
$output = [IO.Path]::GetFullPath($OutputDir)
$runtime = Join-Path (Split-Path $output -Parent) 'runtime'
foreach ($dir in @($output, (Join-Path $runtime 'recovery'), (Join-Path $runtime 'tmp'), (Join-Path $output 'logs'))) {
  New-Item -ItemType Directory -Force -Path $dir | Out-Null
}
$batchExe = Join-Path $root 'bin\win64\comsolbatch.exe'
$outMph = Join-Path $output 'method_smoke_batch.mph'
$args = @('-recoverydir', (Join-Path $runtime 'recovery'), '-tmpdir', (Join-Path $runtime 'tmp'), '-inputfile', $template, '-methodcall', $MethodCallTag, '-outputfile', $outMph)
$commandText = '"' + $batchExe + '" ' + (($args | ForEach-Object { '"' + $_ + '"' }) -join ' ')
$commandText | Tee-Object -FilePath (Join-Path $output 'logs\method_comsolbatch_command.txt')
& $batchExe @args 1> (Join-Path $output 'logs\method_stdout.txt') 2> (Join-Path $output 'logs\method_stderr.txt')
$exitCode = $LASTEXITCODE
"exit_code=$exitCode" | Tee-Object -FilePath (Join-Path $output 'logs\method_exit_code.txt')
Get-ChildItem -LiteralPath $outMph, (Join-Path $output 'method_smoke.mphtxt') -ErrorAction SilentlyContinue | Select-Object FullName, Length, LastWriteTime | Tee-Object -FilePath (Join-Path $output 'logs\method_output_files.txt')
if ($exitCode -ne 0) { exit $exitCode }
if (-not (Test-Path -LiteralPath (Join-Path $output 'method_smoke.mphtxt'))) { throw 'method_smoke.mphtxt was not created.' }
if ((Get-Item -LiteralPath (Join-Path $output 'method_smoke.mphtxt')).Length -le 0) { throw 'method_smoke.mphtxt is empty.' }
