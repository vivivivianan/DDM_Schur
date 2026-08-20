param(
    [string]$Root = 'E:\tsv_pdn4_m3_cn_rom_dt_sweep_v2'
)
$ErrorActionPreference = 'Stop'
$Project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$Exe = Join-Path $Project 'build\Release\SIPGHeatDDM3D.exe'
$Config = 'E:\tsv_pdn4_ddm32\ddm_input.txt'
$ModelCache = Join-Path $Root 'm3_euclidean_model'
$Waveform = 'E:\tsv_pdn4_ddm32_augmented_direct_cn\tsv_constant_heat_sources.csv'
New-Item -ItemType Directory -Force -Path $Root | Out-Null
$env:OMP_NUM_THREADS='8'; $env:MKL_NUM_THREADS='16'; $env:MKL_DYNAMIC='FALSE'; $env:SIPG_SOLVER_WORKERS='8'

function Invoke-M3CnRom([string]$Name, [double]$Dt, [bool]$BuildModel) {
    $out = Join-Path $Root $Name
    New-Item -ItemType Directory -Force -Path $out | Out-Null
    $argumentList = @(
        '--transient','--config',$Config,
        '--mor-transient-generate','--mor-transient-method','local-block-arnoldi',
        '--mor-arnoldi-moments','3','--mor-construction-traces','global-fom',
        '--mor-interface-rank','19','--mor-arnoldi-rank-tolerance','1e-6',
        '--mor-arnoldi-second-moment-energy','1','--mor-arnoldi-second-moment-max-columns','0',
        '--mor-basis-orthogonalization','euclidean',
        '--mor-transient-dt',$Dt.ToString('R',[Globalization.CultureInfo]::InvariantCulture),
        '--mor-transient-t-end','2e-7','--mor-transient-integrator','crank-nicolson',
        '--mor-transient-input',$Waveform,
        '--mor-transient-initial-mode','uniform','--mor-transient-initial-temperature','293.15',
        '--mor-transient-production','--mor-interface-krylov','augmented-direct',
        '--local-mor-matrix-free-threshold','0','--schur-local-solve-threads','8',
        '--schur-local-pardiso-threads','16','--no-mor-native-reduced-history',
        '--no-mor-full-residual-fallback','--mor-full-residual-tolerance','1.0',
        '--mor-transient-output','max-temperature','--output-dir',$out,'--fast-run'
    )
    if ($BuildModel) { $argumentList += @('--mor-transient-save',$ModelCache) }
    else { $argumentList += @('--mor-transient-load',$ModelCache) }
    $argumentList | Set-Content -LiteralPath (Join-Path $out 'command_args.txt')
    & $Exe @argumentList *>> (Join-Path $out 'run.log')
    if ($LASTEXITCODE -ne 0) { throw "M3 CN ROM $Name failed with exit code $LASTEXITCODE" }
    $byTime = Join-Path $out 'local_dynamic_schur_accuracy_by_time.csv'
    if (-not (Test-Path -LiteralPath $byTime)) { throw "Missing time history: $byTime" }
}

Invoke-M3CnRom 'cn_dt_1ns' 1e-9 $true
Invoke-M3CnRom 'cn_dt_0p5ns' 5e-10 $false
Invoke-M3CnRom 'cn_dt_0p25ns' 2.5e-10 $false
