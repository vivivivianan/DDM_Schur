[CmdletBinding()]
param(
    [string]$Config = 'E:\tsv_pdn4_ddm32\ddm_input.txt',
    [string]$Waveform = 'E:\tsv_pdn4_ddm32_basis_enrichment\residual_base_m3\tsv_constant_heat_sources.csv',
    [string]$Root = 'E:\tsv_pdn4_fom_only_time_integrator_diagnosis'
)
$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe = Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
New-Item -ItemType Directory -Force -Path $Root | Out-Null
$env:OMP_NUM_THREADS='8'; $env:MKL_NUM_THREADS='16'; $env:MKL_DYNAMIC='FALSE'; $env:SIPG_SOLVER_WORKERS='8'

function Invoke-Diagnosis([string]$Name, [string]$Integrator, [double]$Dt) {
    $out = Join-Path $Root $Name
    New-Item -ItemType Directory -Force -Path $out | Out-Null
    $args = @(
        '--transient','--config',$Config,
        '--mor-transient-generate','--mor-transient-method','local-block-arnoldi','--mor-transient-fom-only',
        '--mor-transient-dt',$Dt.ToString('R',[Globalization.CultureInfo]::InvariantCulture),
        '--mor-transient-t-end','2e-7','--mor-transient-integrator',$Integrator,
        '--mor-transient-input',$Waveform,
        '--mor-transient-output','max-temperature',
        '--mor-transient-initial-mode','uniform','--mor-transient-initial-temperature','293.15',
        '--output-dir',$out,'--fast-run')
    & $exe @args *>&1 | Tee-Object -FilePath (Join-Path $out 'run.log')
    if ($LASTEXITCODE -ne 0) { throw "$Name failed with exit code $LASTEXITCODE" }
}

Invoke-Diagnosis 'cn_dt_1ns' 'crank-nicolson' 1e-9
Invoke-Diagnosis 'cn_dt_0p5ns' 'crank-nicolson' 5e-10
Invoke-Diagnosis 'cn_dt_0p25ns' 'crank-nicolson' 2.5e-10
Invoke-Diagnosis 'be_dt_1ns' 'backward-euler' 1e-9
