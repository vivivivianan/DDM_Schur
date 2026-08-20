[CmdletBinding()]
param(
    [string]$Config = 'E:\tsv_pdn4_ddm32\ddm_input.txt',
    [string]$Waveform = 'E:\tsv_pdn4_ddm32_basis_enrichment\residual_base_m3\tsv_constant_heat_sources.csv',
    [string]$Root = 'E:\tsv_pdn4_metric_orchestration'
)
$ErrorActionPreference='Stop'
$exe = Join-Path (Resolve-Path (Join-Path $PSScriptRoot '..')) 'build\Release\SIPGHeatDDM3D.exe'
New-Item -ItemType Directory -Force -Path $Root | Out-Null
function Invoke-Metric([string]$metric) {
  $out=Join-Path $Root $metric; New-Item -ItemType Directory -Force -Path $out|Out-Null
  $args=@('--transient','--config',$Config,'--mor-transient-generate','--mor-transient-method','local-block-arnoldi','--mor-arnoldi-moments','3','--mor-basis-orthogonalization',$metric,'--mor-construction-traces','global-fom','--mor-interface-rank','19','--mor-arnoldi-rank-tolerance','1e-6','--mor-basis','standard','--mor-transient-dt','1e-9','--mor-transient-t-end','2e-7','--mor-transient-integrator','crank-nicolson','--mor-transient-input',$Waveform,'--mor-transient-output','full-field','--mor-transient-initial-mode','uniform','--mor-transient-initial-temperature','293.15','--mor-transient-production','--mor-interface-krylov','augmented-direct','--local-mor-matrix-free-threshold','0','--schur-local-solve-threads','8','--schur-local-pardiso-threads','16','--no-mor-native-reduced-history','--no-mor-full-residual-fallback','--mor-full-residual-tolerance','1','--mor-transient-save',(Join-Path $out 'local_dynamic_model'),'--output-dir',$out,'--fast-run')
  & $exe @args 2>&1 | Tee-Object -FilePath (Join-Path $out 'run.log')
  if($LASTEXITCODE -ne 0){throw "$metric failed with exit code $LASTEXITCODE"}
}
$env:OMP_NUM_THREADS='8';$env:MKL_NUM_THREADS='16';$env:MKL_DYNAMIC='FALSE';$env:SIPG_SOLVER_WORKERS='8'
Invoke-Metric 'mass'
Invoke-Metric 'k-energy'
