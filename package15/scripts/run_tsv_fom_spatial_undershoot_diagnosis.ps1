param([string]$Root = 'E:\tsv_pdn4_fom_spatial_undershoot_diagnosis')
$ErrorActionPreference='Stop'
$project=(Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exe=Join-Path $project 'build\Release\SIPGHeatDDM3D.exe'
$config='E:\tsv_pdn4_ddm32\ddm_input.txt'
$waveform='E:\tsv_pdn4_ddm32_augmented_direct_cn\tsv_constant_heat_sources.csv'
New-Item -ItemType Directory -Force -Path $Root|Out-Null
$env:OMP_NUM_THREADS='8';$env:MKL_NUM_THREADS='16';$env:MKL_DYNAMIC='FALSE';$env:SIPG_SOLVER_WORKERS='8'
function Invoke-Fom([string]$Name,[double]$InternalPenalty,[string]$DirichletMethod,[double]$BoundaryPenalty){
  $out=Join-Path $Root $Name;New-Item -ItemType Directory -Force -Path $out|Out-Null
  $a=@('--transient','--config',$config,'--mor-transient-generate','--mor-transient-method','local-block-arnoldi','--mor-transient-fom-only',
       '--mor-transient-dt','1e-9','--mor-transient-t-end','6e-8','--mor-transient-integrator','crank-nicolson',
       '--mor-transient-input',$waveform,'--mor-transient-initial-mode','uniform','--mor-transient-initial-temperature','293.15',
       '--penalty-factor',$InternalPenalty.ToString('R',[Globalization.CultureInfo]::InvariantCulture),
       '--dirichlet-method',$DirichletMethod,'--nitsche-penalty-factor',$BoundaryPenalty.ToString('R',[Globalization.CultureInfo]::InvariantCulture),
       '--output-dir',$out,'--fast-run')
  $a|Set-Content -LiteralPath (Join-Path $out 'command_args.txt')
  & $exe @a *>> (Join-Path $out 'run.log')
  if($LASTEXITCODE -ne 0){throw "FOM case $Name failed with exit code $LASTEXITCODE"}
}
# Experiment B: strong Dirichlet is the actual application baseline; only the
# interior SIPG penalty changes.
Invoke-Fom 'penalty_25_strong' 25 'strong' 15
Invoke-Fom 'penalty_50_strong_baseline' 50 'strong' 15
Invoke-Fom 'penalty_100_strong' 100 'strong' 15
Invoke-Fom 'penalty_200_strong' 200 'strong' 15
# Experiment C: the production configuration is strong Dirichlet, hence these
# are deliberately separate Nitsche diagnostics.  Internal SIPG remains 50.
Invoke-Fom 'nitsche_boundary_7p5' 50 'nitsche' 7.5
Invoke-Fom 'nitsche_boundary_15' 50 'nitsche' 15
Invoke-Fom 'nitsche_boundary_30' 50 'nitsche' 30
Invoke-Fom 'nitsche_boundary_60' 50 'nitsche' 60
