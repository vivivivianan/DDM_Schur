param([string]$Root = 'E:\tsv_pdn4_m3_cn_rom_dt_sweep_v2')
$ErrorActionPreference='Stop'
$cases=@(
 [pscustomobject]@{folder='cn_dt_1ns';dt_ns=1.0},
 [pscustomobject]@{folder='cn_dt_0p5ns';dt_ns=0.5},
 [pscustomobject]@{folder='cn_dt_0p25ns';dt_ns=0.25}
)
$summary=foreach($c in $cases){
  $path=Join-Path $Root ($c.folder+'\local_dynamic_schur_accuracy_by_time.csv')
  $r=Import-Csv $path
  $r | Select-Object step,time_s,@{n='tmin_k';e={[double]$_.local_minimum_k}},@{n='tmax_k';e={[double]$_.local_maximum_k}} |
    Export-Csv (Join-Path $Root ($c.folder+'\rom_minmax_by_time.csv')) -NoTypeInformation
  $mn=$r|Sort-Object {[double]$_.local_minimum_k}|Select-Object -First 1
  $mx=$r|Sort-Object {[double]$_.local_maximum_k} -Descending|Select-Object -First 1
  $last=$r|Select-Object -Last 1
  [pscustomobject]@{method='M3 Local Block-Arnoldi ROM + augmented-direct';integrator='crank-nicolson';dt_ns=$c.dt_ns;samples=$r.Count;rank_total=2784;final_tmin_k=[double]$last.local_minimum_k;final_tmax_k=[double]$last.local_maximum_k;global_tmin_k=[double]$mn.local_minimum_k;global_tmin_time_ns=([double]$mn.time_s*1e9);global_tmax_k=[double]$mx.local_maximum_k;global_tmax_time_ns=([double]$mx.time_s*1e9)}
}
$summary|Export-Csv (Join-Path $Root 'm3_cn_rom_dt_comparison.csv') -NoTypeInformation
$summary|Format-Table -AutoSize
