param([string]$Root='E:\tsv_pdn4_fom_spatial_undershoot_diagnosis')
$ErrorActionPreference='Stop'
$cases=Get-ChildItem -LiteralPath $Root -Directory | Where-Object {$_.Name -match '^(penalty_|nitsche_)'} | Sort-Object Name
$summary=foreach($case in $cases){
  $rows=Import-Csv (Join-Path $case.FullName 'fom_minmax_by_time.csv')
  $global=$rows|Sort-Object {[double]$_.tmin_k}|Select-Object -First 1
  $last=$rows|Select-Object -Last 1
  $runtime=(Import-Csv (Join-Path $case.FullName 'fom_time_integrator_summary.csv'))
  $kind=if($case.Name -like 'penalty_*'){'interior_sipg'}else{'boundary_nitsche'}
  $parameter=if($case.Name -match 'penalty_(\d+)'){$Matches[1]}else{($case.Name -replace 'nitsche_boundary_','')}
  [pscustomobject]@{case=$case.Name;type=$kind;parameter=$parameter;Tmin_60ns=[double]$last.tmin_k;global_Tmin=[double]$global.tmin_k;global_Tmin_time_ns=([double]$global.time_s*1e9);Tmin_x=[double]$global.tmin_x_m;Tmin_y=[double]$global.tmin_y_m;Tmin_z=[double]$global.tmin_z_m;subdomain=[int]$global.tmin_subdomain;element_id=[int]$global.tmin_element;boundary_id=[int]$global.tmin_boundary_entity;is_interface=[int]$global.tmin_is_interface;runtime_s=[double]$runtime.total_seconds}
}
$summary|Export-Csv (Join-Path $Root 'fom_spatial_undershoot_diagnosis.csv') -NoTypeInformation

$baseline=Join-Path $Root 'penalty_50_strong_baseline\fom_nodes_below_290k_at_global_min.csv'
$cold=Import-Csv $baseline
$coldSummary=[pscustomobject]@{
  case='penalty_50_strong_baseline'; nodes_below_290k=$cold.Count;
  subdomain_distribution=(($cold|Group-Object subdomain|Sort-Object Name|ForEach-Object {"$($_.Name):$($_.Count)"}) -join ';');
  z_min_m=($cold|Measure-Object z_m -Minimum).Minimum; z_max_m=($cold|Measure-Object z_m -Maximum).Maximum;
  external_boundary_nodes=(($cold|Where-Object {[int]$_.boundary_entity -ge 0}).Count);
  interface_nodes=(($cold|Where-Object {[int]$_.is_interface -eq 1}).Count)
}
$coldSummary|Export-Csv (Join-Path $Root 'baseline_cold_node_distribution.csv') -NoTypeInformation
$summary|Format-Table -AutoSize
$coldSummary|Format-List
