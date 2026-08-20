param([string]$Root = 'E:\tsv_pdn4_fom_only_time_integrator_diagnosis')
$ErrorActionPreference = 'Stop'
$cases = @(
    [pscustomobject]@{ name='CN dt=1 ns'; folder='cn_dt_1ns'; integrator='crank-nicolson'; dt_ns=1.0 },
    [pscustomobject]@{ name='CN dt=0.5 ns'; folder='cn_dt_0p5ns'; integrator='crank-nicolson'; dt_ns=0.5 },
    [pscustomobject]@{ name='CN dt=0.25 ns'; folder='cn_dt_0p25ns'; integrator='crank-nicolson'; dt_ns=0.25 },
    [pscustomobject]@{ name='BE dt=1 ns'; folder='be_dt_1ns'; integrator='backward-euler'; dt_ns=1.0 }
)
$summary = foreach ($c in $cases) {
    $rows = Import-Csv (Join-Path $Root ($c.folder + '\fom_minmax_by_time.csv'))
    $minRow = $rows | Sort-Object {[double]$_.tmin_k} | Select-Object -First 1
    $maxRow = $rows | Sort-Object {[double]$_.tmax_k} -Descending | Select-Object -First 1
    $last = $rows | Select-Object -Last 1
    [pscustomobject]@{
      method=$c.name; integrator=$c.integrator; dt_ns=$c.dt_ns; samples=$rows.Count
      final_tmin_k=[double]$last.tmin_k; final_tmax_k=[double]$last.tmax_k
      global_tmin_k=[double]$minRow.tmin_k; global_tmin_time_ns=([double]$minRow.time_s*1e9)
      global_tmax_k=[double]$maxRow.tmax_k; global_tmax_time_ns=([double]$maxRow.time_s*1e9)
    }
}
$summary | Export-Csv (Join-Path $Root 'fom_time_integrator_comparison.csv') -NoTypeInformation

# Intersect the four trajectories at integer ns for direct like-for-like comparison.
$common = @{}
foreach ($c in $cases) {
  $rows=Import-Csv (Join-Path $Root ($c.folder+'\fom_minmax_by_time.csv'))
  foreach($r in $rows){
    $ns=[math]::Round(([double]$r.time_s)*1e9, 9)
    if([math]::Abs($ns-[math]::Round($ns)) -lt 1e-7){
      if(-not $common.ContainsKey($ns)){$common[$ns]=@{time_ns=$ns}}
      $common[$ns][($c.folder+'_tmin_k')]=[double]$r.tmin_k
      $common[$ns][($c.folder+'_tmax_k')]=[double]$r.tmax_k
    }
  }
}
$common.Values | Sort-Object time_ns | ForEach-Object {[pscustomobject]$_} |
  Export-Csv (Join-Path $Root 'fom_minmax_common_1ns.csv') -NoTypeInformation
$summary | Format-Table -AutoSize
