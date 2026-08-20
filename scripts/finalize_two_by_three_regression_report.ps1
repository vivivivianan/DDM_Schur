param([string]$OutputDir=(Join-Path $PSScriptRoot '..\data\generated\two_by_three_comsol61'))
$ErrorActionPreference='Stop'
$out=[IO.Path]::GetFullPath($OutputDir); $report=Join-Path $out 'two_by_three_regression_report.txt'; $log=Join-Path $out 'cpp_ddm_regression_run.log'
if(!(Test-Path -LiteralPath $report) -or !(Test-Path -LiteralPath $log)){throw 'Missing generator report or C++ regression log.'}
$t=Get-Content -LiteralPath $log -Raw
function Need([string]$p,[string]$what){if($script:t -notmatch $p){throw "Regression evidence missing: $what"}}
Need 'Domains: 3' 'three MPHTXT meshes loaded'; Need 'COMSOL domain material mappings: 2' 'domain-material mappings';
Need 'Dirichlet conditions: 4' 'all Dirichlet faces'; Need 'Convection/Robin conditions: 3' 'all convection faces'; Need 'Matched inward heat-flux boundary faces: [1-9]' 'heat-flux faces';
for($i=0;$i -lt 3;$i++){Need ("interface {0}: left_area=[0-9.e+-]+ m\^2, right_area=[0-9.e+-]+ m\^2, matched_overlap=[0-9.e+-]+ m\^2, overlap_ratio=1" -f $i) "interface $i precheck"}
Need 'interface_assembly=[0-9.]+' 'interface assembly'; Need '\[Schwarz-Precond-FGMRES\] step 100' 'transient step 100'; Need 'Schwarz-Precond-FGMRES \[success\]' 'solver success'
$step=[regex]::Match($t,'\[Schwarz-Precond-FGMRES\] step 100\s+time=\s*100 s\s+Tmin=\s*([^\s]+)\s+Tmax=\s*([^\s]+)\s+Tavg=\s*([^\s]+)\s+interface_avg_jump=([^\s]+)\s+iterations=(\d+)\s+rel_residual=([^\s]+)')
if(!$step.Success){throw 'Could not parse final transient result.'}
$sum=[regex]::Match($t,'total_iterations=(\d+), final_rel_residual=([^,]+)')
$append=@"

C++ regression execution:
input automatic generation: PASS
C++ input load: PASS
mesh read / 6 domain_material mappings: PASS
physical boundary-condition mapping: PASS
interface precheck (3/3): PASS
nonconforming overlap interface assembly: PASS
transient solve to 100 s: PASS
final step iterations=$($step.Groups[5].Value); final relative residual=$($step.Groups[6].Value)
total iterations=$($sum.Groups[1].Value); Tmin=$($step.Groups[1].Value); Tmax=$($step.Groups[2].Value); Tavg=$($step.Groups[3].Value); interface jump=$($step.Groups[4].Value)

=== TWO-BY-THREE REGRESSION PASSED ===
"@
$base=(Get-Content -LiteralPath $report -Raw) -replace '(?s)\r?\nC\+\+ regression execution:.*$',''
Set-Content -LiteralPath $report -Value ($base.TrimEnd()+"`r`n"+$append)
