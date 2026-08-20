param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('Generate', 'InitialConditions', 'Waveforms')]
    [string]$Mode,
    [Parameter(Mandatory = $true)]
    [string]$Exe,
    [Parameter(Mandatory = $true)]
    [string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [System.IO.Path]::GetFullPath($Root)
$buildRoot = [System.IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

function Assert-Less([double]$Value, [double]$Limit, [string]$Message) {
    if (-not ($Value -lt $Limit)) { throw "$Message (value=$Value, limit=$Limit)" }
}

function Invoke-Case([string]$Name, [string]$Waveform, [string]$InitialMode) {
    $output = Join-Path $rootFull $Name
    & $Exe --transient --config configs\two_cube_parametric_h.txt `
        --mor-transient-generate --mor-transient-method local-block-arnoldi `
        --mor-arnoldi-moments 3 --mor-transient-dt 0.1 --mor-transient-t-end 0.5 `
        --mor-transient-waveform $Waveform --mor-transient-initial-mode $InitialMode `
        --output-dir $output --fast-run
    if ($LASTEXITCODE -ne 0) { throw "Local dynamic Schur case $Name failed." }
    $summary = Import-Csv (Join-Path $output 'local_dynamic_schur_summary.csv')
    if ($summary.status -ne 'success' -or [int]$summary.subdomains -ne 2) {
        throw "Local dynamic Schur case $Name did not report success with two subdomains."
    }
    Assert-Less ([double]$summary.space_time_relative_l2) 1.0e-6 `
        "$Name space-time accuracy failed"
    Assert-Less ([double]$summary.maximum_absolute_k) 0.01 `
        "$Name maximum error failed"
    Assert-Less ([double]$summary.maximum_reduced_residual) 1.0e-8 `
        "$Name reduced residual failed"
    if ([int]$summary.dynamic_schur_symbolic_calls -ne 1 -or
        [int]$summary.dynamic_schur_numerical_calls -ne 1) {
        throw "$Name rebuilt the fixed-dt Dynamic Schur factor."
    }
    return $summary
}

Push-Location $project
try {
    switch ($Mode) {
        'Generate' {
            $summary = Invoke-Case 'ambient_step' 'single_step' 'ambient'
            if ([int]$summary.full_interface_dofs -le 0 -or [int]$summary.total_local_rank -le 0) {
                throw 'Local transient model lost the full interface or local bases.'
            }
            $rank = Import-Csv (Join-Path $rootFull 'ambient_step\local_block_arnoldi_rank.csv')
            $domains = @($rank | Select-Object -ExpandProperty subdomain -Unique)
            if ($domains.Count -ne 2) { throw 'Expected two independently generated local bases.' }
            foreach ($domain in $domains) {
                $last = @($rank | Where-Object subdomain -eq $domain)[-1]
                Assert-Less ([double]$last.orthogonality_error) 1.0e-10 `
                    "Subdomain $domain basis lost orthogonality"
                if ([int]$last.symbolic_calls -ne 1 -or [int]$last.numerical_calls -ne 1) {
                    throw "Subdomain $domain did not reuse its local PARDISO factor."
                }
            }
        }
        'InitialConditions' {
            Invoke-Case 'uniform_step' 'single_step' 'uniform' | Out-Null
            $steady = Invoke-Case 'steady_dc' 'multi_step' 'steady'
            Assert-Less ([double]$steady.space_time_relative_l2) 1.0e-9 `
                'Steady-state initial condition lost DC consistency'
        }
        'Waveforms' {
            foreach ($case in @(
                @('pulse', 'rectangular_pulse'),
                @('multilevel', 'piecewise_multilevel'),
                @('duty', 'variable_duty_cycle'),
                @('async', 'asynchronous_hotspots'))) {
                Invoke-Case $case[0] $case[1] 'ambient' | Out-Null
            }
        }
    }
} finally {
    Pop-Location
}
