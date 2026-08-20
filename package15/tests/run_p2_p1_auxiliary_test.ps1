param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [IO.Path]::GetFullPath($Root)
$buildRoot = [IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null

Push-Location $project
try {
    & $Exe --steady --config configs\two_cube_schur.txt `
        --solvers p2-p1-aux-pcg,direct --direct-mode spd `
        --pcg-tolerance 1e-12 --max-pcg-iterations 1000 --output-dir $rootFull
    if ($LASTEXITCODE -ne 0) { throw 'P2/P1 auxiliary-space run failed.' }
    $aux = Import-Csv (Join-Path $rootFull 'p2_p1_auxiliary_summary.csv')
    if ($aux.status -ne 'success') { throw "Auxiliary status: $($aux.status)" }
    if ([int]$aux.p_rows_one + [int]$aux.p_rows_two -ne [int]$aux.p2_dofs) {
        throw 'P has unclassified or multiply classified rows.'
    }
    foreach ($field in @('constant_error','linear_x_error','linear_y_error','linear_z_error',
                         'a2_symmetry_error','ac_symmetry_error','preconditioner_symmetry_error')) {
        if ([double]$aux.$field -ge 1e-12) { throw "$field failed: $($aux.$field)" }
    }
    if ([double]$aux.ac_positive_energy -le 0 -or
        [double]$aux.preconditioner_positive_energy -le 0) {
        throw 'Auxiliary coarse matrix or preconditioner is not positive.'
    }
    $comparison = Import-Csv (Join-Path $rootFull 'solver_comparison.csv') |
        Where-Object solver -eq 'P2-P1-Auxiliary-PCG'
    if ([double]$comparison.relative_l2_diff_vs_global -ge 1e-9 -or
        [double]$comparison.max_abs_diff_vs_global -ge 1e-5) {
        throw 'Auxiliary PCG/direct accuracy gate failed.'
    }
    $residual = Import-Csv (Join-Path $rootFull 'rram_true_residual_check.csv') |
        Where-Object solver -eq 'P2-P1-Auxiliary-PCG'
    if ([double]$residual.true_relative_residual -ge 1e-10) {
        throw 'Auxiliary PCG true residual gate failed.'
    }
} finally {
    Pop-Location
}
