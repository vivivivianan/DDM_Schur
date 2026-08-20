param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('TwoCube','Serialization','Provenance','Steklov','Pcg','Fgmres','Generalized','TargetSolvers','Refinement','ResidualKrylov','Hybrid','HistoryCompression','GlobalRandomized','ProjectionDiagnosis','ProductionBasisStop','TenTransfer','Randomized','RandomizedSerialization','RandomizedReproducibility')]
    [string]$Mode,
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [System.IO.Path]::GetFullPath($Root)
$buildRoot = [System.IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith(
        $buildRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null
$generated = Join-Path $rootFull 'two_cube'
$randomizedGenerated = Join-Path $rootFull 'randomized_two_cube'

function Assert-Less(
    [double]$Value, [double]$Limit, [string]$Message) {
    if (-not ($Value -lt $Limit)) {
        throw "$Message (value=$Value, limit=$Limit)"
    }
}

function Invoke-OptimalCase(
    [string]$Config, [string]$Output, [string[]]$Extra) {
    & $Exe --transient --config $Config --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --port-basis-method optimal-transfer --mor-arnoldi-moments 2 `
        --optimal-port-rank 6 --optimal-port-inner-solver direct `
        --optimal-port-ablation original-mandatory-trace `
        --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
        --mor-transient-waveform single_step `
        --mor-transient-initial-mode ambient --output-dir $Output `
        --fast-run @Extra
    if ($LASTEXITCODE -ne 0) {
        throw "Optimal-port case failed: $Output"
    }
}

function Invoke-RandomizedCase(
    [string]$Output, [string[]]$Extra) {
    & $Exe --transient --config configs\two_cube_parametric_h.txt `
        --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --port-basis-method randomized-transfer `
        --mor-arnoldi-moments 2 --randomized-port-rank 8 `
        --randomized-oversampling 5 `
        --randomized-power-iterations 1 `
        --randomized-seed 12345 `
        --optimal-port-source-mode generalized-dynamic `
        --optimal-port-inner-solver woodbury-exact `
        --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
        --mor-transient-waveform single_step `
        --mor-transient-initial-mode ambient `
        --output-dir $Output --fast-run @Extra
    if ($LASTEXITCODE -ne 0) {
        throw "Randomized-transfer case failed: $Output"
    }
}

function Invoke-HybridCase(
    [string]$Output,
    [string[]]$Extra = @()) {
    & $Exe --transient --config configs\two_cube_parametric_h.txt `
        --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi `
        --port-basis-method hybrid-randomized `
        --mor-arnoldi-moments 2 --randomized-port-rank 8 `
        --randomized-oversampling 5 `
        --randomized-power-iterations 1 `
        --randomized-seed 12345 `
        --residual-krylov-max-rank 4 `
        --residual-krylov-max-sweeps 2 `
        --residual-krylov-tol 1e-4 `
        --residual-krylov-block-size 4 `
        --residual-krylov-probe-mode operator-geometry `
        --optimal-port-source-mode generalized-dynamic `
        --optimal-port-inner-solver woodbury-exact `
        --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
        --mor-transient-waveform single_step `
        --mor-transient-initial-mode ambient `
        --output-dir $Output --fast-run @Extra
    if ($LASTEXITCODE -ne 0) {
        throw "Hybrid randomized-port case failed: $Output"
    }
}

Push-Location $project
try {
    switch ($Mode) {
        'TwoCube' {
            Invoke-OptimalCase 'configs\two_cube_parametric_h.txt' `
                $generated @(
                    '--mor-transient-save',
                    (Join-Path $generated 'model'))
            $summary = Import-Csv (
                Join-Path $generated 'local_dynamic_schur_summary.csv')
            $rank = @(Import-Csv (
                Join-Path $generated 'optimal_port_rank_by_interface.csv'))
            $operator = Import-Csv (
                Join-Path $generated 'optimal_port_operator_diagnostics.csv')
            if ($summary.status -ne 'success' -or
                $summary.port_basis_method -ne 'optimal-transfer' -or
                [int]$summary.port_snapshot_used -ne 0 -or
                [int]$summary.port_fom_used_for_basis -ne 0 -or
                $rank.Count -ne 1 -or
                [int]$rank[0].target_rows -ne
                    [int]$summary.full_interface_dofs) {
                throw 'Two-cube optimal-port structure/provenance failed.'
            }
            Assert-Less ([double]$summary.space_time_relative_l2) 1e-6 `
                'Two-cube optimal-port temperature accuracy failed'
            Assert-Less ([double]$summary.maximum_absolute_k) 0.01 `
                'Two-cube optimal-port maximum error failed'
            Assert-Less ([double]$rank[0].orthogonality_error) 1e-10 `
                'Two-cube optimal-port Gram orthogonality failed'
            Assert-Less ([double]$operator.schur_relative_asymmetry) 1e-10 `
                'Two-cube target Schur symmetry failed'
        }
        'Serialization' {
            $reload = Join-Path $rootFull 'reload'
            Invoke-OptimalCase 'configs\two_cube_parametric_h.txt' `
                $reload @(
                    '--mor-transient-load',
                    (Join-Path $generated 'model'))
            $left = Import-Csv (
                Join-Path $generated 'local_dynamic_schur_summary.csv')
            $right = Import-Csv (
                Join-Path $reload 'local_dynamic_schur_summary.csv')
            # QRCP preserves the serialized port space while the local
            # descriptor/FOM solves are rebuilt.  Their last-bit rounding
            # differs across those independent solves, so use a tolerance
            # that remains far below the accuracy and residual gates.
            Assert-Less ([Math]::Abs(
                [double]$left.space_time_relative_l2 -
                [double]$right.space_time_relative_l2)) 1e-14 `
                'Reload changed optimal-port relative L2'
            Assert-Less ([Math]::Abs(
                [double]$left.maximum_absolute_k -
                [double]$right.maximum_absolute_k)) 3e-12 `
                'Reload changed optimal-port maximum error'
        }
        'Provenance' {
            $output = Join-Path $rootFull 'provenance_mismatch'
            & $Exe --transient --config configs\two_cube_parametric_h.txt `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method optimal-transfer `
                --mor-arnoldi-moments 2 --optimal-port-rank 7 `
                --optimal-port-inner-solver direct `
                --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
                --mor-transient-waveform single_step `
                --mor-transient-initial-mode ambient `
                --mor-transient-load (Join-Path $generated 'model') `
                --output-dir $output --fast-run
            if ($LASTEXITCODE -eq 0) {
                throw 'M8 cache accepted a mismatched rank configuration.'
            }
        }
        'Steklov' {
            $output = Join-Path $rootFull 'steklov'
            & $Exe --transient --config configs\two_cube_parametric_h.txt `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method steklov-schur `
                --mor-arnoldi-moments 2 --local-port-rank 16 `
                --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
                --mor-transient-waveform single_step `
                --mor-transient-initial-mode ambient `
                --output-dir $output --fast-run
            if ($LASTEXITCODE -ne 0) {
                throw 'Steklov baseline failed.'
            }
            $summary = Import-Csv (
                Join-Path $output 'local_dynamic_schur_summary.csv')
            $rank = Import-Csv (
                Join-Path $output 'local_port_rank_by_interface.csv')
            if ($summary.port_basis_method -ne 'steklov-schur' -or
                [int]$summary.port_snapshot_used -ne 0 -or
                [int]$rank.converged_transfer_rank -le 0) {
                throw 'Steklov baseline did not retain spectral modes.'
            }
            Assert-Less ([double]$rank.orthogonality_error) 1e-10 `
                'Steklov Gram orthogonality failed'
        }
        'Pcg' {
            $output = Join-Path $rootFull 'pcg'
            Invoke-OptimalCase 'configs\two_cube_parametric_h.txt' `
                $output @(
                    '--optimal-port-inner-solver','pcg',
                    '--optimal-port-inner-tol','1e-10',
                    '--optimal-port-inner-max-iters','1000')
            $inner = Import-Csv (
                Join-Path $output 'optimal_port_inner_solver.csv')
            if ($inner.actual_solver -ne 'matrix-free-pcg' -or
                $inner.status -ne 'success' -or
                [int]$inner.total_iterations -le 0) {
                throw 'Matrix-free S_tt PCG path was not exercised.'
            }
            Assert-Less ([double]$inner.final_relative_residual) 1e-9 `
                'Matrix-free S_tt PCG residual failed'
        }
        'Fgmres' {
            $output = Join-Path $rootFull 'fgmres'
            Invoke-OptimalCase 'configs\two_cube_parametric_h.txt' `
                $output @(
                    '--optimal-port-inner-solver','fgmres',
                    '--optimal-port-inner-tol','1e-10',
                    '--optimal-port-inner-max-iters','1000')
            $inner = Import-Csv (
                Join-Path $output 'optimal_port_inner_solver.csv')
            if ($inner.actual_solver -ne 'matrix-free-fgmres' -or
                $inner.status -ne 'success' -or
                [int]$inner.total_iterations -le 0) {
                throw 'Matrix-free S_tt FGMRES path was not exercised.'
            }
            Assert-Less ([double]$inner.final_relative_residual) 1e-9 `
                'Matrix-free S_tt FGMRES residual failed'
        }
        'Generalized' {
            $traceInput = Join-Path $rootFull 'trace_plus_input'
            Invoke-OptimalCase 'configs\two_cube_parametric_h.txt' `
                $traceInput @(
                    '--optimal-port-rank','5',
                    '--optimal-port-eigen-tol','1e-4',
                    '--optimal-port-source-mode','trace-plus-input',
                    '--optimal-port-ablation','transfer-only')
            $traceInputRank = Import-Csv (
                Join-Path $traceInput 'optimal_port_rank_by_interface.csv')
            if ([int]$traceInputRank.trace_source_rows -ne 0 -or
                [int]$traceInputRank.input_source_rows -le 0 -or
                [int]$traceInputRank.boundary_source_rows -ne 0 -or
                [int]$traceInputRank.history_source_rows -ne 0) {
                throw 'trace-plus-input source composition is invalid.'
            }
            $output = Join-Path $rootFull 'generalized'
            Invoke-OptimalCase 'configs\two_cube_parametric_h.txt' `
                $output @(
                    '--optimal-port-rank','7',
                    '--optimal-port-eigen-tol','1e-4',
                    '--optimal-port-source-mode','generalized-dynamic',
                    '--optimal-port-ablation',
                    'constant-geometry-generalized')
            $summary = Import-Csv (
                Join-Path $output 'local_dynamic_schur_summary.csv')
            $rank = Import-Csv (
                Join-Path $output 'optimal_port_rank_by_interface.csv')
            $operator = Import-Csv (
                Join-Path $output 'optimal_port_operator_diagnostics.csv')
            if ([int]$rank.trace_source_rows -ne 0 -or
                [int]$rank.input_source_rows -le 0 -or
                [int]$rank.boundary_source_rows -le 0 -or
                [int]$rank.history_source_rows -le 0 -or
                [int]$rank.converged_transfer_rank -le 0 -or
                [int]$operator.eigen_converged -ne 1 -or
                $summary.status -ne 'success' -or
                $summary.optimal_port_source_mode -ne
                    'generalized-dynamic' -or
                [int]$summary.port_snapshot_used -ne 0) {
                throw 'Two-cube generalized source is empty or mislabeled.'
            }
            Assert-Less ([double]$operator.adjoint_relative_error) 1e-10 `
                'Generalized weighted-adjoint consistency failed'
            Assert-Less (
                [double]$operator.explicit_column_reference_error) 1e-10 `
                'Generalized explicit-column consistency failed'
        }
        'TargetSolvers' {
            $output = Join-Path $rootFull 'target_solvers'
            & $Exe --transient --config configs\two_cube_parametric_h.txt `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method optimal-transfer `
                --mor-arnoldi-moments 2 --optimal-port-rank 5 `
                --optimal-port-eigen-max-iters 200 `
                --optimal-port-eigen-tol 1e-8 `
                --optimal-port-inner-tol 1e-12 `
                --optimal-port-inner-max-iters 1000 `
                --optimal-port-source-mode generalized-dynamic `
                --optimal-port-ablation constant-geometry-generalized `
                --optimal-port-target-solver-comparison `
                --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
                --mor-transient-waveform single_step `
                --mor-transient-initial-mode ambient `
                --output-dir $output --fast-run
            if ($LASTEXITCODE -ne 0) {
                throw 'Target-solver comparison execution failed.'
            }
            $rows = @(Import-Csv (
                Join-Path $output 'milestone8_target_solver_comparison.csv'))
            if ($rows.Count -ne 3 -or
                @($rows | Where-Object { $_.status -ne 'pass' }).Count -ne 0) {
                throw 'Target-solver comparison gate failed.'
            }
            foreach ($row in $rows) {
                Assert-Less ([double]$row.relative_solution_difference) 1e-10 `
                    "$($row.solver) solution consistency failed"
                Assert-Less ([double]$row.target_solve_residual) 1e-9 `
                    "$($row.solver) target residual failed"
                Assert-Less ([double]$row.weighted_adjoint_error) 1e-10 `
                    "$($row.solver) weighted adjoint failed"
                Assert-Less ([double]$row.eigenvalue_relative_difference) 1e-8 `
                    "$($row.solver) eigenvalue consistency failed"
                Assert-Less ([double]$row.basis_projector_relative_difference) 1e-8 `
                    "$($row.solver) port subspace consistency failed"
            }
        }
        'Refinement' {
            $validationRows = @()
            foreach ($case in @(
                    [pscustomobject]@{
                        Name = 'two-cube'
                        Config = 'configs\two_cube_parametric_h.txt'
                    },
                    [pscustomobject]@{
                        Name = 'ten-cube'
                        Config = 'configs\ten_cube_parametric_h.txt'
                    })) {
                $output = Join-Path $rootFull (
                    'refinement_' + $case.Name)
                & $Exe --transient --config $case.Config `
                    --mor-transient-generate `
                    --mor-transient-method local-port-block-arnoldi `
                    --port-basis-method optimal-transfer `
                    --mor-arnoldi-moments 1 --optimal-port-rank 6 `
                    --optimal-port-inner-solver woodbury-exact `
                    --optimal-port-inner-tol 1e-10 `
                    --optimal-port-inner-max-iters 1000 `
                    --optimal-port-inner-refinement-max-iters 3 `
                    --optimal-port-inner-refinement-tol 1e-10 `
                    --optimal-port-eigen-max-iters 200 `
                    --optimal-port-eigen-tol 1e-8 `
                    --optimal-port-source-mode generalized-dynamic `
                    --optimal-port-ablation constant-geometry-generalized `
                    --optimal-port-refinement-validation `
                    --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
                    --mor-transient-waveform single_step `
                    --mor-transient-initial-mode ambient `
                    --output-dir $output --fast-run
                if ($LASTEXITCODE -ne 0) {
                    throw "$($case.Name) refinement validation failed."
                }
                if (Test-Path (
                        Join-Path $output 'local_dynamic_schur_summary.csv')) {
                    throw "$($case.Name) refinement validation advanced a transient."
                }
                $rows = @(Import-Csv (
                    Join-Path $output (
                        'milestone8_woodbury_refinement_validation.csv')))
                if ($rows.Count -ne 2 -or
                    @($rows | Where-Object {
                        $_.status -ne 'pass' }).Count -ne 0) {
                    throw "$($case.Name) refinement comparison gate failed."
                }
                $enabled = $rows | Where-Object {
                    [int]$_.refinement_enabled -eq 1 }
                Assert-Less (
                    [double]$enabled.relative_solution_difference_vs_disabled) `
                    1e-10 "$($case.Name) refinement changed the solution"
                Assert-Less ([double]$enabled.target_solve_residual) 1e-9 `
                    "$($case.Name) refined target residual failed"
                Assert-Less ([double]$enabled.weighted_adjoint_error) 1e-10 `
                    "$($case.Name) refined weighted adjoint failed"
                Assert-Less (
                    [double]$enabled.eigenvalue_relative_difference_vs_disabled) `
                    1e-8 "$($case.Name) refinement changed eigenvalues"
                Assert-Less (
                    [double]$enabled.port_subspace_projector_difference_vs_disabled) `
                    1e-8 "$($case.Name) refinement changed the port subspace"
                $validationRows += $rows
            }
            if ($validationRows.Count -ne 4) {
                throw 'Refinement validation did not cover two- and ten-cube.'
            }
        }
        'ResidualKrylov' {
            $output = Join-Path $rootFull 'residual_krylov'
            & $Exe --transient --config configs\two_cube_parametric_h.txt `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method residual-krylov `
                --mor-arnoldi-moments 2 `
                --residual-krylov-max-rank 1 `
                --residual-krylov-max-sweeps 2 `
                --residual-krylov-tol 1e-4 `
                --residual-krylov-block-size 1 `
                --residual-krylov-inner-solver woodbury-exact `
                --optimal-port-inner-solver woodbury-exact `
                --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
                --mor-transient-waveform single_step `
                --mor-transient-initial-mode ambient `
                --output-dir $output --fast-run
            if ($LASTEXITCODE -ne 0) {
                throw 'Residual-Krylov two-cube smoke failed.'
            }
            $summary = Import-Csv (
                Join-Path $output 'local_dynamic_schur_summary.csv')
            $rows = @(Import-Csv (
                Join-Path $output `
                    'residual_krylov_interface_diagnostics.csv'))
            if ($rows.Count -ne 1 -or
                [int]$rows[0].accepted_enrichment_rank -gt 1 -or
                [int]$rows[0].snapshot_used -ne 0 -or
                [int]$rows[0].pod_used -ne 0 -or
                [int]$rows[0].svd_used -ne 0 -or
                $summary.port_basis_method -ne 'residual-krylov') {
                throw 'Residual-Krylov rank/provenance failed.'
            }
            Assert-Less ([double]$rows[0].final_max_probe_residual) 1e-4 `
                'Residual-Krylov tolerance gate failed'
            Assert-Less ([double]$rows[0].target_residual) 1e-9 `
                'Residual-Krylov target solve residual failed'
            Assert-Less ([double]$rows[0].weighted_adjoint_error) 1e-8 `
                'Residual-Krylov weighted adjoint failed'
            Assert-Less ([double]$summary.space_time_relative_l2) 1e-6 `
                'Residual-Krylov temperature accuracy failed'
        }
        'Randomized' {
            Invoke-RandomizedCase $randomizedGenerated @(
                '--randomized-port-compare-optimal',
                '--mor-transient-save',
                (Join-Path $randomizedGenerated 'model'))
            $summary = Import-Csv (
                Join-Path $randomizedGenerated `
                    'local_dynamic_schur_summary.csv')
            $row = Import-Csv (
                Join-Path $randomizedGenerated `
                    'randomized_transfer_interface_diagnostics.csv')
            $comparison = Import-Csv (
                Join-Path $randomizedGenerated `
                    'randomized_port_subspace_comparison.csv')
            if ($summary.status -ne 'success' -or
                $summary.port_basis_method -ne 'randomized-transfer' -or
                [int]$row.accepted_rank -ne 8 -or
                [int]$row.snapshot_used -ne 0 -or
                [int]$row.fom_used_for_basis -ne 0 -or
                [int]$row.target_solve_count -ge
                    [int]$comparison.optimal_target_solve_count) {
                throw 'Randomized-transfer rank/provenance/cost failed.'
            }
            Assert-Less ([double]$row.orthogonality_error) 1e-10 `
                'Randomized weighted orthogonality failed'
            Assert-Less ([double]$row.weighted_adjoint_error) 1e-10 `
                'Randomized weighted adjoint failed'
            Assert-Less ([double]$row.target_residual) 1e-9 `
                'Randomized target residual failed'
            Assert-Less ([double]$summary.space_time_relative_l2) 1e-6 `
                'Randomized temperature accuracy failed'
        }
        'Hybrid' {
            $output = Join-Path $rootFull 'hybrid_two_cube'
            Invoke-HybridCase $output
            $summary = Import-Csv (
                Join-Path $output 'local_dynamic_schur_summary.csv')
            $randomized = Import-Csv (
                Join-Path $output `
                    'randomized_transfer_interface_diagnostics.csv')
            $residual = Import-Csv (
                Join-Path $output `
                    'residual_krylov_interface_diagnostics.csv')
            $rank = Import-Csv (
                Join-Path $output 'local_port_rank_by_interface.csv')
            $composedRank =
                [int]$residual.mandatory_rank_total +
                [int]$residual.accepted_randomized_rank +
                [int]$residual.accepted_enrichment_rank
            if ($summary.status -ne 'success' -or
                $summary.port_basis_method -ne 'hybrid-randomized' -or
                [int]$summary.port_snapshot_used -ne 0 -or
                [int]$summary.port_fom_used_for_basis -ne 0 -or
                [int]$randomized.snapshot_used -ne 0 -or
                [int]$randomized.fom_used_for_basis -ne 0 -or
                [int]$residual.snapshot_used -ne 0 -or
                [int]$residual.pod_used -ne 0 -or
                [int]$residual.svd_used -ne 0 -or
                [int]$residual.mandatory_rank_total -le 0 -or
                [int]$residual.requested_randomized_rank -ne 8 -or
                [int]$rank.total_port_rank -ne $composedRank) {
                throw 'Hybrid port rank/provenance composition failed.'
            }
            Assert-Less ([double]$residual.target_residual) 1e-9 `
                'Hybrid target residual failed'
            Assert-Less ([double]$residual.weighted_adjoint_error) 1e-8 `
                'Hybrid weighted adjoint failed'
            Assert-Less ([double]$summary.space_time_relative_l2) 1e-6 `
                'Hybrid temperature accuracy failed'
        }
        'HistoryCompression' {
            $output = Join-Path $rootFull 'history_compression'
            $model = Join-Path $output 'model.bin'
            Invoke-HybridCase $output @(
                '--history-compression-method',
                'deterministic-rrqr',
                '--history-compression-rank', '16',
                '--history-compression-tolerance', '1e-12',
                '--mor-transient-save', $model)
            $summary = Import-Csv (
                Join-Path $output 'local_dynamic_schur_summary.csv')
            $rows = @(Import-Csv (
                Join-Path $output `
                    'residual_krylov_interface_diagnostics.csv'))
            if ($rows.Count -ne 1) {
                throw 'History compression did not return one physical port.'
            }
            $row = $rows[0]
            if ($summary.status -ne 'success' -or
                $summary.port_basis_method -ne 'hybrid-randomized' -or
                $row.history_compression_method -ne
                    'deterministic-rrqr' -or
                [int]$row.requested_history_rank -ne 16 -or
                [int]$row.compressed_history_rank -gt 16 -or
                [int]$row.compressed_history_rank -gt
                    [int]$row.active_history_channels -or
                [int]$row.deflated_history_channels -ne
                    ([int]$row.raw_history_channels -
                        [int]$row.compressed_history_rank) -or
                [int]$row.history_target_rhs -ne
                    [int]$row.compressed_history_rank -or
                [uint64]$row.history_compression_fingerprint -eq 0 -or
                [int]$row.snapshot_used -ne 0 -or
                [int]$row.pod_used -ne 0 -or
                [int]$row.svd_used -ne 0) {
                throw 'History compression rank/provenance accounting failed.'
            }
            if ([int]$row.active_history_channels -gt 16 -and
                [int]$row.history_target_rhs -ge
                    [int]$row.raw_history_channels) {
                throw 'History compression did not reduce target RHS count.'
            }
            Assert-Less ([double]$row.target_residual) 1e-9 `
                'History compression target residual failed'
            Assert-Less ([double]$row.weighted_adjoint_error) 1e-8 `
                'History compression weighted adjoint failed'
            Assert-Less ([double]$summary.space_time_relative_l2) 1e-6 `
                'History compression temperature accuracy failed'

            $reload = Join-Path $rootFull 'history_compression_reload'
            Invoke-HybridCase $reload @(
                '--history-compression-method',
                'deterministic-rrqr',
                '--history-compression-rank', '16',
                '--history-compression-tolerance', '1e-12',
                '--mor-transient-load', $model)
            $reloadSummary = Import-Csv (
                Join-Path $reload 'local_dynamic_schur_summary.csv')
            Assert-Less ([Math]::Abs(
                [double]$summary.space_time_relative_l2 -
                [double]$reloadSummary.space_time_relative_l2)) 1e-14 `
                'Version-7 compressed basis reload changed temperature'

            $repeat = Join-Path $rootFull 'history_compression_repeat'
            Invoke-HybridCase $repeat @(
                '--history-compression-method',
                'deterministic-rrqr',
                '--history-compression-rank', '16',
                '--history-compression-tolerance', '1e-12')
            $leftFingerprints = @(
                Import-Csv (
                    Join-Path $output `
                        'residual_krylov_interface_diagnostics.csv') |
                Sort-Object {[int]$_.interface_id} |
                ForEach-Object {$_.history_compression_fingerprint})
            $rightFingerprints = @(
                Import-Csv (
                    Join-Path $repeat `
                        'residual_krylov_interface_diagnostics.csv') |
                Sort-Object {[int]$_.interface_id} |
                ForEach-Object {$_.history_compression_fingerprint})
            if (($leftFingerprints -join ';') -ne
                ($rightFingerprints -join ';')) {
                throw 'History compression selection is not deterministic.'
            }
        }
        'GlobalRandomized' {
            $modelOutput = Join-Path $rootFull (
                'global_randomized_model')
            $model = Join-Path $modelOutput 'model'
            Invoke-HybridCase $modelOutput @(
                '--history-compression-method',
                'deterministic-rrqr',
                '--history-compression-rank', '64',
                '--history-compression-tolerance', '1e-12',
                '--mor-transient-save', $model)
            $output = Join-Path $rootFull (
                'global_randomized_rank10')
            Invoke-HybridCase $output @(
                '--history-compression-method',
                'deterministic-rrqr',
                '--history-compression-rank', '64',
                '--history-compression-tolerance', '1e-12',
                '--global-randomized-schur',
                '--global-randomized-rank', '10',
                '--global-randomized-composition', 'global-only',
                '--global-randomized-seed', '12345',
                '--global-randomized-inner-max-iters', '1000',
                '--global-randomized-inner-tol', '1e-10',
                '--mor-transient-load', $model)
            $diagnostics = Import-Csv (
                Join-Path $output (
                    'milestone8_global_randomized_diagnostics.csv'))
            $single = Import-Csv (
                Join-Path $output (
                    'milestone8_global_randomized_single_step.csv'))
            if ($diagnostics.status -ne 'passed' -or
                [int]$diagnostics.global_port_rank -ne 10 -or
                [int]$diagnostics.global_rhs_count -ne 10 -or
                [int]$diagnostics.pardiso_phase33_calls -ne 0 -or
                [int]$diagnostics.snapshot_used -ne 0 -or
                [int]$diagnostics.fom_used_for_basis -ne 0 -or
                [int]$diagnostics.pod_used -ne 0 -or
                [int]$diagnostics.svd_used -ne 0 -or
                [int]$single.active_port_rank -ne 10 -or
                $single.composition -ne 'global-only') {
                throw 'Global randomized rank/cost/provenance failed.'
            }
            Assert-Less (
                [double]$diagnostics.orthogonality_error) 1e-10 `
                'Global randomized weighted orthogonality failed'
            Assert-Less (
                [double]$diagnostics.schur_residual) 1e-8 `
                'Global randomized Schur residual failed'
            Assert-Less (
                [double]$diagnostics.target_solve_residual) 1e-10 `
                'Global randomized target solve residual failed'
        }
        'ProjectionDiagnosis' {
            $modelOutput = Join-Path $rootFull 'projection_model'
            $model = Join-Path $modelOutput 'model'
            Invoke-HybridCase $modelOutput @(
                '--history-compression-method',
                'deterministic-rrqr',
                '--history-compression-rank', '64',
                '--history-compression-tolerance', '1e-12',
                '--mor-transient-save', $model)
            $output = Join-Path $rootFull 'projection_diagnosis'
            Invoke-HybridCase $output @(
                '--history-compression-method',
                'deterministic-rrqr',
                '--history-compression-rank', '64',
                '--history-compression-tolerance', '1e-12',
                '--milestone8-flux-operator-audit',
                '--projection-interface-ids', '0',
                '--mor-transient-load', $model)
            $summary = Import-Csv (
                Join-Path $output `
                    'milestone8_projection_summary.csv')
            $temperature = Import-Csv (
                Join-Path $output `
                    'milestone8_projection_temperature.csv')
            $flux = Import-Csv (
                Join-Path $output `
                    'milestone8_projection_flux.csv')
            $forcing = Import-Csv (
                Join-Path $output `
                    'milestone8_projection_forcing.csv')
            $linearity = @(
                Import-Csv (Join-Path $output `
                    'milestone8_flux_operator_linearity.csv'))
            $affine = @(
                Import-Csv (Join-Path $output `
                    'milestone8_flux_affine_decomposition.csv'))
            $orientation = Import-Csv (
                Join-Path $output `
                    'milestone8_flux_orientation_audit.csv')
            if ($summary.status -ne 'passed' -or
                [int]$summary.requested_interfaces -ne 1 -or
                [int]$summary.completed_interfaces -ne 1 -or
                [int]$summary.transient_steps -ne 1 -or
                [int]$summary.full_field_read -ne 0 -or
                [int]$summary.snapshot_used -ne 0 -or
                [int]$summary.fom_used_for_basis -ne 0 -or
                [int]$summary.pod_used -ne 0 -or
                [int]$summary.svd_used -ne 0 -or
                [int]$temperature.energy_projection_rank -le 0 -or
                $linearity.Count -ne 2 -or
                @($linearity | Where-Object status -ne 'passed').Count -ne 0 -or
                $affine.Count -ne 2 -or
                @($affine | Where-Object status -ne 'passed').Count -ne 0 -or
                $orientation.status -ne 'passed' -or
                [int]$temperature.energy_projection_rank -gt
                    [int]$temperature.local_port_rank) {
                throw 'M8.12 projection diagnosis provenance/rank failed.'
            }
            Assert-Less (
                [double]$temperature.temperature_projection_error) 1e-6 `
                'M8.12 temperature energy projection failed'
            Assert-Less (
                [double]$temperature.galerkin_orthogonality) 1e-8 `
                'M8.12 energy Galerkin orthogonality failed'
            Assert-Less (
                [double]$flux.flux_projection_error) 1e-5 `
                'M8.12 flux projection failed'
            Assert-Less (
                [double]$forcing.forcing_projection_error) 1e-6 `
                'M8.12 forcing projection failed'
            Assert-Less (
                [double]$forcing.target_solve_residual) 1e-9 `
                'M8.12 exact target solve residual failed'
        }
        'ProductionBasisStop' {
            $output = Join-Path $rootFull 'production_basis_stop'
            Invoke-HybridCase $output @(
                '--randomized-port-rank', '16',
                '--history-compression-method',
                'deterministic-rrqr',
                '--history-compression-rank', '64',
                '--history-compression-tolerance', '1e-12',
                '--milestone8-production-basis-only',
                '--mor-transient-save', $output)
            $summary = Import-Csv (
                Join-Path $output `
                    'milestone8_production_basis_summary.csv')
            if ($summary.status -ne 'success' -or
                [int]$summary.physical_interfaces -ne 1 -or
                [int]$summary.history_rank -ne 64 -or
                [int]$summary.randomized_rank -ne 16 -or
                [int]$summary.residual_rank -ne 4 -or
                [int]$summary.total_port_rank -le 0 -or
                [int]$summary.target_solve_count -le 0 -or
                [int]$summary.snapshot_used -ne 0 -or
                [int]$summary.fom_used_for_basis -ne 0 -or
                [int]$summary.full_field_read -ne 0 -or
                [int]$summary.transient_advanced -ne 0 -or
                -not (Test-Path -LiteralPath (
                    Join-Path $output 'local_port_basis.bin'))) {
                throw 'M8 production basis-only stop/provenance failed.'
            }
            if (Test-Path -LiteralPath (
                    Join-Path $output `
                        'local_dynamic_schur_summary.csv')) {
                throw 'Production basis-only mode advanced the transient.'
            }
        }
        'RandomizedSerialization' {
            $reload = Join-Path $rootFull 'randomized_reload'
            Invoke-RandomizedCase $reload @(
                '--mor-transient-load',
                (Join-Path $randomizedGenerated 'model'))
            $left = Import-Csv (
                Join-Path $randomizedGenerated `
                    'local_dynamic_schur_summary.csv')
            $right = Import-Csv (
                Join-Path $reload 'local_dynamic_schur_summary.csv')
            Assert-Less ([Math]::Abs(
                [double]$left.space_time_relative_l2 -
                [double]$right.space_time_relative_l2)) 1e-14 `
                'Randomized reload changed relative L2'
            Assert-Less ([Math]::Abs(
                [double]$left.maximum_absolute_k -
                [double]$right.maximum_absolute_k)) 1e-11 `
                'Randomized reload changed maximum error'
        }
        'RandomizedReproducibility' {
            $repeat = Join-Path $rootFull 'randomized_repeat'
            Invoke-RandomizedCase $repeat @()
            $left = Import-Csv (
                Join-Path $randomizedGenerated `
                    'local_port_rank_by_interface.csv')
            $right = Import-Csv (
                Join-Path $repeat `
                    'local_port_rank_by_interface.csv')
            if ($left.fingerprint -ne $right.fingerprint) {
                throw 'Fixed seed did not reproduce the exact port basis.'
            }
        }
        'TenTransfer' {
            $output = Join-Path $rootFull 'ten_transfer'
            & $Exe --transient --config configs\ten_cube_parametric_h.txt `
                --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --port-basis-method optimal-transfer `
                --mor-arnoldi-moments 1 --optimal-port-rank 8 `
                --optimal-port-inner-solver direct `
                --optimal-port-eigen-max-iters 8 `
                --optimal-port-eigen-tol 1e-4 `
                --mor-transient-dt 0.1 --mor-transient-t-end 0.1 `
                --mor-transient-waveform single_step `
                --mor-transient-initial-mode ambient `
                --output-dir $output --fast-run
            if ($LASTEXITCODE -ne 0) {
                throw 'Ten-cube transfer smoke failed.'
            }
            $summary = Import-Csv (
                Join-Path $output 'local_dynamic_schur_summary.csv')
            $rank = @(Import-Csv (
                Join-Path $output 'optimal_port_rank_by_interface.csv'))
            $operator = @(Import-Csv (
                Join-Path $output 'optimal_port_operator_diagnostics.csv'))
            $target = ($rank | Measure-Object target_rows -Sum).Sum
            $transfer = (
                $rank | Measure-Object converged_transfer_rank -Sum).Sum
            if ($rank.Count -ne 9 -or
                [int]$target -ne [int]$summary.full_interface_dofs -or
                [int]$transfer -le 0 -or
                @($operator | Where-Object {
                    [int]$_.eigen_converged -ne 1 }).Count -ne 0) {
                throw 'Ten-cube transfer ownership/eigen solve failed.'
            }
            Assert-Less ([double](($operator |
                Measure-Object adjoint_relative_error -Maximum).Maximum)) `
                1e-10 'Transfer weighted-adjoint consistency failed'
            Assert-Less ([double](($operator |
                Measure-Object explicit_column_reference_error -Maximum).Maximum)) `
                1e-10 'Transfer explicit-column reference failed'
        }
    }
} finally {
    Pop-Location
}
