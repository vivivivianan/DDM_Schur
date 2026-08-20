param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('Generate','Projection','Assembly','Enrichment','Flux','Corrected',
        'Serialization','Fingerprint','TenCube','TemplateReuse','Rram','Chiplet')]
    [string]$Mode,
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Root,
    [string]$Config = ''
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$rootFull = [System.IO.Path]::GetFullPath($Root)
$buildRoot = [System.IO.Path]::GetFullPath((Join-Path $project 'build'))
if (-not $rootFull.StartsWith($buildRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "CTest output root must remain under build: $rootFull"
}
New-Item -ItemType Directory -Force -Path $rootFull | Out-Null
$generated = Join-Path $rootFull 'generated'

function Assert-Less([double]$Value, [double]$Limit, [string]$Message) {
    if (-not ($Value -lt $Limit)) { throw "$Message (value=$Value, limit=$Limit)" }
}

function Invoke-PortCase([string]$Output, [string]$Config, [string[]]$Extra) {
    & $Exe --transient --config $Config --mor-transient-generate `
        --mor-transient-method local-port-block-arnoldi --mor-arnoldi-moments 2 `
        --local-port-rank 16 --mor-transient-dt 0.1 --mor-transient-t-end 0.2 `
        --mor-transient-waveform rectangular_pulse --mor-transient-initial-mode ambient `
        --output-dir $Output --fast-run @Extra
    if ($LASTEXITCODE -ne 0) { throw "Local port case failed: $Output" }
}

function Summary() {
    return Import-Csv (Join-Path $generated 'local_dynamic_schur_summary.csv')
}

function Assert-PortPartition([string]$Directory, $SummaryRow) {
    $ports = @(Import-Csv (Join-Path $Directory 'local_port_rank_by_interface.csv'))
    $covered = ($ports | Measure-Object full_interface_rows -Sum).Sum
    if ($ports.Count -le 0 -or [int]$covered -ne [int]$SummaryRow.full_interface_dofs) {
        throw "Physical ports are not a disjoint, complete SIPG trace partition " +
            "(covered=$covered, interface=$($SummaryRow.full_interface_dofs))."
    }
}

Push-Location $project
try {
    switch ($Mode) {
        'Generate' {
            Invoke-PortCase $generated 'configs\two_cube_parametric_h.txt' @(
                '--local-port-enrichment-rounds','1','--local-port-corrected',
                '--mor-transient-save',(Join-Path $generated 'model'))
            $summary = Summary
            if ($summary.status -ne 'success' -or
                $summary.interface_solver -ne 'port-reduced-dense') {
                throw 'Port-reduced two-cube solve did not report success.'
            }
            Assert-PortPartition $generated $summary
            Assert-Less ([double]$summary.space_time_relative_l2) 1e-6 `
                'Port two-cube temperature accuracy failed'
            Assert-Less ([double]$summary.maximum_absolute_k) 0.01 `
                'Port two-cube maximum error failed'
        }
        'Projection' {
            $summary = Summary
            $rank = @(Import-Csv (Join-Path $generated 'local_port_rank_by_interface.csv'))
            if ($rank.Count -ne 1 -or [int]$rank[0].left_subdomain -eq
                [int]$rank[0].right_subdomain) {
                throw 'Nonmatching physical port did not retain both interface sides.'
            }
            if ([int]$rank[0].full_interface_rows -ne [int]$summary.full_interface_dofs) {
                throw 'Two-sided physical port does not cover the full SIPG interface.'
            }
            Assert-Less ([double]$rank[0].orthogonality_error) 1e-10 `
                'Port projection lost orthogonality'
        }
        'Assembly' {
            $summary = Summary
            if ([int]$summary.port_dimension -le 0 -or
                [int]$summary.port_dimension -ge [int]$summary.full_interface_dofs -or
                [int]$summary.dynamic_schur_numerical_calls -ne 1) {
                throw 'Projected Dynamic Schur was not assembled/factored once.'
            }
            Assert-Less ([double]$summary.maximum_interface_relative_residual) 0.01 `
                'Projected Dynamic Schur residual is too large'
        }
        'Enrichment' {
            $summary = Summary
            $history = @(Import-Csv (Join-Path $generated 'local_port_enrichment_history.csv'))
            if ([int]$summary.enrichment_added_rank -le 0 -or $history.Count -ne 2) {
                throw 'Residual-driven local enrichment did not add per-subdomain modes.'
            }
            Assert-Less ([double]$summary.maximum_full_residual) 0.01 `
                'Residual enrichment did not control the full residual'
        }
        'Flux' {
            $summary = Summary
            Assert-Less ([double]$summary.maximum_fom_rom_flux_relative_l2) 0.01 `
                'Flux-aware enrichment regression failed'
            $flux = @(Import-Csv (Join-Path $generated 'local_dynamic_schur_interface_flux.csv'))
            if ($flux.Count -le 0) { throw 'Detailed FOM/ROM flux output is missing.' }
        }
        'Corrected' {
            $corrected = Import-Csv (Join-Path $generated 'local_port_corrected_accuracy.csv')
            if ($corrected.status -ne 'success') { throw 'Corrected port solve failed.' }
            Assert-Less ([double]$corrected.relative_l2) 1e-8 `
                'Corrected temperature error failed'
            Assert-Less ([double]$corrected.maximum_full_residual) 1e-8 `
                'Corrected full residual failed'
        }
        'Serialization' {
            $source = Join-Path $rootFull 'serialization_source'
            $reload = Join-Path $rootFull 'serialization_reload'
            Invoke-PortCase $source 'configs\two_cube_parametric_h.txt' @(
                '--mor-transient-save',(Join-Path $source 'model'))
            Invoke-PortCase $reload 'configs\two_cube_parametric_h.txt' @(
                '--mor-transient-load',(Join-Path $source 'model'))
            $left = Import-Csv (Join-Path $source 'local_dynamic_schur_final_temperature.csv')
            $right = Import-Csv (Join-Path $reload 'local_dynamic_schur_final_temperature.csv')
            if ($left.Count -ne $right.Count) { throw 'Reloaded port field size differs.' }
            [double]$maximum = 0.0
            for ($row = 0; $row -lt $left.Count; ++$row) {
                $maximum = [Math]::Max($maximum, [Math]::Abs(
                    [double]$left[$row].temperature_k - [double]$right[$row].temperature_k))
            }
            Assert-Less $maximum 1e-9 'Serialized/reloaded port model changed the solution'
        }
        'Fingerprint' {
            $model = Join-Path $generated 'model\local_port_basis.bin'
            $bad = Join-Path $rootFull 'fingerprint_bad'
            & $Exe --transient --config configs\ten_cube_parametric_h.txt `
                --mor-transient-generate --mor-transient-method local-port-block-arnoldi `
                --mor-arnoldi-moments 1 --mor-transient-dt 0.1 `
                --mor-transient-t-end 0.1 --mor-transient-waveform single_step `
                --mor-transient-load $model --output-dir $bad --fast-run
            if ($LASTEXITCODE -eq 0) {
                throw 'Mismatched two-cube port fingerprint/order was not rejected.'
            }
        }
        'TenCube' {
            $ten = Join-Path $rootFull 'ten_cube'
            Invoke-PortCase $ten 'configs\ten_cube_parametric_h.txt' @()
            $summary = Import-Csv (Join-Path $ten 'local_dynamic_schur_summary.csv')
            $rank = @(Import-Csv (Join-Path $ten 'local_port_rank_by_interface.csv'))
            Assert-PortPartition $ten $summary
            if ($rank.Count -ne 9 -or [int]$summary.port_dimension -ge 500) {
                throw 'Ten-cube did not create nine compact independent ports.'
            }
            Assert-Less ([double]$summary.space_time_relative_l2) 1e-5 `
                'Ten-cube port regression failed'
        }
        'TemplateReuse' {
            $ten = Join-Path $rootFull 'template_reuse'
            Invoke-PortCase $ten 'configs\ten_cube_parametric_h.txt' @(
                '--local-port-rank','4','--local-port-temperature-weight','0',
                '--local-port-flux-weight','0','--local-port-residual-weight','0')
            $rank = @(Import-Csv (Join-Path $ten 'local_port_rank_by_interface.csv'))
            $unique = @($rank | Select-Object -ExpandProperty fingerprint -Unique)
            if ($rank.Count -ne 9 -or $unique.Count -le 0) {
                throw 'Port template fingerprint diagnostics are incomplete.'
            }
            foreach ($row in $rank) {
                if ([int]$row.template_reused -eq 1) {
                    $prototype = @($rank | Where-Object {
                        $_.template_id -eq $row.template_id -and
                        [int]$_.template_reused -eq 0 })
                    if ($prototype.Count -ne 1 -or
                        $prototype[0].fingerprint -ne $row.fingerprint) {
                        throw 'Port template reuse bypassed the exact fingerprint gate.'
                    }
                }
            }
        }
        'Rram' {
            if (-not (Test-Path -LiteralPath $Config)) {
                throw "RRAM integration config is missing: $Config"
            }
            $rram = Join-Path $rootFull 'rram26_integration'
            & $Exe --transient --config $Config --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --mor-arnoldi-moments 1 --mor-interface-rank 0 `
                --local-port-rank 16 --local-port-corrected `
                --mor-transient-dt 0.01 --mor-transient-t-end 0.01 `
                --mor-transient-waveform mixed_frequency `
                --mor-transient-initial-mode ambient --output-dir $rram --fast-run
            if ($LASTEXITCODE -ne 0) { throw 'RRAM port integration run failed.' }
            $summary = Import-Csv (Join-Path $rram 'local_dynamic_schur_summary.csv')
            $corrected = Import-Csv (Join-Path $rram 'local_port_corrected_accuracy.csv')
            Assert-PortPartition $rram $summary
            if ([int]$summary.port_dimension -le 0 -or
                [int]$summary.port_dimension -ge [int]$summary.full_interface_dofs -or
                $corrected.status -ne 'success') {
                throw 'RRAM integration compression/corrected-residual gate failed.'
            }
        }
        'Chiplet' {
            if (-not (Test-Path -LiteralPath $Config)) {
                throw "Chiplet integration config is missing: $Config"
            }
            $chiplet = Join-Path $rootFull 'chiplet_flux_integration'
            & $Exe --transient --config $Config --mor-transient-generate `
                --mor-transient-method local-port-block-arnoldi `
                --mor-arnoldi-moments 4 --mor-interface-rank 0 `
                --local-port-rank 200 --local-port-energy-tolerance 1e-20 `
                --local-port-enrichment-rounds 5 --local-port-corrected `
                --mor-transient-dt 0.1 --mor-transient-t-end 1.0 `
                --mor-transient-waveform asynchronous_hotspots `
                --mor-transient-initial-mode ambient --output-dir $chiplet --fast-run
            if ($LASTEXITCODE -ne 0) { throw 'Chiplet port integration run failed.' }
            $summary = Import-Csv (Join-Path $chiplet 'local_dynamic_schur_summary.csv')
            Assert-PortPartition $chiplet $summary
            Assert-Less ([double]$summary.space_time_relative_l2) 1e-4 `
                'Chiplet port temperature regression failed'
            Assert-Less ([double]$summary.maximum_fom_rom_flux_relative_l2) 0.05 `
                'Chiplet port flux regression failed'
            Assert-Less ([double]$summary.maximum_full_residual) 0.05 `
                'Chiplet port full-residual regression failed'
        }
    }
} finally {
    Pop-Location
}
