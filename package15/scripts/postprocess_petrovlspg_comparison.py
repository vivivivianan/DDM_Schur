import csv
from pathlib import Path

ROOT = Path(r"E:\tsv_pdn4_petrov_lspg")
CASES = {
    # The metric experiment owns the final Galerkin field/timing; this
    # separate same-discretization FOM run owns its FOM accuracy columns.
    "Galerkin M3": (
        Path(r"E:\tsv_pdn4_metric_experiment\euclidean"),
        Path(r"E:\tsv_pdn4_ddm32_augmented_direct_cn\m3_full_200_fom"),
    ),
    "CN Petrov-LSPG": (ROOT / "petrov_200ns", ROOT / "petrov_200ns"),
}

def last_row(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))[-1]

def temperature_stats(path):
    minimum = float("inf")
    maximum = float("-inf")
    below = 0
    # The exported field labels use zero-based subdomain indices: 0=SD1,
    # 31=SD32. This is independent of COMSOL domain numbering.
    sd_min = {0: float("inf"), 31: float("inf")}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            temperature = float(row["temperature_k"])
            minimum = min(minimum, temperature)
            maximum = max(maximum, temperature)
            below += temperature < 293.15
            subdomain = int(row["subdomain"])
            if subdomain in sd_min:
                sd_min[subdomain] = min(sd_min[subdomain], temperature)
    return minimum, maximum, below, sd_min[0], sd_min[31]

rows = []
for method, (directory, fom_directory) in CASES.items():
    summary = last_row(directory / "local_dynamic_schur_summary.csv")
    fom_summary = last_row(fom_directory / "local_dynamic_schur_summary.csv")
    minimum, maximum, below, sd1_minimum, sd32_minimum = temperature_stats(
        directory / "local_dynamic_schur_final_temperature.csv")
    accuracy = last_row(fom_directory / "local_dynamic_schur_accuracy_by_time.csv")
    rows.append({
        "method": method,
        "rank": summary["total_local_rank"],
        "tmin_200ns_k": f"{minimum:.15g}",
        "tmax_200ns_k": f"{maximum:.15g}",
        "nodes_below_293_15k": below,
        "sd1_tmin_200ns_k": f"{sd1_minimum:.15g}",
        "sd32_tmin_200ns_k": f"{sd32_minimum:.15g}",
        "space_time_l2_vs_fom": fom_summary["space_time_relative_l2"],
        "max_point_error_vs_fom_k": fom_summary["maximum_absolute_k"],
        "final_time_l2_vs_fom": accuracy["relative_l2"],
        "final_time_max_error_vs_fom_k": accuracy["maximum_absolute_k"],
        "basis_setup_s": summary["local_basis_setup_seconds"],
        "solve_200steps_s": summary["time_stepping_seconds"],
        "total_s": summary["total_seconds"],
        "status": summary["status"],
    })

target = ROOT / "petrov_lspg_comparison.csv"
with target.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(target)
for row in rows:
    print(row)
