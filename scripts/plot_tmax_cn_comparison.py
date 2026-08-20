"""Plot a read-only COMSOL-CN Tmax history against the DDM-CN run log."""

import csv
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main(case_dir: Path, comsol_file: Path | None = None, ddm_log: Path | None = None) -> None:
    comsol_file = comsol_file or case_dir / "comsol_cn_tmax_history.csv"
    ddm_log = ddm_log or case_dir / "cpp_run_cn_200ns_1p5um_rerun.log"
    if not comsol_file.is_file() or not ddm_log.is_file():
        raise FileNotFoundError("Expected COMSOL Tmax CSV and final DDM-CN run log.")

    raw_comsol = []
    with comsol_file.open(newline="") as handle:
        for row in csv.DictReader(handle):
            raw_comsol.append((float(row["time_s"]) * 1e9, float(row["tmax_K"])))
    # COMSOL can retain a duplicate final stored value. Keep its last instance.
    # Times are reported in ns here; a 1e-6 ns key tolerates COMSOL's
    # floating-point endpoint duplication without merging distinct 1 ns steps.
    comsol_by_time = {round(t, 6): value for t, value in raw_comsol}
    comsol_t = np.array(sorted(comsol_by_time), dtype=float)
    comsol_max = np.array([comsol_by_time[t] for t in comsol_t], dtype=float)

    pattern = re.compile(
        r"\[Schwarz-Precond-FGMRES\] step\s+(\d+)\s+time=\s*([\deE+\-.]+) s"
        r"\s+Tmin=.*?Tmax=\s*([\deE+\-.]+)"
    )
    text = ddm_log.read_text(encoding="utf-16")
    ddm = [(0.0, 293.15)]
    for match in pattern.finditer(text):
        ddm.append((float(match.group(2)) * 1e9, float(match.group(3))))
    ddm_by_time = {round(t, 6): value for t, value in ddm}
    ddm_t = np.array(sorted(ddm_by_time), dtype=float)
    ddm_max = np.array([ddm_by_time[t] for t in ddm_t], dtype=float)

    comsol_at_ddm = np.interp(ddm_t, comsol_t, comsol_max)
    delta = ddm_max - comsol_at_ddm
    table = case_dir / "tmax_cn_comparison.csv"
    with table.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_ns", "comsol_cn_tmax_K", "ddm_cn_tmax_K", "ddm_minus_comsol_K"])
        for row in zip(ddm_t, comsol_at_ddm, ddm_max, delta):
            writer.writerow(row)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=180)
    ax.plot(comsol_t, comsol_max, color="#1f77b4", linewidth=2.5,
            label="COMSOL CN (E: dense reference)")
    ddm_label = "DDM CN (32 subdomains, 1.5 µm, P2)" if "ddm32" in str(case_dir).lower() else "DDM CN (1.5 µm, P2)"
    ax.plot(ddm_t, ddm_max, color="#d62728", marker="o", markersize=4.2,
            linewidth=1.7, label=ddm_label)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Maximum temperature, Tmax (K)")
    ax.set_title("Tmax transient comparison: COMSOL CN vs DDM CN")
    ax.set_xlim(float(comsol_t.min()), float(comsol_t.max()))
    ax.legend(loc="upper left", frameon=True)
    ax.grid(True, alpha=0.28)
    fig.tight_layout()
    fig.savefig(case_dir / "tmax_cn_comparison.png", bbox_inches="tight")

    summary = case_dir / "tmax_cn_comparison_summary.txt"
    summary.write_text(
        "COMSOL points (deduplicated): %d\nDDM reported points: %d\n"
        "Final time: %.12g ns\nCOMSOL Tmax: %.12g K\nDDM Tmax: %.12g K\n"
        "DDM-COMSOL final difference: %.12g K\n"
        "Maximum absolute difference at DDM reported times: %.12g K\n"
        % (len(comsol_t), len(ddm_t), ddm_t[-1], comsol_at_ddm[-1], ddm_max[-1],
           delta[-1], np.max(np.abs(delta))),
        encoding="utf-8",
    )


if __name__ == "__main__":
    if len(sys.argv) not in (2, 4):
        raise SystemExit("usage: plot_tmax_cn_comparison.py <output_case_dir> [<comsol_history.csv> <ddm_run.log>]")
    main(Path(sys.argv[1]), *(Path(arg) for arg in sys.argv[2:]))
