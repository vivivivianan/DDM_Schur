import csv
import re
from pathlib import Path
import matplotlib.pyplot as plt

case = Path(r"D:\AI agent\codex\DDM_Schur\data\generated\tsv_pdn4_ddm_refined_1um")
comsol_t, comsol_tmax = [], []
with (case / "comsol_dense_tmax_history.csv").open(newline="") as f:
    for row in csv.DictReader(f):
        comsol_t.append(float(row["time_s"]) * 1e9)
        comsol_tmax.append(float(row["tmax_K"]))

pattern = re.compile(r"step\s+(\d+)\s+time=\s*([\deE+\-.]+) s\s+Tmin=.*?Tmax=\s*([\deE+\-.]+)")
ddm_t, ddm_tmax = [], []
for line in (case / "cpp_run_200ns_refined_1um.log").read_text(encoding="utf-16").splitlines():
    match = pattern.search(line)
    if match:
        ddm_t.append(float(match.group(2)) * 1e9)
        ddm_tmax.append(float(match.group(3)))

if len(comsol_t) != 201:
    raise RuntimeError(f"Expected 201 COMSOL times, got {len(comsol_t)}")
if len(ddm_t) < 2:
    raise RuntimeError("DDM Tmax history was not found in solver log")

with (case / "tmax_time_comparison.csv").open("w", newline="") as f:
    w = csv.writer(f); w.writerow(["series", "time_ns", "tmax_K"])
    w.writerows(("COMSOL dense", t, T) for t, T in zip(comsol_t, comsol_tmax))
    w.writerows(("DDM 1 um", t, T) for t, T in zip(ddm_t, ddm_tmax))

plt.style.use("seaborn-v0_8-whitegrid")
fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=180)
ax.plot(comsol_t, comsol_tmax, color="#1f77b4", linewidth=2.4, label="COMSOL dense (706,043 DOF)")
ax.plot(ddm_t, ddm_tmax, color="#d62728", marker="o", markersize=4.2, linewidth=1.8, label="DDM 1 μm P2 (2,264,867 DOF)")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Maximum temperature, Tmax (K)")
ax.set_title("Tmax transient comparison: COMSOL vs DDM")
ax.set_xlim(0, 200)
ax.legend(loc="upper left", frameon=True)
ax.grid(True, alpha=0.28)
fig.tight_layout()
fig.savefig(case / "tmax_time_comparison.png", bbox_inches="tight")
