"""Read-only audit of SIPG interfaces and physical-boundary mesh coverage."""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def split_values(line: str):
    return [part.strip() for part in line.split("=", 1)[1].split(",")]


def main(case: Path) -> None:
    output = case / "ddm_output"
    boundary = {}
    with (output / "boundary_summary.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            boundary[(int(row["subdomain"]), int(row["boundary_entity"]))] = (
                int(row["triangles"]), float(row["area_m2"])
            )

    conditions = defaultdict(list)
    for raw in (case / "ddm_input.txt").read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#") or "=" not in raw:
            continue
        key = raw.split("=", 1)[0].strip()
        if key in {"dirichlet", "convection", "heat_flux", "insulation"}:
            values = split_values(raw)
            conditions[key].append((int(values[0]), int(values[1]), values[2:]))

    lines = ["DDM interface and physical-boundary discretization audit", ""]
    lines.append("Physical boundary conditions read from ddm_input.txt")
    lines.append("type, conditions, P2 boundary triangles, area_m2, missing_entities")
    for key in ("dirichlet", "convection", "heat_flux", "insulation"):
        entries = conditions[key]
        triangles = area = missing = 0
        for subdomain, entity, _ in entries:
            item = boundary.get((subdomain, entity))
            if item is None:
                missing += 1
            else:
                triangles += item[0]
                area += item[1]
        lines.append(f"{key},{len(entries)},{triangles},{area:.17g},{missing}")
    lines.append("")

    lines.append("Interfaces from actual SIPG overlap integration")
    lines.append("left,right,entities,triangles_per_side,area_m2,overlap_ratio,normal_dot_range")
    grouped = defaultdict(lambda: {"entities": 0, "left": 0, "right": 0, "area": 0.0})
    with (output / "interface_summary.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["left_subdomain"], row["right_subdomain"])
            group = grouped[key]
            group["entities"] += 1
            group["left"] += int(row["left_boundary_triangles"])
            group["right"] += int(row["right_boundary_triangles"])
            group["area"] += float(row["integration_area_m2"])
    with (output / "interface_build_summary.csv").open(newline="") as handle:
        build = {(row["left_subdomain"], row["right_subdomain"]): row for row in csv.DictReader(handle)}
    for key, group in sorted(grouped.items()):
        row = build[key]
        lines.append(
            f"{key[0]},{key[1]},{group['entities']},{group['left']}/{group['right']},"
            f"{group['area']:.17g},{row['overlap_ratio_min_side']},"
            f"{row['normal_dot_min']}..{row['normal_dot_max']}"
        )
    lines.append("")

    penalty = next(csv.DictReader((output / "interface_penalty_stats.csv").open(newline="")))
    lines.append("SIPG penalty statistics")
    lines.append("face_pairs=%s, eta_F[min/avg/max]=%s/%s/%s, h_F[min/avg/max]=%s/%s/%s" % (
        penalty["face_pair_count"], penalty["min_eta_F"], penalty["avg_eta_F"], penalty["max_eta_F"],
        penalty["min_h_F"], penalty["avg_h_F"], penalty["max_h_F"],
    ))
    text = (case / "cpp_run_cn_200ns_1p5um_rerun.log").read_text(encoding="utf-16")
    final = re.findall(r"step\s+200\s+time=.*?interface_avg_jump=([\deE+\-.]+)", text)
    if final:
        lines.append(f"final interface_avg_jump_K={final[-1]}")

    report = case / "interface_boundary_discretization_audit.txt"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(report)
    print("\n".join(lines))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
