#!/usr/bin/env python3
"""Generate deterministic structured tetrahedral meshes for package15."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


FACE_PATTERNS = (
    (1, 2, 3),
    (0, 3, 2),
    (0, 1, 3),
    (0, 2, 1),
)


@dataclass(frozen=True)
class StructuredMesh:
    vertices: list[tuple[float, float, float]]
    triangles: list[tuple[int, int, int]]
    triangle_entities: list[int]
    tets: list[tuple[int, int, int, int]]
    edges: set[tuple[int, int]]

    @property
    def p2_dofs(self) -> int:
        return len(self.vertices) + len(self.edges)


def vertex_index(i: int, j: int, k: int, ny: int, nz: int) -> int:
    return (i * (ny + 1) + j) * (nz + 1) + k


def signed_six_volume(
    vertices: Sequence[tuple[float, float, float]],
    tet: Sequence[int],
) -> float:
    p0, p1, p2, p3 = (vertices[index] for index in tet)
    a = tuple(p1[d] - p0[d] for d in range(3))
    b = tuple(p2[d] - p0[d] for d in range(3))
    c = tuple(p3[d] - p0[d] for d in range(3))
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def classify_boundary_entity(
    face: Sequence[int],
    vertices: Sequence[tuple[float, float, float]],
    size_mm: Sequence[float],
) -> int:
    tolerance = 1.0e-10 * max(1.0, *size_mm)
    coordinates = [vertices[index] for index in face]
    planes = (
        (0, 0.0, 1),
        (0, size_mm[0], 2),
        (1, 0.0, 3),
        (1, size_mm[1], 4),
        (2, 0.0, 5),
        (2, size_mm[2], 6),
    )
    for axis, value, entity in planes:
        if all(abs(point[axis] - value) <= tolerance for point in coordinates):
            return entity
    raise ValueError(f"Boundary face {tuple(face)} does not lie on the cuboid boundary")


def apply_boundary_patches(
    entity: int,
    face: Sequence[int],
    vertices: Sequence[tuple[float, float, float]],
    patches: Sequence[dict[str, object]],
) -> int:
    for patch in patches:
        if entity != int(patch["base_entity"]):
            continue
        axis = int(patch["axis"])
        breaks = [float(value) for value in patch["breaks_mm"]]
        entities = [int(value) for value in patch["entities"]]
        if len(entities) != len(breaks) + 1:
            raise ValueError("A boundary patch needs one more entity than break")
        values = [vertices[index][axis] for index in face]
        tolerance = 1.0e-10 * max(1.0, *map(abs, values), *map(abs, breaks))
        for split in breaks:
            if min(values) < split - tolerance and max(values) > split + tolerance:
                raise ValueError(
                    f"Boundary face {tuple(face)} crosses unmeshed patch split {split}"
                )
        centroid = sum(values) / len(values)
        segment = sum(centroid > split + tolerance for split in breaks)
        return entities[segment]
    return entity


def build_structured_mesh(
    size_mm: Sequence[float],
    divisions: Sequence[int],
    boundary_patches: Sequence[dict[str, object]] = (),
) -> StructuredMesh:
    if len(size_mm) != 3 or len(divisions) != 3:
        raise ValueError("Cuboid size and divisions must contain exactly three values")
    nx, ny, nz = (int(value) for value in divisions)
    if min(nx, ny, nz) <= 0 or min(size_mm) <= 0.0:
        raise ValueError("Cuboid sizes and mesh divisions must be positive")

    vertices = [
        (
            size_mm[0] * i / nx,
            size_mm[1] * j / ny,
            size_mm[2] * k / nz,
        )
        for i in range(nx + 1)
        for j in range(ny + 1)
        for k in range(nz + 1)
    ]
    tets: list[tuple[int, int, int, int]] = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                v000 = vertex_index(i, j, k, ny, nz)
                v100 = vertex_index(i + 1, j, k, ny, nz)
                v010 = vertex_index(i, j + 1, k, ny, nz)
                v110 = vertex_index(i + 1, j + 1, k, ny, nz)
                v001 = vertex_index(i, j, k + 1, ny, nz)
                v101 = vertex_index(i + 1, j, k + 1, ny, nz)
                v011 = vertex_index(i, j + 1, k + 1, ny, nz)
                v111 = vertex_index(i + 1, j + 1, k + 1, ny, nz)
                tets.extend(
                    (
                        (v000, v100, v110, v111),
                        (v000, v110, v010, v111),
                        (v000, v010, v011, v111),
                        (v000, v011, v001, v111),
                        (v000, v001, v101, v111),
                        (v000, v101, v100, v111),
                    )
                )

    minimum_volume = math.prod(size_mm) / (nx * ny * nz) * 1.0e-12
    for tet in tets:
        if signed_six_volume(vertices, tet) <= minimum_volume:
            raise ValueError(f"Non-positive or degenerate tetrahedron: {tet}")

    unmatched_faces: dict[tuple[int, int, int], tuple[int, int, int]] = {}
    edges: set[tuple[int, int]] = set()
    for tet in tets:
        for first in range(4):
            for second in range(first + 1, 4):
                edge = tuple(sorted((tet[first], tet[second])))
                edges.add(edge)
        for pattern in FACE_PATTERNS:
            oriented = tuple(tet[local] for local in pattern)
            key = tuple(sorted(oriented))
            if key in unmatched_faces:
                del unmatched_faces[key]
            else:
                unmatched_faces[key] = oriented

    triangles: list[tuple[int, int, int]] = []
    triangle_entities: list[int] = []
    for _, face in sorted(unmatched_faces.items()):
        triangles.append(face)
        entity = classify_boundary_entity(face, vertices, size_mm)
        triangle_entities.append(
            apply_boundary_patches(entity, face, vertices, boundary_patches)
        )

    expected_tets = 6 * nx * ny * nz
    expected_triangles = 4 * (nx * ny + nx * nz + ny * nz)
    if len(tets) != expected_tets or len(triangles) != expected_triangles:
        raise ValueError(
            "Structured mesh topology count mismatch: "
            f"tets={len(tets)}/{expected_tets}, "
            f"triangles={len(triangles)}/{expected_triangles}"
        )
    return StructuredMesh(vertices, triangles, triangle_entities, tets, edges)


def format_number(value: float) -> str:
    if value == 0.0:
        return "0"
    return format(value, ".17g")


def write_mphtxt(path: Path, name: str, mesh: StructuredMesh) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="ascii", newline="\n") as output:
        output.write("# Deterministic package15 mesh generated by generate_package15_mesh.py.\n\n")
        output.write("0 1\n1 # number of tags\n# Tags\n")
        output.write(f"{len(name)} {name}\n")
        output.write("1 # number of types\n# Types\n3 obj\n\n")
        output.write("# --------- Object 0 ----------\n\n")
        output.write("0 0 1\n4 Mesh # class\n4 # version\n3 # sdim\n")
        output.write(f"{len(mesh.vertices)} # number of mesh vertices\n")
        output.write("0 # lowest mesh vertex index\n\n# Mesh vertex coordinates\n")
        for vertex in mesh.vertices:
            output.write(" ".join(format_number(value) for value in vertex) + "\n")

        output.write("\n2 # number of element types\n\n# Type #0\n\n")
        output.write("3 tri # type name\n\n3 # number of vertices per element\n")
        output.write(f"{len(mesh.triangles)} # number of elements\n# Elements\n")
        for triangle in mesh.triangles:
            output.write(" ".join(str(value) for value in triangle) + "\n")
        output.write(
            f"\n{len(mesh.triangle_entities)} # number of geometric entity indices\n"
            "# Geometric entity indices\n"
        )
        for entity in mesh.triangle_entities:
            output.write(f"{entity}\n")

        output.write("\n# Type #1\n\n3 tet # type name\n\n")
        output.write("4 # number of vertices per element\n")
        output.write(f"{len(mesh.tets)} # number of elements\n# Elements\n")
        for tet in mesh.tets:
            output.write(" ".join(str(value) for value in tet) + "\n")
        output.write(
            f"\n{len(mesh.tets)} # number of geometric entity indices\n"
            "# Geometric entity indices\n"
        )
        for _ in mesh.tets:
            output.write("1\n")
    temporary.replace(path)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def domain_bounds(
    domain: dict[str, object], templates: dict[str, object]
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    origin = tuple(float(value) for value in domain["origin_mm"])
    template = templates[str(domain["template"])]
    size = tuple(float(value) for value in template["size_mm"])
    return origin, tuple(origin[axis] + size[axis] for axis in range(3))


def face_plane(
    bounds: tuple[Sequence[float], Sequence[float]], entity: int
) -> tuple[int, float, int]:
    # Entity 7 is the right half of the Logic z-min face. It is geometrically
    # coplanar with entity 5 and exists only to satisfy per-interface coverage.
    if entity == 7:
        entity = 5
    axis = (entity - 1) // 2
    side = (entity - 1) % 2
    return axis, bounds[side][axis], side


def configured_contact_area(
    left_bounds: tuple[Sequence[float], Sequence[float]],
    left_entity: int,
    right_bounds: tuple[Sequence[float], Sequence[float]],
    right_entity: int,
) -> float:
    left_axis, left_plane, left_side = face_plane(left_bounds, left_entity)
    right_axis, right_plane, right_side = face_plane(right_bounds, right_entity)
    if left_axis != right_axis or left_side == right_side:
        return 0.0
    tolerance = 1.0e-10 * max(1.0, abs(left_plane), abs(right_plane))
    if abs(left_plane - right_plane) > tolerance:
        return 0.0
    area = 1.0
    for axis in range(3):
        if axis == left_axis:
            continue
        overlap = min(left_bounds[1][axis], right_bounds[1][axis]) - max(
            left_bounds[0][axis], right_bounds[0][axis]
        )
        if overlap <= tolerance:
            return 0.0
        area *= overlap
    return area


def discover_contacts(
    domains: Sequence[dict[str, object]], templates: dict[str, object]
) -> set[tuple[int, int]]:
    bounds = {int(domain["id"]): domain_bounds(domain, templates) for domain in domains}
    contacts: set[tuple[int, int]] = set()
    for left_index, left in enumerate(domains):
        left_id = int(left["id"])
        for right in domains[left_index + 1 :]:
            right_id = int(right["id"])
            found = False
            for left_entity in range(1, 7):
                for right_entity in range(1, 7):
                    if configured_contact_area(
                        bounds[left_id], left_entity, bounds[right_id], right_entity
                    ) > 0.0:
                        contacts.add((min(left_id, right_id), max(left_id, right_id)))
                        found = True
                        break
                if found:
                    break
    return contacts


def validate_geometry(manifest: dict[str, object]) -> None:
    templates = manifest["mesh_templates"]
    domains = manifest["domains"]
    interfaces = manifest["interfaces"]
    ids = [int(domain["id"]) for domain in domains]
    if ids != list(range(15)):
        raise ValueError(f"Package15 domain IDs must be exactly 0..14, got {ids}")
    if len(interfaces) != 18:
        raise ValueError(f"Package15 must define exactly 18 interfaces, got {len(interfaces)}")

    by_id = {int(domain["id"]): domain for domain in domains}
    bounds = {domain_id: domain_bounds(domain, templates) for domain_id, domain in by_id.items()}
    configured_pairs: set[tuple[int, int]] = set()
    for interface in interfaces:
        left = int(interface["left"])
        right = int(interface["right"])
        area = configured_contact_area(
            bounds[left],
            int(interface["left_entity"]),
            bounds[right],
            int(interface["right_entity"]),
        )
        expected = float(interface["area_mm2"])
        if not math.isclose(area, expected, rel_tol=1.0e-12, abs_tol=1.0e-10):
            raise ValueError(
                f"Interface {interface['name']} area is {area} mm^2, expected {expected} mm^2"
            )
        pair = (min(left, right), max(left, right))
        if pair in configured_pairs:
            raise ValueError(f"Duplicate physical interface pair: {pair}")
        configured_pairs.add(pair)

    discovered = discover_contacts(domains, templates)
    if discovered != configured_pairs:
        raise ValueError(
            "Configured contacts do not match cuboid geometry: "
            f"missing={sorted(discovered - configured_pairs)}, "
            f"extra={sorted(configured_pairs - discovered)}"
        )


def generate_profile(
    manifest_path: Path, output_root: Path, profile: str
) -> Path:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    validate_geometry(manifest)
    output_directory = output_root / profile
    output_directory.mkdir(parents=True, exist_ok=True)

    template_results: dict[str, object] = {}
    for name, specification in manifest["mesh_templates"].items():
        size_mm = tuple(float(value) for value in specification["size_mm"])
        divisions = tuple(int(value) for value in specification["profiles"][profile])
        mesh = build_structured_mesh(
            size_mm, divisions, specification.get("boundary_patches", ())
        )
        mesh_path = output_directory / f"{name}.mphtxt"
        checksum = write_mphtxt(mesh_path, f"package15_{profile}_{name}", mesh)
        template_results[name] = {
            "path": str(mesh_path.resolve()),
            "size_mm": list(size_mm),
            "divisions": list(divisions),
            "vertices": len(mesh.vertices),
            "edges": len(mesh.edges),
            "p2_dofs": mesh.p2_dofs,
            "tetrahedra": len(mesh.tets),
            "boundary_triangles": len(mesh.triangles),
            "sha256": checksum,
        }

    domains = manifest["domains"]
    total_p2_dofs = sum(
        int(template_results[str(domain["template"])]["p2_dofs"]) for domain in domains
    )
    total_tetrahedra = sum(
        int(template_results[str(domain["template"])]["tetrahedra"]) for domain in domains
    )
    generated_manifest = {
        "schema_version": 1,
        "geometry_manifest": str(manifest_path.resolve()),
        "profile": profile,
        "domain_count": len(domains),
        "interface_count": len(manifest["interfaces"]),
        "total_expected_p2_dofs": total_p2_dofs,
        "total_tetrahedra": total_tetrahedra,
        "total_interface_area_mm2": sum(
            float(interface["area_mm2"]) for interface in manifest["interfaces"]
        ),
        "templates": template_results,
    }
    generated_path = output_directory / "mesh_manifest.json"
    generated_path.write_text(
        json.dumps(generated_manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    return generated_path


def parse_args() -> argparse.Namespace:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile", choices=("smoke", "medium", "large"), default="smoke"
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=project / "data" / "manifests" / "package15_geometry.json",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=project / "data" / "generated" / "package15",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_args()
    generated = generate_profile(
        arguments.manifest.resolve(), arguments.output_root.resolve(), arguments.profile
    )
    summary = json.loads(generated.read_text(encoding="ascii"))
    print(
        f"package15 {arguments.profile}: domains={summary['domain_count']}, "
        f"interfaces={summary['interface_count']}, "
        f"P2 DOFs={summary['total_expected_p2_dofs']}, "
        f"tetrahedra={summary['total_tetrahedra']}"
    )
    print(f"manifest={generated}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
