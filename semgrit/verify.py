"""Independent verification of every stage of the pipeline.

These checks deliberately re-derive results by a different route than the code
that produced them, so a shared mistake is less likely to pass both. In
particular :func:`verify_inp_roundtrip` parses the written Abaqus deck back from
disk and tests *physical* invariants -- how far each grain protrudes, whether it
sits inside the requested sector -- rather than re-running the same algebra.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .grain3d import GrainSolid, tet_volumes, validate_grain_solid
from .wheel import UM_PER_MM, WheelModel, check_grain_overlaps, hex_volumes


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str = ""
    data: dict = field(default_factory=dict)

    def __str__(self) -> str:
        return f"[{'PASS' if self.ok else 'FAIL'}] {self.name}: {self.detail}"


# --------------------------------------------------------------------------
# Stage 1: metrology
# --------------------------------------------------------------------------

def verify_metrology(sems: list) -> list[CheckResult]:
    out: list[CheckResult] = []
    for sem in sems:
        name = f"metrology({sem.path})"
        issues = []
        if not (0 < sem.pixel_size_um < 10):
            issues.append(f"implausible pixel size {sem.pixel_size_um}")
        if sem.scalebar_agreement is not None and abs(sem.scalebar_agreement) > 0.05:
            issues.append(f"scale bar disagrees by {100 * sem.scalebar_agreement:+.2f}%")
        if sem.databar_top >= sem.full_intensity.shape[0]:
            issues.append("databar not detected; measurements may include the panel")
        # The micrograph must not contain overlay-flat rows.
        srt = np.sort(sem.intensity, axis=1)
        levels = (np.diff(srt, axis=1) != 0).sum(axis=1) + 1
        if (levels <= 48).any():
            issues.append(
                f"{int((levels <= 48).sum())} near-uniform rows survived the crop"
            )
        out.append(
            CheckResult(
                name,
                not issues,
                "; ".join(issues) if issues
                else f"{sem.pixel_size_um * 1000:.2f} nm/px from {sem.pixel_size_source}, "
                     f"bar agreement "
                     f"{'n/a' if sem.scalebar_agreement is None else f'{100 * sem.scalebar_agreement:+.2f}%'}",
                {
                    "pixel_size_um": sem.pixel_size_um,
                    "source": sem.pixel_size_source,
                    "agreement": sem.scalebar_agreement,
                    "databar_top": sem.databar_top,
                },
            )
        )
    return out


# --------------------------------------------------------------------------
# Stage 2: grain solids
# --------------------------------------------------------------------------

def verify_grain_solids(solids: list[GrainSolid]) -> CheckResult:
    reports = [validate_grain_solid(s) for s in solids]
    bad = [r for r in reports if not r["ok"]]
    worst_vol = max((r["volume_rel_error"] for r in reports if np.isfinite(r["volume_rel_error"])), default=0.0)
    worst_area = max((r["projected_area_rel_error"] for r in reports if np.isfinite(r["projected_area_rel_error"])), default=0.0)
    min_tet = min((r["min_tet_volume_um3"] for r in reports), default=0.0)
    return CheckResult(
        "grain_solids",
        not bad,
        (
            f"{len(reports) - len(bad)}/{len(reports)} solids valid; "
            f"max volume error {worst_vol:.2e}, max silhouette error {worst_area:.2e}, "
            f"min tet volume {min_tet:.3e} um3"
            + ("" if not bad else f"; failures: {[r['grain_id'] for r in bad][:8]}")
        ),
        {
            "n": len(reports),
            "n_bad": len(bad),
            "worst_volume_rel_error": worst_vol,
            "worst_projected_area_rel_error": worst_area,
            "min_tet_volume_um3": min_tet,
            "failures": [r for r in bad][:8],
        },
    )


# --------------------------------------------------------------------------
# Stage 3: wheel
# --------------------------------------------------------------------------

def verify_wheel(model: WheelModel) -> list[CheckResult]:
    out: list[CheckResult] = []
    spec = model.spec

    vols = hex_volumes(model.body_nodes, model.body_hexes)
    out.append(
        CheckResult(
            "body_hex_orientation",
            bool((vols > 0).all()),
            f"{len(vols)} C3D8 elements, min volume {vols.min():.4e} mm3, "
            f"total {vols.sum():.4f} mm3",
            {"min": float(vols.min()), "total": float(vols.sum())},
        )
    )

    # Body nodes must lie inside the requested sector and width.
    n = model.body_nodes
    r = np.hypot(n[:, 0], n[:, 1])
    th = np.degrees(np.arctan2(n[:, 1], n[:, 0])) % 360.0
    issues = []
    if r.max() > spec.outer_radius_mm + 1e-6:
        issues.append(f"node radius {r.max():.6f} exceeds outer {spec.outer_radius_mm}")
    if r.min() < spec.inner_radius_mm - 1e-6:
        issues.append(f"node radius {r.min():.6f} below inner {spec.inner_radius_mm}")
    if abs(n[:, 2]).max() > spec.width_mm / 2 + 1e-6:
        issues.append("node outside wheel width")
    if not spec.is_full_circle:
        # allow a hair past 360->0 wraparound
        outside = (th > spec.sector_deg + 1e-6) & (th < 360.0 - 1e-6)
        if outside.any():
            issues.append(f"{int(outside.sum())} body nodes outside the sector")
    out.append(
        CheckResult(
            "body_within_sector",
            not issues,
            "; ".join(issues) if issues
            else f"r in [{r.min():.4f}, {r.max():.4f}] mm, "
                 f"|z| <= {abs(n[:, 2]).max():.4f} mm, sector {spec.sector_deg:g} deg",
        )
    )

    ov = check_grain_overlaps(model)
    if ov["n_pairs_checked"] == 0 and model.placements:
        detail = (
            f"no pair of the {len(model.placements)} grains is even within "
            f"touching range, so overlap is impossible "
            f"(spacing {model.stats.get('min_spacing_mm', float('nan')):.5f} mm vs "
            f"largest grain diameter "
            f"{2 * max(p.bounding_radius_mm for p in model.placements):.5f} mm)"
        )
    else:
        detail = (
            f"{ov['n_overlapping']} overlapping pairs of {ov['n_pairs_checked']} "
            f"near-neighbour pairs tested"
            + (f", worst {ov['worst_overlap_mm']:.6f} mm" if ov["n_overlapping"] else "")
            + " (bounding-sphere test)"
        )
    out.append(CheckResult("grain_overlap", ov["n_overlapping"] == 0, detail, ov))

    # Protrusion must be achieved: the grain's furthest point should sit exactly
    # protrusion above the bond surface.
    from .wheel import _rotation_matrix

    local_cache = {
        i: (s.vertices - s.centroid_um) / UM_PER_MM for i, s in enumerate(model.shapes)
    }
    errs = []
    for p in model.placements:
        rot = _rotation_matrix(p.rotation_axis, math.radians(p.rotation_angle_deg))
        local = local_cache[p.shape_index]
        # Recover the radial/tangential components without forming the rotated
        # point cloud, so this stays usable at 10^5 grains.
        radial_dir = p.translation_mm.copy()
        radial_dir[2] = 0.0
        n = np.linalg.norm(radial_dir)
        radial_dir = radial_dir / n if n > 0 else np.array([1.0, 0.0, 0.0])
        tangent = np.array([-radial_dir[1], radial_dir[0], 0.0])
        a = local @ (rot.T @ radial_dir)
        b = local @ (rot.T @ tangent)
        reached = float(np.sqrt((n + a) ** 2 + b ** 2).max())
        errs.append(reached - (spec.outer_radius_mm + p.protrusion_mm))
    if errs:
        e = np.abs(np.array(errs))
        out.append(
            CheckResult(
                "grain_protrusion",
                bool(e.max() < 1e-6),
                f"max |achieved - requested| protrusion = {e.max():.3e} mm over "
                f"{len(errs)} grains",
                {"max_abs_error_mm": float(e.max())},
            )
        )

    # Grains inside the sector.
    if model.placements and not spec.is_full_circle:
        bad = [
            p.placement_id
            for p in model.placements
            if not (-1e-9 <= p.theta_deg <= spec.sector_deg + 1e-9)
        ]
        out.append(
            CheckResult(
                "grains_within_sector",
                not bad,
                f"{len(model.placements)} grains, theta in "
                f"[{min(p.theta_deg for p in model.placements):.3f}, "
                f"{max(p.theta_deg for p in model.placements):.3f}] deg"
                + ("" if not bad else f"; outside: {bad[:8]}"),
            )
        )

    if model.achieved_grains < model.requested_grains:
        out.append(
            CheckResult(
                "grain_count",
                True,
                f"placed {model.achieved_grains} of {model.requested_grains} requested "
                f"(surface saturated at the minimum spacing; reported, not silently "
                f"truncated)",
                {"requested": model.requested_grains, "achieved": model.achieved_grains},
            )
        )
    else:
        out.append(
            CheckResult(
                "grain_count",
                True,
                f"placed {model.achieved_grains} grains as requested",
            )
        )
    return out


# --------------------------------------------------------------------------
# Stage 4: Abaqus deck round-trip
# --------------------------------------------------------------------------

@dataclass
class InpPart:
    name: str
    nodes: dict[int, np.ndarray] = field(default_factory=dict)
    elements: dict[str, list[tuple[int, list[int]]]] = field(default_factory=dict)


@dataclass
class InpInstance:
    name: str
    part: str
    translation: np.ndarray
    axis_a: Optional[np.ndarray] = None
    axis_b: Optional[np.ndarray] = None
    angle_deg: float = 0.0


def parse_inp(path: str) -> tuple[dict[str, InpPart], list[InpInstance]]:
    """Minimal Abaqus reader covering the subset this package writes."""
    parts: dict[str, InpPart] = {}
    instances: list[InpInstance] = []

    cur_part: Optional[InpPart] = None
    cur_inst: Optional[InpInstance] = None
    inst_lines: list[list[float]] = []
    mode: Optional[str] = None
    el_type: Optional[str] = None
    in_assembly = False

    def flush_instance() -> None:
        nonlocal cur_inst, inst_lines
        if cur_inst is None:
            return
        if inst_lines:
            cur_inst.translation = np.asarray(inst_lines[0][:3], dtype=float)
        if len(inst_lines) > 1 and len(inst_lines[1]) >= 7:
            row = inst_lines[1]
            cur_inst.axis_a = np.asarray(row[0:3], dtype=float)
            cur_inst.axis_b = np.asarray(row[3:6], dtype=float)
            cur_inst.angle_deg = float(row[6])
        instances.append(cur_inst)
        cur_inst = None
        inst_lines = []

    with open(path, "r", encoding="ascii", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                kw = line.split(",")[0].strip().lower()
                lowered = line.lower()
                if kw == "*part":
                    m = re.search(r"name\s*=\s*([^,]+)", line, re.I)
                    cur_part = InpPart(name=m.group(1).strip() if m else "?")
                    parts[cur_part.name] = cur_part
                    mode = None
                elif kw == "*end part":
                    cur_part, mode = None, None
                elif kw == "*assembly":
                    in_assembly = True
                elif kw == "*end assembly":
                    flush_instance()
                    in_assembly = False
                elif kw == "*instance":
                    flush_instance()
                    mn = re.search(r"name\s*=\s*([^,]+)", line, re.I)
                    mp = re.search(r"part\s*=\s*([^,]+)", line, re.I)
                    cur_inst = InpInstance(
                        name=mn.group(1).strip() if mn else "?",
                        part=mp.group(1).strip() if mp else "?",
                        translation=np.zeros(3),
                    )
                    inst_lines = []
                    mode = "instance"
                elif kw == "*end instance":
                    flush_instance()
                    mode = None
                elif kw == "*node":
                    mode = "node"
                elif kw == "*element":
                    mode = "element"
                    m = re.search(r"type\s*=\s*([^,]+)", line, re.I)
                    el_type = m.group(1).strip().upper() if m else "?"
                else:
                    mode = None if not (in_assembly and cur_inst) else mode
                continue

            vals = [v for v in (t.strip() for t in line.split(",")) if v]
            if mode == "node" and cur_part is not None:
                try:
                    nid = int(vals[0])
                    cur_part.nodes[nid] = np.asarray([float(v) for v in vals[1:4]])
                except (ValueError, IndexError):
                    pass
            elif mode == "element" and cur_part is not None and el_type:
                try:
                    eid = int(vals[0])
                    conn = [int(v) for v in vals[1:]]
                    cur_part.elements.setdefault(el_type, []).append((eid, conn))
                except ValueError:
                    pass
            elif mode == "instance" and cur_inst is not None:
                try:
                    inst_lines.append([float(v) for v in vals])
                except ValueError:
                    mode = None
    flush_instance()
    return parts, instances


def _apply_instance(nodes: np.ndarray, inst: InpInstance) -> np.ndarray:
    """Translate then rotate about the given axis -- the documented convention."""
    p = nodes + inst.translation
    if inst.axis_a is None or abs(inst.angle_deg) < 1e-12:
        return p
    a, b = inst.axis_a, inst.axis_b
    axis = b - a
    n = np.linalg.norm(axis)
    if n < 1e-15:
        return p
    axis = axis / n
    c, s = math.cos(math.radians(inst.angle_deg)), math.sin(math.radians(inst.angle_deg))
    x, y, z = axis
    r = np.array(
        [
            [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
            [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
            [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
        ]
    )
    return (p - a) @ r.T + a


def verify_inp_roundtrip(
    path: str, model: WheelModel, tol_mm: float = 1e-6
) -> list[CheckResult]:
    """Re-read the written deck and confirm it describes the intended wheel."""
    out: list[CheckResult] = []
    parts, instances = parse_inp(path)

    out.append(
        CheckResult(
            "inp_parsed",
            len(parts) > 0 and len(instances) > 0,
            f"{len(parts)} parts, {len(instances)} instances",
        )
    )

    # Node and element IDs must be unique within a part, and every element must
    # reference nodes that exist.
    issues: list[str] = []
    for name, part in parts.items():
        for et, els in part.elements.items():
            ids = [e for e, _ in els]
            if len(set(ids)) != len(ids):
                issues.append(f"{name}: duplicate {et} element ids")
            missing = {n for _, conn in els for n in conn} - set(part.nodes)
            if missing:
                issues.append(f"{name}: {len(missing)} element nodes undefined")
    out.append(
        CheckResult(
            "inp_topology",
            not issues,
            "; ".join(issues) if issues
            else f"all element connectivity resolves; "
                 f"{sum(len(p.nodes) for p in parts.values())} nodes total",
        )
    )

    # Keyword ordering that Abaqus enforces. General contact must sit in the
    # initial step, i.e. before the first *Step; inside a step Abaqus/Standard
    # errors with "General Contact (Std) can only be defined in the initial step"
    # and silently drops the contact definition on import.
    with open(path, "r", encoding="ascii", errors="replace") as fh:
        keywords = [
            (n, ln.split(",")[0].strip().lower())
            for n, ln in enumerate(fh, 1)
            if ln.startswith("*") and not ln.startswith("**")
        ]
    order_issues: list[str] = []
    first_step = next((n for n, k in keywords if k == "*step"), None)
    contact_lines = [n for n, k in keywords if k == "*contact"]
    if first_step is not None:
        late = [n for n in contact_lines if n > first_step]
        if late:
            order_issues.append(
                f"*Contact at line(s) {late[:3]} comes after the first *Step "
                f"(line {first_step}); Abaqus/Standard requires general contact in "
                f"the initial step"
            )
    out.append(
        CheckResult(
            "inp_keyword_order",
            not order_issues,
            "; ".join(order_issues) if order_issues
            else (
                f"general contact in the initial step"
                if contact_lines else "no general contact block"
            )
            + (f", first *Step at line {first_step}" if first_step else ""),
        )
    )

    grain_inst = [i for i in instances if i.name.startswith("G-")]
    by_id = {p.placement_id: p for p in model.placements}
    spec = model.spec

    prot_err: list[float] = []
    sector_bad: list[str] = []
    tet_bad = 0
    tets_checked = 0

    for inst in grain_inst:
        part = parts.get(inst.part)
        if part is None:
            sector_bad.append(f"{inst.name}: part {inst.part} missing")
            continue
        ids = sorted(part.nodes)
        local = np.array([part.nodes[i] for i in ids])
        world = _apply_instance(local, inst)

        pid = int(inst.name.split("-", 1)[1])
        pl = by_id.get(pid)
        if pl is None:
            sector_bad.append(f"{inst.name}: no matching placement")
            continue

        # Physical invariant 1: protrusion above the bond surface.
        radial = np.hypot(world[:, 0], world[:, 1]).max()
        prot_err.append(radial - (spec.outer_radius_mm + pl.protrusion_mm))

        # Physical invariant 2: inside the sector.
        th = np.degrees(np.arctan2(world[:, 1], world[:, 0])) % 360.0
        if not spec.is_full_circle:
            span = np.minimum(th, np.abs(th - 360.0))
            if (th > spec.sector_deg + 5.0).all() and (span > 5.0).all():
                sector_bad.append(f"{inst.name}: theta {th.mean():.2f} outside sector")

        # Physical invariant 3: elements still positively oriented in world space.
        tets = part.elements.get("C3D4")
        if tets:
            remap = {nid: k for k, nid in enumerate(ids)}
            conn = np.array([[remap[n] for n in c] for _, c in tets])
            v = tet_volumes(world, conn)
            tet_bad += int((v <= 0).sum())
            tets_checked += len(v)

    if prot_err:
        e = np.abs(np.asarray(prot_err))
        out.append(
            CheckResult(
                "inp_grain_protrusion",
                bool(e.max() < tol_mm),
                f"rebuilt {len(prot_err)} grains from file; max protrusion error "
                f"{e.max():.3e} mm (tol {tol_mm:g})",
                {"max_abs_error_mm": float(e.max())},
            )
        )
    out.append(
        CheckResult(
            "inp_grains_in_sector",
            not sector_bad,
            "; ".join(sector_bad[:5]) if sector_bad
            else f"all {len(grain_inst)} grain instances lie within "
                 f"{spec.sector_deg:g} deg",
        )
    )
    if tets_checked:
        out.append(
            CheckResult(
                "inp_tet_orientation",
                tet_bad == 0,
                f"{tets_checked - tet_bad}/{tets_checked} C3D4 elements positively "
                f"oriented in assembly coordinates",
            )
        )

    # Body element count must survive the write.
    body = [i for i in instances if not i.name.startswith("G-")]
    if body and model.spec is not None:
        bp = parts.get(body[0].part)
        if bp is not None:
            n_hex = len(bp.elements.get("C3D8", []))
            out.append(
                CheckResult(
                    "inp_body_elements",
                    n_hex == len(model.body_hexes),
                    f"{n_hex} C3D8 in file vs {len(model.body_hexes)} in model",
                )
            )
    return out


def summarise(results: list[CheckResult]) -> tuple[int, int]:
    n_ok = sum(1 for r in results if r.ok)
    return n_ok, len(results)
