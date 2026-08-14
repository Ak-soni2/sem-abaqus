"""Standalone verification suite for the semgrit pipeline.

Run with:  python verify_all.py [image_glob]

Covers unit-level maths, integration across all SEM images, wheel assembly at
several sector angles, Abaqus export round-trips, and determinism. Exits non-zero
if anything fails, so it can gate a commit.
"""

from __future__ import annotations

import glob
import hashlib
import math
import os
import shutil
import sys
import tempfile
import warnings

warnings.filterwarnings("ignore")

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from semgrit.abaqus import AbaqusExportOptions, write_inp
from semgrit.grain3d import (
    HeightModel,
    LoftProfile,
    _signed_area,
    build_grain_library,
    convex_dihedral_angles,
    fillet_polyline,
    round_outline_corners,
    prism_to_tets,
    surface_from_tets,
    tet_volumes,
    triangulate_clean,
    validate_grain_solid,
)
from semgrit.measure import grain_statistics, measure_all, measure_grain
from semgrit.metrology import SemImage, load_sem_image, parse_length_um, snap_to_nice
from semgrit.segment import Segmentation, SegmentationParams, segment_grains
from semgrit.step import (
    StepExportOptions,
    StepWriter,
    check_step_solids,
    read_step_faceted_solids,
    signed_volume as step_signed_volume,
    surface_from_hexes,
    wheel_triangles,
    write_binary_stl,
    write_wheel_step,
)
from semgrit.verify import (
    verify_grain_solids,
    verify_inp_roundtrip,
    verify_metrology,
    verify_wheel,
)
from semgrit.wheel import (
    GrainPopulationSpec,
    WheelSpec,
    build_rim_mesh,
    build_wheel,
    hex_volumes,
    jittered_grid_2d,
    matrix_to_axis_angle,
    _rotation_matrix,
)

PASS, FAIL = [], []


def check(name: str, ok: bool, detail: str = "") -> None:
    (PASS if ok else FAIL).append(name)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


def section(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


# --------------------------------------------------------------------------
def unit_tests() -> None:
    section("1. UNIT TESTS (analytic ground truth)")

    # --- unit parsing ---
    cases = [
        ("Image Pixel Size = 29.30 nm", 0.02930),
        ("Width = 30.00 \xb5m", 30.0),
        ("Pixel Size = 1.5 mm", 1500.0),
        ("Foo = 2.0 pm", 2e-6),
    ]
    ok = all(
        v is not None and abs(v - exp) < 1e-9 * max(exp, 1)
        for text, exp in cases
        for v in [parse_length_um(text)]
    )
    check("length/unit parsing", ok, f"{len(cases)} formats incl. non-ASCII micron sign")

    snaps = [(1.992, 2.0), (0.996, 1.0), (4.9, 5.0), (0.0203, 0.02), (3.3, None)]
    ok = all(snap_to_nice(v) == exp for v, exp in snaps)
    check("scale-label 1-2-5 snapping", ok, "including rejection of a non-1-2-5 value")

    # --- prism -> tets ---
    V = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]], float)
    tets = np.array(prism_to_tets([0, 1, 2], [3, 4, 5]))
    vols = tet_volumes(V, tets)
    check("prism splits into 3 tets summing to prism volume",
          abs(abs(vols).sum() - 0.5) < 1e-12, f"sum|v|={abs(vols).sum():.6f} exact=0.5")

    rng = np.random.default_rng(0)
    P = rng.random((60000, 3))
    pts = P[P[:, 0] + P[:, 1] <= 1]
    cnt = np.zeros(len(pts), int)
    for t in tets:
        A = np.array([V[t[1]] - V[t[0]], V[t[2]] - V[t[0]], V[t[3]] - V[t[0]]]).T
        lam = np.linalg.solve(A, (pts - V[t[0]]).T).T
        cnt += ((lam >= -1e-9).all(1) & (lam.sum(1) <= 1 + 1e-9)).astype(int)
    check("prism decomposition is an exact partition",
          bool((cnt == 1).all()), f"{np.mean(cnt == 1):.4f} of points in exactly one tet")

    # --- conformity of shared quad faces ---
    def faces_of(ts):
        fd = ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))
        return {tuple(sorted((t[a], t[b], t[c]))) for t in ts for a, b, c in fd}
    shared = faces_of(prism_to_tets([0, 1, 2], [4, 5, 6])) & faces_of(
        prism_to_tets([1, 2, 3], [5, 6, 7]))
    check("adjacent prisms agree on the shared-face diagonal", len(shared) == 2,
          f"{len(shared)} shared triangles (2 = conforming)")

    # --- surface extraction closes ---
    faces = surface_from_tets(tets)
    e = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    _, counts = np.unique(np.sort(e, axis=1), axis=0, return_counts=True)
    check("boundary surface of a tet mesh is closed", bool((counts == 2).all()),
          f"{len(faces)} faces, all edges shared twice")

    # --- axis-angle round trip ---
    worst = 0.0
    rng = np.random.default_rng(1)
    for _ in range(2000):
        ax = rng.normal(size=3); ax /= np.linalg.norm(ax)
        ang = rng.uniform(0, 180)
        R = _rotation_matrix(ax, math.radians(ang))
        a2, g2 = matrix_to_axis_angle(R)
        worst = max(worst, float(np.abs(R - _rotation_matrix(a2, math.radians(g2))).max()))
    for ax in ([1, 0, 0], [0, 0, 1], [1, 1, 1], [0.3, -0.7, 0.2]):
        ax = np.array(ax, float); ax /= np.linalg.norm(ax)
        R = _rotation_matrix(ax, math.pi)
        a2, g2 = matrix_to_axis_angle(R)
        worst = max(worst, float(np.abs(R - _rotation_matrix(a2, math.radians(g2))).max()))
    check("axis-angle decomposition round-trips", worst < 1e-6,
          f"max matrix error {worst:.2e} over 2000 random + 4 exact-180 cases")

    # --- jittered grid spacing guarantee ---
    allok = True
    for w, h, n, d in [(26.18, 10, 10472, 0.02), (314.16, 10, 40000, 0.02), (50, 10, 500, 0.8)]:
        p, ach = jittered_grid_2d(w, h, n, d, np.random.default_rng(5))
        if len(p) != n:
            allok = False
        q = p if len(p) <= 3000 else p[np.random.default_rng(2).choice(len(p), 3000, replace=False)]
        dx = np.abs(q[:, None, 0] - q[None, :, 0]); dx = np.minimum(dx, w - dx)
        dy = q[:, None, 1] - q[None, :, 1]
        dist = np.sqrt(dx ** 2 + dy ** 2); np.fill_diagonal(dist, np.inf)
        if dist.min() < d - 1e-9:
            allok = False
    check("jittered grid: exact count and guaranteed min spacing", allok)

    # --- measurement against analytic shapes ---
    ps = 0.05
    H = W = 400
    yy, xx = np.mgrid[0:H, 0:W]

    def measure(mask):
        import cv2
        img = (mask * 200).astype(np.uint8)
        sem = SemImage(path="synthetic", intensity=img, full_intensity=img,
                       pixel_size_um=ps, pixel_size_source="override", databar_top=H)
        dist = (cv2.distanceTransform(mask.astype(np.uint8), cv2.DIST_L2, 5) * ps).astype(np.float32)
        seg = Segmentation(labels=mask.astype(np.int32), foreground=mask, distance_um=dist,
                           n_seeds=1, pixel_size_um=ps, params=SegmentationParams(),
                           border_labels=set())
        return measure_grain(1, 1, seg, sem)

    R0 = 100.0
    m = measure(((xx - 200) ** 2 + (yy - 200) ** 2) <= R0 ** 2)
    errs = {
        "area": abs(m.area_um2 - math.pi * (R0 * ps) ** 2) / (math.pi * (R0 * ps) ** 2),
        "eqdiam": abs(m.equivalent_diameter_um - 2 * R0 * ps) / (2 * R0 * ps),
        "feret_max": abs(m.feret_max_um - 2 * R0 * ps) / (2 * R0 * ps),
        "perimeter": abs(m.perimeter_um - 2 * math.pi * R0 * ps) / (2 * math.pi * R0 * ps),
        "circularity": abs(m.circularity - 1.0),
    }
    check("circle measurements within 1%", max(errs.values()) < 0.01,
          ", ".join(f"{k} {100 * v:+.2f}%" for k, v in errs.items()))

    a, b = 150.0, 60.0
    m = measure((((xx - 200) / a) ** 2 + ((yy - 200) / b) ** 2) <= 1)
    e2 = {
        "feret_max": abs(m.feret_max_um - 2 * a * ps) / (2 * a * ps),
        "feret_min": abs(m.feret_min_um - 2 * b * ps) / (2 * b * ps),
        "aspect": abs(m.aspect_ratio - a / b) / (a / b),
        "major": abs(m.major_axis_um - 2 * a * ps) / (2 * a * ps),
    }
    check("ellipse Feret/aspect/axes within 1%", max(e2.values()) < 0.01,
          ", ".join(f"{k} {100 * v:+.2f}%" for k, v in e2.items()))

    sq = (np.abs(xx - 200) < 100) & (np.abs(yy - 200) < 100)
    con = sq.copy(); con[200:300, 200:300] = False
    m = measure(con)
    # hull of the L clips the notch diagonally: (199^2 - 100^2/2)
    exact_sol = (199 ** 2 - 100 ** 2) / (199 ** 2 - 100 ** 2 / 2)
    check("solidity of a concave L-shape", abs(m.solidity - exact_sol) / exact_sol < 0.01,
          f"measured {m.solidity:.4f} vs exact {exact_sol:.4f}")

    m2 = measure(sq)
    check("square corner angle = 90 deg", abs(m2.min_corner_angle_deg - 90) < 1.0,
          f"{m2.min_corner_angle_deg:.2f} deg")

    # --- convexity classification, and the rounding primitives ---
    cube_v0 = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    cq = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
          (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    ct = []
    for q in cq:
        ct += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
    d = convex_dihedral_angles(cube_v0, np.array(ct))
    check("convex-edge detection on a cube", (np.abs(d - 90) < 1e-6).sum() == 12,
          f"{(np.abs(d - 90) < 1e-6).sum()} edges at 90 deg (expect exactly 12)")

    # A square blunted by r must lose exactly the four corner offcuts:
    # area = s^2 - (4 - pi) r^2.
    side = 10.0
    sq = np.array([[0, 0], [side, 0], [side, side], [0, side]], float)
    ok, detail = True, []
    for r in (0.5, 1.0, 2.0):
        rounded, got = round_outline_corners(sq, r, arc_segments=12)
        exact = side * side - (4 - math.pi) * r * r
        area = abs(_signed_area(rounded))
        err = abs(area - exact) / exact
        detail.append(f"r={r}: {100 * err:.2f}%")
        if err > 0.01 or abs(got - r) / r > 0.02:
            ok = False
    check("outline blunting removes the exact analytic corner area", ok,
          "square of side 10, area should be s^2-(4-pi)r^2 -> " + ", ".join(detail))

    # A filleted right angle: the arc must sit at the requested radius from the
    # corner's bisector, i.e. arc points are r from the fillet centre.
    poly = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)]
    pts, got_r = fillet_polyline(poly, 2.0, arc_segments=8)
    centre = np.array([8.0, 2.0])          # for a 90 deg corner at (10,0), r=2
    arc = pts[1:-1]
    radii = np.linalg.norm(arc - centre, axis=1)
    check("profile fillet places arc points at the requested radius",
          abs(got_r - 2.0) < 1e-9 and np.abs(radii - 2.0).max() < 1e-9,
          f"radius {got_r:.6f}, arc deviation {np.abs(radii - 2.0).max():.2e}")

    # --- STEP writer: exact volume round-trip through the file ---
    cube_v = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    cube_f = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
              (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    tet_v = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    tet_f = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
    tmp = tempfile.mkdtemp(prefix="semgrit_step_")
    try:
        results = []
        for name, v, f, exact in (("cube", cube_v, cube_f, 1.0),
                                  ("tet", tet_v, tet_f, 1.0 / 6.0)):
            p = os.path.join(tmp, f"{name}.step")
            w = StepWriter(p, name)
            w.add_faceted_solid(v, f, name.upper())
            w.finalize()
            got = read_step_faceted_solids(p)
            audit = check_step_solids(p)
            results.append(
                len(got) == 1
                and abs(got[0]["volume"] - exact) / exact < 1e-9
                and audit["ok"]
            )
        check("STEP write/read round-trip reproduces exact volumes", all(results),
              "unit cube (1.0) and tetrahedron (1/6), both closed and outward-oriented")

        # An inward-wound solid must be caught, so the audit is not vacuous.
        p = os.path.join(tmp, "bad.step")
        w = StepWriter(p, "bad")
        w.add_faceted_solid(cube_v, [f[::-1] for f in cube_f], "INSIDEOUT")
        w.finalize()
        bad_audit = check_step_solids(p)
        check("STEP audit rejects an inward-wound solid", not bad_audit["ok"],
              "; ".join(bad_audit["issues"]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # --- hex boundary extraction closes ---
    sp = WheelSpec(diameter_mm=50, width_mm=5, sector_deg=40, rim_depth_mm=1,
                   radial_divisions=2, axial_divisions=3)
    nodes, hxs, _ = build_rim_mesh(sp)
    quads = surface_from_hexes(hxs)
    vol_surf = step_signed_volume(nodes, [list(q) for q in quads])
    vol_hex = float(hex_volumes(nodes, hxs).sum())
    edges: dict = {}
    for q in quads:
        for k in range(4):
            e = (int(q[k]), int(q[(k + 1) % 4]))
            edges[e] = edges.get(e, 0) + 1
    closed = all(c == 1 and edges.get((e[1], e[0]), 0) == 1 for e, c in edges.items())
    check("hex-mesh boundary surface is closed and matches the hex volume",
          closed and abs(vol_surf - vol_hex) / vol_hex < 1e-12,
          f"{len(quads)} boundary quads, surface volume {vol_surf:.6f} vs "
          f"hex volume {vol_hex:.6f} mm3")

    # --- rim mesh volume vs analytic ---
    allok = True
    details = []
    for sec in (360.0, 180.0, 90.0, 30.0, 25.0, 7.5):
        sp = WheelSpec(diameter_mm=100, width_mm=10, sector_deg=sec, rim_depth_mm=3,
                       radial_divisions=3, axial_divisions=4,
                       circumferential_divisions_per_deg=1.0)
        n, hx, sets = build_rim_mesh(sp)
        v = hex_volumes(n, hx)
        nt = sp.circumferential_divisions()
        dth = sp.sector_rad / nt
        exact = 0.5 * (sp.outer_radius_mm ** 2 - sp.inner_radius_mm ** 2) * math.sin(dth) * nt * sp.width_mm
        if not (v > 0).all() or abs(v.sum() - exact) / exact > 1e-10:
            allok = False
        has_cuts = ("SECTOR_FACE_START" in sets)
        if has_cuts == sp.is_full_circle:
            allok = False
        details.append(f"{sec:g}deg")
    check("rim hex mesh: positive volumes, matches analytic, correct cut-face sets",
          allok, "sectors " + ", ".join(details))


# --------------------------------------------------------------------------
def integration(images: list[str]) -> tuple[list, list]:
    section(f"2. INTEGRATION ({len(images)} SEM images)")
    sems, all_solids = [], []
    seg_cache = []
    for path in images:
        sem = load_sem_image(path)
        seg = segment_grains(sem)
        grains = measure_all(seg, sem)
        solids, reports = build_grain_library(grains, seg, sem)
        bad = [r for r in reports if not r.get("ok")]
        sems.append(sem)
        all_solids.extend(solids)
        seg_cache.append((sem, seg, grains))
        st = grain_statistics(grains, sem)
        print(f"    {os.path.basename(path):16s} {sem.pixel_size_um * 1000:6.2f} nm/px  "
              f"{seg.n_grains:4d} grains ({st['n_grains_used']:3d} interior)  "
              f"d50={st['equivalent_diameter_um']['d50']:5.2f}um  "
              f"{len(solids):3d} solids, {len(bad)} invalid")

    # --- edge-radius blunting, on the first image ---
    sem0, seg0, gr0 = seg_cache[0]
    sharp, sharp_rep = build_grain_library(
        gr0, seg0, sem0, profile=LoftProfile(edge_radius_um=0.0), max_vertices=64
    )
    r_req = 0.35
    blunt, blunt_rep = build_grain_library(
        gr0, seg0, sem0,
        profile=LoftProfile(edge_radius_um=r_req, arc_segments=1),
        max_vertices=12,
    )
    if sharp and blunt:
        def sharp_count(sset, thr=60.0):
            d = [convex_dihedral_angles(s.vertices, s.faces) for s in sset]
            return float(np.mean([(x > thr).sum() for x in d if len(x)]))

        s_before, s_after = sharp_count(sharp), sharp_count(blunt)
        got = float(np.mean([s.edge_radius_inplane_um for s in blunt]))
        all_valid = all(r.get("ok") for r in blunt_rep)
        vol_ok = max(r["volume_rel_error"] for r in blunt_rep if r.get("ok")) < 1e-9
        check(
            "edge rounding blunts the geometry and stays valid",
            s_after < 0.4 * s_before and all_valid and vol_ok
            and abs(got - r_req) / r_req < 0.10,
            f"edges >60 deg per grain {s_before:.1f} -> {s_after:.1f}; "
            f"in-plane radius {got:.3f} um vs {r_req} requested; "
            f"{sum(1 for r in blunt_rep if r.get('ok'))}/{len(blunt_rep)} solids valid; "
            f"volume still exact",
        )
        mono = [r for r in blunt_rep if not r.get("ok")
                and any("monotonic" in i for i in r.get("issues", []))]
        check("blunted profiles remain monotonic (no overhangs)", not mono,
              f"{len(blunt_rep) - len(mono)}/{len(blunt_rep)} grains monotonic about "
              f"their widest cross-section")

    mres = verify_metrology(sems)
    check("metrology on every image", all(r.ok for r in mres),
          f"{sum(r.ok for r in mres)}/{len(mres)} images calibrated and cross-checked")

    # every image must crop the databar identically for this instrument
    tops = {s.databar_top for s in sems}
    check("databar cropped consistently", len(tops) == 1, f"top row {tops}")

    gres = verify_grain_solids(all_solids)
    check("all grain solids geometrically valid", gres.ok, gres.detail)

    # physical plausibility of the size distribution
    eq = np.array([g.equivalent_diameter_um for sem, seg, gr in seg_cache
                   for g in gr if not g.touches_border])
    check("grain sizes physically plausible (1-30 um)",
          bool(eq.min() > 0.5 and eq.max() < 30.0),
          f"{len(eq)} interior grains, range {eq.min():.2f}-{eq.max():.2f} um, "
          f"d50 {np.median(eq):.2f} um")

    # The original pipeline's calibration bug would show up as ~15-30x smaller.
    check("calibration sanity vs original notebook",
          bool(np.median(eq) > 1.0),
          f"d50 {np.median(eq):.2f} um; the original 1019-px scale bar would have "
          f"given ~{np.median(eq) / 14.9:.3f} um")
    return sems, all_solids


# --------------------------------------------------------------------------
def wheel_tests(solids: list) -> None:
    section("3. WHEEL ASSEMBLY + ABAQUS EXPORT")
    tmp = tempfile.mkdtemp(prefix="semgrit_verify_")
    try:
        for sector, elem, parts in [
            (360.0, "R3D3", "shared"),
            (180.0, "C3D4", "shared"),
            (30.0, "R3D3", "shared"),
            (25.0, "C3D4", "shared"),
            (10.0, "C3D4", "baked"),
        ]:
            spec = WheelSpec(diameter_mm=50, width_mm=5, sector_deg=sector, rim_depth_mm=1.0)
            pop = GrainPopulationSpec(areal_density_per_mm2=12, seed=4242)
            model = build_wheel(spec, solids, pop)
            wres = verify_wheel(model)
            path = os.path.join(tmp, f"w{sector:g}_{elem}_{parts}.inp")
            write_inp(path, model, AbaqusExportOptions(grain_element=elem, grain_parts=parts))
            ires = verify_inp_roundtrip(path, model)
            allr = wres + ires
            ok = all(r.ok for r in allr)
            size = os.path.getsize(path) / 1e6
            check(f"sector {sector:g} deg / {elem} / {parts}", ok,
                  f"{model.achieved_grains} grains, {size:.2f} MB, "
                  f"{sum(r.ok for r in allr)}/{len(allr)} checks"
                  + ("" if ok else "  FAILED: " + "; ".join(r.name for r in allr if not r.ok)))

        # --- STEP export of a real wheel: geometry must survive the file ---
        spec = WheelSpec(diameter_mm=100, width_mm=10, sector_deg=3, rim_depth_mm=1)
        pop = GrainPopulationSpec(areal_density_per_mm2=40, seed=20260728)
        model = build_wheel(spec, solids, pop)
        spath = os.path.join(tmp, "cad.step")
        sinfo = write_wheel_step(spath, model, StepExportOptions(max_grains=60))
        audit = check_step_solids(spath)
        got = read_step_faceted_solids(spath)

        rim = [s for s in got if s["name"] == "BOND_RIM"]
        fe_vol = float(hex_volumes(model.body_nodes, model.body_hexes).sum())
        rim_err = abs(rim[0]["volume"] - fe_vol) / fe_vol if rim else float("inf")

        by_id = {p.placement_id: p for p in model.placements}
        v_err, r_err = [], []
        for s in got:
            if not s["name"].startswith("GRAIN_"):
                continue
            p = by_id[int(s["name"].split("_")[1])]
            src = model.shapes[p.shape_index]
            expect = src.mesh_volume_um3 / 1e9
            v_err.append(abs(s["volume"] - expect) / expect)
            radial = float(np.hypot(s["vertices"][:, 0], s["vertices"][:, 1]).max())
            r_err.append(abs(radial - (spec.outer_radius_mm + p.protrusion_mm)))

        ok = (audit["ok"] and rim_err < 1e-9
              and max(v_err, default=1) < 1e-5 and max(r_err, default=1) < 1e-6)
        check("STEP of a real wheel matches the FE model geometry", ok,
              f"{audit['n_solids']} solids / {audit['n_faces']} faces; "
              f"rim volume err {rim_err:.2e}, grain volume err "
              f"{max(v_err, default=0):.2e}, protrusion err "
              f"{max(r_err, default=0):.2e} mm")

        # --- STL fallback ---
        tris = wheel_triangles(model, max_grains=100)
        tinfo = write_binary_stl(os.path.join(tmp, "cad.stl"), tris)
        check("binary STL export", tinfo["n_triangles"] > 0,
              f"{tinfo['n_triangles']} triangles, {tinfo['size_bytes'] / 1e6:.2f} MB")

        # --- saturated packing: exercises the overlap-rejection path, which the
        # --- sparse cases above never reach ---
        biggest = 2 * max(s.bounding_radius_um for s in solids) / 1000.0
        spec = WheelSpec(diameter_mm=20, width_mm=1, sector_deg=20, rim_depth_mm=0.5)
        # ask for a density whose cell size is below the largest grain diameter
        dense = 1.5 / (biggest ** 2)
        pop = GrainPopulationSpec(areal_density_per_mm2=dense, seed=99)
        model = build_wheel(spec, solids, pop)
        wres = verify_wheel(model)
        ov = [r for r in wres if r.name == "grain_overlap"][0]
        saturated = model.achieved_grains < model.requested_grains
        check("saturated packing rejects overlaps instead of interpenetrating",
              all(r.ok for r in wres) and saturated,
              f"requested {model.requested_grains}, placed {model.achieved_grains} "
              f"({ov.data.get('n_pairs_checked', 0)} near-neighbour pairs tested, "
              f"{ov.data.get('n_overlapping', 0)} overlapping); "
              f"largest grain {biggest * 1000:.2f} um")

        # --- determinism ---
        def build_and_hash(seed: int) -> str:
            spec = WheelSpec(diameter_mm=50, width_mm=5, sector_deg=20, rim_depth_mm=1.0)
            pop = GrainPopulationSpec(areal_density_per_mm2=12, seed=seed)
            m = build_wheel(spec, solids, pop)
            p = os.path.join(tmp, f"det{seed}.inp")
            write_inp(p, m, AbaqusExportOptions(grain_element="R3D3"))
            return hashlib.sha256(open(p, "rb").read()).hexdigest()

        h1, h1b, h2 = build_and_hash(1), build_and_hash(1), build_and_hash(2)
        check("same seed reproduces a byte-identical deck", h1 == h1b, h1[:16])
        check("different seed changes the deck", h1 != h2, f"{h1[:12]} vs {h2[:12]}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
def main() -> int:
    pattern = sys.argv[1] if len(sys.argv) > 1 else "*.tif"
    images = sorted(glob.glob(pattern))
    print(f"semgrit verification suite -- {len(images)} images matching {pattern!r}")
    unit_tests()
    if images:
        sems, solids = integration(images)
        if solids:
            wheel_tests(solids)
    else:
        print("\n(no images found; skipping integration and wheel stages)")

    section("SUMMARY")
    print(f"  passed: {len(PASS)}")
    print(f"  failed: {len(FAIL)}")
    for f in FAIL:
        print(f"    FAILED: {f}")
    return 0 if not FAIL else 1


if __name__ == "__main__":
    raise SystemExit(main())
