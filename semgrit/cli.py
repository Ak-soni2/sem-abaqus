"""Command line interface.

Two subcommands:

``analyze``
    Calibrate, segment and measure SEM images; write per-grain tables, overlays
    and a reusable grain library.

``wheel``
    Assemble measured grains onto a wheel or sector and export an Abaqus deck.

Both are batch-capable and fully seeded, replacing the notebook's interactive
``files.upload()`` / ``input()`` prompts which could only ever handle one image
and could not be reproduced.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import sys
import time
from dataclasses import asdict
from typing import Optional, Sequence

import numpy as np

from .abaqus import (
    MATERIALS,
    AbaqusExportOptions,
    write_cae_import_script,
    write_grain_stl,
    write_inp,
    write_placement_csv,
)
from .grain3d import HeightModel, LoftProfile, build_grain_library
from .measure import grain_statistics, measure_all, write_csv, write_json
from .metrology import MetrologyError, load_sem_image
from .segment import SegmentationParams, segment_grains
from .verify import (
    summarise,
    verify_grain_solids,
    verify_inp_roundtrip,
    verify_metrology,
    verify_wheel,
)
from .wheel import GrainPopulationSpec, WheelSpec, build_wheel

LIBRARY_NAME = "grain_library.pkl"


def _expand(patterns: Sequence[str]) -> list[str]:
    files: list[str] = []
    for p in patterns:
        hits = sorted(glob.glob(p))
        if hits:
            files.extend(hits)
        elif os.path.exists(p):
            files.append(p)
        else:
            print(f"warning: no files match {p!r}", file=sys.stderr)
    seen, out = set(), []
    for f in files:
        if f not in seen:
            seen.add(f)
            out.append(f)
    return out


def _save_overlay(path: str, sem, seg) -> None:
    import cv2
    from skimage.color import label2rgb

    ov = (label2rgb(seg.labels, image=sem.intensity, bg_label=0, alpha=0.35) * 255).astype(
        np.uint8
    )
    for lid in seg.label_ids:
        m = (seg.labels == lid).astype(np.uint8)
        cs, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        colour = (255, 60, 60) if lid in seg.border_labels else (0, 255, 0)
        cv2.drawContours(ov, cs, -1, colour, 1)
    cv2.imwrite(path, cv2.cvtColor(ov, cv2.COLOR_RGB2BGR))


def cmd_analyze(args: argparse.Namespace) -> int:
    files = _expand(args.images)
    if not files:
        print("error: no input images", file=sys.stderr)
        return 2
    os.makedirs(args.out, exist_ok=True)

    seg_params = SegmentationParams(
        min_grain_um=args.min_grain_um,
        h_maxima_um=args.h_maxima_um,
        gradient_weight=args.gradient_weight,
        min_edge_strength=args.min_edge_strength,
        min_area_um2=args.min_area_um2,
        threshold_method=args.threshold,
    )
    height_model = HeightModel(
        mean_ratio=args.thickness_ratio,
        std_ratio=args.thickness_std,
        seed=args.seed,
    )
    profile = LoftProfile(
        base_scale=args.base_scale,
        top_scale=args.top_scale,
        mid_height_fraction=args.mid_height,
        edge_radius_um=args.edge_radius_um,
        arc_segments=args.arc_segments,
    )

    all_grains, all_solids, all_sems = [], [], []
    per_image: list[dict] = []
    t0 = time.time()

    for path in files:
        try:
            sem = load_sem_image(path, pixel_size_um=args.pixel_size_um)
        except MetrologyError as exc:
            print(f"FAIL {path}: {exc}", file=sys.stderr)
            per_image.append({"image": path, "error": str(exc)})
            continue

        seg = segment_grains(sem, seg_params)
        grains = measure_all(seg, sem)
        stats = grain_statistics(grains, sem, interior_only=not args.include_border)
        solids, reports = build_grain_library(
            grains,
            seg,
            sem,
            height_model=height_model,
            profile=profile,
            simplify_um=args.simplify_um,
            max_vertices=args.max_vertices,
            interior_only=not args.include_border,
        )

        base = os.path.splitext(os.path.basename(path))[0]
        write_csv(grains, os.path.join(args.out, f"{base}_grains.csv"))
        write_json(
            {
                "statistics": stats,
                "calibration": {
                    "pixel_size_um": sem.pixel_size_um,
                    "source": sem.pixel_size_source,
                    "scalebar_length_px": sem.scale_bar.length_px if sem.scale_bar else None,
                    "scalebar_label_um": sem.scale_bar.snapped_label_um if sem.scale_bar else None,
                    "scalebar_agreement": sem.scalebar_agreement,
                    "magnification": sem.magnification,
                    "stage_tilt_deg": sem.stage_tilt_deg,
                    "databar_top_row": sem.databar_top,
                    "warnings": sem.warnings,
                },
                "segmentation": {
                    "n_seeds": seg.n_seeds,
                    "n_grains": seg.n_grains,
                    "rejected": seg.rejected,
                    "thresholds": list(seg.threshold_values),
                    "params": asdict(seg_params),
                },
                "solids": {
                    "n_built": len(solids),
                    "n_valid": sum(1 for r in reports if r.get("ok")),
                    "failures": [r for r in reports if not r.get("ok")],
                },
            },
            os.path.join(args.out, f"{base}_report.json"),
        )
        if not args.no_overlay:
            _save_overlay(os.path.join(args.out, f"{base}_segmentation.png"), sem, seg)

        n_valid = sum(1 for r in reports if r.get("ok"))
        print(
            f"{path}: {sem.pixel_size_um * 1000:.2f} nm/px "
            f"({sem.pixel_size_source}, bar "
            f"{'n/a' if sem.scalebar_agreement is None else f'{100 * sem.scalebar_agreement:+.2f}%'}) "
            f"| {seg.n_grains} grains ({stats['n_grains_used']} interior) "
            f"| d50 {stats['equivalent_diameter_um']['d50']:.2f} um "
            f"| {len(solids)} solids, {n_valid} valid"
        )
        if solids:
            edge = _edge_summary(solids)
            print(
                f"    geometry: {edge['faces']:.0f} faces / {edge['tets']:.0f} tets per grain"
                + (
                    f" | edge radius {args.edge_radius_um:.3f} um requested -> "
                    f"{edge['inplane']:.3f} in-plane / {edge['meridional']:.3f} "
                    f"circumferential achieved"
                    if args.edge_radius_um > 0
                    else " | edges left sharp (--edge-radius-um 0)"
                )
            )
            print(
                f"    sharpness: {edge['gt60']:.1f} edges >60 deg, {edge['gt30']:.1f} >30 deg "
                f"per grain (max {edge['max']:.0f} deg); blunting removed "
                f"{100 * edge['loss']:.1f}% of silhouette area"
            )
            print(
                f"    tet quality (C3D4 only): min volume {edge['min_tet_volume_um3']:.2e} um3, "
                f"median quality {edge['median_tet_quality']:.3f}, "
                f"{100 * edge['frac_slivers']:.1f}% slivers "
                f"-- irrelevant for R3D3 rigid grains, which use the surface only"
            )
        for w in sem.warnings:
            print(f"    warning: {w}")

        all_grains.extend(grains)
        all_solids.extend(solids)
        all_sems.append(sem)
        per_image.append(
            {"image": path, "statistics": stats, "n_solids": len(solids), "n_valid": n_valid}
        )

    if not all_solids:
        print("error: no valid grain solids were produced", file=sys.stderr)
        return 1

    # Pooled statistics across every image.
    pooled = _pooled_statistics(all_grains, all_sems, not args.include_border)
    write_json(
        {"pooled": pooled, "per_image": per_image, "seed": args.seed},
        os.path.join(args.out, "summary.json"),
    )

    lib_path = os.path.join(args.out, LIBRARY_NAME)
    with open(lib_path, "wb") as fh:
        pickle.dump(
            {
                "solids": all_solids,
                "pooled": pooled,
                "seed": args.seed,
                "height_model": asdict(height_model),
                "loft_profile": asdict(profile),
            },
            fh,
        )

    if args.stl:
        stl_dir = os.path.join(args.out, "stl")
        os.makedirs(stl_dir, exist_ok=True)
        for s in all_solids:
            base = os.path.splitext(os.path.basename(s.source_image))[0]
            write_grain_stl(os.path.join(stl_dir, f"{base}_grain{s.grain_id}.stl"), s)
        print(f"wrote {len(all_solids)} per-grain STL files -> {stl_dir}")

    if args.step:
        from .step import check_step_solids, write_grains_step

        step_path = os.path.join(args.out, "grains.step")
        sinfo = write_grains_step(
            step_path, all_solids, max_grains=args.step_max_grains, laid_out=True
        )
        print(f"wrote {step_path} ({sinfo['size_bytes'] / 1e6:.2f} MB): "
              f"{sinfo['n_grain_bodies']} grain solid bodies laid out on a grid, "
              f"in microns -- open in SOLIDWORKS via File > Open")
        for w in sinfo["warnings"]:
            print(f"  warning: {w}")
        audit = check_step_solids(step_path)
        print(f"  [{'PASS' if audit['ok'] else 'FAIL'}] STEP audit: "
              f"{audit['n_solids']} solids, {audit['n_faces']} faces"
              + ("" if audit["ok"] else f" -- {'; '.join(audit['issues'])}"))

    print(
        f"\n{len(all_solids)} grain solids from {len(all_sems)} images in "
        f"{time.time() - t0:.1f}s"
    )
    print(f"pooled d50 = {pooled['equivalent_diameter_um']['d50']:.3f} um, "
          f"areal density = {pooled['areal_density_per_mm2']:.0f} /mm2")
    print(f"grain library -> {lib_path}")

    if args.verify:
        print("\n--- verification ---")
        results = verify_metrology(all_sems) + [verify_grain_solids(all_solids)]
        for r in results:
            print(f"  {r}")
        ok, total = summarise(results)
        print(f"  {ok}/{total} checks passed")
        return 0 if ok == total else 1
    return 0


def _edge_summary(solids: list) -> dict:
    """Geometry, sharpness and mesh-quality averages over a grain set."""
    from .grain3d import convex_dihedral_angles, tet_quality

    q = [tet_quality(s.vertices, s.tets) for s in solids]
    dih = [convex_dihedral_angles(s.vertices, s.faces) for s in solids]
    dih = [d for d in dih if len(d)]
    outline = np.array(
        [abs(_signed_area_2d(s.outline_um)) for s in solids], dtype=float
    )
    widest = np.array([_widest(s) for s in solids], dtype=float)
    loss = np.mean((outline - widest) / np.where(outline > 0, outline, np.nan))
    return {
        "faces": float(np.mean([len(s.faces) for s in solids])),
        "tets": float(np.mean([s.n_tets for s in solids])),
        "inplane": float(np.mean([s.edge_radius_inplane_um for s in solids])),
        "meridional": float(np.mean([s.edge_radius_meridional_um for s in solids])),
        "max": float(max(d.max() for d in dih)) if dih else float("nan"),
        "gt60": float(np.mean([(d > 60).sum() for d in dih])) if dih else float("nan"),
        "gt30": float(np.mean([(d > 30).sum() for d in dih])) if dih else float("nan"),
        "loss": float(loss) if np.isfinite(loss) else 0.0,
        "min_tet_volume_um3": float(min(x["min_volume_um3"] for x in q)) if q else float("nan"),
        "min_tet_quality": float(min(x["min_quality"] for x in q)) if q else float("nan"),
        "median_tet_quality": float(np.mean([x["median_quality"] for x in q])) if q else float("nan"),
        "frac_slivers": float(np.mean([x["frac_quality_below_0p01"] for x in q])) if q else float("nan"),
    }


def _signed_area_2d(ring: np.ndarray) -> float:
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _widest(solid) -> float:
    zs = np.unique(np.round(solid.vertices[:, 2], 9))
    best = 0.0
    for z in zs:
        pts = solid.vertices[np.isclose(solid.vertices[:, 2], z)][:, :2]
        if len(pts) >= 3:
            best = max(best, abs(_signed_area_2d(pts)))
    return best


def _pooled_statistics(grains: list, sems: list, interior_only: bool) -> dict:
    used = [g for g in grains if not (interior_only and g.touches_border)]
    total_area = sum(s.field_area_um2 for s in sems)

    def stat(vals: list[float]) -> dict:
        v = np.array([x for x in vals if np.isfinite(x)], dtype=float)
        if v.size == 0:
            return {"n": 0}
        return {
            "n": int(v.size),
            "mean": float(v.mean()),
            "std": float(v.std(ddof=1)) if v.size > 1 else 0.0,
            "min": float(v.min()),
            "max": float(v.max()),
            "d10": float(np.percentile(v, 10)),
            "d50": float(np.percentile(v, 50)),
            "d90": float(np.percentile(v, 90)),
        }

    return {
        "n_images": len(sems),
        "n_grains_total": len(grains),
        "n_grains_used": len(used),
        "interior_only": interior_only,
        "total_field_area_um2": total_area,
        "areal_density_per_mm2": len(used) / (total_area / 1e6) if total_area > 0 else float("nan"),
        "equivalent_diameter_um": stat([g.equivalent_diameter_um for g in used]),
        "feret_max_um": stat([g.feret_max_um for g in used]),
        "feret_min_um": stat([g.feret_min_um for g in used]),
        "aspect_ratio": stat([g.aspect_ratio for g in used]),
        "solidity": stat([g.solidity for g in used]),
        "circularity": stat([g.circularity for g in used]),
        "min_corner_angle_deg": stat([g.min_corner_angle_deg for g in used]),
    }


def cmd_wheel(args: argparse.Namespace) -> int:
    lib_path = args.library
    if os.path.isdir(lib_path):
        lib_path = os.path.join(lib_path, LIBRARY_NAME)
    if not os.path.exists(lib_path):
        print(f"error: grain library not found: {lib_path}", file=sys.stderr)
        print("run 'analyze' first", file=sys.stderr)
        return 2
    with open(lib_path, "rb") as fh:
        lib = pickle.load(fh)
    solids = lib["solids"]
    if not solids:
        print("error: grain library is empty", file=sys.stderr)
        return 1

    os.makedirs(args.out, exist_ok=True)

    spec = WheelSpec(
        diameter_mm=args.diameter,
        width_mm=args.width,
        sector_deg=args.sector,
        rim_depth_mm=args.rim_depth,
        hub_diameter_mm=args.hub_diameter,
        radial_divisions=args.radial_divisions,
        axial_divisions=args.axial_divisions,
        circumferential_divisions_per_deg=args.circ_divisions_per_deg,
    )
    population = GrainPopulationSpec(
        areal_density_per_mm2=args.areal_density,
        concentration=args.concentration,
        volume_fraction=args.volume_fraction,
        max_grains=args.max_grains,
        protrusion_mean=args.protrusion,
        protrusion_std=args.protrusion_std,
        max_tilt_deg=args.max_tilt,
        spacing_factor=args.spacing_factor,
        seed=args.seed,
    )

    t0 = time.time()
    model = build_wheel(spec, solids, population)

    print(f"wheel: D={spec.diameter_mm} mm, W={spec.width_mm} mm, "
          f"sector={spec.sector_deg:g} deg, rim={model.stats['rim_depth_mm']:.3f} mm")
    print(f"  surface area          : {spec.surface_area_mm2:.3f} mm2")
    print(f"  areal density         : requested "
          f"{model.stats['areal_density_per_mm2']:.1f}/mm2, achieved "
          f"{model.stats['achieved_areal_density_per_mm2']:.1f}/mm2")
    print(f"  grains                : {model.achieved_grains} "
          f"(from {len({p.shape_index for p in model.placements})} distinct measured shapes)")
    print(f"  bond mesh             : {model.n_body_nodes} nodes, "
          f"{model.n_body_elements} C3D8")
    print(f"  grain elements         : {model.total_grain_tets()} C3D4 "
          f"/ {model.total_grain_faces()} R3D3 if rigid")
    for w in model.warnings:
        print(f"  warning: {w}")

    opts = AbaqusExportOptions(
        grain_element=args.grain_element,
        grain_parts=args.grain_parts,
        grain_material=args.grain_material,
        bond_material=args.bond_material,
        include_body=not args.no_body,
        # geometry-only leaves out the step, the contact and the boundary
        # conditions, so the deck imports with nothing to object to and the
        # analyst adds the workpiece and step themselves.
        include_step_template=not (args.no_step or args.geometry_only),
        include_contact=not args.geometry_only,
    )
    inp_path = os.path.join(args.out, args.name + ".inp")
    summary = write_inp(
        inp_path,
        model,
        opts,
        provenance={
            "grain_source_images": len({s.source_image for s in solids}),
            "grain_shapes_available": len(solids),
            "analysis_seed": lib.get("seed"),
            "placement_seed": population.seed,
            "height_model": lib.get("height_model"),
        },
    )
    write_placement_csv(os.path.join(args.out, args.name + "_placements.csv"), model)

    # CAE has no dependable menu route for loading an .inp *with* its assembly, so
    # ship the script that does it.
    script_path = os.path.join(args.out, args.name + "_import_into_cae.py")
    write_cae_import_script(script_path, inp_path, model_name=args.name[:38])
    print(f"wrote {script_path}")
    print("  -> in Abaqus/CAE: File > Run Script... and pick that file. It loads the")
    print("     deck WITH the assembly; File > Import > Part would drop every grit.")
    write_json(
        {"wheel": model.stats, "export": summary, "warnings": model.warnings},
        os.path.join(args.out, args.name + "_report.json"),
    )

    size_mb = os.path.getsize(inp_path) / 1e6
    print(f"\nwrote {inp_path} ({size_mb:.2f} MB) in {time.time() - t0:.1f}s")
    print(f"  {summary['n_grain_parts']} grain parts, {summary['n_instances']} instances, "
          f"{summary['n_grain_elements']} grain elements "
          f"({summary['grain_element']}), {summary['n_body_elements']} bond elements")

    # ---- CAD exports -----------------------------------------------------
    if args.step:
        from .step import (
            StepExportOptions,
            check_step_solids,
            estimate_step_size_bytes,
            write_wheel_step,
        )

        step_path = os.path.join(args.out, args.name + ".step")
        n_cad = (
            len(model.placements)
            if args.step_max_grains <= 0
            else min(args.step_max_grains, len(model.placements))
        )
        faces_each = (
            float(np.mean([len(s.faces) for s in model.shapes])) if model.shapes else 0.0
        )
        est = estimate_step_size_bytes(n_cad, faces_each)
        print(f"\nSTEP export: {n_cad} grain bodies at ~{faces_each:.0f} faces each, "
              f"estimated {est / 1e6:.0f} MB")
        if est > 500e6:
            print(f"  NOTE: a file this large will be slow or impossible to open in "
                  f"CAD. It is written because you asked for every grain; the .inp "
                  f"already contains all {len(model.placements)} regardless.")
        t1 = time.time()
        sinfo = write_wheel_step(
            step_path,
            model,
            StepExportOptions(
                max_grains=args.step_max_grains,
                include_body=not args.no_body,
                name=args.name,
            ),
        )
        print(f"\nwrote {step_path} ({sinfo['size_bytes'] / 1e6:.2f} MB) in "
              f"{time.time() - t1:.1f}s  -- open in SOLIDWORKS via File > Open")
        print(f"  {sinfo['n_solids']} solid bodies "
              f"({sinfo['n_grain_bodies']} grains + bond rim), "
              f"{sinfo['n_entities']} STEP entities")
        for w in sinfo["warnings"]:
            print(f"  warning: {w}")
        if args.verify:
            audit = check_step_solids(step_path)
            status = "PASS" if audit["ok"] else "FAIL"
            print(f"  [{status}] STEP audit: {audit['n_solids']} solids, "
                  f"{audit['n_faces']} faces, all closed and outward-oriented"
                  if audit["ok"] else
                  f"  [FAIL] STEP audit: {'; '.join(audit['issues'])}")

    if args.stl:
        from .step import wheel_triangles, write_binary_stl

        stl_path = os.path.join(args.out, args.name + ".stl")
        tris = wheel_triangles(model, max_grains=args.stl_max_grains)
        tinfo = write_binary_stl(stl_path, tris)
        print(f"\nwrote {stl_path} ({tinfo['size_bytes'] / 1e6:.2f} MB, "
              f"{tinfo['n_triangles']} triangles)")
        if args.stl_max_grains and len(model.placements) > args.stl_max_grains:
            print(f"  warning: STL limited to {args.stl_max_grains} of "
                  f"{len(model.placements)} grains")

    if args.verify:
        print("\n--- verification ---")
        results = verify_wheel(model) + verify_inp_roundtrip(inp_path, model)
        for r in results:
            print(f"  {r}")
        ok, total = summarise(results)
        print(f"  {ok}/{total} checks passed")
        return 0 if ok == total else 1
    return 0


def cmd_grind(args: argparse.Namespace) -> int:
    """Build a runnable Abaqus/Explicit multi-grit scratch simulation."""
    import math

    from .grinding_sim import (
        ProcessSpec,
        StoneMaterial,
        WorkpieceSpec,
        write_grinding_sim,
    )
    from .wheel import GrainPopulationSpec, WheelSpec, build_wheel

    lib_path = args.library
    if os.path.isdir(lib_path):
        lib_path = os.path.join(lib_path, LIBRARY_NAME)
    if not os.path.exists(lib_path):
        print(f"error: grain library not found: {lib_path}", file=sys.stderr)
        return 2
    with open(lib_path, "rb") as fh:
        lib = pickle.load(fh)
    solids = lib["solids"]
    os.makedirs(args.out, exist_ok=True)

    wp = WorkpieceSpec(
        length_mm=args.workpiece_length,
        width_mm=args.workpiece_width,
        depth_mm=args.workpiece_depth,
        element_size_mm=args.element_size,
    )
    pr = ProcessSpec(
        wheel_speed_m_s=args.wheel_speed,
        engagement_um=args.engagement_um,
        approach_gap_um=args.approach_gap_um,
        duration_s=args.duration,
    )
    st = StoneMaterial(
        youngs_modulus_mpa=args.stone_E,
        poisson_ratio=args.stone_nu,
        density_kg_m3=args.stone_density,
        compressive_strength_mpa=args.stone_strength,
        failure_strain=args.stone_failure_strain,
    )

    # Grit patch: long enough that grits sweep the whole workpiece, same width.
    arc = args.grit_arc_mm
    spec = WheelSpec(
        diameter_mm=args.diameter,
        width_mm=wp.width_mm,
        sector_deg=math.degrees(arc / (args.diameter / 2.0)),
        rim_depth_mm=args.rim_depth,
        circumferential_divisions_per_deg=max(200.0, 20.0 / max(arc, 1e-9)),
    )
    pop = (
        GrainPopulationSpec(areal_density_per_mm2=args.areal_density, seed=args.seed)
        if args.areal_density
        else GrainPopulationSpec(concentration=args.concentration, seed=args.seed)
    )
    model = build_wheel(spec, solids, pop)
    print(f"grit patch: {arc:.4f} mm arc x {wp.width_mm} mm "
          f"({spec.sector_deg:.4f} deg of a {args.diameter} mm wheel)")
    print(f"  grits: {model.achieved_grains} at "
          f"{model.stats['achieved_areal_density_per_mm2']:.0f}/mm2")
    for w in model.warnings:
        print(f"  warning: {w}")

    path = os.path.join(args.out, args.name + ".inp")
    info = write_grinding_sim(path, model, wp, pr, st, grain_material=args.grain_material)
    print(f"\nwrote {path} ({info['size_bytes'] / 1e6:.2f} MB)")
    print(f"  workpiece      : {wp.length_mm} x {wp.width_mm} x {wp.depth_mm} mm, "
          f"{wp.element_size_mm * 1000:g} um elements -> "
          f"{info['n_workpiece_elements']:,} C3D8R")
    print(f"  grits          : {info['n_grits']} rigid, "
          f"{info['n_grit_elements']:,} R3D3, {info['n_grit_parts']} distinct shapes")
    print(f"  engagement     : {info['engagement_um']:g} um of a "
          f"{info['max_protrusion_um']:.2f} um tallest protrusion")
    print(f"  step time      : {info['duration_s']:.3e} s at "
          f"{info['wheel_speed_m_s']:g} m/s")
    print(f"  stable dt      : {info['dt_stable_s']:.3e} s  -> "
          f"~{info['n_increments']:,.0f} increments")
    write_json(info, os.path.join(args.out, args.name + "_report.json"))

    from .abaqus import write_cae_import_script

    write_cae_import_script(
        os.path.join(args.out, args.name + "_import_into_cae.py"), path,
        model_name=args.name[:38])
    print(f"  CAE script     : {args.name}_import_into_cae.py")
    print("\nJH-2 note: Abaqus/Explicit has no native Johnson-Holmquist keyword. The")
    print("JH-2 constants are written into the deck as comments in VUMAT order; the")
    print("ACTIVE model is Drucker-Prager with shear failure so it runs as-is.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="semgrit",
        description="SEM abrasive-grain measurement and Abaqus grinding-wheel model generation",
    )
    sub = p.add_subparsers(dest="command", required=True)

    a = sub.add_parser("analyze", help="measure grains in SEM images")
    a.add_argument("images", nargs="+", help="image files or globs")
    a.add_argument("-o", "--out", default="results", help="output directory")
    a.add_argument("--pixel-size-um", type=float, default=None,
                   help="override calibration (default: read from SEM metadata)")
    a.add_argument("--threshold", default="multiotsu", choices=["otsu", "multiotsu"])
    a.add_argument("--min-grain-um", type=float, default=0.9)
    a.add_argument("--h-maxima-um", type=float, default=0.12)
    a.add_argument("--gradient-weight", type=float, default=1.0)
    a.add_argument("--min-edge-strength", type=float, default=1.5)
    a.add_argument("--min-area-um2", type=float, default=0.7)
    a.add_argument("--simplify-um", type=float, default=0.10,
                   help="outline simplification tolerance")
    a.add_argument("--max-vertices", type=int, default=64)
    a.add_argument("--thickness-ratio", type=float, default=0.70,
                   help="mean grain height / minimum Feret width")
    a.add_argument("--thickness-std", type=float, default=0.12)
    a.add_argument("--base-scale", type=float, default=0.70)
    a.add_argument("--top-scale", type=float, default=0.30)
    a.add_argument("--mid-height", type=float, default=0.42)
    a.add_argument("--edge-radius-um", type=float, default=0.0,
                   help="cutting edge radius in microns. 0 leaves knife edges, which "
                        "are stress singularities in FEA. A useful starting point is "
                        "~10%% of the measured d50 grain size")
    a.add_argument("--arc-segments", type=int, default=3,
                   help="facets per rounded edge; higher is smoother but adds elements")
    a.add_argument("--include-border", action="store_true",
                   help="include border-truncated grains (biases size stats low)")
    a.add_argument("--stl", action="store_true", help="also write per-grain STL files")
    a.add_argument("--step", action="store_true",
                   help="also write grains.step: individual grain solids laid out on "
                        "a grid, for inspection in SOLIDWORKS")
    a.add_argument("--step-max-grains", type=int, default=200)
    a.add_argument("--no-overlay", action="store_true")
    a.add_argument("--seed", type=int, default=20260728)
    a.add_argument("--verify", action="store_true", help="run verification checks")
    a.set_defaults(func=cmd_analyze)

    w = sub.add_parser("wheel", help="build an Abaqus grinding-wheel model")
    w.add_argument("library", help="results directory or grain_library.pkl from 'analyze'")
    w.add_argument("-o", "--out", default="wheel", help="output directory")
    w.add_argument("--name", default="grinding_wheel")
    w.add_argument("--diameter", type=float, required=True, help="wheel diameter (mm)")
    w.add_argument("--width", type=float, required=True, help="wheel width (mm)")
    w.add_argument("--sector", type=float, default=360.0,
                   help="angular sector in degrees: 360 full, 180 half, 30, 25, ...")
    w.add_argument("--rim-depth", type=float, default=None,
                   help="model only this radial depth (mm); omit for a solid body")
    w.add_argument("--hub-diameter", type=float, default=0.0)
    w.add_argument("--radial-divisions", type=int, default=4)
    w.add_argument("--axial-divisions", type=int, default=8)
    w.add_argument("--circ-divisions-per-deg", type=float, default=1.0)

    g = w.add_mutually_exclusive_group()
    g.add_argument("--areal-density", type=float, default=None,
                   help="grains per mm2 of wheel surface")
    g.add_argument("--concentration", type=float, default=None,
                   help="abrasive concentration number (C100 = 25 vol%%)")
    g.add_argument("--volume-fraction", type=float, default=None)

    w.add_argument("--max-grains", type=int, default=200_000)
    w.add_argument("--protrusion", type=float, default=0.55,
                   help="mean protruding fraction of grain height")
    w.add_argument("--protrusion-std", type=float, default=0.12)
    w.add_argument("--max-tilt", type=float, default=35.0)
    w.add_argument("--spacing-factor", type=float, default=1.05)
    w.add_argument("--grain-element", default="C3D4", choices=["C3D4", "R3D3"])
    w.add_argument("--grain-parts", default="shared", choices=["shared", "baked"])
    w.add_argument("--grain-material", default="diamond", choices=sorted(MATERIALS))
    w.add_argument("--bond-material", default="vitrified_bond", choices=sorted(MATERIALS))
    w.add_argument("--no-body", action="store_true", help="grains only, no bond mesh")
    w.add_argument("--no-step", action="store_true",
                   help="omit the Abaqus template *Step block")
    w.add_argument("--geometry-only", action="store_true",
                   help="wheel geometry, sets and surfaces only: no *Step, no contact, "
                        "no boundary conditions. Imports into CAE with nothing to "
                        "object to; add your own workpiece and step")
    w.add_argument("--step", action="store_true",
                   help="also write a STEP CAD file (opens in SOLIDWORKS as solid bodies)")
    w.add_argument("--step-max-grains", type=int, default=200,
                   help="grain bodies in the STEP file. **0 = every grain**, at ~430 "
                        "bytes per face (a full wheel runs to gigabytes and will not "
                        "open in CAD). Default 200")
    w.add_argument("--stl", action="store_true",
                   help="also write a binary STL (mesh body, opens anywhere)")
    w.add_argument("--stl-max-grains", type=int, default=2000)
    w.add_argument("--seed", type=int, default=20260728)
    w.add_argument("--verify", action="store_true")
    w.set_defaults(func=cmd_wheel)

    g = sub.add_parser("grind",
                       help="build a runnable Abaqus/Explicit grinding simulation")
    g.add_argument("library", help="results directory or grain_library.pkl")
    g.add_argument("-o", "--out", default="grind")
    g.add_argument("--name", default="grind_stone")
    g.add_argument("--diameter", type=float, default=200.0, help="wheel diameter (mm)")
    g.add_argument("--rim-depth", type=float, default=0.02)
    g.add_argument("--grit-arc-mm", type=float, default=0.15,
                   help="length of wheel arc carrying grits")
    g.add_argument("--areal-density", type=float, default=None,
                   help="grits/mm2; omit to use --concentration")
    g.add_argument("--concentration", type=float, default=100.0,
                   help="abrasive concentration number (C100 = 25 vol%%)")
    g.add_argument("--workpiece-length", type=float, default=0.10,
                   help="mm, along the grinding direction")
    g.add_argument("--workpiece-width", type=float, default=0.06, help="mm, axial")
    g.add_argument("--workpiece-depth", type=float, default=0.025, help="mm")
    g.add_argument("--element-size", type=float, default=0.0015,
                   help="workpiece element size (mm). Must be finer than a grit")
    g.add_argument("--wheel-speed", type=float, default=30.0, help="m/s")
    g.add_argument("--engagement-um", type=float, default=2.0,
                   help="how deep the tallest grit tips bite")
    g.add_argument("--approach-gap-um", type=float, default=5.0)
    g.add_argument("--duration", type=float, default=None,
                   help="step time (s); default is one full traverse")
    g.add_argument("--stone-E", type=float, default=50000.0, help="MPa")
    g.add_argument("--stone-nu", type=float, default=0.25)
    g.add_argument("--stone-density", type=float, default=2650.0, help="kg/m3")
    g.add_argument("--stone-strength", type=float, default=150.0, help="MPa")
    g.add_argument("--stone-failure-strain", type=float, default=0.10)
    g.add_argument("--grain-material", default="diamond", choices=sorted(MATERIALS))
    g.add_argument("--seed", type=int, default=20260728)
    g.set_defaults(func=cmd_grind)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "wheel" and (
        args.areal_density is None
        and args.concentration is None
        and args.volume_fraction is None
    ):
        args.concentration = 100.0
        print("note: no grain density given; using concentration C100 (25 vol%)")
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
