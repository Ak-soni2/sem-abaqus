"""The pipeline as a handful of calls, so a simple path and a full one share one body.

The notebook grew two audiences. Someone who wants to control every knob needs the
seventeen cells and their ninety widgets. Someone who wants a deck needs about seven
numbers. Writing the second as its own cell would mean two implementations of measure →
build → verify → package, and two implementations drift: this project has already been
bitten twice by a preview that disagreed with the build it claimed to predict.

So the orchestration moves here and both paths call it. The advanced cells pass every
widget; :func:`simple_params` fills the rest from one curated profile.

Nothing here decides geometry. Placement, tangency and the deck itself stay in
:mod:`semgrit.build_deck` and :mod:`semgrit.rigid_wheel`.
"""

from __future__ import annotations

import os
import pickle
import shutil
import subprocess
import sys
from typing import Optional, Sequence

import numpy as np

from .build_deck import DeckParams

# Grit population, as a user thinks of it rather than as DeckParams spells it.
GRIT_KINDS = ("a fixed number", "single grain", "concentration", "grains per mm2")

# Named workpiece blocks. The small one is the default the decks were validated with;
# the others are the sizes actually asked for during this project.
WORKPIECE_SIZES = {
    "small  48 x 15 x 6 um": (0.048, 0.015, 0.006),
    "medium  100 x 40 x 20 um": (0.100, 0.040, 0.020),
    "large  200 x 200 x 200 um": (0.200, 0.200, 0.200),
}

# Grain-meshing settings simple mode uses. They match what the Advanced cells offer by
# default, which matters twice over: the two paths then build the *same* grain library,
# and the second one to run gets it from the cache instead of measuring again.
SIMPLE_MEASURE = dict(simplify_um=0.1, max_vertices=64, interior_only=True)

# Everything simple mode does not ask about. These are the values the validated decks
# were built with, not invented defaults.
SIMPLE_PROFILE = dict(
    rim_depth_mm=0.012, width_mm=0.030,
    shell_circumferential_divisions=200, shell_axial_divisions=6,
    shell_radial_divisions=1, bond_density_kg_m3=2700.0,
    inset_grit_band=True,
    protrusion_mean=0.55, protrusion_std=0.12,
    protrusion_min=0.25, protrusion_max=0.85,
    max_tilt_deg=35.0, spacing_factor=1.05, seed=20260731,
    wp_element_size_mm=0.0003,
    wp_material="STONE", wp_density_kg_m3=2650.0,
    wp_youngs_modulus_mpa=50_000.0, wp_poisson_ratio=0.25,
    surface_speed_mm_s=30_000.0, travel_margin_mm=0.006, cores=8,
)


class QuickError(RuntimeError):
    pass


def _cache_key(images, pixel_size_um, seg_params, height_model, profile,
               simplify_um, max_vertices, interior_only) -> str:
    """Identifies a measurement run by its inputs, so a repeat can be skipped."""
    import dataclasses
    import hashlib

    def described(o):
        try:
            return repr(sorted(dataclasses.asdict(o).items()))
        except TypeError:
            return repr(o)

    bits = []
    for p in images:
        try:
            st = os.stat(p)
            bits.append("%s|%d|%d" % (os.path.basename(p), st.st_size,
                                      int(st.st_mtime)))
        except OSError:
            bits.append(p)
    bits += [repr(pixel_size_um), described(seg_params), described(height_model),
             described(profile), repr(simplify_um), repr(max_vertices),
             repr(bool(interior_only))]
    return hashlib.sha1("\n".join(bits).encode("utf-8")).hexdigest()


def _restage(images, pixel_size_um, seg_params, solids, grains, log):
    """Rebuild the per-image stage records for a cache hit.

    Only the two cheap steps are repeated -- load the image, segment it -- and the
    grains and solids are taken from the cached library rather than rebuilt, so
    the records the figures read are the ones the cache already describes. The
    per-grain validation reports are not cached, so they come back empty and the
    verification figure says so rather than inventing them.
    """
    import time

    from .measure import measure_all
    from .metrology import load_sem_image
    from .segment import segment_grains

    t0 = time.time()
    out = []
    by_image = {}
    for s in solids:
        by_image.setdefault(s.source_image, []).append(s)
    for path in images:
        base = os.path.splitext(os.path.basename(path))[0]
        try:
            sem = load_sem_image(path, pixel_size_um=(pixel_size_um or None))
            stages = {}
            seg = segment_grains(sem, seg_params, stages=stages)
            g = [x for x in grains if x.source_image == path] or measure_all(seg, sem)
            out.append({"path": path, "name": base, "sem": sem, "seg": seg,
                        "stages": stages, "grains": g,
                        "solids": by_image.get(path, []), "reports": []})
        except Exception as exc:                              # noqa: BLE001
            log("%-22s could not re-derive stages: %s" % (base, exc))
    if out:
        log("%-22s stages re-derived in %.1f s" % ("", time.time() - t0))
    return out


def measure_images(images: Sequence[str], outdir: str, *,
                   pixel_size_um: float = 0.0,
                   seg_params=None, height_model=None, profile=None,
                   simplify_um: float = 0.0, max_vertices: int = 0,
                   interior_only: bool = True, cache: bool = True,
                   keep_stages: bool = False, log=print) -> dict:
    """SEM images -> a library of meshed grain solids, with the report printed.

    Returns ``{"solids", "grains", "images", "cached", "per_image"}`` and writes
    the per-image CSVs and ``grain_library.pkl`` into ``outdir``.

    ``keep_stages=True`` additionally returns, per image, the ``SemImage``, the
    ``Segmentation``, every segmentation intermediate (see
    :data:`semgrit.segment.STAGE_KEYS`) and the per-grain validation reports --
    everything needed to *show* what the pipeline did to the picture rather than
    just report how many grains came out. These are megabytes of arrays, so they
    are never written to the cache; asking for them forces a fresh measurement.

    ``pixel_size_um = 0`` reads the scale from the image metadata, which is the only
    trustworthy source: the drawn scale bar is rounded for display, and getting the
    scale wrong scales every grain -- and so the whole wheel -- with it.

    Segmentation is the slowest step in the notebook, and changing a *wheel* setting
    and re-running should not pay for it again. ``cache`` reuses the stored library
    when the images and every measurement setting are unchanged; touch any of them and
    it re-measures. Set it to False to force a fresh run.
    """
    import time

    from .grain3d import HeightModel, LoftProfile, build_grain_library
    from .measure import measure_all, write_csv
    from .metrology import load_sem_image
    from .segment import SegmentationParams, segment_grains

    if seg_params is None:
        seg_params = SegmentationParams()
    if height_model is None:
        height_model = HeightModel()
    if profile is None:
        profile = LoftProfile()

    os.makedirs(outdir, exist_ok=True)
    pkl = os.path.join(outdir, "grain_library.pkl")
    key = _cache_key(images, pixel_size_um, seg_params, height_model, profile,
                     simplify_um, max_vertices, interior_only)
    if cache and os.path.exists(pkl):
        try:
            with open(pkl, "rb") as fh:
                held = pickle.load(fh)
            if held.get("key") == key and held.get("solids"):
                log("grain library : reusing %d solids measured earlier "
                    "(same images, same settings)" % len(held["solids"]))
                log("                delete %s or set cache=False to re-measure"
                    % os.path.basename(pkl))
                per = []
                if keep_stages:
                    # The stages are tens of megabytes of arrays and are never
                    # pickled. Re-deriving them is only the image load plus the
                    # segmentation -- under a second an image -- while the slow
                    # part of a measurement, lofting and validating every grain
                    # into a solid, still comes from the cache. So asking for the
                    # figures does not throw the library away.
                    log("                re-deriving the segmentation stages for "
                        "the figures (the solids stay cached)")
                    per = _restage(images, pixel_size_um, seg_params,
                                   held["solids"], held.get("grains", []), log)
                return {"solids": held["solids"], "grains": held.get("grains", []),
                        "images": list(images), "cached": True, "per_image": per}
        except Exception:                      # a stale or truncated cache is not fatal
            pass

    t0 = time.time()
    solids, grains, per_image = [], [], []
    for path in images:
        sem = load_sem_image(path, pixel_size_um=(pixel_size_um or None))
        stages = {} if keep_stages else None
        seg = segment_grains(sem, seg_params, stages=stages)
        g = measure_all(seg, sem)
        s, reports = build_grain_library(
            g, seg, sem, height_model=height_model, profile=profile,
            simplify_um=simplify_um, max_vertices=max_vertices,
            interior_only=interior_only)
        base = os.path.splitext(os.path.basename(path))[0]
        write_csv(g, os.path.join(outdir, base + "_grains.csv"))
        bad = [r for r in reports if not r.get("ok", True)]
        if bad:
            import collections as _c
            for msg, n in _c.Counter(
                    "; ".join(r.get("issues", ["?"])) for r in bad).most_common(3):
                log("%-22s %4d rejected: %s" % ("", n, msg[:100]))
        log("%-22s %8.5f um/px (%s)  %4d grains -> %4d solids%s"
            % (base, sem.pixel_size_um, sem.pixel_size_source, len(g), len(s),
               "  [%d rejected]" % len(bad) if bad else ""))
        if sem.scalebar_agreement is not None:
            # a *relative difference*: 0 is perfect. Printing (a - 1) once turned a
            # healthy +0.4% into an alarming -99.6%.
            log("%-22s scale-bar cross-check: %+.2f%%"
                % ("", 100 * sem.scalebar_agreement))
        for w in sem.warnings:
            log("%-22s warning: %s" % ("", w))
        if keep_stages:
            per_image.append({"path": path, "name": base, "sem": sem, "seg": seg,
                              "stages": stages, "grains": g, "solids": s,
                              "reports": reports})
        solids.extend(s)
        grains.extend(g)

    if not solids:
        raise QuickError(
            "no grain solids were built. The reasons are printed above -- read them "
            "before changing anything. If they all say the same thing it is a missing "
            "dependency or a bad outline, NOT a segmentation setting: %d grains were "
            "found perfectly well." % len(grains))

    with open(pkl, "wb") as fh:
        pickle.dump({"solids": solids, "grains": grains, "key": key}, fh)
    log("%-22s measured in %.1f s" % ("", time.time() - t0))
    return {"solids": solids, "grains": grains, "images": list(images),
            "cached": False, "per_image": per_image}


def library_summary(solids: Sequence, log=print) -> dict:
    """The size distribution of the measured grains, printed and returned."""
    d = np.sort([s.height_um for s in solids])
    w = np.sort([max(s.extent_um()[:2]) for s in solids])
    out = {"n": len(solids),
           "height_um": [float(np.percentile(d, q)) for q in (10, 50, 90)],
           "width_um": [float(np.percentile(w, q)) for q in (10, 50, 90)],
           "height_max_um": float(d.max()), "width_max_um": float(w.max()),
           "faces": [int(min(len(s.faces) for s in solids)),
                     int(max(len(s.faces) for s in solids))]}
    log("grain library : %d solids" % out["n"])
    log("  height  um  : d10 %.2f  d50 %.2f  d90 %.2f  max %.2f"
        % (*out["height_um"], out["height_max_um"]))
    log("  width   um  : d10 %.2f  d50 %.2f  d90 %.2f  max %.2f"
        % (*out["width_um"], out["width_max_um"]))
    log("  facets      : %d..%d per grain" % tuple(out["faces"]))
    return out


def simple_params(*, diameter_mm: float, slice_mm: float, grit_kind: str,
                  grit_value: float, workpiece_mm, wp_position: str,
                  standoff_um: float, run_ready: bool, cae_deck: bool,
                  cad: bool, name: str = "wheel",
                  analysis=None) -> DeckParams:
    """The seven simple choices, mapped onto the full parameter set.

    Everything not named here comes from :data:`SIMPLE_PROFILE`, which is the
    configuration the two Abaqus-validated decks were built with.
    """
    from .analysis import AnalysisParams
    from .rigid_wheel import WP_POSITIONS

    if grit_kind not in GRIT_KINDS:
        raise QuickError("grit_kind must be one of %s, not %r"
                         % (", ".join(GRIT_KINDS), grit_kind))
    if wp_position not in WP_POSITIONS:
        raise QuickError("wp_position must be one of %s, not %r"
                         % (", ".join(WP_POSITIONS), wp_position))
    if diameter_mm <= 0 or slice_mm <= 0:
        raise QuickError("wheel diameter and slice length must both be positive")
    L, W, D = (float(x) for x in workpiece_mm)
    if min(L, W, D) <= 0:
        raise QuickError("every workpiece dimension must be positive: got %s"
                         % (workpiece_mm,))
    if slice_mm < L:
        raise QuickError(
            "a %.4f mm slice is shorter than the %.4f mm workpiece, so the block "
            "would hang off both ends of the wheel. Lengthen the slice."
            % (slice_mm, L))

    kw = dict(SIMPLE_PROFILE)
    kw.update(diameter_mm=float(diameter_mm), sector_mode="arc",
              arc_length_mm=float(slice_mm),
              wp_length_mm=L, wp_width_mm=W, wp_depth_mm=D,
              wp_position=wp_position, clearance_um=float(standoff_um),
              name=name)
    if grit_kind == "single grain":
        kw.update(grit_mode="single", single_grain_index=-1)
    elif grit_kind == "a fixed number":
        kw.update(grit_mode="count", grit_count=max(1, int(grit_value)))
    elif grit_kind == "concentration":
        kw.update(grit_mode="concentration", concentration=float(grit_value))
    else:
        kw.update(grit_mode="areal_density", areal_density_per_mm2=float(grit_value))

    if run_ready:
        # depth_of_cut_um = 0 asks the build to close the standoff and then cut 85% of
        # the grain protrusion, which is the only choice that is safe for any wheel.
        kw["analysis"] = analysis or AnalysisParams(enabled=True, depth_of_cut_um=0.0)
        kw["also_write_cae_deck"] = bool(cae_deck)
    kw.update(write_step=bool(cad), write_stl=bool(cad))
    return DeckParams(**kw)


def verify_decks(work: str, decks: Sequence[str], log=print) -> bool:
    """Run both independent deck verifiers over every deck. True if all pass."""
    ok = True
    for inp in decks:
        log("#" * 78)
        log("# %s" % os.path.basename(inp))
        log("#" * 78)
        for v in ("verify_rigid_deck.py", "verify_rigid_deck2.py"):
            script = os.path.join(work, v)
            if not os.path.exists(script):
                log("  %s not found next to the notebook payload - skipped" % v)
                ok = False
                continue
            r = subprocess.run([sys.executable, script, inp],
                               capture_output=True, text=True)
            log(r.stdout)
            if r.stderr.strip():
                log(r.stderr)
            ok = ok and r.returncode == 0
    log("=" * 78)
    log(("ALL %d DECK(S) GOOD" % len(decks)) if ok else
        "DECK FAILED VERIFICATION - do not run it; read the FAIL lines above")
    log("=" * 78)
    return ok


def bundle(work: str, dirs: Sequence[str], name: str, log=print) -> str:
    """Collect the outputs into one zip and return its path."""
    out = os.path.join(work, name + "_bundle")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    import glob as _g

    deck_dir = next((d for d in dirs if os.path.basename(d) == "2_abaqus"), None)
    if deck_dir and os.path.isdir(deck_dir):
        # the viewer's glTF files live beside the notebook, not in the deck folder
        for f in _g.glob(os.path.join(work, "*.glb")):
            shutil.copy(f, os.path.join(deck_dir, os.path.basename(f)))
    for src in dirs:
        if os.path.isdir(src):
            shutil.copytree(src, os.path.join(out, os.path.basename(src)))
    zip_path = shutil.make_archive(out, "zip", out)
    log("%.2f MB -> %s" % (os.path.getsize(zip_path) / 1e6, zip_path))
    return zip_path
