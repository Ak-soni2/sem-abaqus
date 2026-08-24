"""Build the three trajectory jobs.

    python _make_trajectory_jobs.py            # design, build and report all three
    python _make_trajectory_jobs.py --plan     # design and report only, write nothing
    python _make_trajectory_jobs.py --job 2    # one job

The three
---------
1. ARC_26UM   one abrasive on the drawn trajectory: depth 0 at entry, 26 um at
              mid-span, back to 0 at exit. This is the SHAPE-verification job.
              26 um is about 490 x SiC's dc, so removal is brittle essentially
              everywhere and no transition is expected. What it proves is that
              the grit really does follow the curve in the figure.

2. ARC_80NM   the same arc, scaled so the maximum depth is 80 nm instead of
              26 um. 80 nm against SiC's dc of 52.92 nm is 1.51 x, so the cut
              starts ductile, crosses into brittle, and crosses back -- the
              transition appears TWICE in one pass, symmetrically about the
              centre. This is the picture the model exists to produce.

3. RAMP_80NM  a straight ramp from 0 to 80 nm over the same length. One crossing
              instead of two, and depth linear in position, so the transition
              station can be read off the figure and checked against dc by hand.
              The control for job 2.

Why the block is 2.28 mm and not the 17 mm on the figure
--------------------------------------------------------
The figure's axes are 17 mm across by 26 um deep -- a vertical exaggeration of
about 650:1. Taken literally, a 17 mm chord with 26 um of sagitta needs a tip
radius of 1389 mm, i.e. a 2.8 m wheel. The wheel in this project is 50 mm
diameter, and over 17 mm of arc a 25 mm radius would cut 1445 um deep, 56 x the
figure.

So the SHAPE is reproduced exactly and the LENGTH follows from the real wheel:
for a sagitta D on radius R the chord is L = sqrt(8 R D), giving 2.28 mm for
26 um and 126.5 um for 80 nm. Both are honest arcs of the actual 25 mm wheel and
neither invents a radius. The figure's 0..17 is a normalised axis.

How the depth is prescribed
---------------------------
``vumat_grind.for`` reads h(u) = H0 + HG*u - u^2/(2*RTIP) from the material
card. That is already a line plus a parabola, so no new physics is needed:

  * an ARC is HG = 0 -- no radial infeed -- with the parabola alone taking h
    from 0 up to the peak and back down. The depth is the wheel's own curvature,
    which is exactly what the figure draws.
  * a RAMP is the opposite: the linear infeed term does the work and the
    parabola is negligible over so short a block.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
R_MM = 25.0                     # wheel radius, the real one
MATERIAL = "silicon_carbide"    # every job, as asked
ELEMENTS_PER_DC = 5.0           # the project standard, from _make_run_packages
RAMP_BOW_FRACTION = 0.025       # job 3's parabola, as a fraction of its peak
PROTRUSION_MEAN_FRACTION = 0.55  # DeckParams protrusion_mean: the share of a
                                 # grain's height that stands proud of the bond

JOBS = {
    1: dict(name="ARC_26UM", peak_um=26.0, kind="arc",
            note="trajectory verification: 26 um sits far above dc, so removal "
                 "is brittle throughout. Proves the grit follows the drawn "
                 "curve; it is not a transition experiment."),
    2: dict(name="ARC_80NM", peak_um=0.080, kind="arc",
            note="the transition job: 80 nm straddles dc = 52.92 nm, so the "
                 "pass runs ductile - brittle - ductile with two crossings."),
    3: dict(name="RAMP_80NM", peak_um=0.080, kind="ramp",
            note="linear control: one crossing, depth linear in position, so "
                 "the transition station is directly checkable against dc."),
}


def arc_chord_mm(peak_mm: float, radius_mm: float = R_MM) -> float:
    """Chord whose sagitta on ``radius_mm`` is exactly ``peak_mm``.

    Sagitta D = L^2/(8R), inverted. This is what makes the drawn shape a real
    arc of the real wheel rather than a curve fitted to the figure's axes.
    """
    return math.sqrt(8.0 * radius_mm * peak_mm)


def scale_solids(solids, factor: float):
    """A geometrically similar copy of the grain library, ``factor`` times bigger.

    Job 1 asks for a 26 um scallop, and the measured B4C grains are about 7 um
    tall. A grit can only cut as deep as it stands proud of the bond, so a 26 um
    cut from a 7 um grain is not a mesh problem or a seating problem -- it is
    physically impossible, and ``verify_rigid_deck`` says so: the bond rim would
    be driven 22 um into the workpiece.

    Scaling is the honest fix and was explicitly permitted. Every length scales
    together -- vertices, outline, height -- so the grain keeps its measured
    SHAPE exactly: the same aspect ratio, the same corner angles, the same
    concave cutting features. Only its size changes, and the report says by how
    much. Nothing about the segmentation or the measurement is re-run or
    re-interpreted; this is the same grain, larger.
    """
    import copy

    out = []
    for g in solids:
        c = copy.deepcopy(g)
        c.vertices = c.vertices * factor
        c.outline_um = c.outline_um * factor
        c.centroid_um = c.centroid_um * factor
        c.height_um = c.height_um * factor
        c.analytic_volume_um3 = c.analytic_volume_um3 * factor ** 3
        for a in ("edge_radius_requested_um", "edge_radius_inplane_um",
                  "edge_radius_meridional_um"):
            if getattr(c, a, 0.0):
                setattr(c, a, getattr(c, a) * factor)
        out.append(c)
    return out


def design(job: int, dc_mm: float) -> dict:
    """Every number for one job, derived rather than typed."""
    spec = JOBS[job]
    peak = spec["peak_um"] / 1000.0

    if job == 3:
        # A ramp has to be STRAIGHT, and h(u) always carries the wheel's
        # parabola. Over the arc's own 126 um chord that parabola is 80 nm --
        # the whole ramp -- so a ramp of that length would bow to 120 nm at
        # mid-span and be no kind of line. RAMP_BOW_FRACTION caps the sagitta
        # over the ramp at a fixed fraction of the peak; 2.5 % is straight to
        # four figures and to the eye. L = sqrt(8 R f D) inverts that.
        length = arc_chord_mm(RAMP_BOW_FRACTION * peak)
    else:
        length = arc_chord_mm(peak)

    if job == 1:
        # dc is irrelevant here -- 26 um is ~490 x dc -- so the mesh is sized
        # for the SHAPE. Resolving a 53 nm dc across a 2.28 mm block would be
        # two million elements spent on a question this job does not ask, and
        # the groove is 26 um deep, so a 3 um element still draws it with nine
        # elements through the cut.
        el_lat, el_depth = 0.004, 0.003
        surface, band = 0.030, 0.0
        width = 0.140          # wider than the scaled grain, which is ~134 um
        depth = 1.35 * peak
    else:
        # The DEPTH element is what resolves dc and what sets the stable
        # increment, so it is fixed at dc/5 and everything else is trimmed
        # around it. Along the scratch the transition is a station, not a
        # gradient, so 0.35 um there costs nothing; across the face only the
        # groove lane needs to be fine, which is what the band is for.
        el_lat = 0.00035
        el_depth = dc_mm / ELEMENTS_PER_DC
        surface = max(4.0 * peak, 0.0002)
        band, width = 0.004, 0.005
        # Five times the cut is ample: a 0.8 um block under an 80 nm cut is
        # already ten times deeper than anything that happens in it.
        depth = max(5.0 * peak, 0.0004)

    return dict(job=job, name=spec["name"], kind=spec["kind"],
                note=spec["note"], peak_um=spec["peak_um"],
                length_mm=length, width_mm=width, depth_mm=depth,
                el_lat_mm=el_lat, el_depth_mm=el_depth,
                surface_mm=surface, band_mm=band,
                arc=spec["kind"] == "arc")


def build(d: dict, solids, write: bool):
    """Plan, and optionally write, one job. Returns (params, plan, info)."""
    from semgrit import materials
    from semgrit.analysis import AnalysisParams
    from semgrit.build_deck import DeckParams, build_deck, plan_deck
    from semgrit.hybrid import HYBRID_DEPVAR

    peak_mm = d["peak_um"] / 1000.0

    # The trajectory is the SPECIFICATION here, so H0 and HG are prescribed
    # rather than left to fall out of the seating. h(u) = H0 + HG*u - u^2/(2R):
    #
    #   arc  : H0 = peak, HG = 0. The parabola alone carries h from 0 at the
    #          entry up to the peak at mid-span and back to 0 at the exit --
    #          which is precisely the drawn curve.
    #   ramp : h must go 0 -> peak across the block, so HG = peak / L with H0
    #          placing the zero at the entry edge. The parabola is still there
    #          (it is the real wheel) but over 126 um it contributes 80 pm, so
    #          the profile is a straight line to four figures.
    if d["arc"]:
        h0, hg = peak_mm, 0.0
    else:
        # Solve for the two card constants that put h = 0 at the entry edge and
        # h = peak at the exit edge WITH the parabola present, rather than
        # pretending it is not there:
        #     h(-hb) = H0 - HG*hb - hb^2/(2R) = 0
        #     h(+hb) = H0 + HG*hb - hb^2/(2R) = peak
        # subtracting gives HG = peak / L, adding gives H0 = peak/2 + hb^2/(2R).
        L = d["length_mm"]
        hb = 0.5 * L
        hg = peak_mm / L
        h0 = 0.5 * peak_mm + hb * hb / (2.0 * R_MM)

    hp = materials.hybrid_params(MATERIAL, h_source=0, dc_form=2,
                                 h0_override_mm=h0, hg_override=hg)

    # depth_of_cut_um still drives the STEP: how far the wheel feeds in over the
    # pass, and therefore how long the step lasts. It no longer sets the profile.
    ae_um = max(d["peak_um"], 0.001)

    p = DeckParams(
        name=d["name"].lower(),
        diameter_mm=2 * R_MM,
        sector_mode="arc",
        arc_length_mm=max(4.0 * d["length_mm"], 0.5),
        # A scaled grain needs a wheel that can carry it: the dressed face must
        # be wider than the grain, and the rim deeper than it protrudes, or the
        # build refuses -- rightly. Both follow the scale rather than being
        # retyped, so the wheel and the grit stay consistent by construction.
        rim_depth_mm=max(0.012, 3.0 * d["peak_um"] / 1000.0),
        width_mm=max(1.5 * d["width_mm"], 0.030,
                     3.0 * d.get("grain_across_mm", 0.0)),
        include_bond=False,            # one grit and the block, nothing else
        grit_mode="single",
        single_grain_index=-1,         # the largest measured grain
        single_grit_offset_mm=0.0,     # centred, so the arc is symmetric
        include_workpiece=True,
        wp_length_mm=d["length_mm"],
        wp_width_mm=d["width_mm"],
        wp_depth_mm=d["depth_mm"],
        wp_element_size_length_mm=d["el_lat_mm"],
        wp_element_size_width_mm=d["el_lat_mm"],
        wp_element_size_depth_mm=d["el_depth_mm"],
        wp_surface_layer_mm=d["surface_mm"],
        wp_depth_growth=1.45,
        wp_width_band_mm=d["band_mm"],
        wp_width_growth=1.35,
        clearance_um=0.0,
        wp_position="centred",
        surface_speed_mm_s=30_000.0,   # 30 m/s at r = 25 mm
        cores=8,
        analysis=AnalysisParams(
            enabled=True, depth_of_cut_um=ae_um,
            material_model="hybrid", hybrid=hp,
            n_depvar=HYBRID_DEPVAR, element_deletion=True),
    )
    materials.apply(p, MATERIAL)
    plan = plan_deck(p, solids)
    info = None
    if write:
        out = os.path.join(HERE, "RUN_TRAJECTORY", d["name"])
        os.makedirs(out, exist_ok=True)
        info = build_deck(p, solids, out)
        for f in ("vumat_grind.for", "vumat_jh2.for"):
            src = os.path.join(HERE, f)
            if os.path.exists(src):
                import shutil
                shutil.copy2(src, os.path.join(out, f))
    return p, plan, info


def h_curve(field, dc_mm, half_mm, n=2001):
    """h(u) across the block, every dc crossing, and the ductile fraction.

    The package's own ``_transition_station`` bisects, which assumes h is
    monotone in u -- true for a ramp, false for an arc, which crosses dc twice.
    A bisection on the arc would silently report one crossing or none, so the
    sampled curve is scanned for every sign change instead and each one is then
    refined. This is the only reason this helper exists rather than reusing the
    package function.
    """
    us = [-half_mm + 2.0 * half_mm * i / (n - 1) for i in range(n)]
    hs = [field.h_at(u) for u in us]
    cross = []
    for i in range(1, n):
        a, b = hs[i - 1] - dc_mm, hs[i] - dc_mm
        if a == 0.0:
            cross.append(us[i - 1])
        elif a * b < 0.0:
            lo, hi = us[i - 1], us[i]
            for _ in range(90):
                mid = 0.5 * (lo + hi)
                if (field.h_at(mid) - dc_mm) * a > 0.0:
                    lo = mid
                else:
                    hi = mid
            cross.append(0.5 * (lo + hi))
    duct = sum(1 for h in hs if h < dc_mm) / float(n)
    return us, hs, cross, duct


def main(argv):
    from semgrit.hybrid import plan_hybrid
    from semgrit.materials import MATERIALS
    from semgrit.quick import measure_images

    write = "--plan" not in argv
    only = int(argv[argv.index("--job") + 1]) if "--job" in argv else None

    dc_mm = MATERIALS[MATERIAL].dc_nm(2) / 1e6
    print("=" * 78)
    print("THREE TRAJECTORY JOBS  -  %s" % MATERIALS[MATERIAL].label)
    print("=" * 78)
    print("  wheel radius        : %.1f mm  (diameter %.0f mm)" % (R_MM, 2 * R_MM))
    print("  dc (Bifano, form 2) : %.3f nm" % (dc_mm * 1e6))
    print("  depth element       : dc / %.0f = %.3f nm"
          % (ELEMENTS_PER_DC, dc_mm / ELEMENTS_PER_DC * 1e6))
    print("  80 nm / dc          : %.2f   -> the transition lands inside the cut"
          % (0.080 / (dc_mm * 1e3)))
    print()

    got = measure_images([os.path.join(HERE, "B4C_15.tif")],
                         os.path.join(HERE, "_traj_meas"), log=lambda *a: None)
    solids = got["solids"]
    tallest = max(s.height_um for s in solids)
    print("  grain library       : %d verified solids from B4C_15.tif"
          % len(solids))
    print("  tallest grain       : %.2f um" % tallest)
    print()

    summary, bad = [], 0
    for job in sorted(JOBS):
        if only and job != only:
            continue
        d = design(job, dc_mm)
        print("-" * 78)
        print("JOB %d  %s   (%s trajectory)" % (job, d["name"], d["kind"]))
        print("-" * 78)
        print("  %s" % d["note"])
        print()
        t0 = time.time()
        # Two passes. RTIP is the GRIT TIP radius -- the ground radius plus the
        # grain's own protrusion -- not the nominal 25 mm, and the chord that
        # closes the arc to zero depends on it. Plan once to learn RTIP, resize
        # the block to that radius, then plan again. Without this the ends of
        # the arc miss zero by about (RTIP-R)/RTIP of the peak: 4 nm on job 1,
        # and proportionally the same fraction on job 2.
        # A grit cuts only as deep as it protrudes. Job 1's 26 um scallop needs
        # a grain that stands 26 um proud of the bond, and the measured ones are
        # 7 um tall, so the library is scaled for that job -- shape preserved,
        # size changed, and the factor reported. Jobs 2 and 3 cut 80 nm, which
        # any measured grain clears by two orders of magnitude, so they use the
        # grains exactly as measured.
        lib = solids
        d["grain_scale"] = 1.0
        need_um = d["peak_um"] / PROTRUSION_MEAN_FRACTION
        if need_um > tallest:
            d["grain_scale"] = round(1.15 * need_um / tallest, 3)
            lib = scale_solids(solids, d["grain_scale"])
            d["grain_across_mm"] = max(float(max(g.extent_um()[:2]))
                                       for g in lib) / 1000.0
            print("  grains scaled x%.3f so a %.1f um protrusion is reachable"
                  % (d["grain_scale"], d["peak_um"]))
            print("  (shape preserved exactly; only the size changes)")
            print()
        p, plan, _unused = build(d, lib, False)
        fld0, _ = plan_hybrid(plan, p.analysis.hybrid)
        if d["arc"]:
            d = dict(d, length_mm=arc_chord_mm(d["peak_um"] / 1000.0,
                                               fld0.rtip_mm))
        p, plan, info = build(d, lib, write)
        fld, dcm = plan_hybrid(plan, p.analysis.hybrid)
        half = d["length_mm"] / 2.0
        us, hs, cross, duct = h_curve(fld, dcm, half)
        h_max = max(hs)

        print("  THE TRAJECTORY")
        print("    block            : %.5f mm long x %.2f um wide x %.4f um deep"
              % (d["length_mm"], d["width_mm"] * 1000, d["depth_mm"] * 1000))
        print("    peak depth wanted: %.4f um" % d["peak_um"])
        print("    h at entry       : %.5f um" % (hs[0] * 1000))
        print("    h at mid-span    : %.5f um" % (hs[len(hs) // 2] * 1000))
        print("    h at exit        : %.5f um" % (hs[-1] * 1000))
        print("    h max            : %.5f um" % (h_max * 1000))
        print("    card constants   : H0 %+.6e  HG %+.6e  RTIP %.4f mm"
              % (fld.h0_mm, fld.hg, fld.rtip_mm))
        print()
        print("  THE TRANSITION")
        print("    dc               : %.3f nm" % (dcm * 1e6))
        print("    dc crossings     : %d%s" % (len(cross), "" if cross else "  (none)"))
        for c in cross:
            print("        u = %+.3f um   (%.1f %% along the pass)"
                  % (c * 1000, 100 * (c + half) / (2 * half)))
        print("    ductile fraction : %.1f %% of the pass" % (100 * duct))
        print()
        print("  THE MESH")
        print("    elements         : %s" % format(plan["n_workpiece_elements"], ","))
        el = plan["element_um"]
        print("    surface element  : %.4f x %.4f x %.5f um  (cut x axial x depth)"
              % (el[0], el[1], el[2]))
        print("    across dc        : %.2f elements" % (dcm * 1000.0 / max(el[2], 1e-12)))
        lo = min(v for v in el[:3] if v > 0)
        print("    aspect ratio     : %.1f : 1" % (max(el[:3]) / lo))
        cost = plan.get("cost") or {}
        print("    stable increment : %.3e s" % (cost.get("stable_dt_s") or 0))
        print("    increments       : %s" % format(int(cost.get("increments") or 0), ","))
        print("    runtime estimate : %.2f h on 8 cores"
              % ((cost.get("est_hours") or {}).get("8", 0.0)))
        if info:
            print()
            print("    WROTE %s  (%.1f MB in %.0f s)"
                  % (os.path.basename(info["path"]),
                     info["size_bytes"] / 1e6, time.time() - t0))
        for w in plan.get("warnings", []):
            print("    warning: %s" % w)

        # --- the assertions that make this a gate, not a report -------------
        checks = []
        if d["arc"]:
            checks.append(("peak depth is the arc sagitta, to 1 %",
                           abs(h_max * 1000 - d["peak_um"]) <= 0.01 * d["peak_um"]))
            checks.append(("depth returns to ~0 at both ends",
                           hs[0] * 1e6 < 1.0 and hs[-1] * 1e6 < 1.0))
            checks.append(("the arc is symmetric about the centre",
                           abs(hs[0] - hs[-1]) * 1e6 < 1.0))
        else:
            checks.append(("depth is monotone along the pass",
                           all(hs[i] >= hs[i - 1] - 1e-12 for i in range(1, len(hs)))
                           or all(hs[i] <= hs[i - 1] + 1e-12
                                  for i in range(1, len(hs)))))
            checks.append(("peak depth reaches the asked-for value, to 5 %",
                           abs(h_max * 1000 - d["peak_um"]) <= 0.05 * d["peak_um"]))
        if job == 1:
            checks.append(("no transition expected: brittle throughout",
                           duct < 0.02))
        else:
            checks.append(("the transition is inside the pass",
                           len(cross) >= 1))
            checks.append(("dc is resolved by >= 4 elements",
                           dcm * 1000.0 / max(el[2], 1e-12) >= 4.0))
            checks.append(("both regimes are actually present",
                           0.05 < duct < 0.95))
        if d["arc"] and job == 2:
            checks.append(("an arc crosses dc twice", len(cross) == 2))
        print()
        print("  CHECKS")
        for label, ok in checks:
            print("    [%s] %s" % ("PASS" if ok else "FAIL", label))
            if not ok:
                bad += 1
        print()

        summary.append(dict(
            job=job, name=d["name"], kind=d["kind"], peak_um=d["peak_um"],
            grain_scale=d.get("grain_scale", 1.0),
            length_mm=d["length_mm"], width_mm=d["width_mm"],
            depth_mm=d["depth_mm"], h_max_um=h_max * 1000,
            h_entry_um=hs[0] * 1000, h_exit_um=hs[-1] * 1000,
            crossings_um=[c * 1000 for c in cross],
            ductile_fraction=duct, dc_nm=dcm * 1e6,
            elements=plan["n_workpiece_elements"],
            element_um=list(el),
            elements_across_dc=dcm * 1000.0 / max(el[2], 1e-12),
            stable_dt_s=cost.get("stable_dt_s"),
            increments=cost.get("increments"),
            est_hours_8=(cost.get("est_hours") or {}).get("8"),
            H0_mm=fld.h0_mm, HG=fld.hg, RTIP_mm=fld.rtip_mm))

    if write and summary:
        out = os.path.join(HERE, "RUN_TRAJECTORY")
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "SUMMARY.json"), "w") as fh:
            json.dump({"material": MATERIAL,
                       "material_label": MATERIALS[MATERIAL].label,
                       "wheel_radius_mm": R_MM, "dc_nm": dc_mm * 1e6,
                       "elements_per_dc": ELEMENTS_PER_DC,
                       "jobs": summary}, fh, indent=1)
        print("wrote RUN_TRAJECTORY/SUMMARY.json")

    print("=" * 78)
    print("ALL CHECKS PASSED" if not bad else "%d CHECK(S) FAILED" % bad)
    print("=" * 78)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
