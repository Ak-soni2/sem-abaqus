"""Verification gate for the swept-envelope chip-thickness engine.

The engine replaces four hard-coded constants with a field computed for every
element, so it has to earn the same trust those constants had. Four independent
things are checked:

**It reproduces the model it replaces.** For one grit the swept field must agree
with the closed-form wedge that ``verify_hybrid_deck.py`` already validates. That
is the whole argument for trusting it at 700 grits: the general method must give
the special case's answer.

**It agrees with brute force.** The sweep is vectorised, chunked and adaptively
sampled. A slow, direct, unchunked evaluation of the same kinematics must reach
the same penetration.

**Shadowing is right.** Synthetic two-grit cases with known geometry, where the
answer can be written down: a grit following a deeper one takes nothing, a grit
following a shallower one takes the difference.

**The deck carries it correctly.** The field is read back out of the written
``.inp`` and compared bit for bit, and the compiled VUMAT is asked which branch
it picks from that card at that field value.

    python verify_envelope.py [-v] [--library path/to/grain_library.pkl]

With no ``--library`` it uses this project's reference library; failing that it
searches for one. The synthetic checks -- surface sampling and shadowing on
geometry whose answer can be written down -- need no library at all and run
either way.
"""

from __future__ import annotations

import math
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "_envchk")


def find_library(argv) -> str:
    """The grain library to build the reference deck from.

    Explicit argument first, then this project's own reference library, then a
    search -- so the gate works in a Colab runtime where the only library is the
    one the notebook just measured.
    """
    for i, a in enumerate(argv):
        if a == "--library" and i + 1 < len(argv):
            return argv[i + 1]
        if a.startswith("--library="):
            return a.split("=", 1)[1]
    ref = os.path.join(HERE, "WHEEL_FIXED", "1_measurements",
                       "grain_library.pkl")
    if os.path.exists(ref):
        return ref
    for root, _d, files in os.walk(os.getcwd()):
        if "grain_library.pkl" in files:
            return os.path.join(root, "grain_library.pkl")
    return ""


LIBRARY = ""

VERBOSE = "-v" in sys.argv or "--verbose" in sys.argv
_PASS = 0
_FAIL: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    global _PASS
    if ok:
        _PASS += 1
        if VERBOSE:
            print("  PASS  %-54s %s" % (name, detail))
    else:
        _FAIL.append(name)
        print("  FAIL  %-54s %s" % (name, detail))


def close(name: str, got: float, want: float, rtol=1e-9, atol=0.0) -> None:
    err = abs(got - want)
    check(name, err <= atol + rtol * abs(want),
          "got %.10g want %.10g (err %.3g)" % (got, want, err))


# --------------------------------------------------------------------------
# synthetic geometry, where the answer can be written down
# --------------------------------------------------------------------------

def _prism(a_top: float, b_c: float, z_c: float, half: float) -> tuple:
    """A tiny axis-aligned box in the block frame, as vertices + triangles.

    ``a_top`` is its outermost radial coordinate, so a box with
    ``a_top = r_ground + d`` reaches depth ``d`` into the block. Simple enough
    that the chip thickness it should produce is arithmetic.
    """
    a0 = a_top - 2 * half
    v = np.array([[a0, b_c - half, z_c - half], [a0, b_c + half, z_c - half],
                  [a0, b_c + half, z_c + half], [a0, b_c - half, z_c + half],
                  [a_top, b_c - half, z_c - half],
                  [a_top, b_c + half, z_c - half],
                  [a_top, b_c + half, z_c + half],
                  [a_top, b_c - half, z_c + half]], dtype=np.float64)
    f = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
                  [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                  [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]],
                 dtype=np.int64)
    return v, f


def _synthetic(depths, stations, r_ground=25.0, half=0.004):
    """A place dict with boxes reaching the given depths at the given stations."""
    frames, faces = [], []
    for d, b in zip(depths, stations):
        v, f = _prism(r_ground + d, b, 0.0, half)
        frames.append(v)
        faces.append(f)
    return {"frames": frames, "faces": faces, "r_ground": r_ground,
            "theta_c": 0.0, "gov": 0, "baked": frames,
            "protrusion_um": [0.0] * len(frames)}


# --------------------------------------------------------------------------

def main() -> int:
    from semgrit.analysis import wheel_motion
    from semgrit.build_deck import build_deck, workpiece_of
    from semgrit.hybrid import chip_field
    from semgrit.rigid_wheel import place_workpiece
    from semgrit.wheel_workpiece import WorkpieceBlock
    from semgrit_multi.build import MultiParams, build_multi
    from semgrit_multi.envelope import (EnvelopeParams, nodal_field,
                                        sample_grit_surface, sweep_envelope)
    from semgrit_multi.fieldinject import (InjectError, inject_field,
                                           read_field)

    print("=" * 78)
    print("swept-envelope chip-thickness verification")
    print("=" * 78)

    # ------------------------------------------------------------------
    print("1. surface sampling")
    v, f = _prism(25.01, 0.0, 0.0, 0.004)
    for sub, least in ((1, 8), (2, 20), (3, 40)):
        pts = sample_grit_surface(v, f, sub)
        check("facet_subdivision=%d gives at least %d points" % (sub, least),
              len(pts) >= least, "%d points" % len(pts))
        # every sample must lie on the body it came from
        check("subdivision=%d samples stay inside the box" % sub,
              bool((pts[:, 0] >= v[:, 0].min() - 1e-12).all()
                   and (pts[:, 0] <= v[:, 0].max() + 1e-12).all()))
    check("vertices are all present at subdivision 1",
          len(np.unique(np.round(sample_grit_surface(v, f, 1), 12),
                        axis=0)) == 8)

    # ------------------------------------------------------------------
    print("2. shadowing, on geometry whose answer is arithmetic")
    wp = WorkpieceBlock(length_mm=0.048, width_mm=0.015, depth_mm=0.006,
                        element_size_mm=0.0003)
    # No infeed and no rotation would mean nothing sweeps, so give it a tiny
    # rotation and zero infeed: then each box cuts exactly its own depth.
    step = 1.0e-6
    omega = -1200.0
    mot = {"vr3": omega, "radial_speed_mm_s": 0.0, "omega_rad_s": 1200.0}
    ep = EnvelopeParams(depth_resolution_mm=1e-8, facet_subdivision=2)

    # (a) one box, depth 1 um
    pl = _synthetic([0.001], [0.0])
    env = sweep_envelope(pl, mot, wp, step_time_s=step, params=ep)
    close("one box of depth 1 um gives h = 1 um", env.h_elem[env.cut].max(),
          0.001, rtol=2e-3)
    close("and removes 1 um of depth", env.depth_removed.max(), 0.001,
          rtol=2e-3)

    # (b) deep box first, then a shallower one over the SAME station range.
    #     Same starting station means the same crossing window, so the sort is
    #     a tie broken by index and grit 0 leads deterministically. Offsetting
    #     them in station instead would not separate them at all: every box
    #     sweeps the whole block, so a station offset is only a time offset.
    pl = _synthetic([0.001, 0.0005], [0.0, 0.0])
    env2 = sweep_envelope(pl, mot, wp, step_time_s=step, params=ep)
    order = env2.grit_order
    check("the leading grit is swept first", order[0] == 0, str(order))
    h_first = env2.per_grit_h[0][0]
    h_second = env2.per_grit_h[1][0]
    close("the deeper leading grit takes its full depth", h_first, 0.001,
          rtol=2e-3)
    check("the shallower following grit takes nothing", h_second <= 1e-9,
          "h = %.4g mm" % h_second)
    close("the groove is still 1 um deep", env2.depth_removed.max(), 0.001,
          rtol=2e-3)

    # (c) shallow first, then deeper: the second takes only the difference.
    pl = _synthetic([0.0005, 0.001], [0.0, 0.0])
    env3 = sweep_envelope(pl, mot, wp, step_time_s=step, params=ep)
    close("the shallow leading grit takes 0.5 um",
          env3.per_grit_h[0][0], 0.0005, rtol=3e-3)
    close("the deeper following grit takes only the extra 0.5 um",
          env3.per_grit_h[1][0], 0.0005, rtol=3e-3)
    close("together they remove 1 um", env3.depth_removed.max(), 0.001,
          rtol=2e-3)
    check("removing the same total by two grits gives the same groove",
          abs(env3.depth_removed.max() - env2.depth_removed.max()) < 2e-6)

    # (d) two grits on different stations do not interact at all
    pl = _synthetic([0.001, 0.001], [0.010, -0.010])
    env4 = sweep_envelope(pl, mot, wp, step_time_s=step, params=ep)
    check("grits on separate stations both cut fully",
          min(env4.per_grit_h[0][0], env4.per_grit_h[1][0]) > 0.0009,
          "%.4g and %.4g mm" % (env4.per_grit_h[0][0], env4.per_grit_h[1][0]))

    # ------------------------------------------------------------------
    print("3. the real deck: agreement with brute force")
    lib = find_library(sys.argv)
    if not lib or not os.path.exists(lib):
        print("   no grain library found, so the deck checks are SKIPPED.")
        print("   Pass one with --library <grain_library.pkl>.")
        print("=" * 78)
        if _FAIL:
            print("%d passed, %d FAILED" % (_PASS, len(_FAIL)))
            for f in _FAIL:
                print("   - " + f)
            return 1
        print("%d CHECKS PASSED, deck checks SKIPPED (no grain library)"
              % _PASS)
        return 0
    print("   grain library: %s" % lib)
    with open(lib, "rb") as fh:
        solids = pickle.load(fh)["solids"]

    mp = MultiParams(name="env_single", grit_mode="single")
    dp = mp.deck_params()
    info, model = build_deck(dp, solids, None, return_model=True, dry_run=True)
    wpb = workpiece_of(dp)
    place = place_workpiece(model, wpb, dp.clearance_um, dp.wp_position,
                           dp.wp_position_deg, True)
    st = float(info["cost"]["step_time_s"])
    import dataclasses
    an = dataclasses.replace(dp.analysis, depth_of_cut_um=0.0)
    # resolve the automatic depth the same way the writer does
    from semgrit.build_deck import _auto_depth
    an = dataclasses.replace(
        an, depth_of_cut_um=_auto_depth(dp, place["clearance_um"]))
    motion = wheel_motion(an, place["theta_c"], dp.surface_speed_mm_s,
                          dp.outer_radius_mm, st)
    env = sweep_envelope(place, motion, wpb, step_time_s=st)

    # Brute force: the same kinematics, unchunked, 200k uniform samples over
    # the whole step, vertices and facet samples alike.
    cloud = sample_grit_surface(place["frames"][0], place["faces"][0], 2)
    hl, hw = wpb.length_mm / 2.0, wpb.width_mm / 2.0
    best = -np.inf
    for t in np.linspace(0.0, st, 200_001):
        al = motion["vr3"] * t
        a = (cloud[:, 0] * math.cos(al) - cloud[:, 1] * math.sin(al)
             + motion["radial_speed_mm_s"] * t)
        b = cloud[:, 0] * math.sin(al) + cloud[:, 1] * math.cos(al)
        ok = (b >= -hl) & (b < hl) & (np.abs(cloud[:, 2]) < hw)
        if ok.any():
            best = max(best, float((a[ok] - place["r_ground"]).max()))
    close("deepest penetration matches a 200k-sample brute force",
          env.depth_removed.max(), best, rtol=0.0,
          atol=3.0 * 2.0e-7)          # the stated 0.2 nm depth resolution

    # ------------------------------------------------------------------
    print("4. the single-grit closed form it replaces")
    fld = chip_field(place, motion, wpb, rotation_reversed=False)
    nl, nw, nd = wpb.divisions()
    u_edges = env.u_edges
    worst = 0.0
    worst_at = None
    n_cmp = 0
    for i in range(nl):
        col = env.h_elem[i, :, :][env.cut[i, :, :]]
        if col.size == 0:
            continue
        # h grows as u falls, so the deepest cut inside a station bin is at its
        # low-u edge -- that is what the bin's maximum records.
        want = fld.h_at(u_edges[i])
        got = float(col.max())
        n_cmp += 1
        if abs(got - want) > worst:
            worst, worst_at = abs(got - want), (i, got, want)
    check("the swept field reproduces the closed-form wedge at every cut "
          "station", worst < 5.0e-6,
          "%d stations compared, worst %.4f nm at station %s"
          % (n_cmp, worst * 1e6, None if worst_at is None else worst_at[0]))
    print("      worst station disagreement: %.4f nm (%d stations)"
          % (worst * 1e6, n_cmp))

    # ------------------------------------------------------------------
    print("5. internal consistency")
    check("every cut element has a positive chip thickness",
          bool((env.h_elem[env.cut] > 0).all()))
    check("no NaN survives the fill", bool(np.isfinite(env.h_elem).all()))
    check("depth_removed is the envelope of the per-grit depths",
          env.depth_removed.max() <= env.h_elem.max() + 1e-12)
    check("the cut mask is contiguous from the surface downwards",
          bool(np.all(np.diff(env.cut.astype(np.int8), axis=2) <= 0)),
          "an element cut below an uncut one would mean material was removed "
          "from under a lid")
    vol = float(env.depth_removed.sum() * (wpb.length_mm / nl)
                * (wpb.width_mm / nw))
    close("removed volume matches the depth map", env.stats[
        "removed_volume_mm3"], vol, rtol=1e-12)

    # ------------------------------------------------------------------
    print("6. element field -> nodal field")
    nod = nodal_field(env, wpb)
    check("one value per node", nod.size == (nl + 1) * (nw + 1) * (nd + 1),
          "%d values, %d nodes" % (nod.size,
                                   (nl + 1) * (nw + 1) * (nd + 1)))
    check("nodal values are bounded by the element values",
          nod.min() >= env.h_elem.min() - 1e-15
          and nod.max() <= env.h_elem.max() + 1e-15)
    # Averaging a LINEAR field over the eight nodes of a brick returns its
    # centroid value exactly, which is what the integration point will see.
    lin = np.add.outer(np.add.outer(np.arange(nl) * 3.0,
                                    np.arange(nw) * 5.0),
                       np.arange(nd) * 7.0) + 11.0
    from semgrit_multi.envelope import ChipEnvelope
    fake = ChipEnvelope(h_elem=lin, depth_removed=np.zeros((nl, nw)),
                        n_grits_engaged=0, grit_order=[], per_grit_h={},
                        u_edges=env.u_edges, z_edges=env.z_edges,
                        depth_edges=env.depth_edges,
                        cut=np.ones(lin.shape, dtype=bool))
    nl2 = nodal_field(fake, wpb).reshape(nl + 1, nw + 1, nd + 1)
    interp = np.zeros_like(lin)
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                interp += nl2[di:di + nl, dj:dj + nw, dk:dk + nd]
    interp /= 8.0
    err = np.abs(interp - lin)
    # Interior: exact. A boundary node has fewer elements around it, so its
    # average is biased inward -- unavoidable in any nodal representation, and
    # bounded by half the field's change over one element.
    close("the integration point recovers a linear field exactly in the "
          "interior", float(err[1:-1, 1:-1, 1:-1].max()), 0.0, atol=1e-12)
    # A boundary node averages fewer elements, so it reports its own element
    # rather than the extrapolated nodal value. The bias is bounded by half the
    # field's total variation across one element -- and it is measured, not
    # assumed: (3 + 5 + 7)/4 = 3.75 at the corner.
    bound = 0.5 * (3.0 + 5.0 + 7.0)
    check("and the boundary bias is under half the field's variation over one "
          "element", float(err.max()) <= bound + 1e-12,
          "worst %.4f, bound %.4f" % (err.max(), bound))

    # ------------------------------------------------------------------
    print("7. injection into the deck")
    res = build_multi(mp, solids, OUT, log=(print if VERBOSE else
                                            (lambda *a: None)))
    deck = res["path"]
    back = read_field(deck, 1)
    check("every field value is written", len(back) == nod.size,
          "%d read, %d written" % (len(back), nod.size))
    worst = max(abs(back[i + 1] - v) for i, v in enumerate(nod))
    check("the field round-trips out of the deck EXACTLY", worst == 0.0,
          "max |dh| = %.3g mm" % worst)

    lines = open(deck, encoding="ascii").readlines()
    i_ic = next(i for i, ln in enumerate(lines)
                if ln.lower().startswith("*initial conditions"))
    i_ea = next(i for i, ln in enumerate(lines)
                if ln.lower().startswith("*end assembly"))
    i_st = next(i for i, ln in enumerate(lines)
                if ln.lower().startswith("*step"))
    check("the field sits after *End Assembly and before *Step",
          i_ea < i_ic < i_st, "%d < %d < %d" % (i_ea, i_ic, i_st))
    check("the un-injected deck is kept alongside it",
          res["plain_path"] and os.path.exists(res["plain_path"]))
    plain = open(res["plain_path"], encoding="ascii").readlines()
    added = len(lines) - len(plain)
    check("injection adds only the field block", added == res["injected"][
        "n_lines_added"], "%d lines added" % added)
    k = next(i for i, (x, y) in enumerate(zip(plain, lines)) if x != y)
    check("and changes nothing else",
          plain == lines[:k] + lines[k + added:],
          "the deck is the original with %d lines inserted at %d" % (added, k))

    # refusals
    for tag, kwargs, want in (
            ("a field of the wrong length", dict(values=nod[:-1]), "nodes"),
            ("negative chip thickness", dict(values=-np.abs(nod)), "negative"),
            ("a non-finite value", dict(values=np.where(
                np.arange(nod.size) == 3, np.nan, nod)), "non-finite")):
        try:
            inject_field(res["plain_path"],
                         os.path.join(OUT, "_reject.inp"), **kwargs)
            check("refuses " + tag, False, "it did not")
        except InjectError as exc:
            check("refuses " + tag, want in str(exc), str(exc)[:60])
    # and a deck whose card does not read field variable 1
    from semgrit.build_deck import hybrid_single_grit
    from semgrit.hybrid import HybridParams, kic_from_mpa_sqrt_m
    hp0 = HybridParams(enabled=True, kic=kic_from_mpa_sqrt_m(0.3), h_source=0)
    info0 = build_deck(hybrid_single_grit(hybrid=hp0, name="env_coords"),
                       solids, OUT)
    try:
        inject_field(info0["path"], os.path.join(OUT, "_reject2.inp"), nod)
        check("refuses a deck that does not read field variable 1", False,
              "it did not")
    except InjectError as exc:
        check("refuses a deck that does not read field variable 1",
              "PROPS(56)" in str(exc), str(exc)[:70])

    # ------------------------------------------------------------------
    # The card-level checks. verify_hybrid_deck.py owns these for a deck whose
    # chip thickness comes from the four wedge constants; it cannot be used on
    # this one, because it would compare the VUMAT against a wedge the VUMAT is
    # not reading. So the parts that still apply live here.
    print("7b. the material card of an injected deck")
    from semgrit.analysis import JH2_SANDSTONE
    from semgrit.hybrid import (HYBRID_DELETE_SDV, HYBRID_DEPVAR,
                                N_HYBRID_PROPS, critical_depth_mm)
    from verify_hybrid_deck import check_user_material_format, parse_deck
    dk = parse_deck(deck)
    check_user_material_format("*User Material data lines carry 8 values each",
                               deck)
    check("*User Material declares %d constants" % N_HYBRID_PROPS,
          dk["n_declared"] == N_HYBRID_PROPS, str(dk["n_declared"]))
    check("%d constants are written" % N_HYBRID_PROPS,
          len(dk["props"]) == N_HYBRID_PROPS, str(len(dk["props"])))
    check("*Depvar is %d, deleting on SDV%d" % (HYBRID_DEPVAR,
                                                HYBRID_DELETE_SDV),
          dk["depvar"] == HYBRID_DEPVAR
          and dk["delete_sdv"] == HYBRID_DELETE_SDV,
          "%s / %s" % (dk["depvar"], dk["delete_sdv"]))
    check("props 1-17 are the unmodified JH-2 card",
          all(abs(dk["props"][i] - JH2_SANDSTONE[i]) < 1e-12
              for i in range(17)))
    close("dc in the card matches lambda_c, H, E, Kc as written",
          dk["props"][46],
          critical_depth_mm(dk["props"][47], dk["props"][48], dk["props"][27],
                            dk["props"][49], int(round(dk["props"][50]))),
          rtol=0.0, atol=0.0)
    close("dc matches the build report", dk["props"][46], res["dc_mm"],
          rtol=0.0, atol=0.0)
    hdr = chr(10).join(dk["header"])
    check("the header names the field and says where it came from",
          "CHIP THICKNESS FIELD" in hdr and "swept from" in hdr)
    check("the header warns the JC constants are placeholders",
          "PLACEHOLDER" in hdr.upper())

    # ------------------------------------------------------------------
    print("8. the compiled VUMAT agrees with the injected field")
    try:
        import verify_vumat_grind as vg
        drv = vg.Driver(vg.find_gfortran(),
                        os.path.join(HERE, "vumat_grind.for"), "grind.exe")
    except SystemExit as exc:
        print("   SKIPPED: %s" % str(exc).split("\n")[0])
        drv = None
    if drv is not None:
        from semgrit.hybrid import HYBRID_DEPVAR
        from verify_hybrid_deck import parse_deck
        props = parse_deck(deck)["props"]
        close("the card asks for field variable 1", props[55], 1.0, atol=0.0)
        dc = res["dc_mm"]
        seg = [(4, (0.0, 0.0, 0.0, 1.0e-9, 0.0, 0.0))]
        picks = np.unique(np.concatenate([
            np.linspace(0, nod.size - 1, 12).astype(int),
            np.argsort(np.abs(nod - dc))[:6]]))
        bad = []
        worst_h = 0.0
        for n in picks:
            hv = float(nod[n])
            r = drv.run(props, seg, nstatev=HYBRID_DEPVAR, fields=(hv,),
                        nout=4)
            worst_h = max(worst_h, abs(r[-1]["sdv"][13] - hv))
            got = int(round(r[-1]["sdv"][12]))
            want = 1 if hv < dc else 2
            if got != want:
                bad.append((int(n), hv, got, want))
        check("the VUMAT reads the field value exactly", worst_h == 0.0,
              "max |dh| = %.3g" % worst_h)
        check("and picks the branch the field implies at %d nodes"
              % len(picks), not bad, str(bad[:2]))

    # ------------------------------------------------------------------
    print("9. determinism")
    r2 = build_multi(mp, solids, os.path.join(OUT, "again"),
                     log=(lambda *a: None))
    h1 = np.load(os.path.join(OUT, mp.name + "_h_elem.npy"))
    h2 = np.load(os.path.join(OUT, "again", mp.name + "_h_elem.npy"))
    check("the same inputs give a bit-identical field",
          bool((h1 == h2).all()),
          "max |dh| = %.3g" % float(np.abs(h1 - h2).max()))

    # ------------------------------------------------------------------
    print("9b. tip_paths agrees with the sweep it illustrates")
    from semgrit_multi.envelope import tip_paths
    tp = tip_paths(place, motion, wpb, step_time_s=st)
    check("a tip path is returned for the grit that crosses", len(tp) == 1,
          str(sorted(tp)))
    p0 = tp[0]
    # The path is swept a little past the block on purpose, so the plot shows
    # the approach; the sweep only records what lands inside. Compare like with
    # like.
    hlb = wpb.length_mm / 2.0
    ins = (p0[:, 1] >= -hlb) & (p0[:, 1] < hlb)
    # tip_paths draws 400 samples for a picture; the sweep uses 14,000 for a
    # number. So the tolerance is the PATH's own sampling blur -- the depth it
    # gains between consecutive samples -- not an invented figure.
    blur = float(np.abs(np.diff(p0[:, 2])).max())
    close("inside the block the tip path reaches the depth the sweep recorded",
          float(p0[ins, 2].max()), float(env.depth_removed.max()), rtol=0.0,
          atol=2.0 * blur)
    print("      tip-path sampling blur %.3f nm" % (blur * 1e6))
    check("the path is swept past the block edge, for context in the plot",
          float(np.abs(p0[:, 1]).max()) > hlb)
    check("depth increases along the path, as the infeed says it must",
          bool(np.all(np.diff(p0[:, 2]) > -1e-12)))

    # ------------------------------------------------------------------
    print("9c. custom measured trajectories")
    from semgrit_multi.trajectory import (Trajectory, TrajectoryError,
                                          from_csv, from_points,
                                          sweep_trajectory)
    hl = wpb.length_mm / 2.0
    # A path whose answer is arithmetic: a straight ramp to a known depth.
    d_max = 4.0e-4
    uu = np.linspace(hl, -hl, 500)
    dd = np.linspace(0.0, d_max, 500)
    tr = from_points(np.column_stack([uu, dd]), columns="u,depth")
    check("a two-column table is read as u and depth", tr.n == 500)
    e_tr = sweep_trajectory(place, wpb, tr, step_time_s=st)
    close("the swept depth is the depth the path prescribed",
          float(e_tr.depth_removed.max()), d_max, rtol=0.0, atol=3.0 * 2.0e-7)
    check("and it is NOT the ideal sweep's depth, so the path really was used",
          abs(e_tr.depth_removed.max() - env.depth_removed.max()) > 1e-6,
          "measured %.4f um vs ideal %.4f um"
          % (e_tr.depth_removed.max() * 1000, env.depth_removed.max() * 1000))
    check("the trajectory field still classifies against dc",
          e_tr.split(res["dc_mm"])["n_cut"] > 0)

    # Units are the classic mistake: a profile in microns is 1000x too deep.
    tr_um = from_points(np.column_stack([uu * 1000.0, dd * 1000.0]),
                        columns="u,depth", scale_mm=1e-3)
    close("scale_mm converts a table given in microns", float(tr_um.depth.max()),
          d_max, rtol=1e-12)
    tr_neg = from_points(np.column_stack([uu, -dd]), columns="u,depth",
                         depth_sign=-1.0)
    close("depth_sign flips a profilometer's negative groove",
          float(tr_neg.depth.max()), d_max, rtol=1e-12)
    check("a negative-depth path is called out in the notes",
          any("negative depth" in n for n in
              from_points(np.column_stack([uu, -dd]),
                          columns="u,depth").notes))
    # round trip through a file
    csv_path = os.path.join(OUT, "traj.csv")
    with open(csv_path, "w") as fh:
        fh.write("u_mm,depth_mm" + chr(10))
        for a_, b_ in zip(uu, dd):
            fh.write("%.9g,%.9g%s" % (a_, b_, chr(10)))
    tr_csv = from_csv(csv_path, columns="u,depth")
    check("a CSV with a header round-trips", tr_csv.n == 500
          and abs(float(tr_csv.depth.max()) - d_max) < 1e-15)
    for tag, fn in (
            ("a one-column table", lambda: from_points(uu.reshape(-1, 1))),
            ("a single sample", lambda: from_points(np.zeros((1, 2)),
                                                    columns="u,depth")),
            ("a non-finite value", lambda: from_points(
                np.array([[0.0, 0.0], [1.0, np.nan]]), columns="u,depth"))):
        try:
            fn()
            check("refuses " + tag, False, "it did not")
        except TrajectoryError:
            check("refuses " + tag, True)
    try:
        from_points(np.column_stack([uu * 1000.0, dd]),
                    columns="u,depth").clipped_to_block(wpb)
        check("refuses a path that misses the block", False, "it did not")
    except TrajectoryError as exc:
        check("refuses a path that misses the block", "millimetres" in str(exc))
    # resampling and retiming must not move the geometry
    close("resampling preserves the deepest point",
          float(tr.resampled(97).depth.max()), d_max, rtol=1e-12)
    close("retiming does not touch the geometry",
          float(tr.retimed(0.0, 1e-6).depth.max()), d_max, rtol=0.0, atol=0.0)

    # ------------------------------------------------------------------
    print("10. the SWMODE rewrite for vumat_grind2.for")
    from semgrit_multi.swmode import (DEPVAR_GRIND2, N_PROPS_GRIND2,
                                      SwModeError, set_energy_mode)
    e_deck = os.path.join(OUT, "energy_mode.inp")
    inf = set_energy_mode(deck, e_deck, 1, 0.0)
    de = parse_deck(e_deck)
    check("the card now declares %d constants" % N_PROPS_GRIND2,
          de["n_declared"] == N_PROPS_GRIND2, str(de["n_declared"]))
    check("and %d are written" % N_PROPS_GRIND2,
          len(de["props"]) == N_PROPS_GRIND2, str(len(de["props"])))
    check("the first 56 are untouched",
          de["props"][:56] == dk["props"],
          "worst %.3g" % max((abs(a - b) for a, b in
                              zip(de["props"][:56], dk["props"])),
                             default=0.0))
    close("PROPS(57) is SWMODE", de["props"][56], 1.0, atol=0.0)
    close("PROPS(58) is PSI", de["props"][57], 0.0, atol=0.0)
    check("*Depvar grew to %d, keeping delete=%d" % (DEPVAR_GRIND2,
                                                     de["delete_sdv"] or 0),
          de["depvar"] == DEPVAR_GRIND2 and de["delete_sdv"] == 12,
          "%s / %s" % (de["depvar"], de["delete_sdv"]))
    # Nothing but the card and the note may change.
    a_lines = open(deck, encoding="ascii").readlines()
    b_lines = open(e_deck, encoding="ascii").readlines()
    a_body = [x for x in a_lines if not x.startswith("**")]
    b_body = [x for x in b_lines if not x.startswith("**")]
    # A zip comparison is useless here: inserting one line shifts every line
    # after it. Count the actual edit operations instead.
    import difflib
    ins = rep = dele = 0
    for op, i1, i2, j1, j2 in difflib.SequenceMatcher(
            None, a_body, b_body, autojunk=False).get_opcodes():
        if op == "insert":
            ins += j2 - j1
        elif op == "delete":
            dele += i2 - i1
        elif op == "replace":
            rep += max(i2 - i1, j2 - j1)
    # A line COUNT is the wrong assertion here, and it was hiding a real bug.
    # The rewriter used to append a data line holding SWMODE and PSI, which was
    # safe only because 56 = 7*8 exactly; at 57 constants that leaves a card of
    # 8,8,8,8,8,8,8,1,2 and Abaqus reads *User Material EIGHT to a line. The
    # rewriter now repacks the whole block, so the edit is 3 replaced rather
    # than 1 inserted + 2 replaced. What actually matters is that every changed
    # line belongs to the *User Material block or is the *Depvar data line --
    # i.e. that nothing structural moved -- so assert that instead.
    # A set, not `x not in a_body`. The bodies are ~500,000 lines each, so the
    # list form is 2.5e11 comparisons and the gate simply never returns -- which
    # is how it was first written here and it hung the whole build.
    a_set = set(a_body)
    changed = [x for x in b_body if x not in a_set]
    def _is_card_or_depvar(t):
        t = t.strip()
        if t.lower().startswith("*user material"):
            return True
        # a constants data line: only numbers and commas
        try:
            [float(v) for v in t.split(",") if v.strip()]
            return True
        except ValueError:
            return False
    check("outside the comments, only the card and *Depvar data lines change",
          dele == 0 and all(_is_card_or_depvar(x) for x in changed),
          "%d inserted, %d replaced, %d deleted; %d changed lines, all in the "
          "card: %s" % (ins, rep, dele, len(changed),
                        all(_is_card_or_depvar(x) for x in changed)))
    check("no keyword line other than *User Material is touched",
          not [x for x in changed
               if x.startswith("*")
               and not x.lower().startswith("*user material")],
          str([x.strip() for x in changed if x.startswith("*")][:3]))
    check_user_material_format(
        "the rewritten card still carries 8 values per line", e_deck)
    check("the note says which subroutine to submit with",
          any("vumat_grind2.for" in x for x in b_lines if x.startswith("**")))
    for tag, bad in (("a card that is already %d constants"
                  % N_PROPS_GRIND2, e_deck),):
        try:
            set_energy_mode(bad, os.path.join(OUT, "_r.inp"), 1, 0.0)
            check("refuses " + tag, False, "it did not")
        except SwModeError as exc:
            check("refuses " + tag, "already" in str(exc), str(exc)[:60])
    try:
        set_energy_mode(deck, os.path.join(OUT, "_r.inp"), 7, 0.0)
        check("refuses an invalid SWMODE", False, "it did not")
    except SwModeError:
        check("refuses an invalid SWMODE", True)

    # And vumat_grind2.for must actually read it back.
    if drv is not None:
        g2 = vg.Driver(vg.find_gfortran(),
                       os.path.join(HERE, "vumat_grind2.for"), "grind2.exe")
        hv = float(nod[int(np.argmax(nod))])
        r = g2.run(de["props"], [(20, (0, 0, 0, 1e-9, 0, 0))],
                   nstatev=DEPVAR_GRIND2, fields=(hv,), nout=20)
        check("vumat_grind2.for reads the rewritten card and starts ductile",
              int(round(r[-1]["sdv"][12])) == 1,
              "SWMODE 1 must ignore h = %.4g mm >= dc" % hv)
        check("and reports an energy ratio slot",
              len(r[-1]["sdv"]) == DEPVAR_GRIND2)

    print("=" * 78)
    if _FAIL:
        print("%d passed, %d FAILED" % (_PASS, len(_FAIL)))
        for f in _FAIL:
            print("   - " + f)
        return 1
    print("ALL %d CHECKS PASSED" % _PASS)
    return 0


if __name__ == "__main__":
    sys.exit(main())
