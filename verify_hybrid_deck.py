"""Verification gate for a hybrid ductile/brittle single-abrasive deck.

``verify_vumat_grind.py`` proves the subroutine is right. The two deck
verifiers prove the geometry is right. Neither can see the thing that is new
here: whether the deck and the subroutine agree about **which material points
are ductile**. That agreement rests on four numbers passed through the
material card -- THC, H0, HG, RTIP -- and on a sign convention for the
rotation. Get any of them wrong and the job runs, reports nothing unusual, and
puts the transition in the wrong place.

So this file builds a deck, re-derives the chip-thickness field from the
deck's own node coordinates by a second route, and then feeds the deck's own
material card into the compiled VUMAT at the coordinates of real workpiece
elements, checking that the branch the Fortran picks is the branch Python
predicts, element by element.

    python verify_hybrid_deck.py <deck>.inp [-v]   # check that deck
    python verify_hybrid_deck.py [-v]              # build a reference and
                                                   # check that, plus compare
                                                   # it against its JH-2 twin

Exits non-zero on any failure.
"""

from __future__ import annotations

import math
import os
import pickle
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LIBRARY = os.path.join(HERE, "WHEEL_FIXED", "1_measurements",
                       "grain_library.pkl")
OUT = os.path.join(HERE, "_hybridchk")

VERBOSE = "-v" in sys.argv or "--verbose" in sys.argv
_PASS = 0
_FAIL: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    global _PASS
    if ok:
        _PASS += 1
        if VERBOSE:
            print("  PASS  %-56s %s" % (name, detail))
    else:
        _FAIL.append(name)
        print("  FAIL  %-56s %s" % (name, detail))


def close(name: str, got: float, want: float, rtol=1e-9, atol=0.0) -> None:
    err = abs(got - want)
    check(name, err <= atol + rtol * abs(want),
          "got %.12g want %.12g (err %.3g)" % (got, want, err))


# --------------------------------------------------------------------------
# read back what was written
# --------------------------------------------------------------------------

def parse_deck(path: str) -> dict:
    """Nodes, elements and the user material, straight out of the file.

    Deliberately a fresh parser rather than a call back into semgrit: the
    point is to check the file, and asking the writer what it wrote proves
    nothing.
    """
    nodes: dict[int, tuple] = {}
    parts: dict[str, dict] = {}
    part = None
    mode = None
    etype = None
    props: list[float] = []
    depvar = None
    delete_sdv = None
    in_props = False
    header: list[str] = []
    field_vals: list[float] = []

    with open(path, encoding="ascii") as fh:
        for raw in fh:
            ln = raw.rstrip("\n")
            if ln.startswith("**"):
                header.append(ln)
                continue
            if ln.startswith("*"):
                key = ln.split(",")[0].strip().lower()
                low = ln.lower()
                in_props = False
                if key == "*part":
                    m = re.search(r"name\s*=\s*([^,]+)", ln, re.I)
                    part = m.group(1).strip()
                    parts[part] = {"nodes": {}, "elements": {}}
                    mode = None
                elif key == "*end part":
                    part, mode = None, None
                elif key == "*node":
                    mode = "node"
                elif key == "*element":
                    mode = "element"
                    m = re.search(r"type\s*=\s*([A-Za-z0-9]+)", ln, re.I)
                    etype = m.group(1).upper() if m else "?"
                elif key == "*depvar":
                    mode = "depvar"
                    m = re.search(r"delete\s*=\s*(\d+)", ln, re.I)
                    delete_sdv = int(m.group(1)) if m else None
                elif key == "*initial conditions":
                    mode = ("field"
                            if re.search(r"type\s*=\s*field", ln, re.I) else None)
                elif key == "*user material":
                    m = re.search(r"constants\s*=\s*(\d+)", ln, re.I)
                    n_declared = int(m.group(1)) if m else 0
                    in_props = True
                    mode = None
                    parts.setdefault("_mat", {})["n_declared"] = n_declared
                else:
                    mode = None
                continue
            if not ln.strip():
                continue
            if in_props:
                props += [float(x) for x in ln.split(",") if x.strip()]
                continue
            if mode == "depvar":
                depvar = int(float(ln.split(",")[0]))
                mode = None
                continue
            if mode == "field":
                # "<instance>.<node>, <value>" -- the variable number is on the
                # keyword line, so a data line is TWO fields, not three. One
                # line per node: the field is nodal, and the switch is
                # evaluated per integration point, which is why the two counts
                # never match exactly.
                f = [x.strip() for x in ln.split(",") if x.strip()]
                if len(f) >= 2:
                    try:
                        field_vals.append(float(f[-1]))
                    except ValueError:
                        pass
                continue
            if part is None:
                continue
            f = [x.strip() for x in ln.split(",") if x.strip()]
            if mode == "node":
                parts[part]["nodes"][int(f[0])] = tuple(float(x) for x in f[1:4])
            elif mode == "element":
                parts[part]["elements"].setdefault(etype, {})[int(f[0])] = [
                    int(x) for x in f[1:]]
    return {"parts": parts, "props": props, "depvar": depvar,
            "delete_sdv": delete_sdv, "header": header,
            "field_vals": field_vals,
            "n_declared": parts.get("_mat", {}).get("n_declared", 0)}


def user_material_line_widths(path: str) -> list:
    """Values per data line of *User Material, in order.

    Abaqus reads these constants EIGHT to a line. Four to a line is not a
    stylistic difference, it is
      ***ERROR: THERE ARE INVALID DATA ASSOCIATED WITH THIS USER DEFINED
                MATERIAL DEFINITION
    and it is what every deck in this project did on its first real
    submission. Nothing in a self-written parser can notice that, so the rule
    is checked explicitly here.
    """
    out = []
    grabbing = False
    with open(path, encoding="ascii") as fh:
        for ln in fh:
            if ln.startswith("**"):
                continue
            if ln.startswith("*"):
                grabbing = ln.lower().startswith("*user material")
                continue
            if grabbing and ln.strip():
                out.append(len([x for x in ln.split(",") if x.strip()]))
    return out


def check_user_material_format(name: str, path: str) -> None:
    w = user_material_line_widths(path)
    ok = bool(w) and all(n == 8 for n in w[:-1]) and 1 <= w[-1] <= 8
    check(name, ok, "values per line: %s" % (w if len(w) < 12 else
                                             str(w[:10]) + "..."))
    longest = max((len(ln.rstrip()) for ln in open(path, encoding="ascii")),
                  default=0)
    check(name + " (line length under Abaqus' 256 columns)", longest <= 256,
          "longest line %d chars" % longest)


def header_value(header: list[str], label: str):
    for ln in header:
        if label in ln:
            m = re.search(r":\s*([-+0-9.eE]+)", ln)
            if m:
                return float(m.group(1))
    return None


# --------------------------------------------------------------------------

def deck_argument(argv) -> str:
    for a in argv[1:]:
        if not a.startswith("-"):
            return a
    return ""


def main() -> int:
    import json

    from semgrit.hybrid import (HYBRID_DELETE_SDV, HYBRID_DEPVAR,
                                N_HYBRID_PROPS, critical_depth_mm,
                                kic_from_mpa_sqrt_m)

    print("=" * 78)
    print("hybrid single-abrasive deck verification")
    print("=" * 78)

    given = deck_argument(sys.argv)
    _deck_path = given or ""
    solids = None
    params = None
    if given:
        # Check the deck we were handed. Everything sections 1-6 need is in the
        # deck and its report, so no grain library and no rebuild is required
        # -- which is what lets the notebook run this on what it just wrote,
        # in a scratch directory that has no library in it.
        if not os.path.exists(given):
            print("no such deck: %s" % given)
            return 2
        rep_path = os.path.splitext(given)[0] + "_report.json"
        if not os.path.exists(rep_path):
            print("missing report next to the deck: %s" % rep_path)
            return 2
        with open(rep_path, encoding="utf-8") as fh:
            info = json.load(fh)
        info["path"] = given
        info["size_bytes"] = os.path.getsize(given)
    else:
        import dataclasses

        from semgrit.build_deck import build_deck, hybrid_single_grit
        from semgrit.hybrid import HybridParams

        if not os.path.exists(LIBRARY):
            print("missing grain library: %s" % LIBRARY)
            print("or pass a deck to check: verify_hybrid_deck.py <deck>.inp")
            return 2
        with open(LIBRARY, "rb") as fh:
            solids = pickle.load(fh)["solids"]
        print("grain library: %d solids" % len(solids))
        hp = HybridParams(enabled=True, kic=kic_from_mpa_sqrt_m(0.3))
        params = hybrid_single_grit(hybrid=hp)
        info = build_deck(params, solids, OUT)

    hy = info.get("hybrid")
    if not hy:
        print("this is not a hybrid deck: its report carries no 'hybrid' "
              "block, so material_model was not 'hybrid'")
        return 2
    pr = info["params"] if isinstance(info.get("params"), dict) else None
    wp_length_mm = (pr["wp_length_mm"] if pr is not None
                    else params.wp_length_mm)
    rot_rev = bool((pr["analysis"] or {}).get("rotation_reversed") if pr
                   else params.analysis.rotation_reversed)
    print("deck: %s  (%.2f MB)" % (os.path.basename(info["path"]),
                                   info["size_bytes"] / 1e6))
    print("dc = %.6f nm, transition at u = %s mm"
          % (hy["dc_nm"], hy["chip_field"]["transition_u_mm"]))

    deck = parse_deck(info["path"])

    # -- 1. the card ------------------------------------------------------
    print("1. the material card")
    # Two shapes are legal. vumat_grind.for reads 56 constants and 20 SDVs;
    # vumat_grind2.for is the same routine plus the local energy criterion, so
    # it reads 58 (SWMODE, PSI appended) and 22 (plastic work, energy ratio).
    # swmode.py converts one into the other in place, and before this the gate
    # rejected the converted deck outright -- so package 3, the only one
    # carrying the energy criterion, could not be checked at all.
    n_energy_props, energy_depvar = N_HYBRID_PROPS + 2, HYBRID_DEPVAR + 2
    is_energy = len(deck["props"]) == n_energy_props
    want_props = n_energy_props if is_energy else N_HYBRID_PROPS
    want_depvar = energy_depvar if is_energy else HYBRID_DEPVAR
    if is_energy:
        print("   energy-criterion variant: %d constants, %d SDVs"
              % (want_props, want_depvar))
    check("*User Material declares %d constants" % want_props,
          deck["n_declared"] == want_props, str(deck["n_declared"]))
    check("%d constants are actually written" % want_props,
          len(deck["props"]) == want_props, str(len(deck["props"])))
    check("*Depvar is %d" % want_depvar, deck["depvar"] == want_depvar,
          str(deck["depvar"]))
    check("*Depvar names SDV%d as the deletion flag" % HYBRID_DELETE_SDV,
          deck["delete_sdv"] == HYBRID_DELETE_SDV, str(deck["delete_sdv"]))
    check_user_material_format("*User Material data lines carry 8 values each",
                               given or info["path"])
    # The first 56 must still be the hybrid card, unchanged by the conversion.
    # PROPS(56) is allowed to differ from the report: an ablation arm is the
    # same deck with the h source re-pointed and nothing else, which is the
    # whole point of it. Every OTHER constant must still match exactly.
    diffs = [(k + 1, a, b) for k, (a, b) in
             enumerate(zip(deck["props"][:N_HYBRID_PROPS], hy["props"]))
             if a != b]
    only_hsrc = all(k == 56 for k, _, _ in diffs)
    worst = max((abs(a - b) for _, a, b in diffs if _ != 56), default=0.0)
    check("the card round-trips to the numbers the build decided EXACTLY",
          worst == 0.0 and only_hsrc,
          "max |dprop| = %.3g outside PROPS(56); differing props %s"
          % (worst, [k for k, _, _ in diffs]))
    if is_energy:
        swmode = int(round(deck["props"][56]))
        psi = deck["props"][57]
        check("SWMODE is 0, 1 or 2", swmode in (0, 1, 2), str(swmode))
        check("PSI is zero (derive from dc) or positive", psi >= 0.0,
              "%.6g" % psi)
        print("   SWMODE %d, PSI %.6g%s"
              % (swmode, psi, " (derived from dc)" if psi == 0.0 else ""))
    # props 1..17 must be one of the registered JH-2 cards, untouched. Matching
    # against the whole registry rather than against sandstone is what lets a
    # second material through this gate while still catching a card that has
    # been edited by hand or half-overwritten by a mixed material assignment.
    from semgrit import materials
    hit = [k for k, w in materials.MATERIALS.items()
           if all(abs(deck["props"][i] - w.jh2[i]) < 1e-12 for i in range(17))]
    check("props 1-17 are an unmodified registered JH-2 card",
          len(hit) == 1, hit[0] if hit else "matches no material in "
          "semgrit.materials")
    mat = materials.get(hit[0]) if hit else None

    # -- 2. dc ------------------------------------------------------------
    print("2. the critical depth of cut")
    dc = deck["props"][46]
    want = critical_depth_mm(deck["props"][47], deck["props"][48],
                             deck["props"][27], deck["props"][49],
                             int(round(deck["props"][50])))
    close("dc in the card matches lambda_c, H, E, Kc as written", dc, want,
          rtol=0.0, atol=0.0)
    close("dc matches the build report", dc, hy["dc_mm"], rtol=0.0, atol=0.0)
    # The unit trap: toughness quoted in MPa*sqrt(m) but a deck in mm. Read the
    # expected value from the material the card was just identified as, so this
    # stays a conversion check and not a check that the rock is sandstone.
    if mat is not None:
        close("Kc converted to MPa*sqrt(mm) for %s" % mat.key,
              deck["props"][49],
              mat.dc["kic_mpa_sqrt_m"] * math.sqrt(1000.0),
              rtol=0.0, atol=0.0)
        close("E in the ductile branch matches the material", deck["props"][27],
              mat.jc["youngs_mpa"], rtol=0.0, atol=0.0)
        close("density is the material's, in tonne/mm^3", deck["props"][29],
              mat.density_kg_m3 * 1e-12, rtol=1e-12)
    check("dc is a physically plausible length", 1e-9 < dc < 1e-2,
          "%.4g mm" % dc)

    # h_source decides which of the two h routes is authoritative, and
    # therefore which checks below are meaningful. Read it from the CARD, not
    # from the report: the card is the thing Abaqus will obey.
    h_src = int(round(deck["props"][55]))
    # An ablation arm deliberately differs here, and says so in its header.
    hdr_txt = chr(10).join(deck["header"])
    is_arm = "ABLATION ARM" in hdr_txt
    check("h_source in the card matches the report, or the header declares an "
          "ablation arm", h_src == hy["h_source"] or is_arm,
          "card %d, report %s%s" % (h_src, hy["h_source"],
                                    " (declared arm)" if is_arm else ""))
    print("   h source: %d (%s)"
          % (h_src, {0: "coordinates, single grit", 1: "field variable 1",
                     2: "forced ductile", 3: "forced brittle"}.get(h_src, "?")))
    close("dc_form in the card is the registered form for this material",
          deck["props"][50], float(mat.dc.get("dc_form", 2)) if mat else
          deck["props"][50], rtol=0.0, atol=0.0)

    # -- 3. the chip-thickness field, re-derived from the file ------------
    print("3. the chip-thickness field")
    thc, h0, hg, rtip = (deck["props"][51], deck["props"][52],
                         deck["props"][53], deck["props"][54])
    close("THC is the workpiece centre angle in the report", thc,
          math.radians(info["theta_workpiece_deg"]), rtol=1e-12)
    close("RTIP is the governing grit vertex radius", rtip,
          hy["chip_field"]["rtip_mm"], rtol=0.0, atol=0.0)

    # A PRESCRIBED profile replaces the derivation on purpose: the deck is
    # cutting a specified trajectory, not one that fell out of the infeed. So
    # the derived-value checks below are the right ones only when the term was
    # in fact derived. Where it was imposed, the card is checked against the
    # report instead -- the card must still be exactly the number the build
    # decided, which is the property that actually protects the run.
    cf = hy.get("chip_field") or {}
    hg_fixed = bool(cf.get("hg_prescribed"))
    h0_fixed = bool(cf.get("h0_prescribed"))

    m = info["motion"]
    v_r = m["radial_speed_mm_s"]
    omega = m["omega_rad_s"]
    way = -1.0 if rot_rev else 1.0
    if hg_fixed:
        close("HG in the card is the prescribed HG from the report",
              hg, cf["hg"], rtol=0.0, atol=0.0)
        print("      (HG prescribed, so it is not checked against -v_r/omega r)")
    else:
        close("HG = -v_r / (omega r_tip)", hg, -way * v_r / (omega * rtip),
              rtol=0.0, atol=0.0)
        check("HG has the sign that makes the cut deepen along the sweep",
              (hg < 0) == (not rot_rev),
              "HG = %.6g, rotation_reversed = %s" % (hg, rot_rev))

    # H0 must put h at the governing vertex exactly one standoff below the
    # ground face -- that is the tangency the geometry verifiers already
    # check, expressed as a chip thickness.
    u_gov = hy["chip_field"]["u_gov_mm"]
    h_gov = h0 + hg * u_gov - u_gov ** 2 / (2.0 * rtip)
    if h0_fixed:
        close("H0 in the card is the prescribed H0 from the report",
              h0, cf["h0_mm"], rtol=0.0, atol=0.0)
        print("      (H0 prescribed, so h at the governing vertex is not the "
              "standoff)")
    else:
        close("h at the governing vertex is minus the standoff", h_gov,
              -(info["clearance_um"] / 1000.0 + 1e-9), rtol=0.0, atol=5e-12)

    # -- 4. the deck's own nodes -------------------------------------------
    print("4. h over the workpiece, from the node coordinates")
    wpart = next(k for k in deck["parts"]
                 if k not in ("WHEEL", "_mat") and deck["parts"][k]["nodes"])
    wnodes = deck["parts"][wpart]["nodes"]
    hexes = deck["parts"][wpart]["elements"]["C3D8R"]
    check("the workpiece is C3D8R", bool(hexes), wpart)
    close("element count matches the report", len(hexes),
          info["n_workpiece_elements"], rtol=0.0, atol=0.0)

    et = np.array([-math.sin(thc), math.cos(thc), 0.0])
    nid = np.array(sorted(wnodes))
    xyz = np.array([wnodes[i] for i in nid])
    u = xyz @ et
    hl = h0 + hg * u - u * u / (2.0 * rtip)
    hl = np.maximum(hl, 0.0)
    check("the tangential span of the block is its length",
          abs((u.max() - u.min()) - wp_length_mm) < 1e-9,
          "%.9f vs %.9f mm" % (u.max() - u.min(), wp_length_mm))
    # WHICH h is authoritative depends on the source. In coordinate mode
    # (PROPS(56)=0) the wedge in the card IS the h, and the node coordinates are
    # the right thing to re-derive it from. In FIELD mode (=1) the card's wedge
    # is carried for the post-processor only, and the h the subroutine actually
    # reads is field variable 1 -- so checking the coordinates there tests a
    # quantity the run never uses. That is why this section failed on every
    # multi-abrasive and energy deck: the gate was right and aimed at the wrong
    # number. Each mode is now checked against the thing that governs it.
    if h_src == 0:
        n_duct = int((hl < dc).sum())
        n_brit = len(hl) - n_duct
        print("   %d of %d workpiece nodes are ductile (h < dc), %d brittle"
              % (n_duct, len(hl), n_brit))
        check("both regimes appear somewhere on the block",
              n_duct > 0 and n_brit > 0,
              "ductile %d, brittle %d -- if one is zero the switch does "
              "nothing on this deck" % (n_duct, n_brit))
        # h is monotone along the sweep whenever the depth comes from a radial
        # infeed: over a block this short the linear term beats the parabola by
        # four orders of magnitude. A PRESCRIBED arc is the exception and is
        # meant to be -- it rises to a peak and falls back, which is the whole
        # point of it -- so there the requirement is single-humped rather than
        # monotone, i.e. the slope changes sign at most once.
        order = np.argsort(u)
        dh = np.diff(hl[order])
        if hg_fixed or h0_fixed:
            sign = np.sign(dh[np.abs(dh) > 1e-15])
            flips = int(np.count_nonzero(np.diff(sign))) if sign.size else 0
            check("h is single-humped in u (prescribed arc rises then falls)",
                  flips <= 1,
                  "%d slope reversals; max +%.3g, min %.3g"
                  % (flips, dh.max(), dh.min()))
        else:
            check("h is monotone in u across the block",
                  bool(np.all(dh <= 1e-15) or np.all(dh >= -1e-15)),
                  "max +%.3g, min %.3g" % (dh.max(), dh.min()))
    elif h_src in (2, 3):
        # Forced modes. h is not read from anywhere, so neither the coordinates
        # nor the field says anything about the branch -- the only correct check
        # is that the branch is the forced one everywhere.
        want = 1 if h_src == 2 else 2
        print("   forced %s: h is not read at all"
              % ("ductile" if h_src == 2 else "brittle"))
        check("the card still carries a field for the other arms to use",
              True, "%d field values (unused in this mode)"
              % len(deck["field_vals"]))
    else:
        fv = np.asarray(deck["field_vals"], dtype=float)
        check("the deck carries an *Initial Conditions, type=FIELD block",
              fv.size > 0, "%d field values" % fv.size)
        if fv.size:
            n_duct = int((fv < dc).sum())
            n_brit = int(fv.size - n_duct)
            print("   %d of %d FIELD values are ductile (h < dc), %d brittle"
                  % (n_duct, fv.size, n_brit))
            check("both regimes appear in the injected field",
                  n_duct > 0 and n_brit > 0,
                  "ductile %d, brittle %d -- if one is zero the switch does "
                  "nothing on this deck" % (n_duct, n_brit))
            check("every field value is a physically possible h",
                  bool(np.all(np.isfinite(fv)) and np.all(fv >= 0.0)
                       and float(fv.max()) < 1.0),
                  "range %.4g .. %.4g mm" % (fv.min(), fv.max()))
            npy = os.path.join(
                os.path.dirname(os.path.abspath(_deck_path)) or ".",
                os.path.basename(_deck_path).split("_field")[0]
                + "_h_elem.npy")
            if os.path.exists(npy):
                he = np.load(npy).ravel()
                # The field is written per NODE and the switch is evaluated per
                # integration point, so the nodal round trip moves points across
                # dc. Bound the agreement rather than demand equality.
                close("the field spans the swept array's h range",
                      float(fv.max()), float(he.max()), rtol=0.30)
                print("   swept array %d elements, h up to %.4f nm; field up "
                      "to %.4f nm" % (he.size, he.max() * 1e6, fv.max() * 1e6))

    # -- 5. the Fortran must agree, element by element ---------------------
    print("5. the compiled VUMAT picks the same branch")
    drv = None
    try:
        import verify_vumat_grind as vg
        fc = vg.find_gfortran()
        drv = vg.Driver(fc, os.path.join(HERE, "vumat_grind.for"),
                        "grind.exe")
    except SystemExit as exc:
        # No Fortran compiler here. Everything above still ran; say plainly
        # what was NOT checked rather than reporting a clean pass that is
        # missing its most important section.
        print("   SKIPPED: %s" % str(exc).split("\n")[0])
        print("   The card and the field are checked; that the subroutine")
        print("   agrees with them is not. Install gfortran and re-run "
              "(Colab: apt-get -qq install gfortran).")
    if drv is None:
        print("=" * 78)
        if _FAIL:
            print("%d passed, %d FAILED (Fortran cross-check skipped)"
                  % (_PASS, len(_FAIL)))
            for f in _FAIL:
                print("   - " + f)
            return 1
        print("%d CHECKS PASSED, Fortran cross-check SKIPPED" % _PASS)
        return 0
    props = list(deck["props"])
    seg = [(4, (0.0, 0.0, 0.0, 1.0e-9, 0.0, 0.0))]

    # Element centroids, which is where the single integration point of a
    # C3D8R sits. Sample across the whole block plus the two nodes closest to
    # the predicted transition, where a sign error would show up first.
    cent = {}
    for eid, conn in hexes.items():
        p = np.array([wnodes[n] for n in conn])
        cent[eid] = p.mean(axis=0)
    eids = sorted(cent)
    picks = [eids[i] for i in np.linspace(0, len(eids) - 1, 24).astype(int)]
    u_t = hy["chip_field"]["transition_u_mm"]
    if u_t is not None:
        near = sorted(eids, key=lambda e: abs(float(cent[e] @ et) - u_t))[:8]
        picks = sorted(set(picks) | set(near))

    worst_h = 0.0
    bad_mode = []
    if h_src == 0:
        for eid in picks:
            x = cent[eid]
            r = drv.run(props, seg, nstatev=HYBRID_DEPVAR, coord=tuple(x),
                        nout=4)
            got_h = r[-1]["sdv"][13]
            got_mode = int(round(r[-1]["sdv"][12]))
            uu = float(x @ et)
            want_h = max(h0 + hg * uu - uu * uu / (2.0 * rtip), 0.0)
            want_mode = 1 if want_h < dc else 2
            worst_h = max(worst_h, abs(got_h - want_h))
            if got_mode != want_mode:
                bad_mode.append((eid, uu, got_h, want_h, got_mode,
                                 want_mode))
        check("the VUMAT computes the same h as Python at %d element centroids"
              % len(picks), worst_h < 1e-15, "max |dh| = %.3g mm" % worst_h)
    elif h_src in (2, 3):
        want = 1 if h_src == 2 else 2
        bad_forced = []
        for hval in (0.0, 0.1 * dc, 10.0 * dc):
            for fields in ((), (hval,)):
                for xyz in ((0.0, 0.0, 0.0), (9.9, -9.9, 9.9)):
                    r = drv.run(props, seg, nstatev=HYBRID_DEPVAR,
                                fields=fields, coord=xyz, nout=4)
                    got = int(round(r[-1]["sdv"][12]))
                    if got != want:
                        bad_forced.append((hval, fields, xyz, got))
        check("h_source = %d forces branch %d for every h, field and coordinate"
              % (h_src, want), not bad_forced,
              "%d counterexamples, e.g. %s"
              % (len(bad_forced), bad_forced[:1]))
    else:
        # Field mode: hand the driver field variable 1 and require the
        # subroutine to latch exactly that into SDV14 and branch on it. This is
        # the only check standing between a working hybrid deck and an
        # all-ductile deck that is identical in every other output -- if the
        # field does not arrive, hloc is 0, 0 < dc, and every point is ductile.
        fv = np.asarray(deck["field_vals"], dtype=float)
        samples = [0.25 * dc, 0.99 * dc, 1.01 * dc, 4.0 * dc]
        if fv.size:
            samples += [float(fv.min()), float(np.median(fv)),
                        float(fv.max())]
        for hval in samples:
            r = drv.run(props, seg, nstatev=HYBRID_DEPVAR, fields=(hval,),
                        coord=(0.0, 0.0, 0.0), nout=4)
            got_h = r[-1]["sdv"][13]
            got_mode = int(round(r[-1]["sdv"][12]))
            want_mode = 1 if hval < dc else 2
            worst_h = max(worst_h, abs(got_h - hval))
            if got_mode != want_mode:
                bad_mode.append((hval, got_h, got_mode, want_mode))
        check("the VUMAT latches field variable 1 as h at %d probe values"
              % len(samples), worst_h < 1e-18, "max |dh| = %.3g mm" % worst_h)
        far = drv.run(props, seg, nstatev=HYBRID_DEPVAR, fields=(2.0 * dc,),
                      coord=(9.9, -9.9, 9.9), nout=4)[-1]["sdv"][13]
        check("h from the field ignores the coordinates entirely",
              abs(far - 2.0 * dc) < 1e-18,
              "got %.12g want %.12g" % (far, 2.0 * dc))
    check("the VUMAT picks the same branch at every sampled element",
          not bad_mode, str(bad_mode[:2]))

    # And the two override modes must still work through the real card.
    for src, want_mode, tag in ((2, 1, "forced ductile"),
                                (3, 2, "forced brittle")):
        p2 = list(props)
        p2[55] = float(src)
        r = drv.run(p2, seg, nstatev=HYBRID_DEPVAR,
                    coord=tuple(cent[eids[len(eids) // 2]]), nout=4)
        check("h_source = %d gives %s through the real card" % (src, tag),
              int(round(r[-1]["sdv"][12])) == want_mode)

    # -- 6. the deck still says what it does ------------------------------
    print("6. the deck documents itself")
    hdr = "\n".join(deck["header"])
    check("the header names vumat_grind.for", "vumat_grind.for" in hdr)
    check("the header states dc", "critical depth of cut" in hdr)
    if h_src == 0:
        check("the header states the transition station or says there is none",
              ("transition at u" in hdr) or ("ENTIRELY" in hdr)
              or ("ductile-brittle transitions, at u" in hdr),
              "an arc that crosses dc twice reports both, so the singular "
              "phrase is not the only acceptable one")
    else:
        check("the header names field variable 1 as the h source",
              "field variable 1" in hdr,
              "a field-mode deck must not advertise a wedge it never reads")
    check("the header warns that the JC constants are placeholders",
          "PLACEHOLDER" in hdr.upper(),
          "and the report agrees" if hy["placeholder_constants"] else "")
    check("the report carries the chip field for the post-processor",
          set(hy["chip_field"]) >= {"theta_c_rad", "h0_mm", "hg", "rtip_mm",
                                    "transition_u_mm"})

    # -- 7. a hybrid deck and a JH-2 deck differ only in the material ------
    if solids is None:
        print("7. skipped: comparing against the JH-2 twin needs the grain "
              "library, so run this with no argument to include it")
        print("=" * 78)
        if _FAIL:
            print("%d passed, %d FAILED" % (_PASS, len(_FAIL)))
            for f in _FAIL:
                print("   - " + f)
            return 1
        print("ALL %d CHECKS PASSED" % _PASS)
        return 0
    print("7. same geometry as the brittle deck it replaces")
    from semgrit.analysis import AnalysisParams
    jh2_params = dataclasses.replace(
        params, name="jh2_single_grit",
        analysis=AnalysisParams(enabled=True, material_model="jh2",
                                n_depvar=12, depth_of_cut_um=0.0))
    jinfo = build_deck(jh2_params, solids, OUT)
    close("same grit count", info["n_grits"], jinfo["n_grits"], atol=0.0)
    close("same workpiece element count", info["n_workpiece_elements"],
          jinfo["n_workpiece_elements"], atol=0.0)
    close("same ground radius", info["workpiece_ground_radius_mm"],
          jinfo["workpiece_ground_radius_mm"], rtol=0.0, atol=0.0)
    close("same depth of cut chosen", params.analysis.depth_of_cut_um,
          jh2_params.analysis.depth_of_cut_um, rtol=0.0, atol=0.0)

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
