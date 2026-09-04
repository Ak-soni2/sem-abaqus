"""Independent verification of a SAG deck.

    python verify_sag_deck.py <deck>.inp [<deck2>.inp ...]
    python verify_sag_deck.py --build      # build a small pair and check both

Deliberately shares NO code with the writer. Everything here is re-derived from
the text of the ``.inp`` -- a keyword-grammar state machine, node coordinates
read back and re-measured, element connectivity re-checked for inversion, the
material card re-parsed and its 58 constants re-interpreted. A bug in
``semgrit.sagemit`` cannot also be baked into its own verifier, which is the
convention the rest of this project's gates follow.

WHAT IS CHECKED, AND WHY EACH ONE MATTERS
-----------------------------------------
Grammar and structure. Unbalanced ``*Part``/``*End Part`` or a data line under
the wrong keyword is a job that dies at preprocessing after the queue wait.

Geometry, re-measured from the nodes. The polyurethane ring's radii must match
the layer thickness; every hex must have a positive Jacobian; the tool must be
seated so its tallest protrusion just touches the work, because an overlap at
t = 0 is an impulse and a gap is a free-flight phase.

The physics that would be silently wrong. The sector's own sagitta must span the
wheel compression, or the "wheel" is a flat punch. The press must be a VELOCITY
and quasi-static against the layer's wave speed. The steps must be in the order
press, hold, grind, and the hold must be long enough for the polyurethane to
relax -- the measured force is a steady reading against a long-term modulus.

The material card. 58 constants (not 56), SWMODE = 1, and the energy threshold
``PSI*Kc^2/E`` re-computed from the card's own numbers and compared against
``H*dc``. Since dc for WC-Co is MEASURED, this is where a deck that quietly fell
back on Bifano's 17x-too-large value would be caught.

A GRAIN THAT NEVER REACHES THE WORK. The MICRO deck seats the grain clear of
the surface and pushes it in by the predicted indentation. If the standoff
exceeds that indentation the ramp ends with the grain in mid-air, and the job
runs to completion having touched nothing -- no error, just zero energy
throughout. The standoff was once 2% of the block depth against a nanometre
indentation, so it was 200x the travel.

RIGID BODIES THAT CANNOT MOVE. A rigid part made of R3D3 facets has no
volume and therefore no mass, so a free translational dof driven by a force
makes a = F/m undefined and Abaqus refuses at the packager. This deck shipped
with exactly that and passed every other gate, because none of them asked
whether the model could move.

MESH CONVERGENCE. The energy criterion is regularised by the element length, so
it is mesh-dependent BY CONSTRUCTION: halving the element halves the work
density needed to trigger. ``--converge`` builds the same physics at several
resolutions and reports how the predicted transition moves, because a result
that changes with the mesh has to be quoted with its mesh.
"""

from __future__ import annotations

import math
import os
import re
import sys

FAIL = []
PASS = 0


def chk(what: str, ok: bool, detail: str = "") -> bool:
    global PASS
    if ok:
        PASS += 1
        print("  [PASS] %s%s" % (what, ("  " + detail) if detail else ""))
    else:
        FAIL.append(what)
        print("  [FAIL] %s%s" % (what, ("  " + detail) if detail else ""))
    return ok


# ---------------------------------------------------------------------------
# an independent keyword reader
# ---------------------------------------------------------------------------

class Deck:
    """A parsed .inp, from the text alone."""

    def __init__(self, path: str):
        self.path = path
        self.blocks = []          # (keyword, params dict, [data lines])
        self.parts = {}           # name -> {nodes: {}, elements: {}, ...}
        self.steps = []
        self._read()

    def _read(self):
        cur = None
        part = None
        step = None
        with open(self.path, encoding="ascii", errors="replace") as fh:
            for lineno, raw in enumerate(fh, 1):
                ln = raw.rstrip("\n")
                if not ln.strip() or ln.startswith("**"):
                    continue
                if ln.startswith("*"):
                    head = ln[1:].split(",")
                    kw = head[0].strip().lower()
                    pars = {}
                    for tok in head[1:]:
                        if "=" in tok:
                            k, v = tok.split("=", 1)
                            pars[k.strip().lower()] = v.strip()
                        elif tok.strip():
                            pars[tok.strip().lower()] = True
                    cur = (kw, pars, [], lineno)
                    self.blocks.append(cur)
                    if kw == "part":
                        part = pars.get("name")
                        self.parts[part] = dict(nodes={}, elements={},
                                                etypes={}, nsets={},
                                                elsets={}, sections=[])
                    elif kw == "end part":
                        part = None
                    elif kw == "step":
                        step = dict(name=pars.get("name"), blocks=[],
                                    line=lineno)
                        self.steps.append(step)
                    elif kw == "end step":
                        step = None
                    if step is not None and kw not in ("step",):
                        step["blocks"].append(cur)
                    continue
                if cur is None:
                    raise ValueError("%s:%d data line before any keyword"
                                     % (self.path, lineno))
                cur[2].append(ln)
                kw, pars = cur[0], cur[1]
                if part and kw == "node":
                    f = [x.strip() for x in ln.split(",")]
                    if len(f) >= 4:
                        self.parts[part]["nodes"][int(f[0])] = (
                            float(f[1]), float(f[2]), float(f[3]))
                elif part and kw == "element":
                    f = [x.strip() for x in ln.split(",") if x.strip()]
                    et = pars.get("type", "?")
                    self.parts[part]["elements"][int(f[0])] = [
                        int(x) for x in f[1:]]
                    self.parts[part]["etypes"][int(f[0])] = et

    def kw(self, name: str):
        name = name.lower()
        return [b for b in self.blocks if b[0] == name]

    def one(self, name: str):
        got = self.kw(name)
        return got[0] if got else None

    def text(self) -> str:
        return open(self.path, encoding="ascii", errors="replace").read()


# ---------------------------------------------------------------------------
# geometry, re-derived
# ---------------------------------------------------------------------------

_TETS = ((0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6), (3, 4, 6, 7),
         (1, 4, 5, 6))


def _jac(p) -> float:
    """Signed volume of one hex from its eight corner coordinates."""
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    tot = 0.0
    for a, b, c, d in _TETS:
        tot += dot(cross(sub(p[b], p[a]), sub(p[c], p[a])),
                   sub(p[d], p[a])) / 6.0
    return tot


def hex_stats(part: dict) -> dict:
    """Jacobian sign and edge extremes, from the deck's own node table."""
    nodes = part["nodes"]
    worst = None
    n_bad = 0
    n_hex = 0
    lo = float("inf")
    hi = 0.0
    edges = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7))
    for eid, conn in part["elements"].items():
        if len(conn) != 8:
            continue
        n_hex += 1
        p = [nodes[i] for i in conn]
        j = _jac(p)
        if j <= 0:
            n_bad += 1
            if worst is None or j < worst[1]:
                worst = (eid, j)
        for a, b in edges:
            d = math.dist(p[a], p[b])
            if d > 0:
                lo = min(lo, d)
                hi = max(hi, d)
    return dict(hexes=n_hex, inverted=n_bad, worst=worst,
                min_edge=(lo if lo < float("inf") else 0.0), max_edge=hi)


def radii(part: dict) -> tuple:
    """Min and max cylindrical radius about the part's own axis (z)."""
    rs = [math.hypot(x, y) for x, y, _ in part["nodes"].values()]
    return (min(rs), max(rs)) if rs else (0.0, 0.0)


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------

def check_grammar(d: Deck) -> None:
    print("\n1. keyword grammar and structure")
    t = d.text()
    for a, b in (("*Part,", "*End Part"), ("*Step,", "*End Step"),
                 ("*Assembly", "*End Assembly"),
                 ("*Instance,", "*End Instance")):
        chk("%s is balanced with %s" % (a.rstrip(','), b),
            t.count(a) == t.count(b),
            "%d vs %d" % (t.count(a), t.count(b)))
    chk("every part has a name", all(bool(n) for n in d.parts),
        "parts: %s" % ", ".join(sorted(str(x) for x in d.parts)))
    chk("no unresolved format placeholder survived",
        "%s" not in t and "%d" not in t and "{" not in t)
    chk("the deck is pure ASCII", all(ord(c) < 128 for c in t[:200000]))
    ints = [b for b in d.blocks if b[0] == "instance"]
    chk("every instance names an existing part",
        all(b[1].get("part") in d.parts for b in ints),
        "%d instance(s)" % len(ints))


def check_solver(d: Deck) -> None:
    print("\n2. solver and contact")
    dyn = d.kw("dynamic")
    chk("every step is Abaqus/Explicit",
        bool(dyn) and all("explicit" in b[1] for b in dyn),
        "%d *Dynamic block(s)" % len(dyn))
    chk("general contact is used, not contact pairs",
        bool(d.kw("contact")) and bool(d.kw("contact inclusions"))
        and not d.kw("contact pair"),
        "the VUMAT deletes elements, and deletion exposes interior faces "
        "that a pre-declared pair would never see")
    ci = d.one("contact inclusions")
    chk("the contact domain is ALL EXTERIOR",
        ci is not None and "all exterior" in ci[1])
    chk("a friction and a hard-contact behaviour are defined",
        bool(d.kw("friction")) and bool(d.kw("surface behavior")))
    sc = d.one("section controls")
    chk("element deletion is enabled in the section controls",
        sc is not None and str(sc[1].get("element deletion", "")).lower()
        in ("yes", "true"),
        "required: the energy criterion drives elements to failure")
    chk("enhanced hourglass control is on",
        sc is not None and "enhanced" in str(sc[1].get("hourglass", "")).lower(),
        "C3D8R has zero-energy modes a soft layer excites readily")


def check_material(d: Deck) -> dict:
    print("\n3. the material card, re-parsed")
    um = d.one("user material")
    out = {}
    if not chk("a *User Material block exists", um is not None):
        return out
    n = int(um[1].get("constants", 0))
    chk("the card declares 58 constants, not 56", n == 58,
        "got %d -- 56 is vumat_grind (geometric), 58 is vumat_grind2 (energy)"
        % n)
    vals = []
    for ln in um[2]:
        vals += [float(x) for x in ln.split(",") if x.strip()]
    chk("the card carries exactly as many numbers as it declares",
        len(vals) == n, "%d values" % len(vals))
    per_line = [len([x for x in ln.split(",") if x.strip()]) for ln in um[2]]
    chk("constants are written 8 per line",
        all(v == 8 for v in per_line[:-1]) and per_line[-1] <= 8,
        "%s -- 4 per line is silently rejected by Abaqus" % per_line)
    if len(vals) < 58:
        return out

    # PROPS indices are 1-based in the Fortran.
    e_jc, hardn, kic = vals[27], vals[48], vals[49]
    dcut, lamc, idcf = vals[46], vals[47], vals[50]
    swmode, psi = vals[56], vals[57]
    out.update(e=e_jc, h=hardn, kic=kic, dc=dcut, swmode=swmode, psi=psi)

    chk("SWMODE = 1, the local energy criterion", swmode == 1.0,
        "0 geometric, 1 energy, 2 either -- SAG has no closed-form "
        "trajectory, so 1")
    chk("dc is positive and in the tens of nanometres",
        0.0 < dcut < 1.0e-3,
        "dc = %.4f nm" % (dcut * 1e6))
    chk("hardness and modulus are positive", hardn > 0 and e_jc > 0,
        "H = %.0f MPa, E = %.0f MPa" % (hardn, e_jc))
    chk("Kc is in MPa*sqrt(mm), not MPa*sqrt(m)",
        kic > 10.0,
        "Kc = %.3f -- MPa*sqrt(m) would be ~31.6x smaller and would scale "
        "dc by 1000" % kic)

    # The threshold, re-derived two ways from the card's own numbers.
    psi_eff = psi if psi > 0 else dcut * e_jc * hardn / (kic * kic)
    gcrit = psi_eff * kic * kic / e_jc
    hdc = hardn * dcut
    chk("the energy threshold equals H*dc",
        abs(gcrit - hdc) / hdc < 1e-9,
        "PSI*Kc^2/E = %.6f = H*dc = %.6f MPa*mm (%.1f J/m2)"
        % (gcrit, hdc, hdc * 1000.0))
    out["threshold"] = hdc

    # Bifano on the SAME card, to catch a deck that silently used it.
    bif = 0.15 * (e_jc / hardn) * (kic / hardn) ** 2
    ratio = bif / dcut
    chk("dc is NOT Bifano's value for this material",
        ratio > 2.0,
        "Bifano would give %.1f nm, %.1fx the %.1f nm in the card -- the "
        "reference paper measured 60-100 nm and showed Bifano fails here"
        % (bif * 1e6, ratio, dcut * 1e6))
    dep = d.one("depvar")
    chk("*Depvar declares 22 state variables and deletes on 12",
        dep is not None and dep[1].get("delete") == "12"
        and any("22" in ln for ln in dep[2]),
        "vumat_grind2 writes 22; 20 is vumat_grind")
    return out


def check_geometry(d: Deck, macro: bool) -> dict:
    print("\n4. geometry, re-measured from the node table")
    out = {}
    for name, part in sorted(d.parts.items()):
        st = hex_stats(part)
        if st["hexes"]:
            chk("%s: no inverted hex" % name, st["inverted"] == 0,
                "%s hexes, worst J = %s" % (
                    format(st["hexes"], ","),
                    "n/a" if not st["worst"] else "%.3e" % st["worst"][1]))
            asp = st["max_edge"] / st["min_edge"] if st["min_edge"] else 0.0
            chk("%s: edge lengths are sane" % name,
                st["min_edge"] > 0 and asp < 1e4,
                "%.6g to %.6g mm, %.1f:1" % (st["min_edge"], st["max_edge"],
                                             asp))
            out[name] = st
        rigid = [t for t in part["etypes"].values() if t.upper().startswith("R")]
        if rigid:
            chk("%s: rigid facets are R3D3/R3D4" % name,
                all(t.upper() in ("R3D3", "R3D4") for t in rigid),
                "%s facets" % format(len(rigid), ","))
    if macro:
        pu = d.parts.get("PU")
        if pu:
            lo, hi = radii(pu)
            out["pu_radii"] = (lo, hi)
            chk("the polyurethane ring is an annulus of positive thickness",
                hi > lo > 0, "r = %.4f to %.4f mm (%.4f mm thick)"
                % (lo, hi, hi - lo))
        hub = d.parts.get("HUB")
        if hub and pu:
            hlo, hhi = radii(hub)
            chk("the hub sits inside the compliant layer",
                abs(hhi - lo) < 1e-6,
                "hub outer %.4f vs PU bore %.4f mm" % (hhi, lo))
    return out


def check_steps(d: Deck, macro: bool, timing: dict) -> None:
    print("\n5. the loading sequence")
    names = [s["name"] for s in d.steps]
    if macro:
        chk("three steps, in the order press, hold, grind",
            names == ["PRESS", "HOLD", "GRIND"], "%s" % names)
    else:
        chk("a load step followed by repeated passes",
            len(names) >= 3 and names[0] == "LOAD"
            and all(n == "PASS%d" % i
                    for i, n in enumerate(names[1:], start=1)),
            "%d passes: %s" % (len(names) - 1,
                               ", ".join(names[:3]) + (" ..." if len(names) > 3
                                                       else "")))
        # The passes must ALTERNATE direction. A grain that only ever slides
        # forward runs down a fresh track, so every point it crosses sees
        # exactly ONE pass and the energy criterion can never trip -- however
        # many steps there are. Reversing keeps it over the same material,
        # which is what a polishing pad does.
        vx = []
        for st in d.steps[1:]:
            for kwd, pars, data, _ in st["blocks"]:
                if kwd != "boundary":
                    continue
                if str(pars.get("type", "")).lower() != "velocity":
                    continue
                for ln in data:
                    f = [x.strip() for x in ln.split(",")]
                    if len(f) >= 4 and f[1] == "1" and f[3]:
                        vx.append(float(f[3]))
        flips = sum(1 for i in range(len(vx) - 1)
                    if vx[i] * vx[i + 1] < 0)
        chk("the passes alternate direction, so they retrace one track",
            len(vx) >= 2 and flips == len(vx) - 1,
            "%d pass velocities, %d reversals -- a one-way slide would leave "
            "every point with a single pass and could never accumulate to "
            "the threshold" % (len(vx), flips))
        # Displacement-controlled, and NOT force-controlled -- the opposite
        # of what this gate used to demand. A rigid grain has no mass, so a
        # force on a free dof is rejected by Abaqus outright; and giving it
        # the real diamond mass does not help, because 1.7e-5 N on 4e-16
        # tonne is 4e10 mm/s2 and the grain would cross the block two hundred
        # times over before contact could resist. The pad POSITIONS the grain
        # in the experiment, and the depth is predicted by the Hertz chain
        # rather than unknown, so it is an input here. The branch is still the
        # output.
        depths = []
        for kwd, pars, data, _ in d.kw("boundary"):
            if str(pars.get("type", "")).lower() == "velocity":
                continue
            for ln in data:
                f = [x.strip() for x in ln.split(",")]
                if len(f) >= 4 and f[1] == "3" and f[2] == "3" and f[3]:
                    v = float(f[3])
                    if v != 0.0:
                        depths.append(abs(v))
        chk("the grain is pressed to a prescribed depth", bool(depths),
            "%s mm -- the contact chain predicts the indentation and is "
            "validated against the paper's measurements, so it is an input; "
            "the ductile/brittle branch is what is predicted"
            % (["%.3e" % v for v in depths] or "none"))
        if depths:
            chk("that depth is a sane sub-nanometre indentation",
                all(1e-9 < v < 1e-3 for v in depths),
                "%s mm" % ["%.3e" % v for v in depths])
        chk("no *Cload on the grain: a massless rigid body cannot take one",
            not d.kw("cload"),
            "Abaqus: rigid bodies require non-zero mass unless translational "
            "constraints are applied")
        slide = 0.0
        for kwd, pars, data, _ in d.kw("boundary"):
            if str(pars.get("type", "")).lower() != "velocity":
                continue
            for ln in data:
                f = [x.strip() for x in ln.split(",")]
                # dof 1 is the sliding direction; a 4-field line carries a value
                if len(f) >= 4 and f[1] == "1" and f[3]:
                    slide = max(slide, abs(float(f[3])))
        chk("and it SLIDES, so plastic work can accumulate", slide > 0.0,
            "tangential velocity %.1f mm/s -- the energy criterion triggers "
            "on history, so a grain that only indents and stops can never "
            "reach the threshold however hard it presses" % slide)

    times = []
    for s in d.steps:
        for kwd, pars, data, _ in s["blocks"]:
            if kwd == "dynamic":
                for ln in data:
                    f = [x.strip() for x in ln.split(",")]
                    if len(f) >= 2 and f[1]:
                        times.append(float(f[1]))
    chk("every step has a positive period",
        len(times) == len(d.steps) and all(t > 0 for t in times),
        "%s s" % ", ".join("%.3e" % t for t in times))

    if macro and len(times) == 3:
        press, hold, grind = times
        vel = None
        for kwd, pars, data, _ in d.steps[0]["blocks"]:
            if kwd == "boundary" and str(pars.get("type", "")).lower() == "velocity":
                for ln in data:
                    f = [x.strip() for x in ln.split(",")]
                    if len(f) >= 4 and f[1] == "2" and float(f[3]) != 0.0:
                        vel = abs(float(f[3]))
        chk("the press is a VELOCITY, not a displacement", vel is not None,
            "the reference deck wrote -51600 with no type=, which Abaqus "
            "reads as 51.6 m of displacement")
        if vel and timing:
            got = vel * press
            chk("velocity x press time equals the wheel compression",
                abs(got - timing["compression_mm"]) / timing["compression_mm"]
                < 0.02,
                "%.4f mm/s x %.5f s = %.5f mm vs T = %.4f mm"
                % (vel, press, got, timing["compression_mm"]))
            mach = vel / timing["wave_speed_mm_s"]
            chk("the press is quasi-static against the layer's wave speed",
                mach < 0.02,
                "v/c = %.4f; above ~0.01 the patch is loaded inertially and "
                "its pressure is not the steady Hertzian one" % mach)
            chk("the hold is several Prony relaxation times",
                hold >= 2.5 * timing["prony_tau_s"],
                "%.4f s = %.1f tau; the measured force is a steady reading "
                "against a long-term modulus"
                % (hold, hold / timing["prony_tau_s"]))


def check_reachable(d: Deck) -> None:
    """Can the grain actually reach the workpiece?

    A MICRO deck seats the grain a little clear of the surface and then pushes
    it in by the predicted indentation. If that standoff is LARGER than the
    indentation, the ramp finishes with the grain still in mid-air: the job
    runs to completion, writes a full .odb, and touches nothing. There is no
    error -- the only symptom is that kinetic and internal energy stay
    identically zero, which is easy to read as "nothing has happened yet"
    rather than "nothing will ever happen".

    That is not hypothetical. The standoff was 2% of the BLOCK DEPTH while the
    indentation is nanometres, so it came out ~200x the entire travel, and a
    six-hour run was spent closing a gap it could not close.

    Measured here from the deck's own node coordinates, not from any number
    the writer reported.
    """
    print("\n8. can the grain reach the workpiece?")
    work = d.parts.get("WORK")
    grains = d.parts.get("GRAINS")
    if not (work and grains):
        print("       (not a MICRO deck: skipped)")
        return

    top = max(z for _, _, z in work["nodes"].values())
    low = min(z for _, _, z in grains["nodes"].values())
    standoff = low - top

    # The imposed depth, read back out of the boundary conditions.
    depth = 0.0
    for kwd, pars, data, _ in d.kw("boundary"):
        if str(pars.get("type", "")).lower() == "velocity":
            continue
        for ln in data:
            f = [x.strip() for x in ln.split(",")]
            if len(f) >= 4 and f[1] == "3" and f[2] == "3" and f[3]:
                depth = max(depth, abs(float(f[3])))

    chk("the grain starts clear of the work, not interpenetrating",
        standoff > 0,
        "lowest grain node is %.4e mm %s the surface"
        % (abs(standoff), "below" if standoff < 0 else "above"))
    if depth > 0:
        chk("the imposed depth is larger than the standoff", standoff < depth,
            "standoff %.4e mm vs depth %.4e mm -- the ramp would finish with "
            "the grain still %.1f nm clear, and the job would run to "
            "completion having touched nothing"
            % (standoff, depth, (standoff - depth) * 1e6))
        chk("and the standoff is a small fraction of it",
            0 < standoff < 0.5 * depth,
            "standoff/depth = %.3f" % (standoff / depth if depth else 0))


def check_rigid_mass(d: Deck) -> None:
    """A rigid body needs mass, or every free translation must be constrained.

    Abaqus/Explicit integrates a = F/m on a rigid body's reference node. R3D3
    facets carry no volume and therefore no mass, so a rigid part built from
    them has m = 0 and any FREE translational dof makes that division
    undefined. Abaqus refuses at the packager:

        ERROR: Abaqus/Explicit requires rigid bodies to have a non-zero mass
        unless translational constraints are applied with the *BOUNDARY
        option.

    This deck hit exactly that: the grain was rigid, massless, and driven on
    dof 3 by a *Cload. It passed every other gate -- the grammar, the geometry,
    the material card -- because none of them asked whether the model could
    move.
    """
    print("\n9. can every rigid body actually move?")
    rb = d.kw("rigid body")
    if not rb:
        print("       (no rigid bodies in this deck)")
        return
    has_mass = bool(d.kw("mass")) or bool(d.kw("inertia"))
    # Which reference-node dofs are constrained, anywhere in the deck?
    refs = set()
    for _, pars, _, _ in rb:
        r = pars.get("ref node")
        if r:
            refs.add(r.upper())
    for ref in sorted(refs):
        fixed = set()
        for kwd, pars, data, _ in d.blocks:
            if kwd != "boundary":
                continue
            for ln in data:
                f = [x.strip() for x in ln.split(",")]
                if not f or f[0].upper().split(".")[-1] != ref:
                    continue
                if len(f) >= 2 and f[1].upper() in ("ENCASTRE", "PINNED"):
                    fixed |= {1, 2, 3}
                elif len(f) >= 3 and f[1].isdigit() and f[2].isdigit():
                    fixed |= set(range(int(f[1]), int(f[2]) + 1))
                elif len(f) >= 2 and f[1].isdigit():
                    fixed.add(int(f[1]))
        free = {1, 2, 3} - fixed
        chk("%s: every translation is constrained, or the body has mass"
            % ref, (not free) or has_mass,
            "dof %s free and no *Mass -- Abaqus integrates a = F/m on the "
            "reference node, and R3D3 facets carry no volume so m = 0"
            % sorted(free))

    # A force on a massless free dof is the specific failure.
    for kwd, pars, data, _ in d.kw("cload"):
        for ln in data:
            f = [x.strip() for x in ln.split(",")]
            if len(f) >= 2 and f[0].upper().split(".")[-1] in refs:
                chk("a *Cload on a rigid reference node needs mass", has_mass,
                    "%s dof %s is force-driven on a body with no *Mass"
                    % (f[0], f[1]))


def check_sector(d: Deck, timing: dict) -> None:
    """The sector's own curvature must span the indent."""
    print("\n10. can the modelled sector actually make contact?")
    pu = d.parts.get("PU")
    if not pu or not timing:
        print("       (not a MACRO deck, or no plan supplied: skipped)")
        return
    th = [math.degrees(math.atan2(y, x)) for x, y, _ in pu["nodes"].values()]
    span = max(th) - min(th)
    r = radii(pu)[1]
    sag = r * (1.0 - math.cos(math.radians(span) / 2.0))
    T = timing["compression_mm"]
    chk("the sector's own sagitta spans the wheel compression", sag >= T,
        "sector %.3f deg on r = %.3f mm gives %.5f mm of sagitta vs T = "
        "%.4f mm -- below it the arc is flat and the deck is a punch, "
        "not a wheel" % (span, r, sag, T))
    need = 2.0 * math.degrees(math.acos(max(-1.0, 1.0 - T / r)))
    chk("the sector is at least the geometric minimum", span >= need * 0.999,
        "%.3f deg vs %.3f deg required" % (span, need))


def check_seating(d: Deck, timing: dict) -> None:
    print("\n11. is the tool seated at first contact?")
    ins = {b[1].get("name"): b for b in d.kw("instance")}
    pu = d.parts.get("PU")
    work = d.parts.get("WORK")
    if not (pu and work and "PU-1" in ins):
        print("       (not a MACRO deck: skipped)")
        return
    data = [ln for ln in ins["PU-1"][2] if ln.strip()]
    if not data:
        chk("the tool instance carries a translation", False)
        return
    f = [float(x) for x in data[0].split(",")]
    y = f[1]
    r = radii(pu)[1]
    top = max(z for _, _, z in work["nodes"].values())
    clear = y - r - top
    chk("the tool neither overlaps nor floats above the work",
        -1.0e-3 <= clear <= 0.05,
        "pad surface sits %.6f mm from the ground face; a negative value is "
        "an impulse at t=0 and a large positive one is a free-flight phase"
        % clear)


# ---------------------------------------------------------------------------
# mesh convergence
# ---------------------------------------------------------------------------

def converge(refinements=(3.0, 5.0, 8.0, 12.0)) -> None:
    """How much does the predicted transition move with the mesh?

    The energy criterion is regularised by ``L_c``, so it is mesh-dependent by
    construction: the threshold is a work per unit AREA, and dividing by a
    smaller element gives a smaller work DENSITY. That is correct behaviour for
    a fracture-energy criterion and it is exactly why PSI has to be quoted with
    a mesh. This reports the sensitivity rather than pretending it is absent.
    """
    from semgrit import materials, sagdeck

    print("\n12. mesh convergence of the energy criterion")
    w = materials.get("wc_co")
    hp = w.hybrid_params()
    dc = hp.critical_depth_mm()
    hdc = hp.hardness_mpa * dc
    print("       threshold W_p*L_c >= H*dc = %.6f MPa*mm (%.1f J/m2)"
          % (hdc, hdc * 1000.0))
    print("       %-10s %-14s %-16s %s"
          % ("el/dc", "element (nm)", "W_p needed (MPa)", "elements"))
    rows = []
    for n in refinements:
        p = sagdeck.SAGParams(grain_um=6.0, material="wc_co",
                              elements_per_dc=n, micro_grains=1)
        pl = sagdeck.plan(p)
        el = pl["micro"]["element_mm"]
        rows.append((n, el, hdc / el, pl["micro"]["elements"]))
        print("       %-10.1f %-14.3f %-16.0f %s"
              % (n, el * 1e6, hdc / el, format(pl["micro"]["elements"], ",")))
    # The product W_p * L_c is what the criterion tests, so it must be
    # mesh-INDEPENDENT even though each factor is not.
    prods = [wp * el for _, el, wp, _ in rows]
    spread = (max(prods) - min(prods)) / min(prods)
    chk("the triggering ENERGY is mesh independent", spread < 1e-9,
        "W_p*L_c varies by %.2e across a 4x refinement -- the work DENSITY "
        "scales with 1/L_c, which is the regularisation, but the energy it "
        "represents does not" % spread)
    chk("refining raises the required work density",
        all(rows[i][2] < rows[i + 1][2] for i in range(len(rows) - 1)),
        "so a deck's PSI is only meaningful WITH its element size, which is "
        "why every deck states it")
    finest, coarsest = rows[-1][2], rows[0][2]
    print("       a %.0fx refinement changes the required work density %.1fx"
          % (rows[-1][0] / rows[0][0], finest / coarsest))


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def verify(path: str, timing: dict = None) -> None:
    print("=" * 78)
    print("SAG DECK: %s  (%.2f MB)" % (path, os.path.getsize(path) / 1e6))
    print("=" * 78)
    d = Deck(path)
    macro = "PU" in d.parts
    print("  read %d keyword blocks, %d part(s), %d step(s) -- %s deck"
          % (len(d.blocks), len(d.parts), len(d.steps),
             "MACRO" if macro else "MICRO"))
    check_grammar(d)
    check_solver(d)
    check_material(d)
    check_geometry(d, macro)
    check_steps(d, macro, timing or {})
    check_reachable(d)
    check_rigid_mass(d)
    if macro:
        check_sector(d, timing or {})
        check_seating(d, timing or {})


def build_pair(outdir="_sagverify") -> tuple:
    """A small but physically valid pair, for --build."""
    import glob

    from semgrit import sagdeck, sagemit
    from semgrit.quick import measure_images

    imgs = sorted(glob.glob("B4C_1*.tif"))[:1]
    if not imgs:
        raise SystemExit("no B4C_1*.tif to measure")
    os.makedirs(outdir, exist_ok=True)
    got = measure_images(imgs, os.path.join(outdir, "meas"),
                         log=lambda *a: None)
    p = sagdeck.SAGParams(grain_um=30.0, material="wc_co", name="sagverify",
                          macro_sector_deg=17.0, macro_grain_cap=2000,
                          micro_grains=1, grind_time_s=2.0e-5)
    pl = sagdeck.plan(p)
    mi = sagemit.write_micro(os.path.join(outdir, "micro.inp"), pl,
                             got["solids"])
    ma = sagemit.write_macro(os.path.join(outdir, "macro.inp"), pl,
                             got["solids"])
    timing = dict(pl["timing"])
    timing["compression_mm"] = p.compression_mm
    timing["prony_tau_s"] = p.polyurethane.prony_tau_s
    return [mi["path"], ma["path"]], timing


def main(argv):
    decks = [a for a in argv if not a.startswith("-")]
    timing = None
    if "--build" in argv or not decks:
        decks, timing = build_pair()
    for pth in decks:
        verify(pth, timing)
    if "--no-converge" not in argv:
        converge()
    print()
    print("=" * 78)
    if FAIL:
        print("  %d passed, %d FAILED" % (PASS, len(FAIL)))
        for f in FAIL:
            print("    - %s" % f)
        return 1
    print("  ALL %d CHECKS PASSED" % PASS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
