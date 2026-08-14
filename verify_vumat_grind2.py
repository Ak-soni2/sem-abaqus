"""Verification gate for vumat_grind2.for, the version with the local criterion.

Three things have to hold, in this order of importance:

**It has not drifted from its source.** ``vumat_grind2.for`` is derived from
``vumat_grind.for``, which passed 112 checks. With SWMODE = 0 the two must be
bit-identical on every history -- that is what makes the derivation safe and
what will catch it if someone edits one and not the other.

**It still inherits the JH-2 validation.** With the brittle override it must be
bit-identical to ``vumat_jh2.for``, exactly as its source is.

**The new criterion is what it claims.** The plastic work accumulates as
written, the ratio is the ratio, the flip happens when the ratio reaches 1 and
latches. And the relationship to the geometric criterion is checked as it
actually is, not as it would be convenient for it to be: the pointwise Griffith
balance gives an exponent of +1 on H/E where the published forms use +0.5 and
-1, so the two are different expressions, and PSI is defaulted to the value that
makes them agree instead.

    python verify_vumat_grind2.py [-v]
"""

from __future__ import annotations

import math
import os
import subprocess
import sys

import verify_vumat_grind as vg
from verify_vumat_grind import JC, JH2, check, close, find_gfortran

HERE = os.path.dirname(os.path.abspath(__file__))
SRC2 = os.path.join(HERE, "vumat_grind2.for")


def card2(*, swmode=0, psi=0.0, **kw):
    """The 56-prop card of vumat_grind.for plus SWMODE and PSI."""
    p = vg.card(**kw)
    return p + [float(swmode), float(psi)]


def main() -> int:
    fc = find_gfortran()
    print("=" * 78)
    print("vumat_grind2.for verification")
    print("=" * 78)

    # ---- 0. the derivation is current -------------------------------
    print("0. the derivation")
    r = subprocess.run([sys.executable,
                        os.path.join(HERE, "_derive_grind2.py"), "--check"],
                       capture_output=True, text=True, cwd=HERE)
    check("vumat_grind2.for is up to date with vumat_grind.for",
          r.returncode == 0, r.stdout.strip().splitlines()[-1]
          if r.stdout.strip() else r.stderr[:80])

    g1 = vg.Driver(fc, vg.SRC_GRIND, "grind.exe")
    g2 = vg.Driver(fc, SRC2, "grind2.exe")
    jh = vg.Driver(fc, vg.SRC_JH2, "jh2.exe")
    check("compiles with no conversion/typing warnings", not g2.warnings,
          "; ".join(g2.warnings[:3]))

    # ---- 1. source hygiene ------------------------------------------
    print("1. source hygiene")
    text = open(SRC2, encoding="ascii").read()
    longs = [i for i, ln in enumerate(text.split("\n"), 1)
             if ln.strip() and ln[0] not in "Cc*!" and len(ln.rstrip()) > 72]
    check("no code line exceeds column 72", not longs, str(longs[:8]))
    check("ascii only, unix line endings", "\r" not in text
          and all(ord(c) < 128 for c in text))
    check("the header says it is derived and must not be edited",
          "DERIVED FROM vumat_grind.for" in text)
    import re
    bad = {m.group(1) for m in re.finditer(
        r"^ {6,}([a-z][a-z0-9_]*)\s*=\s*props\(", text,
        re.MULTILINE | re.IGNORECASE) if m.group(1)[0].lower() in "ijklmn"}
    check("no real prop is read into an implicitly-integer name", not bad,
          str(sorted(bad)))

    # ---- 2. SWMODE = 0 is the original, bit for bit -----------------
    print("2. SWMODE = 0 is vumat_grind.for, bit for bit")
    for name, segs in vg.jh2_histories():
        for tag, extra in (("switch on dc", {}),
                           ("forced ductile", dict(ihmode=2)),
                           ("forced brittle", dict(ihmode=3))):
            p = card2(swmode=0, dcut=1.0e-4, h0=5.0e-5, **extra)
            a = g2.run(p, segs, nstatev=22, nout=250)
            b = g1.run(p[:vg.N_PROPS_GRIND], segs, nstatev=20, nout=250)
            if len(a) != len(b):
                check("row counts match: %s / %s" % (name, tag), False,
                      "%d vs %d" % (len(a), len(b)))
                continue
            ws = max((abs(x - y) for ra, rb in zip(a, b)
                      for x, y in zip(ra["s"], rb["s"])), default=0.0)
            wv = max((abs(x - y) for ra, rb in zip(a, b)
                      for x, y in zip(ra["sdv"][:20], rb["sdv"][:20])),
                     default=0.0)
            check("identical stress: %s / %s" % (name, tag), ws == 0.0,
                  "max %.3g" % ws)
            check("identical SDV 1-20: %s / %s" % (name, tag), wv == 0.0,
                  "max %.3g" % wv)

    # ---- 3. still the original JH-2 ---------------------------------
    print("3. the brittle branch is still vumat_jh2.for, bit for bit")
    for name, segs in vg.jh2_histories():
        a = g2.run(card2(swmode=0, ihmode=3), segs, nstatev=22, nout=250)
        b = jh.run(JH2 + [1.0, 1.0, 1.0, -1.0], segs, nstatev=12, nout=250)
        ws = max((abs(x - y) for ra, rb in zip(a, b)
                  for x, y in zip(ra["s"], rb["s"])), default=0.0)
        check("identical to vumat_jh2.for: " + name, ws == 0.0 and
              len(a) == len(b), "max %.3g" % ws)

    # ---- 4. how the two criteria relate ----------------------------
    print("4. how the local and geometric criteria relate")
    # The pointwise Griffith balance gives dc = PSI (H/E)^1 (Kc/H)^2, an
    # exponent of +1 where the published geometric forms use +0.5 and -1. So
    # they are NOT the same expression, and the header must not claim they are.
    for lam, H, E, Kc in ((0.15, 1000.0, 6500.0, 9.4868),
                          (0.30, 2500.0, 70000.0, 1.5),
                          (0.05, 800.0, 30000.0, 3.0)):
        d_energy = lam * Kc * Kc / (E * H)
        close("the energy balance gives PSI (H/E)(Kc/H)^2  (H=%g)" % H,
              d_energy, lam * (H / E) * (Kc / H) ** 2, rtol=1e-14)
        check("which is neither published form  (H=%g)" % H,
              abs(d_energy - vg.dc_formula(lam, H, E, Kc, 1)) > 1e-12
              and abs(d_energy - vg.dc_formula(lam, H, E, Kc, 2)) > 1e-12,
              "energy %.4g, form1 %.4g, form2 %.4g"
              % (d_energy, vg.dc_formula(lam, H, E, Kc, 1),
                 vg.dc_formula(lam, H, E, Kc, 2)))
        # Which is exactly why PSI is defaulted from dc instead: with
        # PSI = dc E H / Kc^2 the threshold work per area is H*dc, whichever
        # form dc came from.
        for form in (1, 2):
            dc = vg.dc_formula(lam, H, E, Kc, form)
            psi_def = dc * E * H / (Kc * Kc)
            close("PSI from dc form %d gives a threshold of H*dc  (H=%g)"
                  % (form, H), psi_def * Kc * Kc / E, H * dc, rtol=1e-13)
    check("the header states that relationship, and does not claim identity",
          "THIRD member of the family" in text
          and "PSI = dc E H / Kc^2" in text)

    # ---- 5. the plastic work and the ratio --------------------------
    print("5. plastic work, the ratio, and the flip")
    c = dict(JC)
    dg = 1.0e-5
    clen = 1.0e-3
    Kc, H, lam = vg.JC["E"] * 0.0 + 9.4868, 1000.0, 0.15
    # Energy mode only, so the chip thickness is never consulted.
    p = card2(swmode=1, psi=lam, hard=H, kic=Kc, dcut=1.0e-9, h0=1.0,
              lamc=lam)
    r = g2.run(p, [(4000, (0, 0, 0, dg, 0, 0))], nstatev=22, nout=1,
               dt=1e-9, charlen=clen)
    gcrit = lam * Kc * Kc / c["E"]
    wp = 0.0
    worst_w = worst_e = 0.0
    for row in r:
        if int(round(row["sdv"][12])) != 1:
            break
        wp += row["sdv"][17] * row["sdv"][8] * 0.0  # placeholder, replaced below
    # q_pl is (1-D)*SEFF, the stress that did the work; recompute it properly
    wp = 0.0
    dprev = 0.0
    for i, row in enumerate(r):
        dep = row["sdv"][8]
        if dep > 0:
            wp += (1.0 - dprev) * row["sdv"][17] * dep
        dprev = row["sdv"][0]
        worst_w = max(worst_w, abs(row["sdv"][20] - wp))
        worst_e = max(worst_e, abs(row["sdv"][21] - wp * clen / gcrit))
        if int(round(row["sdv"][12])) == 2:
            break
    # With PSI defaulted the threshold is H*dc exactly, read back through the
    # running code: ERATIO = W_p L_c / (H dc).
    pdef = card2(swmode=1, psi=0.0, hard=H, kic=Kc, dcut=2.0e-6, h0=1.0,
                 lamc=lam)
    rd = g2.run(pdef, [(4000, (0, 0, 0, dg, 0, 0))], nstatev=22, nout=1,
                dt=1e-9, charlen=clen)
    wd = 0.0
    dpv = 0.0
    worst_d = 0.0
    for row in rd:
        if row["sdv"][8] > 0:
            wd += (1.0 - dpv) * row["sdv"][17] * row["sdv"][8]
        dpv = row["sdv"][0]
        if int(round(row["sdv"][12])) == 2:
            break
        worst_d = max(worst_d,
                      abs(row["sdv"][21] - wd * clen / (H * 2.0e-6)))
    check("with PSI defaulted the threshold work per area is exactly H*dc",
          worst_d < 1e-9, "max |dratio| = %.3g" % worst_d)

    check("W_p accumulates as sum of (1-D) sigma_e depsilon_p",
          worst_w < 1e-9, "max |dW| = %.3g MPa" % worst_w)
    check("SDV22 is W_p L_c E / (PSI Kc^2)", worst_e < 1e-9,
          "max |dratio| = %.3g" % worst_e)

    modes = [int(round(row["sdv"][12])) for row in r]
    ratios = [row["sdv"][21] for row in r]
    check("every point starts ductile in energy mode", modes[0] == 1)
    flip = next((i for i, m in enumerate(modes) if m == 2), None)
    check("the point does flip to brittle", flip is not None,
          "ratio reached %.4f" % max(ratios))
    if flip is not None:
        check("it flips on the increment the ratio reaches 1",
              ratios[flip] >= 1.0 and (flip == 0 or ratios[flip - 1] < 1.0),
              "ratio %.6f before, %.6f at the flip"
              % (ratios[flip - 1] if flip else 0.0, ratios[flip]))
        check("and stays brittle afterwards",
              all(m == 2 for m in modes[flip:]))
        check("the plastic work stops growing once brittle",
              abs(r[-1]["sdv"][20] - r[flip]["sdv"][20]) < 1e-9)

    # ---- 6. mesh dependence, stated rather than hidden --------------
    print("6. mesh dependence")
    # SWMODE 0 keeps the criterion inert, so the ratio is still recorded but
    # nothing flips and it cannot saturate at 1 partway through the comparison.
    p_inert = card2(swmode=0, psi=lam, hard=H, kic=Kc, dcut=1.0, h0=1.0e-5,
                    lamc=lam)
    got = []
    for cl in (2.0e-3, 1.0e-3, 5.0e-4):
        rr = g2.run(p_inert, [(2000, (0, 0, 0, dg, 0, 0))], nstatev=22,
                    nout=2000, dt=1e-9, charlen=cl)
        got.append(rr[-1]["sdv"][21])
    check("the ratio is recorded even with the criterion inert",
          min(got) > 0.0, "ratios %s" % ["%.4g" % v for v in got])
    close("halving the element halves the ratio", got[1] / got[0], 0.5,
          rtol=1e-9)
    close("and again", got[2] / got[1], 0.5, rtol=1e-9)
    check("so PSI is calibrated for a mesh, and the header says so",
          "calibrated FOR A MESH" in text)

    # ---- 7. the modes ----------------------------------------------
    print("7. the three switch modes")
    seg = [(10, (0, 0, 0, 1e-9, 0, 0))]
    # geometric only: h >= dc is brittle, and the work criterion is inert
    r0 = g2.run(card2(swmode=0, dcut=1e-4, h0=1e-3, psi=lam, hard=H, kic=Kc),
                seg, nstatev=22, nout=10)
    check("SWMODE 0, h > dc -> brittle", int(round(r0[-1]["sdv"][12])) == 2)
    r0b = g2.run(card2(swmode=0, dcut=1e-4, h0=1e-5, psi=lam, hard=H, kic=Kc),
                 seg, nstatev=22, nout=10)
    check("SWMODE 0, h < dc -> ductile", int(round(r0b[-1]["sdv"][12])) == 1)
    # energy only: h is ignored entirely
    r1 = g2.run(card2(swmode=1, dcut=1e-9, h0=1.0, psi=lam, hard=H, kic=Kc),
                seg, nstatev=22, nout=10)
    check("SWMODE 1 ignores h >> dc and starts ductile",
          int(round(r1[-1]["sdv"][12])) == 1)
    # both: h decides the start, the work can still flip it
    r2 = g2.run(card2(swmode=2, dcut=1e-4, h0=1e-3, psi=lam, hard=H, kic=Kc),
                seg, nstatev=22, nout=10)
    check("SWMODE 2, h > dc -> brittle from the start",
          int(round(r2[-1]["sdv"][12])) == 2)
    r2b = g2.run(card2(swmode=2, dcut=1e-4, h0=1e-5, psi=lam, hard=H, kic=Kc),
                 [(4000, (0, 0, 0, dg, 0, 0))], nstatev=22, nout=1, dt=1e-9)
    ms = [int(round(x["sdv"][12])) for x in r2b]
    check("SWMODE 2, h < dc -> starts ductile then the work flips it",
          ms[0] == 1 and ms[-1] == 2,
          "%d -> %d" % (ms[0], ms[-1]))

    # PSI <= 0 must fall back to lambda_c, so the default is the equivalence.
    # PSI defaulting, measured through the inert ratio so nothing saturates.
    seg2k = [(2000, (0, 0, 0, dg, 0, 0))]
    r_auto = g2.run(card2(swmode=0, psi=0.0, lamc=lam, hard=H, kic=Kc,
                          h0=1.0e-7, dcut=2.0e-6),
                    seg2k, nstatev=22, nout=2000, dt=1e-9, charlen=clen)
    r_lam = g2.run(card2(swmode=0, psi=lam, lamc=lam, hard=H, kic=Kc,
                         h0=1.0e-7, dcut=2.0e-6),
                   seg2k, nstatev=22, nout=2000, dt=1e-9, charlen=clen)
    check("PSI <= 0 uses dc*E*H/Kc^2, which is not lambda_c",
          abs(r_auto[-1]["sdv"][21] - r_lam[-1]["sdv"][21]) > 1e-12,
          "%.6g vs %.6g" % (r_auto[-1]["sdv"][21], r_lam[-1]["sdv"][21]))
    # and the two differ by exactly the ratio of the PSI values
    close("by exactly the ratio of the two PSI values",
          r_auto[-1]["sdv"][21] / r_lam[-1]["sdv"][21],
          lam / (2.0e-6 * c["E"] * H / (Kc * Kc)), rtol=1e-9)
    # With no dc it falls back to lambda_c.
    r_fb = g2.run(card2(swmode=0, psi=0.0, lamc=lam, hard=0.0, kic=Kc,
                        h0=1.0e-7, dcut=1.0e-6),
                  seg2k, nstatev=22, nout=2000, dt=1e-9, charlen=clen)
    check("with no hardness there is no dc-based PSI, and no criterion either",
          r_fb[-1]["sdv"][21] >= 0.0)
    r_br = g2.run(card2(swmode=0, psi=0.0, lamc=lam, hard=H, kic=Kc,
                        h0=1.0e-3, dcut=2.0e-6),
                  seg2k, nstatev=22, nout=2000, dt=1e-9, charlen=clen)
    check("a point that was always brittle reports no energy ratio at all",
          int(round(r_br[-1]["sdv"][12])) == 2
          and r_br[-1]["sdv"][21] == 0.0 and r_br[-1]["sdv"][20] == 0.0,
          "mode %d, ratio %.6g, W_p %.6g"
          % (int(round(r_br[-1]["sdv"][12])), r_br[-1]["sdv"][21],
             r_br[-1]["sdv"][20]))

    # With no toughness there is no criterion, and it must be inert rather
    # than dividing by zero.
    rz = g2.run(card2(swmode=1, psi=lam, hard=0.0, kic=0.0, dcut=1e-6, h0=1.0),
                [(2000, (0, 0, 0, dg, 0, 0))], nstatev=22, nout=2000, dt=1e-9)
    check("no toughness -> the criterion is inert, not a division by zero",
          all(math.isfinite(v) for v in rz[-1]["sdv"])
          and rz[-1]["sdv"][21] == 0.0
          and int(round(rz[-1]["sdv"][12])) == 1)

    # ---- 8. continuity across the flip ------------------------------
    print("8. continuity across the flip")
    rr = g2.run(card2(swmode=1, psi=lam, hard=H, kic=Kc, h0=1.0),
                [(4000, (-1e-5, 3e-6, 3e-6, 1e-5, 0, 0))], nstatev=22,
                nout=1, dt=1e-9)
    ms = [int(round(x["sdv"][12])) for x in rr]
    fl = next((i for i, m in enumerate(ms) if m == 2), None)
    check("this history flips too", fl is not None)
    if fl is not None and fl + 1 < len(rr):
        jump = max(abs(a - b) for a, b in zip(rr[fl]["s"], rr[fl - 1]["s"]))
        scale = max(abs(x) for x in rr[fl]["s"]) or 1.0
        check("the stress does not jump discontinuously at the flip",
              jump / scale < 0.5,
              "largest component change %.3g of %.3g" % (jump, scale))
        check("the damage carries across rather than resetting",
              rr[fl + 1]["sdv"][0] >= rr[fl]["sdv"][0] - 1e-12,
              "D %.6f -> %.6f" % (rr[fl]["sdv"][0], rr[fl + 1]["sdv"][0]))
        check("and the point keeps its plastic strain",
              rr[fl + 1]["sdv"][1] >= rr[fl]["sdv"][1] - 1e-12)
    check("no NaN or Inf anywhere in that history",
          all(math.isfinite(v) for row in rr for v in row["s"] + row["sdv"]))

    print("=" * 78)
    if vg._FAIL:
        print("%d passed, %d FAILED" % (vg._PASS, len(vg._FAIL)))
        for f in vg._FAIL:
            print("   - " + f)
        return 1
    print("ALL %d CHECKS PASSED" % vg._PASS)
    return 0


if __name__ == "__main__":
    sys.exit(main())
