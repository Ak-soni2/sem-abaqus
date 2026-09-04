"""Build the smallest Abaqus job that can answer the one open question.

    python _make_probe.py [material]        # default sandstone

Abaqus has never run any deck in this project. The only submission attempt died
at preprocessing on ``*User Material constants=56`` written four to a line, which
is fixed -- but nothing has reached increment 1, so every statement about this
model is a statement about code rather than about a result.

Two things are genuinely unknown, and both are settled by ONE eight-node job that
finishes in seconds:

1. **Does field variable 1 reach the VUMAT?** The multi-abrasive and energy decks
   carry the chip thickness as ``*Initial Conditions, type=FIELD, variable=1`` and
   read it at ``PROPS(56)=1``. If it does not arrive, ``hloc`` is 0, ``0 < dc``,
   and the deck runs 100% DUCTILE while still producing a plausible chip and a
   clean .sta. There is no output that distinguishes that from success, so it has
   to be probed directly: SDV14 must come back as the number that was injected.

2. **Does element deletion fire?** ``*Depvar, delete=12`` plus
   ``*Section Controls ELEMENT DELETION=YES`` plus the VUMAT setting SDV12 to
   zero is a three-way agreement no gate can check.

The probe writes eight one-element blocks, each with a different injected h that
straddles dc, so one job reports the branch the subroutine picks across the whole
transition. It uses the material card of a real shipped deck verbatim -- read out
of the .inp, not rebuilt -- so a card that Abaqus rejects here is a card Abaqus
would reject there.

    cd PROBE && run.bat          (or: abaqus job=probe input=probe.inp
                                       user=vumat_grind.for double=both cpus=1)

Then read probe.dat / probe.sta and PROBE/EXPECTED.md. Nothing in the roadmap has
a defined value until this returns.
"""

from __future__ import annotations

import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# One element per h, as a multiple of dc. Straddles the transition so a single
# job maps the whole switch, including the two values either side of it where a
# sign error or an off-by-one comparison would show up first.
H_OVER_DC = (0.10, 0.50, 0.90, 0.99, 1.01, 1.10, 2.00, 5.00)

# 0.3 um cube. The real decks use 0.3 x 1.5 x 0.03 um, but the probe is not
# measuring a force -- a cube keeps the stable increment obvious and makes the
# prescribed strain trivially readable.
L_MM = 3.0e-4


def read_card(inp_path: str):
    """The *User Material constants and *Depvar out of a real deck.

    Reads the file rather than rebuilding from semgrit on purpose: the point of
    the probe is to put the exact bytes Abaqus will see in front of Abaqus.
    """
    props: list[float] = []
    depvar = delete_sdv = None
    grabbing = False
    mode = None
    with open(inp_path, encoding="ascii") as fh:
        for ln in fh:
            if ln.startswith("**"):
                continue
            if ln.startswith("*"):
                low = ln.lower()
                grabbing = low.startswith("*user material")
                mode = "depvar" if low.startswith("*depvar") else None
                if mode == "depvar":
                    import re
                    m = re.search(r"delete\s*=\s*(\d+)", ln, re.I)
                    delete_sdv = int(m.group(1)) if m else None
                continue
            if grabbing and ln.strip():
                props += [float(x) for x in ln.split(",") if x.strip()]
            elif mode == "depvar" and ln.strip():
                depvar = int(float(ln.split(",")[0]))
                mode = None
    if not props:
        raise SystemExit("no *User Material found in %s" % inp_path)
    return props, depvar, delete_sdv


def predict(props, dc_mm):
    """Run the SAME card through the compiled VUMAT with the same field values.

    So the shipped EXPECTED.md carries numbers that have already been checked
    outside Abaqus, not numbers someone hopes are right. If Abaqus disagrees
    with this table the difference is Abaqus' field-variable plumbing, because
    everything else has been held identical.
    """
    import verify_vumat_grind as V
    drv = V.Driver(V.find_gfortran(), V.SRC_GRIND, "probechk.exe")
    seg = [(400, (0.0, -0.0005, 0.0, 0.0, 0.0, 0.0))]
    rows = []
    for f in H_OVER_DC:
        h = f * dc_mm
        s = drv.run(props, seg, nstatev=20, fields=(h,), coord=(0.0, 0.0, 0.0),
                    nout=2)[-1]["sdv"]
        rows.append({"f": f, "h": h, "sdv14": s[13], "sdv13": int(round(s[12])),
                     "sdv19": s[18], "err": abs(s[13] - h)})
    worst = max(r["err"] for r in rows)
    bad = [r for r in rows if r["sdv13"] != (1 if r["f"] < 1.0 else 2)]
    if worst > 0.0 or bad:
        raise SystemExit("the driver does not reproduce the probe's own "
                         "expectations: worst |dh| %.3g, %d branch mismatches"
                         % (worst, len(bad)))
    return rows


def write_probe(folder: str, props, depvar, delete_sdv, dc_mm, material,
                src_deck, rows=None):
    os.makedirs(folder, exist_ok=True)
    n = len(H_OVER_DC)
    path = os.path.join(folder, "probe.inp")
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        w = fh.write
        w("*Heading\n")
        w("** FIELD-VARIABLE PROBE for vumat_grind.for\n")
        w("** Eight C3D8R cubes of %g mm, one per injected chip thickness.\n"
          % L_MM)
        w("** Material card copied verbatim from %s\n"
          % os.path.basename(src_deck))
        w("** material: %s     dc = %.6e mm (%.4f nm)\n"
          % (material, dc_mm, dc_mm * 1e6))
        w("**\n")
        w("** WHAT THIS ANSWERS\n")
        w("**   SDV14 must come back as the h injected into each element.\n")
        w("**   SDV13 must be 1 (ductile) where h < dc and 2 (brittle) above.\n")
        w("**   If every SDV14 is zero, field variable 1 is NOT reaching the\n")
        w("**   VUMAT, and every field-carrying deck in this project is\n")
        w("**   running 100%% ductile while looking entirely healthy.\n")
        w("**\n")
        w("** element   h/dc        h [mm]\n")
        for i, f in enumerate(H_OVER_DC, start=1):
            w("**   %-8d %-11.2f %.9e\n" % (i, f, f * dc_mm))
        w("**\n")

        # ---- one part, n disjoint cubes ---------------------------------
        w("*Part, name=P\n")
        w("*Node\n")
        nid = 0
        for e in range(n):
            x0 = e * 2.0 * L_MM          # spaced apart: no shared nodes
            for (dx, dy, dz) in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                                 (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)):
                nid += 1
                w("%d, %.12e, %.12e, %.12e\n"
                  % (nid, x0 + dx * L_MM, dy * L_MM, dz * L_MM))
        w("*Element, type=C3D8R\n")
        for e in range(n):
            b = e * 8
            w("%d, %s\n" % (e + 1, ", ".join(str(b + k) for k in range(1, 9))))
        w("*Elset, elset=ALL, generate\n1, %d, 1\n" % n)
        w("*Solid Section, elset=ALL, material=PROBE, controls=EC-1\n,\n")
        w("*End Part\n")
        w("**\n*Assembly, name=A\n*Instance, name=P-1, part=P\n*End Instance\n")

        # Node sets, one per element, so the field can differ element to element.
        for e in range(n):
            b = e * 8
            w("*Nset, nset=NS%d, instance=P-1\n" % (e + 1))
            w("%s\n" % ", ".join(str(b + k) for k in range(1, 9)))
        # Fully constrained base, prescribed compression on the top face.
        w("*Nset, nset=BOT, instance=P-1\n")
        w("%s\n" % ", ".join(str(e * 8 + k) for e in range(n)
                             for k in (1, 2, 3, 4)))
        w("*Nset, nset=TOP, instance=P-1\n")
        w("%s\n" % ", ".join(str(e * 8 + k) for e in range(n)
                             for k in (5, 6, 7, 8)))
        w("*End Assembly\n")

        w("**\n*Material, name=PROBE\n")
        w("*Density\n%.8e,\n" % props[29])
        if delete_sdv:
            w("*Depvar, delete=%d\n%d,\n" % (delete_sdv, depvar))
        else:
            w("*Depvar\n%d,\n" % depvar)
        # EIGHT per line. Same rule that killed the first submission.
        w("*User Material, constants=%d\n" % len(props))
        for i in range(0, len(props), 8):
            w(", ".join(repr(float(v)) for v in props[i:i + 8]) + "\n")
        w("**\n*Section Controls, name=EC-1, hourglass=ENHANCED, "
          "ELEMENT DELETION=YES\n1., 1., 1.\n")

        # ---- the field: one value per element, on its own node set ------
        w("**\n** The whole point. variable=1 is what PROPS(56)=1 reads.\n")
        # repr(), not %.12e. The value being injected is compared against dc a
        # few nanometres away, and %g-style formatting does not round-trip a
        # double -- which is the same trap the material card was bitten by, and
        # which the assert at the bottom of this file caught here. The
        # production writer (semgrit_multi/fieldinject) already uses %r.
        w("*Initial Conditions, type=FIELD, variable=1\n")
        for e, f in enumerate(H_OVER_DC, start=1):
            w("NS%d, %r\n" % (e, float(f * dc_mm)))

        # ---- one short explicit step, prescribed compression ------------
        # 20% nominal strain at 1e6 /s: enough plastic strain for the brittle
        # branch to damage and the ductile branch to harden, short enough to
        # finish in seconds.
        step_t = 2.0e-7
        vel = -0.20 * L_MM / step_t
        w("**\n*Step, name=S1, nlgeom=YES\n")
        w("*Dynamic, Explicit\n, %.9e\n" % step_t)
        w("*Bulk Viscosity\n0.06, 1.2\n")
        w("*Boundary\nBOT, 1, 3\n")
        w("*Boundary, type=VELOCITY\nTOP, 2, 2, %.9e\nTOP, 1, 1\nTOP, 3, 3\n"
          % vel)
        w("*Output, field, number interval=20\n")
        w("*Node Output\nU, V\n")
        w("*Element Output, directions=YES\nS, SDV, STATUS\n")
        w("*Output, history, time interval=%.9e\n" % (step_t / 200.0))
        w("*Energy Output\nALLIE, ALLKE, ALLAE, ALLPD\n")
        w("*End Step\n")

    # ---- how to run it, and what to accept -------------------------------
    # The subroutine has to sit BESIDE the deck: Abaqus resolves user= against
    # the working directory, so a probe folder without it aborts at compile --
    # and the probe is the gate that proves the field variable reaches the
    # material points at all, so it failing silently costs two 5 h runs.
    shutil.copy(os.path.join(HERE, "vumat_grind.for"), folder)

    line = ("abaqus job=probe input=probe.inp user=vumat_grind.for "
            "double=both cpus=1")
    with open(os.path.join(folder, "run.bat"), "w", newline="") as fh:
        fh.write("@echo off\r\n")
        fh.write("cd /d \"%~dp0\"\r\n")
        fh.write("rem Step 0: does the Fortran toolchain work at all?\r\n")
        fh.write("abaqus verify -user_exp\r\n")
        fh.write("rem Step 1: preprocessing only. Seconds. Reads the card.\r\n")
        fh.write(line + " datacheck\r\n")
        fh.write("if errorlevel 1 exit /b 1\r\n")
        fh.write("rem Step 2: solve it.\r\n")
        fh.write(line + " interactive\r\n")
    with open(os.path.join(folder, "run.sh"), "w", newline="\n") as fh:
        fh.write("#!/bin/sh\nset -e\ncd \"$(dirname \"$0\")\"\n")
        fh.write("abaqus verify -user_exp\n")
        fh.write(line + " datacheck\n" + line + " interactive\n")

    exp = ["# What to accept", "",
           "Material: **%s**, dc = **%.4f nm**, card copied from `%s`."
           % (material, dc_mm * 1e6, os.path.basename(src_deck)), "",
           "Every number below was produced by running **this card with these",
           "field values** through the compiled VUMAT outside Abaqus, so the",
           "table is a verified prediction rather than a hope. Any disagreement",
           "is therefore Abaqus' field-variable plumbing and nothing else.", "",
           "| element | h/dc | h [mm] | SDV14 must be | SDV13 must be | SDV19 |",
           "|---|---|---|---|---|---|"]
    for i, f in enumerate(H_OVER_DC, start=1):
        r = (rows or [{}])[i - 1] if rows else {}
        exp.append("| %d | %.2f | %.6e | %.6e | %d (%s) | %s |"
                   % (i, f, f * dc_mm, f * dc_mm, 1 if f < 1.0 else 2,
                      "ductile" if f < 1.0 else "brittle",
                      ("%.4f" % r["sdv19"]) if r else "-"))
    if rows:
        exp += ["",
                "SDV19 is the strain-gradient amplification. It is above 1 in "
                "the ductile",
                "branch and rises as h falls (%.4f at h = %.2f dc, %.4f at "
                "h = %.2f dc);"
                % (rows[0]["sdv19"], rows[0]["f"], rows[3]["sdv19"],
                   rows[3]["f"]),
                "it is exactly 0 in the brittle branch, which never evaluates "
                "it."]
    exp += ["",
            "## The three outcomes",
            "",
            "**All eight SDV14 equal the injected h, and SDV13 flips between",
            "element 4 and element 5.** The field route works. Every",
            "field-carrying deck in the project is a genuine hybrid deck, and",
            "the roadmap's Phase 2 onwards is aimed at code that executes.",
            "",
            "**Every SDV14 is zero.** Field variable 1 is not arriving. Then",
            "`hloc` is 0, `0 < dc`, and `RUN_ME*/2_multi_abrasive` and",
            "`3_energy_criterion` are running **100% ductile** -- with a",
            "plausible chip, a clean `.sta` and nothing in any output to say",
            "so. Fix that before spending 5.3 h per material on anything.",
            "",
            "**The job dies at preprocessing.** Read the `.dat`. Compare",
            "against `error/single_abrasive.dat`, which is the same failure",
            "this project already had once (`*User Material` four values to a",
            "line instead of eight).",
            "",
            "## Also worth reading off the same job",
            "",
            "* `STATUS` -- if nothing ever deletes, the three-way agreement",
            "  between `*Depvar, delete=%d`, `ELEMENT DELETION=YES` and the"
            % (delete_sdv or 12),
            "  VUMAT zeroing SDV%d is broken, and no chip will ever separate"
            % (delete_sdv or 12),
            "  in the real decks either.",
            "* `SDV13` on elements 4 and 5 (h = 0.99 dc and 1.01 dc) -- an",
            "  off-by-one or a `<=` where a `<` belongs shows up only there.",
            "* `SDV19`, the strain-gradient amplification. At h = 0.1 dc it",
            "  should be well above 1; at h = 5 dc close to 1.",
            "* `ALLAE/ALLIE` -- hourglass energy on a single-element test",
            "  should be negligible. If it is not, the enhanced hourglass",
            "  control is not doing its job at this aspect ratio.",
            ""]
    with open(os.path.join(folder, "EXPECTED.md"), "w", encoding="utf-8",
              newline="\n") as fh:
        fh.write("\n".join(exp))
    return path, line


def main(argv) -> int:
    material = argv[1] if len(argv) > 1 else "sandstone"
    src = {"sandstone": os.path.join(HERE, "RUN_ME", "2_multi_abrasive",
                                     "multi_abrasive_field.inp"),
           "silicon_carbide": os.path.join(HERE, "RUN_ME_SIC",
                                           "2_multi_abrasive",
                                           "multi_abrasive_field.inp")}[material]
    if not os.path.exists(src):
        print("no such deck: %s  (run _make_run_packages.py first)" % src)
        return 2
    props, depvar, delete_sdv = read_card(src)
    dc_mm = props[46]
    folder = os.path.join(HERE, "PROBE_" + material)
    shutil.rmtree(folder, ignore_errors=True)
    rows = predict(props, dc_mm)
    path, line = write_probe(folder, props, depvar, delete_sdv, dc_mm,
                             material, src, rows)

    # The card must be the deck's card, byte for byte in value.
    from verify_hybrid_deck import (check_user_material_format, parse_deck)
    d = parse_deck(path)
    assert len(d["props"]) == len(props), "prop count changed"
    assert max(abs(a - b) for a, b in zip(d["props"], props)) == 0.0, \
        "the probe card does not round-trip to the deck card"
    fv = d["field_vals"]
    assert len(fv) == len(H_OVER_DC), "expected %d field values, got %d" % (
        len(H_OVER_DC), len(fv))
    want = [f * dc_mm for f in H_OVER_DC]
    assert max(abs(a - b) for a, b in zip(sorted(fv), sorted(want))) == 0.0, \
        "the injected field does not round-trip"
    check_user_material_format("probe *User Material is 8 per line", path)

    print("=" * 78)
    print("PROBE_%s/  %s  (%d constants, %d SDVs, delete=%s)"
          % (material, os.path.basename(path), len(props), depvar, delete_sdv))
    print("dc = %.6e mm (%.4f nm); %d elements, h/dc = %s"
          % (dc_mm, dc_mm * 1e6, len(H_OVER_DC),
             ", ".join("%g" % f for f in H_OVER_DC)))
    print("card and field both round-trip out of the written file EXACTLY")
    print("compiled VUMAT reproduces all %d expectations: SDV14 exact to 0, "
          "SDV13 flips" % len(rows))
    print("between h/dc = %.2f and %.2f, SDV19 %.4f -> 1 -> 0 across the "
          "transition" % (H_OVER_DC[3], H_OVER_DC[4], rows[0]["sdv19"]))
    print("=" * 78)
    print("  cd PROBE_%s" % material)
    print("  %s datacheck      <- seconds, reads the card" % line)
    print("  %s interactive    <- solves it" % line)
    print()
    print("Then read PROBE_%s/EXPECTED.md. Accept only if SDV14 equals the"
          % material)
    print("injected h and SDV13 flips between elements 4 and 5.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
