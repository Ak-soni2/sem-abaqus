"""Turn one hybrid deck into the ablation set the model's claim rests on.

The whole argument of this project is "the hybrid law predicts something that
neither pure law predicts". That is only demonstrable against the two pure laws
run on the SAME geometry, and nothing shipped ever produced them -- the three
packages compared a single-grit geometric deck, a multi-grit field deck and an
energy deck, which confounds grit count and h source with the criterion.

``PROPS(56)`` already selects the source:

    0  h from the coordinates      1  h from field variable 1
    2  FORCE DUCTILE everywhere    3  FORCE BRITTLE everywhere

and mode 3 is proven bit-identical to ``vumat_jh2.for`` -- exact equality over
21,000 increments across five strain histories, not a tolerance. So the ablation
is a ONE-TOKEN edit of an existing deck: rewrite the last value of the last
``*User Material`` data line. No rebuild, no re-mesh, no re-sweep. Everything
else in the file stays byte for byte, which is exactly what makes the comparison
attributable to the switch.

    from semgrit_multi.ablate import write_arms
    write_arms("multi_abrasive_field.inp", "ARMS")

Not ``swmode.set_energy_mode``: that rewrites PROPS(58)/(59), which is a
different question (which CRITERION), not this one (which LAW).

ONE CAVEAT, AND IT IS LOAD-BEARING. On a FIELD-CARRYING deck the forced arms
(2 and 3) currently read the field as a strength multiplier as well. The
heterogeneity guard at ``vumat_grind.for`` reads

    if (nfieldv .ge. 1 .and. ihmode .ne. 1) then
      if (fieldOld(km,1) .gt. z0) het = fieldOld(km,1)

which excludes only mode 1, so modes 2 and 3 fall through and multiply the JH-2
intact strength by a chip thickness in mm -- order 1e-4. Measured through the
driver on the shipped card: sigis 0.0505 -> 5.40e-06, a 9355x collapse, silent.
The comment three lines above that test says the feature is "not available when
field variable 1 is being used to carry the chip thickness", so the code does not
match its own stated intent.

This is a question about the constitutive routine and it has NOT been changed
here. Until it is resolved, the FORCED-BRITTLE arm on a field-carrying deck is
not a pure-JH-2 baseline -- it is JH-2 with the strength scaled by h. Either
resolve the guard first, or build the forced arms from a deck with no
``*Initial Conditions`` field (``h_source=0``), where the question does not
arise. The geometric-hybrid arm is unaffected.
"""

from __future__ import annotations

import os
import re
from typing import Optional

H_SOURCE_PROP = 56

ARMS = {
    1: ("geometric_hybrid",
        "the hybrid law, h from the swept field. The deck as shipped."),
    2: ("forced_ductile",
        "Johnson-Cook + SGE everywhere. The switch is disabled; every point "
        "is ductile regardless of its chip thickness."),
    3: ("forced_brittle",
        "Johnson-Holmquist II everywhere. Bit-identical to vumat_jh2.for, so "
        "this arm is the published brittle model on this geometry."),
}

_USER_MAT = re.compile(r"^(\*User Material,\s*constants\s*=\s*)(\d+)\s*$", re.I)


class AblateError(RuntimeError):
    pass


def set_h_source(deck_in: str, deck_out: str, h_source: int, *,
                 comment: Optional[list] = None) -> dict:
    """Copy a deck with ``PROPS(56)`` set to ``h_source``. One value changes."""
    if h_source not in (0, 1, 2, 3):
        raise AblateError("h_source must be 0, 1, 2 or 3, not %r" % (h_source,))
    if not os.path.exists(deck_in):
        raise AblateError("no such deck: %s" % deck_in)
    with open(deck_in, encoding="ascii") as fh:
        lines = fh.readlines()

    i_mat = n_declared = None
    for i, ln in enumerate(lines):
        m = _USER_MAT.match(ln.rstrip("\n"))
        if m:
            if i_mat is not None:
                raise AblateError("more than one *User Material block")
            i_mat, n_declared = i, int(m.group(2))
    if i_mat is None:
        raise AblateError("no *User Material block: this is not a VUMAT deck")

    j = i_mat + 1
    vals: list[str] = []
    while j < len(lines) and not lines[j].startswith("*"):
        vals += [x.strip() for x in lines[j].split(",") if x.strip()]
        j += 1
    if len(vals) != n_declared:
        raise AblateError("the card declares %d constants but %d are written"
                          % (n_declared, len(vals)))
    if n_declared < H_SOURCE_PROP:
        raise AblateError("the card has only %d constants, so it has no "
                          "PROPS(%d) to set" % (n_declared, H_SOURCE_PROP))

    was = int(round(float(vals[H_SOURCE_PROP - 1])))
    vals[H_SOURCE_PROP - 1] = repr(float(h_source))

    # Repack EIGHT to a line. Abaqus reads *User Material constants eight to a
    # line and rejects anything else outright; re-emitting rather than editing
    # in place means the rule holds whatever the constant count becomes.
    block = ["*User Material, constants=%d\n" % n_declared]
    block += [", ".join(vals[k:k + 8]) + "\n" for k in range(0, len(vals), 8)]

    name, why = ARMS.get(h_source, ("h_source_%d" % h_source, ""))
    note = ["**",
            "** ---------------- ABLATION ARM: %s ----------------" % name,
            "** PROPS(%d) h source %d -> %d." % (H_SOURCE_PROP, was, h_source),
            "** %s" % why,
            "**",
            "** ONE VALUE differs from the deck this was copied from. Same mesh,",
            "** same grains, same seating, same field, same everything else --",
            "** so any difference in the result is the constitutive law and",
            "** nothing else. That is what makes this an ablation rather than a",
            "** second model.",
            "**"]
    if h_source in (2, 3):
        note[-1:] = [
            "** WARNING, read before using this arm as a baseline.",
            "** vumat_grind.for's heterogeneity guard reads",
            "**     nfieldv .ge. 1 .and. ihmode .ne. 1",
            "** which excludes only mode 1. This deck CARRIES a field, so in",
            "** mode %d that guard passes and the JH-2 intact strength is" % h_source,
            "** multiplied by the chip thickness in mm -- order 1e-4, measured",
            "** at 9355x on the shipped card. So this arm is NOT the pure law;",
            "** it is the pure law with the strength scaled by h.",
            "**",
            "** This is a constitutive question and has deliberately not been",
            "** changed. To get a clean baseline, build the forced arms from a",
            "** deck with no *Initial Conditions field (h_source = 0).",
            "**"]
    for ln in (comment or []):
        note.insert(-1, "** " + ln)

    out = (lines[:i_mat] + [x + "\n" for x in note] + block + lines[j:])
    with open(deck_out, "w", encoding="ascii", newline="\n") as fh:
        fh.writelines(out)

    return {"path": deck_out, "arm": name, "why": why,
            "h_source": h_source, "h_source_was": was,
            "n_props": n_declared,
            "size_bytes": os.path.getsize(deck_out),
            "subroutine": ("vumat_grind2.for" if n_declared > 57
                           else "vumat_grind.for")}


def write_arms(deck_in: str, outdir: str, *,
               arms=(1, 2, 3), cores: int = 8) -> dict:
    """Write one folder per ablation arm, each with its own run command."""
    os.makedirs(outdir, exist_ok=True)
    made = {}
    for h in arms:
        name, why = ARMS[h]
        folder = os.path.join(outdir, "%d_%s" % (h, name))
        os.makedirs(folder, exist_ok=True)
        job = name
        dst = os.path.join(folder, job + ".inp")
        info = set_h_source(deck_in, dst, h)
        base = ("abaqus job=%s input=%s user=%s double=both"
                % (job, os.path.basename(dst), info["subroutine"]))
        line = base + " cpus=%d interactive" % cores
        # datacheck is serial work: cpus=1 needs 5 licence tokens against the 12
        # that cpus=8 reserves, and it is not a second faster for spending them.
        check = base + " cpus=1 datacheck"
        with open(os.path.join(folder, "run.bat"), "w", newline="") as fh:
            fh.write("@echo off\r\ncd /d \"%~dp0\"\r\n")
            # `call`, because on Windows abaqus is abaqus.bat and
            # running one batch file from another without it transfers
            # control permanently -- the caller never resumes, so the
            # solve never runs and nothing says so.
            fh.write("call abaqus verify -user_exp\r\n")
            fh.write("call " + check + "\r\n")
            fh.write("if errorlevel 1 exit /b 1\r\n")
            fh.write("call " + line + "\r\n")
        with open(os.path.join(folder, "run.sh"), "w", newline="\n") as fh:
            fh.write("#!/bin/sh\nset -e\ncd \"$(dirname \"$0\")\"\n")
            fh.write("abaqus verify -user_exp\n")
            fh.write(check + "\n")
            fh.write(line + "\n")
        # The arms are the comparison the whole model rests on, and until now they
        # were the only decks shipped with no way to read their .odb -- the README
        # documented a post step that did not exist here.
        from semgrit.odbpost import write_odb_postprocess_script
        write_odb_postprocess_script(
            os.path.join(folder, job + "_postprocess_odb.py"))
        info["size_bytes"] = os.path.getsize(dst)
        info["command"] = line
        # The folder is NUMBERED ("3_forced_brittle"), not bare. The README was
        # generated from the arm name alone, so its copy-paste block cd'd into
        # directories that do not exist.
        info["folder"] = os.path.basename(folder)
        made[name] = info
    return made


def demo() -> None:
    """Self-check against a shipped deck: python -m semgrit_multi.ablate"""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    src = os.path.join(here, "RUN_ME", "2_multi_abrasive",
                       "multi_abrasive_field.inp")
    if not os.path.exists(src):
        print("no shipped deck to check against: %s" % src)
        return
    import tempfile
    tmp = tempfile.mkdtemp(prefix="ablate_")
    base = open(src, encoding="ascii").read().splitlines()
    for h in (1, 2, 3):
        out = os.path.join(tmp, "a%d.inp" % h)
        info = set_h_source(src, out, h)
        got = open(out, encoding="ascii").read().splitlines()
        # Outside the ** comments and the card block, the file must be identical.
        a = [x for x in base if not x.startswith("**")]
        b = [x for x in got if not x.startswith("**")]
        assert len(a) == len(b), "line count changed: %d -> %d" % (len(a), len(b))
        diff = [(x, y) for x, y in zip(a, b) if x != y]
        # Only card data lines may differ, and only in ONE value.
        for x, y in diff:
            xa = [v.strip() for v in x.split(",")]
            ya = [v.strip() for v in y.split(",")]
            assert len(xa) == len(ya), "a card line changed width"
            assert sum(1 for p, q in zip(xa, ya) if p != q) == 1, \
                "more than one value changed on a line"
        assert len(diff) <= 1, "more than one line changed: %d" % len(diff)
        # And the value that changed must be PROPS(56).
        from verify_hybrid_deck import parse_deck
        pd = parse_deck(out)
        assert int(round(pd["props"][H_SOURCE_PROP - 1])) == h
        ref = parse_deck(src)
        for k, (p, q) in enumerate(zip(ref["props"], pd["props"])):
            if k != H_SOURCE_PROP - 1:
                assert p == q, "PROPS(%d) changed too" % (k + 1)
        print("  arm %d %-18s one value differs, %d card lines changed"
              % (h, info["arm"], len(diff)))
    print("ablate: ok -- the arms differ from the source deck by exactly "
          "PROPS(%d)" % H_SOURCE_PROP)


if __name__ == "__main__":
    demo()
