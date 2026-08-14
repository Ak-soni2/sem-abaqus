"""Switch a hybrid deck over to vumat_grind2.for's local energy criterion.

``semgrit.hybrid`` writes a 56-constant card with ``*Depvar 20``, which is what
``vumat_grind.for`` reads. ``vumat_grind2.for`` reads two more constants --
SWMODE and PSI -- and writes two more state variables. Rather than teach the
verified writer about a subroutine it does not know, this module rewrites the
card of an already-built deck, the same way :mod:`semgrit_multi.fieldinject`
adds the chip-thickness field.

The alternative was to print instructions and let the card be edited by hand.
That is a bad trade on a deck this size: the two constants go on a specific line,
``constants=`` has to change with them, ``*Depvar`` has to change too, and
getting any of it wrong produces a job that either aborts at preprocessing or --
worse -- runs with SWMODE read out of whatever the 57th number happens to be.

Everything else in the deck is left byte for byte, which is what makes a
three-way comparison meaningful: same mesh, same grains, same seating, same
field, only the criterion differs.
"""

from __future__ import annotations

import os
import re
from typing import Optional

SWMODE_PROP = 57
PSI_PROP = 58
N_PROPS_GRIND2 = 58
DEPVAR_GRIND2 = 22

_USER_MAT = re.compile(r"^(\*User Material,\s*constants\s*=\s*)(\d+)\s*$",
                       re.I)
_DEPVAR = re.compile(r"^\*Depvar(\s*,\s*delete\s*=\s*(\d+))?\s*$", re.I)

SWMODE_TEXT = {
    0: "geometric only: h vs dc. Identical to vumat_grind.for.",
    1: "energy only: W_p L_c vs PSI Kc^2/E. The chip thickness is not read.",
    2: "both: brittle if either criterion says so.",
}


class SwModeError(RuntimeError):
    pass


def set_energy_mode(deck_in: str, deck_out: str, swmode: int,
                    psi: float = 0.0, *,
                    comment: Optional[list] = None) -> dict:
    """Rewrite a 56-constant hybrid card as a 58-constant vumat_grind2 card."""
    if swmode not in (0, 1, 2):
        raise SwModeError("swmode must be 0, 1 or 2, not %r" % (swmode,))
    if not os.path.exists(deck_in):
        raise SwModeError("no such deck: %s" % deck_in)
    with open(deck_in, encoding="ascii") as fh:
        lines = fh.readlines()

    i_mat = None
    n_declared = 0
    for i, ln in enumerate(lines):
        m = _USER_MAT.match(ln.rstrip("\n"))
        if m:
            if i_mat is not None:
                raise SwModeError("more than one *User Material block")
            i_mat, n_declared = i, int(m.group(2))
    if i_mat is None:
        raise SwModeError("no *User Material block: this is not a VUMAT deck")
    if n_declared == N_PROPS_GRIND2:
        raise SwModeError("this deck already carries %d constants; rewriting it "
                          "again would append two more" % N_PROPS_GRIND2)
    if n_declared != 56:
        raise SwModeError(
            "expected the 56-constant card vumat_grind.for reads, found %d. "
            "Refusing to guess where SWMODE and PSI belong." % n_declared)

    # Find the end of the constants block: data lines until the next keyword.
    j = i_mat + 1
    vals: list[str] = []
    while j < len(lines) and not lines[j].startswith("*"):
        vals += [x.strip() for x in lines[j].split(",") if x.strip()]
        j += 1
    n_read = len(vals)
    if n_read != 56:
        raise SwModeError("the card declares 56 constants but %d are written"
                          % n_read)

    # *Depvar must grow too, or SDV21 and SDV22 are silently dropped and the
    # energy ratio has nowhere to live.
    i_dep = None
    delete_sdv = None
    for i, ln in enumerate(lines):
        m = _DEPVAR.match(ln.rstrip("\n"))
        if m:
            i_dep = i
            delete_sdv = m.group(2)
    if i_dep is None:
        raise SwModeError("no *Depvar block")
    if i_dep + 1 >= len(lines) or lines[i_dep + 1].startswith("*"):
        raise SwModeError("*Depvar has no data line")
    old_depvar = int(float(lines[i_dep + 1].split(",")[0]))
    if old_depvar > DEPVAR_GRIND2:
        raise SwModeError("*Depvar is already %d, more than vumat_grind2.for "
                          "writes" % old_depvar)

    # Re-emit the WHOLE constants block packed eight to a line, rather than
    # appending a ninth line holding the two new values.
    #
    # Appending happens to work at 56 constants: 56 = 7*8 exactly, so the
    # appended line of 2 leaves no short line in the middle. At any count that is
    # NOT a multiple of 8 it gives e.g. 8,8,8,8,8,8,8,1,2 -- and Abaqus reads
    # *User Material constants EIGHT to a
    # line, so it takes that 1 as a full line, mis-aligns every value after it
    # and rejects the deck with
    #   ***ERROR: THERE ARE INVALID DATA ASSOCIATED WITH THIS USER DEFINED
    #             MATERIAL DEFINITION
    # which is exactly the failure this project already shipped once. Repacking
    # removes the class of bug rather than this instance of it.
    packed = list(vals) + [repr(float(swmode)), repr(float(psi))]
    if len(packed) != N_PROPS_GRIND2:
        raise SwModeError("built %d constants, vumat_grind2.for reads %d"
                          % (len(packed), N_PROPS_GRIND2))
    block = ["*User Material, constants=%d\n" % N_PROPS_GRIND2]
    block += [", ".join(packed[i:i + 8]) + "\n"
              for i in range(0, len(packed), 8)]

    out = lines[:i_mat] + block + lines[j:]
    shift = len(block) - (j - i_mat)
    i_dep2 = i_dep + (shift if i_dep >= j else 0)
    out[i_dep2 + 1] = "%d,\n" % DEPVAR_GRIND2

    note = ["**",
            "** ---------------- LOCAL ENERGY CRITERION ----------------",
            "** Rewritten for vumat_grind2.for by semgrit_multi.swmode.",
            "**   SWMODE = %d  %s" % (swmode, SWMODE_TEXT[swmode]),
            "**   PSI    = %g%s" % (psi, "  (0 = derive it from dc, giving a"
                                    " threshold of H*dc)" if psi <= 0 else ""),
            "** constants 56 -> %d, *Depvar %d -> %d (SDV21 W_p, SDV22 ratio)."
            % (N_PROPS_GRIND2, old_depvar, DEPVAR_GRIND2),
            "** Submit with:  user=vumat_grind2.for",
            "**"]
    for ln in (comment or []):
        note.insert(-1, "** " + ln)
    ins = i_mat
    out[ins:ins] = [x + "\n" for x in note]

    with open(deck_out, "w", encoding="ascii", newline="\n") as fh:
        fh.writelines(out)
    return {
        "path": deck_out,
        "size_bytes": os.path.getsize(deck_out),
        "swmode": swmode,
        "psi": psi,
        "n_props": N_PROPS_GRIND2,
        "n_depvar": DEPVAR_GRIND2,
        "depvar_was": old_depvar,
        "delete_sdv": int(delete_sdv) if delete_sdv else None,
        "subroutine": "vumat_grind2.for",
    }
