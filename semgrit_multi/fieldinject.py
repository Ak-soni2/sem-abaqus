"""Write a per-node chip-thickness field into an already-built deck.

``vumat_grind.for`` reads field variable 1 as the chip thickness when
PROPS(56) = 1, and that path is covered by ``verify_hybrid_deck.py``. So the
whole multi-abrasive extension needs no change to the subroutine and no change
to ``semgrit``: build the deck exactly as before with ``h_source = 1``, then add
one ``*Initial Conditions, type=FIELD`` block to it.

Injecting rather than teaching the writer a new keyword is deliberate. The
writer produces two decks that Abaqus has already accepted and that
``_check_presets.py`` compares against byte for byte; adding an output path to
it would put those at risk for no benefit. Injection is additive, reversible,
and works on any deck the existing pipeline emits.

Where the block goes
--------------------
``*Initial Conditions`` is a model-level keyword: it must come after
``*End Assembly`` and before the first ``*Step``. Both anchors are unambiguous
in these decks, and this module refuses to write anything if it cannot find
exactly one of each.
"""

from __future__ import annotations

import os
import re
from typing import Optional, Sequence

import numpy as np


class InjectError(RuntimeError):
    pass


_STEP = re.compile(r"^\*Step\b", re.I)
_END_ASSEMBLY = re.compile(r"^\*End Assembly\b", re.I)
_USER_MATERIAL = re.compile(r"^\*User Material,\s*constants\s*=\s*(\d+)", re.I)

# PROPS index (1-based) of the h source in vumat_grind.for, and the value that
# means "read field variable 1".
H_SOURCE_PROP = 56
H_SOURCE_FIELD = 1


def _read_props(lines: Sequence[str]) -> tuple:
    """(first_value_line, n_lines, values) of the *User Material block."""
    for i, ln in enumerate(lines):
        m = _USER_MATERIAL.match(ln)
        if not m:
            continue
        n = int(m.group(1))
        vals: list[float] = []
        j = i + 1
        while j < len(lines) and len(vals) < n:
            if lines[j].startswith("*"):
                break
            vals += [float(x) for x in lines[j].split(",") if x.strip()]
            j += 1
        if len(vals) != n:
            raise InjectError("*User Material declares %d constants but %d "
                              "were readable" % (n, len(vals)))
        return i + 1, j - (i + 1), vals
    raise InjectError("no *User Material block: this is not a VUMAT deck")


def inject_field(deck_in: str, deck_out: str, values: np.ndarray,
                 *, instance: str = "WP-1", variable: int = 1,
                 per_line: int = 1, comment: Optional[Sequence[str]] = None,
                 require_h_source: bool = True) -> dict:
    """Copy ``deck_in`` to ``deck_out`` with a nodal field variable added.

    ``values[i]`` is the field value at node ``i + 1`` of ``instance``, which is
    the order :func:`semgrit.wheel_workpiece.build_block_mesh` numbers them and
    the order the deck writes them.
    """
    if not os.path.exists(deck_in):
        raise InjectError("no such deck: %s" % deck_in)
    v = np.asarray(values, dtype=np.float64).ravel()
    if v.size == 0:
        raise InjectError("no field values given")
    if not np.isfinite(v).all():
        raise InjectError("the field carries %d non-finite values"
                          % int((~np.isfinite(v)).sum()))
    if (v < 0).any():
        raise InjectError("the chip thickness cannot be negative; %d values are"
                          % int((v < 0).sum()))

    with open(deck_in, encoding="ascii") as fh:
        lines = fh.readlines()

    steps = [i for i, ln in enumerate(lines) if _STEP.match(ln)]
    ends = [i for i, ln in enumerate(lines) if _END_ASSEMBLY.match(ln)]
    if len(ends) != 1:
        raise InjectError("expected exactly one *End Assembly, found %d"
                          % len(ends))
    if len(steps) != 1:
        raise InjectError(
            "expected exactly one *Step, found %d. A field written before the "
            "wrong step would be applied at the wrong time." % len(steps))
    at = steps[0]
    if at <= ends[0]:
        raise InjectError("*Step precedes *End Assembly; the deck is malformed")

    # The card has to be asking for a field, or the values are inert and the
    # run silently uses whatever PROPS(53..55) happen to say.
    _first, _n, props = _read_props(lines)
    if require_h_source:
        if len(props) < H_SOURCE_PROP:
            raise InjectError(
                "the card has only %d constants, so it is not vumat_grind.for "
                "and has no h source to check" % len(props))
        got = int(round(props[H_SOURCE_PROP - 1]))
        if got != H_SOURCE_FIELD:
            raise InjectError(
                "PROPS(%d) is %d, not %d: this deck's material card does not "
                "read the chip thickness from field variable 1, so injecting "
                "one would do nothing. Build it with HybridParams(h_source=1)."
                % (H_SOURCE_PROP, got, H_SOURCE_FIELD))

    # Count the instance's nodes, so a field of the wrong length is caught here
    # rather than by Abaqus halfway through preprocessing.
    n_nodes = _count_instance_nodes(lines, instance)
    if n_nodes and v.size != n_nodes:
        raise InjectError(
            "the field has %d values but %s has %d nodes"
            % (v.size, instance, n_nodes))

    block: list[str] = ["**\n"]
    block.append("** ---------------- CHIP THICKNESS FIELD ----------------\n")
    block.append("** Undeformed chip thickness at every workpiece node, mm,\n")
    block.append("** swept from the grit trajectories before the run. Read by\n")
    block.append("** vumat_grind.for as field variable 1 because PROPS(%d) = %d.\n"
                 % (H_SOURCE_PROP, H_SOURCE_FIELD))
    for ln in (comment or []):
        block.append("** %s\n" % ln)
    block.append("**\n")
    block.append("*Initial Conditions, type=FIELD, variable=%d\n" % variable)
    if per_line <= 1:
        block += ["%s.%d, %r\n" % (instance, i + 1, float(x))
                  for i, x in enumerate(v)]
    else:
        for i in range(0, v.size, per_line):
            chunk = v[i:i + per_line]
            block.append(", ".join("%s.%d, %r" % (instance, i + j + 1,
                                                  float(x))
                                   for j, x in enumerate(chunk)) + "\n")

    out = lines[:at] + block + lines[at:]
    with open(deck_out, "w", encoding="ascii", newline="\n") as fh:
        fh.writelines(out)

    return {
        "path": deck_out,
        "size_bytes": os.path.getsize(deck_out),
        "n_values": int(v.size),
        "instance": instance,
        "variable": variable,
        "inserted_before_line": at,
        "n_lines_added": len(block),
        "h_min_mm": float(v.min()),
        "h_max_mm": float(v.max()),
        "h_mean_mm": float(v.mean()),
    }


def _count_instance_nodes(lines: Sequence[str], instance: str) -> int:
    """Nodes in the part that ``instance`` instantiates, or 0 if unfindable."""
    want = None
    for ln in lines:
        if ln.lower().startswith("*instance"):
            m = re.search(r"name\s*=\s*([^,]+),\s*part\s*=\s*([^,\s]+)", ln,
                          re.I)
            if m and m.group(1).strip().upper() == instance.upper():
                want = m.group(2).strip().upper()
                break
    if want is None:
        return 0
    part = None
    mode = None
    n = 0
    for ln in lines:
        if ln.startswith("*"):
            key = ln.split(",")[0].strip().lower()
            if key == "*part":
                m = re.search(r"name\s*=\s*([^,\s]+)", ln, re.I)
                part = m.group(1).strip().upper() if m else None
                mode = None
            elif key == "*end part":
                part, mode = None, None
            elif key == "*node":
                mode = "node"
            else:
                mode = None
            continue
        if part == want and mode == "node" and ln.strip():
            n += 1
    return n


def read_field(deck: str, variable: int = 1) -> dict:
    """Read an injected field back out of a deck, for verification."""
    out: dict[int, float] = {}
    grabbing = False
    pat = re.compile(r"^\*Initial Conditions,\s*type\s*=\s*FIELD"
                     r"(?:,\s*variable\s*=\s*(\d+))?", re.I)
    with open(deck, encoding="ascii") as fh:
        for ln in fh:
            if ln.startswith("*"):
                m = pat.match(ln)
                grabbing = bool(m) and int(m.group(1) or 1) == variable
                continue
            if not grabbing or not ln.strip() or ln.startswith("**"):
                continue
            f = [x.strip() for x in ln.split(",") if x.strip()]
            for i in range(0, len(f) - 1, 2):
                node = f[i].split(".")[-1]
                out[int(node)] = float(f[i + 1])
    return out
