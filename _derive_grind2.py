"""Derive vumat_grind2.for from vumat_grind.for by anchored surgery.

Why a derivation and not a hand-written copy. ``vumat_grind.for`` is a thousand
lines that passed 112 checks, including bit-identity with ``vumat_jh2.for``.
Retyping it to add one criterion would put every one of those lines at risk of a
transcription slip that no reviewer would spot. So the shared body is copied
byte for byte and only the named anchors are touched -- and
``verify_vumat_grind2.py`` then proves the result is still bit-identical to
``vumat_grind.for`` wherever the new criterion is switched off, which is the
check that keeps the two from drifting.

Run it again after any fix to ``vumat_grind.for`` to carry that fix across:

    python _derive_grind2.py            # writes vumat_grind2.for
    python _derive_grind2.py --check    # exits 1 if it is out of date

What it adds
------------
PROPS 57 SWMODE and 58 PSI, SDV 21 and 22, and a second, purely local switch:

    W_p * L_c  >=  PSI * Kc^2 / E     ->  brittle

``W_p`` is accumulated plastic work per unit volume, ``L_c`` is the element's
own characteristic length, which Abaqus hands the VUMAT. No coordinates, no
kinematics, no field variable, no grit count.

It is the same statement Bifano's dc came from, written pointwise instead of
geometrically -- but the exponents do not line up, and pretending they do would
be wrong. The pointwise balance, plastic work per volume ~ H over a depth d
against the fracture energy Kc^2/E, gives

    d_c(energy) = PSI * Kc^2 / (E H) = PSI * (H/E)^1 * (Kc/H)^2

an exponent of +1 on H/E. The two published geometric forms use +0.5 and -1. So
the local criterion is a THIRD member of the same family, not either of them.

Rather than leave that as a discrepancy, PSI is defaulted so the local criterion
trips at exactly the dc the deck already chose:

    PSI = dc * E * H / Kc^2     =>     W_p * L_c >= H * dc

which reads as plainly as it should: brittle once the plastic work per unit area
exceeds what it costs to plastically remove a layer of thickness dc at flow
stress H. Whichever dc form the deck used, the two criteria then agree by
construction, and PSI can still be overridden to calibrate against data.
"""

import io
import os
import sys

SRC = "vumat_grind.for"
DST = "vumat_grind2.for"


def edit(text, anchor, replacement, count=1):
    n = text.count(anchor)
    if n != count:
        raise SystemExit("anchor found %d times, expected %d:\n%r"
                         % (n, count, anchor[:90]))
    return text.replace(anchor, replacement)


HEADER_OLD = """C=======================================================================
C  VUMAT_GRIND.FOR
C
C  Hybrid ductile / brittle constitutive model for SINGLE-GRIT grinding,
C  Abaqus/Explicit VUMAT.
C
C  One material point follows ONE of two laws, chosen once, from the
C  undeformed chip thickness h that the grit takes at that point's
C  station along the scratch:
"""

HEADER_NEW = """C=======================================================================
C  VUMAT_GRIND2.FOR
C
C  Hybrid ductile / brittle constitutive model for grinding with ANY
C  number of abrasives, Abaqus/Explicit VUMAT.
C
C  DERIVED FROM vumat_grind.for BY _derive_grind2.py. Do not edit the
C  shared body here: edit vumat_grind.for and re-run the derivation, or
C  the two will drift. verify_vumat_grind2.py proves this file is still
C  bit-identical to vumat_grind.for wherever SWMODE = 0.
C
C  It adds a SECOND, purely local criterion that needs no geometry at
C  all, so the switch works for one abrasive or seven hundred, for a
C  plunge or a traverse, and for a second pass over the same groove:
C
C        W_p * L_c  >=  PSI * Kc^2 / E     ->  brittle
C
C  W_p is the accumulated plastic work per unit volume and L_c the
C  element's own characteristic length, which Abaqus hands the VUMAT.
C
C  This is the statement Bifano's critical depth came from, written
C  pointwise -- but the exponents do NOT line up and it would be wrong to
C  claim they do. Plastic work per volume of order H over a depth d,
C  against the fracture energy Kc^2/E, gives
C
C        d_c(energy) = PSI Kc^2/(E H) = PSI (H/E)^1 (Kc/H)^2
C
C  an exponent of +1 on H/E, where the two published geometric forms use
C  +0.5 and -1. The local criterion is a THIRD member of the family.
C
C  So PSI is not defaulted to lambda_c. It is defaulted to the value that
C  makes the local criterion trip at exactly the dc the deck already
C  chose:
C
C        PSI = dc E H / Kc^2      giving      W_p L_c >= H dc
C
C  which reads as plainly as it should: brittle once the plastic work per
C  unit area exceeds the cost of plastically removing a layer of
C  thickness dc at flow stress H. Whichever dc form the card carries, the
C  two criteria then agree by construction. Override PSI to calibrate it
C  against scratch or nanoindentation data instead.
C
C  Unlike the geometric switch, this one triggers on HISTORY, so a point
C  starts ductile and turns brittle as the cut deepens under it -- which
C  is the physical transition, rather than being told in advance where it
C  happens. Once triggered it latches.
C
C  SWMODE, PROPS(57):
C     0  geometric only: h vs dc. Identical to vumat_grind.for.
C     1  energy only: every point starts ductile and flips when the work
C        criterion is met. Needs no chip thickness at all.
C     2  both: brittle if either criterion says so.
C
C  On a flip the JH-2 branch inherits the Johnson-Cook damage as its own.
C  The two damages are different mechanisms with the same meaning -- the
C  fraction of the way to failure -- and carrying it is the only
C  continuation that does not either forgive or double-count the damage
C  the point already has.
C
C  ONE POINT OF PHYSICS TO WATCH. The energy criterion is regularised by
C  L_c, so it is mesh-dependent by construction, exactly as every
C  energy-based failure criterion is. Halving the element size halves the
C  work density needed to trigger. That is the correct behaviour for a
C  fracture-energy criterion and it means PSI is calibrated FOR A MESH.
C  State the element size alongside PSI.
C
C-----------------------------------------------------------------------
C  ORIGINAL HEADER OF vumat_grind.for FOLLOWS
C=======================================================================
C  VUMAT_GRIND.FOR
C
C  Hybrid ductile / brittle constitutive model for SINGLE-GRIT grinding,
C  Abaqus/Explicit VUMAT.
C
C  One material point follows ONE of two laws, chosen once, from the
C  undeformed chip thickness h that the grit takes at that point's
C  station along the scratch:
"""

PROPS_DOC_OLD = """C    56  IHMODE  0 coords, 1 field variable 1, 2 all ductile,
C                3 all brittle
"""

PROPS_DOC_NEW = """C    56  IHMODE  0 coords, 1 field variable 1, 2 all ductile,
C                3 all brittle
C
C  57..58  the local energy criterion (vumat_grind2 only):
C    57  SWMODE  0 geometric only (this file then equals vumat_grind.for),
C                1 energy only, 2 both
C    58  PSI     calibration constant in W_p L_c >= PSI Kc^2/E.
C                <=0 -> dc*E*H/Kc^2, which makes the local criterion trip
C                at the same dc the geometric one uses; if H or Kc is
C                absent it falls back to LAMC, PROPS(48).
"""

SDV_DOC_OLD = """C   20  INIT    1 once the point has been initialised
"""

SDV_DOC_NEW = """C   20  INIT    1 once the point has been initialised
C   21  WPLAS   accumulated plastic work per unit volume  [stress]
C   22  ERATIO  W_p L_c E / Kc^2, the energy criterion's own ratio.
C               Reaching PSI is what flips the point. Plot it to see the
C               transition coming.
"""

READ_OLD = """      if (nprops .ge. 56) ihmode = int(props(56) + half)
"""

READ_NEW = """      if (nprops .ge. 56) ihmode = int(props(56) + half)
C
C     The local energy criterion. Off unless asked for, so a card written
C     for vumat_grind.for behaves identically here.
      iswm = 0
      if (nprops .ge. 57) iswm = int(props(57) + half)
      if (iswm .lt. 0 .or. iswm .gt. 2) iswm = 0
      psi = z0
      if (nprops .ge. 58) psi = props(58)
      if (psi .le. z0) then
C       Default: the value that makes the local criterion trip at the same
C       critical depth the geometric one uses, so the two agree instead of
C       being two different thresholds with one name.
        if (dcut .gt. z0 .and. hardn .gt. z0 .and. rkic .gt. z0
     1      .and. ejc .gt. z0) then
          psi = dcut*ejc*hardn/(rkic*rkic)
        else
          psi = rlamc
        endif
      endif
C     Fracture energy per unit area, Kc^2/E, times PSI: the work per area
C     a point has to do before it is allowed to fracture. With the default
C     PSI this is exactly H*dc.
      gcrit = z0
      if (rkic .gt. z0 .and. ejc .gt. z0) gcrit = psi*rkic*rkic/ejc
"""

MODE_OLD = """          imode = 1
          if (hloc .ge. dcut) imode = 2
          if (ihmode .eq. 2) imode = 1
"""

MODE_NEW = """          imode = 1
          if (hloc .ge. dcut) imode = 2
C         Energy only: start every point ductile and let the work
C         criterion decide. The chip thickness is then never consulted,
C         which is the whole point of the mode.
          if (iswm .eq. 1) imode = 1
          if (ihmode .eq. 2) imode = 1
"""

INIT_OLD = """          if (nstatev .ge. 20) stateNew(km,20) = z1
"""

INIT_NEW = """          if (nstatev .ge. 20) stateNew(km,20) = z1
          if (nstatev .ge. 21) stateNew(km,21) = z0
          if (nstatev .ge. 22) stateNew(km,22) = z0
"""

WORK_OLD = """        if (nstatev .ge. 19) stateNew(km,19) = fsge
"""

WORK_NEW = """        if (nstatev .ge. 19) stateNew(km,19) = fsge
C
C       The local energy criterion. Plastic work per unit volume times the
C       element's characteristic length is work per unit area; once it
C       reaches the fracture energy the point is allowed to crack, and
C       from the next increment it follows JH-2.
C
C       Deliberately evaluated AFTER this increment's stress: the point
C       has already done this increment's work under the ductile law, and
C       re-solving it under the brittle one would need the increment
C       repeated. One increment of lag at dt ~ 1e-10 s is nothing.
        wpold = z0
        if (nstatev .ge. 21) wpold = stateOld(km,21)
        if (wpold .lt. z0) wpold = z0
        wpnew = wpold + qpl*dep
        clen = charLength(km)
        if (clen .le. z0) clen = z1
        eratio = z0
        if (gcrit .gt. z0) eratio = wpnew*clen/gcrit
        if (nstatev .ge. 21) stateNew(km,21) = wpnew
        if (nstatev .ge. 22) stateNew(km,22) = eratio
        if (iswm .ne. 0 .and. eratio .ge. z1 .and. nstatev .ge. 13) then
          stateNew(km,13) = z2
        endif
"""

DECL_OLD = """      integer itcut, idcf, ihmode, imode
"""

DECL_NEW = """      integer itcut, idcf, ihmode, imode, iswm
"""


def derive(src_text: str) -> str:
    t = src_text
    t = edit(t, HEADER_OLD, HEADER_NEW)
    t = edit(t, PROPS_DOC_OLD, PROPS_DOC_NEW)
    t = edit(t, SDV_DOC_OLD, SDV_DOC_NEW)
    t = edit(t, DECL_OLD, DECL_NEW)
    t = edit(t, READ_OLD, READ_NEW)
    t = edit(t, MODE_OLD, MODE_NEW)
    t = edit(t, INIT_OLD, INIT_NEW)
    t = edit(t, WORK_OLD, WORK_NEW)
    t = t.replace("C  STATEV  (20; *Depvar, delete=12)",
                  "C  STATEV  (22; *Depvar, delete=12)")
    # The brittle branch is left entirely alone. Overwriting ERATIO there with
    # a magic 1 was tempting and wrong: a point that was ALWAYS brittle should
    # read 0, meaning the energy criterion never had anything to say about it,
    # and a point that FLIPPED should keep the ratio it flipped at, which is
    # more informative than a constant. Both fall out of carrying the state
    # forward, which the shared do-45 copy already does.
    return t


def check_columns(text: str) -> list:
    bad = []
    for i, ln in enumerate(text.split("\n"), start=1):
        if not ln.strip() or ln[0] in "Cc*!":
            continue
        if len(ln.rstrip()) > 72:
            bad.append(i)
    return bad


def main(argv) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    src = io.open(os.path.join(here, SRC), encoding="ascii").read()
    out = derive(src)
    bad = check_columns(out)
    if bad:
        print("derived source has %d lines past column 72: %s"
              % (len(bad), bad[:8]))
        return 1
    dst = os.path.join(here, DST)
    if "--check" in argv:
        if not os.path.exists(dst):
            print("%s does not exist" % DST)
            return 1
        cur = io.open(dst, encoding="ascii").read()
        if cur != out:
            print("%s is OUT OF DATE with respect to %s" % (DST, SRC))
            return 1
        print("%s is up to date with %s" % (DST, SRC))
        return 0
    io.open(dst, "w", encoding="ascii", newline="\n").write(out)
    print("wrote %s  (%d lines, %d added over %s)"
          % (DST, out.count("\n"), out.count("\n") - src.count("\n"), SRC))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
