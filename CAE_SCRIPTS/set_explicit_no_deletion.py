"""Force every element in the model to the Explicit library with element deletion OFF.

Run from Abaqus/CAE:   File -> Run Script...   and pick this file.

Doing it by hand means opening the Element Type dialog once per part. With ~100 grain
parts that is not practical, and it is easy to miss one -- a single part left on the
Standard library will stop an Explicit job.

What it changes, for every part in every model (or just MODELS_TO_DO if you set it):

  * element library  -> EXPLICIT
  * element deletion -> OFF
  * full-integration hexes -> reduced integration (C3D8 -> C3D8R), because
    Abaqus/Explicit has no full-integration C3D8

Rigid elements (R3D3, R3D4) also get the Explicit library. It makes no difference to
the input deck -- `*Element, type=R3D3` is written identically either way, because
rigid elements are not library-specific and Abaqus/Explicit accepts them as they are.
It is done purely so the Element Type dialog reads "Explicit" on every part and there
is nothing left to doubt. Element deletion is not applied to them: a rigid body has no
stress and cannot fail.

Steps are only *reported*, never modified. Changing a step's procedure means deleting
and recreating it, which silently discards every load, BC and interaction attached to
it -- far too destructive to do without asking.
"""

import abaqusConstants as AC
import mesh
from abaqus import mdb

# Leave empty to process every model in the database.
MODELS_TO_DO = []

# Deformable element codes, mapped to their Abaqus/Explicit equivalent.
# Explicit has no full-integration C3D8 or C3D6, so those drop to reduced integration.
EXPLICIT_EQUIVALENT = {
    "C3D8": "C3D8R", "C3D8R": "C3D8R", "C3D8I": "C3D8R", "C3D8H": "C3D8R",
    "C3D6": "C3D6", "C3D4": "C3D4", "C3D4H": "C3D4",
    "C3D10": "C3D10M", "C3D10M": "C3D10M", "C3D10H": "C3D10M",
    "S3": "S3R", "S3R": "S3R", "S4": "S4R", "S4R": "S4R",
    "CPE4": "CPE4R", "CPE4R": "CPE4R", "CPS4": "CPS4R", "CPS4R": "CPS4R",
}
RIGID_CODES = ("R3D3", "R3D4", "RAX2", "R2D2")


def element_types_in(part):
    """Distinct element type strings present in a part."""
    found = set()
    try:
        for e in part.elements:
            found.add(str(e.type))
    except Exception as exc:
        print("      could not read elements: %s" % exc)
    return found


def apply_to_part(part):
    """Returns (changed_codes, rigid_codes, errors)."""
    present = element_types_in(part)
    if not present:
        return [], [], []

    rigid = sorted(t for t in present if t in RIGID_CODES)
    deformable = sorted(t for t in present if t not in RIGID_CODES)

    targets = []          # (from, to, code, is_rigid)
    unknown = []
    for t in deformable:
        tgt = EXPLICIT_EQUIVALENT.get(t)
        code = getattr(AC, tgt, None) if tgt else None
        if code is None:
            unknown.append(t)
            continue
        targets.append((t, tgt, code, False))
    for t in rigid:
        code = getattr(AC, t, None)
        if code is None:
            unknown.append(t)
            continue
        targets.append((t, t, code, True))

    if not targets:
        return [], rigid, unknown

    elem_types = []
    for _, _, code, is_rigid in targets:
        if is_rigid:
            # No elemDeletion on rigid elements: a rigid body has no stress and
            # cannot fail, and passing it is rejected.
            elem_types.append(mesh.ElemType(elemCode=code, elemLibrary=AC.EXPLICIT))
            continue
        try:
            elem_types.append(mesh.ElemType(
                elemCode=code, elemLibrary=AC.EXPLICIT, elemDeletion=AC.OFF))
        except Exception:
            # Older releases reject elemDeletion on some codes; set the library at least.
            elem_types.append(mesh.ElemType(elemCode=code, elemLibrary=AC.EXPLICIT))

    errors = []
    try:
        part.setElementType(regions=(part.elements,), elemTypes=tuple(elem_types))
    except Exception as exc:
        # Geometry-based parts want cells rather than an element sequence.
        try:
            part.setElementType(regions=(part.cells,), elemTypes=tuple(elem_types))
        except Exception as exc2:
            errors.append("%s / %s" % (exc, exc2))

    labels = []
    for a, b, _, is_rigid in targets:
        if is_rigid:
            labels.append("%s (rigid, library only)" % a)
        else:
            labels.append("%s -> %s" % (a, b))
    return labels, rigid, errors + unknown


def main():
    names = MODELS_TO_DO if MODELS_TO_DO else sorted(mdb.models.keys())
    print("=" * 74)
    print("Setting every element to EXPLICIT with element deletion OFF")
    print("=" * 74)

    for mname in names:
        if mname not in mdb.models:
            print("model %r not found, skipping" % mname)
            continue
        m = mdb.models[mname]
        print("")
        print("MODEL %r  (%d parts)" % (mname, len(m.parts)))

        n_done = n_untouched = 0
        conversions = {}
        problems = []
        for pname in sorted(m.parts.keys()):
            changed, rigid, errs = apply_to_part(m.parts[pname])
            if errs:
                problems.append((pname, errs))
            if changed:
                n_done += 1
                for c in changed:
                    conversions[c] = conversions.get(c, 0) + 1
            else:
                n_untouched += 1

        print("  parts set to the Explicit library : %d of %d" % (n_done, len(m.parts)))
        if n_untouched:
            print("  parts with nothing to set        : %d" % n_untouched)
        for k in sorted(conversions):
            print("    %-18s on %d part(s)" % (k, conversions[k]))
        if problems:
            print("  parts that could not be set      : %d" % len(problems))
            for pname, errs in problems[:6]:
                print("    %s : %s" % (pname, errs[0]))

        # Steps are reported, never touched.
        analysis_steps = [s for s in m.steps.keys() if s != "Initial"]
        if analysis_steps:
            print("  steps:")
            wrong = []
            for sname in analysis_steps:
                proc = type(m.steps[sname]).__name__
                ok = "ExplicitDynamicsStep" in proc
                if not ok:
                    wrong.append(sname)
                print("    %-22s %-28s %s" % (
                    sname, proc, "OK" if ok else "<-- NOT Explicit"))
            # Only warn when there is actually something wrong.
            if wrong:
                print("")
                print("  %s is not Explicit. A step's procedure cannot be changed in"
                      % ", ".join(wrong))
                print("  place: doing so means deleting and recreating it, which discards")
                print("  the loads, BCs and interactions attached to it. Delete and")
                print("  recreate it as Dynamic, Explicit yourself, then reapply them.")
            else:
                print("    all analysis steps are already Explicit - nothing to do")

    print("")
    print("=" * 74)
    print("Done. Check the counts above, then re-open Element Type on one part to")
    print("confirm it reads Explicit with element deletion = No.")
    print("=" * 74)


main()
