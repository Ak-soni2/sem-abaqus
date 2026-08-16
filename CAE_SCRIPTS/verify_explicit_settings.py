"""Confirm the model really is Explicit with element deletion off.

Run from Abaqus/CAE:   File -> Run Script...

Rather than trusting the Element Type dialog part by part, this writes the model out
as an input file and reads back what Abaqus actually recorded. The input deck is the
thing the solver consumes, so it is the only definitive check.

Reports, per model:
  * every *Element type present, and how many elements use it
  * any *Section Controls block and whether ELEMENT DELETION is on
  * the step procedure
"""

import os
import re

from abaqus import mdb

MODEL = None          # None = every model that has parts
OUT_DIR = None        # None = Abaqus's current working directory


def check(model_name):
    m = mdb.models[model_name]
    if not m.parts:
        return
    print("")
    print("=" * 74)
    print("MODEL %r" % model_name)
    print("=" * 74)

    job = "_typecheck_%s" % re.sub(r"[^A-Za-z0-9_]", "_", model_name)[:30]
    if job in mdb.jobs:
        del mdb.jobs[job]
    mdb.Job(name=job, model=model_name)
    print("writing %s.inp ..." % job)
    mdb.jobs[job].writeInput(consistencyChecking=OFF)

    path = job + ".inp"
    if OUT_DIR:
        path = os.path.join(OUT_DIR, path)
    if not os.path.isfile(path):
        print("could not find %s -- check Abaqus's working directory" % path)
        return

    el_counts = {}
    controls = []
    steps = []
    cur = None
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("**"):
                continue
            if s.startswith("*"):
                low = s.lower()
                if low.startswith("*element,"):
                    mt = re.search(r"type\s*=\s*([A-Za-z0-9]+)", s, re.I)
                    cur = mt.group(1).upper() if mt else "?"
                    el_counts.setdefault(cur, 0)
                elif low.startswith("*section controls"):
                    controls.append(s)
                    cur = None
                elif low.startswith("*step"):
                    steps.append(s)
                    cur = None
                elif low.startswith("*dynamic") or low.startswith("*static"):
                    steps.append("   " + s)
                    cur = None
                else:
                    cur = None
            elif cur and s:
                el_counts[cur] += 1

    print("")
    print("  element types written to the deck:")
    rigid = ("R3D3", "R3D4", "RAX2", "R2D2")
    problems = []
    for t in sorted(el_counts):
        kind = "rigid" if t in rigid else "deformable"
        # Explicit has no full-integration C3D8; its presence means Standard.
        flag = ""
        if t in ("C3D8", "C3D6", "C3D10", "S4", "S3", "CPE4", "CPS4"):
            flag = "  <-- full integration, this is a STANDARD element"
            problems.append(t)
        print("    %-10s %10d elements   %-11s%s" % (t, el_counts[t], kind, flag))

    print("")
    if controls:
        print("  *Section Controls found:")
        for c in controls:
            on = "element deletion=yes" in c.lower()
            print("    %s   %s" % (c, "<-- DELETION IS ON" if on else "(deletion off)"))
            if on:
                problems.append("section controls")
    else:
        print("  no *Section Controls block -> element deletion is not enabled")

    print("")
    print("  step definition:")
    for s in steps:
        print("    %s" % s)
    if not any("*dynamic" in s.lower() and "explicit" in s.lower() for s in steps):
        print("    <-- no '*Dynamic, Explicit' found")
        problems.append("step")

    print("")
    if problems:
        print("  RESULT: not clean -- see the flags above (%s)" % ", ".join(
            sorted(set(str(p) for p in problems))))
    else:
        print("  RESULT: clean. Explicit elements, no element deletion, Explicit step.")


names = [MODEL] if MODEL else sorted(mdb.models.keys())
for n in names:
    if n in mdb.models:
        check(n)
print("")
print("done")
