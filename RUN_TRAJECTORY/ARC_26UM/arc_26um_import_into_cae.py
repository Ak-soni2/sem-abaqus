"""Load a semgrit wheel deck into Abaqus/CAE as a complete model.

Run from Abaqus/CAE:   File -> Run Script...   and pick this file.

Why a script: CAE's File -> Import -> Part reads only the *Part blocks and skips
the *Assembly, so the grains arrive unplaced; and File -> Import -> Model does not
accept .inp. ModelFromInputFile reads both.

The .inp is located relative to this script, so the pair can be copied to any
folder or machine without editing paths.
"""
import os

from abaqus import mdb, session
from abaqusConstants import *   # noqa: F401,F403  (CAE scripts rely on these names)

INP_NAME = "arc_26um.inp"
MODEL = "arc_26um"

# Extra places to look, in case the script and the deck get separated.
SEARCH_ROOTS = [r"D:\temp", r"C:\temp"]


def find_inp():
    cands = []
    try:
        cands.append(os.path.dirname(os.path.abspath(__file__)))
    except NameError:          # some CAE versions do not set __file__
        pass
    cands.append(os.getcwd())
    for root in SEARCH_ROOTS:
        cands.append(root)
        if os.path.isdir(root):
            try:
                for entry in os.listdir(root):
                    sub = os.path.join(root, entry)
                    if os.path.isdir(sub):
                        cands.append(sub)
                        for entry2 in os.listdir(sub):
                            sub2 = os.path.join(sub, entry2)
                            if os.path.isdir(sub2):
                                cands.append(sub2)
            except OSError:
                pass
    seen = set()
    for c in cands:
        if not c or c in seen:
            continue
        seen.add(c)
        p = os.path.join(c, INP_NAME)
        if os.path.isfile(p):
            return p
    raise IOError(
        "could not find %s. Looked in:\n  %s\nPut this script next to the .inp, "
        "or edit SEARCH_ROOTS at the top." % (INP_NAME, "\n  ".join(sorted(seen)))
    )


INP = find_inp()

if MODEL in mdb.models:
    del mdb.models[MODEL]

print("reading " + INP)
mdb.ModelFromInputFile(name=MODEL, inputFileName=INP)

m = mdb.models[MODEL]
a = m.rootAssembly
print("parts      : %d" % len(m.parts))
print("instances  : %d" % len(a.instances))
grits = len([k for k in a.instances.keys() if k.startswith("G-")])
print("grit instances placed on the wheel : %d" % grits)

# Show the assembled model rather than a single part.
try:
    vp = session.viewports[session.currentViewportName]
    vp.setValues(displayedObject=a)
    vp.view.fitView()
except Exception as exc:              # viewport naming differs between versions
    print("could not set the viewport automatically: %s" % exc)
    print("switch Module to Assembly by hand instead")

print("")
print("If 'grit instances placed' matches the grain count in the report JSON,")
print("every grit is in the model. Zoom onto the outer diameter to see them:")
print("they are only a few microns across.")
