"""Read the .odb with Abaqus' own Python and write the numbers out.

    abaqus python postprocess_odb.py

Runs INSIDE Abaqus, which has no matplotlib -- so this writes CSV and JSON
only, and plotting is left to the host Python. That is not a limitation being
worked around; it is the reason nothing in this project plots from inside
Abaqus.

The result is SDV13, the branch each material point took: 1 ductile, 2
brittle. What matters is not the final frame but HOW IT EVOLVES -- the pass at
which points start flipping. So the brittle fraction is written for every
frame of every step.
"""
import json
import os

from odbAccess import openOdb

JOB = "sagv2_6um"
DUCTILE, BRITTLE = 1, 2


def main():
    odb = openOdb(JOB + ".odb", readOnly=True)
    inst = odb.rootAssembly.instances["WORK-1"]
    n_el = len(inst.elements)

    rows = [("step", "frame", "time", "n_active", "n_ductile", "n_brittle",
             "brittle_fraction", "peak_mises", "peak_peeq")]
    per_step = {}
    for sname in odb.steps.keys():
        step = odb.steps[sname]
        for f, frame in enumerate(step.frames):
            fo = frame.fieldOutputs
            if "SDV13" not in fo:
                continue
            sdv = fo["SDV13"].getSubset(region=inst).values
            status = (fo["SDV12"].getSubset(region=inst).values
                      if "SDV12" in fo else None)
            nd = nb = na = 0
            for i, v in enumerate(sdv):
                if status is not None and status[i].data < 0.5:
                    continue
                na += 1
                if v.data >= 1.5:
                    nb += 1
                elif v.data >= 0.5:
                    nd += 1
            mis = 0.0
            if "S" in fo:
                for v in fo["S"].getSubset(region=inst).values:
                    if v.mises > mis:
                        mis = v.mises
            pk = 0.0
            if "PEEQ" in fo:
                for v in fo["PEEQ"].getSubset(region=inst).values:
                    if v.data > pk:
                        pk = v.data
            frac = (float(nb) / na) if na else 0.0
            rows.append((sname, f, frame.frameValue, na, nd, nb, frac,
                         mis, pk))
            per_step[sname] = frac
    odb.close()

    with open(JOB + "_sdv13.csv", "w") as fh:
        for r in rows:
            fh.write(",".join(str(x) for x in r) + "\n")

    passes = [k for k in per_step if k.startswith("PASS")]
    passes.sort(key=lambda s: int(s[4:]))
    first_brittle = None
    for k in passes:
        if per_step[k] > 0.0:
            first_brittle = int(k[4:])
            break
    out = dict(job=JOB, elements=n_el,
               brittle_fraction_by_step=per_step,
               passes=len(passes),
               first_pass_with_brittle=first_brittle,
               final_brittle_fraction=(per_step[passes[-1]] if passes
                                       else None))
    with open(JOB + "_summary.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print("")
    print("  %s: %d elements" % (JOB, n_el))
    print("  brittle fraction by pass:")
    for k in passes:
        print("    %-8s %.4f" % (k, per_step[k]))
    if first_brittle is None:
        print("")
        print("  NO BRITTLE ELEMENTS in any pass.")
        print("  For the 6 um pad that AGREES with the paper (pure ductile).")
        print("  For the 30 um pad it does NOT: the paper sees fracture there,")
        print("  so the energy criterion is accumulating too slowly.")
    else:
        print("")
        print("  first brittle at PASS%d of %d" % (first_brittle,
                                                    len(passes)))
    print("")
    print("  wrote %s_sdv13.csv and %s_summary.json" % (JOB, JOB))


if __name__ == "__main__":
    main()
