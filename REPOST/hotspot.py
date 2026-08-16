"""Find the peak-stress element and say WHY it is at that stress.

    abaqus python hotspot.py <job>.odb
    abaqus python hotspot.py --table          (after the .odb runs, pure json)

WHAT THIS IS FOR. single_abrasive2_sic reported a peak Mises of 40 GPa, which
is 2.8x SiC's HEL. Two mechanisms in vumat_grind.for can produce that and they
are not distinguishable from a contour plot:

  BRITTLE, and legitimate.  JH-2's intact surface sig = A(P*+T*)^N * SIGHEL is
      NOT capped -- only the fractured branch gets SFMAX (vumat_grind2.for
      line 1009). Reaching 40 GPa needs P = 35,736 MPa, i.e. 6.06 x PHEL and
      2.47 x HEL. High, but a shock under a grit is where you would find it.

  DUCTILE, and a size-effect artefact.  The strain-gradient term uses
      eta = 4 ep / h with h floored at one Burgers vector, and for SiC
      sgec = r' b (M alpha G)^2 = 18,772 against 5.84 for sandstone -- 3,214x.
      40 GPa needs only h ~ 2.5 nm at ep = 0.05. The single-abrasive deck is
      the only one with PROPS(56) = 0, where hloc = h0 + hg*u - u^2/2r clamps
      to zero across the rubbing zone, so nanometre h is freely available.
      At the Burgers floor the same term reaches 110 GPa at ep = 0.05.

SDV19 is FSGE, the amplification SEFF/SJC, and it settles it outright:
FSGE ~ 1 means the size effect is idle and the stress is JH-2's; FSGE ~ 4.8
means the size effect built the whole thing.

The postprocessor reads SDVs at the LAST frame only, by which time the hot
element is usually deleted. This walks every frame instead, which is the only
reason it exists as a separate script.
"""

from __future__ import print_function

import glob
import json
import os
import sys

# Only probe() needs Abaqus. --table is plain json, so it stays runnable with
# a normal interpreter after you copy the results back off the solve machine.
try:
    from odbAccess import openOdb
    from abaqusConstants import MISES
except ImportError:
    openOdb = MISES = None

# SDVs worth printing at the hot spot, from the vumat_grind2.for header.
SDV = [("SDV1", "D        damage 0..1"),
       ("SDV2", "EPBAR    equiv plastic strain"),
       ("SDV3", "P        pressure (+ve compression)  MPa"),
       ("SDV4", "Q        von Mises                   MPa"),
       ("SDV13", "MODE     1 = ductile JC+SGE, 2 = brittle JH-2"),
       ("SDV14", "HLOC     chip thickness this point got   mm"),
       ("SDV15", "DCLOC    dc actually used                mm"),
       ("SDV17", "SJC      Johnson-Cook flow stress   MPa"),
       ("SDV18", "SEFF     after SGE, before damage   MPa"),
       ("SDV19", "FSGE     SEFF/SJC  <-- THE ANSWER"),
       ("SDV22", "ERATIO   energy criterion ratio")]


def hot_frame(step):
    """(max mises, element label, frame index) over every frame in the step."""
    best = (-1.0, None, None)
    for i, fr in enumerate(step.frames):
        if "S" not in fr.fieldOutputs:
            continue
        for v in fr.fieldOutputs["S"].getScalarField(invariant=MISES).values:
            if v.data > best[0]:
                best = (float(v.data), v.elementLabel, i)
    return best


def probe(odb_path):
    if openOdb is None:
        print("run this with Abaqus' Python:  abaqus python hotspot.py <job>.odb")
        raise SystemExit(2)
    odb = openOdb(odb_path, readOnly=True)
    try:
        step = odb.steps[list(odb.steps.keys())[-1]]
        smax, elem, iframe = hot_frame(step)
        if elem is None:
            print("no S field in %s" % odb_path)
            return None
        fr = step.frames[iframe]
        print("=" * 70)
        print("%s" % os.path.basename(odb_path))
        print("  peak Mises %.4g MPa at element %d, frame %d/%d (t = %.4g s)"
              % (smax, elem, iframe, len(step.frames) - 1, fr.frameValue))
        got = {}
        for key, desc in SDV:
            if key not in fr.fieldOutputs:
                continue
            for v in fr.fieldOutputs[key].values:
                if v.elementLabel == elem:
                    got[key] = float(v.data)
                    print("    %-6s %-44s %14.6g" % (key, desc, got[key]))
                    break
        print("  " + "-" * 66)
        mode = int(round(got.get("SDV13", 0.0)))
        fsge = got.get("SDV19")
        if mode == 2:
            print("  VERDICT: BRITTLE branch. The stress is JH-2's intact")
            print("           surface, uncapped. Cross-check SDV3 against the")
            print("           35,736 MPa that 40 GPa implies for SiC.")
        elif mode == 1 and fsge is not None:
            print("  VERDICT: DUCTILE branch, size effect amplifying x%.3g." % fsge)
            if fsge > 2.0:
                print("           This stress is BUILT BY THE SGE TERM, not by")
                print("           JH-2. h at this point was %.4g mm."
                      % got.get("SDV14", float("nan")))
                print("           Do not quote the peak stress as a material")
                print("           strength -- it is a gradient-hardening value.")
            else:
                print("           Size effect is near idle; the stress is the")
                print("           Johnson-Cook flow stress, not an artefact.")
        else:
            print("  VERDICT: inconclusive -- SDV13/SDV19 absent. The deck")
            print("           requests SDV output, so check the VUMAT ran.")
        out = {"odb": os.path.basename(odb_path), "peak_mises_MPa": smax,
               "element": elem, "frame": iframe, "time_s": fr.frameValue,
               "sdv": got}
        base = os.path.splitext(odb_path)[0]
        with open(base + "_hotspot.json", "w") as fh:
            json.dump(out, fh, indent=2, sort_keys=True)
        print("  wrote %s_hotspot.json" % os.path.basename(base))
        return out
    finally:
        odb.close()


def table():
    """One comparison across every *_summary.json in the folder."""
    rows = []
    for p in sorted(glob.glob("*_summary.json")):
        d = json.load(open(p))
        f, e, r = d.get("forces", {}), d.get("energy", {}), d.get("removal", {})
        s = d.get("sdv", {})
        hs = p.replace("_summary.json", "_hotspot.json")
        h = json.load(open(hs)) if os.path.exists(hs) else {}
        ie = e.get("ALLIE_final") or 0.0
        rows.append((
            d.get("odb", p)[:26],
            f.get("peak_magnitude_N"), ie,
            e.get("artificial_fraction"),
            r.get("elements_deleted"),
            s.get("n_ductile_alive"), s.get("n_brittle_alive"),
            d.get("specific_energy_J_mm3"),
            h.get("peak_mises_MPa"),
            (h.get("sdv") or {}).get("SDV19")))
    if not rows:
        print("no *_summary.json here -- run postprocess_odb.py first")
        return
    hdr = ("job", "peakF_N", "ALLIE_mJ", "AE/IE", "deleted", "duct", "britt",
           "J/mm3", "maxMises", "FSGE")
    print("%-26s %10s %11s %7s %8s %6s %6s %9s %10s %7s" % hdr)
    for r in rows:
        cells = [r[0]]
        for v in r[1:]:
            cells.append("-" if v is None else
                         ("%.4g" % v if isinstance(v, float) else str(v)))
        print("%-26s %10s %11s %7s %8s %6s %6s %9s %10s %7s" % tuple(cells))


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--help"):
        print(__doc__)
    elif args[0] == "--table":
        table()
    else:
        for a in args:
            probe(a)
