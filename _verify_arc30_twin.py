"""Verify the arc30 configuration at a size this machine can parse.

The 30 degree deck is 602 MB and 36,410 grains; verify_rigid_deck.py reads the whole
file into a Python string and then builds a dict of node arrays, which needs several GB
and dies with MemoryError here. That is a limit of the checker on this laptop, not a
statement about the deck.

So: build the *same configuration* -- same diameter, same rim, same full-face dressing
at C100, same protrusion statistics, same mesh, same standoff and infeed fraction --
with a shorter arc, and put that through both verifiers. Every code path the 30 degree
deck used is exercised; only the number of grains differs.

Run the verifiers on the full deck on the HPC before submitting. They ship in the
bundle and take a few minutes with enough RAM.
"""
import dataclasses
import os
import pickle
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from semgrit.build_deck import build_deck, plan_deck
from _build_arc30 import PARAMS as FULL

OUT = "_arc30_twin"


def main():
    solids = pickle.load(open("WHEEL_FIXED/1_measurements/grain_library.pkl",
                              "rb"))["solids"]
    sieve = [s for s in solids if 2.0 <= s.height_um <= 6.6]
    os.makedirs(OUT, exist_ok=True)

    # Scale the arc and the sweep together, so the block still stays on the dressed
    # band for the whole pass. Only the number of grains changes; the infeed direction,
    # the tangency, the entry placement, the sets and the contact are size-independent.
    twin = dataclasses.replace(
        FULL, name="arc30_twin", sector_deg=8.0, travel_mm=0.25,
        write_step=False, write_stl=False, also_write_cae_deck=True)
    # same infeed rule as the full deck: 95% of the face-to-bond gap
    plan = plan_deck(twin, sieve)
    twin.analysis.depth_of_cut_um = round(0.95 * plan["depth_ceiling_um"], 4)
    print("twin: %.1f mm dia, %.2f deg, arc %.3f mm, %d grains, ae %.3f um"
          % (twin.diameter_mm, twin.sector_deg, plan["arc_length_mm"],
             plan["n_grits"], twin.analysis.depth_of_cut_um))
    print("full: %.1f mm dia, %.2f deg, arc 4.189 mm, 36410 grains, travel %.2f mm"
          % (FULL.diameter_mm, FULL.sector_deg, FULL.travel_mm))
    print("both place the block with: %s" % FULL.wp_position)
    same = [f.name for f in dataclasses.fields(twin)
            if f.name not in ("name", "sector_deg", "write_step", "write_stl",
                              "analysis", "also_write_cae_deck")
            and f.name != "travel_mm"
            and getattr(twin, f.name) != getattr(FULL, f.name)]
    print("settings that differ beyond the arc length and CAD flags: %s"
          % (same or "none"))

    info = build_deck(twin, sieve, OUT)
    print("built %s  %.1f MB" % (os.path.basename(info["path"]),
                                 info["size_bytes"] / 1e6))
    ok = True
    for deck in [info["path"], info.get("cae_deck")]:
        if not deck:
            continue
        for v in ("verify_rigid_deck.py", "verify_rigid_deck2.py"):
            r = subprocess.run([sys.executable, v, deck], capture_output=True,
                               text=True)
            fails = [l.strip() for l in r.stdout.splitlines() if "[FAIL]" in l]
            tot = [l.strip() for l in r.stdout.splitlines() if "TOTAL" in l]
            print("  %-24s %-22s %s   %s"
                  % (os.path.basename(deck), v,
                     "PASS" if r.returncode == 0 else "FAIL",
                     tot[-1] if tot else ""))
            for l in fails:
                print("      " + l)
            ok = ok and r.returncode == 0
    print()
    print("TWIN VERIFIES" if ok else "TWIN FAILED - the configuration is wrong")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
