"""A 30 degree grinding wheel arc, whole face dressed at real concentration.

Three things were asked for, and one of them fights the other two, so the reasoning is
recorded here.

1. **30 degrees, clearly an arc.** For a sector, sagitta/arc = theta/8, so a 30 degree
   segment is 6.5% deep relative to its chord *at any radius* -- it always reads as an
   arc. What the radius changes is the arc's *length*, and therefore the grain count.
2. **Shorter outer radius.** Taken as the lever for keeping a 30 degree arc affordable:
   D = 16 mm gives a 4.19 mm arc instead of the 13.1 mm a 50 mm wheel would give, and a
   1.2 mm rim on an 8 mm radius reads as a fat annular segment rather than a gentle
   banana.
3. **The workpiece as close as possible.** It is already as close as geometry allows:
   the tallest grain under the block is *exactly* tangent to it, to within the write
   guard of 0.001 nm, and both verifiers check that. It cannot be closer without
   initial overclosure, which Abaqus/Explicit answers with a spurious force spike -- the
   depth of cut is the right way in, and here it is set to 95% of the face-to-bond gap
   rather than the usual 85%.

   The apparent gap is not standoff, it is **curvature**: a flat block sits b^2/2R above
   a round wheel at distance b from the tangent point. Shrinking the radius makes that
   *worse*, not better -- at R = 8 mm the block is 10 um clear only 0.4 mm from the
   contact. That is why the block here is short: 0.70 mm, so its ends are 7.7 um off
   rather than hundreds. Wanting a long flat block hugging a small wheel is the one
   thing the geometry will not give.
"""

import os
import pickle
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from semgrit.analysis import AnalysisParams
from semgrit.build_deck import (CONTACT_OVERHEAD, PARALLEL_EFFICIENCY, RATE_PER_CORE,
                                DeckParams, build_deck, plan_deck)
from semgrit.preview import summary_text

OUT = "ARC30_CUTS"
NAME = "grinding_wheel_arc30_cuts"

PARAMS = DeckParams(
    name=NAME,
    diameter_mm=16.0, sector_mode="angle", sector_deg=30.0,
    rim_depth_mm=1.20, width_mm=0.60,
    shell_circumferential_divisions=300, shell_axial_divisions=10,
    shell_radial_divisions=1, bond_density_kg_m3=2700.0,
    # whole rim face dressed, edge to edge, at C100
    grit_mode="concentration", concentration=100.0,
    grit_arc_window_mm=0.0, grit_width_window_mm=0.0,
    inset_grit_band=True,
    protrusion_mean=0.85, protrusion_std=0.03,
    protrusion_min=0.80, protrusion_max=0.90,
    max_tilt_deg=35.0, spacing_factor=1.05, seed=20260731,
    # short block, so curvature does not lift its ends away from the wheel
    wp_length_mm=0.70, wp_width_mm=0.55, wp_depth_mm=0.25,
    wp_element_size_mm=0.002,
    wp_element_size_length_mm=0.002, wp_element_size_width_mm=0.004,
    wp_surface_layer_mm=0.030, wp_element_size_depth_mm=0.0015,
    wp_depth_growth=1.3, wp_max_depth_element_mm=0.020,
    wp_material="STONE", wp_density_kg_m3=2650.0,
    wp_youngs_modulus_mpa=50_000.0, wp_poisson_ratio=0.25,
    clearance_um=0.0, wp_position="first grit at entry",
    surface_speed_mm_s=30_000.0, travel_mm=1.50, travel_margin_mm=0.006, cores=64,
    analysis=AnalysisParams(enabled=True, depth_of_cut_um=0.0,
                            mass_scaling_factor=10.0, field_frames=120,
                            restart_intervals=20, history_intervals=400),
    also_write_cae_deck=True,
    write_step=False, write_stl=False,
)


def main():
    solids = pickle.load(open("WHEEL_FIXED/1_measurements/grain_library.pkl",
                              "rb"))["solids"]
    sieve = [s for s in solids if 2.0 <= s.height_um <= 6.6]
    print("grain library : %d measured, %d in the 2.0-6.6 um sieve"
          % (len(solids), len(sieve)))
    os.makedirs(OUT, exist_ok=True)

    print()
    print("=" * 78)
    print("PLAN  (tens of thousands of grains -- this takes a few minutes)")
    print("=" * 78)
    t0 = time.time()
    plan = plan_deck(PARAMS, sieve)
    print("planned in %.0f s" % (time.time() - t0))
    print()
    print(summary_text(plan))

    arc, w = plan["arc_length_mm"], PARAMS.width_mm
    print()
    print("COVERAGE   %.3f x %.3f mm dressed of a %.3f x %.3f mm face = %.1f%%"
          % (plan["grit_band_arc_mm"], plan["grit_band_width_mm"], arc, w,
             100.0 * plan["grit_band_arc_mm"] * plan["grit_band_width_mm"] / (arc * w)))
    # Push the infeed as deep as is safe, since the ask was "as close as possible".
    ceil = plan["depth_ceiling_um"]
    PARAMS.analysis.depth_of_cut_um = round(0.95 * ceil, 4)
    sw = plan["swept_clearances_um"]
    act = sum(1 for x in sw if x <= PARAMS.analysis.depth_of_cut_um)
    R = PARAMS.diameter_mm / 2.0
    print("GAP        standoff 0: the tallest grain under the block is tangent to it.")
    print("           curvature lifts the block %.1f um at +-%.2f mm from the contact"
          % (( (PARAMS.wp_length_mm / 2) ** 2 / (2 * R)) * 1000,
             PARAMS.wp_length_mm / 2))
    print("INFEED     %.3f um = 95%% of the %.3f um face-to-bond gap -> %d of %d "
          "grains cut (%.1f%%)"
          % (PARAMS.analysis.depth_of_cut_um, ceil, act, len(sw),
             100.0 * act / max(len(sw), 1)))
    per_s = RATE_PER_CORE * 64 * PARALLEL_EFFICIENCY / CONTACT_OVERHEAD
    print("RUN        %.2f h on 64 cores, %.1f h on 16"
          % (plan["cost"]["element_increments"] / per_s / 3600.0,
             plan["cost"]["element_increments"]
             / (RATE_PER_CORE * 16 * PARALLEL_EFFICIENCY / CONTACT_OVERHEAD) / 3600.0))

    print()
    print("=" * 78)
    print("BUILDING  (%.0f MB)" % plan["estimated_mb"])
    print("=" * 78)
    t0 = time.time()
    info = build_deck(PARAMS, sieve, OUT)
    print("wrote %s  %.1f MB in %.0f s"
          % (os.path.basename(info["path"]), info["size_bytes"] / 1e6,
             time.time() - t0))

    print()
    print("=" * 78)
    print("VERIFYING")
    print("=" * 78)
    ok = True
    for deck in [info["path"]] + ([info["cae_deck"]] if info.get("cae_deck") else []):
        for v in ("verify_rigid_deck.py", "verify_rigid_deck2.py"):
            r = subprocess.run([sys.executable, v, deck], capture_output=True,
                               text=True)
            fails = [l.strip() for l in r.stdout.splitlines() if "[FAIL]" in l]
            print("%-32s %-22s %s" % (os.path.basename(deck), v,
                                      "PASS" if r.returncode == 0 else "FAIL"))
            for l in fails:
                print("    " + l)
            ok = ok and r.returncode == 0
    print()
    print("ALL VERIFIERS PASS" if ok else "VERIFICATION FAILED - do not run it")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
