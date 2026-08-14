"""Build the presentation deck: a dressed wheel sector at real grit density.

Design decisions and why, so they can be defended in the talk:

* **Grit density is C100 concentration**, which for this measured grain library works
  out at ~14,900 grains/mm2 -- a mean spacing of under two grain widths. That is a
  genuinely dressed abrasive surface, not a sparse cartoon.
* **The grit is sieved to 2.0-6.6 um tall.** A real wheel is made from graded abrasive
  of one mesh size; the raw SEM library spans 0.5-6.5 um because segmentation catches
  every fragment. Sieving is more faithful, and it is also what makes the grain tips
  land close enough together to all reach the work.
* **The dressed strip is sized to the contact arc, not to the block.** Over a 1.6 mm
  block on a 25 mm radius the rim falls away 5 um at the ends -- far more than the
  ~2 um spread in grain protrusion -- so who can touch is set by wheel curvature. The
  strip is kept inside the arc that a 5 um infeed closes, which is what puts two
  thirds of the grains into cut at once. A wider strip would only add grains that
  never touch.
* **The workpiece is millimetre scale** (1.60 x 0.70 x 0.35 mm) and thick enough not
  to behave like a membrane.
* **Zero standoff**: the tallest grain under the block is exactly tangent to it, and
  the 5 um infeed then drives the rest in.
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

OUT = "PRESENT_FINAL"
NAME = "grinding_wheel_dressed_10deg"

PARAMS = DeckParams(
    name=NAME,
    # ---- the wheel: a small sector, but a solid wedge of one ----
    diameter_mm=50.0, sector_mode="angle", sector_deg=10.0,
    rim_depth_mm=3.0, width_mm=1.5,
    shell_circumferential_divisions=240, shell_axial_divisions=12,
    shell_radial_divisions=1, bond_density_kg_m3=2700.0,
    # ---- the abrasive: real concentration, on the strip that can cut ----
    grit_mode="concentration", concentration=100.0,
    grit_arc_window_mm=0.95, grit_width_window_mm=0.62,
    inset_grit_band=True,
    protrusion_mean=0.85, protrusion_std=0.03,
    protrusion_min=0.80, protrusion_max=0.90,
    max_tilt_deg=35.0, spacing_factor=1.05, seed=20260731,
    # ---- the workpiece: millimetre scale, and thick ----
    wp_length_mm=1.60, wp_width_mm=0.70, wp_depth_mm=0.35,
    wp_element_size_mm=0.002,
    wp_element_size_length_mm=0.002,      # along the cut
    wp_element_size_width_mm=0.005,       # across the face
    wp_surface_layer_mm=0.030,            # fine skin where the chips form
    wp_element_size_depth_mm=0.0015,
    wp_depth_growth=1.3, wp_max_depth_element_mm=0.020,
    wp_material="STONE", wp_density_kg_m3=2650.0,
    wp_youngs_modulus_mpa=50_000.0, wp_poisson_ratio=0.25,
    clearance_um=0.0, wp_position="centred",
    # ---- kinematics ----
    surface_speed_mm_s=30_000.0, travel_mm=0.30, travel_margin_mm=0.006, cores=64,
    analysis=AnalysisParams(enabled=True, depth_of_cut_um=5.0,
                            mass_scaling_factor=10.0, field_frames=120,
                            restart_intervals=20, history_intervals=400),
    also_write_cae_deck=True,
    write_step=True, write_stl=True, step_max_grains=1200, stl_max_grains=1200,
)


def main():
    solids = pickle.load(open("WHEEL_FIXED/1_measurements/grain_library.pkl",
                              "rb"))["solids"]
    # Sieve the library the way a wheel manufacturer grades abrasive.
    sieve = [s for s in solids if 2.0 <= s.height_um <= 6.6]
    print("grain library : %d measured, %d in the 2.0-6.6 um sieve" % (len(solids),
                                                                       len(sieve)))
    os.makedirs(OUT, exist_ok=True)

    print()
    print("=" * 78)
    print("PLAN")
    print("=" * 78)
    plan = plan_deck(PARAMS, sieve)
    print(summary_text(plan))
    sw = plan["swept_clearances_um"]
    act = sum(1 for x in sw if x <= plan["depth_of_cut_um"])
    per_s = RATE_PER_CORE * 64 * PARALLEL_EFFICIENCY / CONTACT_OVERHEAD
    print()
    print("ENGAGEMENT %d of %d grains in the swept band cut at %.2f um (%.1f%%)"
          % (act, len(sw), plan["depth_of_cut_um"], 100.0 * act / max(len(sw), 1)))
    print("RUN        %.2f h on 64 cores, %.1f h on 16"
          % (plan["cost"]["element_increments"] / per_s / 3600.0,
             plan["cost"]["element_increments"]
             / (RATE_PER_CORE * 16 * PARALLEL_EFFICIENCY / CONTACT_OVERHEAD) / 3600.0))

    print()
    print("=" * 78)
    print("BUILDING  (%.0f MB, this takes a few minutes)" % plan["estimated_mb"])
    print("=" * 78)
    t0 = time.time()
    info = build_deck(PARAMS, sieve, OUT)
    print("wrote %s  %.1f MB in %.0f s"
          % (os.path.basename(info["path"]), info["size_bytes"] / 1e6,
             time.time() - t0))
    for k in ("cae_deck", "postprocess_script"):
        if info.get(k):
            print("      %s" % os.path.basename(info[k]))

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
            print("%-34s %-22s %s" % (os.path.basename(deck), v,
                                      "PASS" if r.returncode == 0 else "FAIL"))
            for l in fails:
                print("    " + l)
            ok = ok and r.returncode == 0
    print()
    print("ALL VERIFIERS PASS" if ok else "VERIFICATION FAILED - do not run it")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
