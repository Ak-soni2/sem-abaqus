"""A wheel segment whose ENTIRE rim face is dressed, at real grit concentration.

The previous deck dressed only the strip the block sweeps, which is cheap and
physically sufficient -- but in CAD it reads as a bare block with a postage stamp of
abrasive on it, because that is what it is. This one dresses the whole periphery, and
pays for it by shrinking the segment: a smaller angle and a narrower face, so the full
surface can carry C100 concentration without the deck becoming unopenable.

  arc 1.20 mm x face 0.80 mm, every square micron of it dressed
  2.75 deg of a 50 mm wheel, rim 0.70 mm deep

Honest consequence, worth saying out loud in the talk: dressing the whole face does not
put every grain in contact. Over a 1.2 mm arc on a 25 mm radius the rim falls away
7 um at the ends, against grain protrusions of about 5 um, so contact is confined to
the middle of the arc. That is exactly what a real wheel does -- the contact arc is a
property of the geometry, not of the dressing -- and the fraction of grains cutting at
any instant is reported below rather than glossed over.
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

OUT = "PRESENT_FULLWHEEL"
NAME = "grinding_wheel_full_face"

PARAMS = DeckParams(
    name=NAME,
    # A small segment, so the whole of it can be dressed
    diameter_mm=50.0, sector_mode="angle", sector_deg=2.75,
    rim_depth_mm=0.70, width_mm=0.80,
    shell_circumferential_divisions=240, shell_axial_divisions=10,
    shell_radial_divisions=1, bond_density_kg_m3=2700.0,
    # 0 and 0 mean: dress the entire rim face, edge to edge
    grit_mode="concentration", concentration=100.0,
    grit_arc_window_mm=0.0, grit_width_window_mm=0.0,
    inset_grit_band=True,
    protrusion_mean=0.85, protrusion_std=0.03,
    protrusion_min=0.80, protrusion_max=0.90,
    max_tilt_deg=35.0, spacing_factor=1.05, seed=20260731,
    # workpiece: as wide as the dressed face, so nothing is wasted off the side
    wp_length_mm=0.90, wp_width_mm=0.70, wp_depth_mm=0.30,
    wp_element_size_mm=0.002,
    wp_element_size_length_mm=0.002, wp_element_size_width_mm=0.004,
    wp_surface_layer_mm=0.030, wp_element_size_depth_mm=0.0015,
    wp_depth_growth=1.3, wp_max_depth_element_mm=0.020,
    wp_material="STONE", wp_density_kg_m3=2650.0,
    wp_youngs_modulus_mpa=50_000.0, wp_poisson_ratio=0.25,
    clearance_um=0.0, wp_position="centred",
    surface_speed_mm_s=30_000.0, travel_mm=0.25, travel_margin_mm=0.006, cores=64,
    analysis=AnalysisParams(enabled=True, depth_of_cut_um=0.0,
                            mass_scaling_factor=10.0, field_frames=120,
                            restart_intervals=20, history_intervals=400),
    also_write_cae_deck=True,
    write_step=True, write_stl=True, step_max_grains=0, stl_max_grains=0,
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
    print("PLAN")
    print("=" * 78)
    plan = plan_deck(PARAMS, sieve)
    print(summary_text(plan))
    arc, w = plan["arc_length_mm"], PARAMS.width_mm
    print()
    print("COVERAGE   dressed band %.3f x %.3f mm of a %.3f x %.3f mm rim face "
          "= %.1f%% of it"
          % (plan["grit_band_arc_mm"], plan["grit_band_width_mm"], arc, w,
             100.0 * plan["grit_band_arc_mm"] * plan["grit_band_width_mm"]
             / (arc * w)))
    sw = plan["swept_clearances_um"]
    act = sum(1 for x in sw if x <= plan["depth_of_cut_um"])
    per_s = RATE_PER_CORE * 64 * PARALLEL_EFFICIENCY / CONTACT_OVERHEAD
    print("ENGAGEMENT %d of %d grains cut at %.2f um (%.1f%%) -- the rest are outside "
          "the contact arc" % (act, len(sw), plan["depth_of_cut_um"],
                               100.0 * act / max(len(sw), 1)))
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
