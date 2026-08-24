"""Build the CAD viewer on a real deck and check the page holds together.

    python _cadviewer_check.py            # a small deck, fast
    python _cadviewer_check.py --dense    # a 30 deg wheel, ~5000 grits
    python _cadviewer_check.py --open     # write the html and open it in a browser

The viewer is ~1700 lines of JavaScript inside a Python string, built once per
notebook run and exercised only by a human looking at it. Nothing else in this
repo can catch a typo in it, and a syntax error shows up as a blank grey box in
front of whoever is being shown the model.

So this is the cheap gate: build the page for real, then assert the things that
break silently -- every control the template promises is present, every id the
script reaches for exists, the percent-escaping survived, and the JS has
balanced braces. It is not a browser and does not pretend to be; it catches the
class of mistake that comes from editing a large string by hand.
"""
from __future__ import annotations

import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def deck(dense: bool):
    """A real plan, from real measured grains."""
    from semgrit import materials
    from semgrit.analysis import AnalysisParams
    from semgrit.build_deck import DeckParams, plan_deck
    from semgrit.hybrid import HYBRID_DEPVAR
    from semgrit.quick import measure_images

    imgs = ["B4C_15.tif"] + (["B4C_16.tif", "DIAMOND_07.tif"] if dense else [])
    got = measure_images([os.path.join(HERE, i) for i in imgs],
                         os.path.join(HERE, "_cadcheck_meas"),
                         log=lambda *a: None)
    solids = got["solids"]

    hp = materials.hybrid_params("sandstone", h_source=0, dc_form=2)
    dc = hp.critical_depth_mm()
    if dense:
        geom = dict(sector_mode="angle", sector_deg=30.0, rim_depth_mm=1.0,
                    width_mm=10.0, grit_mode="areal_density",
                    areal_density_per_mm2=40.0)
    else:
        geom = dict(sector_mode="arc", arc_length_mm=2.0, rim_depth_mm=0.012,
                    width_mm=0.030, grit_mode="single", single_grain_index=-1,
                    single_grit_offset_mm=0.015)
    p = DeckParams(
        name="cadcheck", diameter_mm=50.0, include_bond=True,
        include_workpiece=True, wp_length_mm=0.048, wp_width_mm=0.020,
        wp_depth_mm=0.006, wp_element_size_length_mm=0.0003,
        wp_element_size_width_mm=0.0003, wp_element_size_depth_mm=dc / 5.0,
        clearance_um=0.0, wp_position="centred", surface_speed_mm_s=30_000.0,
        cores=8,
        analysis=AnalysisParams(enabled=True, depth_of_cut_um=0.20,
                                material_model="hybrid", hybrid=hp,
                                n_depvar=HYBRID_DEPVAR, element_deletion=True),
        **geom)
    materials.apply(p, "sandstone")
    return plan_deck(p, solids), solids


# Every control the template advertises. If one of these ids stops being emitted
# the feature is gone from the UI and nothing else would notice.
REQUIRED_IDS = [
    "cadwrap", "cadcanvas", "cadtools", "cadinfo", "cadstat", "cadhelp",
    "cadedges", "cadortho", "cadspin", "cadaxis", "cadcut", "cadflip",
    "cadclear", "cadtree", "cadbclist", "cadeditlist", "cadshot", "cadfit",
    "caddrag", "cadpick",
    # added by the upgrade
    "cadfull", "cadgrip", "cadfirst", "cadkeys", "cadkeysbtn", "cadwhere",
    "cadcolour", "cadlegend", "cadlegbar", "cadleglo", "cadleghi",
    "cadcap", "cadexplode", "cadexplodetxt",
]

# Ids the script reaches for that are created at runtime rather than in the HTML.
RUNTIME_IDS = ["cadbandrange", "cadbandtrack", "cadbandpin", "cadbandtext",
               "cadapply", "cadcopy", "caddl", "cadreset", "cadeditstate",
               "cadeditjson"]


def check(html: str, meta: dict, info: dict, tag: str) -> int:
    bad = 0

    def fail(msg):
        nonlocal bad
        bad += 1
        print("  FAIL  %s" % msg)

    for i in REQUIRED_IDS:
        if ('id="%s"' % i) not in html:
            fail("control missing from the page: #%s" % i)
    for i in RUNTIME_IDS:
        if i not in html:
            fail("runtime id never referenced: #%s" % i)

    # A literal percent in the TEMPLATE that lost its double becomes a stray format
    # spec. Checking the rendered page cannot see that -- '120% 100%' is correct
    # output -- so check the template source instead, where every literal percent
    # must be doubled and every real placeholder is %(name)s.
    from semgrit import cadviewer as _cv
    tpl = _cv._TEMPLATE
    for m in re.finditer(r"(?<!%)%(?!%|\()", tpl):
        frag = tpl[max(0, m.start() - 45):m.start() + 15].replace("\n", " ")
        fail("undoubled %% in _TEMPLATE near: ...%s..." % frag)
        break

    # Parse the module script with node. Brace counting cannot be made correct --
    # regex literals and template strings defeat it -- and a real parse is both
    # simpler and stricter. Skipped, with a note, where node is not installed.
    js = html.split('<script type="module">', 1)
    if len(js) != 2:
        fail("the module script is not in the page")
    else:
        body = js[1].rsplit("</script>", 1)[0]
        import shutil
        import subprocess
        import tempfile
        node = shutil.which("node")
        if not node:
            print("       (node not found: the JS was not parsed)")
        else:
            fd, path = tempfile.mkstemp(suffix=".mjs")
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as fh:
                    fh.write(body)
                r = subprocess.run([node, "--check", path],
                                   capture_output=True, text=True)
                if r.returncode != 0:
                    fail("the viewer JavaScript does not parse:\n%s"
                         % (r.stderr or r.stdout)[:900])
            finally:
                os.unlink(path)

    if "data:model/gltf-binary;base64," not in html:
        fail("the model is not embedded")

    print("  %-13s %6.2f MB page  %6.2f MB glb  %s tris  %d/%d grits (%d full)  %s"
          % (tag, len(html) / 1e6, info["bytes"] / 1e6,
             format(info["triangles"], ","), meta["grits_drawn"],
             meta["grits_total"], meta["grits_full_detail"],
             "OK" if not bad else "%d PROBLEM(S)" % bad))
    for n in meta.get("notes", []):
        print("                note: %s" % n[:96])
    return bad


def main(argv):
    dense = "--dense" in argv
    from semgrit.cadviewer import build as build_cad

    print("building a %s deck ..." % ("dense 30 deg wheel" if dense else "single-grit"))
    t0 = time.time()
    plan, _ = deck(dense)
    print("  plan: %d grits, %.1f s" % (plan["n_grits"], time.time() - t0))

    bad = 0
    pages = {}
    for mode in ("contact", "wheel", "whole wheel"):
        glb = os.path.join(HERE, "_cadcheck_%s.glb" % mode.replace(" ", "_"))
        t0 = time.time()
        try:
            html, meta, info = build_cad(plan, glb, mode=mode, max_grits=0,
                                         height=720)
        except Exception as exc:                                  # noqa: BLE001
            print("  FAIL  %s: %s" % (mode, exc))
            bad += 1
            continue
        bad += check(html, meta, info, mode)
        pages[mode] = html

    if "--open" in argv and pages:
        mode = "whole wheel" if "whole wheel" in pages else list(pages)[0]
        out = os.path.join(HERE, "_cadcheck_view.html")
        with open(out, "w", encoding="utf-8") as fh:
            fh.write("<!doctype html><meta charset=utf-8>" + pages[mode])
        print("\nwrote %s" % out)
        import webbrowser
        webbrowser.open("file:///" + out.replace("\\", "/"))

    print()
    print("ALL OK" if not bad else "%d PROBLEM(S)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
