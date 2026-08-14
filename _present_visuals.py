"""Figures for the talk: the 4-panel preview, the glTF, and two rendered views."""
import os
import pickle
import re
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from semgrit.build_deck import plan_deck
from semgrit.cadviewer import build as build_cad_view
from semgrit.preview import preview
from _build_presentation import NAME, OUT, PARAMS

EDGE = [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe"]


def shoot(page, png, button=None, wait=6000):
    exe = next((b for b in EDGE if os.path.exists(b)), None)
    if not exe:
        print("no browser found - skipping %s" % os.path.basename(png))
        return False
    prof = tempfile.mkdtemp(prefix="presshot_")
    try:
        url = "file:///" + os.path.abspath(page).replace("\\", "/").replace(" ", "%20")
        subprocess.run([exe, "--headless=new", "--disable-gpu",
                        "--enable-unsafe-swiftshader", "--no-sandbox",
                        "--user-data-dir=" + prof, "--window-size=1600,1000",
                        "--virtual-time-budget=%d" % (wait + 12000),
                        "--screenshot=" + os.path.abspath(png), url],
                       capture_output=True, timeout=300)
    finally:
        # The probe used to leave these behind; 74 of them had piled up.
        import shutil
        shutil.rmtree(prof, ignore_errors=True)
    return os.path.exists(png)


def main():
    solids = pickle.load(open("WHEEL_FIXED/1_measurements/grain_library.pkl",
                              "rb"))["solids"]
    sieve = [s for s in solids if 2.0 <= s.height_um <= 6.6]
    plan = plan_deck(PARAMS, sieve)

    fig = preview(plan)
    p1 = os.path.join(OUT, NAME + "_preview.png")
    fig.savefig(p1, dpi=140)
    plt.close(fig)
    print("wrote %s" % os.path.basename(p1))

    glb = os.path.join(OUT, NAME + "_cad.glb")
    html, meta, info = build_cad_view(plan, glb, mode="whole wheel", max_grits=0,
                                      height=940, max_inline_mb=24.0)
    print("wrote %s  %.1f MB, %s triangles, %d grains in full detail + %d as boxes"
          % (os.path.basename(glb), info["bytes"] / 1e6,
             format(info["triangles"], ","), len(meta["grains"]),
             len(meta["grains_far"])))
    for n in meta["notes"]:
        print("   note: %s" % n)

    for tag, btn in (("wheel", "wheelview"), ("contact", "contact")):
        pose = ("<script>setTimeout(function(){document.querySelector("
                "'#cadtools button[data-v=\"%s\"]').click();},6000);</script>" % btn)
        page = os.path.join(OUT, "_view_%s.html" % tag)
        with open(page, "w", encoding="utf-8", newline="\n") as fh:
            fh.write('<!doctype html><meta charset="utf-8">'
                     '<body style="margin:0">' + html + pose + "</body>")
        png = os.path.join(OUT, "%s_view_%s.png" % (NAME, tag))
        print("  %s -> %s" % (tag, "ok" if shoot(page, png, btn) else "failed"))
        os.remove(page)
    return 0


if __name__ == "__main__":
    sys.exit(main())
