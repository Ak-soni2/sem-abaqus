"""Render the CAD viewer in a real browser, screenshot it, and report its console.

    python _cadviewer_shot.py                 # default deck, iso view
    python _cadviewer_shot.py --dense         # ~5000 grits
    python _cadviewer_shot.py --act contact   # click a control first

The static check (`_cadviewer_check.py`) proves the page parses. It cannot tell
whether the thing actually draws, whether a panel overflows its column, or
whether a control throws the moment it is clicked -- and those are exactly the
failures a viewer has. This drives headless Chromium instead: load the page,
wait for the model, run an optional interaction, then save a PNG and print every
console message and page error.

Needs `playwright` and its chromium (`python -m playwright install chromium`).
It is a development tool, not part of the pipeline or the notebook.
"""
from __future__ import annotations

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

# What to do after the model loads, before the shot. Each is a small script run
# in the page; they exercise the controls the way a user would.
ACTIONS = {
    "none": "",
    "contact": "document.querySelector('[data-v=contact]').click()",
    "wheel": "document.querySelector('[data-v=wheelview]').click()",
    "face": "document.querySelector('[data-v=face]').click()",
    "section": ("document.getElementById('cadaxis').selectedIndex=1;"
                "document.getElementById('cadaxis').onchange();"
                "document.getElementById('cadcut').value=420;"
                "document.getElementById('cadcut').oninput();"),
    "ortho": ("document.getElementById('cadortho').checked=true;"
              "document.getElementById('cadortho').onchange();"),
    "noedges": ("const b=document.getElementById('cadedges');b.checked=false;"
                "b.onchange({target:b});"),
    "keys": "document.getElementById('cadkeys').style.display='block'",
    # The first-run hint dismisses itself on the first pointerdown anywhere in the
    # wrapper, and merely moving the mouse over the page during a screenshot run is
    # enough to lose it. Re-show it rather than chase the timing.
    "hint": ("const f=document.getElementById('cadfirst');"
             "if(f) f.style.display='block';"),
    "collapse": ("document.querySelectorAll('#cadtools .cadhd').forEach((h,i)=>"
                 "{if(i>2)h.click();})"),
    "colour": ("const s=document.getElementById('cadcolour');"
               "s.value='protrusion_um';s.onchange();"
               # hide the bond rim: it is a solid wall in front of the grains
               "document.querySelector('#cadtree input').click();"
               "document.getElementById('cadfit').click();"
               "const f=document.getElementById('cadfirst');if(f)f.remove();"),
    "colour_face": ("const s=document.getElementById('cadcolour');"
                    "s.value='protrusion_um';s.onchange();"
                    "document.querySelector('[data-v=contact]').click();"
                    "document.querySelector('[data-v=face]').click();"
                    "const f=document.getElementById('cadfirst');if(f)f.remove();"),
    "colour_engage": ("const s=document.getElementById('cadcolour');"
                      "s.value='engage';s.onchange();"),
    "band": ("const r=document.getElementById('cadbandrange');"
             "if(r){r.value=String(Math.round(Number(r.max)*0.45));"
             "r.oninput();r.onchange();}"),
}


def main(argv):
    from playwright.sync_api import sync_playwright

    from _cadviewer_check import deck
    from semgrit.cadviewer import build as build_cad

    dense = "--dense" in argv
    mode = "whole wheel"
    if "--mode" in argv:
        mode = argv[argv.index("--mode") + 1]
    act = argv[argv.index("--act") + 1] if "--act" in argv else "none"
    if act not in ACTIONS:
        raise SystemExit("--act must be one of: %s" % ", ".join(ACTIONS))

    plan, _ = deck(dense)
    glb = os.path.join(HERE, "_cadshot.glb")
    html, meta, info = build_cad(plan, glb, mode=mode, max_grits=0, height=760)
    page_path = os.path.join(HERE, "_cadshot.html")
    with open(page_path, "w", encoding="utf-8") as fh:
        fh.write("<!doctype html><meta charset=utf-8>"
                 "<body style='margin:0;background:#fff'>" + html)

    out = os.path.join(HERE, "_cadshot_%s_%s.png"
                       % (mode.replace(" ", "_"), act))
    msgs, errors = [], []
    with sync_playwright() as p:
        br = p.chromium.launch(args=["--use-gl=swiftshader",
                                     "--enable-unsafe-swiftshader"])
        # A fresh context every run, so localStorage is empty. Otherwise the
        # first-run hint and the collapsed-section memory persist between shots
        # and the screenshot stops showing what a new user actually sees.
        ctx = br.new_context(viewport={"width": 1560, "height": 820},
                             device_scale_factor=2)
        pg = ctx.new_page()
        pg.on("console", lambda m: msgs.append("%s: %s" % (m.type, m.text)))
        pg.on("pageerror", lambda e: errors.append(str(e)))
        pg.goto("file:///" + page_path.replace("\\", "/"))

        # The status line stops saying "loading" once the model is in.
        try:
            pg.wait_for_function(
                "() => {const s=document.getElementById('cadstat');"
                "return s && !/loading/i.test(s.textContent);}", timeout=60_000)
        except Exception:
            errors.append("the model never finished loading")
        pg.wait_for_timeout(1200)

        if ACTIONS[act]:
            try:
                pg.evaluate("() => {%s}" % ACTIONS[act])
            except Exception as exc:                              # noqa: BLE001
                errors.append("action %r threw: %s" % (act, exc))
            pg.wait_for_timeout(900)

        stat = pg.text_content("#cadstat") or ""
        pg.screenshot(path=out)
        br.close()

    print("mode=%s  action=%s" % (mode, act))
    print("  grits %d/%d (%d full)   page %.2f MB   tris %s"
          % (meta["grits_drawn"], meta["grits_total"], meta["grits_full_detail"],
             len(html) / 1e6, format(info["triangles"], ",")))
    print("  status line: %s" % " ".join(stat.split())[:150])
    bad = [m for m in msgs if m.startswith(("error", "warning"))]
    for m in bad[:12]:
        print("  console %s" % m[:150])
    for e in errors:
        print("  PAGE ERROR %s" % e[:400])
    print("  wrote %s" % os.path.relpath(out, HERE))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
