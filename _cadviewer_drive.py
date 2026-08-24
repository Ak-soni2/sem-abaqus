"""Click every control in the CAD viewer and assert nothing throws.

    python _cadviewer_drive.py            # small deck
    python _cadviewer_drive.py --dense    # ~5000 grits

`_cadviewer_check.py` proves the page parses; `_cadviewer_shot.py` shows what one
state looks like. Neither exercises the controls, and a viewer's failures are
overwhelmingly "this button throws once the model is loaded" -- a renamed id, a
handler bound before the thing it binds to exists, a null on a deck that happens
to have no workpiece.

So this drives the whole panel in a real browser and fails on the first page
error. It is the gate to run after editing the template.

Needs playwright and its chromium. Development tool; not part of the pipeline.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# (label, javascript). Each runs in the page after the model has loaded. Anything
# that throws is reported with its label, so a failure names the control.
STEPS = [
    ("view iso", "document.querySelector('[data-v=iso]').click()"),
    ("view front", "document.querySelector('[data-v=front]').click()"),
    ("view top", "document.querySelector('[data-v=top]').click()"),
    ("view right", "document.querySelector('[data-v=right]').click()"),
    ("view face", "document.querySelector('[data-v=face]').click()"),
    ("view axial", "document.querySelector('[data-v=axial]').click()"),
    ("view wheel", "document.querySelector('[data-v=wheelview]').click()"),
    ("view contact", "document.querySelector('[data-v=contact]').click()"),
    ("fit", "document.getElementById('cadfit').click()"),

    ("edges off", "const b=document.getElementById('cadedges');"
                  "b.checked=false;b.onchange({target:b})"),
    ("edges on", "const b=document.getElementById('cadedges');"
                 "b.checked=true;b.onchange({target:b})"),
    # ortho/spin handlers read e.target.checked, so a synthetic call passes one.
    ("ortho on", "const b=document.getElementById('cadortho');"
                 "b.checked=true;b.onchange({target:b})"),
    ("ortho off", "const b=document.getElementById('cadortho');"
                  "b.checked=false;b.onchange({target:b})"),
    ("spin on", "const b=document.getElementById('cadspin');b.checked=true"),
    ("spin off", "const b=document.getElementById('cadspin');b.checked=false"),

    ("section X", "const a=document.getElementById('cadaxis');"
                  "a.selectedIndex=1;a.onchange()"),
    ("section slide", "const s=document.getElementById('cadcut');"
                      "s.value=300;s.oninput()"),
    ("section flip", "const f=document.getElementById('cadflip');"
                     "f.checked=true;f.onchange()"),
    ("cap off", "const c=document.getElementById('cadcap');"
                "c.checked=false;c.onchange()"),
    ("cap on", "const c=document.getElementById('cadcap');"
               "c.checked=true;c.onchange()"),
    ("section Z", "const a=document.getElementById('cadaxis');"
                  "a.selectedIndex=2;a.onchange()"),
    ("section off", "const a=document.getElementById('cadaxis');"
                    "a.selectedIndex=0;a.onchange()"),

    ("colour protrusion", "const s=document.getElementById('cadcolour');"
                          "s.value='protrusion_um';s.onchange()"),
    ("colour height", "const s=document.getElementById('cadcolour');"
                      "s.value='height_um';s.onchange()"),
    ("colour width", "const s=document.getElementById('cadcolour');"
                     "s.value='width_um';s.onchange()"),
    ("colour volume", "const s=document.getElementById('cadcolour');"
                      "s.value='volume_um3';s.onchange()"),
    ("colour engage", "const s=document.getElementById('cadcolour');"
                      "s.value='engage';s.onchange()"),
    ("colour off", "const s=document.getElementById('cadcolour');"
                   "s.value='none';s.onchange()"),

    ("explode 50", "const r=document.getElementById('cadexplode');"
                   "r.value=50;r.oninput()"),
    ("explode 100", "const r=document.getElementById('cadexplode');"
                    "r.value=100;r.oninput()"),
    ("explode 0", "const r=document.getElementById('cadexplode');"
                  "r.value=0;r.oninput()"),

    ("parts toggle off", "document.querySelectorAll('#cadtree input')"
                         ".forEach(b=>{b.checked=false;b.onchange({target:b})})"),
    ("parts toggle on", "document.querySelectorAll('#cadtree input')"
                        ".forEach(b=>{b.checked=true;b.onchange({target:b})})"),
    ("bc toggle off", "document.querySelectorAll('#cadbclist input')"
                      ".forEach(b=>{b.checked=false;b.onchange({target:b})})"),
    ("bc toggle on", "document.querySelectorAll('#cadbclist input')"
                     ".forEach(b=>{b.checked=true;b.onchange({target:b})})"),
    ("bc details open", "const d=document.querySelector('#cadbclist details');"
                        "if(d)d.open=true"),

    ("keys open", "document.getElementById('cadkeysbtn').click()"),
    ("keys close", "document.getElementById('cadkeysbtn').click()"),
    ("collapse all", "document.querySelectorAll('#cadtools .cadhd')"
                     ".forEach(h=>h.click())"),
    ("expand all", "document.querySelectorAll('#cadtools .cadhd')"
                   ".forEach(h=>h.click())"),

    ("depth band drag", "const r=document.getElementById('cadbandrange');"
                        "if(r){r.value=String(Math.round(Number(r.max)*0.5));"
                        "r.oninput();r.onchange();}"),
    ("edit a field", "const e=document.getElementById('cadedit_depth_of_cut_um');"
                     "if(e){e.value=0.25;e.onchange();}"),
    ("copy json", "document.getElementById('cadcopy').click()"),
    ("reset edits", "document.getElementById('cadreset').click()"),
    ("arm drag", "document.getElementById('caddrag').click()"),
    ("disarm drag", "document.getElementById('caddrag').click()"),
    ("clear measure", "document.getElementById('cadclear').click()"),
    ("screenshot", "document.getElementById('cadshot').click()"),

    ("api parts", "return document.getElementById('cadwrap').cadAPI.parts()"),
    ("api settings", "return document.getElementById('cadwrap').cadAPI.settings()"),
    ("api bc", "return document.getElementById('cadwrap').cadAPI.bc()"),
]


def main(argv):
    from playwright.sync_api import sync_playwright

    from _cadviewer_check import deck
    from semgrit.cadviewer import build as build_cad

    dense = "--dense" in argv
    modes = ["contact", "wheel", "whole wheel"]
    plan, _ = deck(dense)
    print("deck: %d grits" % plan["n_grits"])

    total_bad = 0
    with sync_playwright() as p:
        br = p.chromium.launch(args=["--use-gl=swiftshader",
                                     "--enable-unsafe-swiftshader"])
        for mode in modes:
            glb = os.path.join(HERE, "_caddrive.glb")
            html, meta, _ = build_cad(plan, glb, mode=mode, max_grits=0,
                                      height=720)
            page = os.path.join(HERE, "_caddrive.html")
            with open(page, "w", encoding="utf-8") as fh:
                fh.write("<!doctype html><meta charset=utf-8>"
                         "<body style='margin:0'>" + html)

            ctx = br.new_context(viewport={"width": 1400, "height": 820},
                                 accept_downloads=True)
            # navigator.clipboard is denied on a file:// origin, and the viewer
            # already handles that by falling back to selecting the textarea.
            # Silence it here so a harness limitation is not read as a defect.
            ctx.add_init_script(
                "Object.defineProperty(navigator,'clipboard',"
                "{get:()=>({writeText:async()=>{}})});")
            pg = ctx.new_page()
            errs = []
            pg.on("pageerror", lambda e: errs.append(str(e)))
            pg.goto("file:///" + page.replace("\\", "/"))
            try:
                pg.wait_for_function(
                    "() => {const s=document.getElementById('cadstat');"
                    "return s && !/loading/i.test(s.textContent);}",
                    timeout=90_000)
            except Exception:
                print("  %-12s FAIL the model never loaded" % mode)
                total_bad += 1
                ctx.close()
                continue
            pg.wait_for_timeout(700)

            bad = []
            for label, js in STEPS:
                before = len(errs)
                try:
                    pg.evaluate("() => {%s}" % js)
                except Exception as exc:                          # noqa: BLE001
                    bad.append("%s: %s" % (label, str(exc).split("\n")[0][:110]))
                    continue
                pg.wait_for_timeout(45)
                if len(errs) > before:
                    bad.append("%s: %s" % (label, errs[-1][:110]))
            ctx.close()

            print("  %-12s %2d steps, %s" % (mode, len(STEPS),
                                             "all OK" if not bad
                                             else "%d FAILED" % len(bad)))
            for b in bad:
                print("      FAIL %s" % b)
            total_bad += len(bad)
        br.close()

    print()
    print("ALL OK" if not total_bad else "%d PROBLEM(S)" % total_bad)
    return 1 if total_bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
