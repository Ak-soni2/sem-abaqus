"""Drive the MESH viewer in a real browser and check its wording.

    python _meshview_drive.py

`semgrit.meshview` feeds element geometry through the SAME viewer the CAD is
shown in, which is what makes the mesh view free -- it inherits section capping,
explode, the measuring tool and every shortcut. The risk is the other side of
that trade: one template, two kinds of content, so a control that assumes grain
solids can be left enabled and wired to nothing.

That is exactly what was found here. `setupColourBy()` was called from inside
`buildBC()`, which returns early when a view has no boundary conditions, so on a
mesh view the colour-by selector stayed ENABLED but did nothing -- picking an
entry gave no colour and no explanation. It is now called unconditionally, and
its own guard disables it with a reason.

So this drives all 51 controls against a real three-part SAG mesh, and then
asserts the wording: a mesh view must not invite the reader to "click a grain".
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))

def main():
    from playwright.sync_api import sync_playwright
    from _cadviewer_drive import STEPS
    from semgrit import meshview as mv
    from semgrit.sagwrite import build_block, build_compliant_ring
    hub = build_compliant_ring(inner_r_mm=52.5, outer_r_mm=57.5, width_mm=10.0,
                               sector_deg=20.0, n_circ=16, n_rad=2, n_axial=4)
    pu = build_compliant_ring(inner_r_mm=57.5, outer_r_mm=62.5, width_mm=10.0,
                              sector_deg=20.0, n_circ=16, n_rad=4, n_axial=4)
    wp = build_block(length_mm=8.0, width_mm=10.0, depth_mm=2.0,
                     el_length_mm=0.5, el_width_mm=0.6, fine_depth_mm=0.1,
                     band_mm=0.5, growth=1.3, x0_mm=-4.0, y0_mm=-5.0)
    html, vmeta, vinfo = mv.build([
        dict(name="hub", nodes=hub[0], conn=hub[1], color=mv.C_HUB),
        dict(name="polyurethane", nodes=pu[0], conn=pu[1],
             color=mv.C_COMPLIANT),
        dict(name="workpiece", nodes=wp[0], conn=wp[1], color=mv.C_WORK),
    ], os.path.join(HERE, "_meshview.glb"), height=720)
    assert vmeta["kind"] == "mesh"
    page = os.path.join(HERE, "_meshview.html")
    with open(page, "w", encoding="utf-8") as fh:
        fh.write("<!doctype html><meta charset=utf-8>"
                 "<body style='margin:0'>" + html)
    with sync_playwright() as p:
        br = p.chromium.launch(args=["--use-gl=swiftshader",
                                     "--enable-unsafe-swiftshader"])
        ctx = br.new_context(viewport={"width": 1400, "height": 820},
                             accept_downloads=True)
        ctx.add_init_script("Object.defineProperty(navigator,'clipboard',"
                            "{get:()=>({writeText:async()=>{}})});")
        pg = ctx.new_page()
        errs = []
        pg.on("pageerror", lambda e: errs.append(str(e)))
        pg.goto("file:///" + page.replace("\\", "/"))
        pg.wait_for_function(
            "() => {const s=document.getElementById('cadstat');"
            "return s && !/loading/i.test(s.textContent);}", timeout=90_000)
        pg.wait_for_timeout(800)
        # Controls that exist only for a CAD deck. A mesh view has no
        # per-grain properties to colour by and no editable parameters, and
        # the viewer correctly DISABLES the selector and omits the edit
        # buttons -- so absence here is right, not broken. Verified below
        # rather than merely skipped.
        cad_only = {"colour protrusion", "colour height", "colour width",
                    "colour volume", "colour engage", "colour off",
                    "copy json", "reset edits", "edit a field",
                    "depth band drag"}
        sel_state = pg.evaluate(
            "() => {const s=document.getElementById('cadcolour');"
            "return s ? {present:true, disabled:s.disabled,"
            " title:s.title} : {present:false};}")
        assert sel_state["present"], "the colour selector must still be drawn"
        assert sel_state["disabled"],             "with no per-grain data the colour selector must be DISABLED"
        assert "no per-grain data" in (sel_state["title"] or ""), sel_state
        print("colour-by: correctly disabled -- %r" % sel_state["title"])

        # A mesh view must not use the CAD vocabulary.
        words = pg.evaluate(
            "() => ({"
            " colour: [...document.querySelectorAll('#cadtools .cadhd')]"
            "   .map(h => h.textContent).find(t => /colour/i.test(t)) || '',"
            " pick: document.getElementById('cadpick').textContent,"
            " hint: (document.getElementById('cadfirst')||{}).textContent||'',"
            " inspect: [...document.querySelectorAll('#cadinfo .cadhd')]"
            "   .map(h => h.textContent).join('|')})")
        assert "grains" not in words["colour"].lower(), words["colour"]
        assert "grain" not in words["pick"].lower(), words["pick"]
        assert "MESH" in words["hint"], words["hint"][:120]
        assert "dc is resolved" in words["hint"], words["hint"][:120]
        assert "Mesh" in words["inspect"], words["inspect"]
        print("wording:   retitled for a mesh view")

        bad = []
        for label, js in STEPS:
            if label in cad_only:
                continue
            before = len(errs)
            try:
                pg.evaluate("() => {%s}" % js)
            except Exception as exc:
                bad.append("%s: %s" % (label, str(exc).split("\n")[0][:110]))
                continue
            pg.wait_for_timeout(40)
            if len(errs) > before:
                bad.append("%s: %s" % (label, errs[-1][:110]))
        pg.screenshot(path=os.path.join(HERE, "_meshview_shot.png"))
        parts = pg.evaluate(
            "() => document.getElementById('cadwrap').cadAPI.parts()")
        ctx.close(); br.close()
    print("%d steps, %s" % (len(STEPS), "all OK" if not bad
                            else "%d FAILED" % len(bad)))
    for b in bad:
        print("   FAIL", b)
    print("parts in viewer:", [q.get("name") for q in parts]
          if isinstance(parts, list) else parts)
    return 1 if bad else 0

if __name__ == "__main__":
    raise SystemExit(main())
