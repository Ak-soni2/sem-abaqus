"""Drive the three.js CAD viewer in a real headless browser and check what it does.

Static checks cannot tell you that the viewer *works*: whether three.js loads inside a
sandboxed output frame, whether the .glb parses, whether the section plane actually
cuts, and -- the one that bit already -- whether clicking a grain identifies it, given
that GLTFLoader renames every node it loads.

So this loads the generated page in headless Edge/Chrome, waits for the model, then
drives it with real events -- view buttons, wheel-zoom, pointer clicks, the tree
checkboxes, the section slider -- and reads back what the panels say. Everything it
reports is what a user would see.

Zooming before sampling matters: a grain is a few microns across on a 1 mm arc, so from
a whole-wheel view a click grid lands on it about once in a thousand tries. That is a
limit of the probe, not of the viewer, so the probe zooms in first.
"""
from __future__ import annotations

import json
import math
import os
import pickle
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

BROWSERS = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
]

PROBE = r"""
<script>
const OUT = document.createElement('div');
OUT.id = 'probeout'; OUT.style.display = 'none';
document.body.appendChild(OUT);
const sleep = ms => new Promise(r => setTimeout(r, ms));
const $ = s => document.querySelector(s);

// Two animation frames, not a timer. The viewer draws in requestAnimationFrame, so a
// setTimeout can return before anything has been redrawn -- reading pixels then
// compares two copies of the same stale frame and every "did it change?" check
// silently answers no.
// ...but headless Chrome does not always service rAF, so never block on it: race each
// wait against a timer. RAF counts how often the frame path actually won, which tells
// us whether the pixel comparisons below mean anything.
let RAF = 0, TMO = 0;
const frame = () => new Promise(r => {
  let done = false;
  const fin = which => { if (!done) { done = true; which(); r(); } };
  const t = setTimeout(() => fin(() => TMO++), 150);
  requestAnimationFrame(() => requestAnimationFrame(
    () => fin(() => { RAF++; clearTimeout(t); })));
});
async function settle(n) {
  for (let i = 0; i < (n || 3); i++) {
    const api = $('#cadwrap').cadAPI;
    if (api) api.draw();          // step damping even if no animation frame arrives
    await frame();
  }
}

function click(x, y, shift) {
  const c = $('#cadcanvas'), r = c.getBoundingClientRect();
  const o = {clientX: r.left + x, clientY: r.top + y, bubbles: true, button: 0,
             pointerId: 1, pointerType: 'mouse', isPrimary: true, shiftKey: !!shift};
  c.dispatchEvent(new PointerEvent('pointerdown', o));
  // Release it. OrbitControls tracks pressed pointers, and a pointerdown with no
  // matching pointerup leaves a phantom in its list; after enough of them it decides
  // a multi-touch gesture is under way and pans the camera on its own.
  c.dispatchEvent(new PointerEvent('pointerup', o));
  document.dispatchEvent(new PointerEvent('pointerup', o));
  return $('#cadpick').innerText.replace(/\s+/g, ' ').trim();
}

// Records the screen positions that actually landed on the model, so later steps can
// aim at geometry instead of guessing coordinates and reporting a miss as a failure.
function sweep(seen, nx, ny, hits) {
  const c = $('#cadcanvas'), W = c.clientWidth, H = c.clientHeight;
  for (let i = 1; i < nx; i++) {
    for (let j = 1; j < ny; j++) {
      const x = W * i / nx, y = H * j / ny;
      const t = click(x, y, false);
      if (t && t !== 'click a grain') {
        seen[t] = (seen[t] || 0) + 1;
        if (hits) hits.push([x, y]);
      }
    }
  }
}

function zoom(n) {
  const c = $('#cadcanvas');
  for (let k = 0; k < n; k++) {
    c.dispatchEvent(new WheelEvent('wheel',
      {deltaY: -120, bubbles: true, cancelable: true}));
  }
}

(async () => {
  const err = [];
  // OrbitControls calls setPointerCapture, which throws for synthetic pointer events
  // that carry no real pointer id. That is an artefact of driving the page from a
  // script, not a viewer fault, so record distinct messages only.
  addEventListener('error', e => {
    const m = String(e.message);
    if (m.indexOf('setPointerCapture') >= 0) return;
    if (err.indexOf(m) < 0) err.push(m);
  });
  $('#cadcanvas').addEventListener('webglcontextlost',
                                   () => err.push('WEBGL CONTEXT LOST'));
  // Compare whole images, not their compressed length: two different frames can
  // encode to the same number of bytes and a length test then reports "no change".
  // Force the redraw through the viewer's own hook rather than hoping an animation
  // frame has run -- headless Chrome services rAF erratically, and reading the canvas
  // before it redraws compares two copies of the same stale frame.
  const shot = () => {
    const api = $('#cadwrap').cadAPI;
    return api ? api.png() : $('#cadcanvas').toDataURL();
  };
  const stages = {};
  // Publish after every stage. The browser dumps the DOM when its virtual-time budget
  // runs out, which can land mid-run; writing as we go means a truncated run still
  // reports what it managed to check instead of nothing at all.
  const R = {picks: {}, isolated: {}, stages: stages, errors: err, tree: []};
  const save = () => { OUT.textContent = JSON.stringify(R); };
  save();
  try {
    // Wait for the model, do not guess at a delay. three.js comes from a CDN and the
    // .glb has to parse; on a loaded machine that takes longer than any fixed sleep,
    // and every check after it then measures a page that had not finished loading.
    for (let i = 0; i < 200; i++) {
      if ($('#cadwrap').cadAPI &&
          $('#cadstat').innerText.indexOf('rendered locally') >= 0) break;
      await sleep(150);
    }
    stages.ready = !!$('#cadwrap').cadAPI;
    await settle(6);
    const c = $('#cadcanvas'), W = c.clientWidth, H = c.clientHeight;
    R.canvas = [W, H];
    R.status = $('#cadstat').innerText;
    R.tree = Array.from(document.querySelectorAll('#cadtree label'))
                  .map(l => l.innerText.trim());
    save();
    const seen = R.picks;
    sweep(seen, 44, 30);                                 // default iso view
    stages.iso = Object.keys(seen).length;
    save();

    $('#cadtools button[data-v="face"]').click();        // straight at the face
    await settle(4);
    sweep(seen, 44, 30);
    stages.face = Object.keys(seen).length;
    save();

    zoom(14);                                            // as a user would
    await settle(24);                                    // damping needs frames
    const hits = [];
    sweep(seen, 52, 36, hits);
    stages.zoomed = Object.keys(seen).length;
    stages.hits = hits.length;
    stages.blank_after_zoom = shot().length < 3000;
    save();

    // measurement, between two points the sweep just proved are on the model, and
    // without moving the camera in between
    R.measurement = '';
    R.components = '';
    for (let k = 0; k + 1 < hits.length && !R.measurement; k++) {
      $('#cadclear').click();
      click(hits[k][0], hits[k][1], true);
      click(hits[hits.length - 1 - k][0], hits[hits.length - 1 - k][1], true);
      const h = $('#cadhelp').innerText;
      if (h.indexOf('distance') === 0) {
        R.measurement = h;
        R.components = $('#cadwrap').cadAPI.measure();
        R.panel = $('#cadpick').innerText.replace(/\s+/g, ' ').trim();
      }
    }
    $('#cadclear').click();
    save();

    // the two new standard views must move the camera somewhere different
    const shots = {};
    for (const v of ['wheelview', 'contact']) {
      $('#cadtools button[data-v="' + v + '"]').click();
      await settle(24);
      shots[v] = shot();
    }
    R.wheel_vs_contact_differ = shots.wheelview !== shots.contact;
    R.wheel_blank = shots.wheelview.length < 3000;
    R.contact_blank = shots.contact.length < 3000;
    save();

    // From here on, pin the camera to a known view so each check is comparable.
    $('#cadtools button[data-v="face"]').click();
    await settle(24);

    // section plane: it must change the rendered pixels
    const before = shot();
    $('#cadaxis').value = '1'; $('#cadaxis').dispatchEvent(new Event('change'));
    $('#cadcut').value = '680'; $('#cadcut').dispatchEvent(new Event('input'));
    await settle(4);
    R.section_changed = before !== shot();
    $('#cadaxis').value = '-1'; $('#cadaxis').dispatchEvent(new Event('change'));
    await settle(4);
    save();

    // edges toggle must change the picture too
    const e0 = shot();
    $('#cadedges').click();
    await settle(4);
    R.edges_changed = e0 !== shot();
    $('#cadedges').click();
    await settle(4);
    save();

    // the guardrails and a drag, with no kernel at all
    R.guard = null; R.drag = null;
    try {
      const api = $('#cadwrap').cadAPI;
      if (api && api.guard && api.guard().window) {
        const g0 = api.guard();
        const len0 = api.settings().wp_length_mm;
        api.setSetting('clearance_um', 0.6);
        const gs = api.guard();
        api.setSetting('depth_of_cut_um', 99.0);
        const warnHigh = api.guard().warnings;
        api.setSetting('clearance_um', 0);
        api.setSetting('depth_of_cut_um', api.settings().depth_of_cut_um);
        api.setSetting('wp_length_mm', 99.0);
        const warnFit = api.guard().warnings;
        const mInvalid = api.blockMatrix();
        api.setSetting('wp_length_mm', len0);
        api.reset();
        api.arm(true);
        const th = api.dragAlongArc(25.0);
        const staleAfterDrag = api.guard().stale;
        api.nudge('arc', 1, 5);
        const thNudge = api.settings().wp_position_deg;
        const pb = api.paramBlock();
        api.reset();
        // Disarm: an armed drag deliberately captures clicks on the block, so leaving
        // it on would hijack the picking sweeps that follow.
        api.arm(false);
        R.guard = {window0: g0.window, shifted: gs.window,
                   warn_high: warnHigh, warn_fit: warnFit,
                   matrix_while_invalid: mInvalid,
                   in_force: g0.in_force, resolved_deg: g0.resolved_deg,
                   arc_mm: g0.arc_mm, nudge_um: g0.nudge_um};
        R.drag = {theta: th, theta_nudged: thNudge, stale: staleAfterDrag,
                  radius_mm: (api.meta.edit || {}).ground_radius_mm,
                  param_block: pb, after_reset: api.edited()};
      }
    } catch (e) { R.guard_error = String(e); }
    save();

    // The Apply reply, in the shapes a real Colab runtime produces. Colab formats the
    // callback's return value into a mimetype bundle, so this is where a working edit
    // used to be reported as "Python refused it, no reason given".
    R.reply = null;
    try {
      const api = $('#cadwrap').cadAPI;
      if (api && api.readReply) {
        const jsonShape = {data: {'application/json': {ok: true, message: 'done'},
                                  'text/plain': 'CADREPLY {"ok": true}'}};
        const textOnly = {data: {'text/plain':
          'CADREPLY {"ok": false, "error": "depth of cut 9 um exceeds 3.8"}'}};
        const reprOnly = {data: {'text/plain': "{'ok': True, 'message': 'x'}"}};
        R.reply = {
          json: api.readReply(jsonShape),
          text: api.readReply(textOnly),
          unreadable: api.readReply(reprOnly),
          empty: api.readReply(undefined),
          str_json: api.readReply({data: {'application/json': '{"ok": true}'}})};

        // ...and end to end: a fake kernel of the real shape, through the real button.
        let sent = null;
        window.google = {colab: {kernel: {invokeFunction: async (name, args) => {
          sent = {name: name, settings: args[0]};
          return {data: {'application/json': {ok: false,
                  error: 'the block is 0.05 mm < the 0.20 mm arc'},
                  'text/plain': 'CADREPLY {"ok": false}'}};
        }}}};
        R.transport_name = api.transportName();
        api.setSetting('clearance_um', 0.3);
        api.clickApply();
        await settle(6);
        R.reply.sent = sent;
        R.reply.state_refused = api.applyState();
        R.reply.text_refused = api.applyText();

        // A reply nothing can read must not be called a refusal.
        window.google.colab.kernel.invokeFunction = async () =>
          ({data: {'text/plain': "{'ok': True}"}});
        api.clickApply();
        await settle(6);
        R.reply.state_unreadable = api.applyState();
        R.reply.text_unreadable = api.applyText();
        delete window.google;
        api.reset();
      }
    } catch (e) { R.reply_error = String(e); }
    save();

    // model tree: hide everything but the workpiece, then fit and click it
    const boxes = Array.from(document.querySelectorAll('#cadtree input'));
    await settle(4);
    const shownPixels = shot();
    boxes.forEach((b, i) => {
      if (R.tree[i] !== 'workpiece' && b.checked) { b.click(); }
    });
    await settle(4);
    R.tree_changed = shownPixels !== shot();
    R.parts_hidden = $('#cadwrap').cadAPI.parts();
    save();
    $('#cadfit').click();          // must now frame the workpiece alone
    await settle(24);
    sweep(R.isolated, 34, 24);
    boxes.forEach(b => { if (!b.checked) b.click(); });
    await settle(4);
    R.parts_after = $('#cadwrap').cadAPI.parts();
    R.done = true;
    R.raf = RAF; R.timeouts = TMO;
    save();
  } catch (e) {
    R.fatal = String(e && e.stack || e);
    save();
  }
  document.title = 'PROBE-DONE';
})();
</script>
"""


def browser() -> str:
    for b in BROWSERS:
        if os.path.exists(b):
            return b
    raise SystemExit("no Chrome or Edge found; cannot run the browser probe")


def _run_once(page_path: str, budget_ms: int) -> dict:
    import shutil
    prof = tempfile.mkdtemp(prefix="cadprobe_")
    url = "file:///" + page_path.replace("\\", "/").replace(" ", "%20")
    cmd = [browser(), "--headless=new", "--disable-gpu", "--enable-unsafe-swiftshader",
           "--no-sandbox", "--window-size=1200,800", "--user-data-dir=" + prof,
           "--virtual-time-budget=%d" % budget_ms, "--dump-dom", url]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                           encoding="utf-8", errors="replace")
    finally:
        # A browser profile is tens of MB; one per run per mode had been piling up
        # until the temp directory was the reason the probe failed.
        shutil.rmtree(prof, ignore_errors=True)
    m = re.search(r'<div id="probeout"[^>]*>(.*?)</div>', r.stdout, re.S)
    if not m:
        raise SystemExit("the probe never finished; the viewer likely threw.\n"
                         + r.stdout[-2500:])
    # The results come back through --dump-dom, which serialises the text node holding
    # them, so every '<' the page captured arrives as '&lt;'. Undo exactly that, and undo
    # '&amp;' last so a literal '&lt;' in the data survives as itself.
    txt = (m.group(1).replace("&lt;", "<").replace("&gt;", ">")
           .replace("&amp;", "&"))
    got = json.loads(txt)
    if "fatal" in got:
        raise SystemExit("the probe itself threw:\n" + got["fatal"])
    return got


def run(page_path: str) -> dict:
    """Drive the page, and insist the script reached its end.

    The page saves its results as it goes, so a run that exceeds the virtual-time budget
    dumps a *partial* result. That read as two ordinary failures in whichever checks came
    last -- picking and the automation hook -- which is a misleading way to report "the
    browser ran out of time on a loaded machine". So: retry once with a longer budget,
    and if it still stops short, say exactly that.
    """
    got = _run_once(page_path, 150_000)
    if not got.get("done"):
        print("      (probe stopped short of the end; retrying with a longer budget)")
        got = _run_once(page_path, 400_000)
    return got


FAIL = []


def chk(name, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                           (": " + detail) if detail else ""))
    if not ok:
        FAIL.append(name)


def main() -> int:
    from semgrit.analysis import AnalysisParams
    from semgrit.build_deck import DeckParams, plan_deck
    from semgrit import cadviewer

    work = os.path.join(HERE, "_cadv")
    os.makedirs(work, exist_ok=True)
    solids = pickle.load(open(os.path.join(
        HERE, "WHEEL_FIXED/1_measurements/grain_library.pkl"), "rb"))["solids"]
    p = DeckParams(name="probe", sector_mode="arc", arc_length_mm=1.0,
                   grit_mode="count", grit_count=60,
                   analysis=AnalysisParams(enabled=True))
    plan = plan_deck(p, solids)

    for mode in ("contact", "wheel", "whole wheel"):
        print("\n%s mode" % mode.upper())
        html, meta, info = cadviewer.build(
            plan, os.path.join(work, "%s.glb" % mode), mode=mode,
            max_grits=10 ** 6, height=760)
        page = os.path.join(work, "probe_%s.html" % mode)
        with open(page, "w", encoding="utf-8", newline="\n") as fh:
            fh.write('<!doctype html><meta charset="utf-8">'
                     '<body style="margin:0">' + html + PROBE + "</body>")
        got = run(page)

        chk("the probe script ran every stage to the end", bool(got.get("done")),
            "last stage %s" % (got.get("stages") or "none"))
        chk("three.js loads and the model parses in a real browser",
            "rendered locally" in got.get("status", ""), got.get("status", "")[:80])
        want_tree = ["bond rim", "abrasive grits", "workpiece"]
        if mode == "whole wheel":
            want_tree = ["bond rim", "whole wheel (context)",
                         "contact marker (pointer)", "abrasive grits", "workpiece"]
        chk("the model tree lists every part",
            got.get("tree", []) == want_tree, str(got.get("tree", [])))
        chk("the Face view direction is a unit vector along the contact normal",
            abs(sum(x * x for x in meta.get("face_dir", [0])) - 1) < 1e-9,
            str([round(x, 4) for x in meta.get("face_dir", [])]))

        grains = {}
        for txt in got.get("picks", {}):
            m = re.match(r"grain (\d+)", txt)
            if m:
                grains[int(m.group(1))] = txt
        # Three distinct grains is enough to prove the triangle-range mapping works
        # across grain boundaries; the ceiling here is how often a coarse click grid
        # lands on a few-micron grain, not the viewer.
        chk("clicking picks out individual grains", len(grains) >= 3,
            "%d distinct grains identified by real pointer events" % len(grains))
        ids = {g["id"] for g in meta["grains"]}
        chk("every picked id is a real grain of this deck",
            bool(grains) and set(grains) <= ids, "picked %s" % sorted(grains)[:10])

        bad = []
        for gid, txt in grains.items():
            rec = next(g for g in meta["grains"] if g["id"] == gid)
            for label, k in (("protrusion", "protrusion_um"), ("height", "height_um"),
                             ("width", "width_um"), ("volume", "volume_um3")):
                m = re.search(label + r" ([-\d.]+)", txt)
                # The panel prints three decimals -- nanometre resolution on a micron
                # -- so compare at what is shown, not at the stored precision.
                if not m or abs(float(m.group(1)) - rec[k]) > 5e-4:
                    bad.append("grain %s %s: %s vs %s"
                               % (gid, label, m and m.group(1), rec[k]))
        chk("the numbers on screen are the deck's own numbers", not bad, str(bad[:3]))
        chk("clicking the bond names it and gives its radius",
            any(t.startswith("bond rim") for t in got.get("picks", {})),
            "%d distinct results" % len(got.get("picks", {})))
        # Context parts may legitimately be hidden by the zoom rule, so compare the
        # names and the context flags, not the visibility.
        chk("the automation hook reports every part and flags the context",
            isinstance(got.get("parts_after"), list)
            and [p["name"] for p in got["parts_after"]] == want_tree
            and sum(1 for p in got["parts_after"] if p["context"])
            == (2 if mode == "whole wheel" else 0),
            str([(p["name"], p["context"]) for p in got.get("parts_after", [])]))
        chk("hiding parts in the tree isolates the workpiece",
            any(t.startswith("workpiece") for t in got.get("isolated", {}))
            and not any(t.startswith("grain ") for t in got.get("isolated", {})),
            str(list(got.get("isolated", {}))[:3]))
        print("      stages: %s  raf/timeout %s/%s  js errors: %s"
              % (got.get("stages"), got.get("raf"), got.get("timeouts"),
                 got.get("errors") or "none"))
        chk("the tree toggle changes what is drawn", bool(got.get("tree_changed")))
        chk("the section plane changes what is drawn", bool(got.get("section_changed")))
        chk("the edges toggle changes what is drawn", bool(got.get("edges_changed")))
        chk("shift-click measures a distance",
            got.get("measurement", "").startswith("distance"),
            got.get("measurement", "").encode("ascii", "replace").decode())
        comp = got.get("components") or {}
        want = ["dx_mm", "dy_mm", "dz_mm", "radial_mm", "along_arc_mm",
                "across_face_mm", "total_mm"]
        chk("it reports components as well as the total",
            all(k in comp for k in want), str(sorted(comp)))
        panel = got.get("panel", "")
        chk("and the panel shows them",
            all(t in panel for t in ("X", "Y", "Z", "radial", "along arc",
                                     "across face")),
            panel.encode("ascii", "replace").decode()[:100])
        if all(k in comp for k in want):
            tot = comp["total_mm"]
            g = (comp["dx_mm"] ** 2 + comp["dy_mm"] ** 2 + comp["dz_mm"] ** 2) ** 0.5
            f = (comp["radial_mm"] ** 2 + comp["along_arc_mm"] ** 2
                 + comp["across_face_mm"] ** 2) ** 0.5
            # Both triplets are orthonormal decompositions of the same vector, so
            # either failing to reproduce the length means the frame is wrong.
            chk("the components reproduce the distance in both frames",
                abs(g - tot) < 1e-9 and abs(f - tot) < 1e-9,
                "global %.9f, contact %.9f, total %.9f" % (g, f, tot))
        g = got.get("guard")
        if g:
            chk("the viewer knows the depth-of-cut window",
                g["window0"]["hi"] > g["window0"]["lo"],
                "%.4f .. %.4f um" % (g["window0"]["lo"], g["window0"]["hi"]))
            chk("a standoff shifts both ends of the window by exactly itself",
                abs((g["shifted"]["lo"] - g["window0"]["lo"]) - 0.6) < 1e-9
                and abs((g["shifted"]["hi"] - g["window0"]["hi"]) - 0.6) < 1e-9,
                "lo %+.6f hi %+.6f" % (g["shifted"]["lo"] - g["window0"]["lo"],
                                       g["shifted"]["hi"] - g["window0"]["hi"]))
            chk("too deep is warned about before Apply",
                any("bond rim into the work" in w for w in g["warn_high"]),
                (g["warn_high"] or ["none"])[0][:62])
            chk("a block that would not fit is warned about before Apply",
                any("hang off both ends" in w for w in g["warn_fit"]),
                (g["warn_fit"] or ["none"])[0][:62])
            chk("a state that would be refused is not previewed",
                g["matrix_while_invalid"][0] == 1
                and g["matrix_while_invalid"][5] == 1
                and g["matrix_while_invalid"][10] == 1,
                "block transform stayed identity")
            chk("the wheel extent is reported resolved, not raw",
                g["in_force"] in ("arc_length_mm", "sector_deg")
                and g["resolved_deg"] > 0,
                "%s in force, %.4f deg, %.4f mm arc"
                % (g["in_force"], g["resolved_deg"], g["arc_mm"]))
        d = got.get("drag")
        if d:
            want = math.degrees(0.025 / d["radius_mm"])
            chk("a 25 um drag becomes the angle that displacement subtends",
                abs(d["theta"] - want) < 1e-5,
                "browser %.6f vs python %.6f deg" % (d["theta"], want))
            step = math.degrees(0.0005 / d["radius_mm"])
            chk("five 0.1 um nudges add exactly that much again",
                abs((d["theta_nudged"] - d["theta"]) - step) < 1e-5,
                "%.6f vs %.6f deg" % (d["theta_nudged"] - d["theta"], step))
            chk("dragging marks the seating stale rather than guessing it",
                bool(d["stale"]))
            chk("dragging offers the widget lines to paste back",
                "WP_POSITION_DEG" in d["param_block"],
                d["param_block"].replace("\n", " . ")[:66])
            chk("reset clears every edit", d["after_reset"] == {})
        if got.get("guard_error"):
            chk("the guardrail drive did not throw", False, got["guard_error"][:70])
        rp = got.get("reply")
        if rp:
            chk("a Colab application/json reply is read as the applied result",
                rp["json"].get("ok") is True and rp["json"].get("message") == "done")
            chk("and one delivered only as text is read exactly, not guessed at",
                rp["text"].get("ok") is False
                and "exceeds" in (rp["text"].get("error") or ""),
                (rp["text"].get("error") or "")[:52])
            chk("a JSON string under that mimetype is parsed too",
                rp["str_json"].get("ok") is True)
            # The old code turned anything it did not recognise into {ok:false} with no
            # message, which is how a successful Apply came back as a refusal.
            chk("a reply nothing can parse is flagged unreadable, not refused",
                rp["unreadable"].get("unreadable") is True
                and "ok" not in rp["unreadable"],
                str(rp["unreadable"])[:60])
            chk("and it keeps the raw text so the reason is not lost",
                "ok" in (rp["unreadable"].get("raw") or ""),
                (rp["unreadable"].get("raw") or "")[:52])
            chk("no reply at all is unreadable rather than a refusal",
                rp["empty"].get("unreadable") is True)
            chk("the viewer finds a Colab kernel and calls cad.commit with the settings",
                got.get("transport_name") == "colab"
                and (rp.get("sent") or {}).get("name") == "cad.commit"
                and abs(((rp.get("sent") or {}).get("settings")
                         or {}).get("clearance_um", -1) - 0.3) < 1e-12,
                "%s -> %s" % (got.get("transport_name"),
                              (rp.get("sent") or {}).get("name")))
            chk("a real refusal reaches the panel with its reason intact",
                rp["state_refused"] == "failed"
                and "refused" in rp["text_refused"]
                and "0.05 mm < the 0.20 mm arc" in rp["text_refused"],
                rp["text_refused"][:72])
            chk("an unreadable reply does not accuse Python of refusing",
                rp["state_unreadable"] == "unknown"
                and "refused" not in rp["text_unreadable"]
                and "could not be read" in rp["text_unreadable"],
                rp["text_unreadable"][:72])
        if got.get("reply_error"):
            chk("the reply drive did not throw", False, got["reply_error"][:70])
        chk("Wheel and Contact frame different things",
            bool(got.get("wheel_vs_contact_differ"))
            and not got.get("wheel_blank") and not got.get("contact_blank"),
            "wheel blank=%s contact blank=%s"
            % (got.get("wheel_blank"), got.get("contact_blank")))

    print("\n" + "=" * 74)
    print("BROWSER PROBE: %d failure(s)%s"
          % (len(FAIL), "" if not FAIL else " -> " + str(sorted(set(FAIL)))))
    print("=" * 74)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
