"""A CAD viewer, built on three.js, that runs inside the notebook with no account.

``<model-viewer>`` renders the model beautifully but it is a product viewer: it has no
section plane, no edges, no model tree and no way to interrogate a part. Those are the
things that make SolidWorks or the Autodesk viewer feel like an engineering tool rather
than a turntable, and they are all client-side features -- none of them is the reason
those products charge.

So this builds them directly:

* **shaded with edges** -- the single strongest "this is CAD" cue. Feature edges above
  a threshold angle, drawn over a physically-shaded surface.
* **section plane** on any axis, with a slider, so you can cut through the contact and
  see how deep each grain sits in the bond.
* **standard views** -- iso, front, top, right -- plus zoom-to-fit, and an
  orthographic/perspective toggle, because a CAD user expects ortho.
* **model tree** with per-part visibility and isolate.
* **click a grain** to read its id, protrusion, height, width and volume. The grits are
  one merged mesh for speed, so picking maps the triangle index back to a grain through
  the ranges recorded when the mesh was built.
* studio lighting, contact shadow, gradient background, screenshot.

Everything is a single self-contained HTML string with the model embedded as a data
URI. Nothing is uploaded and no key is needed.
"""

from __future__ import annotations

import base64
import json
import os
from typing import Optional

THREE = "0.160.0"


def viewer_html(glb_path: str, meta: dict, height: int = 720,
                max_inline_mb: float = 24.0) -> str:
    """Self-contained HTML for the viewer, with the .glb inlined."""
    mb = os.path.getsize(glb_path) / 1e6
    if mb > max_inline_mb:
        raise ValueError(
            "%.1f MB .glb is too large to inline (cap %.1f MB). Raise "
            "max_inline_mb, use mode='contact', or open the .glb in Blender."
            % (mb, max_inline_mb))
    b64 = base64.b64encode(open(glb_path, "rb").read()).decode()
    return _TEMPLATE % {
        "three": THREE,
        "h": height,
        "b64": b64,
        "meta": json.dumps(meta),
        "mb": mb,
    }


# The template is one string so the whole viewer travels as a single artefact. %(x)s
# placeholders are filled above; every literal percent in the CSS/JS is doubled.
_TEMPLATE = r"""
<div id="cadwrap" style="position:relative;width:100%%;height:%(h)dpx;
     border:1px solid #b9c0c8;border-radius:6px;overflow:hidden;
     background:radial-gradient(120%% 100%% at 50%% 0%%,#fdfefe 0%%,#e8eef4 45%%,#b9c6d4 100%%);
     font-family:system-ui,'Segoe UI',sans-serif;color:#22282e">
  <canvas id="cadcanvas" style="display:block;width:100%%;height:100%%"></canvas>

  <div id="cadtools" style="position:absolute;top:10px;left:10px;background:#ffffffe8;
       border:1px solid #b9c0c8;border-radius:6px;padding:8px 10px;font-size:12px;
       line-height:1.65;width:206px;max-height:calc(100%% - 46px);overflow-y:auto;
       overflow-x:hidden;overscroll-behavior:contain;
       box-shadow:0 2px 8px #0002;backdrop-filter:blur(4px)">
    <div class="cadhd">View</div>
    <div style="display:flex;gap:3px;flex-wrap:wrap">
      <button data-v="iso">Iso</button><button data-v="front">Front</button>
      <button data-v="top">Top</button><button data-v="right">Right</button>
      <button data-v="face" title="straight at the dressed face">Face</button>
      <button data-v="axial" title="along the wheel axis">Axial</button>
      <button data-v="wheelview" title="the whole wheel, where the curvature shows"
              >Wheel</button>
      <button data-v="contact" title="zoom to the grains on the workpiece"
              >Contact</button>
      <button id="cadfit" title="zoom to fit (F)">Fit</button>
      <button id="caddrag" title="drag the workpiece (G). Shift-drag for standoff, arrow keys nudge">Drag block</button>
    </div>
    <label><input type="checkbox" id="cadedges" checked> shaded with edges</label>
    <label><input type="checkbox" id="cadortho"> orthographic</label>
    <label><input type="checkbox" id="cadspin"> spin</label>

    <div class="cadhd">Section plane</div>
    <select id="cadaxis" style="width:100%%">
      <option value="-1">off</option>
      <option value="0">normal to X</option>
      <option value="1">normal to Z &nbsp;(wheel axis)</option>
      <option value="2">normal to Y</option>
    </select>
    <input type="range" id="cadcut" min="0" max="1000" value="500"
           style="width:100%%;margin:4px 0 0">
    <label><input type="checkbox" id="cadflip"> flip side</label>
    <label><input type="checkbox" id="cadcap" checked> cap the cut face</label>

    <div class="cadhd">Explode</div>
    <input type="range" id="cadexplode" min="0" max="100" value="0"
           style="width:100%%" title="pull the parts apart along the radius">
    <div id="cadexplodetxt" style="font-size:10px;color:#5a636c">assembled</div>

    <div class="cadhd">Colour the grains by</div>
    <select id="cadcolour" style="width:100%%">
      <option value="none">part colour</option>
      <option value="protrusion_um">protrusion above the bond</option>
      <option value="height_um">grain height</option>
      <option value="width_um">grain width</option>
      <option value="volume_um3">grain volume</option>
      <option value="engage">engages the block (yes / no)</option>
    </select>
    <div id="cadlegend" style="margin-top:4px;display:none">
      <div id="cadlegbar" style="height:10px;border:1px solid #a9b2bb;
           border-radius:2px"></div>
      <div style="display:flex;justify-content:space-between;font-size:10px;
           color:#5a636c"><span id="cadleglo"></span><span id="cadleghi"></span></div>
    </div>

    <div class="cadhd">Measure</div>
    <div style="color:#5a636c">shift-click two points</div>
    <button id="cadclear">clear</button>

    <div class="cadhd">Parts</div>
    <div id="cadtree"></div>

    <div class="cadhd">Boundary conditions</div>
    <div id="cadbclist" style="font-size:11px"></div>

    <div class="cadhd">Edit and rebuild</div>
    <div id="cadeditlist" style="font-size:11px"></div>
    <button id="cadshot" style="margin-top:8px">Save PNG</button>
  </div>

  <div id="cadinfo" style="position:absolute;top:10px;right:10px;background:#ffffffe8;
       border:1px solid #b9c0c8;border-radius:6px;padding:8px 10px;font-size:12px;
       line-height:1.6;width:224px;max-height:calc(100%% - 66px);overflow-y:auto;
       overflow-x:hidden;overscroll-behavior:contain;
       box-shadow:0 2px 8px #0002;backdrop-filter:blur(4px)">
    <div class="cadhd" style="margin-top:0">Inspect</div>
    <div id="cadpick" style="color:#5a636c">click a grain</div>
  </div>

  <div id="cadstat" style="position:absolute;bottom:9px;left:11px;right:210px;
       font-size:11px;color:#39424b;white-space:nowrap;overflow:hidden;
       text-overflow:ellipsis;
       font-family:ui-monospace,Consolas,monospace">loading three.js &hellip;</div>
  <div id="cadhelp" style="position:absolute;bottom:9px;right:64px;font-size:11px;
       color:#5a636c;max-width:44%%;text-align:right;white-space:nowrap;
       overflow:hidden;text-overflow:ellipsis"
       >drag orbit &middot; right-drag pan &middot; scroll zoom</div>

  <!-- Viewport size. A notebook output cell is a fixed box, and 720 px is a
       compromise between "big enough to work in" and "does not push the next cell
       off the screen". Rather than pick a bigger compromise, let the user choose:
       fullscreen for inspection, a drag handle for a permanent nudge. -->
  <div style="position:absolute;bottom:26px;right:11px;z-index:5;display:flex;gap:4px">
    <button id="cadkeysbtn" title="keyboard shortcuts (?)">?</button>
    <button id="cadfull"
            title="fullscreen (double-click the canvas, or Esc to leave)"
            >&#x26F6;</button>
  </div>
  <div id="cadgrip" title="drag to resize the viewer"
       style="position:absolute;left:0;right:0;bottom:0;height:7px;cursor:ns-resize;
              background:linear-gradient(#0000,#0000001a)"></div>

  <!-- First run: the viewer opens with no indication that anything is clickable.
       Three lines, dismissed for good on the first click. -->
  <div id="cadfirst" style="position:absolute;left:50%%;top:50%%;
       transform:translate(-50%%,-50%%);background:#12365ee8;color:#fff;
       padding:13px 17px;border-radius:8px;font-size:12.5px;line-height:1.75;
       box-shadow:0 6px 24px #0005;z-index:6;pointer-events:none;text-align:left">
    <b style="font-size:13.5px">This is the model the .inp contains</b><br>
    <span style="opacity:.92">
    &#x25CF;&nbsp; <b>click a grain</b> to read its size and protrusion<br>
    &#x25CF;&nbsp; <b>G</b> to drag the workpiece, <b>F</b> to fit, <b>?</b> for all keys<br>
    &#x25CF;&nbsp; <b>section plane</b> and <b>colour by</b> are on the left</span>
    <div style="opacity:.66;margin-top:5px;font-size:11px">click anywhere to dismiss</div>
  </div>

  <!-- "You are here". At wheel scale the deck is a 2 mm slice on a 50 mm disc and
       the contact is far below one pixel, so the opening view is an empty grey
       circle unless something points at it. The 3-D marker cannot be made bigger
       without becoming larger than what it points at, so the label lives in screen
       space instead: always legible, never part of the geometry. -->
  <div id="cadwhere" title="jump to the contact (C)"
       style="position:absolute;left:0;top:0;z-index:4;display:none;cursor:pointer;
       transform:translate(-50%%,-50%%)">
    <div style="position:relative">
      <div style="width:34px;height:34px;border:2px solid #f2731a;border-radius:50%%;
                  box-shadow:0 0 0 2px #ffffffcc,0 0 9px #f2731a55"></div>
      <div id="cadwherelab" style="position:absolute;left:44px;top:3px;
           white-space:nowrap;background:#f2731aE8;color:#fff;font-size:11px;
           font-weight:600;padding:2px 7px;border-radius:4px;
           box-shadow:0 1px 5px #0003">the cut &middot; press C</div>
    </div>
  </div>

  <!-- Keyboard reference, toggled with ? -->
  <div id="cadkeys" style="position:absolute;left:50%%;top:50%%;
       transform:translate(-50%%,-50%%);background:#ffffffF2;border:1px solid #b9c0c8;
       border-radius:8px;padding:12px 16px;font-size:12px;line-height:1.8;
       box-shadow:0 6px 24px #0004;z-index:7;display:none">
    <b>Keyboard</b><br>
    <table style="font-size:11.5px;border-spacing:9px 1px">
      <tr><td><kbd>F</kbd></td><td>zoom to fit</td>
          <td><kbd>G</kbd></td><td>arm / disarm block drag</td></tr>
      <tr><td><kbd>&larr;</kbd><kbd>&rarr;</kbd></td><td>nudge along the arc</td>
          <td><kbd>&uarr;</kbd><kbd>&darr;</kbd></td><td>nudge the standoff</td></tr>
      <tr><td><kbd>1</kbd>..<kbd>4</kbd></td><td>iso / front / top / right</td>
          <td><kbd>5</kbd> <kbd>6</kbd></td><td>face / axial</td></tr>
      <tr><td><kbd>W</kbd></td><td>wheel view</td>
          <td><kbd>C</kbd></td><td>contact view</td></tr>
      <tr><td><kbd>E</kbd></td><td>edges on / off</td>
          <td><kbd>O</kbd></td><td>orthographic</td></tr>
      <tr><td><kbd>X</kbd></td><td>section: cycle axis</td>
          <td><kbd>shift</kbd>+click</td><td>measure two points</td></tr>
      <tr><td><kbd>Esc</kbd></td><td>cancel a drag</td>
          <td><kbd>?</kbd></td><td>this list</td></tr>
    </table>
  </div>

  <style>
    #cadwrap button{font:inherit;font-size:11px;padding:2px 7px;cursor:pointer;
      border:1px solid #a9b2bb;border-radius:4px;background:
      linear-gradient(#ffffff,#e9eef3)}
    #cadwrap button:hover{background:linear-gradient(#ffffff,#d7e3ee)}
    #cadwrap button.on{background:linear-gradient(#cfe4f7,#a9cdec);border-color:#5b8fbf}
    #cadwrap label{display:block;cursor:pointer;user-select:none}
    #cadwrap .cadhd{font-weight:600;margin:9px 0 3px;padding-top:6px;
      border-top:1px solid #dfe4e9;cursor:pointer;user-select:none;
      display:flex;align-items:center;gap:4px}
    #cadwrap .cadhd:first-child{margin-top:0;padding-top:0;border-top:none}
    /* The left column holds ten sections in 198 px. Collapsing is what makes it
       navigable rather than a scroll well; the caret is the affordance. */
    #cadwrap .cadhd .car{display:inline-block;width:8px;font-size:9px;color:#7a838c;
      transition:transform .12s}
    #cadwrap .cadhd.shut .car{transform:rotate(-90deg)}
    #cadwrap .cadsec.shut{display:none}
    #cadwrap kbd{background:#eef1f4;border:1px solid #c3cad1;border-bottom-width:2px;
      border-radius:3px;padding:0 4px;font:inherit;font-size:10.5px;
      font-family:ui-monospace,Consolas,monospace}
    #cadwrap:fullscreen{border-radius:0;height:100%%!important}
    #cadwrap:fullscreen #cadgrip{display:none}
    #cadtree label{white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
    #cadtree .sw{display:inline-block;width:9px;height:9px;border:1px solid #7a838c;
      border-radius:2px;margin-right:4px;vertical-align:-1px}
  </style>
</div>

<script type="importmap">
{"imports":{"three":"https://unpkg.com/three@%(three)s/build/three.module.js",
 "three/addons/":"https://unpkg.com/three@%(three)s/examples/jsm/"}}
</script>

<script type="module">
const STAT = document.getElementById('cadstat');
const PICK = document.getElementById('cadpick');
try {
  const THREE = await import('three');
  const { OrbitControls } = await import('three/addons/controls/OrbitControls.js');
  const { GLTFLoader }    = await import('three/addons/loaders/GLTFLoader.js');
  const { RoomEnvironment } =
        await import('three/addons/environments/RoomEnvironment.js');

  const META = %(meta)s;
  const wrap = document.getElementById('cadwrap');
  const canvas = document.getElementById('cadcanvas');

  // stencil: true is required by the section cap, which fills the cut face using a
  // stencil pass. WebGL gives a stencil buffer by default, but asking is explicit
  // and costs nothing.
  const renderer = new THREE.WebGLRenderer({canvas, antialias:true, alpha:true,
                                            stencil:true,
                                            preserveDrawingBuffer:true});
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.localClippingEnabled = true;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.02;

  const scene = new THREE.Scene();
  const pmrem = new THREE.PMREMGenerator(renderer);
  scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
  scene.add(new THREE.HemisphereLight(0xffffff, 0x8e97a3, 1.25));
  const keyL = new THREE.DirectionalLight(0xffffff, 1.9); keyL.position.set(1, 1.5, 1);
  const fill = new THREE.DirectionalLight(0xffffff, 0.6); fill.position.set(-1, 0.3, -1);
  const rim = new THREE.DirectionalLight(0xffffff, 0.45); rim.position.set(0, -1, 0.4);
  scene.add(keyL, fill, rim);

  const persp = new THREE.PerspectiveCamera(32, 1, 1e-4, 1e5);
  const ortho = new THREE.OrthographicCamera(-1, 1, 1, -1, -1e4, 1e5);
  let cam = persp;
  const controls = new OrbitControls(cam, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.12;

  // ---- section plane ------------------------------------------------------
  // A plane keeps the half-space where dot(n,p)+constant > 0, so "off" is just a
  // constant large enough that nothing is ever behind it.
  const clip = new THREE.Plane(new THREE.Vector3(0, -1, 0), 1e9);
  const axisSel  = document.getElementById('cadaxis');
  const cutSlide = document.getElementById('cadcut');
  const flipBox  = document.getElementById('cadflip');

  const parts = [];                    // {mesh, edges, name, pretty}
  let root = null, box = null, deckBox = null, sphere = null, gritMesh = null;
  let hilite = null;

  // GLTFLoader sanitises node names (spaces become underscores), so keep the
  // pretty label separately and match on the sanitised form. Getting this wrong is
  // silent: picking would simply never recognise a grain.
  // GLTFLoader sanitises node names (spaces become underscores and dots go), so keep
  // the pretty label separately and match on the sanitised form. Getting this wrong is
  // silent: picking would simply never recognise a grain.
  const sane = s => String(s || 'part').replace(/[\s.:]/g, '_');
  // Python's messages compare numbers, so they contain '<' and '>'. Inserted raw, an
  // error like "0.05 mm < the 0.20 mm arc" swallows the rest of itself as a tag.
  const esc = s => String(s === undefined || s === null ? '' : s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const FAR_KEY = sane(META.far_part || 'abrasive grits (far, simplified)');
  const GHOST_KEY = sane(META.ghost_part || 'whole wheel (context, not in the deck)');
  const MARK_KEY = sane(META.mark_part || 'contact marker (not in the deck)');
  const PRETTY = {bond_rim:'bond rim', abrasive_grits:'abrasive grits',
                  workpiece:'workpiece'};
  PRETTY[FAR_KEY] = 'grits far off (boxes)';
  PRETTY[GHOST_KEY] = 'whole wheel (context)';
  PRETTY[MARK_KEY] = 'contact marker (pointer)';
  const SWATCH = {bond_rim:'#b8bcc4', abrasive_grits:'#29a05c',
                  workpiece:'#2e70b5'};
  SWATCH[FAR_KEY] = '#6b9478';
  SWATCH[GHOST_KEY] = '#8b939b';
  SWATCH[MARK_KEY] = '#f2731a';

  // Feature edges are what make this read as CAD, but EdgesGeometry over hundreds of
  // thousands of triangles is slow to build and draws a wire ball nobody can see
  // through. Above this many triangles the grit meshes start with edges off; the bond
  // and the workpiece -- which carry the CAD look -- always get them.
  const EDGE_TRI_BUDGET = 120000;
  let edgesSuppressed = 0;
  let farMesh = null, ghostMesh = null;

  new GLTFLoader().load('data:model/gltf-binary;base64,%(b64)s', (gltf) => {
    root = gltf.scene;
    scene.add(root);
    box = new THREE.Box3().setFromObject(root);
    sphere = box.getBoundingSphere(new THREE.Sphere());

    root.traverse(o => {
      if (!o.isMesh) return;
      const key = sane(o.name);
      const ntri = o.geometry.index ? o.geometry.index.count / 3
                                    : o.geometry.attributes.position.count / 3;
      const isGrit = key === 'abrasive_grits' || key === FAR_KEY;
      const isGhost = key === GHOST_KEY || key === MARK_KEY;
      o.material.side = THREE.DoubleSide;      // a cut face must not look hollow
      o.material.clippingPlanes = [clip];
      o.material.polygonOffset = true;         // keep edges off the surface
      o.material.polygonOffsetFactor = 1;
      o.material.polygonOffsetUnits = 1;
      o.material.needsUpdate = true;
      // The ghost is context, not the deck: never let it swallow a click, and keep it
      // behind everything else.
      // Context, not the deck: never let it swallow a click, and do not let a
      // translucent shell write depth over the geometry it is meant to sit around.
      if (isGhost) {
        o.raycast = () => {};
        o.renderOrder = -1;
        o.material.depthWrite = false;
      }

      let edges = null;
      const heavy = isGrit && ntri > EDGE_TRI_BUDGET;
      if (heavy) { edgesSuppressed += ntri; }
      // 24 deg threshold: silhouettes and the bond's structure, without every facet.
      const em = new THREE.LineBasicMaterial({color:0x243039, transparent:true,
                                              opacity:isGhost ? 0.25 : 0.6});
      em.clippingPlanes = [clip];
      // An empty placeholder must NOT go through EdgesGeometry: its constructor reads
      // geometry.attributes.position.count, so on an empty BufferGeometry it throws
      // "Cannot read properties of undefined (reading 'count')" -- and because
      // GLTFLoader routes exceptions from onLoad into onError, the whole model then
      // reports as unloadable. An empty geometry is fine to *render*, just not to
      // derive edges from.
      edges = new THREE.LineSegments(
        heavy ? new THREE.BufferGeometry() : new THREE.EdgesGeometry(o.geometry, 24),
        em);
      edges.renderOrder = 1;
      edges.userData.heavy = heavy;
      edges.userData.owner = o;
      o.add(edges);

      if (key === 'abrasive_grits') gritMesh = o;
      if (key === FAR_KEY) farMesh = o;
      if (isGhost) ghostMesh = o;
      parts.push({mesh:o, edges, key, isGhost, wanted:true,
                  pretty:PRETTY[key] || o.name, ntri:ntri});
    });

    const tree = document.getElementById('cadtree');
    parts.forEach((p, i) => {
      const lb = document.createElement('label');
      lb.innerHTML = '<input type="checkbox" checked> <span class="sw" style="' +
                     'background:' + (SWATCH[p.key] || '#999') + '"></span>' + p.pretty;
      lb.title = Math.round(p.ntri).toLocaleString() + ' triangles' +
                 (p.isGhost ? ' - context only, not part of the deck' : '');
      lb.querySelector('input').onchange = e => {
        p.wanted = e.target.checked;
        p.mesh.visible = p.wanted;      // the loop re-applies the context rule
      };
      tree.appendChild(lb);
    });

    buildBC();

    // NOT inside buildBC(): that returns early when the view has no
    // boundary conditions, and it also used to be nested in the engage
    // block. Either way a view with no per-grain data -- a mesh view, or
    // a deck whose grains all miss the work -- never reached it, so the
    // colour-by selector stayed ENABLED but wired to nothing: picking an
    // entry did nothing and explained nothing. setupColourBy's own guard
    // disables it and says why, so it has to be reached to do that.
    setupColourBy();
    // A mesh view is fed element geometry, not grain solids, so the panel's
    // grain wording is wrong there. Retitled here rather than forking the
    // template: one template, two vocabularies, and the labels stay honest
    // about what a click will actually tell you.
    if (META.kind === 'mesh') {
      const hd = [...document.querySelectorAll('#cadtools .cadhd')]
        .find(h => /colour the grains by/i.test(h.textContent));
      if (hd) hd.textContent = 'Colour by';
      const pick = document.getElementById('cadpick');
      if (pick) pick.textContent = 'click an element face';
      const first = document.getElementById('cadfirst');
      if (first) {
        const b = first.querySelector('b');
        if (b) b.textContent = 'This is the MESH the .inp contains';
        first.innerHTML = first.innerHTML
          .replace(/<b>click a grain<\/b> to read its size and protrusion/,
                   '<b>element edges</b> are drawn: check that dc is resolved')
          .replace(/<b>colour by<\/b>/, '<b>explode</b>');
      }
      const insp = [...document.querySelectorAll('#cadinfo .cadhd')]
        .find(h => /inspect/i.test(h.textContent));
      if (insp) insp.textContent = 'Mesh';
    }
    buildEditPanel();
    makeCollapsible();

    // The section plane must range over the deck, not over the ghost: a slider that
    // spans 50 mm of context moves in steps far bigger than the whole model and never
    // appears to cut anything.
    deckBox = new THREE.Box3();
    parts.forEach(p => { if (!p.isGhost) deckBox.expandByObject(p.mesh); });
    if (deckBox.isEmpty()) deckBox = box.clone();

    // In whole-wheel mode open on the wheel: fitting the deck would frame a 2 mm
    // slice and leave the 50 mm disc off screen, which is the opposite of the point.
    if (META.ghost) {
      frameAt(new THREE.Vector3(0, 0, 0), META.wheel_radius_mm || sphere.radius,
              AXIS.clone().addScaledVector(FACE, 0.25));
    } else {
      frame(ISO);
    }
    applyClip();
    // Report the deck's own extent. The ghost and the marker are pointers, and
    // letting a 50 mm pointer set the number would misdescribe the model.
    const d = deckBox.getSize(new THREE.Vector3());
    const far = (META.grains_far || []).length;
    // glTF is Y-up, the model is Z-up: report in the model's own frame.
    STAT.textContent = 'extent  X ' + d.x.toFixed(3) + '  Y ' + d.z.toFixed(3) +
      '  Z ' + d.y.toFixed(3) + ' mm   |   ' + META.grits_drawn + ' of ' +
      META.grits_total + ' grits' + (far ? ' (' + far + ' as boxes)' : '') +
      '   |   %(mb).2f MB, rendered locally' +
      (edgesSuppressed ? '   |   edges off on the grits (' +
       Math.round(edgesSuppressed).toLocaleString() + ' triangles) - tick to force'
       : '');
  }, undefined, e => { STAT.textContent = 'could not load the model: ' + e; });

  // ---- camera -------------------------------------------------------------
  function visibleSphere() {
    // Fit what is on screen, not what was loaded: with the bond hidden, a fit that
    // still framed the bond would leave the workpiece a speck in the middle.
    const b = new THREE.Box3();
    let any = false;
    parts.forEach(p => {
      if (p.mesh.visible && !p.isGhost) { b.expandByObject(p.mesh); any = true; }
    });
    // ...unless the deck itself is hidden, in which case fit whatever is left
    if (!any) {
      parts.forEach(p => { if (p.mesh.visible) { b.expandByObject(p.mesh); any = true; } });
    }
    return (any ? b : box).getBoundingSphere(new THREE.Sphere());
  }
  function frame(dir) {
    if (!box) return;
    const s0 = visibleSphere();
    const r = s0.radius || 1, c = s0.center;
    const dist = 1.08 * r / Math.sin(THREE.MathUtils.degToRad(persp.fov * 0.5));
    persp.near = Math.max(r / 5000, 1e-6); persp.far = dist * 40;
    const s = r * 1.08;
    ortho.top = s; ortho.bottom = -s; ortho.near = -dist * 20; ortho.far = dist * 40;
    cam.position.copy(c).add(dir.clone().normalize().multiplyScalar(dist));
    controls.target.copy(c);
    resize();                       // recomputes the ortho frustum from the aspect
    controls.update();
  }
  const ISO = new THREE.Vector3(1, 0.72, 1);
  // Face and Axial are in the wheel's own frame: looking back along the contact
  // normal shows the dressed surface with every grain in view, which is the view a
  // grinding paper prints. The global Front/Top/Right are kept for orientation.
  const FACE = new THREE.Vector3().fromArray(META.face_dir || [1, 0, 0]);
  const AXIS = new THREE.Vector3().fromArray(META.axis_dir || [0, 1, 0]);
  const ARC = new THREE.Vector3().fromArray(META.arc_dir || [0, 0, 1]);
  const CONTACT = new THREE.Vector3().fromArray(META.contact_centre || [0, 0, 0]);

  // Zoom to a chosen point and radius rather than to the whole model. This is what
  // makes one viewer serve both scales: the whole 50 mm wheel, where the curvature is
  // the thing to see, and the few-micron contact, where the grains are.
  function frameAt(centre, radius, dir) {
    const r = Math.max(radius, 1e-5);
    const dist = 1.15 * r / Math.sin(THREE.MathUtils.degToRad(persp.fov * 0.5));
    persp.near = Math.max(r / 5000, 1e-6); persp.far = dist * 200 + sphere.radius * 4;
    ortho.top = r * 1.15; ortho.bottom = -r * 1.15;
    ortho.near = -dist * 20 - sphere.radius * 4; ortho.far = dist * 40 +
                                                            sphere.radius * 4;
    cam.position.copy(centre).add(dir.clone().normalize().multiplyScalar(dist));
    controls.target.copy(centre);
    resize();
    controls.update();
  }
  document.querySelectorAll('#cadtools button[data-v]').forEach(b => {
    b.onclick = () => {
      if (b.dataset.v === 'contact') {
        frameAt(CONTACT, META.contact_radius_mm || 0.05,
                FACE.clone().addScaledVector(AXIS, 0.35).addScaledVector(ARC, 0.35));
        return;
      }
      if (b.dataset.v === 'wheelview') {
        // Straight down the axis at the whole disc: the one view where a 50 mm wheel
        // reads as a wheel.
        frameAt(new THREE.Vector3(0, 0, 0), META.wheel_radius_mm || sphere.radius,
                AXIS.clone().addScaledVector(FACE, 0.25));
        return;
      }
      frame({
        iso:   ISO,
        front: new THREE.Vector3(0, 0, 1),
        top:   new THREE.Vector3(0, 1, 1e-4),
        right: new THREE.Vector3(1, 0, 0),
        face:  FACE,
        axial: AXIS.clone().addScaledVector(FACE, 1e-3)}[b.dataset.v]);
    };
  });
  document.getElementById('cadfit').onclick = () => frame(ISO);

  const clickView = v => {
    const b = document.querySelector('#cadtools button[data-v="' + v + '"]');
    if (b) b.click();
  };

  addEventListener('keydown', e => {
    // Never steal a key from a field the user is typing in: the number inputs in
    // the edit panel are full of digits, and 'e' is legal in a float literal.
    const t = e.target, tag = t && t.tagName;
    if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA') {
      if (e.key === 'Escape') t.blur();
      return;
    }
    if (e.ctrlKey || e.metaKey || e.altKey) return;
    const k = e.key;
    if (k === 'f' || k === 'F') frame(ISO);
    if (k === 'g' || k === 'G') setDragArmed(!dragArm);
    if (k === '?' || (k === '/' && e.shiftKey)) toggleKeys();
    if (k === '1') clickView('iso');
    if (k === '2') clickView('front');
    if (k === '3') clickView('top');
    if (k === '4') clickView('right');
    if (k === '5') clickView('face');
    if (k === '6') clickView('axial');
    if (k === 'w' || k === 'W') clickView('wheelview');
    if (k === 'c' || k === 'C') clickView('contact');
    if (k === 'e' || k === 'E') {
      const b = document.getElementById('cadedges');
      b.checked = !b.checked; b.onchange({target: b});
    }
    if (k === 'o' || k === 'O') {
      const b = document.getElementById('cadortho');
      b.checked = !b.checked;
      b.onchange({target: b});      // the handler reads e.target.checked
    }
    if (k === 'x' || k === 'X') {
      // Cycle off -> X -> Z -> Y -> off, so one key walks every section plane.
      const s = document.getElementById('cadaxis');
      s.selectedIndex = (s.selectedIndex + 1) %% s.options.length;
      s.onchange();
    }
    if (k === 'Escape') {
      if (document.getElementById('cadkeys').style.display === 'block') toggleKeys();
      else cancelDrag();
    }
    if (!dragArm || !EDIT || !EDIT.block) return;
    const nk = {ArrowLeft: ['arc', -1], ArrowRight: ['arc', 1],
                ArrowUp: ['rad', 1], ArrowDown: ['rad', -1]}[k];
    if (nk) { nudge(nk[0], nk[1]); e.preventDefault(); }
  });

  function toggleKeys() {
    const p = document.getElementById('cadkeys');
    p.style.display = (p.style.display === 'block') ? 'none' : 'block';
  }
  document.getElementById('cadkeysbtn').onclick = toggleKeys;
  document.getElementById('cadkeys').onclick = toggleKeys;
  document.getElementById('cadwhere').onclick = () => clickView('contact');

  // ---- viewport size -------------------------------------------------------
  const fullBtn = document.getElementById('cadfull');
  fullBtn.onclick = () => {
    if (document.fullscreenElement) document.exitFullscreen();
    else if (wrap.requestFullscreen) wrap.requestFullscreen();
  };
  // A notebook iframe may refuse fullscreen (no allow="fullscreen"). Say so once
  // rather than leaving a button that silently does nothing.
  wrap.addEventListener('fullscreenerror', () => {
    fullBtn.title = 'this notebook frame does not permit fullscreen; ' +
                    'drag the bottom edge instead';
    fullBtn.style.opacity = '0.45';
  });
  document.addEventListener('fullscreenchange', () => {
    fullBtn.innerHTML = document.fullscreenElement ? '&#x2715;' : '&#x26F6;';
    setTimeout(resize, 60);
  });
  canvas.addEventListener('dblclick', () => fullBtn.onclick());

  // Drag the bottom edge to resize. Persisted per notebook so the choice survives a
  // re-run of the cell, which is what makes it feel like a setting and not a fidget.
  const HKEY = 'semgrit.cadviewer.height';
  try {
    const kept = parseInt(localStorage.getItem(HKEY) || '', 10);
    if (kept >= 260 && kept <= 4000) wrap.style.height = kept + 'px';
  } catch (err) { /* storage blocked in a sandboxed frame; the default stands */ }
  (() => {
    const grip = document.getElementById('cadgrip');
    let y0 = 0, h0 = 0;
    grip.addEventListener('pointerdown', ev => {
      y0 = ev.clientY; h0 = wrap.clientHeight;
      grip.setPointerCapture(ev.pointerId);
      ev.preventDefault();
    });
    grip.addEventListener('pointermove', ev => {
      if (!grip.hasPointerCapture(ev.pointerId)) return;
      const h = Math.max(260, Math.min(4000, h0 + (ev.clientY - y0)));
      wrap.style.height = h + 'px';
      resize();
    });
    grip.addEventListener('pointerup', ev => {
      grip.releasePointerCapture(ev.pointerId);
      try { localStorage.setItem(HKEY, String(wrap.clientHeight)); } catch (err) {}
    });
  })();

  // ---- collapsible panel sections -----------------------------------------
  // Everything after a .cadhd up to the next .cadhd is that section's body. Called
  // AFTER the tree, BC and edit panels are populated -- they are built inside the
  // GLTF callback, so running this at module scope would bind empty sections.
  function makeCollapsible() {
    document.querySelectorAll('#cadtools .cadhd').forEach(hd => {
      if (hd.dataset.coll) return;                  // idempotent: safe to re-run
      hd.dataset.coll = '1';
      const car = document.createElement('span');
      car.className = 'car'; car.textContent = '▼';
      hd.insertBefore(car, hd.firstChild);
      const body = [];
      for (let n = hd.nextElementSibling; n && !n.classList.contains('cadhd');
           n = n.nextElementSibling) { n.classList.add('cadsec'); body.push(n); }
      const key = 'semgrit.cadviewer.shut.' + hd.textContent.trim();
      const set = shut => {
        hd.classList.toggle('shut', shut);
        body.forEach(n => n.classList.toggle('shut', shut));
        try { localStorage.setItem(key, shut ? '1' : '0'); } catch (err) {}
      };
      let shut = false;
      try { shut = localStorage.getItem(key) === '1'; } catch (err) {}
      if (shut) set(true);
      hd.onclick = () => set(!hd.classList.contains('shut'));
    });
  }

  // ---- first-run hint ------------------------------------------------------
  (() => {
    const f = document.getElementById('cadfirst');
    if (!f) return;
    let seen = false;
    try { seen = localStorage.getItem('semgrit.cadviewer.seen') === '1'; } catch (e) {}
    if (seen) { f.remove(); return; }
    const go = () => {
      f.remove();
      try { localStorage.setItem('semgrit.cadviewer.seen', '1'); } catch (e) {}
    };
    setTimeout(go, 12000);
    wrap.addEventListener('pointerdown', go, {once: true});
  })();

  document.getElementById('cadedges').onchange = e => {
    parts.forEach(p => {
      // A heavy mesh started with an empty edge buffer. Build it on demand the first
      // time edges are forced on, so the cost is paid only if it is asked for.
      if (e.target.checked && p.edges.userData.heavy) {
        p.edges.geometry.dispose();
        p.edges.geometry = new THREE.EdgesGeometry(p.mesh.geometry, 24);
        p.edges.userData.heavy = false;
      }
      p.edges.visible = e.target.checked;
    });
  };
  document.getElementById('cadortho').onchange = e => {
    const t = controls.target.clone(), o = cam.position.clone();
    cam = e.target.checked ? ortho : persp;
    controls.object = cam;
    cam.position.copy(o); controls.target.copy(t);
    frame(o.sub(t));
  };
  let spin = false;
  document.getElementById('cadspin').onchange = e => { spin = e.target.checked; };

  function applyClip() {
    const a = parseInt(axisSel.value, 10);
    if (a < 0 || !deckBox) { clip.constant = 1e9; return; }
    const s = flipBox.checked ? 1 : -1;
    const n = [new THREE.Vector3(s, 0, 0), new THREE.Vector3(0, s, 0),
               new THREE.Vector3(0, 0, s)][a];
    const lo = deckBox.min.getComponent(a), hi = deckBox.max.getComponent(a);
    const at = lo + (hi - lo) * (cutSlide.value / 1000);
    clip.setFromNormalAndCoplanarPoint(n, deckBox.min.clone().setComponent(a, at));
  }
  axisSel.onchange = () => { applyClip(); syncCap(); };
  cutSlide.oninput = () => { applyClip(); syncCap(); };
  flipBox.onchange = () => { applyClip(); syncCap(); };

  // ---- capping the cut face ------------------------------------------------
  // Clipping alone leaves the section hollow: you look through the cut and see the
  // inside of the far surface, which reads as a broken model rather than a
  // section. The standard fix is a stencil pass -- draw the back faces of the
  // clipped geometry into the stencil buffer, the front faces out of it, and what
  // remains set is exactly the solid interior, which a full-screen quad then
  // fills. Cheap, exact, and it needs no CSG.
  let capPlane = null, capStencil = [];
  const CAP_COL = 0xcfd6dd;

  function buildCap() {
    if (capPlane) return;
    const g = new THREE.PlaneGeometry(1, 1);
    const m = new THREE.MeshStandardMaterial({
      color: CAP_COL, metalness: 0.0, roughness: 0.85,
      side: THREE.DoubleSide,
      stencilWrite: true, stencilRef: 0,
      stencilFunc: THREE.NotEqualStencilFunc,
      stencilFail: THREE.ReplaceStencilOp,
      stencilZFail: THREE.ReplaceStencilOp,
      stencilZPass: THREE.ReplaceStencilOp});
    capPlane = new THREE.Mesh(g, m);
    capPlane.renderOrder = 3;
    capPlane.visible = false;
    // The cap must not be clipped by the very plane it fills, and it must not be
    // pickable -- a click on the section should still reach the geometry behind.
    capPlane.raycast = () => {};
    scene.add(capPlane);

    // One stencil pair per part, sharing the part's geometry.
    parts.forEach(p => {
      if (p.isGhost) return;
      const mk = (side, zpass) => {
        const mm = new THREE.MeshBasicMaterial();
        mm.depthWrite = false; mm.depthTest = false;
        mm.colorWrite = false; mm.stencilWrite = true;
        mm.stencilFunc = THREE.AlwaysStencilFunc;
        mm.side = side;
        mm.clippingPlanes = [clip];
        mm.stencilFail = THREE.KeepStencilOp;
        mm.stencilZFail = zpass;
        mm.stencilZPass = zpass;
        const o = new THREE.Mesh(p.mesh.geometry, mm);
        o.renderOrder = 1;
        o.raycast = () => {};
        capStencil.push({obj: o, src: p});
        scene.add(o);
        return o;
      };
      mk(THREE.BackSide, THREE.IncrementWrapStencilOp);
      mk(THREE.FrontSide, THREE.DecrementWrapStencilOp);
    });
  }

  function syncCap() {
    const on = document.getElementById('cadcap').checked
               && parseInt(axisSel.value, 10) >= 0;
    if (on) buildCap();
    if (!capPlane) return;
    capPlane.visible = on;
    // A stencil helper follows the visibility of the part it stands for, or a
    // hidden part would still punch a hole in the cap.
    capStencil.forEach(s => { s.obj.visible = on && s.src.mesh.visible; });
    if (!on) return;
    // Sit the quad on the clip plane, big enough to cover the model from any angle.
    const n = clip.normal, d = -clip.constant;
    capPlane.position.copy(n).multiplyScalar(d);
    capPlane.lookAt(capPlane.position.clone().add(n));
    const s = (sphere ? sphere.radius : 1) * 4;
    capPlane.scale.set(s, s, 1);
  }
  document.getElementById('cadcap').onchange = syncCap;

  // ---- explode -------------------------------------------------------------
  // Pull the parts apart along the radius, which is the axis they are stacked on:
  // bond, then grits, then the workpiece outside them. Standard CAD, and here it
  // is the quickest way to see that the grits really do stand proud of the bond
  // and that the block sits clear of both.
  const EXPLODE_ORDER = {bond_rim: -1.0, abrasive_grits: 0.35, workpiece: 1.0};
  function applyExplode() {
    const r = document.getElementById('cadexplode');
    const t = Number(r.value) / 100.0;
    const span = deckBox ? deckBox.getSize(new THREE.Vector3()).length() : 1;
    parts.forEach(p => {
      let k = EXPLODE_ORDER[p.name];
      if (k === undefined) k = p.isGhost ? 0.0 : 0.35;
      p.mesh.position.copy(FACE).multiplyScalar(t * k * span * 0.55);
    });
    // The stencil helpers shadow their parts, so they have to move with them or
    // the cap is punched by geometry that is no longer there.
    capStencil.forEach(s2 => { s2.obj.position.copy(s2.src.mesh.position); });
    const txt = document.getElementById('cadexplodetxt');
    if (txt) {
      txt.textContent = t <= 0 ? 'assembled'
        : 'exploded ' + (t * span * 0.55 * 1000).toFixed(0) + ' um along the radius';
    }
  }
  document.getElementById('cadexplode').oninput = applyExplode;

  // ---- pick a grain, and measure -----------------------------------------
  const ray = new THREE.Raycaster();
  const mpts = [];
  let mline = null;
  const dots = new THREE.Group(); scene.add(dots);

  function clearMeasure() {
    mpts.length = 0;
    lastMeasure = null;
    dots.clear();
    if (mline) { scene.remove(mline); mline = null; }
    document.getElementById('cadhelp').textContent =
      'drag orbit · right-drag pan · scroll zoom';
  }
  document.getElementById('cadclear').onclick = clearMeasure;

  function grainAt(list, faceIndex) {
    // triangle index -> grain, through the ranges recorded when the mesh was merged
    let lo = 0, hi = (list || []).length - 1;
    while (lo <= hi) {
      const m = (lo + hi) >> 1, g = list[m];
      if (faceIndex < g.tri0) hi = m - 1;
      else if (faceIndex >= g.tri0 + g.ntri) lo = m + 1;
      else return g;
    }
    return null;
  }

  function highlight(g, mesh) {
    if (hilite) { hilite.parent.remove(hilite); hilite.geometry.dispose();
                  hilite = null; }
    if (!g || !mesh) return;
    const idx = mesh.geometry.index.array;
    const sub = new THREE.BufferGeometry();
    sub.setAttribute('position', mesh.geometry.attributes.position);
    sub.setIndex(Array.from(idx.subarray(g.tri0 * 3, (g.tri0 + g.ntri) * 3)));
    hilite = new THREE.Mesh(sub, new THREE.MeshBasicMaterial({
      color:0xff9500, side:THREE.DoubleSide, clippingPlanes:[clip],
      polygonOffset:true, polygonOffsetFactor:-2, polygonOffsetUnits:-2}));
    hilite.renderOrder = 2;
    mesh.add(hilite);
  }

  // glTF is Y-up; the deck is Z-up. Report in the deck's frame or the numbers mean
  // nothing next to the .inp.
  const toModel = v => new THREE.Vector3(v.x, -v.z, v.y);
  let lastMeasure = null;
  const fmt = (mm) => (mm * 1000).toFixed(3) + ' &micro;m';

  renderer.domElement.addEventListener('pointermove', ev => {
    if (dragging) { moveDrag(ev); ev.preventDefault(); }
  });
  renderer.domElement.addEventListener('pointerup', () => { if (dragging) endDrag(); });
  renderer.domElement.addEventListener('pointerleave', () => {
    if (dragging) endDrag();
  });

  renderer.domElement.addEventListener('pointerdown', ev => {
    if (!root || ev.button !== 0) return;
    if (dragArm && !ev.altKey) {
      // Only start a drag if the pointer is actually on the block, so an armed mode
      // does not hijack every orbit.
      const r0 = renderer.domElement.getBoundingClientRect();
      ray.setFromCamera(new THREE.Vector2(
        ((ev.clientX - r0.left) / r0.width) * 2 - 1,
        -((ev.clientY - r0.top) / r0.height) * 2 + 1), cam);
      const wpm = parts.find(p => p.key === 'workpiece');
      if (wpm && wpm.mesh.visible && ray.intersectObject(wpm.mesh, false).length) {
        if (beginDrag(ev)) { ev.preventDefault(); return; }
      }
    }
    const r = renderer.domElement.getBoundingClientRect();
    ray.setFromCamera(new THREE.Vector2(
      ((ev.clientX - r.left) / r.width) * 2 - 1,
      -((ev.clientY - r.top) / r.height) * 2 + 1), cam);
    const all = ray.intersectObjects(
      parts.filter(p => p.mesh.visible).map(p => p.mesh), false);
    // A see-through part should not swallow the click. From the Face view the
    // workpiece sits between the camera and every grain, so prefer the first opaque
    // surface and fall back to the transparent one only when there is nothing behind.
    const hits = all.length
      ? [all.find(h => !h.object.material.transparent) || all[0]]
      : [];

    if (ev.shiftKey) {                                   // measurement
      if (!hits.length) return;
      mpts.push(hits[0].point.clone());
      const d = new THREE.Mesh(
        new THREE.SphereGeometry((sphere.radius || 1) * 0.006, 12, 8),
        new THREE.MeshBasicMaterial({color:0xd11a2a}));
      d.position.copy(hits[0].point); dots.add(d);
      if (mpts.length === 2) {
        if (mline) scene.remove(mline);
        mline = new THREE.Line(
          new THREE.BufferGeometry().setFromPoints(mpts),
          new THREE.LineBasicMaterial({color:0xd11a2a}));
        scene.add(mline);
        const L = mpts[0].distanceTo(mpts[1]);
        const dv = new THREE.Vector3().subVectors(mpts[1], mpts[0]);
        const dm = toModel(dv);
        // On a curved part the useful directions are not the global axes but the
        // contact frame: radially into the wheel, along the arc, across the face.
        const dr = dv.dot(FACE), da = dv.dot(ARC), dz = dv.dot(AXIS);
        // Also expose it as numbers: the panel is for reading, this is for scripting
        // and for the test that checks the components really add up to the distance.
        lastMeasure = {total_mm:L, dx_mm:dm.x, dy_mm:dm.y, dz_mm:dm.z,
                       radial_mm:dr, along_arc_mm:da, across_face_mm:dz};
        PICK.innerHTML =
          '<b>distance ' + (L * 1000).toFixed(3) + ' &micro;m</b>' +
          '<div style="color:#5a636c">' + L.toFixed(6) + ' mm</div>' +
          '<div class="cadhd">components, global</div>' +
          row('&Delta;X', dm.x * 1000, '&micro;m') +
          row('&Delta;Y', dm.y * 1000, '&micro;m') +
          row('&Delta;Z', dm.z * 1000, '&micro;m') +
          '<div class="cadhd">components, contact frame</div>' +
          row('radial', dr * 1000, '&micro;m') +
          row('along arc', da * 1000, '&micro;m') +
          row('across face', dz * 1000, '&micro;m');
        document.getElementById('cadhelp').textContent =
          'distance ' + (L * 1000).toFixed(3) + ' µm — see the panel for components';
        mpts.length = 0;
      } else {
        document.getElementById('cadhelp').textContent =
          'first point set — shift-click a second';
      }
      return;
    }

    if (!hits.length) { highlight(null); PICK.innerHTML =
        '<span style="color:#5a636c">click a grain</span>'; return; }
    const h = hits[0];
    const key = sane(h.object.name);
    const isGrit = key === 'abrasive_grits' || key === FAR_KEY;
    if (!isGrit) {
      highlight(null);
      const wp = META.workpiece_mm;
      PICK.innerHTML = '<b>' + (PRETTY[key] || key) + '</b>' + (
        key === 'workpiece' && wp
          ? '<br>' + wp[0] + ' &times; ' + wp[1] + ' &times; ' + wp[2] + ' mm'
          : key === 'bond_rim'
            ? '<br>outer radius ' + META.outer_radius_mm + ' mm' +
              '<br>bond clearance ' + META.bond_clearance_um.toFixed(3) + ' &micro;m'
            : '');
      return;
    }
    const far = key === FAR_KEY;
    const g = grainAt(far ? META.grains_far : META.grains, h.faceIndex);
    highlight(g, h.object);
    PICK.innerHTML = g
      ? '<b>grain ' + g.id + '</b>' +
        (far ? '<div style="color:#8a5a00">drawn as a box to save space &mdash; the '
             + 'numbers below are still the real grain</div>' : '') +
        row('protrusion', g.protrusion_um, '&micro;m') +
        row('height',     g.height_um,     '&micro;m') +
        row('width',      g.width_um,      '&micro;m') +
        row('volume',     g.volume_um3,    '&micro;m&sup3;') +
        row('along arc',  g.b_um,          '&micro;m') +
        row('across face', g.z_um,         '&micro;m')
      : '<b>abrasive grits</b>';
  });
  function row(k, v, u) {
    const t = (typeof v === 'number') ? v.toFixed(3) : v;
    return '<div style="display:flex;justify-content:space-between;gap:8px"><span>' +
           k + '</span><span><b>' + t + '</b> ' + u + '</span></div>';
  }

  // ---- boundary conditions ------------------------------------------------
  // Every symbol here is placed from numbers computed in Python (semgrit/bcspec.py):
  // an anchor, a unit direction, a magnitude. No trig, no e_r reconstruction, no
  // deciding which grains engage. The browser is a draughtsman, not a solver -- which
  // is the whole reason the drawn BCs cannot drift from the ones in the deck.
  const BCSPEC = META.bc || {has_analysis: false, items: [], notes: []};
  const bcGroup = new THREE.Group();
  const bcByKind = {};
  const screenScaled = [];        // objects kept a constant size on screen
  const CONE = new THREE.ConeGeometry(0.5, 1.4, 12);
  const BALL = new THREE.SphereGeometry(0.55, 16, 12);
  const UP = new THREE.Vector3(0, 1, 0);
  const BCCOL = {encastre: 0xff8c1a, velocity: 0x1668d6, rotation: 0x8a3fc0,
                 refnode: 0x18a05a, contact: 0xe03a3a};

  function glyphLabel(text, colour) {
    const c = document.createElement('canvas');
    c.width = 256; c.height = 64;
    const x = c.getContext('2d');
    x.font = 'bold 34px system-ui';
    x.fillStyle = colour;
    x.textAlign = 'center'; x.textBaseline = 'middle';
    x.fillText(text, 128, 34);
    const s = new THREE.Sprite(new THREE.SpriteMaterial(
      {map: new THREE.CanvasTexture(c), depthTest: false, transparent: true}));
    s.scale.set(4.0, 1.0, 1.0);
    return s;
  }

  function patchMesh(quad, colour, opacity) {
    const q = quad.map(a => new THREE.Vector3().fromArray(a));
    const geo = new THREE.BufferGeometry().setFromPoints(
      [q[0], q[1], q[2], q[0], q[2], q[3]]);
    geo.computeVertexNormals();
    const mat = new THREE.MeshBasicMaterial({
      color: colour, transparent: true, opacity: opacity,
      side: THREE.DoubleSide, depthWrite: false, clippingPlanes: [clip]});
    return new THREE.Mesh(geo, mat);
  }

  function quadNormal(quad) {
    const q = quad.map(a => new THREE.Vector3().fromArray(a));
    return new THREE.Vector3().subVectors(q[1], q[0])
      .cross(new THREE.Vector3().subVectors(q[3], q[0])).normalize();
  }

  function addScaled(obj, base) {
    obj.userData.baseScale = base;
    screenScaled.push(obj);
    return obj;
  }

  function buildBC() {
    if (!root || !BCSPEC.items.length) return;
    root.add(bcGroup);
    BCSPEC.items.forEach((it, i) => {
      const g = new THREE.Group();
      const col = BCCOL[it.kind] || 0x666666;
      const hex = '#' + col.toString(16).padStart(6, '0');

      if (it.kind === 'encastre') {
        // The whole held face, plus a handful of pin symbols on it. The face patch is
        // what tells the truth about the set's extent; the pins are only legibility,
        // and the panel carries the real node count so the sampling cannot mislead.
        g.add(patchMesh(it.quad, col, 0.26));
        const n = quadNormal(it.quad);
        (it.anchors || []).forEach(a => {
          const cone = new THREE.Mesh(CONE, new THREE.MeshBasicMaterial({color: col}));
          cone.position.fromArray(a);
          cone.quaternion.setFromUnitVectors(UP, n);
          g.add(addScaled(cone, 1.0));
        });
      } else if (it.kind === 'contact') {
        g.add(patchMesh(it.quad, col, 0.22));
      } else if (it.kind === 'velocity') {
        const at = new THREE.Vector3().fromArray(it.at);
        const dir = new THREE.Vector3().fromArray(it.dir).normalize();
        // Drawn from the anchor BACK along the feed direction, so the head sits on the
        // contact and points the way the wheel is actually driven.
        const shaft = new THREE.Group();
        shaft.position.copy(at);
        const head = new THREE.Mesh(CONE, new THREE.MeshBasicMaterial({color: col}));
        head.quaternion.setFromUnitVectors(UP, dir);
        shaft.add(head);
        const tail = new THREE.Line(
          new THREE.BufferGeometry().setFromPoints(
            [dir.clone().multiplyScalar(-6.0), new THREE.Vector3()]),
          new THREE.LineBasicMaterial({color: col}));
        shaft.add(tail);
        const lab = glyphLabel(it.label, hex);
        lab.position.copy(dir.clone().multiplyScalar(-8.5));
        shaft.add(lab);
        g.add(addScaled(shaft, 1.0));
      } else if (it.kind === 'rotation') {
        // A curved arrow about the wheel axis at a fraction of the rim radius, with the
        // head on the end the surface travels TOWARD. sign < 0 means decreasing theta.
        const R = (it.radius_mm || 1) * 0.55;
        const ax = new THREE.Vector3().fromArray(it.axis).normalize();
        const arc = new THREE.Mesh(
          new THREE.TorusGeometry(R, R * 0.012, 8, 64, Math.PI * 1.15),
          new THREE.MeshBasicMaterial({color: col}));
        const spin = new THREE.Group();
        // The torus lies in its own XY plane; rotate that plane onto the axis.
        spin.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), ax);
        spin.add(arc);
        const end = Math.PI * 1.15 * (it.sign < 0 ? 0 : 1);
        const head = new THREE.Mesh(CONE, new THREE.MeshBasicMaterial({color: col}));
        head.position.set(R * Math.cos(end), R * Math.sin(end), 0);
        // tangent at that end, reversed when the sense is reversed
        const tang = new THREE.Vector3(-Math.sin(end), Math.cos(end), 0)
          .multiplyScalar(it.sign < 0 ? -1 : 1);
        head.quaternion.setFromUnitVectors(UP, tang);
        spin.add(addScaled(head, 1.6));
        const mid = Math.PI * 1.15 * 0.5;
        const lab = glyphLabel(it.label, hex);
        lab.position.set(R * 1.14 * Math.cos(mid), R * 1.14 * Math.sin(mid), 0);
        spin.add(addScaled(lab, 1.0));
        // Remember the head so the sense can be checked numerically rather than by eye.
        rotCheck = {at: head.position.clone(), tang: tang.clone(),
                    axisLocal: new THREE.Vector3(0, 0, 1), group: spin,
                    sign: it.sign, magnitude: it.magnitude};
        g.add(spin);
      } else if (it.kind === 'refnode') {
        const at = new THREE.Vector3().fromArray(it.at);
        const ball = new THREE.Mesh(BALL, new THREE.MeshBasicMaterial({color: col}));
        ball.position.copy(at);
        g.add(addScaled(ball, 1.0));
        // A leader to the contact, because the reference node is on the axis and the
        // contact is microns wide: zoomed in, the wheel's only BC is otherwise off
        // screen and the model looks unconstrained.
        if (it.leader) {
          g.add(new THREE.Line(
            new THREE.BufferGeometry().setFromPoints(
              [at, new THREE.Vector3().fromArray(it.leader)]),
            new THREE.LineDashedMaterial({color: col, dashSize: 0.4,
                                          gapSize: 0.25})));
          g.children[g.children.length - 1].computeLineDistances();
        }
        const lab = glyphLabel(it.label, hex);
        lab.position.copy(at);
        g.add(addScaled(lab, 1.0));
      }
      g.userData.kind = it.kind;
      g.userData.index = i;
      bcGroup.add(g);
      (bcByKind[it.kind] = bcByKind[it.kind] || []).push(g);
    });

    // Light the grains the deck says can engage. Reuses the existing index-only
    // sub-geometry trick, so this is the merged grit mesh's own triangles -- the very
    // facets that go into ES_GRITS_ENGAGE.
    const eng = (META.grains || []).filter(x => x.engage);
    if (eng.length && gritMesh) {
      const idx = gritMesh.geometry.index.array;
      const list = [];
      eng.forEach(x => {
        for (let k = x.tri0 * 3; k < (x.tri0 + x.ntri) * 3; k++) list.push(idx[k]);
      });
      const sub = new THREE.BufferGeometry();
      sub.setAttribute('position', gritMesh.geometry.attributes.position);
      sub.setIndex(list);
      engageMesh = new THREE.Mesh(sub, new THREE.MeshBasicMaterial({
        color: 0xe03a3a, side: THREE.DoubleSide, clippingPlanes: [clip],
        transparent: true, opacity: 0.85,
        polygonOffset: true, polygonOffsetFactor: -3, polygonOffsetUnits: -3}));
      engageMesh.renderOrder = 2;
      engageMesh.visible = false;
      gritMesh.add(engageMesh);
    }

    // the panel
    const box = document.getElementById('cadbclist');
    if (!BCSPEC.has_analysis) {
      box.innerHTML = '<div style="color:#8a5a00">' +
        (BCSPEC.notes[0] || 'no boundary conditions in this deck') + '</div>';
      return;
    }
    const seen = [];
    BCSPEC.items.forEach(it => {
      if (seen.indexOf(it.kind) >= 0) return;
      seen.push(it.kind);
      const id = 'cadbc_' + it.kind;
      const lb = document.createElement('label');
      lb.innerHTML = '<input type="checkbox" id="' + id + '" checked> ' +
        '<span class="sw" style="background:' + hexOf(it.kind) + '"></span>' +
        nameOf(it.kind);
      lb.title = it.detail;
      box.appendChild(lb);
      lb.querySelector('input').onchange = e => {
        (bcByKind[it.kind] || []).forEach(gg => { gg.visible = e.target.checked; });
        if (it.kind === 'contact' && engageMesh) {
          engageMesh.visible = e.target.checked && showEngage;
        }
      };
    });
    const eb = document.createElement('label');
    eb.innerHTML = '<input type="checkbox" id="cadbc_eng"> highlight engaging grits';
    box.appendChild(eb);
    eb.querySelector('input').onchange = e => {
      showEngage = e.target.checked;
      if (engageMesh) engageMesh.visible = showEngage;
    };
    // The detail is one full sentence per boundary condition, read out of the deck.
    // Left open it is the single thing that overflows the 206 px column -- five
    // items push the panel past the bottom of the viewer and the last one is cut
    // off mid-word. Fold it behind a disclosure: the checkboxes above are the
    // control, this is the evidence, and evidence can be asked for.
    //
    // A warning note is different: it is the reason a run would be wasted, so it
    // stays outside the fold and stays red.
    const det = document.createElement('details');
    det.style.cssText = 'margin-top:5px;color:#5a636c;font-size:11px';
    const sum = document.createElement('summary');
    sum.textContent = 'what each one writes (' + BCSPEC.items.length + ')';
    sum.style.cssText = 'cursor:pointer;user-select:none;color:#39424b';
    det.appendChild(sum);
    const body = document.createElement('div');
    body.style.cssText = 'max-height:190px;overflow:auto;margin-top:3px';
    body.innerHTML = BCSPEC.items.map(it =>
      '<div style="margin-bottom:3px"><b style="color:' + hexOf(it.kind) + '">' +
      it.label + '</b> ' + esc(it.detail) + '</div>').join('');
    det.appendChild(body);
    box.appendChild(det);
    if (BCSPEC.notes.length) {
      const n = document.createElement('div');
      n.style.cssText = 'color:#b00;font-size:11px;margin-top:4px';
      n.textContent = BCSPEC.notes.join(' ');
      box.appendChild(n);
    }
  }
  let engageMesh = null, showEngage = false, rotCheck = null;

  // ---- colour the grains by a measured property ----------------------------
  // Every grain's numbers are already in META, and the grits are ONE merged mesh
  // with a triangle range per grain -- the same ranges picking and the engage
  // highlight use. So a per-vertex colour attribute painted over those ranges
  // turns the wheel into a map of the dressing, which is otherwise only visible
  // one grain at a time by clicking. Vertices are not shared between grains in
  // the merged mesh, so writing per grain cannot bleed into its neighbour.

  // Viridis, sampled. Perceptually uniform and legible in greyscale, unlike the
  // rainbow a viewer usually reaches for; matches what the Python figures use.
  const VIRIDIS = [[0.267,0.005,0.329],[0.283,0.141,0.458],[0.254,0.265,0.530],
                   [0.207,0.372,0.553],[0.164,0.471,0.558],[0.128,0.567,0.551],
                   [0.135,0.659,0.518],[0.267,0.749,0.441],[0.478,0.821,0.318],
                   [0.741,0.873,0.150],[0.993,0.906,0.144]];
  function viridis(t) {
    t = Math.max(0, Math.min(1, t));
    const x = t * (VIRIDIS.length - 1), i = Math.min(Math.floor(x),
                                                     VIRIDIS.length - 2);
    const f = x - i, a = VIRIDIS[i], b = VIRIDIS[i + 1];
    return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]),
            a[2] + f * (b[2] - a[2])];
  }

  function setupColourBy() {
    const sel = document.getElementById('cadcolour');
    if (!sel) return;
    if (!gritMesh || !(META.grains || []).length) {
      sel.disabled = true;
      sel.title = 'no per-grain data in this view';
      return;
    }
    sel.onchange = () => paintColour(sel.value);
  }

  function paintColour(key) {
    const leg = document.getElementById('cadlegend');
    // Both grit meshes, so a dense wheel -- where most grains are box proxies in
    // the far mesh -- colours all the way out instead of only near the contact.
    const targets = [[gritMesh, META.grains || []],
                     [farMesh, META.grains_far || []]].filter(t => t[0]);
    if (!targets.length) return;

    if (key === 'none') {
      targets.forEach(([m]) => {
        const b = m.userData.baseColour;
        if (b !== undefined) m.material.color.setHex(b);
        m.material.vertexColors = false;
        m.material.needsUpdate = true;
      });
      leg.style.display = 'none';
      return;
    }

    // One scale across BOTH meshes: colouring each to its own range would make a
    // near grain and a far grain of the same size different colours, which is
    // exactly the comparison the map exists to support.
    const all = (META.grains || []).concat(META.grains_far || []);
    const grainsAll = all;

    // Boolean and continuous are different questions and deserve different
    // encodings: a two-value viridis ramp reads as "a bit more of something".
    const boolean = (key === 'engage');
    let lo = Infinity, hi = -Infinity;
    if (!boolean) {
      grainsAll.forEach(g => {
        const v = g[key];
        if (typeof v === 'number' && isFinite(v)) {
          if (v < lo) lo = v;
          if (v > hi) hi = v;
        }
      });
      if (!isFinite(lo)) { lo = 0; hi = 1; }
      if (hi - lo < 1e-12) hi = lo + 1e-12;
    }

    targets.forEach(([mesh, grains]) => {
      const geom = mesh.geometry, mat = mesh.material;
      if (mesh.userData.baseColour === undefined) {
        mesh.userData.baseColour = mat.color.getHex();
      }
      const n = geom.attributes.position.count;
      let attr = geom.getAttribute('color');
      if (!attr || attr.count !== n) {
        attr = new THREE.BufferAttribute(new Float32Array(n * 3), 3);
        geom.setAttribute('color', attr);
      }
      const idx = geom.index ? geom.index.array : null;
      const arr = attr.array;
      arr.fill(0.72);                     // anything not covered stays neutral grey
      grains.forEach(g => {
        let c;
        if (boolean) c = g.engage ? [0.88, 0.23, 0.23] : [0.62, 0.66, 0.70];
        else {
          const v = g[key];
          c = (typeof v === 'number' && isFinite(v))
            ? viridis((v - lo) / (hi - lo)) : [0.72, 0.72, 0.72];
        }
        for (let k = g.tri0 * 3; k < (g.tri0 + g.ntri) * 3; k++) {
          const vi = idx ? idx[k] : k;
          arr[vi * 3] = c[0]; arr[vi * 3 + 1] = c[1]; arr[vi * 3 + 2] = c[2];
        }
      });
      attr.needsUpdate = true;
      mat.vertexColors = true;
      mat.color.setHex(0xffffff);         // white, or the tint multiplies through
      mat.needsUpdate = true;
    });

    // Legend
    leg.style.display = 'block';
    const bar = document.getElementById('cadlegbar');
    if (boolean) {
      bar.style.background = 'linear-gradient(90deg,#9ea9b3 0 50%%,#e03a3a 50%% 100%%)';
      document.getElementById('cadleglo').textContent = 'no';
      document.getElementById('cadleghi').textContent = 'engages';
    } else {
      const stops = [];
      for (let i = 0; i <= 10; i++) {
        const c = viridis(i / 10);
        stops.push('rgb(' + Math.round(c[0] * 255) + ',' + Math.round(c[1] * 255)
                   + ',' + Math.round(c[2] * 255) + ') ' + (i * 10) + '%%');
      }
      bar.style.background = 'linear-gradient(90deg,' + stops.join(',') + ')';
      const unit = key === 'volume_um3' ? ' um3' : ' um';
      const fmt = v => (Math.abs(v) >= 100 ? v.toFixed(0)
                        : Math.abs(v) >= 1 ? v.toFixed(2) : v.toPrecision(2));
      document.getElementById('cadleglo').textContent = fmt(lo) + unit;
      document.getElementById('cadleghi').textContent = fmt(hi) + unit;
    }
  }

  function hexOf(k) {
    return '#' + (BCCOL[k] || 0x666666).toString(16).padStart(6, '0');
  }
  function nameOf(k) {
    return {encastre: 'held faces (ENCASTRE)', velocity: 'infeed velocity',
            rotation: 'wheel rotation', refnode: 'reference node',
            contact: 'contact surfaces'}[k] || k;
  }

  // Keep the symbols a constant size on screen, from the same span the context rule
  // uses, so they read at 50 mm wheel scale and at 5 um grain scale alike.
  function scaleGlyphs() {
    if (!screenScaled.length) return;
    const fov = THREE.MathUtils.degToRad(persp.fov);
    const span = (cam === ortho)
      ? ortho.top * 2
      : 2 * cam.position.distanceTo(controls.target) * Math.tan(fov / 2);
    const k = span * 0.018;
    screenScaled.forEach(o => o.scale.setScalar(k * (o.userData.baseScale || 1)));
  }

  // ---- editing -------------------------------------------------------------
  // The browser collects numbers and previews them by moving the box it was given. It
  // does not decide anything. The seating of the block is settled by a facet-clipping
  // tangency in Python (rigid_wheel.ground_radius) that cannot be reproduced here, so
  // a live edit is labelled a PREVIEW and the authoritative numbers come back when it
  // is applied. Every commit goes through one Python function, editable.apply.
  const EDIT = META.edit || null;
  const EDITED = {};                       // key -> value, only what the user changed
  let wpEdit = null;                       // the group the preview transform acts on
  let applyState = 'idle';

  // Colab formats a callback's return value into a mimetype bundle before the page sees
  // it, so the reply has to be read, not assumed. This used to be
  //   (res.data && res.data['application/json']) || {ok: false}
  // which turned any reply it did not recognise into a refusal with no reason -- and a
  // plain Python dict formats to text/plain only, so that was *every* reply, successful
  // ones included. Python now sends both mimetypes (editable.CommitReply); this reads
  // either, and when it can read neither it says so instead of blaming Python.
  const REPLY_SENTINEL = 'CADREPLY ';

  function readReply(res) {
    const d = (res && res.data) || (res && res['data']) || null;
    let j = d ? d['application/json'] : null;
    if (j === null || j === undefined) j = res && res['application/json'];
    if (typeof j === 'string') {
      try { j = JSON.parse(j); } catch (e) { j = null; }
    }
    if (j && typeof j === 'object') return j;
    // text/plain, but exact: Python repr()s the reply as the sentinel plus real JSON.
    const t = (d && d['text/plain']) || (res && res['text/plain']) || '';
    const at = String(t).indexOf(REPLY_SENTINEL);
    if (at >= 0) {
      try { return JSON.parse(String(t).slice(at + REPLY_SENTINEL.length)); }
      catch (e) {}
    }
    return {unreadable: true, raw: String(t || (res === undefined ? '' : res))};
  }

  function transport() {
    // Resolved at call time, never captured: a page reopened without a live runtime
    // still shows the viewer, and must fall back to export rather than fail mutely.
    if (window.cadTransport) return window.cadTransport;
    const k = (window.google && google.colab && google.colab.kernel) || null;
    if (!k) return null;
    return {
      name: 'colab',
      commit: async (settings) => readReply(
        await k.invokeFunction('cad.commit', [settings], {}))};
  }

  function settingsNow() {
    return Object.assign({}, (EDIT && EDIT.settings) || {}, EDITED);
  }

  function liveKeys() {
    return (EDIT ? EDIT.fields : []).filter(f => f.tier === 'live').map(f => f.key);
  }

  function pendingRebuild() {
    const rb = (EDIT ? EDIT.fields : []).filter(f => f.tier === 'rebuild')
      .map(f => f.key);
    return Object.keys(EDITED).filter(k => rb.indexOf(k) >= 0);
  }

  // Group the workpiece and the glyphs anchored to it, so one transform moves them
  // together and nothing drifts apart during a preview.
  function makeEditGroup() {
    if (!root || wpEdit) return;
    wpEdit = new THREE.Group();
    wpEdit.matrixAutoUpdate = false;
    root.add(wpEdit);
    const wpPart = parts.find(p => p.key === 'workpiece');
    if (wpPart) wpEdit.attach(wpPart.mesh);
    bcGroup.children.forEach(g => {
      if (g.userData.kind === 'encastre' || g.userData.kind === 'contact') {
        wpEdit.attach(g);
      }
    });
  }

  // The preview transform, built only from vectors Python supplied.
  function applyLive() {
    if (!wpEdit || !EDIT || !EDIT.block) return;
    // Never preview a state Python would refuse. A block typed four millimetres long on
    // a one-millimetre slice is an 83x scale on the drawn box: it swamps the view and
    // tells the user nothing the warning does not say better. Hold the last valid
    // transform and let the warning speak.
    if (fitWarnings().length) { paintEditState(); paintBand(); return; }
    const B = EDIT.basis, s0 = EDIT.settings, s = settingsNow();
    const er = new THREE.Vector3().fromArray(B.radial).normalize();
    const et = new THREE.Vector3().fromArray(B.arc).normalize();
    const ez = new THREE.Vector3().fromArray(B.axial).normalize();
    const c = new THREE.Vector3().fromArray(EDIT.block.centre);

    // scale about the block's own centre, along the supplied axes
    const sd = (s.wp_depth_mm || 1) / (s0.wp_depth_mm || 1);
    const sl = (s.wp_length_mm || 1) / (s0.wp_length_mm || 1);
    const sw = (s.wp_width_mm || 1) / (s0.wp_width_mm || 1);
    const R = new THREE.Matrix4().makeBasis(er, et, ez);
    const Rt = R.clone().transpose();
    const S = new THREE.Matrix4().makeScale(sd, sl, sw);
    const M = new THREE.Matrix4()
      .multiply(new THREE.Matrix4().makeTranslation(c.x, c.y, c.z))
      .multiply(R).multiply(S).multiply(Rt)
      .multiply(new THREE.Matrix4().makeTranslation(-c.x, -c.y, -c.z));

    // standoff: straight out along the radial direction
    const ds = ((s.clearance_um || 0) - (s0.clearance_um || 0)) / 1000.0;
    if (ds) {
      M.premultiply(new THREE.Matrix4().makeTranslation(
        er.x * ds, er.y * ds, er.z * ds));
    }
    // position along the arc: a rotation about the wheel axis, which passes through
    // the origin. Only meaningful for 'custom angle'; the named modes are decided in
    // Python, so changing the angle also switches the mode (and the panel says so).
    const dth = ((s.wp_position_deg || 0) - (s0.wp_position_deg || 0))
      * Math.PI / 180.0;
    if (dth) {
      M.premultiply(new THREE.Matrix4().makeRotationAxis(ez, dth));
    }
    wpEdit.matrix.copy(M);
    wpEdit.matrixWorldNeedsUpdate = true;

    // the rotation glyph follows the sense immediately -- it is the one BC a user can
    // flip without any re-placement
    const rev = !!s.rotation_reversed, rev0 = !!s0.rotation_reversed;
    if (rotCheck && rev !== rev0 !== rotCheck.flipped) {
      rotCheck.group.scale.z = rev !== rev0 ? -1 : 1;
      rotCheck.flipped = rev !== rev0;
    }
    paintEditState();
    paintBand();
  }

  function paintEditState() {
    const n = Object.keys(EDITED).length;
    const rb = pendingRebuild();
    const box = document.getElementById('cadeditstate');
    if (!box) return;
    if (!n) {
      box.innerHTML = '<span style="color:#5a636c">no edits</span>';
      return;
    }
    const warn = fitWarnings();
    box.innerHTML =
      '<b>' + n + ' edit' + (n > 1 ? 's' : '') + ' pending</b>' +
      (warn.length
        ? '<div style="color:#b00;margin:2px 0"><b>this will be refused:</b><br>'
          + warn.join('<br>') + '</div>'
        : '') +
      (seatingStale()
        ? '<div style="color:#8a5a00"><s>standoff, engaging grains and the depth window</s>'
          + ' recomputed on Apply - the block has moved, and only Python re-seats it.</div>'
        : '') +
      '<div style="color:#8a5a00">the box you see is a PREVIEW: Python re-seats it on '
      + 'the grains when you apply, so the standoff and depth-of-cut numbers will be '
      + 'recomputed.</div>' +
      (rb.length
        ? '<div style="color:#b00">' + rb.join(', ') + ' change which grains exist, so '
          + 'they need a full rebuild and are not previewed at all.</div>'
        : '') +
      (EDITED.wp_position_deg !== undefined
        ? '<div style="color:#8a5a00">moving along the arc sets "where on the arc" to '
          + '<b>custom angle</b>, overriding the named placement.</div>'
        : '');
  }

  function buildEditPanel() {
    const host = document.getElementById('cadeditlist');
    if (!host) return;
    if (!EDIT) {
      host.innerHTML = '<div style="color:#5a636c">this view was built without ' +
        'parameters, so nothing can be edited from here</div>';
      return;
    }
    makeEditGroup();
    ['live', 'rebuild'].forEach(tier => {
      const hd = document.createElement('div');
      hd.style.cssText = 'font-weight:600;margin:5px 0 2px;color:'
        + (tier === 'live' ? '#1668d6' : '#b0631a');
      hd.textContent = tier === 'live' ? 'previewed live' : 'needs a rebuild';
      host.appendChild(hd);
      EDIT.fields.filter(f => f.tier === tier).forEach(f => {
        const row = document.createElement('div');
        row.style.cssText = 'display:flex;justify-content:space-between;gap:4px;' +
          'align-items:center';
        const lab = document.createElement('span');
        lab.textContent = f.label + (f.unit ? ' (' + f.unit + ')' : '');
        lab.title = f.note || '';
        lab.style.cssText = 'overflow:hidden;text-overflow:ellipsis;white-space:nowrap';
        row.appendChild(lab);
        const v = EDIT.settings[f.key];
        let inp;
        if (f.choices && f.choices.length) {
          inp = document.createElement('select');
          f.choices.forEach(c => {
            const o = document.createElement('option');
            o.value = c; o.textContent = c;
            if (c === v) o.selected = true;
            inp.appendChild(o);
          });
        } else if (typeof v === 'boolean') {
          inp = document.createElement('input');
          inp.type = 'checkbox'; inp.checked = v;
        } else {
          inp = document.createElement('input');
          inp.type = 'number'; inp.value = v;
          inp.step = 'any';
          if (f.lo !== null && f.lo !== undefined) inp.min = f.lo;
          if (f.hi !== null && f.hi !== undefined) inp.max = f.hi;
        }
        inp.id = 'cadedit_' + f.key;
        inp.style.cssText = 'width:74px;font:inherit;font-size:11px';
        inp.onchange = () => {
          const raw = (inp.type === 'checkbox') ? inp.checked
            : (f.choices && f.choices.length ? inp.value : parseFloat(inp.value));
          if (raw === EDIT.settings[f.key]) delete EDITED[f.key];
          else EDITED[f.key] = raw;
          if (f.key === 'wp_position_deg' && EDITED.wp_position_deg !== undefined) {
            EDITED.wp_position = 'custom angle';
            const sel = document.getElementById('cadedit_wp_position');
            if (sel) sel.value = 'custom angle';
          }
          applyLive();
        };
        row.appendChild(inp);
        host.appendChild(row);
        if (f.key === 'depth_of_cut_um') {
          host.appendChild(depthBand());
        }
        // The wheel extent is given one way or the other. Grey the one that is not in
        // force and say which it is, because the raw sector_deg on an arc-mode wheel is
        // whatever the default happened to be -- it read 30 deg for a 2.29 deg arc.
        if (f.key === 'sector_deg' || f.key === 'arc_length_mm') {
          const live = (EDIT.settings._sector_mode === 'angle')
            ? 'sector_deg' : 'arc_length_mm';
          if (f.key !== live) {
            row.style.opacity = '0.45';
            lab.title = 'not in force: the extent is set as '
              + (live === 'angle' ? 'an angle' : 'an arc length')
              + '. Editing this switches to it.';
          } else {
            const nb = document.createElement('div');
            nb.style.cssText = 'font-size:10px;color:#5a636c;margin:-2px 0 3px';
            nb.textContent = 'in force: ' + EDIT.arc_length_mm.toFixed(4) + ' mm arc = '
              + EDIT.sector_resolved_deg.toFixed(4) + ' deg';
            host.appendChild(nb);
          }
        }
      });
    });
    const bar = document.createElement('div');
    bar.style.cssText = 'display:flex;gap:3px;flex-wrap:wrap;margin-top:6px';
    bar.innerHTML =
      '<button id="cadapply">Apply</button>' +
      '<button id="cadcopy">Copy JSON</button>' +
      '<button id="caddl">Download</button>' +
      '<button id="cadreset">Reset</button>';
    host.appendChild(bar);
    const st = document.createElement('div');
    st.id = 'cadeditstate';
    st.style.cssText = 'margin-top:4px;font-size:11px';
    host.appendChild(st);
    const out = document.createElement('textarea');
    out.id = 'cadeditjson';
    out.readOnly = true;
    out.style.cssText = 'width:100%%;height:52px;margin-top:4px;font:11px ui-monospace,' +
      'monospace;display:none';
    host.appendChild(out);
    paintEditState();
    paintBand();

    const db = document.getElementById('caddrag');
    if (db) db.onclick = () => setDragArmed(!dragArm);
    document.getElementById('cadreset').onclick = () => {
      Object.keys(EDITED).forEach(k => delete EDITED[k]);
      EDIT.fields.forEach(f => {
        const el = document.getElementById('cadedit_' + f.key);
        if (!el) return;
        if (el.type === 'checkbox') el.checked = EDIT.settings[f.key];
        else el.value = EDIT.settings[f.key];
      });
      if (wpEdit) { wpEdit.matrix.identity(); wpEdit.matrixWorldNeedsUpdate = true; }
      if (rotCheck) { rotCheck.group.scale.z = 1; rotCheck.flipped = false; }
      paintEditState();
    };
    const payload = () => JSON.stringify(
      {settings: settingsNow(), edited: EDITED,
       note: 'made in the CAD viewer; apply with semgrit.editable.apply'}, null, 1);
    // The notebook's form widgets cannot be written to from Python, so after an edit
    // they still show the old values. These are the lines to paste back, in the
    // widgets' own names and units -- the map comes from editable.Field, because two of
    // them differ (WHEEL_WIDTH_MM, and SURFACE_SPEED_M_S in m/s).
    const paramBlock = () => {
      const s = settingsNow(), out = [];
      EDIT.fields.forEach(f => {
        if (!f.widget) return;
        const v = s[f.key];
        if (EDIT.settings[f.key] === v) return;
        let txt;
        if (typeof v === 'boolean') txt = v ? 'True' : 'False';
        else if (f.choices && f.choices.length) txt = '"' + v + '"';
        else txt = String(Number((v * (f.widget_scale || 1)).toPrecision(10)));
        out.push(f.widget + ' = ' + txt);
      });
      if (out.some(l => l.indexOf('SECTOR_DEG ') === 0)) out.push('SECTOR_MODE = "angle"');
      else if (out.some(l => l.indexOf('ARC_LENGTH_MM ') === 0))
        out.push('SECTOR_MODE = "arc"');
      return out.join('\n');
    };
    document.getElementById('cadcopy').onclick = () => {
      const t = document.getElementById('cadeditjson');
      t.style.display = 'block';
      const pb = paramBlock();
      t.value = payload() + (pb ? '\n\n# paste these into the form cells above:\n'
                                  + pb : '');
      t.select();
      try { navigator.clipboard.writeText(t.value); } catch (e) {}
    };
    document.getElementById('caddl').onclick = () => {
      const a = document.createElement('a');
      a.download = 'viewer_settings.json';
      a.href = 'data:application/json;charset=utf-8,'
        + encodeURIComponent(payload());
      a.click();
    };
    document.getElementById('cadapply').onclick = async () => {
      const t = transport();
      const box = document.getElementById('cadeditstate');
      if (!Object.keys(EDITED).length) {
        box.innerHTML = '<span style="color:#5a636c">nothing to apply</span>';
        return;
      }
      if (!t) {
        // No live kernel: this is the normal case outside Colab and after a reopen.
        // Say so plainly and hand over the settings instead of failing silently.
        document.getElementById('cadcopy').click();
        box.innerHTML = '<b style="color:#8a5a00">no live Python kernel here.</b>' +
          '<div>The settings are in the box below and in viewer_settings.json - run ' +
          'the REBUILD cell to apply them.</div>';
        return;
      }
      applyState = 'busy';
      box.innerHTML = '<b>applying through ' + t.name + ' ...</b>';
      try {
        const res = await t.commit(settingsNow()) || {};
        applyState = res.ok ? 'done' : 'failed';
        if (res.ok) {
          box.innerHTML = '<b style="color:#1a7a3a">applied.</b><div>' +
            esc(res.message || 'rebuilt; re-run this cell to see the new geometry.') +
            '</div>';
        } else if (res.unreadable) {
          // Do not call this a refusal. The edit may well have been applied; what failed
          // is reading the answer, and the raw text usually still carries it.
          applyState = 'unknown';
          box.innerHTML = '<b style="color:#8a5a00">the reply could not be read.</b>' +
            '<div>The edit may or may not have been applied. Re-run this cell to see ' +
            'the current state, or use Copy JSON and the REBUILD cell.</div>' +
            (res.raw ? '<div style="color:#5a636c;word-break:break-all">' +
                       esc(res.raw).slice(0, 400) + '</div>' : '');
        } else {
          box.innerHTML = '<b style="color:#b00">Python refused it.</b><div>' +
            esc(res.error || 'no reason given') + '</div>';
        }
      } catch (err) {
        applyState = 'failed';
        box.innerHTML = '<b style="color:#b00">the commit failed.</b><div>' + esc(err) +
          '</div><div>Use Copy JSON and the REBUILD cell instead.</div>';
      }
    };
  }

  // ---- guardrails, dragging, and staleness ---------------------------------
  // The numbers here are all supplied by Python. The only arithmetic done in the
  // browser is comparing them, shifting the depth window by a standoff delta (exact --
  // a standoff lifts the ground face, so both ends move together, which plan_deck's own
  // check already proves), and turning a mouse displacement into an angle by dividing
  // by a supplied radius. Nothing decides where the block seats.
  const GUARD = {win: null, stale: false};

  function windowNow() {
    if (!EDIT || EDIT.first_contact_um === null
        || EDIT.first_contact_um === undefined) return null;
    const s0 = EDIT.settings, s = settingsNow();
    const ds = (s.clearance_um || 0) - (s0.clearance_um || 0);
    return {lo: EDIT.first_contact_um + ds, hi: EDIT.depth_ceiling_um + ds,
            shifted: ds !== 0};
  }

  // Anything derived from the seating goes stale the moment the block is re-placed,
  // because only Python can re-seat it. Stale is shown, never guessed.
  function seatingStale() {
    const keys = (EDIT ? EDIT.fields : []).filter(f => f.seating).map(f => f.key);
    return Object.keys(EDITED).some(k => keys.indexOf(k) >= 0);
  }

  function fitWarnings() {
    if (!EDIT) return [];
    const s = settingsNow(), out = [];
    const arc = (EDITED.arc_length_mm !== undefined)
      ? EDITED.arc_length_mm : EDIT.arc_length_mm;
    if (s.wp_length_mm > arc) {
      out.push('the block is ' + s.wp_length_mm.toFixed(4) + ' mm long but the wheel '
        + 'slice is only ' + arc.toFixed(4) + ' mm of arc, so it would hang off both '
        + 'ends');
    }
    if (s.wp_width_mm > s.width_mm) {
      out.push('the block is ' + s.wp_width_mm.toFixed(4) + ' mm wide but the wheel '
        + 'face is only ' + s.width_mm.toFixed(4) + ' mm, so it would overhang the '
        + 'sides');
    }
    const w = windowNow();
    if (w && !seatingStale()) {
      const ae = s.depth_of_cut_um;
      if (ae > 0 && ae <= w.lo) {
        out.push('a ' + ae.toFixed(3) + ' um depth of cut never reaches the work: the '
          + 'nearest grain is ' + w.lo.toFixed(3) + ' um clear, so the wheel would turn '
          + 'for the whole step and touch nothing');
      } else if (ae >= w.hi) {
        out.push('a ' + ae.toFixed(3) + ' um depth of cut drives the bond rim into the '
          + 'work: the face-to-bond gap is ' + w.hi.toFixed(3) + ' um');
      }
    }
    return out;
  }

  function depthBand() {
    const wrap = document.createElement('div');
    wrap.style.cssText = 'margin:2px 0 4px';
    wrap.innerHTML =
      '<div id="cadbandtrack" style="position:relative;height:11px;border:1px solid ' +
      '#a9b2bb;border-radius:3px;overflow:hidden;background:#eceff2">' +
      '<div id="cadbandlo" style="position:absolute;left:0;top:0;bottom:0;' +
      'background:#f3c9c9"></div>' +
      '<div id="cadbandok" style="position:absolute;top:0;bottom:0;' +
      'background:#bfe3c8"></div>' +
      '<div id="cadbandhi" style="position:absolute;right:0;top:0;bottom:0;' +
      'background:#f3c9c9"></div>' +
      '<div id="cadbandpin" style="position:absolute;top:-2px;width:2px;height:15px;' +
      'background:#12365e"></div></div>' +
      '<input type="range" id="cadbandrange" min="0" max="1000" value="0" ' +
      'style="width:100%%;margin:2px 0 0" title="drag to set the depth of cut">' +
      '<div id="cadbandtext" style="font-size:10px;color:#5a636c"></div>';
    // The slider was painted from paintBand but had no handler, so dragging it moved
    // the thumb, changed nothing, and then snapped back on the next repaint. Drive
    // the same field the number input drives, through the same EDITED bookkeeping.
    const r = wrap.querySelector('#cadbandrange');
    const push = (commit) => {
      const um = Number(r.value) / 1000.0;
      const key = 'depth_of_cut_um';
      if (um === EDIT.settings[key]) delete EDITED[key];
      else EDITED[key] = um;
      const num = document.getElementById('cadedit_' + key);
      if (num) num.value = String(Number(um.toPrecision(10)));
      // Repaint the band and the state on every move so the pin tracks the thumb,
      // but only re-preview the geometry on release: applyLive() re-seats the drawn
      // block and is far too heavy to run on every pixel of a drag.
      paintBand();
      if (commit) applyLive(); else paintEditState();
    };
    r.oninput = () => push(false);
    r.onchange = () => push(true);
    return wrap;
  }

  function paintBand() {
    const t = document.getElementById('cadbandtrack');
    if (!t) return;
    const w = windowNow();
    const txt = document.getElementById('cadbandtext');
    if (!w) { t.parentNode.style.display = 'none'; return; }
    const full = Math.max(w.hi * 1.15, 1e-6);
    const pc = v => Math.max(0, Math.min(100, 100 * v / full));
    const lo = pc(Math.max(w.lo, 0)), hi = pc(w.hi);
    document.getElementById('cadbandlo').style.width = lo + '%%';
    document.getElementById('cadbandok').style.left = lo + '%%';
    document.getElementById('cadbandok').style.width = Math.max(0, hi - lo) + '%%';
    document.getElementById('cadbandhi').style.left = hi + '%%';
    const ae = settingsNow().depth_of_cut_um;
    document.getElementById('cadbandpin').style.left = pc(ae) + '%%';
    const stale = seatingStale();
    // 0 means "let the build choose". Say what it will choose, or the field reads as
    // "no cut at all" -- the same trap the sector field had.
    const auto = (ae === 0 && META.depth_of_cut_um)
      ? ' (0 = automatic, which is ' + META.depth_of_cut_um.toFixed(3) + ' um)' : '';
    txt.innerHTML = stale
      ? '<span style="color:#8a5a00">window stale until Apply: the block has moved, so '
        + 'Python must re-seat it</span>'
      : 'usable between <b>' + w.lo.toFixed(3) + '</b> and <b>' + w.hi.toFixed(3)
        + '</b> um' + auto
        + (w.shifted ? ' &middot; shifted with the standoff' : '');
    t.style.opacity = stale ? '0.35' : '1';
    const r = document.getElementById('cadbandrange');
    r.max = String(Math.round(full * 1000));
    // Writing value back while the user is dragging fights the thumb: paintBand runs
    // on every input event, so the slider would jump back to the last committed
    // number mid-drag. The thumb is already where the user put it.
    if (document.activeElement !== r) r.value = String(Math.round(ae * 1000));
  }

  // ---- dragging the block --------------------------------------------------
  let dragArm = false, dragging = null, nudgeUm = 0.1;

  function setDragArmed(on) {
    dragArm = !!on;
    const b = document.getElementById('caddrag');
    if (b) b.className = dragArm ? 'on' : '';
    const h = document.getElementById('cadhelp');
    if (h) {
      h.textContent = dragArm
        ? 'drag the block along the arc - shift-drag for standoff - arrows nudge - Esc cancels'
        : 'drag orbit · right-drag pan · scroll zoom';
    }
  }

  // Screen motion -> a displacement in the block's own plane. The plane is defined by
  // the supplied radial normal through the supplied block centre; the result is
  // resolved onto the supplied arc and radial vectors. No frame is reconstructed here.
  function planePoint(ev) {
    if (!EDIT || !EDIT.block) return null;
    const r = renderer.domElement.getBoundingClientRect();
    ray.setFromCamera(new THREE.Vector2(
      ((ev.clientX - r.left) / r.width) * 2 - 1,
      -((ev.clientY - r.top) / r.height) * 2 + 1), cam);
    const n = new THREE.Vector3().fromArray(EDIT.basis.radial).normalize();
    const c = new THREE.Vector3().fromArray(EDIT.block.centre);
    const pl = new THREE.Plane().setFromNormalAndCoplanarPoint(n, c);
    const hit = new THREE.Vector3();
    return ray.ray.intersectPlane(pl, hit) ? hit : null;
  }

  function beginDrag(ev) {
    const p = planePoint(ev);
    if (!p) return false;
    dragging = {
      from: p.clone(), shift: !!ev.shiftKey,
      th0: settingsNow().wp_position_deg || 0,
      cl0: settingsNow().clearance_um || 0,
      before: Object.assign({}, EDITED)};
    controls.enabled = false;
    return true;
  }

  function moveDrag(ev) {
    if (!dragging) return;
    const p = planePoint(ev);
    if (!p) return;
    const d = new THREE.Vector3().subVectors(p, dragging.from);
    if (dragging.shift) {
      const er = new THREE.Vector3().fromArray(EDIT.basis.radial).normalize();
      setField('clearance_um',
               Math.max(0, dragging.cl0 + d.dot(er) * 1000.0));
    } else {
      const et = new THREE.Vector3().fromArray(EDIT.basis.arc).normalize();
      // arc length -> angle, by the radius Python supplied. A chord for an arc, which
      // over a fraction of a millimetre is sub-nanometre and only ever feeds the
      // preview: the committed number is the field's.
      const dth = (d.dot(et) / (EDIT.ground_radius_mm || 1)) * 180.0 / Math.PI;
      setField('wp_position_deg', dragging.th0 + dth);
    }
    showDragReadout();
  }

  function endDrag() {
    dragging = null;
    controls.enabled = true;
  }

  function cancelDrag() {
    if (!dragging) return;
    const keep = dragging.before;
    Object.keys(EDITED).forEach(k => delete EDITED[k]);
    Object.keys(keep).forEach(k => { EDITED[k] = keep[k]; });
    syncInputs();
    endDrag();
    applyLive();
  }

  function showDragReadout() {
    const s = settingsNow(), s0 = EDIT.settings;
    const dth = (s.wp_position_deg || 0) - (s0.wp_position_deg || 0);
    const along = dth * Math.PI / 180.0 * (EDIT.ground_radius_mm || 1) * 1000.0;
    const h = document.getElementById('cadhelp');
    if (h) {
      h.textContent = 'along arc ' + (along >= 0 ? '+' : '') + along.toFixed(2)
        + ' um (theta ' + (s.wp_position_deg || 0).toFixed(4) + ' deg) · standoff '
        + (s.clearance_um || 0).toFixed(3) + ' um';
    }
  }

  function setField(key, value) {
    const el = document.getElementById('cadedit_' + key);
    if (!el) return;
    el.value = (Math.round(value * 1e6) / 1e6);
    el.onchange();
  }

  function syncInputs() {
    (EDIT ? EDIT.fields : []).forEach(f => {
      const el = document.getElementById('cadedit_' + f.key);
      if (!el) return;
      const v = (EDITED[f.key] !== undefined) ? EDITED[f.key] : EDIT.settings[f.key];
      if (el.type === 'checkbox') el.checked = !!v; else el.value = v;
    });
  }

  function nudge(axis, sign) {
    const s = settingsNow();
    if (axis === 'arc') {
      const dth = (sign * nudgeUm / 1000.0 / (EDIT.ground_radius_mm || 1))
        * 180.0 / Math.PI;
      setField('wp_position_deg', (s.wp_position_deg || 0) + dth);
    } else {
      setField('clearance_um',
               Math.max(0, (s.clearance_um || 0) + sign * nudgeUm));
    }
    showDragReadout();
  }

  // ---- orientation triad --------------------------------------------------
  // Labelled in the model's Z-up frame, not glTF's Y-up, so the axes mean what the
  // deck means by them.
  const gizScene = new THREE.Scene();
  const gizCam = new THREE.OrthographicCamera(-1.7, 1.7, 1.7, -1.7, 0.1, 20);
  function label(text, color) {
    const c = document.createElement('canvas'); c.width = c.height = 64;
    const x = c.getContext('2d');
    x.fillStyle = color; x.font = 'bold 44px system-ui';
    x.textAlign = 'center'; x.textBaseline = 'middle'; x.fillText(text, 32, 34);
    const s = new THREE.Sprite(new THREE.SpriteMaterial(
      {map:new THREE.CanvasTexture(c), depthTest:false}));
    s.scale.setScalar(0.62);
    return s;
  }
  // (glTF direction, name, colour) -- glTF +Y is the wheel axis Z, glTF -Z is Y
  [[new THREE.Vector3(1, 0, 0), 'X', '#c0392b'],
   [new THREE.Vector3(0, 0, -1), 'Y', '#1e8449'],
   [new THREE.Vector3(0, 1, 0), 'Z', '#1f5fa8']].forEach(([d, n, col]) => {
    const c = new THREE.Color(col);
    gizScene.add(new THREE.ArrowHelper(d, new THREE.Vector3(), 1.0, c, 0.3, 0.18));
    const s = label(n, col); s.position.copy(d).multiplyScalar(1.42);
    gizScene.add(s);
  });

  function resize() {
    const w = wrap.clientWidth, h = wrap.clientHeight;
    renderer.setSize(w, h, false);
    persp.aspect = w / h; persp.updateProjectionMatrix();
    const t = ortho.top;
    ortho.left = -t * w / h; ortho.right = t * w / h; ortho.updateProjectionMatrix();
  }
  new ResizeObserver(resize).observe(wrap); resize();

  document.getElementById('cadshot').onclick = () => {
    render();
    const a = document.createElement('a');
    a.download = 'wheel_view.png';
    a.href = renderer.domElement.toDataURL('image/png');
    a.click();
  };

  // The ghost disc and the contact marker are pointers sized for a 50 mm view. Zoomed
  // into the contact they are metres wide and swamp everything, so they fade out once
  // the view is much smaller than the wheel. The checkbox still wins: unticking hides
  // them for good, ticking only offers them back at a scale where they mean something.
  function applyContextRule() {
    const fov = THREE.MathUtils.degToRad(persp.fov);
    const span = (cam === ortho)
      ? ortho.top * 2
      : 2 * cam.position.distanceTo(controls.target) * Math.tan(fov / 2);
    // Keyed to the marker's own size, not the wheel's: that is what decides when a
    // pointer stops being a hint and starts being an obstruction.
    const mk = META.mark_size_mm || (0.012 * (META.wheel_radius_mm || 1));
    const on = span > 12.0 * mk;
    parts.forEach(p => {
      p.mesh.visible = p.wanted && (!p.isGhost || on);
    });
    placeWhere(on);
  }

  // The "you are here" callout follows the same rule as the ghost: it is a hint for
  // when the contact is too small to see, and an obstruction once you have arrived.
  function placeWhere(far) {
    const el = document.getElementById('cadwhere');
    if (!el) return;
    if (!far || !CONTACT) { el.style.display = 'none'; return; }
    const v = CONTACT.clone().project(cam);
    // Behind the camera, or off screen: nothing useful to point at.
    if (v.z > 1 || Math.abs(v.x) > 1.06 || Math.abs(v.y) > 1.06) {
      el.style.display = 'none';
      return;
    }
    el.style.display = 'block';
    el.style.left = ((v.x * 0.5 + 0.5) * wrap.clientWidth) + 'px';
    el.style.top = ((-v.y * 0.5 + 0.5) * wrap.clientHeight) + 'px';
  }

  function render() {
    applyContextRule();
    scaleGlyphs();
    renderer.setViewport(0, 0, wrap.clientWidth, wrap.clientHeight);
    renderer.setScissorTest(false);
    renderer.autoClear = true;
    renderer.render(scene, cam);
    // triad, drawn last into its own corner viewport
    const s = 92;
    renderer.autoClear = false;
    renderer.setViewport(wrap.clientWidth - s - 10, 10, s, s);
    renderer.setScissor(wrap.clientWidth - s - 10, 10, s, s);
    renderer.setScissorTest(true);
    renderer.clearDepth();
    gizCam.position.copy(cam.position).sub(controls.target).normalize()
          .multiplyScalar(6);
    gizCam.up.copy(cam.up); gizCam.lookAt(0, 0, 0);
    renderer.render(gizScene, gizCam);
    renderer.setScissorTest(false);
    renderer.autoClear = true;
  }

  // A small hook for automation: force a redraw and read the state without waiting on
  // the animation loop. Used to script screenshots, and to test the viewer.
  wrap.cadAPI = {
    draw: () => { controls.update(); render(); },
    png: () => { controls.update(); render();
                 return renderer.domElement.toDataURL('image/png'); },
    parts: () => parts.map(p => ({name: p.pretty, visible: p.mesh.visible,
                                  context: !!p.isGhost})),
    measure: () => lastMeasure,
    // What the viewer believes the deck constrains, and what it is drawing for each.
    // The probe asserts these against the .inp, so a glyph that does not correspond to
    // a keyword in the deck is a test failure rather than a decoration.
    bc: () => ({
      has_analysis: !!BCSPEC.has_analysis,
      notes: BCSPEC.notes || [],
      items: (BCSPEC.items || []).map((it, i) => ({
        kind: it.kind, keyword: it.keyword, set: it.set, label: it.label,
        detail: it.detail,
        drawn: !!(bcByKind[it.kind] || []).length,
        visible: ((bcByKind[it.kind] || [])[0] || {}).visible !== false})),
      engaging_drawn: (META.grains || []).filter(x => x.engage).length,
      engage_highlight: !!(engageMesh && engageMesh.visible),
      // The drawn arrow's own sense, measured: the moment of the arrowhead's direction
      // about the wheel axis. Its sign must match the sign of VR3 in the deck.
      rotation_drawn_sign: rotCheck
        ? Math.sign(new THREE.Vector3().crossVectors(rotCheck.at, rotCheck.tang)
                    .dot(rotCheck.axisLocal))
        : 0,
      rotation_spec_sign: rotCheck ? Math.sign(rotCheck.sign) : 0}),
    // The edit surface, so a headless test can drive an edit and read back exactly
    // what would be handed to Python -- no kernel required.
    settings: () => settingsNow(),
    edited: () => Object.assign({}, EDITED),
    setSetting: (k, v) => {
      const el = document.getElementById('cadedit_' + k);
      if (!el) return false;
      if (el.type === 'checkbox') el.checked = !!v; else el.value = v;
      el.onchange();
      return true;
    },
    applyState: () => applyState,
    clickApply: () => { const b = document.getElementById('cadapply');
                        if (b) b.click(); },
    reset: () => { const b = document.getElementById('cadreset');
                   if (b) b.click(); },
    blockMatrix: () => (wpEdit ? wpEdit.matrix.elements.slice() : null),
    // The guardrail and drag surface, so a headless test can exercise both with no
    // kernel: the window the viewer believes in, what it would refuse, and a drag
    // expressed in the same numbers a mouse would have produced.
    guard: () => ({window: windowNow(), warnings: fitWarnings(),
                   stale: seatingStale(),
                   armed: dragArm, nudge_um: nudgeUm,
                   in_force: (EDIT && EDIT.settings._sector_mode === 'angle')
                     ? 'sector_deg' : 'arc_length_mm',
                   resolved_deg: EDIT ? EDIT.sector_resolved_deg : null,
                   arc_mm: EDIT ? EDIT.arc_length_mm : null}),
    arm: (on) => setDragArmed(on),
    nudge: (axis, sign, n) => { for (let i = 0; i < (n || 1); i++) nudge(axis, sign); },
    dragAlongArc: (um) => {
      // The same path a mouse takes: an arc displacement turned into an angle by the
      // supplied radius, so the test drives what a user drives.
      const s = settingsNow();
      const dth = (um / 1000.0 / (EDIT.ground_radius_mm || 1)) * 180.0 / Math.PI;
      setField('wp_position_deg', (s.wp_position_deg || 0) + dth);
      return settingsNow().wp_position_deg;
    },
    paramBlock: () => {
      const t = document.getElementById('cadeditjson');
      document.getElementById('cadcopy').click();
      const v = t.value;
      const i = v.indexOf('# paste these');
      return i < 0 ? '' : v.slice(v.indexOf('\n', i) + 1);
    },
    // The Apply reply, so a test can hand the viewer the mimetype bundles a real Colab
    // runtime produces -- the one path that had no test, and the one that was broken.
    readReply: (res) => readReply(res),
    transportName: () => { const t = transport(); return t ? t.name : null; },
    applyText: () => {
      const b = document.getElementById('cadeditstate');
      return b ? (b.textContent || '') : '';
    },
    meta: META};

  (function loop() {
    requestAnimationFrame(loop);
    if (spin && root) root.rotation.y += 0.0035;
    controls.update();
    render();
  })();
} catch (err) {
  STAT.textContent = '';
  PICK.innerHTML = '<b style="color:#b00">viewer failed to start</b>' +
                   '<pre style="white-space:pre-wrap;font-size:10px">' + err + '</pre>';
}
</script>
"""


def build(plan: dict, glb_path: str, mode: str = "contact", max_grits: int = 0,
          window_um: float = 0.0, height: int = 720,
          max_inline_mb: float = 24.0) -> tuple:
    """Write the .glb and return (html, meta, glb_info)."""
    from .glb import parts_from_plan, write_glb

    parts, meta = parts_from_plan(plan, mode=mode, max_grits=max_grits,
                                  window_um=window_um, with_meta=True,
                                  budget_mb=max_inline_mb)
    info = write_glb(glb_path, parts)
    return viewer_html(glb_path, meta, height=height,
                       max_inline_mb=max_inline_mb), meta, info
