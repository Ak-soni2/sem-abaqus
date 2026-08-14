"""Write the model as a glTF binary (.glb), and view it with no account and no key.

Why this rather than a cloud viewer
-----------------------------------
Autodesk's viewer is billed and needs the geometry uploaded and translated. Sketchfab
and the rest are the same trade with a smaller bill. None of that buys anything here,
because the triangles are already in hand -- the only thing missing was a renderer that
shades them properly.

glTF is the answer: an open standard, and Google's ``<model-viewer>`` web component
renders it from a CDN with physically-based lighting, shadows and orbit controls. No
key, no upload, no translation, and the same file opens in Blender, Windows 3D Viewer,
PowerPoint and Sketchfab if you ever want them.

The writer is deliberately dependency-free. A GLB is a 12-byte header, a JSON chunk and
a binary chunk; writing those directly keeps the view provably identical to the deck
rather than routed through a converter.

Axis convention: the model is Z-up (wheel axis on Z) and glTF is Y-up, so positions are
emitted as (x, z, -y). Distances are unchanged, so the shape and proportions are exact.
"""

from __future__ import annotations

import base64
import json
import math
import struct
from typing import Optional, Sequence

import numpy as np

# glTF component and target constants, spelled out rather than magic numbers.
_FLOAT, _UINT32 = 5126, 5125
_ARRAY_BUFFER, _ELEMENT_ARRAY_BUFFER = 34962, 34963


def _pad4(b: bytes, fill: bytes) -> bytes:
    return b + fill * ((4 - len(b) % 4) % 4)


def write_glb(path: str, parts: Sequence[dict]) -> dict:
    """Write ``parts`` to a .glb.

    Each part is ``{"name", "vertices" (N,3), "faces" (M,3), "color" (r,g,b,a)}``.
    Returns a summary. Colours are base-colour factors; alpha < 1 turns on blending so
    the workpiece can be seen through.
    """
    bin_parts: list[bytes] = []
    views, accessors, meshes, nodes, materials = [], [], [], [], []
    offset = 0

    for p in parts:
        v = np.asarray(p["vertices"], dtype=np.float64)
        f = np.asarray(p["faces"], dtype=np.int64)
        if len(v) == 0 or len(f) == 0:
            continue
        # Drop vertices no triangle refers to. A trimmed part -- the bond rim cut down
        # to the contact window -- keeps the whole rim's node array, and a viewer that
        # frames the camera from POSITION min/max would then fit to geometry it is not
        # drawing. Compacting also shrinks the file.
        used = np.unique(f)
        if len(used) < len(v):
            remap = np.full(len(v), -1, dtype=np.int64)
            remap[used] = np.arange(len(used))
            v = v[used]
            f = remap[f]
        # Z-up -> Y-up
        pos = np.column_stack([v[:, 0], v[:, 2], -v[:, 1]]).astype(np.float32)
        idx = f.astype(np.uint32).ravel()

        pb = pos.tobytes()
        ib = idx.tobytes()
        pb = _pad4(pb, b"\0")
        ib = _pad4(ib, b"\0")

        views.append({"buffer": 0, "byteOffset": offset, "byteLength": len(pb),
                      "target": _ARRAY_BUFFER})
        offset += len(pb)
        views.append({"buffer": 0, "byteOffset": offset, "byteLength": len(ib),
                      "target": _ELEMENT_ARRAY_BUFFER})
        offset += len(ib)
        bin_parts += [pb, ib]

        accessors.append({"bufferView": len(views) - 2, "componentType": _FLOAT,
                          "count": int(len(pos)), "type": "VEC3",
                          # min/max on POSITION is required by the spec and is what
                          # lets a viewer frame the camera without scanning the buffer
                          "min": [float(x) for x in pos.min(axis=0)],
                          "max": [float(x) for x in pos.max(axis=0)]})
        accessors.append({"bufferView": len(views) - 1, "componentType": _UINT32,
                          "count": int(len(idx)), "type": "SCALAR"})

        r, g, b, a = p.get("color", (0.7, 0.7, 0.7, 1.0))
        materials.append({
            "name": p["name"],
            "pbrMetallicRoughness": {"baseColorFactor": [r, g, b, a],
                                     "metallicFactor": p.get("metallic", 0.15),
                                     "roughnessFactor": p.get("roughness", 0.6)},
            "doubleSided": True,
            **({"alphaMode": "BLEND"} if a < 1.0 else {})})
        meshes.append({"name": p["name"],
                       "primitives": [{"attributes": {"POSITION": len(accessors) - 2},
                                       "indices": len(accessors) - 1,
                                       "material": len(materials) - 1}]})
        nodes.append({"mesh": len(meshes) - 1, "name": p["name"]})

    if not meshes:
        raise ValueError("nothing to write: every part was empty")

    blob = b"".join(bin_parts)
    gltf = {"asset": {"version": "2.0", "generator": "semgrit"},
            "scene": 0, "scenes": [{"nodes": list(range(len(nodes)))}],
            "nodes": nodes, "meshes": meshes, "materials": materials,
            "accessors": accessors, "bufferViews": views,
            "buffers": [{"byteLength": len(blob)}]}

    js = _pad4(json.dumps(gltf, separators=(",", ":")).encode("utf-8"), b" ")
    bn = _pad4(blob, b"\0")
    total = 12 + 8 + len(js) + 8 + len(bn)
    with open(path, "wb") as fh:
        fh.write(struct.pack("<III", 0x46546C67, 2, total))     # 'glTF', version 2
        fh.write(struct.pack("<II", len(js), 0x4E4F534A))       # 'JSON'
        fh.write(js)
        fh.write(struct.pack("<II", len(bn), 0x004E4942))       # 'BIN'
        fh.write(bn)
    return {"path": path, "bytes": total, "parts": len(meshes),
            "triangles": sum(m["primitives"][0]["indices"] is not None
                             and accessors[m["primitives"][0]["indices"]]["count"] // 3
                             for m in meshes)}


MODES = ("contact", "wheel", "whole wheel")

# A face is 3 uint32 indices, a vertex 3 float32 positions: 12 bytes either way.
_B_PER_VERTEX = _B_PER_FACE = 12
FAR_PART = "abrasive grits (far, simplified)"
GHOST_PART = "whole wheel (context, not in the deck)"
MARK_PART = "contact marker (not in the deck)"

# 6 quads of a box, wound so the right-hand normal points out.
_BOX_QUADS = np.array([(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1),
                       (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)])
_BOX_TRIS = np.vstack([_BOX_QUADS[:, [0, 1, 2]], _BOX_QUADS[:, [0, 2, 3]]])


def _box_from_bounds(lo, hi, basis) -> np.ndarray:
    """The 8 corners of an axis-aligned box in the contact frame, back in world axes."""
    c = np.array([[lo[0], lo[1], lo[2]], [hi[0], lo[1], lo[2]],
                  [hi[0], hi[1], lo[2]], [lo[0], hi[1], lo[2]],
                  [lo[0], lo[1], hi[2]], [hi[0], lo[1], hi[2]],
                  [hi[0], hi[1], hi[2]], [lo[0], hi[1], hi[2]]])
    return c @ basis.T


def _ghost_wheel(r_in: float, r_out: float, width: float,
                 n: int = 180) -> tuple[np.ndarray, np.ndarray]:
    """A closed annular disc: the complete wheel the modelled sector was cut from.

    Drawn only as context. At a 1 mm arc on a 50 mm wheel the sagitta is 5 um, so the
    slice on its own genuinely looks flat -- the curvature only becomes visible when
    you can see the wheel it came from.
    """
    th = np.linspace(0.0, 2 * math.pi, n, endpoint=False)
    cs, sn = np.cos(th), np.sin(th)
    hz = width / 2.0
    # node order: (ring, theta, z) with ring 0 = bore, 1 = outer
    v = np.empty((2 * n * 2, 3), dtype=np.float64)
    for ri, r in enumerate((r_in, r_out)):
        for zi, z in enumerate((-hz, hz)):
            k = (ri * n + np.arange(n)) * 2 + zi
            v[k, 0] = r * cs
            v[k, 1] = r * sn
            v[k, 2] = z

    def nid(ri, ti, zi):
        return (ri * n + (ti % n)) * 2 + zi

    q = []
    for t in range(n):
        q.append((nid(1, t, 0), nid(1, t + 1, 0), nid(1, t + 1, 1), nid(1, t, 1)))
        if r_in > 0:                       # a solid disc has no bore to draw
            q.append((nid(0, t, 0), nid(0, t, 1), nid(0, t + 1, 1), nid(0, t + 1, 0)))
        q.append((nid(0, t, 0), nid(0, t + 1, 0), nid(1, t + 1, 0), nid(1, t, 0)))
        q.append((nid(0, t, 1), nid(1, t, 1), nid(1, t + 1, 1), nid(0, t + 1, 1)))
    quads = np.asarray(q, dtype=np.int64)
    return v, np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])


def _cost(parts) -> int:
    """Encoded byte size, counting only vertices some face refers to."""
    n = 0
    for p in parts:
        f = np.asarray(p["faces"])
        n += _B_PER_FACE * len(f) + _B_PER_VERTEX * len(np.unique(f))
    return n


def parts_from_plan(plan: dict, mode: str = "contact", max_grits: int = 0,
                    window_um: float = 0.0, with_meta: bool = False,
                    budget_mb: float = 0.0):
    """Bond, grits and workpiece as glTF parts, from the same geometry the deck uses.

    ``mode`` is ``contact`` (the patch under the block), ``wheel`` (the modelled
    sector) or ``whole wheel`` (the sector plus a translucent ghost of the complete
    wheel, so the curvature is visible).

    ``max_grits = 0`` draws every grain, which is the point -- the old default of 400
    was a guess, and a 712-grain wheel is only a 1.6 MB glTF. ``budget_mb`` caps the
    encoded size: past it, grains far from the contact are replaced by a 12-triangle
    box of the same size and orientation rather than dropped, and only if that is
    still too big are any left out. Either way the meta records what happened; a
    viewer that quietly showed a third of the wheel would be worse than one that
    refused.

    With ``with_meta`` also returns per-grain records -- id, triangle range in the
    merged grit mesh, protrusion, volume -- so a viewer can identify the grain under
    the cursor. The ranges are built in the same pass as the mesh, because a second
    pass that re-derived the ordering could drift out of step with it.
    """
    from .rigid_wheel import build_rim_shell

    if mode not in MODES:
        raise ValueError("mode must be one of %s, not %r" % (", ".join(MODES), mode))

    model, place, wp = plan["_model"], plan["_place"], plan["_wp"]
    R = plan["outer_radius_mm"]
    e_r, e_t, e_z = place["e_r"], place["e_t"], place["e_z"]
    basis = place["basis"]
    r_ground = plan["ground_radius_mm"] or R
    if window_um <= 0:
        window_um = (max(wp.length_mm, wp.width_mm) * 1000.0 * 1.8
                     if wp is not None else 200.0)
    half = window_um / 2000.0
    centre = e_r * r_ground
    out: list[dict] = []
    notes: list[str] = []

    nodes, quads, _ = build_rim_shell(model.spec)
    if mode == "contact":
        # Trim by whole quads, not by triangles. Cutting the two halves of a quad
        # independently leaves a saw-tooth boundary, because their centroids sit on
        # opposite sides of the window edge.
        c = nodes[quads].mean(axis=1) - centre
        quads = quads[(np.abs(c @ e_t) <= half) & (np.abs(c @ e_z) <= half)
                      & (np.abs(c @ e_r) <= max(half, plan["rim_depth_mm"]))]
    tris = np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]]) if len(quads) \
        else np.zeros((0, 3), dtype=np.int64)
    if len(tris):
        out.append({"name": "bond rim", "vertices": nodes, "faces": tris,
                    "color": (0.72, 0.74, 0.77, 1.0), "metallic": 0.25,
                    "roughness": 0.55})

    ghost = mode == "whole wheel" and not plan.get("full_wheel", False)
    if ghost:
        # The whole disc, not the rim ring: the modelled sector's inner radius is
        # only 12 um in from the surface, so a ghost drawn between those two radii is
        # a hairline circle and shows no wheel at all. The bore is unknown -- the SEM
        # never saw it -- so draw it solid rather than invent a hub diameter.
        gv, gf = _ghost_wheel(0.0, R, model.spec.width_mm)
        out.append({"name": GHOST_PART, "vertices": gv, "faces": gf,
                    "color": (0.55, 0.60, 0.66, 0.38), "metallic": 0.0,
                    "roughness": 0.9})
        # At wheel scale the 30 um rim and its grains are far below one pixel, so
        # without a marker the view is a blank disc and you cannot tell where to zoom.
        # Deliberately oversized -- it is a pointer, not geometry, and it says so.
        # Sized against the wheel, never against the view window: with a millimetre
        # workpiece that window is a few mm across, and scaling the pointer to it made
        # the marker 4 mm wide -- larger than everything it was meant to point at.
        mr = R * 0.012
        out.append({"name": MARK_PART,
                    "vertices": _box_from_bounds(
                        np.array([r_ground - mr, -mr, -mr]),
                        np.array([r_ground + mr, mr, mr]), basis),
                    "faces": _BOX_TRIS, "color": (0.95, 0.45, 0.10, 0.85),
                    "metallic": 0.0, "roughness": 0.5})

    baked, faces, frames = place["baked"], place["faces"], place["frames"]
    order = list(range(len(baked)))
    if mode == "contact":
        order = [i for i in order if abs(frames[i][:, 1].mean()) <= half
                 and abs(frames[i][:, 2].mean()) <= half]
    in_view = len(order)
    # Nearest the block first, so whatever has to be degraded is the furthest away.
    order.sort(key=lambda i: abs(frames[i][:, 1].mean()))
    if max_grits and len(order) > max_grits:
        notes.append("capped at max_grits=%d of %d grains in view; set it to 0 for all"
                     % (max_grits, len(order)))
        order = order[:max_grits]

    # How many can be drawn in full before the budget bites. Everything after that
    # index becomes a box; the split is by distance from the contact, so detail is
    # kept exactly where it is looked at.
    n_full = len(order)
    if budget_mb > 0 and order:
        fixed = _cost(out)
        per = [_B_PER_FACE * len(faces[i]) + _B_PER_VERTEX * len(baked[i])
               for i in order]
        box_cost = _B_PER_FACE * len(_BOX_TRIS) + _B_PER_VERTEX * 8
        allow = budget_mb * 1e6 * 0.75 - fixed        # base64 inflates by 4/3
        total = sum(per) + (0 if wp is None else 12 * (8 + 12))
        if total > allow:
            run, n_full = 0.0, 0
            for k, cst in enumerate(per):
                run += cst - box_cost
                if run + box_cost * len(order) > allow:
                    break
                n_full = k + 1
            notes.append(
                "%d of %d grains drawn in full detail; the %d furthest from the "
                "contact are drawn as boxes to stay inside the %g MB budget"
                % (n_full, len(order), len(order) - n_full, budget_mb))

    # The deck's own engaging set, from the shared helper the writer uses -- not a
    # lookalike recomputed here from a different window.
    _eng = set(plan.get("_engage") or [])
    meta, meta_far = [], []
    if order:
        for lo, hi, dest, name, colour in (
                (0, n_full, meta, "abrasive grits", (0.16, 0.62, 0.36, 1.0)),
                (n_full, len(order), meta_far, FAR_PART, (0.42, 0.58, 0.47, 1.0))):
            if hi <= lo:
                continue
            V, F, off, tri = [], [], 0, 0
            for i in order[lo:hi]:
                pl_i = model.placements[i]
                sh = model.shapes[pl_i.shape_index]
                if dest is meta:
                    v, f = baked[i], np.asarray(faces[i])
                else:
                    fr = frames[i]
                    v = _box_from_bounds(fr.min(axis=0), fr.max(axis=0), basis)
                    f = _BOX_TRIS
                V.append(v)
                F.append(f + off)
                dest.append({
                    "id": int(pl_i.placement_id),
                    "tri0": tri, "ntri": int(len(f)),
                    "proxy": dest is meta_far,
                    "engage": bool(i in _eng),
                    "protrusion_um": round(float(np.hypot(
                        baked[i][:, 0], baked[i][:, 1]).max() - R) * 1000, 4),
                    "volume_um3": round(float(sh.mesh_volume_um3), 3),
                    "height_um": round(float(sh.height_um), 3),
                    "width_um": round(float(max(sh.extent_um()[:2])), 3),
                    "b_um": round(float(frames[i][:, 1].mean()) * 1000, 3),
                    "z_um": round(float(frames[i][:, 2].mean()) * 1000, 3)})
                off += len(v)
                tri += len(f)
            out.append({"name": name, "vertices": np.vstack(V),
                        "faces": np.vstack(F), "color": colour,
                        "metallic": 0.1, "roughness": 0.45})

    if wp is not None:
        hb, hz, d = wp.length_mm / 2.0, wp.width_mm / 2.0, wp.depth_mm
        loc = np.array([[r_ground, -hb, -hz], [r_ground, hb, -hz],
                        [r_ground, hb, hz], [r_ground, -hb, hz],
                        [r_ground + d, -hb, -hz], [r_ground + d, hb, -hz],
                        [r_ground + d, hb, hz], [r_ground + d, -hb, hz]])
        corners = loc[:, 0:1] * e_r + loc[:, 1:2] * e_t + loc[:, 2:3] * e_z
        quad = np.array([(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1),
                         (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)])
        out.append({"name": "workpiece", "vertices": corners,
                    "faces": np.vstack([quad[:, [0, 1, 2]], quad[:, [0, 2, 3]]]),
                    "color": (0.18, 0.44, 0.71, 0.55), "metallic": 0.0,
                    "roughness": 0.8})
    if with_meta:
        # e_r points from the axis out through the contact, so looking back along it
        # frames the dressed face -- the view a grinding paper actually shows. Map it
        # into the viewer's Y-up axes the same way vertices are mapped.
        from .bcspec import to_viewer as to_view
        face, axis, arc = to_view(e_r), to_view(e_z), to_view(e_t)
        ctr = centre
        from .bcspec import build as _bc_build
        _edit = None
        if plan.get("_params") is not None:
            import dataclasses as _dc
            from .editable import FIELDS as _EF
            from .editable import settings_from_params as _sfp
            _edit = {"settings": _sfp(plan["_params"]),
                     "fields": [_dc.asdict(f) for f in _EF],
                     # What the viewer may preview by transforming the drawn box, and
                     # what it must hand back to Python instead. Supplied, not inferred.
                     "basis": {"radial": face, "arc": arc, "axial": axis},
                     "block": ({"centre": [float(x) for x in
                                           to_view(e_r * (r_ground + wp.depth_mm / 2.0))],
                                "half": [wp.depth_mm / 2.0, wp.length_mm / 2.0,
                                         wp.width_mm / 2.0]}
                               if wp is not None else None),
                     "ground_radius_mm": r_ground,
                     # The window a depth of cut has to land in. Both ends are already
                     # computed by plan_deck; without them the viewer edits the one
                     # number that has twice produced a job which ran and ground
                     # nothing, with nothing on screen to say so.
                     "first_contact_um": plan.get("first_contact_um"),
                     "depth_ceiling_um": plan.get("depth_ceiling_um"),
                     "engaging_now": len(plan.get("_engage") or []),
                     "n_grits": plan.get("n_grits"),
                     "arc_length_mm": plan.get("arc_length_mm"),
                     "sector_resolved_deg": plan.get("sector_deg")}
        return out, {"grains": meta, "grains_far": meta_far, "mode": mode,
                     "bc": _bc_build(plan), "edit": _edit,
                     "face_dir": face, "axis_dir": axis, "arc_dir": arc,
                     # where to point the Contact button, in the viewer's Y-up axes
                     "contact_centre": [float(ctr[0]), float(ctr[2]), float(-ctr[1])],
                     "contact_radius_mm": max(window_um / 2000.0, 1e-4),
                     "wheel_radius_mm": R, "ghost": ghost,
                     "mark_size_mm": R * 0.012,
                     "far_part": FAR_PART, "ghost_part": GHOST_PART,
                     "mark_part": MARK_PART,
                     "sector_deg": plan.get("sector_deg", 0.0),
                     "grits_total": len(baked), "grits_in_view": in_view,
                     "grits_drawn": len(order), "grits_full_detail": n_full,
                     "notes": notes,
                     "outer_radius_mm": R,
                     "ground_radius_mm": plan["ground_radius_mm"],
                     "bond_clearance_um": plan["bond_clearance_um"],
                     "depth_of_cut_um": plan.get("depth_of_cut_um") or 0.0,
                     "workpiece_mm": ([wp.length_mm, wp.width_mm, wp.depth_mm]
                                      if wp is not None else None)}
    return out


def model_viewer_html(glb_path: str, height: int = 620,
                      max_inline_mb: float = 12.0) -> str:
    """Embed the .glb with Google's <model-viewer>: free, keyless, CDN-hosted.

    The file goes in as a data URI so nothing has to be served or uploaded. Base64
    inflates by a third, hence the cap -- past it the notebook itself becomes unwieldy.
    """
    import os

    mb = os.path.getsize(glb_path) / 1e6
    if mb > max_inline_mb:
        raise ValueError(
            "%.1f MB .glb is too large to inline (cap %.1f MB). Use mode='contact' or "
            "lower max_grits, or just download the .glb and open it in Blender or "
            "Windows 3D Viewer." % (mb, max_inline_mb))
    uri = ("data:model/gltf-binary;base64,"
           + base64.b64encode(open(glb_path, "rb").read()).decode())
    return """
<script type="module"
  src="https://ajax.googleapis.com/ajax/libs/model-viewer/3.5.0/model-viewer.min.js">
</script>
<model-viewer src="%s" alt="grinding wheel"
  camera-controls touch-action="pan-y" auto-rotate shadow-intensity="1"
  exposure="1.1" environment-image="neutral" ar
  style="width:100%%;height:%dpx;background:#eef1f4;border-radius:6px">
</model-viewer>
<div style="font:12px monospace;color:#555;padding-top:4px">
  drag to orbit &middot; scroll to zoom &middot; %.1f MB glTF, rendered locally --
  nothing uploaded
</div>
""" % (uri, height, mb)
