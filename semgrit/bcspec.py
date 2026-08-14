"""What the deck constrains, described as numbers a viewer can draw.

Why this is a Python module and not JavaScript
----------------------------------------------
The viewer is proven to draw the deck's own triangles -- ``verify_colab`` asserts the
glTF holds exactly ``{bond rim, abrasive grits, workpiece}`` and nothing else. So
boundary-condition symbols cannot be shipped as extra meshes without breaking that
proof. They travel instead as *numbers*: an anchor point, a unit direction, a magnitude,
a set name and the set's true size. The browser places unit shapes at them and does no
geometry of its own.

That division is not fastidiousness. Twice on this project two implementations of one
idea have agreed on the same wrong answer -- the infeed sign, and a deck header that
described the rotation backwards for months. Every direction here is therefore derived
once, from the same basis the writer uses, and mapped Z-up to Y-up by the same single
helper the vertices use.

Nothing is invented: a glyph appears only if the deck actually writes the keyword. A
geometry-only deck reports no boundary conditions at all, because it has none.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np


def to_viewer(v) -> list:
    """Z-up (deck) to Y-up (glTF): (x, z, -y). The one place this sign lives.

    ``glb.write_glb`` applies the same map to vertices. Keeping directions and anchors
    on this single helper is what stops a hand-written second copy from disagreeing with
    it -- which is exactly how the last two sign bugs happened.
    """
    a = np.asarray(v, dtype=float)
    return [float(a[0]), float(a[2]), float(-a[1])]


def _quad(loc, e_r, e_t, e_z) -> list:
    """Four corners given in (radial, tangential, axial), returned in viewer axes."""
    out = []
    for a, b, c in loc:
        p = a * np.asarray(e_r) + b * np.asarray(e_t) + c * np.asarray(e_z)
        out.append(to_viewer(p))
    return out


def _sample(quad, n_side: int = 4) -> list:
    """A small grid of anchor points across a quad, for the symbols drawn on a face.

    A constrained face carries tens of thousands of nodes. Drawing one symbol each is
    both unreadable and slow, so a handful are sampled -- and the caller always pairs
    them with the full face patch and the set's true node count, so the decimation can
    never be mistaken for a smaller boundary condition.
    """
    p0, p1, p2, p3 = (np.asarray(q, dtype=float) for q in quad)
    pts = []
    for i in range(n_side):
        for j in range(n_side):
            u = (i + 0.5) / n_side
            v = (j + 0.5) / n_side
            top = p0 + (p1 - p0) * u
            bot = p3 + (p2 - p3) * u
            pts.append([float(x) for x in (top + (bot - top) * v)])
    return pts


def build(plan: dict) -> dict:
    """Every boundary condition, load and interaction the deck writes, as drawable data.

    Returns ``{"has_analysis", "items": [...], "notes": [...]}``. Each item carries a
    ``kind`` the viewer knows how to draw, the Abaqus keyword it stands for, the set it
    applies to, that set's real size, and the numbers to label it with.
    """
    p = plan.get("_params")
    wp = plan.get("_wp")
    place = plan.get("_place") or {}
    an = getattr(p, "analysis", None) if p is not None else None
    cost = plan.get("cost") or {}
    R = plan["outer_radius_mm"]
    items: list[dict] = []
    notes: list[str] = []

    if an is None or not getattr(an, "enabled", False):
        return {"has_analysis": False, "items": [],
                "notes": ["This deck is geometry only: it contains no step, no boundary "
                          "conditions and no contact. Nothing is drawn because there is "
                          "nothing to draw. Enable the run-ready analysis to see them."]}

    e_r, e_t, e_z = place["e_r"], place["e_t"], place["e_z"]
    r_ground = plan["ground_radius_mm"] or R

    # ---- the workpiece faces that are held ------------------------------
    if wp is not None:
        hb, hz, d = wp.length_mm / 2.0, wp.width_mm / 2.0, wp.depth_mm
        nl, nw, nd = wp.divisions()
        # (name, on?, the four corners in (radial, tangential, axial), node count)
        faces = [
            ("A_WP_BACK_FACE", an.fix_back_face,
             [(r_ground + d, -hb, -hz), (r_ground + d, hb, -hz),
              (r_ground + d, hb, hz), (r_ground + d, -hb, hz)],
             (nl + 1) * (nw + 1), "the face away from the wheel"),
            ("A_WP_END_A", an.fix_ends,
             [(r_ground, -hb, -hz), (r_ground + d, -hb, -hz),
              (r_ground + d, -hb, hz), (r_ground, -hb, hz)],
             (nw + 1) * (nd + 1), "the exit end"),
            ("A_WP_END_B", an.fix_ends,
             [(r_ground, hb, -hz), (r_ground + d, hb, -hz),
              (r_ground + d, hb, hz), (r_ground, hb, hz)],
             (nw + 1) * (nd + 1), "the entry end"),
            ("A_WP_SIDE_A", an.fix_sides,
             [(r_ground, -hb, -hz), (r_ground + d, -hb, -hz),
              (r_ground + d, hb, -hz), (r_ground, hb, -hz)],
             (nl + 1) * (nd + 1), "one side of the face"),
            ("A_WP_SIDE_B", an.fix_sides,
             [(r_ground, -hb, hz), (r_ground + d, -hb, hz),
              (r_ground + d, hb, hz), (r_ground, hb, hz)],
             (nl + 1) * (nd + 1), "the other side"),
        ]
        for name, on, loc, n_nodes, what in faces:
            if not on:
                continue
            q = _quad(loc, e_r, e_t, e_z)
            items.append({
                "kind": "encastre", "keyword": "*Boundary", "value": "ENCASTRE",
                "set": name, "nodes": int(n_nodes), "what": what,
                "quad": q, "anchors": _sample(q),
                "label": "ENCASTRE",
                "detail": "%s, %s nodes, all 6 dof held from t=0"
                          % (name, format(int(n_nodes), ",")),
            })
        if not any(on for _, on, _, _, _ in faces):
            notes.append("no workpiece face is held: the block is free to fly away.")

        # ---- the contact pair -------------------------------------------
        if an.contact_scope != "none":
            eng = plan.get("_engage") or []
            master = ("A_GRITS_ENGAGE_SURF"
                      if an.contact_scope == "engaging" and eng
                      and len(eng) < len(place.get("frames") or []) else "A_GRITS_SURF")
            if an.contact_scope == "all exterior":
                master = "ALL EXTERIOR"
            gq = _quad([(r_ground, -hb, -hz), (r_ground, hb, -hz),
                        (r_ground, hb, hz), (r_ground, -hb, hz)], e_r, e_t, e_z)
            items.append({
                "kind": "contact", "keyword": "*Contact Inclusions",
                "set": master, "slave": "A_WP_GROUND_SURF",
                "quad": gq, "anchors": [],
                "n_engaging": int(len(eng)),
                "n_grits_total": int(len(place.get("frames") or [])),
                "friction": float(an.friction),
                "label": "CONTACT",
                "detail": "%s to A_WP_GROUND_SURF, friction %g%s"
                          % (master, an.friction,
                             ", %s grits able to engage" % format(len(eng), ",")
                             if an.contact_scope == "engaging" else ""),
            })

    # ---- the wheel: one velocity BC on one reference node ----------------
    omega = cost.get("surface_speed_mm_s", p.surface_speed_mm_s) / R if R else 0.0
    ae_mm = (plan.get("depth_of_cut_um") or 0.0) / 1000.0
    t_step = float(cost.get("step_time_s") or 0.0)
    v_r = (ae_mm / t_step) if t_step > 0 else 0.0
    contact_pt = np.asarray(e_r) * r_ground

    items.append({
        "kind": "refnode", "keyword": "*Rigid Body", "set": "A_WHEEL_REF",
        "at": to_viewer((0.0, 0.0, 0.0)),
        "leader": to_viewer(contact_pt),
        "label": "RP",
        "detail": "one rigid body, elset ES_WHEEL_ALL: the bond rim and every grit "
                  "facet, driven by this single node on the axis",
    })
    if v_r > 0:
        items.append({
            "kind": "velocity", "keyword": "*Boundary, type=VELOCITY",
            "set": "A_WHEEL_REF",
            # +e_r: the block sits OUTSIDE the rim, so feeding in means translating the
            # wheel outward. Drawn from the same vector the BC is written from.
            "at": to_viewer(contact_pt), "dir": to_viewer(e_r),
            "magnitude": float(v_r), "unit": "mm/s",
            "label": "INFEED V1,V2",
            "detail": "%.1f mm/s along +e_r = %.4f um of depth of cut over the step"
                      % (v_r, ae_mm * 1000.0),
        })
    else:
        notes.append("the depth of cut is zero, so no infeed is prescribed and the "
                     "wheel will spin without cutting.")
    items.append({
        "kind": "rotation", "keyword": "*Boundary, type=VELOCITY", "set": "A_WHEEL_REF",
        "at": to_viewer((0.0, 0.0, 0.0)),
        "axis": to_viewer(e_z),
        # VR3 = -omega. A positive rotation about +Z carries +X toward +Y, so a negative
        # one drives the surface toward DECREASING theta: grits arrive at the block from
        # its high-theta end. The deck header said the opposite for months.
        "sign": -1.0, "magnitude": float(omega), "unit": "rad/s",
        "rpm": float(omega * 30.0 / math.pi),
        "radius_mm": float(R),
        "surface_speed_mm_s": float(cost.get("surface_speed_mm_s")
                                    or p.surface_speed_mm_s),
        "label": "VR3",
        "detail": "VR3 = -%.1f rad/s (%.0f rpm): the surface travels toward decreasing "
                  "theta, so grits arrive from the high-theta end"
                  % (omega, omega * 30.0 / math.pi),
    })

    held = [it for it in items if it["kind"] == "encastre"]
    return {
        "has_analysis": True,
        "items": items,
        "held_nodes_total": int(sum(it["nodes"] for it in held)),
        "notes": notes,
    }
