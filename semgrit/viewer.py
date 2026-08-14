"""Interactive 3-D view of the model, in the notebook output.

The matplotlib preview answers "are the numbers right". This answers "is that the
object I meant" -- you can orbit it, zoom into the contact and see the actual measured
grains sitting on the bond with the block resting on them.

Built from the same geometry the writer emits: the rim shell quads from
``build_rim_shell``, the baked grit facets from ``place_workpiece``, and the workpiece
block corners. Nothing is re-derived, so the view cannot show a different wheel from
the deck.

Two things force the design:

**Scale.** The wheel is 50 mm and a grit is 3 um -- 4 orders apart. Rendered whole, the
grits are far below one pixel. So ``mode='contact'`` clips to a window around the block,
which is the only zoom at which grains are visible at all, and ``mode='wheel'`` shows
the whole sector for proportion.

**Triangle budget.** A dressed wheel can carry half a million facets; a browser starts
to struggle past ~100k. Grits are therefore capped, nearest the block first, and the
count actually drawn is reported rather than silently truncated.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np


def _quads_to_tris(quads: np.ndarray) -> np.ndarray:
    q = np.asarray(quads)
    return np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])


def _box_tris(corners: np.ndarray) -> np.ndarray:
    """12 triangles of a hexahedron given its 8 corners in C3D8 order."""
    faces = [(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1),
             (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)]
    return _quads_to_tris(np.array(faces))


def view3d(plan: dict, mode: str = "contact", max_grits: int = 400,
           window_um: float = 0.0, show_bond: bool = True,
           show_workpiece: bool = True):
    """Interactive Plotly figure of the wheel, its grits and the workpiece.

    ``mode='contact'`` clips to a window around the block; ``'wheel'`` shows the whole
    sector. Returns the figure and a dict of what was actually drawn.
    """
    import plotly.graph_objects as go

    from .rigid_wheel import build_rim_shell

    model = plan["_model"]
    place = plan["_place"]
    wp = plan["_wp"]
    R = plan["outer_radius_mm"]
    e_r, e_t, e_z = place["e_r"], place["e_t"], place["e_z"]
    r_ground = plan["ground_radius_mm"]

    if window_um <= 0:
        window_um = (max(wp.length_mm, wp.width_mm) * 1000.0 * 1.8
                     if wp is not None else 200.0)
    half = window_um / 2000.0                      # mm, half-window
    centre = e_r * (r_ground if r_ground else R)

    data, drawn = [], {}

    # ---------------- bond rim ----------------
    if show_bond:
        nodes, quads, _ = build_rim_shell(model.spec)
        tris = _quads_to_tris(quads)
        if mode == "contact":
            c = nodes[tris].mean(axis=1)
            d = c - centre
            keep = ((np.abs(d @ e_t) <= half) & (np.abs(d @ e_z) <= half)
                    & (np.abs(d @ e_r) <= max(half, plan["rim_depth_mm"])))
            tris = tris[keep]
        drawn["bond_triangles"] = int(len(tris))
        if len(tris):
            data.append(go.Mesh3d(
                x=nodes[:, 0], y=nodes[:, 1], z=nodes[:, 2],
                i=tris[:, 0], j=tris[:, 1], k=tris[:, 2],
                color="#b9bec4", opacity=1.0, flatshading=True,
                name="bond rim", showlegend=True, hoverinfo="name"))

    # ---------------- grits ----------------
    baked, faces, frames = place["baked"], place["faces"], place["frames"]
    order = list(range(len(baked)))
    if mode == "contact":
        order = [i for i in order
                 if abs(frames[i][:, 1].mean()) <= half
                 and abs(frames[i][:, 2].mean()) <= half]
    # nearest the block centre first, so a cap keeps the ones that matter
    order.sort(key=lambda i: abs(frames[i][:, 1].mean()))
    capped = order[:max_grits]
    drawn["grits_drawn"] = len(capped)
    drawn["grits_in_view"] = len(order)
    drawn["grits_total"] = len(baked)

    if capped:
        V, F, C = [], [], []
        off = 0
        for i in capped:
            v = baked[i]
            V.append(v)
            F.append(np.asarray(faces[i]) + off)
            # colour every vertex of a grain by that grain's protrusion
            prot = float(np.hypot(v[:, 0], v[:, 1]).max()) - R
            C.append(np.full(len(v), prot * 1000.0))
            off += len(v)
        V = np.vstack(V)
        F = np.vstack(F)
        C = np.concatenate(C)
        drawn["grit_triangles"] = int(len(F))
        data.append(go.Mesh3d(
            x=V[:, 0], y=V[:, 1], z=V[:, 2],
            i=F[:, 0], j=F[:, 1], k=F[:, 2],
            intensity=C, colorscale="Viridis", flatshading=True,
            colorbar=dict(title=dict(text="grit protrusion<br>(um)", side="right"),
                          thickness=12, len=0.6),
            name="abrasive grits", showlegend=True, hoverinfo="name"))

    # ---------------- workpiece ----------------
    if show_workpiece and wp is not None:
        hb, hz, d = wp.length_mm / 2.0, wp.width_mm / 2.0, wp.depth_mm
        loc = np.array([[r_ground, -hb, -hz], [r_ground, hb, -hz],
                        [r_ground, hb, hz], [r_ground, -hb, hz],
                        [r_ground + d, -hb, -hz], [r_ground + d, hb, -hz],
                        [r_ground + d, hb, hz], [r_ground + d, -hb, hz]])
        corners = loc[:, 0:1] * e_r + loc[:, 1:2] * e_t + loc[:, 2:3] * e_z
        tris = _box_tris(corners)
        data.append(go.Mesh3d(
            x=corners[:, 0], y=corners[:, 1], z=corners[:, 2],
            i=tris[:, 0], j=tris[:, 1], k=tris[:, 2],
            color="#2f6fb5", opacity=0.45, flatshading=True,
            name="workpiece", showlegend=True, hoverinfo="name"))
        drawn["workpiece_mm"] = (wp.length_mm, wp.width_mm, wp.depth_mm)

    title = ("%s -- %s view, %d of %d grits drawn"
             % (plan.get("title", "model"), mode, drawn.get("grits_drawn", 0),
                drawn.get("grits_total", 0)))
    fig = go.Figure(data=data)
    fig.update_layout(
        title=dict(text=title, font=dict(size=13)),
        # 'data' keeps true proportions: anything else silently distorts a 4-order
        # scale range and the picture stops meaning anything.
        scene=dict(aspectmode="data",
                   xaxis_title="x (mm)", yaxis_title="y (mm)", zaxis_title="z (mm)"),
        margin=dict(l=0, r=0, t=34, b=0), height=620,
        legend=dict(orientation="h", yanchor="bottom", y=0.0))
    if mode == "contact":
        lim = [centre - half * 1.2, centre + half * 1.2]
        rng = np.array(lim)
        fig.update_layout(scene=dict(
            aspectmode="data",
            xaxis=dict(range=[rng[:, 0].min(), rng[:, 0].max()], title="x (mm)"),
            yaxis=dict(range=[rng[:, 1].min(), rng[:, 1].max()], title="y (mm)"),
            zaxis=dict(range=[rng[:, 2].min(), rng[:, 2].max()], title="z (mm)")))
    return fig, drawn
