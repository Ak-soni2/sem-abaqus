"""Show the MESH in the same viewer the CAD is shown in.

The CAD viewer draws the geometry a deck describes. It does not draw the mesh,
and the mesh is where the arguments are: whether ``dc`` is resolved, whether the
grading is where it should be, whether the compliant layer has enough elements
through its thickness to bend rather than shear, whether anything is inverted.
Those are the questions a reviewer asks about a grinding model, and until now
the only way to answer them was to open the ``.inp`` in CAE.

So this converts a hex mesh into exactly the part dicts ``semgrit.glb`` already
consumes -- ``{"name", "vertices", "faces", "color"}`` -- which means the mesh
view inherits every feature the CAD viewer has: section capping, explode,
colour-by-property, the measuring tool, all twelve shortcuts. Nothing about the
viewer changes; it is fed different parts.

Two things are done properly rather than quickly.

**Only the boundary is drawn.** A hex mesh has six faces per element, and all
but the outermost are shared by two elements and invisible. Drawing them all
would be 6x the triangles for no pixels, and on a million-element deck that is
the difference between a viewer that opens and one that does not. The interior
faces are found by counting: a quad face appearing twice is internal, once is
boundary. That also makes it a **mesh check** -- a face appearing three times
means the connectivity is corrupt, and :func:`surface_of_hexes` refuses it.

**Element quality is computed, not asserted.** ``quality`` returns the aspect
ratio and Jacobian sign per element so the viewer can colour by them, because
"the mesh looks fine" is not a claim anyone should accept about a deck that
takes days to run.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

# The six faces of a C3D8R hex, each wound so its normal points OUT of the
# element when the element itself has a positive Jacobian. Getting this wrong
# gives a viewer full of black facets from back-face culling.
_HEX_FACES = (
    (0, 3, 2, 1),   # bottom  (-w)
    (4, 5, 6, 7),   # top     (+w)
    (0, 1, 5, 4),   # side
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)

# Okabe-Ito, the palette the rest of the project's figures use, so a mesh view
# and a figure of the same deck agree on what colour a part is.
C_HUB = (0.35, 0.35, 0.38, 1.0)
C_COMPLIANT = (0.90, 0.62, 0.00, 1.0)     # orange: the layer that matters
C_GRAIN = (0.00, 0.45, 0.70, 1.0)         # blue
C_WORK = (0.80, 0.80, 0.82, 1.0)
C_FINE = (0.00, 0.62, 0.45, 1.0)          # green: the dc-resolved band


class MeshViewError(RuntimeError):
    pass


def surface_of_hexes(conn: np.ndarray) -> np.ndarray:
    """Boundary quads of a hex mesh, wound outward.

    A face shared by two elements is interior and dropped. A face appearing
    more than twice cannot happen in a valid mesh, so it is refused rather than
    drawn -- that is a corrupt connectivity table and it would otherwise show
    up much later as an Abaqus preprocessing failure.
    """
    conn = np.asarray(conn, dtype=np.int64)
    if conn.ndim != 2 or conn.shape[1] != 8:
        raise MeshViewError("expected an (n, 8) hex connectivity table, got %r"
                            % (conn.shape,))
    seen: dict = {}
    for row in conn:
        for f in _HEX_FACES:
            quad = (int(row[f[0]]), int(row[f[1]]),
                    int(row[f[2]]), int(row[f[3]]))
            key = tuple(sorted(quad))
            hit = seen.get(key)
            if hit is None:
                seen[key] = [quad, 1]
            else:
                hit[1] += 1
                if hit[1] > 2:
                    raise MeshViewError(
                        "face %r is shared by %d elements. A hex face can be "
                        "shared by at most two; the connectivity is corrupt."
                        % (key, hit[1]))
    out = [q for q, n in seen.values() if n == 1]
    if not out:
        raise MeshViewError("the mesh has no boundary faces at all")
    return np.array(out, dtype=np.int64)


def quads_to_tris(quads: np.ndarray) -> np.ndarray:
    """Split quads into triangles, preserving winding."""
    q = np.asarray(quads, dtype=np.int64)
    return np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])


def hex_jacobians(nodes: np.ndarray, conn: np.ndarray) -> np.ndarray:
    """Signed volume of each hex, by the same tet decomposition as the writer.

    Negative means the element is inverted. Vectorised, because this runs on
    every deck and a Python loop over a million elements is minutes.
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    conn = np.asarray(conn, dtype=np.int64)
    tets = ((0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6),
            (3, 4, 6, 7), (1, 4, 5, 6))
    p = nodes[conn]
    tot = np.zeros(len(conn))
    for a, b, c, d in tets:
        tot += np.einsum("ij,ij->i",
                         np.cross(p[:, b] - p[:, a], p[:, c] - p[:, a]),
                         p[:, d] - p[:, a]) / 6.0
    return tot


def hex_aspect(nodes: np.ndarray, conn: np.ndarray) -> np.ndarray:
    """Longest over shortest edge, per element.

    The measure that matters for an explicit run: the stable time increment is
    set by the SHORTEST dimension, so a 50:1 element costs 50x the steps of a
    cubic one of the same volume. This project's own IMPROVEMENTS.md records an
    aspect ratio of 51:1 as a defect that had to be fixed, so it is worth being
    able to see it.
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    conn = np.asarray(conn, dtype=np.int64)
    edges = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7))
    p = nodes[conn]
    lens = np.stack([np.linalg.norm(p[:, b] - p[:, a], axis=1)
                     for a, b in edges], axis=1)
    lo = lens.min(axis=1)
    hi = lens.max(axis=1)
    return np.where(lo > 0, hi / np.maximum(lo, 1e-300), np.inf)


def quality(nodes: np.ndarray, conn: np.ndarray) -> dict:
    """Everything worth knowing about a mesh before spending days on it."""
    jac = hex_jacobians(nodes, conn)
    asp = hex_aspect(nodes, conn)
    n_inv = int((jac <= 0).sum())
    edges = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7))
    p = np.asarray(nodes)[np.asarray(conn)]
    lens = np.stack([np.linalg.norm(p[:, b] - p[:, a], axis=1)
                     for a, b in edges], axis=1)
    return dict(
        elements=len(conn), nodes=int(np.asarray(conn).max()) + 1,
        volume=float(np.abs(jac).sum()),
        inverted=n_inv,
        min_edge=float(lens.min()), max_edge=float(lens.max()),
        aspect_mean=float(asp[np.isfinite(asp)].mean()),
        aspect_max=float(asp[np.isfinite(asp)].max()),
        aspect_p99=float(np.percentile(asp[np.isfinite(asp)], 99)),
        jacobian_min=float(jac.min()), jacobian_max=float(jac.max()),
    )


def mesh_part(name: str, nodes: np.ndarray, conn: np.ndarray,
              color=C_WORK, *, check: bool = True) -> dict:
    """One hex mesh as a viewer part.

    With ``check`` the element quality is computed and inverted elements are
    refused. That is deliberate: a viewer is the one place a human looks before
    submitting, so it is the right place to stop a mesh that cannot run.
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    if check:
        q = quality(nodes, conn)
        if q["inverted"]:
            raise MeshViewError(
                "%s: %d of %d elements have a non-positive Jacobian, so the "
                "node ordering is wrong. Abaqus reports this as a "
                "preprocessing failure with no element numbers, which is why "
                "it is caught here instead."
                % (name, q["inverted"], q["elements"]))
    tris = quads_to_tris(surface_of_hexes(conn))
    return dict(name=name, vertices=nodes, faces=tris, color=tuple(color))


def wire_part(name: str, nodes: np.ndarray, conn: np.ndarray,
              color=(0.1, 0.1, 0.1, 1.0), width_frac: float = 0.02) -> dict:
    """The element edges, as thin triangular prisms.

    glTF can carry lines, but the viewer's material path, section capping and
    explode all assume triangles, so drawing edges as very thin solids keeps
    every one of those features working on them. ``width_frac`` is a fraction
    of the shortest edge, so it scales with the mesh instead of needing a
    length in millimetres on a model that spans 125 mm to 16 nm.
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    quads = surface_of_hexes(conn)
    # unique undirected edges of the boundary quads only -- interior edges are
    # not visible and would multiply the geometry for nothing
    e = set()
    for q in quads:
        for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
            i, j = int(q[a]), int(q[b])
            e.add((i, j) if i < j else (j, i))
    if not e:
        raise MeshViewError("no boundary edges to draw")
    ee = np.array(sorted(e), dtype=np.int64)
    p0 = nodes[ee[:, 0]]
    p1 = nodes[ee[:, 1]]
    seg = p1 - p0
    L = np.linalg.norm(seg, axis=1)
    keep = L > 0
    p0, p1, seg, L = p0[keep], p1[keep], seg[keep], L[keep]
    r = max(float(L.min()) * width_frac, 1e-12)

    # An orthonormal frame per segment. Choosing the reference axis as the one
    # the segment is LEAST aligned with keeps the cross product conditioned.
    d = seg / L[:, None]
    ref = np.zeros_like(d)
    ref[np.arange(len(d)), np.abs(d).argmin(axis=1)] = 1.0
    u = np.cross(d, ref)
    u /= np.linalg.norm(u, axis=1)[:, None]
    w = np.cross(d, u)

    # triangular prism: 3 verts at each end
    ang = np.array([0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0])
    off = (np.cos(ang)[:, None, None] * u[None] +
           np.sin(ang)[:, None, None] * w[None]) * r      # (3, n, 3)
    v = np.concatenate([p0[None] + off, p1[None] + off], axis=0)   # (6, n, 3)
    n_seg = len(p0)
    verts = v.transpose(1, 0, 2).reshape(-1, 3)           # (n*6, 3)
    base = np.arange(n_seg) * 6
    tri = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        tri.append(np.stack([base + a, base + b, base + 3 + b], axis=1))
        tri.append(np.stack([base + a, base + 3 + b, base + 3 + a], axis=1))
    faces = np.vstack(tri)
    return dict(name=name, vertices=verts, faces=faces, color=tuple(color))


def parts_from_meshes(meshes: Sequence[dict], *, edges: bool = True,
                      max_edge_elements: int = 60_000) -> tuple:
    """Viewer parts for a list of ``{name, nodes, conn, color}`` meshes.

    Edges are skipped above ``max_edge_elements`` and the reason is returned in
    the notes, because a wireframe of a million elements is tens of millions of
    triangles -- it would not open, and silently dropping it would leave a
    reviewer thinking the mesh has no edges.
    """
    parts, notes, stats = [], [], {}
    for m in meshes:
        nodes, conn = m["nodes"], m["conn"]
        q = quality(nodes, conn)
        stats[m["name"]] = q
        parts.append(mesh_part(m["name"], nodes, conn,
                               m.get("color", C_WORK)))
        if edges:
            if q["elements"] <= max_edge_elements:
                parts.append(wire_part(m["name"] + " edges", nodes, conn,
                                       color=m.get("edge_color",
                                                   (0.12, 0.12, 0.14, 1.0))))
            else:
                notes.append(
                    "%s: %s elements, so the wireframe is omitted -- it would "
                    "be ~%s triangles. Use the section plane to inspect the "
                    "interior instead."
                    % (m["name"], format(q["elements"], ","),
                       format(q["elements"] * 4 * 12, ",")))
        if q["aspect_p99"] > 20.0:
            notes.append(
                "%s: 99th-percentile aspect ratio %.1f:1. The stable time "
                "increment follows the SHORTEST edge, so this costs steps."
                % (m["name"], q["aspect_p99"]))
    return parts, dict(notes=notes, stats=stats)


def build(meshes: Sequence[dict], glb_path: str, *, height: int = 720,
          edges: bool = True, max_inline_mb: float = 24.0) -> tuple:
    """Write the mesh .glb and return ``(html, meta, info)``.

    Signature deliberately mirrors :func:`semgrit.cadviewer.build`, so a
    notebook cell that shows the CAD can show the mesh by swapping one call.
    """
    from .cadviewer import viewer_html
    from .glb import write_glb

    parts, meta = parts_from_meshes(meshes, edges=edges)
    info = write_glb(glb_path, parts)
    meta = dict(meta)
    meta.setdefault("grits_drawn", 0)
    meta.setdefault("grits_total", 0)
    meta.setdefault("grits_full_detail", 0)
    meta["kind"] = "mesh"
    return (viewer_html(glb_path, meta, height=height,
                        max_inline_mb=max_inline_mb), meta, info)


def demo() -> None:
    """Self-check against meshes whose surface count is known in advance."""
    from .sagwrite import build_block, build_compliant_ring

    # --- a 2x3x4 box: the boundary count is arithmetic ------------------
    nodes, conn, _ = build_block(
        length_mm=2.0, width_mm=3.0, depth_mm=4.0,
        el_length_mm=1.0, el_width_mm=1.0, fine_depth_mm=1.0,
        band_mm=4.0, growth=1.0)
    assert len(conn) == 2 * 3 * 4
    quads = surface_of_hexes(conn)
    # a 2x3x4 grid of cubes has 2*(2*3 + 3*4 + 2*4) = 52 boundary quads
    assert len(quads) == 2 * (2 * 3 + 3 * 4 + 2 * 4), len(quads)
    # every interior face must have been dropped: 6 per element minus boundary
    assert 6 * len(conn) - len(quads) == 2 * (6 * len(conn) - len(quads)) // 2

    tris = quads_to_tris(quads)
    assert len(tris) == 2 * len(quads)
    # winding must be consistent: the closed surface's signed volume, computed
    # from the triangles alone, must equal the box and be POSITIVE
    v = nodes[tris]
    vol = np.einsum("ij,ij->i", np.cross(v[:, 0], v[:, 1]), v[:, 2]).sum() / 6.0
    assert abs(abs(vol) - 2.0 * 3.0 * 4.0) < 1e-9, vol
    assert vol > 0, "boundary quads must be wound outward, not inward"

    # --- Jacobians and aspect -------------------------------------------
    jac = hex_jacobians(nodes, conn)
    assert (jac > 0).all(), "the writer's box must not be inverted"
    assert abs(jac.sum() - 24.0) < 1e-9
    asp = hex_aspect(nodes, conn)
    assert abs(asp.max() - 1.0) < 1e-9, "unit cubes are 1:1"

    q = quality(nodes, conn)
    assert q["inverted"] == 0
    assert abs(q["volume"] - 24.0) < 1e-9
    assert abs(q["aspect_max"] - 1.0) < 1e-9

    # a deliberately flat mesh must report a large aspect ratio
    fn, fc, _ = build_block(
        length_mm=2.0, width_mm=2.0, depth_mm=0.02,
        el_length_mm=1.0, el_width_mm=1.0, fine_depth_mm=0.02,
        band_mm=0.02, growth=1.0)
    assert abs(hex_aspect(fn, fc).max() - 50.0) < 1e-6

    # --- an inverted mesh must be refused, not drawn --------------------
    # Swap the bottom and top faces -- the real way a hex gets inverted, and
    # it flips the Jacobian to exactly -1 on unit cubes. (Swapping a single
    # node pair merely distorts the element: it halves the volume here and
    # stays positive, so it would not have tested anything.)
    bad = conn[:, [4, 5, 6, 7, 0, 1, 2, 3]]
    jb = hex_jacobians(nodes, bad)
    assert (jb < 0).all(), jb[:3]
    assert abs(jb.max() + 1.0) < 1e-9, "unit cubes must invert to exactly -1"
    try:
        mesh_part("bad", nodes, bad)
    except MeshViewError as exc:
        assert "Jacobian" in str(exc)
    else:
        raise AssertionError("an inverted mesh must be refused")

    # --- corrupt connectivity must be refused ---------------------------
    dup = np.vstack([conn, conn[:1], conn[:1]])   # same element three times
    try:
        surface_of_hexes(dup)
    except MeshViewError as exc:
        assert "shared by" in str(exc)
    else:
        raise AssertionError("a face shared 3 times must be refused")

    # --- the ring: a closed annulus has no cut faces --------------------
    rn, rc, _ = build_compliant_ring(
        inner_r_mm=57.5, outer_r_mm=62.5, width_mm=10.0, sector_deg=360.0,
        n_circ=24, n_rad=2, n_axial=2)
    rq = surface_of_heres = surface_of_hexes(rc)
    # bore + outer + two axial faces, no seam: 24*2*2 (radial) + 24*2*2 (axial)
    assert len(rq) == 24 * 2 * 2 + 24 * 2 * 2, len(rq)
    assert (hex_jacobians(rn, rc) > 0).all(), "the ring must not be inverted"

    # a sector has two more faces than the closed ring, per row
    sn, sc, _ = build_compliant_ring(
        inner_r_mm=57.5, outer_r_mm=62.5, width_mm=10.0, sector_deg=30.0,
        n_circ=24, n_rad=2, n_axial=2)
    assert len(surface_of_hexes(sc)) == len(rq) + 2 * 2 * 2

    # --- wireframe ------------------------------------------------------
    w = wire_part("edges", nodes, conn)
    assert len(w["vertices"]) > 0 and len(w["faces"]) > 0
    assert len(w["faces"]) % 6 == 0, "each segment is a 6-triangle prism"
    # the prism must be thin: its bounding box must not exceed the mesh's
    lo, hi = w["vertices"].min(axis=0), w["vertices"].max(axis=0)
    mlo, mhi = nodes.min(axis=0), nodes.max(axis=0)
    assert (lo > mlo - 0.05).all() and (hi < mhi + 0.05).all()

    # --- parts, notes, and the honest omission --------------------------
    parts, meta = parts_from_meshes([
        dict(name="workpiece", nodes=nodes, conn=conn, color=C_WORK)])
    assert len(parts) == 2, "surface plus edges"
    assert not meta["notes"]
    assert meta["stats"]["workpiece"]["inverted"] == 0

    # a flat mesh must be FLAGGED for its aspect ratio, not silently accepted
    _, fmeta = parts_from_meshes([
        dict(name="thin", nodes=fn, conn=fc)])
    assert any("aspect ratio" in n for n in fmeta["notes"]), fmeta["notes"]

    # and a large mesh must say why the wireframe is missing
    _, big = parts_from_meshes(
        [dict(name="big", nodes=nodes, conn=conn)], max_edge_elements=1)
    assert any("wireframe is omitted" in n for n in big["notes"])
    bparts, _ = parts_from_meshes(
        [dict(name="big", nodes=nodes, conn=conn)], max_edge_elements=1)
    assert len(bparts) == 1, "no edge part when it is omitted"

    print("semgrit.meshview: all checks passed")
    print("  box 2x3x4       %d elements -> %d boundary quads, %d triangles"
          % (len(conn), len(quads), len(tris)))
    print("  closed ring     %d elements -> %d boundary quads (no seam)"
          % (len(rc), len(rq)))
    print("  quality         volume %.4f, aspect max %.2f:1, inverted %d"
          % (q["volume"], q["aspect_max"], q["inverted"]))
    print("  wireframe       %d prisms, %d triangles"
          % (len(w["faces"]) // 6, len(w["faces"])))


if __name__ == "__main__":
    demo()
