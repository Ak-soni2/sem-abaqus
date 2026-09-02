"""Write the SAG decks: a deformable compliant tool pressed into a workpiece.

This is the part the existing writers cannot do. ``rigid_wheel.py`` emits a
*discrete rigid* wheel -- a surface, no interior, no stiffness -- which is right
when the wheel is three orders stiffer than the cut and wrong here, because in
shape-adaptive grinding the tool's **deformation is the process**. A rigid SAG
tool would make line contact, load one grain, and reproduce conventional
grinding exactly.

So the tool is built as solid elements in three layers:

    hub          C3D8R, stiff, made rigid by *Rigid Body about a reference node
    compliant    C3D8R, hyperelastic + viscoelastic polyurethane  <-- the point
    pad          measured grains from the SEM pipeline, on the outer face

and driven in two steps, which is the physical sequence and what the reference
deck does: press in by the wheel compression T, then rotate.

THE TRANSITION CRITERION IS THE ENERGY ONE
------------------------------------------
``vumat_grind.for``'s geometric switch compares a *prescribed* chip thickness
``h(u)`` against ``dc``. That works for one grit on a known circle. It cannot
work here: with a compliant tool and hundreds of grains there is no closed-form
trajectory, and worse, the load per grain is not knowable in advance -- it is
what the contact solution produces.

``vumat_grind2.for``'s local criterion needs no geometry at all:

    W_p * L_c  >=  PSI * Kc^2 / E    ->  brittle

accumulated plastic work per unit volume, times the element's own length,
against a fracture energy. It triggers on HISTORY, so a point starts ductile and
turns brittle as work accumulates under repeated grain passes -- which is what a
polishing pad actually does to a surface, and it is why a coarse pad (few grains,
each hard-loaded) fractures where a fine pad (many grains, each light) does not.

With PSI left at 0 the subroutine derives it as ``dc*E*H/Kc^2``, which makes the
threshold exactly ``W_p L_c >= H dc``. Since the WC-Co card carries a MEASURED
dc, the energy criterion inherits that measurement and needs no separate
calibration. That is checked in :func:`demo`.

One consequence to state rather than discover later: the criterion is
regularised by the element length, so it is **mesh-dependent by construction**,
as every energy-based failure criterion is. Halving the element halves the work
density needed to trigger. PSI is therefore calibrated FOR A MESH, and every
deck this module writes states its element size next to PSI in the header.
"""

from __future__ import annotations

import math
import os
from typing import Optional, Sequence

import numpy as np

from . import sag, sagdeck

# Abaqus wants no more than 16 values on a data line; 8 is what this project
# uses for material cards because *User Material silently rejects 4.
_PER_LINE = 8


class SAGWriteError(RuntimeError):
    pass


def _fmt(v: float) -> str:
    """A number Abaqus will parse, short enough to keep lines sane.

    17 significant digits round-trip a float exactly, and the geometry here
    spans 125 mm down to 16 nm -- a ratio of 8e6 -- so trimming to 6 digits
    would quantise the fine mesh onto the coarse one.
    """
    return "%.12g" % float(v)


def _node_lines(nodes: np.ndarray, first: int = 1) -> list:
    return ["%d, %s, %s, %s" % (first + i, _fmt(p[0]), _fmt(p[1]), _fmt(p[2]))
            for i, p in enumerate(nodes)]


def _elem_lines(conn: np.ndarray, first: int = 1) -> list:
    return ["%d, %s" % (first + i, ", ".join(str(int(v) + 1) for v in row))
            for i, row in enumerate(conn)]


# ---------------------------------------------------------------------------
# the compliant tool, as solid elements
# ---------------------------------------------------------------------------

def build_compliant_ring(*, inner_r_mm: float, outer_r_mm: float,
                         width_mm: float, sector_deg: float,
                         n_circ: int, n_rad: int, n_axial: int,
                         z0_mm: float = 0.0) -> tuple:
    """Hex mesh of an annular sector, for the compliant layer or the hub.

    Returns ``(nodes, conn, faces)`` where ``faces`` names the node sets a
    caller needs: the bore (tied to the hub), the outer surface (carries the
    pad and contacts the work), and the two cut faces.

    Node ordering is the C3D8R convention -- bottom face counter-clockwise seen
    from outside, then the top face -- so the Jacobian is positive and Abaqus
    does not silently invert elements. That is asserted in :func:`demo` by
    computing the volume and requiring it positive.
    """
    if outer_r_mm <= inner_r_mm:
        raise SAGWriteError("outer radius must exceed inner radius")
    if not 0.0 < sector_deg <= 360.0:
        raise SAGWriteError("sector must be in (0, 360] degrees")
    if min(n_circ, n_rad, n_axial) < 1:
        raise SAGWriteError("every division count must be at least 1")

    full = abs(sector_deg - 360.0) < 1e-9
    n_t_nodes = n_circ if full else n_circ + 1
    th = np.deg2rad(np.linspace(0.0, sector_deg, n_circ + 1)[:n_t_nodes])
    rr = np.linspace(inner_r_mm, outer_r_mm, n_rad + 1)
    zz = np.linspace(z0_mm, z0_mm + width_mm, n_axial + 1)

    # index(i_t, i_r, i_z) -> node number
    def nid(it, ir, iz):
        return (it % n_t_nodes) * (n_rad + 1) * (n_axial + 1) \
            + ir * (n_axial + 1) + iz

    nodes = np.zeros((n_t_nodes * (n_rad + 1) * (n_axial + 1), 3))
    for it in range(n_t_nodes):
        c, s = math.cos(th[it]), math.sin(th[it])
        for ir in range(n_rad + 1):
            for iz in range(n_axial + 1):
                nodes[nid(it, ir, iz)] = (rr[ir] * c, rr[ir] * s, zz[iz])

    conn = []
    for it in range(n_circ):
        for ir in range(n_rad):
            for iz in range(n_axial):
                # bottom face, then top: counter-clockwise about +z
                conn.append([
                    nid(it, ir, iz), nid(it + 1, ir, iz),
                    nid(it + 1, ir + 1, iz), nid(it, ir + 1, iz),
                    nid(it, ir, iz + 1), nid(it + 1, ir, iz + 1),
                    nid(it + 1, ir + 1, iz + 1), nid(it, ir + 1, iz + 1),
                ])
    conn = np.array(conn, dtype=np.int64)

    bore = sorted({nid(it, 0, iz) for it in range(n_t_nodes)
                   for iz in range(n_axial + 1)})
    outer = sorted({nid(it, n_rad, iz) for it in range(n_t_nodes)
                    for iz in range(n_axial + 1)})
    faces = dict(bore=bore, outer=outer)
    if not full:
        faces["cut_lo"] = sorted({nid(0, ir, iz)
                                  for ir in range(n_rad + 1)
                                  for iz in range(n_axial + 1)})
        faces["cut_hi"] = sorted({nid(n_circ, ir, iz)
                                  for ir in range(n_rad + 1)
                                  for iz in range(n_axial + 1)})
    return nodes, conn, faces


def hex_volume(nodes: np.ndarray, conn: np.ndarray) -> float:
    """Total volume of a hex mesh, by decomposing each cell into tetrahedra.

    Exists to catch inverted elements: a negative total means the node ordering
    is wrong, which Abaqus reports as a cryptic preprocessing failure on a
    million-element deck. Cheap here, expensive there.
    """
    tets = ((0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6),
            (3, 4, 6, 7), (1, 4, 5, 6))
    # Vectorised over elements: a Python loop here is minutes on a
    # half-million-element deck, and this runs on every deck written.
    p = nodes[conn]                              # (n_el, 8, 3)
    tot = 0.0
    for a, b, c, d in tets:
        tot += np.abs(np.einsum(
            "ij,ij->i",
            np.cross(p[:, b] - p[:, a], p[:, c] - p[:, a]),
            p[:, d] - p[:, a])).sum() / 6.0
    return float(tot)


def ring_volume_exact(inner_r_mm: float, outer_r_mm: float, width_mm: float,
                      sector_deg: float) -> float:
    """Closed form, to check the mesh against."""
    return (math.pi * (outer_r_mm ** 2 - inner_r_mm ** 2) * width_mm
            * sector_deg / 360.0)


# ---------------------------------------------------------------------------
# the graded workpiece
# ---------------------------------------------------------------------------

def graded_depth_planes(depth_mm: float, fine_mm: float, band_mm: float,
                        growth: float, max_mm: float = 0.0) -> np.ndarray:
    """Depth coordinates: uniform ``fine_mm`` through ``band_mm``, then growing.

    The surface band is where the transition happens and must resolve dc; below
    it nothing needs that resolution, and paying for it would be the difference
    between a deck that runs and one that does not. Same idea as
    ``wp_surface_layer_mm`` / ``wp_depth_growth`` in the rigid pipeline.
    """
    if fine_mm <= 0 or depth_mm <= 0:
        raise SAGWriteError("depth and element size must be positive")
    if growth < 1.0:
        raise SAGWriteError("growth ratio must be >= 1")
    band = min(band_mm, depth_mm)
    z = [0.0]
    while z[-1] < band - 1e-15:
        z.append(min(z[-1] + fine_mm, band))
    step = fine_mm
    cap = max_mm if max_mm > 0 else depth_mm
    while z[-1] < depth_mm - 1e-15:
        step = min(step * growth, cap)
        z.append(min(z[-1] + step, depth_mm))
    return np.array(z)


def build_block(*, length_mm: float, width_mm: float, depth_mm: float,
                el_length_mm: float, el_width_mm: float,
                fine_depth_mm: float, band_mm: float, growth: float,
                max_depth_el_mm: float = 0.0,
                x0_mm: float = 0.0, y0_mm: float = 0.0,
                top_z_mm: float = 0.0) -> tuple:
    """Hex mesh of the workpiece, graded in depth. Top face at ``top_z_mm``."""
    nx = max(int(round(length_mm / el_length_mm)), 1)
    ny = max(int(round(width_mm / el_width_mm)), 1)
    zs = graded_depth_planes(depth_mm, fine_depth_mm, band_mm, growth,
                             max_depth_el_mm)
    xs = np.linspace(x0_mm, x0_mm + length_mm, nx + 1)
    ys = np.linspace(y0_mm, y0_mm + width_mm, ny + 1)
    zc = top_z_mm - zs                      # downwards from the ground face
    nz = len(zc) - 1

    def nid(i, j, k):
        return i * (ny + 1) * (nz + 1) + j * (nz + 1) + k

    nodes = np.zeros(((nx + 1) * (ny + 1) * (nz + 1), 3))
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                nodes[nid(i, j, k)] = (xs[i], ys[j], zc[k])

    conn = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # k increases downwards, so the k+1 face is the lower one and
                # goes FIRST to keep the Jacobian positive.
                conn.append([
                    nid(i, j, k + 1), nid(i + 1, j, k + 1),
                    nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                    nid(i, j, k), nid(i + 1, j, k),
                    nid(i + 1, j + 1, k), nid(i, j + 1, k),
                ])
    conn = np.array(conn, dtype=np.int64)

    top = sorted({nid(i, j, 0) for i in range(nx + 1) for j in range(ny + 1)})
    bottom = sorted({nid(i, j, nz) for i in range(nx + 1)
                     for j in range(ny + 1)})
    sides = sorted({nid(i, j, k)
                    for i in (0, nx) for j in range(ny + 1)
                    for k in range(nz + 1)}
                   | {nid(i, j, k)
                      for j in (0, ny) for i in range(nx + 1)
                      for k in range(nz + 1)})
    return nodes, conn, dict(top=top, bottom=bottom, sides=sides,
                             nx=nx, ny=ny, nz=nz, z_planes=zs)


# ---------------------------------------------------------------------------
# deck assembly
# ---------------------------------------------------------------------------

def _material_block(pl: dict, name: str, psi: float,
                    element_mm: float) -> list:
    """The workpiece material: the 58-constant vumat_grind2 card.

    Written directly at 58 rather than writing 56 and rewriting it, because
    there is no already-built deck here to patch -- ``semgrit_multi.swmode``
    exists to convert an EXISTING deck and is the right tool for that job, not
    this one.
    """
    from . import materials
    from .hybrid import hybrid_props

    w = materials.get(pl["params"].material)
    hp = w.hybrid_params()
    props = list(hybrid_props(w.jh2, hp, None))
    if len(props) != 56:
        raise SAGWriteError("expected 56 hybrid constants, got %d"
                            % len(props))
    # 57 = SWMODE, 58 = PSI. SWMODE 1 is energy-only: no chip thickness is
    # read at all, which is the whole reason it suits SAG.
    props += [1.0, float(psi)]

    dc_mm = hp.critical_depth_mm()
    thresh = hp.hardness_mpa * dc_mm

    out = ["*Material, name=%s" % name,
           "*Density", " %s," % _fmt(hp.density_kg_m3 * 1e-12),
           "*Depvar, delete=12", " 22,",
           "*User Material, constants=58"]
    for i in range(0, len(props), _PER_LINE):
        out.append(", ".join(_fmt(v) for v in props[i:i + _PER_LINE]))
    out += [
        "**",
        "** SWMODE = 1, the LOCAL ENERGY criterion. No chip thickness is read.",
        "**   brittle once  W_p * L_c >= PSI * Kc^2 / E",
        "**   PSI = %s%s" % (_fmt(psi),
                             "  (0 -> the subroutine derives dc*E*H/Kc^2)"
                             if psi <= 0 else ""),
        "**   which is  W_p * L_c >= H * dc = %s MPa*mm = %.1f J/m2"
        % (_fmt(thresh), thresh * 1000.0),
        "**   with dc = %.2f nm, %s." % (dc_mm * 1e6,
                                         "MEASURED" if w.dc_measured
                                         else "computed"),
        "**",
        "** THE CRITERION IS MESH-DEPENDENT BY CONSTRUCTION. It is regularised",
        "** by L_c, so halving the element halves the work density needed to",
        "** trigger -- correct behaviour for a fracture-energy criterion, and",
        "** the reason PSI is calibrated FOR A MESH. This deck's surface",
        "** element is %.4f nm through the depth." % (element_mm * 1e6),
        "**",
        "** It triggers on HISTORY: a point starts ductile and turns brittle as",
        "** plastic work accumulates under repeated grain passes. That is the",
        "** physical transition, rather than being told in advance where it is.",
    ]
    return out


def write_micro_deck(path: str, pl: dict, solids: Sequence, *,
                     psi: float = 0.0) -> dict:
    """The resolved deck: a patch of contact, dc/5 deep, energy criterion.

    Grains are driven by the load per grain the contact solution gives, applied
    as a concentrated force on a rigid grain rather than by prescribing a
    depth -- because the depth is what the model is supposed to predict, and
    prescribing it would assume the answer.
    """
    from . import materials

    p: sagdeck.SAGParams = pl["params"]
    mic = pl["micro"]
    c: sag.SAGContact = pl["contact"]
    w = materials.get(p.material)
    hp = w.hybrid_params()

    if not solids:
        raise SAGWriteError("the grain library is empty")

    side = mic["side_mm"]
    nodes, conn, meta = build_block(
        length_mm=side, width_mm=side, depth_mm=mic["depth_mm"],
        el_length_mm=mic["element_inplane_mm"],
        el_width_mm=mic["element_inplane_mm"],
        fine_depth_mm=mic["element_mm"],
        band_mm=min(mic["depth_mm"], 10.0 * hp.critical_depth_mm()),
        growth=1.25, x0_mm=-0.5 * side, y0_mm=-0.5 * side)

    vol = hex_volume(nodes, conn)
    want = side * side * mic["depth_mm"]
    if abs(vol - want) / want > 1e-6:
        raise SAGWriteError(
            "workpiece mesh volume %.6g mm3 does not match the block %.6g -- "
            "elements are inverted or the grading is wrong" % (vol, want))

    L = sagdeck.micro_header(pl)
    L += ["*Heading", "** SAG MICRO -- %s" % p.name, "*Preprint, echo=NO,"
          " model=NO, history=NO, contact=NO", "**",
          "*Part, name=WORKPIECE", "*Node"]
    L += _node_lines(nodes)
    L += ["*Element, type=C3D8R"]
    L += _elem_lines(conn)
    L += ["*Nset, nset=ALL, generate", " 1, %d, 1" % len(nodes),
          "*Elset, elset=ALL, generate", " 1, %d, 1" % len(conn),
          "*Nset, nset=TOP"]
    L += _pack_ids(meta["top"])
    L += ["*Nset, nset=FIXED"]
    L += _pack_ids(sorted(set(meta["bottom"]) | set(meta["sides"])))
    L += ["*Solid Section, elset=ALL, controls=EC1, material=%s"
          % w.inp_material, " ,", "*End Part", "**"]

    info = dict(
        kind="micro", path=path,
        nodes=len(nodes), elements=len(conn),
        block_mm=(side, side, mic["depth_mm"]),
        element_depth_mm=mic["element_mm"],
        element_inplane_mm=mic["element_inplane_mm"],
        volume_mm3=vol, volume_error=abs(vol - want) / want,
        dc_nm=hp.critical_depth_mm() * 1e6,
        dc_measured=w.dc_measured,
        psi=psi, swmode=1,
        load_per_grain_n=c.load_per_grain_n,
        energy_threshold_mpa_mm=hp.hardness_mpa * hp.critical_depth_mm(),
        elements_per_dc=p.elements_per_dc,
        n_depth_planes=len(meta["z_planes"]),
    )
    return info


def _pack_ids(ids: Sequence[int], per_line: int = 16) -> list:
    out = []
    ids = [int(i) + 1 for i in ids]
    for i in range(0, len(ids), per_line):
        out.append(", ".join(str(v) for v in ids[i:i + per_line]) + ",")
    return out


def demo() -> None:
    """Self-check the geometry and the energy-criterion wiring."""
    # --- the ring mesh must be right, and provably so -------------------
    nodes, conn, faces = build_compliant_ring(
        inner_r_mm=57.5, outer_r_mm=62.5, width_mm=10.0, sector_deg=30.0,
        n_circ=12, n_rad=3, n_axial=4)
    assert len(conn) == 12 * 3 * 4
    v = hex_volume(nodes, conn)
    exact = ring_volume_exact(57.5, 62.5, 10.0, 30.0)
    # a 12-facet chord under-runs the arc, so the mesh is slightly SMALLER
    assert v < exact, "a faceted sector cannot exceed the true annulus"
    assert abs(v - exact) / exact < 0.02, (v, exact)
    # refining must converge toward the closed form
    n2, c2, _ = build_compliant_ring(
        inner_r_mm=57.5, outer_r_mm=62.5, width_mm=10.0, sector_deg=30.0,
        n_circ=48, n_rad=3, n_axial=4)
    assert abs(hex_volume(n2, c2) - exact) / exact < abs(v - exact) / exact

    # every radius must land in the annulus
    r = np.hypot(nodes[:, 0], nodes[:, 1])
    assert r.min() > 57.5 - 1e-9 and r.max() < 62.5 + 1e-9
    # bore and outer sets must sit at the right radii
    assert abs(r[faces["bore"]].max() - 57.5) < 1e-9
    assert abs(r[faces["outer"]].min() - 62.5) < 1e-9
    # a sector has two cut faces; a full ring has none
    assert "cut_lo" in faces and "cut_hi" in faces
    _, _, full_faces = build_compliant_ring(
        inner_r_mm=57.5, outer_r_mm=62.5, width_mm=10.0, sector_deg=360.0,
        n_circ=24, n_rad=2, n_axial=2)
    assert "cut_lo" not in full_faces, "a closed ring must have no seam"

    # --- graded depth planes -------------------------------------------
    z = graded_depth_planes(0.002, 16e-6, 0.0008, 1.25)
    assert abs(z[0]) < 1e-15 and abs(z[-1] - 0.002) < 1e-12
    d = np.diff(z)
    assert (d > 0).all(), "depth planes must increase"
    assert abs(d[0] - 16e-6) < 1e-12, "the surface band must be at dc/5"
    assert d[-1] > d[0], "it must coarsen with depth"
    # Grading must save SOMETHING, but not much here and that is correct: the
    # band is 10*dc of a 20*dc block, so half the depth is legitimately fine --
    # the transition is a near-surface phenomenon and that is where the
    # elements belong. The large saving in this deck is IN-PLANE (grain/20
    # rather than dc/5), which is asserted on the written deck below.
    assert len(z) < 0.002 / 16e-6, "grading must beat a uniform fine mesh"
    assert len(z) > 0.0008 / 16e-6, "the fine band must be fully resolved"

    # --- the block mesh -------------------------------------------------
    bn, bc, bm = build_block(
        length_mm=0.05, width_mm=0.05, depth_mm=0.002,
        el_length_mm=0.0003, el_width_mm=0.0003,
        fine_depth_mm=16e-6, band_mm=0.0008, growth=1.25)
    bv = hex_volume(bn, bc)
    assert abs(bv - 0.05 * 0.05 * 0.002) / (0.05 * 0.05 * 0.002) < 1e-9, \
        "a box mesh must reproduce the box volume exactly"
    assert bm["nz"] + 1 == len(bm["z_planes"])
    # the top face must be the ground surface
    assert abs(bn[bm["top"], 2].max()) < 1e-15

    # --- the energy criterion is wired to the MEASURED dc ---------------
    pl = sagdeck.plan(sagdeck.SAGParams(grain_um=6.0, material="wc_co",
                                        micro_grains=1))
    blk = _material_block(pl, "WCCO", 0.0, pl["micro"]["element_mm"])
    txt = "\n".join(blk)
    assert "constants=58" in txt, "grind2 needs 58 constants, not 56"
    assert "*Depvar, delete=12" in txt
    assert " 22," in txt, "grind2 writes 22 state variables"
    assert "SWMODE = 1" in txt
    assert "MESH-DEPENDENT BY CONSTRUCTION" in txt
    assert "MEASURED" in txt, "the WC-Co dc is measured and must say so"
    # the card body must be 8 values per line -- 4 is silently rejected
    body = [ln for ln in blk if ln and ln[0] not in "*" and "," in ln
            and not ln.startswith(" ") and not ln.startswith("**")]
    assert body, "no material data lines were written"
    assert len(body[0].split(",")) == _PER_LINE, body[0]

    # PSI derived from a MEASURED dc must reproduce H*dc exactly
    from . import materials
    w = materials.get("wc_co")
    hp = w.hybrid_params()
    dc = hp.critical_depth_mm()
    psi = dc * hp.youngs_mpa * hp.hardness_mpa / (hp.kic ** 2)
    assert abs(psi * hp.kic ** 2 / hp.youngs_mpa
               - hp.hardness_mpa * dc) < 1e-12, \
        "PSI must make the energy threshold equal H*dc"
    assert abs(dc * 1e6 - 80.0) < 1e-9, "and dc must be the measurement"

    # --- writing a micro deck ------------------------------------------
    class _S:
        """Smallest thing the writer needs: it only counts them here."""
        height_um = 5.0

    info = write_micro_deck("_unused.inp", pl, [_S()])
    assert info["swmode"] == 1
    assert info["dc_measured"]
    assert info["volume_error"] < 1e-6
    assert info["elements"] > 1000
    assert abs(info["element_depth_mm"] * 1e6 - 80.0 / 5.0) < 1e-9
    assert info["element_inplane_mm"] > info["element_depth_mm"]
    # and the in-plane anisotropy is what makes the deck affordable: at dc/5
    # in-plane the same patch would be ~350x larger
    iso = info["elements"] * (info["element_inplane_mm"]
                              / info["element_depth_mm"]) ** 2
    assert iso / info["elements"] > 100.0, iso / info["elements"]
    assert info["elements"] < 3e6, info["elements"]

    # --- refusals -------------------------------------------------------
    for bad, why in (
        (lambda: build_compliant_ring(inner_r_mm=62.5, outer_r_mm=57.5,
                                      width_mm=10.0, sector_deg=30.0,
                                      n_circ=4, n_rad=1, n_axial=1),
         "inverted radii"),
        (lambda: build_compliant_ring(inner_r_mm=57.5, outer_r_mm=62.5,
                                      width_mm=10.0, sector_deg=0.0,
                                      n_circ=4, n_rad=1, n_axial=1),
         "zero sector"),
        (lambda: build_compliant_ring(inner_r_mm=57.5, outer_r_mm=62.5,
                                      width_mm=10.0, sector_deg=30.0,
                                      n_circ=0, n_rad=1, n_axial=1),
         "zero divisions"),
        (lambda: graded_depth_planes(0.002, 0.0, 0.0008, 1.25),
         "zero element size"),
        (lambda: graded_depth_planes(0.002, 16e-6, 0.0008, 0.5),
         "shrinking growth ratio"),
        (lambda: write_micro_deck("_unused.inp", pl, []),
         "empty grain library"),
    ):
        try:
            bad()
        except SAGWriteError:
            pass
        else:
            raise AssertionError("should have been refused: %s" % why)

    print("semgrit.sagwrite: all checks passed")
    print("  ring mesh   %d elements, volume %.4f mm3 (exact %.4f, %.3f%% low)"
          % (len(conn), v, exact, 100.0 * (exact - v) / exact))
    print("  micro deck  %s elements, %s nodes, %.4f nm depth element"
          % (format(info["elements"], ","), format(info["nodes"], ","),
             info["element_depth_mm"] * 1e6))
    print("  energy      W_p*L_c >= %.4f MPa*mm = %.1f J/m2 (dc = %.1f nm, %s)"
          % (info["energy_threshold_mpa_mm"],
             info["energy_threshold_mpa_mm"] * 1000.0,
             info["dc_nm"], "MEASURED" if info["dc_measured"] else "computed"))


if __name__ == "__main__":
    demo()
