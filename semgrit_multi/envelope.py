"""Undeformed chip thickness for every workpiece element, any number of grits.

The idea
--------
The wheel in these decks is ONE discrete rigid body driven by a prescribed
velocity boundary condition. It cannot deflect, cannot slow down and cannot be
pushed back, so the position of every grit at every instant of the step is
known in closed form before Abaqus is even started. Sweeping that motion and
recording how deep the grit surface reaches gives the undeformed chip thickness
directly -- for one grit or for seven hundred, and without the subroutine
needing to know anything about wheels.

Kinematics, exactly
-------------------
Work in the block frame, columns ``(a, b, z)`` = (radial, tangential, axial),
which is the frame ``rigid_wheel.place_workpiece`` already returns the grit
vertices in. The wheel rotates about Z at the signed rate the deck writes as
VR3, and translates along ``+e_r`` at the infeed speed ``v_r``. Because the
frame itself is a rotation about Z, the rotation looks the same inside it:

    a(t) = a0 cos(w t) - b0 sin(w t) + v_r t
    b(t) = a0 sin(w t) + b0 cos(w t)
    z(t) = z0                                   w = VR3

and the depth of that point below the original ground face at radius
``r_ground`` is simply ``a(t) - r_ground``. No approximation: this is the same
rigid-body motion Abaqus will integrate from the same boundary condition.

Time order, and why it matters
------------------------------
Undeformed chip thickness is depth below the surface *as the grit finds it*,
not depth below the original surface. With one grit those are the same thing.
With several they are not: the second grit over a station only removes what the
first one left. So the sweep runs in time order and carries a running removed
depth ``D(u, z)``:

    grit i takes        h_i = max(0, d_i - D)
    material between    D < depth <= d_i        is removed by grit i,
                        and its chip thickness is h_i
    then                D <- max(D, d_i)

For a single grit ``D`` starts at zero, ``h_1 = d_1``, and every point at a
station shares one chip thickness -- which is exactly what the closed-form
single-grit model does, and what ``verify_envelope.py`` checks.

Elements the grit never reaches keep the chip thickness of the deepest grit at
their station. They are subsurface and are never removed, so the choice only
decides which law governs their (small) damage, and inheriting the surface
regime is the defensible reading.

What this does NOT do
---------------------
It assumes the wheel's motion is prescribed, which is true here and stops being
true if the wheel is ever made deformable or its grits are allowed to wear
during the run. It also ignores elastic deflection of the workpiece, which is
the right call: the quantity the transition criterion is calibrated against is
the *undeformed* chip thickness, a kinematic quantity by definition.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np


class EnvelopeError(RuntimeError):
    pass


@dataclass
class EnvelopeParams:
    """How finely to sweep.

    The sample count is NOT a round number chosen by taste. Depth grows with
    the infeed at ``v_r``, so a time step of ``dt`` blurs the recorded depth by
    ``v_r*dt``. With ``v_r`` of order 1800 mm/s and ``dc`` of order 5 nm, 240
    samples over a 1 us window blurs the depth by 8 nm -- larger than the
    threshold the whole model turns on. So the sample count is derived from a
    stated depth resolution instead, and the resolution is reported.
    """

    depth_resolution_mm: float = 2.0e-7
    """Target blur in the recorded depth, mm. 2e-7 mm = 0.2 nm, which is a
    twenty-fifth of a typical dc."""

    max_time_samples: int = 200_000
    min_time_samples: int = 64
    chunk_samples: int = 2_048
    """Time samples per vectorised block, to bound peak memory."""

    facet_subdivision: int = 2
    """Barycentric samples per facet edge. 1 = vertices only, 2 adds edge
    midpoints and the centroid, 3 subdivides again. Vertices alone
    under-report the envelope wherever a facet spans more than one cell."""

    margin_factor: float = 1.25
    """How far outside the block, in units of the grit's own radius, to keep
    sweeping. Covers the case of a facet whose centre is outside the block but
    whose corner is inside."""

    fill_never_cut: bool = True
    """Give material no grit reaches the chip thickness of the nearest station
    that was cut, instead of leaving it at zero.

    Zero would classify it ductile, which is wrong twice over: it inflates the
    ductile count with material that was never removed at all, and it puts
    subsurface material under a brittle groove on the ductile law, when the
    crack field it will actually see is the brittle one."""

    def validate(self) -> None:
        if self.depth_resolution_mm <= 0:
            raise EnvelopeError("depth_resolution_mm must be positive")
        if self.facet_subdivision < 1:
            raise EnvelopeError("facet_subdivision must be at least 1")
        if self.min_time_samples < 8:
            raise EnvelopeError("min_time_samples must be at least 8")
        if self.chunk_samples < 1:
            raise EnvelopeError("chunk_samples must be at least 1")

    def samples_for(self, window_s: float, v_r: float) -> int:
        """Samples needed to hold the depth blur under the target."""
        if v_r <= 0:
            return self.min_time_samples
        n = int(math.ceil(window_s * v_r / self.depth_resolution_mm)) + 1
        return int(min(max(n, self.min_time_samples), self.max_time_samples))


@dataclass
class ChipEnvelope:
    """The result of a sweep."""

    h_elem: np.ndarray
    """(nl, nw, nd) undeformed chip thickness per workpiece element, mm, in the
    same index order ``wheel_workpiece.build_block_mesh`` numbers them."""
    depth_removed: np.ndarray
    """(nl, nw) total depth the whole grit population removes at each station,
    mm. This is the groove profile the run should end up with."""
    n_grits_engaged: int
    grit_order: list
    """Grit indices in the order they cross the block."""
    per_grit_h: dict
    """grit index -> (max, mean) of its own chip thickness where it cuts, mm."""
    u_edges: np.ndarray
    z_edges: np.ndarray
    depth_edges: np.ndarray
    stats: dict = field(default_factory=dict)

    cut: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    """(nl, nw, nd) True where a grit actually removes this element."""

    def h_at_element(self, i: int, j: int, k: int) -> float:
        return float(self.h_elem[i, j, k])

    def split(self, dc_mm: float) -> dict:
        """How the elements divide, counting "never cut" separately.

        Lumping never-cut material in with the ductile elements makes the
        ductile fraction meaningless -- on a block much wider than one grit,
        most of it is simply untouched. The number worth quoting is the split
        among the elements the grits actually remove.
        """
        h = self.h_elem
        tot = int(h.size)
        cut = (self.cut if self.cut.size == h.size
               else np.ones(h.shape, dtype=bool))
        n_cut = int(cut.sum())
        duct_cut = int((cut & (h < dc_mm)).sum())
        return {
            "dc_mm": dc_mm,
            "n_elements": tot,
            "n_cut": n_cut,
            "n_never_cut": tot - n_cut,
            "n_ductile_of_cut": duct_cut,
            "n_brittle_of_cut": n_cut - duct_cut,
            "ductile_fraction_of_cut": (duct_cut / n_cut) if n_cut else 0.0,
            # the law each element will actually run, cut or not
            "n_ductile_law": int((h < dc_mm).sum()),
            "n_brittle_law": int((h >= dc_mm).sum()),
        }


# --------------------------------------------------------------------------
# surface sampling
# --------------------------------------------------------------------------

def _barycentric(n: int) -> np.ndarray:
    """Barycentric weights for a triangle subdivided ``n`` times per edge."""
    out = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            out.append((i / n, j / n, k / n))
    return np.asarray(out, dtype=np.float64)


def sample_grit_surface(verts: np.ndarray, faces: Sequence,
                        subdivision: int = 2) -> np.ndarray:
    """Point cloud covering a grit's triangulated surface.

    Vertices alone are not enough. A facet is a few tenths of a micron across,
    the station grid is one element wide, and the deepest point of a facet is
    generally in its interior once the facet is tilted -- so a vertex-only
    cloud under-reports the envelope by up to a facet's own relief.
    """
    v = np.asarray(verts, dtype=np.float64)
    f = np.asarray(faces, dtype=np.int64)
    if f.ndim != 2 or f.shape[1] != 3:
        raise EnvelopeError("faces must be (K, 3) triangles")
    w = _barycentric(max(int(subdivision), 1))
    tri = v[f]                                    # (K, 3, 3)
    # (K, S, 3) = sum over the triangle's corners of weight * corner
    pts = np.einsum("sc,kcd->ksd", w, tri).reshape(-1, 3)
    return np.unique(np.round(pts, 12), axis=0)


# --------------------------------------------------------------------------
# the sweep
# --------------------------------------------------------------------------

def _crossing_window(a_c: float, b_c: float, reach: float, omega: float,
                     half_len: float, step_time: float,
                     margin: float) -> Optional[tuple]:
    """When this grit's station is within the block, or None if never.

    The centre crosses ``b = 0`` when ``a_c sin(w t) + b_c cos(w t) = 0``, i.e.
    ``w t = atan2(-b_c, a_c)``. The window around it is the time to travel the
    block half-length plus a margin at the surface speed.
    """
    if omega == 0.0:
        return None
    alpha = math.atan2(-b_c, a_c)
    # atan2 gives the crossing nearest t = 0; the wheel may need one more turn,
    # which for a step far shorter than a revolution never happens, but adding
    # the branch costs nothing and stops a silent miss.
    for turn in (0.0, math.copysign(2.0 * math.pi, omega)):
        t_mid = (alpha + turn) / omega
        span = (half_len + margin) / max(abs(omega) * max(reach, 1e-12), 1e-30)
        t_lo, t_hi = t_mid - span, t_mid + span
        if t_hi < 0.0 or t_lo > step_time:
            continue
        return max(t_lo, 0.0), min(t_hi, step_time)
    return None


def sweep_envelope(place: dict, motion: Optional[dict], wp,
                   *, step_time_s: float, rotation_reversed: bool = False,
                   params: Optional[EnvelopeParams] = None,
                   paths: Optional[dict] = None,
                   log=None) -> ChipEnvelope:
    """Undeformed chip thickness per element, from the deck's own placement.

    ``place`` is :func:`semgrit.rigid_wheel.place_workpiece`'s output and
    ``motion`` is :func:`semgrit.analysis.wheel_motion`'s, so the sweep uses the
    same grit positions and the same infeed the deck writes. Passing anything
    re-derived would let the field disagree with the deck it is written into.

    ``paths`` overrides the kinematics for the grits named in it:
    ``{grit_index: (n, 4) array of [t, u_tip, z_tip, depth_tip]}``. The grit is
    then translated so its own outermost point follows that path, keeping its
    orientation, instead of being carried round by the wheel. That is how a
    trajectory measured off a real groove is replayed -- see
    :mod:`semgrit_multi.trajectory`. Everything downstream is identical, because
    what comes back is the same :class:`ChipEnvelope`.
    """
    params = params or EnvelopeParams()
    params.validate()
    if wp is None:
        raise EnvelopeError("a chip-thickness field needs a workpiece")
    frames = place["frames"]
    faces = place["faces"]
    if not frames:
        raise EnvelopeError("no grits were placed")
    r_ground = float(place["r_ground"])
    if step_time_s <= 0:
        raise EnvelopeError("the step has no duration, so nothing sweeps")

    omega = 0.0 if motion is None else float(motion["vr3"])
    v_r = 0.0 if motion is None else float(motion["radial_speed_mm_s"])
    if omega == 0.0:
        raise EnvelopeError(
            "VR3 is zero: the wheel does not turn, so no grit crosses the "
            "block and the undeformed chip thickness is undefined")
    # The deck's own sign convention, cross-checked rather than assumed.
    if (omega < 0.0) != (not rotation_reversed):
        raise EnvelopeError(
            "VR3 = %g contradicts rotation_reversed = %s; the sweep would run "
            "the grits the wrong way across the block"
            % (omega, rotation_reversed))

    nl, nw, nd = wp.divisions()
    hl, hw = wp.length_mm / 2.0, wp.width_mm / 2.0
    u_edges = np.linspace(-hl, hl, nl + 1)
    z_edges = np.linspace(-hw, hw, nw + 1)
    d_edges = wp.depth_coordinates()
    if len(d_edges) - 1 != nd:
        raise EnvelopeError("depth coordinates disagree with divisions()")
    d_centre = 0.5 * (d_edges[:-1] + d_edges[1:])
    du, dz = (2 * hl) / nl, (2 * hw) / nw

    # --- order the grits by when they cross the block -------------------
    paths = paths or {}
    order = []
    for gi, v in enumerate(frames):
        c = v.mean(axis=0)
        reach = float(np.hypot(v[:, 0], v[:, 1]).max())
        rad = float(np.linalg.norm(v - c, axis=1).max())
        if gi in paths:
            # A measured path carries its own timing, so the analytic crossing
            # window would be answering a question about a motion this grit is
            # not making.
            pt = np.asarray(paths[gi], dtype=np.float64)
            if pt.ndim != 2 or pt.shape[1] != 4 or len(pt) < 2:
                raise EnvelopeError(
                    "paths[%d] must be an (n, 4) array of [t, u, z, depth] "
                    "with at least two samples" % gi)
            win = (float(pt[:, 0].min()), float(pt[:, 0].max()))
            order.append((0.5 * (win[0] + win[1]), gi, win))
            continue
        win = _crossing_window(float(c[0]), float(c[1]), reach, omega, hl,
                               step_time_s, params.margin_factor * rad)
        if win is None:
            continue
        # Also skip grits that miss across the face: no rotation moves them in z.
        if float(np.abs(v[:, 2]).min()) > hw + params.margin_factor * rad:
            continue
        order.append((0.5 * (win[0] + win[1]), gi, win))
    order.sort()
    if not order:
        raise EnvelopeError(
            "no grit crosses the workpiece during the step. Either the step is "
            "too short for the wheel to bring one round, or the grits sit "
            "outside the block's width.")

    # --- sweep each grit, in time order --------------------------------
    removed = np.zeros((nl, nw), dtype=np.float64)
    h_elem = np.full((nl, nw, nd), np.nan, dtype=np.float64)
    cut_mask = np.zeros((nl, nw, nd), dtype=bool)
    best_share = np.zeros((nl, nw, nd), dtype=np.float64)
    last_h = np.zeros((nl, nw), dtype=np.float64)
    per_grit = {}
    n_pts = 0

    n_samples_used = []
    for _t, gi, (t_lo, t_hi) in order:
        cloud = sample_grit_surface(frames[gi], faces[gi],
                                    params.facet_subdivision)
        n_pts += len(cloud)
        a0, b0, z0 = cloud[:, 0], cloud[:, 1], cloud[:, 2]
        pt = paths.get(gi)
        if pt is not None:
            # Reference the cloud to the grit's own outermost point, so putting
            # that point on the measured path puts the whole grit where the real
            # abrasive was.
            pt = np.asarray(pt, dtype=np.float64)
            tip = int(np.argmax(a0))
            ra0 = a0 - a0[tip]
            rb0 = b0 - b0[tip]
            rz0 = z0 - z0[tip]
            ns = len(pt)
        else:
            ns = params.samples_for(t_hi - t_lo, abs(v_r))
        n_samples_used.append(ns)
        d_i = np.zeros((nl, nw), dtype=np.float64)
        touched = np.zeros((nl, nw), dtype=bool)
        # Chunked in time so peak memory stays at chunk_samples x cloud, not
        # samples x cloud: at 0.2 nm depth resolution the sample count runs to
        # tens of thousands.
        for s0 in range(0, ns, params.chunk_samples):
            s1 = min(s0 + params.chunk_samples, ns)
            if pt is not None:
                # Translate only: a measured groove says where the abrasive was,
                # not how it was turned.
                seg = pt[s0:s1]
                depth = (ra0[None, :] + seg[:, 3][:, None]).ravel()
                bb = (rb0[None, :] + seg[:, 1][:, None]).ravel()
                zz = (rz0[None, :] + seg[:, 2][:, None]).ravel()
            else:
                ts = t_lo + (t_hi - t_lo) * (
                    np.arange(s0, s1, dtype=np.float64) / max(ns - 1, 1))
                al = omega * ts
                ca, sa = np.cos(al)[:, None], np.sin(al)[:, None]
                a = a0[None, :] * ca - b0[None, :] * sa + (v_r * ts)[:, None]
                b = a0[None, :] * sa + b0[None, :] * ca
                depth = (a - r_ground).ravel()
                bb = b.ravel()
                zz = np.broadcast_to(z0[None, :], b.shape).ravel()
            keep = ((depth > 0.0) & (bb >= -hl) & (bb < hl)
                    & (zz >= -hw) & (zz < hw))
            if not keep.any():
                continue
            iu = np.clip(((bb[keep] + hl) / du).astype(np.int64), 0, nl - 1)
            iz = np.clip(((zz[keep] + hw) / dz).astype(np.int64), 0, nw - 1)
            np.maximum.at(d_i, (iu, iz), depth[keep])
            touched[iu, iz] = True

        # This grit's own chip thickness: what it removes beyond what is gone.
        h_i = np.maximum(d_i - removed, 0.0)
        cuts = touched & (d_i > removed)
        if cuts.any():
            hv = h_i[cuts]
            per_grit[gi] = (float(hv.max()), float(hv.mean()))
        else:
            per_grit[gi] = (0.0, 0.0)
            continue

        # Assign it to the elements it takes material out of, by the exact
        # overlap of the cut with each element's own depth interval.
        #
        # Testing the element CENTRE instead was wrong in a way that matters:
        # a 1.7 nm cut into a 300 nm element never reaches the centre, so a
        # nanometre-scale infeed -- exactly the ductile regime this model
        # exists to capture -- registered as nothing cut at all. An element the
        # cut enters at all has experienced a cut of that station's thickness,
        # and whether it is removed or ploughed is the constitutive law's
        # business, not the sweep's.
        #
        # Where two grits both take material out of one element, the element
        # keeps the chip thickness of whichever removed the larger share of it.
        for k in range(nd):
            lo, hi = d_edges[k], d_edges[k + 1]
            start = np.maximum(lo, removed)
            end = np.minimum(hi, d_i)
            share = end - start
            sel = share > 0.0
            if not sel.any():
                continue
            better = sel & (share > best_share[:, :, k])
            if better.any():
                h_elem[:, :, k][better] = h_i[better]
                best_share[:, :, k][better] = share[better]
            cut_mask[:, :, k] |= sel
        removed = np.maximum(removed, d_i)
        last_h = np.where(cuts, h_i, last_h)

    # --- elements no grit reaches --------------------------------------
    untouched = np.isnan(h_elem)
    n_never = int(untouched.sum())
    if n_never and params.fill_never_cut:
        h_elem = _fill_from_nearest_cut(h_elem, untouched)
    h_elem = np.where(np.isnan(h_elem), 0.0, h_elem)
    h_elem = np.ascontiguousarray(h_elem)

    cut_h = h_elem[cut_mask] if cut_mask.any() else np.zeros(1)
    stats = {
        "n_grits_total": len(frames),
        "n_grits_engaged": len(order),
        "n_surface_points_swept": int(n_pts),
        "time_samples_min": int(min(n_samples_used)) if n_samples_used else 0,
        "time_samples_max": int(max(n_samples_used)) if n_samples_used else 0,
        "depth_resolution_nm": float(params.depth_resolution_mm * 1e6),
        "facet_subdivision": params.facet_subdivision,
        "station_step_um": float(du * 1000.0),
        "max_depth_removed_um": float(removed.max() * 1000.0),
        "mean_depth_removed_um": float(removed.mean() * 1000.0),
        # over the elements a grit actually removes -- the only ones whose chip
        # thickness means anything
        "h_cut_min_um": float(cut_h.min() * 1000.0),
        "h_cut_max_um": float(cut_h.max() * 1000.0),
        "h_cut_mean_um": float(cut_h.mean() * 1000.0),
        "n_elements": int(h_elem.size),
        "n_elements_cut": int(cut_mask.sum()),
        "n_elements_never_cut": n_never,
        "never_cut_filled_from_nearest": bool(params.fill_never_cut),
        "fill_method": _FILL_METHOD[0],
        "removed_volume_mm3": float(removed.sum() * du * dz),
        # Can the mesh hold the cut at all? The project's own advice is 5-10
        # elements through the deepest cut; below about 1 there is no chip.
        "surface_layer_um": float((d_edges[1] - d_edges[0]) * 1000.0),
        "elements_through_deepest_cut": float(
            removed.max() / max(d_edges[1] - d_edges[0], 1e-30)),
    }
    if log:
        log("envelope: %d of %d grits cross the block, %s surface points, "
            "%s-%s time samples at %.2f nm depth resolution"
            % (stats["n_grits_engaged"], stats["n_grits_total"],
               format(stats["n_surface_points_swept"], ","),
               format(stats["time_samples_min"], ","),
               format(stats["time_samples_max"], ","),
               stats["depth_resolution_nm"]))
        log("          %s of %s elements are cut; their h is %.4f to %.4f nm; "
            "deepest groove %.4f um"
            % (format(stats["n_elements_cut"], ","),
               format(stats["n_elements"], ","),
               stats["h_cut_min_um"] * 1000.0, stats["h_cut_max_um"] * 1000.0,
               stats["max_depth_removed_um"]))

    return ChipEnvelope(h_elem=h_elem, depth_removed=removed,
                        n_grits_engaged=len(order),
                        grit_order=[gi for _t, gi, _w in order],
                        per_grit_h=per_grit, u_edges=u_edges, z_edges=z_edges,
                        depth_edges=np.asarray(d_edges), stats=stats,
                        cut=cut_mask)


_FILL_METHOD = [None]
"""Which branch _fill_from_nearest_cut actually took, for the report."""


def _fill_from_nearest_cut(h: np.ndarray, missing: np.ndarray) -> np.ndarray:
    """Give every uncut element the chip thickness of the nearest cut one.

    A Euclidean distance transform on the element grid, so the value
    propagates sideways across the face and downwards into the depth. That
    puts material under a brittle groove on the brittle law, which is the
    field it will actually see, and it keeps material beside a ductile scratch
    ductile.

    Falls back to nearest-along-depth only if scipy is unavailable, which keeps
    the module importable in a bare Colab runtime.
    """
    out = np.array(h, dtype=np.float64)
    if not missing.any():
        return out
    if missing.all():
        return np.zeros_like(out)
    try:
        from scipy import ndimage
        _d, idx = ndimage.distance_transform_edt(
            missing, return_distances=True, return_indices=True)
        out[missing] = out[tuple(i[missing] for i in idx)]
        _FILL_METHOD[0] = "scipy nearest-neighbour"
        return out
    except ImportError:
        # ImportError ONLY. A bare ``except Exception`` here silently swapped a
        # true nearest-neighbour fill for a depth-only fill on ANY failure, and
        # the two disagree by a factor of six in the resulting field statistics
        # (mean h 33.7 -> 5.6 nm, brittle-law count 5484 -> 914). verify_envelope
        # passes either way, so the swap would never have been noticed.
        #
        # Column-wise forward then backward fill through the depth.
        nl, nw, nd = out.shape
        for k in range(1, nd):
            m = np.isnan(out[:, :, k])
            out[:, :, k][m] = out[:, :, k - 1][m]
        for k in range(nd - 2, -1, -1):
            m = np.isnan(out[:, :, k])
            out[:, :, k][m] = out[:, :, k + 1][m]
        _FILL_METHOD[0] = "depth-only fill (scipy unavailable)"
        return out


def tip_paths(place: dict, motion: Optional[dict], wp, *, step_time_s: float,
              rotation_reversed: bool = False,
              params: Optional[EnvelopeParams] = None,
              paths: Optional[dict] = None, n: int = 400) -> dict:
    """Where each grit's outermost point goes: ``{gi: (n, 3) [t, u, depth]}``.

    Exactly the motion :func:`sweep_envelope` uses, on one point instead of a
    cloud, so the trajectory that gets drawn is the trajectory that was swept.
    Deriving the picture from its own copy of the kinematics is how a plot ends
    up disagreeing with the model it illustrates, and this project has been
    bitten by that before.
    """
    params = params or EnvelopeParams()
    paths = paths or {}
    omega = 0.0 if motion is None else float(motion["vr3"])
    v_r = 0.0 if motion is None else float(motion["radial_speed_mm_s"])
    r_ground = float(place["r_ground"])
    hl = wp.length_mm / 2.0
    out = {}
    for gi, v in enumerate(place["frames"]):
        if gi in paths:
            pt = np.asarray(paths[gi], dtype=np.float64)
            out[gi] = np.column_stack([pt[:, 0], pt[:, 1], pt[:, 3]])
            continue
        if omega == 0.0:
            continue
        tip = int(np.argmax(v[:, 0]))
        a0, b0 = float(v[tip, 0]), float(v[tip, 1])
        c = v.mean(axis=0)
        reach = float(np.hypot(v[:, 0], v[:, 1]).max())
        rad = float(np.linalg.norm(v - c, axis=1).max())
        win = _crossing_window(float(c[0]), float(c[1]), reach, omega, hl,
                               step_time_s, params.margin_factor * rad)
        if win is None:
            continue
        ts = np.linspace(win[0], win[1], n)
        al = omega * ts
        aa = a0 * np.cos(al) - b0 * np.sin(al) + v_r * ts
        bb = a0 * np.sin(al) + b0 * np.cos(al)
        out[gi] = np.column_stack([ts, bb, aa - r_ground])
    return out


# --------------------------------------------------------------------------
# element field -> nodal field
# --------------------------------------------------------------------------

def nodal_field(env: ChipEnvelope, wp) -> np.ndarray:
    """Chip thickness at every workpiece node, mm, in the deck's node order.

    Abaqus interpolates a nodal field variable to the integration point, and a
    C3D8R has one, at the centroid, where the interpolation is the mean of the
    eight nodes. Setting each node to the mean of the elements around it
    therefore hands the integration point back very nearly its own element's
    value: the two differ only by the field's curvature over one element, which
    is nanometres over 0.3 um.

    The one visible consequence is that an element straddling the transition
    gets a blended chip thickness and may fall on either side. The mesh already
    limits the transition to plus or minus one element, so this adds nothing
    that was not there.
    """
    nl, nw, nd = wp.divisions()
    if env.h_elem.shape != (nl, nw, nd):
        raise EnvelopeError("envelope shape %s does not match the block's "
                            "%s elements" % (env.h_elem.shape, (nl, nw, nd)))
    acc = np.zeros((nl + 1, nw + 1, nd + 1), dtype=np.float64)
    cnt = np.zeros((nl + 1, nw + 1, nd + 1), dtype=np.float64)
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                acc[di:di + nl, dj:dj + nw, dk:dk + nd] += env.h_elem
                cnt[di:di + nl, dj:dj + nw, dk:dk + nd] += 1.0
    nodal = acc / np.maximum(cnt, 1.0)
    # build_block_mesh numbers nodes as ((i*(nw+1) + j)*(nd+1) + k)
    return nodal.reshape(-1)
