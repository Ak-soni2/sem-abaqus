"""Replay a measured abrasive path instead of the ideal wheel kinematics.

After a real single-grit experiment you have the groove, and the groove *is* the
measurement. This module turns it into a trajectory the sweep can follow, so the
simulated abrasive goes exactly where the real one went.

Three ways in, all producing the same thing:

``from_points``
    You already have coordinates. Give them as ``(u, depth)`` or
    ``(u, z, depth)`` or ``(t, u, z, depth)``, in mm, and they are used as they
    stand.
``from_csv``
    The same, from a file, with the column names or indices named.
``from_profile_image``
    A scaled image of the groove -- a cross-section, or a top view with a depth
    scale -- traced column by column into a profile. This one is a convenience,
    not a measurement: it traces one boundary with a threshold, so **always look
    at the overlay it returns before believing it**.

Coordinates
-----------
The sweep works in the block frame: ``u`` along the scratch, ``z`` across the
face, ``depth`` into the material from the original ground surface, all in mm
and all zero at the block centre / original surface. A profilometer trace of a
groove is exactly that once it is scaled, with the reference line taken as the
uncut surface.

What "trajectory" means here
----------------------------
The path is the position of the grit's **outermost point** -- its tip. The grit
is translated so that point follows the path, keeping its measured orientation.
It is not rotated along the path: a groove tells you where the abrasive was, not
how it was turned, and inventing a rotation would put shape into the result that
was never measured.

Timing is taken from the path's own ``t`` column when it has one, and otherwise
spread evenly over the step. For a single pass only the order matters, and even
spacing is the honest default.

What this cannot do
-------------------
It replays a kinematic path. It does not make the Abaqus job follow that path:
the wheel is still driven by the velocity boundary condition the deck writes, so
the *simulation* moves the grit on the ideal arc. What the trajectory changes is
the chip-thickness field the constitutive switch reads, which is what decides
ductile against brittle. If you need the simulated grit to physically follow a
measured path as well, that is a boundary-condition problem (``*Amplitude`` on
the reference node, or VDISP) and a different job -- ``deck_amplitudes`` below
writes the tables for it, but wiring them in is left to you, deliberately,
because it changes the mechanics and should be a decision rather than a default.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np


class TrajectoryError(RuntimeError):
    pass


@dataclass
class Trajectory:
    """A measured tip path in the block frame. Columns t, u, z, depth; mm and s."""

    samples: np.ndarray                 # (n, 4) [t, u, z, depth]
    source: str = ""
    notes: list = field(default_factory=list)

    def __post_init__(self):
        a = np.asarray(self.samples, dtype=np.float64)
        if a.ndim != 2 or a.shape[1] != 4:
            raise TrajectoryError("samples must be (n, 4): t, u, z, depth")
        if len(a) < 2:
            raise TrajectoryError("a trajectory needs at least two samples")
        if not np.isfinite(a).all():
            raise TrajectoryError("the trajectory carries non-finite values")
        if np.any(np.diff(a[:, 0]) < 0):
            raise TrajectoryError("the time column must be non-decreasing")
        self.samples = a

    # -- shape -----------------------------------------------------------
    @property
    def n(self) -> int:
        return len(self.samples)

    @property
    def t(self) -> np.ndarray:
        return self.samples[:, 0]

    @property
    def u(self) -> np.ndarray:
        return self.samples[:, 1]

    @property
    def z(self) -> np.ndarray:
        return self.samples[:, 2]

    @property
    def depth(self) -> np.ndarray:
        return self.samples[:, 3]

    def summary(self) -> dict:
        return {
            "source": self.source,
            "n_samples": self.n,
            "t_s": (float(self.t.min()), float(self.t.max())),
            "u_mm": (float(self.u.min()), float(self.u.max())),
            "z_mm": (float(self.z.min()), float(self.z.max())),
            "depth_um": (float(self.depth.min() * 1000.0),
                         float(self.depth.max() * 1000.0)),
            "length_mm": float(np.abs(np.diff(self.u)).sum()),
            "notes": list(self.notes),
        }

    # -- conditioning ----------------------------------------------------
    def retimed(self, t0: float, t1: float) -> "Trajectory":
        """The same path, spread evenly over ``[t0, t1]``."""
        if t1 <= t0:
            raise TrajectoryError("t1 must exceed t0")
        s = self.samples.copy()
        s[:, 0] = np.linspace(t0, t1, len(s))
        return Trajectory(s, self.source, self.notes + [
            "retimed to %.6e..%.6e s" % (t0, t1)])

    def resampled(self, n: int) -> "Trajectory":
        """Linearly resampled to ``n`` points, for a finer or cheaper sweep."""
        if n < 2:
            raise TrajectoryError("n must be at least 2")
        t = self.t
        tt = np.linspace(t[0], t[-1], n)
        out = np.column_stack([tt] + [np.interp(tt, t, self.samples[:, k])
                                      for k in (1, 2, 3)])
        return Trajectory(out, self.source, self.notes + ["resampled to %d" % n])

    def clipped_to_block(self, wp) -> "Trajectory":
        """Drop samples outside the block's footprint, keeping the order."""
        hl, hw = wp.length_mm / 2.0, wp.width_mm / 2.0
        keep = ((self.u >= -hl) & (self.u <= hl)
                & (self.z >= -hw) & (self.z <= hw))
        if keep.sum() < 2:
            raise TrajectoryError(
                "only %d of %d samples lie inside the %.4f x %.4f mm block. "
                "Check the units: this module wants millimetres."
                % (int(keep.sum()), self.n, wp.length_mm, wp.width_mm))
        return Trajectory(self.samples[keep], self.source,
                          self.notes + ["clipped to the block footprint"])


# --------------------------------------------------------------------------
# loaders
# --------------------------------------------------------------------------

def from_points(points, *, columns: str = "auto", scale_mm: float = 1.0,
                depth_sign: float = 1.0, source: str = "points",
                t0: float = 0.0, t1: float = 1.0) -> Trajectory:
    """Coordinates in, trajectory out.

    ``columns`` is ``"u,depth"``, ``"u,z,depth"``, ``"t,u,z,depth"`` or
    ``"auto"``, which picks by width. ``scale_mm`` multiplies the position
    columns, so a table in microns is ``scale_mm=1e-3``. ``depth_sign = -1``
    flips a profile that measures the groove as negative, which most
    profilometers do.
    """
    a = np.asarray(points, dtype=np.float64)
    if a.ndim != 2:
        raise TrajectoryError("points must be a 2-D array")
    if columns == "auto":
        columns = {2: "u,depth", 3: "u,z,depth", 4: "t,u,z,depth"}.get(
            a.shape[1], "")
        if not columns:
            raise TrajectoryError(
                "cannot guess the columns of a %d-wide table; name them"
                % a.shape[1])
    names = [c.strip().lower() for c in columns.split(",")]
    if a.shape[1] != len(names):
        raise TrajectoryError("columns names %d fields, the table has %d"
                              % (len(names), a.shape[1]))
    idx = {nm: i for i, nm in enumerate(names)}
    for req in ("u", "depth"):
        if req not in idx:
            raise TrajectoryError("columns must include %r" % req)
    n = len(a)
    u = a[:, idx["u"]] * scale_mm
    d = a[:, idx["depth"]] * scale_mm * depth_sign
    z = (a[:, idx["z"]] * scale_mm if "z" in idx else np.zeros(n))
    t = (a[:, idx["t"]] if "t" in idx else np.linspace(t0, t1, n))
    notes = []
    if scale_mm != 1.0:
        notes.append("positions scaled by %g into mm" % scale_mm)
    if depth_sign < 0:
        notes.append("depth sign flipped: the groove was given as negative")
    if "t" not in idx:
        notes.append("no time column, so the path is spread evenly over the step")
    if (d < 0).any():
        notes.append("WARNING: %d samples have negative depth, i.e. above the "
                     "uncut surface. They will not cut." % int((d < 0).sum()))
    return Trajectory(np.column_stack([t, u, z, d]), source, notes)


def from_csv(path: str, *, columns: str = "auto", scale_mm: float = 1.0,
             depth_sign: float = 1.0, delimiter: Optional[str] = None,
             skip_header: Optional[int] = None) -> Trajectory:
    """A trajectory from a text table. Header lines are detected and skipped."""
    if not os.path.exists(path):
        raise TrajectoryError("no such file: %s" % path)
    with open(path, encoding="utf-8-sig") as fh:
        raw = [ln for ln in fh.read().splitlines() if ln.strip()]
    if not raw:
        raise TrajectoryError("%s is empty" % path)
    if delimiter is None:
        first = raw[0]
        delimiter = ("\t" if "\t" in first else
                     ";" if ";" in first else
                     "," if "," in first else None)
    if skip_header is None:
        skip_header = 0
        for ln in raw:
            try:
                [float(x) for x in
                 (ln.split(delimiter) if delimiter else ln.split())
                 if x.strip()]
                break
            except ValueError:
                skip_header += 1
    rows = []
    for ln in raw[skip_header:]:
        f = [x for x in (ln.split(delimiter) if delimiter else ln.split())
             if x.strip()]
        try:
            rows.append([float(x) for x in f])
        except ValueError:
            continue
    if len(rows) < 2:
        raise TrajectoryError(
            "%s gave %d usable numeric rows. Check the delimiter and the "
            "header." % (path, len(rows)))
    width = min(len(r) for r in rows)
    arr = np.array([r[:width] for r in rows], dtype=np.float64)
    return from_points(arr, columns=columns, scale_mm=scale_mm,
                       depth_sign=depth_sign,
                       source="%s (%d rows, %d columns)"
                              % (os.path.basename(path), len(arr), width))


def from_profile_image(path: str, *, mm_per_px_x: float, mm_per_px_y: float,
                       surface_row: Optional[int] = None,
                       dark_is_material: bool = True,
                       threshold: Optional[float] = None,
                       smooth_px: int = 3,
                       source: str = "") -> tuple:
    """Trace a groove profile out of a scaled image.

    Returns ``(Trajectory, overlay)`` where ``overlay`` is the image with the
    traced curve drawn on it. **Look at the overlay.** This is one threshold and
    one boundary per column; it is a convenience for getting a shape in quickly,
    not a measurement, and it cannot tell a groove from a scratch, a shadow or a
    speck of dust.

    ``mm_per_px_x`` and ``mm_per_px_y`` are the scales along and into the
    groove -- from the image's own scale bar, in the units this project uses
    everywhere else. ``surface_row`` is the pixel row of the uncut surface;
    without it, the shallowest traced row is taken as the surface, which assumes
    the trace starts outside the groove.
    """
    try:
        import cv2
    except ImportError:  # pragma: no cover
        raise TrajectoryError("tracing an image needs opencv "
                              "(pip install opencv-python-headless)")
    if not os.path.exists(path):
        raise TrajectoryError("no such image: %s" % path)
    img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise TrajectoryError("could not read %s as an image" % path)
    work = img if dark_is_material else 255 - img
    if threshold is None:
        thr, mask = cv2.threshold(work, 0, 255,
                                  cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    else:
        thr = float(threshold)
        _r, mask = cv2.threshold(work, thr, 255, cv2.THRESH_BINARY_INV)
    if smooth_px >= 3:
        k = smooth_px | 1
        mask = cv2.medianBlur(mask, min(k, 9))

    h, w = mask.shape
    rows = np.full(w, -1, dtype=np.int64)
    for x in range(w):
        col = np.nonzero(mask[:, x])[0]
        if col.size:
            rows[x] = int(col[0])          # topmost material pixel
    ok = rows >= 0
    if ok.sum() < 2:
        raise TrajectoryError(
            "the threshold found material in only %d of %d columns. Pass "
            "threshold=, or set dark_is_material=False." % (int(ok.sum()), w))
    xs = np.nonzero(ok)[0]
    ys = rows[ok].astype(np.float64)
    base = float(ys.min()) if surface_row is None else float(surface_row)

    u = (xs - xs.mean()) * mm_per_px_x
    depth = (ys - base) * mm_per_px_y
    traj = from_points(np.column_stack([u, depth]), columns="u,depth",
                       source=source or "%s (%d columns traced)"
                                        % (os.path.basename(path), len(xs)))
    traj.notes.append("threshold %.1f, surface row %.1f, %s is material"
                      % (thr, base, "dark" if dark_is_material else "light"))
    traj.notes.append("CHECK THE OVERLAY: this is one traced boundary, not a "
                      "measurement")

    overlay = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    for x, y in zip(xs, ys.astype(np.int64)):
        cv2.circle(overlay, (int(x), int(y)), 1, (0, 0, 255), -1)
    cv2.line(overlay, (0, int(base)), (w - 1, int(base)), (0, 200, 0), 1)
    return traj, overlay


# --------------------------------------------------------------------------
# using it
# --------------------------------------------------------------------------

def sweep_trajectory(place: dict, wp, traj: Trajectory, *, grit: int = 0,
                     step_time_s: float, params=None, log=None):
    """Chip-thickness field from a measured path, for one grit.

    Returns the same :class:`semgrit_multi.envelope.ChipEnvelope` the ideal
    sweep returns, so the field can be injected, gated and plotted by exactly
    the same code.
    """
    from .envelope import sweep_envelope

    t = traj.retimed(0.0, step_time_s) if (traj.t.max() > step_time_s
                                          or traj.t.max() <= 0) else traj
    t = t.clipped_to_block(wp)
    if grit >= len(place["frames"]):
        raise TrajectoryError("grit %d does not exist; the deck placed %d"
                              % (grit, len(place["frames"])))
    # A rotation of zero and an infeed of zero: the path supplies the motion.
    # sweep_envelope still refuses a non-turning wheel for the grits it has to
    # move itself, so hand it a nominal sense that matches its own convention
    # and let the path override this grit.
    motion = {"vr3": -1.0, "radial_speed_mm_s": 0.0, "omega_rad_s": 1.0}
    return sweep_envelope(
        {**place, "frames": [place["frames"][grit]],
         "faces": [place["faces"][grit]]},
        motion, wp, step_time_s=step_time_s, rotation_reversed=False,
        params=params, paths={0: t.samples}, log=log)


def deck_amplitudes(traj: Trajectory, theta_c: float, r_ground: float) -> str:
    """``*Amplitude`` tables that would drive a reference node along the path.

    Provided because it is the obvious next question and the arithmetic is
    fiddly, NOT wired into any deck: making the simulated wheel follow a
    measured path replaces the rotation-plus-infeed boundary condition with a
    prescribed displacement, which changes the mechanics of the run. That should
    be a decision, not a default.
    """
    e_r = (math.cos(theta_c), math.sin(theta_c))
    e_t = (-math.sin(theta_c), math.cos(theta_c))
    L = ["** Reference-node displacement that follows the measured path.",
         "** Feed these to *Boundary, type=DISPLACEMENT, amplitude=...",
         "** NOT used by any deck this project writes; see the docstring."]
    for comp, name in ((1, "TRAJ_U1"), (2, "TRAJ_U2")):
        L.append("*Amplitude, name=%s, definition=TABULAR" % name)
        vals = []
        for t, u, _z, d in traj.samples:
            # tip at (r_ground + d) along e_r and u along e_t; the wheel has to
            # move by the same amount for the tip to get there
            x = (r_ground + d) * e_r[comp - 1] + u * e_t[comp - 1]
            x0 = r_ground * e_r[comp - 1]
            vals.append("%.9e, %.9e" % (t, x - x0))
        for i in range(0, len(vals), 4):
            L.append(", ".join(vals[i:i + 4]))
    return "\n".join(L) + "\n"
