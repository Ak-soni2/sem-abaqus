"""One figure that shows what the swept field decided.

Three panels, because three questions get asked of this model and each needs a
different view:

* **the groove** -- did the grits remove what you expected, and where;
* **the chip thickness against dc** -- which parts of the cut are on which side
  of the transition, on a log axis because h spans three decades;
* **the map** -- where on the ground face the ductile and brittle regions are,
  which is the picture that goes in the paper next to SDV13 from the run.

Kept out of the notebook cell so it can be exercised without a browser.
"""

from __future__ import annotations

from typing import Optional

import numpy as np


def preview_figure(env, dc_mm: float, wp, *, title: str = "",
                   figsize=(12.0, 8.5)):
    """Groove profile, h against dc, and the ductile/brittle map."""
    import matplotlib
    if matplotlib.get_backend().lower() not in ("agg", "module://ipykernel"
                                                ".pylab.backend_inline"):
        pass
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    nl, nw, nd = wp.divisions()
    u_c = 0.5 * (env.u_edges[:-1] + env.u_edges[1:]) * 1000.0     # um
    z_c = 0.5 * (env.z_edges[:-1] + env.z_edges[1:]) * 1000.0
    dc_nm = dc_mm * 1e6

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 1, height_ratios=(1.0, 1.0, 1.15), hspace=0.42)

    # -- 1. the groove ------------------------------------------------
    ax = fig.add_subplot(gs[0])
    prof = env.depth_removed.max(axis=1) * 1000.0
    mean = env.depth_removed.mean(axis=1) * 1000.0
    ax.plot(u_c, prof, lw=1.6, label="deepest across the face")
    ax.plot(u_c, mean, lw=1.0, ls="--", label="mean across the face")
    ax.invert_yaxis()
    ax.set_xlabel("station along the scratch, u  (um)")
    ax.set_ylabel("depth removed (um)")
    ax.set_title("the groove the grits cut" + (" - " + title if title else ""))
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)

    # -- 2. chip thickness against dc ---------------------------------
    ax = fig.add_subplot(gs[1])
    cut = env.cut if env.cut.size == env.h_elem.size else None
    hn = env.h_elem * 1e6
    if cut is not None and cut.any():
        # per station, the range of chip thickness over the cut elements
        lo = np.full(nl, np.nan)
        hi = np.full(nl, np.nan)
        for i in range(nl):
            v = hn[i][cut[i]]
            if v.size:
                lo[i], hi[i] = v.min(), v.max()
        ok = ~np.isnan(lo)
        if ok.any():
            ax.fill_between(u_c[ok], lo[ok], hi[ok], alpha=0.35,
                            label="range of h over the cut elements")
            ax.plot(u_c[ok], hi[ok], lw=1.2, label="deepest cut at that station")
    ax.axhline(dc_nm, color="crimson", lw=1.4,
               label="dc = %.3f nm" % dc_nm)
    ax.set_yscale("log")
    ax.set_xlabel("station along the scratch, u  (um)")
    ax.set_ylabel("undeformed chip thickness (nm)")
    ax.set_title("chip thickness against the critical depth  "
                 "(below the line is ductile)")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3, which="both")

    # -- 3. the map on the ground face --------------------------------
    ax = fig.add_subplot(gs[2])
    # 0 = never cut, 1 = cut ductilely, 2 = cut brittlely, taken on the
    # shallowest layer, which is the face the run will show.
    surf_h = env.h_elem[:, :, 0]
    surf_cut = (env.cut[:, :, 0] if env.cut.size == env.h_elem.size
                else np.ones_like(surf_h, dtype=bool))
    code = np.where(~surf_cut, 0, np.where(surf_h < dc_mm, 1, 2))
    cmap = ListedColormap(["#e8e8e8", "#2e7d32", "#c62828"])
    im = ax.imshow(code.T, origin="lower", aspect="auto", cmap=cmap,
                   vmin=-0.5, vmax=2.5,
                   extent=(u_c[0], u_c[-1], z_c[0], z_c[-1]))
    cb = fig.colorbar(im, ax=ax, ticks=(0, 1, 2), pad=0.02)
    cb.ax.set_yticklabels(["not cut", "ductile", "brittle"], fontsize=8)
    ax.set_xlabel("station along the scratch, u  (um)")
    ax.set_ylabel("across the face, z  (um)")
    ax.set_title("the ground face: which law each surface element runs")
    return fig


def field_slice_figure(env, dc_mm: float, wp, *, figsize=(11.0, 4.0)):
    """The chip-thickness field through the depth, at mid-face."""
    import matplotlib.pyplot as plt

    nl, nw, nd = wp.divisions()
    j = nw // 2
    u_c = 0.5 * (env.u_edges[:-1] + env.u_edges[1:]) * 1000.0
    d_c = 0.5 * (env.depth_edges[:-1] + env.depth_edges[1:]) * 1000.0
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.pcolormesh(u_c, d_c, (env.h_elem[:, j, :] * 1e6).T,
                       shading="nearest")
    ax.invert_yaxis()
    fig.colorbar(im, ax=ax, label="h (nm)", pad=0.02)
    ax.contour(u_c, d_c, (env.h_elem[:, j, :] * 1e6).T, levels=[dc_mm * 1e6],
               colors="crimson", linewidths=1.4)
    ax.set_xlabel("station along the scratch, u  (um)")
    ax.set_ylabel("depth into the block (um)")
    ax.set_title("chip thickness through the depth at mid-face; "
                 "the red contour is dc")
    return fig


def trajectory_figure(place, motion, wp, dc_mm: float, *, step_time_s: float,
                      rotation_reversed: bool = False, paths=None,
                      env=None, figsize=(12.0, 7.0)):
    """Where the abrasives went, and where that crosses dc.

    Top: depth below the original surface against station, one line per grit,
    with dc drawn across it. Everything under the dc line is ductile, everything
    over it is brittle, so the crossing IS the transition and you can read the
    depth of cut it happens at straight off the axis.

    Bottom: the same paths seen from above, so it is obvious which part of the
    face each grit sweeps and where they overlap.

    The paths come from :func:`semgrit_multi.envelope.tip_paths`, which is the
    sweep's own kinematics -- not a second copy of them.
    """
    import matplotlib.pyplot as plt

    from .envelope import tip_paths

    tp = tip_paths(place, motion, wp, step_time_s=step_time_s,
                   rotation_reversed=rotation_reversed, paths=paths)
    if not tp:
        raise ValueError("no grit moves across the block, so there is no "
                        "trajectory to draw")
    hl, hw = wp.length_mm / 2.0, wp.width_mm / 2.0
    dc_nm = dc_mm * 1e6

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=figsize, height_ratios=(2.0, 1.0),
        gridspec_kw=dict(hspace=0.35))

    deepest = 0.0
    for gi, p in sorted(tp.items()):
        u, d = p[:, 1] * 1000.0, p[:, 2] * 1000.0
        inside = (u >= -hl * 1000.0) & (u <= hl * 1000.0) & (d > 0)
        ax.plot(u, np.maximum(d, 0.0), lw=1.0, alpha=0.85,
                label=("grit %d" % gi) if len(tp) <= 8 else None)
        if inside.any():
            deepest = max(deepest, float(d[inside].max()))
    ax.axhline(dc_nm / 1000.0, color="crimson", lw=1.6,
               label="dc = %.3f nm" % dc_nm)
    ax.axvspan(-hl * 1000.0, hl * 1000.0, color="0.9", zorder=0,
               label="the workpiece")
    top = max(deepest, dc_nm / 1000.0) * 1.15 + 1e-9
    ax.axhspan(0, dc_nm / 1000.0, color="#2e7d32", alpha=0.13, zorder=0)
    ax.axhspan(dc_nm / 1000.0, top, color="#c62828", alpha=0.10, zorder=0)
    ax.text(hl * 1000.0 * 0.98, dc_nm / 1000.0 * 0.5, "ductile", ha="right",
            va="center", fontsize=9, color="#2e7d32")
    if top > dc_nm / 1000.0:
        ax.text(hl * 1000.0 * 0.98, (top + dc_nm / 1000.0) / 2.0, "brittle",
                ha="right", va="center", fontsize=9, color="#c62828")
    ax.set_ylim(0, top)
    ax.invert_yaxis()
    ax.set_xlabel("station along the scratch, u  (um)")
    ax.set_ylabel("tip depth below the original surface  (um)")
    ax.set_title("abrasive trajectories against the critical depth of cut  "
                 "(%d grit%s crossing)" % (len(tp), "" if len(tp) == 1 else "s"))
    if len(tp) <= 8:
        ax.legend(fontsize=8, loc="lower left", ncol=2)
    else:
        ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3)

    for gi, p in sorted(tp.items()):
        z = np.full(len(p), 0.0)
        fr = place["frames"][gi]
        z[:] = float(fr[int(np.argmax(fr[:, 0])), 2]) * 1000.0
        if paths and gi in paths:
            z = np.asarray(paths[gi])[:, 2] * 1000.0
        ax2.plot(p[:, 1] * 1000.0, z, lw=1.0, alpha=0.85)
    ax2.add_patch(plt.Rectangle((-hl * 1000.0, -hw * 1000.0), 2 * hl * 1000.0,
                                2 * hw * 1000.0, fill=False, lw=1.2,
                                color="0.3"))
    ax2.set_xlabel("station along the scratch, u  (um)")
    ax2.set_ylabel("across the face, z  (um)")
    ax2.set_title("the same paths from above, with the workpiece outline")
    ax2.grid(alpha=0.3)
    return fig


def dc_sweep_figure(depths_um, fractions, *, dc_nm: float, chosen_um=None,
                    figsize=(9.0, 4.2)):
    """Ductile share of the cut against depth of cut: how to choose ae.

    ``fractions`` is the ductile fraction of the CUT elements at each depth of
    cut, from a real sweep at each one. The useful reading is the knee: below it
    the pass is all ductile, above it all brittle, and the interesting
    experiments are in between.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize)
    d = np.asarray(depths_um, dtype=float)
    f = np.asarray(fractions, dtype=float) * 100.0
    ax.plot(d, f, "o-", lw=1.6)
    ax.axhline(50.0, color="0.6", ls=":", lw=1.0)
    if chosen_um is not None:
        ax.axvline(chosen_um, color="crimson", lw=1.4,
                   label="this deck: ae = %.4f um" % chosen_um)
        ax.legend(fontsize=8)
    ax.set_xlabel("depth of cut  (um)")
    ax.set_ylabel("ductile share of the cut elements  (%)")
    ax.set_title("where the transition sits, against depth of cut  "
                 "(dc = %.3f nm)" % dc_nm)
    ax.set_ylim(-3, 103)
    ax.grid(alpha=0.3)
    return fig


def trajectory_check_figure(traj, wp, dc_mm: float, overlay=None,
                            figsize=(11.0, 6.0)):
    """A measured trajectory, as read in. Look at this before trusting it.

    Left, or top: the path itself against dc, with the block outline, so a units
    mistake is obvious -- a profile given in microns instead of millimetres
    lands a thousand times too deep and will not fit on the axis.
    Right: the traced overlay, when the path came from an image.
    """
    import matplotlib.pyplot as plt

    n = 2 if overlay is not None else 1
    fig, axes = plt.subplots(1, n, figsize=figsize)
    ax = axes[0] if n == 2 else axes
    hl, hw = wp.length_mm / 2.0, wp.width_mm / 2.0
    ax.plot(traj.u * 1000.0, traj.depth * 1000.0, lw=1.4)
    ax.axhline(dc_mm * 1e6 / 1000.0, color="crimson", lw=1.4,
               label="dc = %.3f nm" % (dc_mm * 1e6))
    ax.axvspan(-hl * 1000.0, hl * 1000.0, color="0.9", zorder=0,
               label="the workpiece")
    ax.invert_yaxis()
    ax.set_xlabel("station u  (um)")
    ax.set_ylabel("depth  (um)")
    ax.set_title("the measured trajectory, as read in")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    if overlay is not None:
        axes[1].imshow(overlay[:, :, ::-1] if overlay.ndim == 3 else overlay)
        axes[1].set_title("traced boundary (red) and the surface line (green)")
        axes[1].set_axis_off()
    return fig
