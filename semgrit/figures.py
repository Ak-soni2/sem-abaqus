"""Show what the pipeline did, at every stage, as figures.

Every other module in this package computes; this one is the only one whose job
is to make the computation *visible*. It exists because the pipeline was
effectively unwatchable: an SEM micrograph went in, two lines of text came out
("45 grains -> 27 solids"), and the eleven image-processing steps in between --
the thresholding, the distance transform, the seeding, the watershed, the
split-retention argument, the outline extraction, the 3-D loft -- were locals
that vanished on return. A reader could accept the grain count or reject it, and
had nothing in between.

So each function here takes the data one stage of the pipeline actually produced
and draws it. Nothing is recomputed and nothing is approximated for display: the
segmentation panels draw the very arrays ``segment_grains`` used, and the grain
panels draw the very triangles that are written into the ``.inp``. A figure that
redrew its own idea of the geometry could agree with a picture and disagree with
the deck, which is the one thing a verification figure must not do.

Colours follow the repo's Okabe-Ito convention (see ``REPOST/plots.py``):
ductile blue, brittle vermillion, both legible to deuteranopes and in greyscale.

    from semgrit.figures import segmentation_stages, grain_gallery
    fig = segmentation_stages(per_image[0])

All functions return a matplotlib ``Figure`` and draw nothing to screen, so the
caller decides between ``plt.show()`` and ``fig.savefig(...)``.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

# Okabe-Ito, matching REPOST/plots.py so one project speaks one visual language.
C_DUCTILE = "#0072B2"
C_BRITTLE = "#D55E00"
C_KEEP = "#009E73"
C_DROP = "#D55E00"
C_WARN = "#B00020"
C_GREY = "#666666"
C_ACCENT = "#CC79A7"
C_SEED = "#F0E442"

# A qualitative cycle for label images. Randomised once with a fixed seed so
# neighbouring grains differ but the figure is reproducible run to run.
_LABEL_SEED = 20260823


def _plt():
    """Import pyplot late, and apply house style once.

    Late because importing this module must not drag matplotlib into an Abaqus
    kernel that has none -- that is exactly the bug REPOST/plots.py exists to
    work around.
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 170, "font.size": 10,
        "axes.titlesize": 11, "axes.labelsize": 10, "axes.grid": False,
        "legend.frameon": False, "figure.constrained_layout.use": True,
    })
    return plt


def _label_cmap(n: int):
    """A shuffled qualitative colormap with label 0 transparent."""
    from matplotlib.colors import ListedColormap
    import matplotlib.pyplot as plt
    base = plt.get_cmap("tab20")(np.linspace(0, 1, 20))
    reps = int(math.ceil(max(n, 1) / 20.0))
    cols = np.tile(base, (reps, 1))[:max(n, 1)]
    rng = np.random.default_rng(_LABEL_SEED)
    rng.shuffle(cols)
    return ListedColormap(np.vstack([[0, 0, 0, 0], cols]))


def _cbar(fig, ax, im, label=None, *, inside=False):
    """A colorbar that does not steal width from its own panel.

    ``fig.colorbar(..., ax=ax)`` shrinks the host axes, so in a grid where only
    some panels have one the images come out different sizes and the figure
    reads as sloppy. An inset keeps every panel identical. ``inside=True`` puts
    it over the image instead of beside it, for the panel whose right-hand
    neighbour would otherwise be crowded.
    """
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    if inside:
        cax = inset_axes(ax, width="3%", height="60%", loc="lower right",
                         bbox_to_anchor=(-0.02, 0.06, 1, 1),
                         bbox_transform=ax.transAxes, borderpad=0)
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=7, colors="black")
        cb.outline.set_edgecolor("black")
    else:
        cax = inset_axes(ax, width="3%", height="80%", loc="center left",
                         bbox_to_anchor=(1.015, 0.0, 1, 1),
                         bbox_transform=ax.transAxes, borderpad=0)
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=7.5)
    if label:
        cb.set_label(label, fontsize=8)
    return cb


def _show_labels(ax, labels, title, *, background=None):
    """Draw a label image over an optional greyscale backdrop."""
    n = int(labels.max())
    if background is not None:
        ax.imshow(background, cmap="gray", interpolation="nearest")
        alpha = np.where(labels > 0, 0.75, 0.0)
    else:
        alpha = None
    ax.imshow(labels, cmap=_label_cmap(n), interpolation="nearest",
              vmin=0, vmax=max(n, 1), alpha=alpha)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    return n


class _TextPanel:
    """A key/value panel on a blank axes, laid out top-down.

    Hand-rolled ``y -= 0.07`` layout was used in three figures here and drifted
    off the bottom of the panel whenever a row was added. This keeps one cursor
    and one line height, so a panel cannot silently overflow.
    """

    def __init__(self, ax, title=None, *, line=0.068, key_x=0.0, val_x=0.55):
        self.ax = ax
        self.y = 1.0
        self.line = line
        self.key_x, self.val_x = key_x, val_x
        ax.axis("off")
        if title:
            ax.text(0.0, self.y, title, fontsize=12, fontweight="bold",
                    va="top", transform=ax.transAxes)
            self.y -= line * 1.55

    def row(self, key, value, *, colour="black"):
        self.ax.text(self.key_x, self.y, key, fontsize=9.5, color=C_GREY,
                     va="top", transform=self.ax.transAxes)
        self.ax.text(self.val_x, self.y, str(value), fontsize=9.5, va="top",
                     family="monospace", color=colour,
                     transform=self.ax.transAxes)
        self.y -= self.line

    def note(self, text, colour="black", *, size=9.5, bold=False, wrap=0):
        import textwrap
        for chunk in (textwrap.wrap(text, wrap) if wrap else [text]):
            self.ax.text(self.key_x, self.y, chunk, fontsize=size, color=colour,
                         va="top", transform=self.ax.transAxes,
                         fontweight="bold" if bold else "normal")
            self.y -= self.line * (size / 9.5)

    def gap(self, n=0.5):
        self.y -= self.line * n


def _scalebar(ax, pixel_size_um, width_px, *, colour="white"):
    """A 1-2-5 scale bar sized to about a fifth of the frame.

    Drawn from the calibration the pipeline actually used, so a wrong pixel size
    is visible in the picture rather than buried in a report field.
    """
    target_um = width_px * pixel_size_um / 5.0
    if target_um <= 0:
        return
    exp = math.floor(math.log10(target_um))
    for mult in (1, 2, 5, 10):
        nice = mult * 10.0 ** exp
        if nice >= target_um:
            break
    n_px = nice / pixel_size_um
    x0 = width_px * 0.04
    y0 = ax.get_ylim()[0] * 0.95 if ax.get_ylim()[0] > 0 else 0
    ax.plot([x0, x0 + n_px], [y0, y0], color=colour, lw=3,
            solid_capstyle="butt")
    label = ("%g um" % nice) if nice >= 1 else ("%g nm" % (nice * 1000))
    ax.text(x0 + n_px / 2.0, y0, label, color=colour, ha="center", va="bottom",
            fontsize=9, fontweight="bold")


# ---------------------------------------------------------------------------
# 1. calibration
# ---------------------------------------------------------------------------

def calibration(rec: dict, figsize=(13, 5.2)):
    """The micrograph, the databar it was cropped from, and the scale check.

    The scale is the single number every later result is proportional to -- the
    original notebook had it wrong by 15-30x and every grain with it -- so it
    gets a panel of its own rather than a line of text.
    """
    plt = _plt()
    sem = rec["sem"]
    fig, axes = plt.subplots(1, 3, figsize=figsize,
                             gridspec_kw={"width_ratios": [1.25, 1.25, 1.15]})

    ax = axes[0]
    ax.imshow(sem.full_intensity, cmap="gray", interpolation="nearest")
    # Everything here is in DATA coordinates (image rows/columns). The databar
    # row is a row index, so mixing it with an axes fraction -- which also runs
    # the other way for an image -- put the label at the top of the frame
    # pointing at a line along the bottom.
    ax.axhline(sem.databar_top, color=C_WARN, lw=1.6, ls="--")
    ax.annotate("databar detected and cropped",
                xy=(sem.full_intensity.shape[1] * 0.5, sem.databar_top),
                xytext=(0, 8), textcoords="offset points", color=C_WARN,
                fontsize=9, fontweight="bold", ha="center", va="bottom",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                          alpha=0.75))
    bar = sem.scale_bar
    if bar is not None:
        y0, x0, y1, x1 = bar.bbox
        ax.add_patch(plt.Rectangle((x0 - 4, y0 - 4), x1 - x0 + 8,
                                   max(y1 - y0, 3) + 8,
                                   fill=False, edgecolor=C_SEED, lw=2))
        # Label above the bar and inside the micrograph, where there is black to
        # read it against -- beside it lands on the Zeiss text.
        ax.annotate("scale bar: %.0f px" % bar.length_px,
                    xy=(x0, y0), xytext=(0, 34), textcoords="offset points",
                    color=C_SEED, fontsize=9, fontweight="bold", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.22", fc="black", ec=C_SEED,
                              lw=0.8, alpha=0.65))
    ax.set_title("as acquired  (%d x %d px)" % sem.full_intensity.shape[::-1])
    ax.set_xticks([]); ax.set_yticks([])

    ax = axes[1]
    ax.imshow(sem.intensity, cmap="gray", interpolation="nearest")
    _scalebar(ax, sem.pixel_size_um, sem.intensity.shape[1])
    ax.set_title("measured region  (%.1f x %.1f um)"
                 % (sem.width_um, sem.height_um))
    ax.set_xticks([]); ax.set_yticks([])

    # The cross-check, as a picture: two independent estimates of one number.
    ax = axes[2]
    panel = _TextPanel(ax, "CALIBRATION")
    panel.row("image", rec["name"])
    panel.row("magnification", sem.magnification or "-")
    panel.row("pixel size", "%.5f um/px" % sem.pixel_size_um)
    panel.row("source", sem.pixel_size_source)
    if bar is not None and bar.pixel_size_um:
        panel.gap()
        panel.row("scale bar", "%.0f px" % bar.length_px)
        panel.row("bar implies", "%.4f um" % (bar.implied_um or 0.0))
        panel.row("snapped to", "%g um" % (bar.snapped_label_um or 0.0))
        panel.row("bar pixel size", "%.5f um/px" % bar.pixel_size_um)
    agree = sem.scalebar_agreement
    panel.gap()
    if agree is None:
        panel.note("no scale-bar cross-check available", C_GREY)
    else:
        ok = abs(agree) <= 0.05
        col = C_KEEP if ok else C_WARN
        panel.note("metadata vs scale bar   %+.2f %%" % (100 * agree), col,
                   size=11, bold=True)
        panel.note("agrees within 5 %: calibration trusted" if ok else
                   "DISAGREES: the run would have stopped", col)
    if sem.warnings:
        panel.gap()
        for w in sem.warnings[:4]:
            panel.note("! " + w, C_WARN, size=8.2, wrap=46)

    fig.suptitle("1 - Calibration: the number every later result scales with",
                 fontsize=13, fontweight="bold")
    return fig


# ---------------------------------------------------------------------------
# 2. segmentation, every stage
# ---------------------------------------------------------------------------

def segmentation_stages(rec: dict, figsize=(15, 11)):
    """All eleven image-processing stages, in the order they were computed.

    Draws the arrays ``segment_grains`` actually used (captured via its
    ``stages`` argument), not a re-run -- so what is shown is what happened.
    """
    plt = _plt()
    st = rec.get("stages") or {}
    if not st:
        raise ValueError("no captured stages: call measure_images(..., "
                         "keep_stages=True)")
    seg, sem = rec["seg"], rec["sem"]
    ps = seg.pixel_size_um

    fig, axes = plt.subplots(3, 4, figsize=figsize)
    A = axes.ravel()

    A[0].imshow(st["raw"], cmap="gray", interpolation="nearest")
    A[0].set_title("1  raw micrograph")
    _scalebar(A[0], ps, st["raw"].shape[1])

    A[1].imshow(st["denoised"], cmap="gray", interpolation="nearest")
    A[1].set_title("2  denoised  (median %.2f um)"
                   % seg.params.median_um)

    ax = A[2]
    ax.hist(st["denoised"].ravel(), bins=128, color=C_GREY)
    for t in seg.threshold_values:
        ax.axvline(t, color=C_BRITTLE, lw=1.6)
        ax.annotate("%.0f" % t, xy=(t, 0.97), xycoords=("data", "axes fraction"),
                    xytext=(3, 0), textcoords="offset points", color=C_BRITTLE,
                    fontsize=9, fontweight="bold", ha="left", va="top")
    ax.set_title("3  %s threshold" % seg.params.threshold_method)
    ax.set_xlabel("grey level   (foreground = above the upper threshold)",
                  fontsize=8.6)
    ax.set_yticks([])
    ax.set_yscale("log")

    A[3].imshow(st["threshold_raw"], cmap="gray", interpolation="nearest")
    A[3].set_title("4  foreground, raw  (%.1f %% of frame)"
                   % (100.0 * st["threshold_raw"].mean()))

    # What morphology changed, as a difference map rather than a second picture
    # that looks the same: added green, removed vermillion.
    before, after = st["threshold_raw"], st["foreground"]
    added, removed = after & ~before, before & ~after
    # A few hundred changed pixels in a 700k-pixel frame is invisible at figure
    # scale, so the changes are dilated purely for display. The caption carries
    # the true counts; the picture is there to say *where* they are.
    import cv2 as _cv
    k = _cv.getStructuringElement(_cv.MORPH_ELLIPSE, (7, 7))
    add_v = _cv.dilate(added.astype(np.uint8), k).astype(bool)
    rem_v = _cv.dilate(removed.astype(np.uint8), k).astype(bool)
    diff = np.zeros(before.shape + (3,), dtype=float)
    diff[after & before] = 0.82
    diff[rem_v] = np.array([0.84, 0.37, 0.0])
    diff[add_v] = np.array([0.0, 0.62, 0.45])
    A[4].imshow(diff, interpolation="nearest")
    A[4].set_title("5  close / open / fill   (+%d  -%d px)"
                   % (int(added.sum()), int(removed.sum())))

    im = A[5].imshow(st["distance_um"], cmap="magma", interpolation="nearest")
    A[5].set_title("6  distance transform  (max %.2f um)"
                   % float(st["distance_um"].max()))
    _cbar(fig, A[5], im, "um")

    seeds = st["seeds"]
    A[6].imshow(st["denoised"], cmap="gray", interpolation="nearest")
    ys, xs = np.nonzero(seeds)
    A[6].scatter(xs, ys, s=1.2, color=C_SEED, marker="s", linewidths=0)
    A[6].set_title("7  h-maxima seeds  (%d, h = %.2f um)"
                   % (seg.n_seeds, seg.params.h_maxima_um))

    im = A[7].imshow(st["gradient"], cmap="cividis", interpolation="nearest")
    A[7].set_title("8  intensity gradient (Sobel)")
    _cbar(fig, A[7], im)

    im = A[8].imshow(st["elevation"], cmap="terrain", interpolation="nearest")
    A[8].set_title("9  elevation = -dist + %.1f x grad"
                   % seg.params.gradient_weight)
    _cbar(fig, A[8], im, inside=True)

    n_raw = _show_labels(A[9], st["watershed_raw"],
                         "10  watershed  (%d regions, deliberately over-split)"
                         % int(st["watershed_raw"].max()))

    # The split-retention argument, drawn as the scatter the rule is a line in.
    ev = st.get("boundary_evidence") or {}
    ax = A[10]
    if ev:
        es = np.array([v["edge_strength"] for v in ev.values()])
        nr = np.array([v["neck_ratio"] for v in ev.values()])
        kp = np.array([v["kept"] for v in ev.values()])
        ax.scatter(es[kp], nr[kp], s=26, color=C_KEEP, label="kept (%d)" % kp.sum(),
                   edgecolor="black", linewidth=0.4, zorder=3)
        ax.scatter(es[~kp], nr[~kp], s=26, color=C_DROP,
                   label="merged back (%d)" % (~kp).sum(),
                   edgecolor="black", linewidth=0.4, zorder=3)
        ax.axvline(seg.params.min_edge_strength, color=C_GREY, ls="--", lw=1.2)
        ax.axhline(seg.params.min_neck_ratio, color=C_GREY, ls="--", lw=1.2)
        ax.set_xlabel("edge strength  (x mean gradient)")
        ax.set_ylabel("neck ratio")
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(alpha=0.25)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "no shared boundaries", ha="center", va="center",
                color=C_GREY, transform=ax.transAxes)
    ax.set_title("11  which splits the image supports")

    n_fin = _show_labels(A[11], seg.labels,
                         "12  final: %d grains  (%d merged, %d too small)"
                         % (seg.n_grains, seg.rejected.get("merged_splits", 0),
                            seg.rejected.get("too_small", 0)),
                         background=st["denoised"])

    for ax in A:
        if ax.images or ax is A[6]:
            ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle("2 - Segmentation, stage by stage:  %d seeds -> %d watershed "
                 "regions -> %d grains" % (seg.n_seeds, n_raw, n_fin),
                 fontsize=13, fontweight="bold")
    return fig


def segmentation_overlay(rec: dict, figsize=(14, 6.4), max_labels=400):
    """The result on top of the image: outlines, ids, and what was excluded.

    Border-truncated grains are the ones the original notebook measured as
    whole, which is why they are drawn in a different colour rather than simply
    dropped without comment.
    """
    plt = _plt()
    seg, sem = rec["seg"], rec["sem"]
    grains = rec["grains"]
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    for ax in axes:
        ax.imshow(sem.intensity, cmap="gray", interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])

    # left: every region, coloured by whether it is used
    from matplotlib.patches import Patch
    import matplotlib.patheffects as pe
    import cv2
    n_int = n_bor = 0
    for g in grains[:max_labels]:
        mask = (seg.labels == g.label).astype(np.uint8)
        cs, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        col = C_BRITTLE if g.touches_border else C_KEEP
        if g.touches_border:
            n_bor += 1
        else:
            n_int += 1
        for c in cs:
            axes[0].plot(c[:, 0, 0], c[:, 0, 1], color=col, lw=1.1)
    axes[0].legend(handles=[Patch(color=C_KEEP, label="interior, measured (%d)" % n_int),
                            Patch(color=C_BRITTLE, label="border-truncated, "
                                                         "excluded (%d)" % n_bor)],
                   loc="upper left", fontsize=9, framealpha=0.85,
                   facecolor="white", frameon=True)
    axes[0].set_title("every grain the segmentation found  (%d)" % len(grains))
    _scalebar(axes[0], seg.pixel_size_um, sem.intensity.shape[1])

    # right: the ones that became solids, numbered, with size annotated
    solids = rec["solids"]
    by_id = {s.grain_id: s for s in solids}
    drawn = 0
    for g in grains:
        s = by_id.get(g.grain_id)
        if s is None:
            continue
        drawn += 1
        if drawn > max_labels:
            break
        mask = (seg.labels == g.label).astype(np.uint8)
        cs, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        for c in cs:
            axes[1].plot(c[:, 0, 0], c[:, 0, 1], color=C_DUCTILE, lw=1.2)
        axes[1].text(g.centroid_x_um / seg.pixel_size_um,
                     g.centroid_y_um / seg.pixel_size_um, str(g.grain_id),
                     color="white", fontsize=7, ha="center", va="center",
                     fontweight="bold",
                     path_effects=[pe.withStroke(linewidth=1.8,
                                                 foreground="black")])
    # The id is the grain's rank by area across EVERY region found, not a
    # renumbering of the survivors -- so the ids are sparse, and saying
    # "numbered by size" without saying over what invites the reader to read
    # a gap as a missing grain.
    axes[1].set_title("%d of them became verified 3-D solids   "
                      "(label = rank by area among all %d)"
                      % (len(solids), len(grains)))

    fig.suptitle("3 - What was measured, and what was deliberately not",
                 fontsize=13, fontweight="bold")
    return fig


# ---------------------------------------------------------------------------
# 3. the measurements themselves
# ---------------------------------------------------------------------------

def measurement_distributions(grains: Sequence, *, interior_only=True,
                              figsize=(14, 7.4)):
    """The six descriptors a grinding model cares about, as distributions.

    25 are measured per grain; these are the ones that change the answer --
    size sets the chip, aspect and circularity set how it cuts, solidity and
    corner angle say how sharp it is.
    """
    plt = _plt()
    g = [x for x in grains if not (interior_only and x.touches_border)]
    if not g:
        g = list(grains)
    fields = [
        ("equivalent_diameter_um", "equivalent diameter", "um", C_DUCTILE),
        ("feret_max_um", "max Feret (length)", "um", C_DUCTILE),
        ("aspect_ratio", "aspect ratio", "-", C_ACCENT),
        ("circularity", "circularity", "-", C_ACCENT),
        ("solidity", "solidity", "-", C_KEEP),
        ("min_corner_angle_deg", "sharpest corner", "deg", C_BRITTLE),
    ]
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    for ax, (key, label, unit, col) in zip(axes.ravel(), fields):
        v = np.array([getattr(x, key) for x in g], dtype=float)
        v = v[np.isfinite(v)]
        if not v.size:
            ax.axis("off")
            continue
        ax.hist(v, bins=min(24, max(6, len(v) // 3)), color=col, alpha=0.85,
                edgecolor="black", linewidth=0.5)
        d50 = float(np.median(v))
        ax.axvline(d50, color="black", ls="--", lw=1.3)
        ax.text(0.97, 0.93, "d50 %.3g %s\nd10 %.3g   d90 %.3g"
                % (d50, unit, np.percentile(v, 10), np.percentile(v, 90)),
                transform=ax.transAxes, ha="right", va="top", fontsize=8.6,
                family="monospace")
        ax.set_xlabel("%s  [%s]" % (label, unit))
        ax.set_ylabel("grains")
        ax.grid(alpha=0.2)
    fig.suptitle("4 - Measured grain population  (n = %d interior grains, "
                 "25 descriptors each)" % len(g),
                 fontsize=13, fontweight="bold")
    return fig


def outline_fidelity(rec: dict, n=6, figsize=(14, 4.6)):
    """Measured outline against its convex hull, worst offenders first.

    The original notebook simplified every grain to a convex hull of 8-10
    vertices, which erases exactly the concave features that do the cutting.
    The shaded gap between the two curves here is what that cost.
    """
    plt = _plt()
    from shapely.geometry import Polygon

    # Ranked by how much the hull would ADD, not by grain size. The claim is
    # that hulling erases concave cutting features, and a panel showing a grain
    # that happens to be convex demonstrates nothing.
    scored = []
    for sol in rec["solids"]:
        poly = Polygon(sol.outline_um)
        hull = poly.convex_hull
        if hull.area > 0 and poly.area > 0:
            scored.append((100.0 * (1.0 - poly.area / hull.area), sol, poly, hull))
    scored.sort(key=lambda t: -t[0])
    scored = scored[:n]
    if not scored:
        raise ValueError("no measurable outlines in this record")

    fig, axes = plt.subplots(1, len(scored), figsize=figsize, squeeze=False)
    axes = axes[0]
    for ax, (lost, sol, poly, hull) in zip(axes, scored):
        hx, hy = hull.exterior.xy
        ring = sol.outline_um
        ax.fill(hx, hy, color=C_BRITTLE, alpha=0.28, label="convex hull")
        ax.plot(hx, hy, color=C_BRITTLE, lw=1.1, ls="--")
        ax.fill(ring[:, 0], ring[:, 1], color=C_DUCTILE, alpha=0.60,
                label="measured outline")
        ax.plot(np.append(ring[:, 0], ring[0, 0]),
                np.append(ring[:, 1], ring[0, 1]), color=C_DUCTILE, lw=1.5)
        ax.set_title("#%d   %d vertices\nhull would add %.1f %%"
                     % (sol.grain_id, len(ring), max(lost, 0.0)), fontsize=9.5)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, fontsize=10, frameon=False)
    fig.suptitle("5 - Real outlines, not convex hulls: the concave features are "
                 "the cutting edges", fontsize=13, fontweight="bold")
    return fig


# ---------------------------------------------------------------------------
# 4. the 3-D solids
# ---------------------------------------------------------------------------

def grain_gallery(solids: Sequence, n=8, figsize=(15, 7.6)):
    """The lofted 3-D grains, drawn from the very triangles the deck writes."""
    plt = _plt()
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    sel = sorted(solids, key=lambda s: -s.mesh_volume_um3)[:n]
    cols = 4
    rows = int(math.ceil(len(sel) / cols)) or 1
    fig = plt.figure(figsize=figsize)
    for i, s in enumerate(sel):
        ax = fig.add_subplot(rows, cols, i + 1, projection="3d")
        v, f = s.vertices, s.faces
        tri = v[f]
        # Flat-shade by facet normal against a fixed light. Without it every
        # facet is one flat blue and a polyhedron reads as a silhouette, which
        # defeats the point of showing that the grain is faceted at all.
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        ln = np.linalg.norm(n, axis=1, keepdims=True)
        n = n / np.where(ln > 0, ln, 1.0)
        light = np.array([0.35, -0.55, 0.76])
        shade = 0.45 + 0.55 * np.clip(np.abs(n @ light), 0, 1)
        from matplotlib.colors import to_rgb
        base = np.array(to_rgb(C_DUCTILE))
        facecols = np.clip(shade[:, None] * base[None, :], 0, 1)
        pc = Poly3DCollection(tri, facecolors=facecols, edgecolor="black",
                              linewidths=0.15)
        ax.add_collection3d(pc)
        c = v.mean(axis=0)
        r = float(np.abs(v - c).max()) * 0.72
        ax.set_xlim(c[0] - r, c[0] + r)
        ax.set_ylim(c[1] - r, c[1] + r)
        ax.set_zlim(c[2] - r, c[2] + r)
        ax.set_box_aspect((1, 1, 1))
        ax.set_title("#%d  %.1f um tall\n%d faces, %d tets"
                     % (s.grain_id, s.height_um, len(f), len(s.tets)),
                     fontsize=9)
        ax.set_axis_off()          # hides ticks, panes AND the spine lines
        ax.view_init(elev=24, azim=-58)
    fig.suptitle("6 - The measured grains as watertight 3-D solids  "
                 "(%d of %d shown; these exact triangles go into the .inp)"
                 % (len(sel), len(solids)), fontsize=13, fontweight="bold")
    return fig


def solid_verification(rec: dict, figsize=(14, 4.6)):
    """Every grain's own verification numbers, as distributions.

    Each solid is checked against closed-form geometry when it is built: mesh
    volume against the analytic prismatoid sum, and maximum projected section
    against the measured outline. Both are asserted per grain; this is what the
    whole population looks like, which is the only way to see that the tolerance
    holds everywhere rather than on the one grain someone spot-checked.
    """
    plt = _plt()
    reports = [r for r in (rec.get("reports") or []) if r.get("ok") is not None]
    good = [r for r in reports if r.get("ok")]
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # A cached grain library carries the solids but not the per-grain validation
    # reports -- they are re-derived only on a fresh measurement. Say so, rather
    # than dividing by an empty set and taking the cell down with it.
    if not reports:
        for ax in axes:
            ax.axis("off")
        axes[1].text(0.5, 0.5,
                     "the per-grain verification reports\n"
                     "are not in the cache\n\n"
                     "they are produced when the grains are measured,\n"
                     "and this run reused a cached library.\n"
                     "Re-measure (cache=False, or delete\n"
                     "grain_library.pkl) to see them.",
                     ha="center", va="center", fontsize=10, color=C_GREY,
                     transform=axes[1].transAxes)
        fig.suptitle("7 - Every solid is verified against closed-form geometry "
                     "before it is allowed into the library",
                     fontsize=13, fontweight="bold")
        return fig

    ax = axes[0]
    v = np.array([abs(r.get("volume_rel_error", 0.0)) for r in good], dtype=float)
    v = v[np.isfinite(v)]
    v = np.maximum(v, 1e-18)
    ax.hist(np.log10(v), bins=22, color=C_KEEP, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("log10 |mesh volume - analytic volume| / analytic")
    ax.set_ylabel("grains")
    ax.set_title("volume closes to %.0e (worst)" % v.max())
    ax.grid(alpha=0.2)

    ax = axes[1]
    a = np.array([abs(r.get("projected_area_rel_error", 0.0)) for r in good],
                 dtype=float)
    a = np.maximum(a[np.isfinite(a)], 1e-18)
    ax.hist(np.log10(a), bins=22, color=C_DUCTILE, edgecolor="black",
            linewidth=0.5)
    ax.set_xlabel("log10 |projected section - measured outline| / outline")
    ax.set_title("silhouette reproduced to %.0e (worst)" % a.max())
    ax.grid(alpha=0.2)

    ax = axes[2]
    ax.axis("off")
    n_bad = len(reports) - len(good)
    issues = {}
    for r in reports:
        if not r.get("ok"):
            for it in r.get("issues", []):
                issues[it[:44]] = issues.get(it[:44], 0) + 1
    ax.text(0.0, 1.0, "PER-GRAIN GATES", fontsize=12, fontweight="bold",
            va="top", transform=ax.transAxes)
    rows = [("solids attempted", len(reports)),
            ("passed every check", len(good)),
            ("rejected", n_bad)]
    y = 0.86
    for k, val in rows:
        ax.text(0.0, y, k, fontsize=10, color=C_GREY, va="top",
                transform=ax.transAxes)
        ax.text(0.72, y, str(val), fontsize=10, family="monospace", va="top",
                color=C_WARN if (k == "rejected" and val) else "black",
                transform=ax.transAxes)
        y -= 0.10
    if issues:
        y -= 0.04
        ax.text(0.0, y, "why they were rejected:", fontsize=9.5, color=C_GREY,
                va="top", transform=ax.transAxes)
        y -= 0.09
        for msg, cnt in sorted(issues.items(), key=lambda kv: -kv[1])[:5]:
            ax.text(0.0, y, "%3d  %s" % (cnt, msg), fontsize=8.4,
                    family="monospace", va="top", transform=ax.transAxes)
            y -= 0.075
    else:
        y -= 0.06
        ax.text(0.0, y, "no grain failed a geometric check", fontsize=9.5,
                color=C_KEEP, va="top", transform=ax.transAxes)

    fig.suptitle("7 - Every solid is verified against closed-form geometry "
                 "before it is allowed into the library",
                 fontsize=13, fontweight="bold")
    return fig


# ---------------------------------------------------------------------------
# 4b. the assembled model
# ---------------------------------------------------------------------------

def assembly(plan: dict, figsize=(14, 5.4)):
    """Wheel, grit and workpiece, at the three scales they live at.

    ``semgrit.preview.preview`` is the general-purpose version of this and stays
    the right tool for a dressed multi-grit band. On a ONE-grit deck two of its
    four panels degenerate -- an engagement curve through a single point, and an
    unwrapped band holding one dot -- and most of the canvas comes out empty.
    This draws the three things that are actually informative at that scale:
    the whole wheel, the grit against the block, and the grit's own profile
    against the depth of cut.

    Everything is read from the plan, which is built by the same placement code
    the writer uses, so this cannot disagree with the deck.
    """
    plt = _plt()
    from matplotlib.patches import Rectangle

    R = plan["outer_radius_mm"]
    frames = plan.get("frames") or []
    wp = plan.get("workpiece")
    r_ground = plan["ground_radius_mm"]
    ae = plan.get("depth_of_cut_um") or 0.0

    fig, axes = plt.subplots(1, 3, figsize=figsize,
                             gridspec_kw={"width_ratios": [1.0, 1.3, 1.0]})

    # --- 1. the whole wheel, to scale ------------------------------------
    ax = axes[0]
    th = np.linspace(0, 2 * math.pi, 400)
    ax.plot(R * np.cos(th), R * np.sin(th), color=C_GREY, lw=1.0, ls="--")
    ax.plot((R - plan["rim_depth_mm"]) * np.cos(th),
            (R - plan["rim_depth_mm"]) * np.sin(th), color=C_GREY, lw=0.7,
            ls=":")
    t_c = math.radians(plan["theta_workpiece_deg"])
    sec = math.radians(plan["sector_deg"])
    ts = np.linspace(t_c - sec / 2, t_c + sec / 2, 60)
    ax.plot(R * np.cos(ts), R * np.sin(ts), color=C_DUCTILE, lw=3.5,
            solid_capstyle="butt")
    ax.plot([0], [0], marker="+", color="black", ms=9)
    ax.annotate("the %.1f mm arc\nthis deck models"
                % plan.get("arc_length_mm", 0.0),
                xy=(R * math.cos(t_c), R * math.sin(t_c)),
                xytext=(-R * 0.75, R * 0.62), fontsize=9, color=C_DUCTILE,
                arrowprops=dict(arrowstyle="->", color=C_DUCTILE, lw=1.2))
    ax.set_aspect("equal")
    ax.set_xlabel("x  [mm]")
    ax.set_ylabel("y  [mm]")
    ax.set_title("wheel, %.0f mm diameter" % (2 * R))
    ax.grid(alpha=0.2)

    # --- 2. the grit against the block, in microns -----------------------
    ax = axes[1]
    if wp is not None and frames:
        hb = wp["length_mm"] / 2.0
        # The grit's LOWER ENVELOPE -- the highest radial reach at each station,
        # which is the surface that can touch. Scattering every vertex instead
        # draws the far side of the grain too, and since this axis is radial
        # height rather than a section, that reads as the grit passing through
        # the block when nothing of the sort is happening.
        top = None
        for k, f in enumerate(frames):
            u = f[:, 1] * 1000.0
            hgt = (f[:, 0] - r_ground) * 1000.0
            # A 3-D grain has many vertices at the same station at different
            # heights, so joining them in station order draws a zigzag through
            # the body of the grain. Bin by station and keep the MAXIMUM: that
            # is the silhouette that can reach the work, which is the only part
            # of the grain this panel is about.
            nb_ = max(int(np.ceil((u.max() - u.min()) / 0.25)), 4)
            edges = np.linspace(u.min(), u.max(), nb_ + 1)
            idx = np.clip(np.digitize(u, edges) - 1, 0, nb_ - 1)
            env_u, env_h = [], []
            for b in range(nb_):
                sel = idx == b
                if sel.any():
                    env_u.append(0.5 * (edges[b] + edges[b + 1]))
                    env_h.append(hgt[sel].max())
            ax.fill_between(env_u, min(hgt.min(), -wp["depth_mm"] * 1000),
                            env_h, color=C_BRITTLE, alpha=0.30, zorder=3)
            ax.plot(env_u, env_h, color=C_BRITTLE, lw=1.8, zorder=4,
                    label="abrasive (rigid)" if k == 0 else None)
            top = hgt.max() if top is None else max(top, hgt.max())

        # Block drawn only as deep as the picture needs: the interesting scale
        # is the depth of cut, which is 0.2 um against a 6 um block, so drawing
        # the whole block makes the cut a hairline.
        lo = min(-wp["depth_mm"] * 1000,
                 min(float((f[:, 0] - r_ground).min() * 1000) for f in frames))
        ax.add_patch(Rectangle((-hb * 1000, lo * 1.05), 2 * hb * 1000,
                               abs(lo * 1.05), facecolor="#E4E4E4",
                               edgecolor="none", zorder=0))
        ax.axhline(0, color=C_DUCTILE, lw=2.2, zorder=4)
        ax.annotate("ground face of the workpiece", xy=(-hb * 1000 * 0.97, 0),
                    xytext=(0, 5), textcoords="offset points", fontsize=9,
                    color=C_DUCTILE, fontweight="bold")
        if ae:
            ax.axhline(-ae, color=C_WARN, lw=1.5, ls="--", zorder=4)
            ax.annotate("after %.2f um infeed" % ae,
                        xy=(hb * 1000 * 0.97, -ae), xytext=(0, -4),
                        textcoords="offset points", fontsize=9, color=C_WARN,
                        ha="right", va="top")
        ax.set_xlim(-hb * 1000 * 1.02, hb * 1000 * 1.02)
        # The interesting scale is the depth of cut against the grit tip, which
        # is a few hundred nanometres. Showing the grain's full 8 um of relief
        # compresses that to a hairline, so the view is clipped to the contact
        # zone and the caption says the grain continues below it.
        deep = max(4.0 * (ae or 0.3), 1.2)
        ax.set_ylim(-deep, max(0.35 * deep, 0.25))
        ax.annotate("zero overclosure at t=0 -- the grain continues %.1f um "
                    "below this view"
                    % (abs(min(float((f[:, 0] - r_ground).min() * 1000)
                               for f in frames)) - deep),
                    xy=(0.5, 0.02), xycoords="axes fraction", ha="center",
                    fontsize=8.2, color=C_GREY, style="italic")
        ax.set_xlabel("along the scratch  [um]")
        ax.set_ylabel("height above the ground face  [um]")
        gap_nm = abs(top or 0.0) * 1000.0
        ax.set_title("the abrasive where it meets the work  "
                     "(%d grit%s, seated %s)"
                     % (len(frames), "" if len(frames) == 1 else "s",
                        "exactly tangent" if gap_nm < 0.01
                        else "%.2f nm clear" % gap_nm))
        ax.legend(loc="upper left", fontsize=8.5, frameon=True,
                  facecolor="white", framealpha=0.9)
        ax.grid(alpha=0.25)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "no workpiece in this deck", ha="center",
                va="center", color=C_GREY, transform=ax.transAxes)

    # --- 3. the numbers ---------------------------------------------------
    panel = _TextPanel(axes[2], "THE MODEL", line=0.077, val_x=0.60)
    panel.row("wheel", "%.0f mm dia" % (2 * R))
    panel.row("arc modelled", "%.2f mm" % plan.get("arc_length_mm", 0.0))
    panel.row("grits", "%d" % plan.get("n_grits", 0))
    panel.row("engaging the block", "%d" % plan.get("n_grits_under_block", 0))
    ph = plan.get("protrusion_um") or {}
    if ph.get("n"):
        panel.row("protrusion", "%.2f um" % ph.get("max", 0.0))
    panel.gap()
    if wp is not None:
        panel.row("workpiece", "%.0f x %.0f x %.0f um"
                  % (wp["length_mm"] * 1000, wp["width_mm"] * 1000,
                     wp["depth_mm"] * 1000))
        panel.row("elements", format(plan.get("n_workpiece_elements", 0), ","))
        el = plan.get("element_um") or (0, 0, 0, 0)
        panel.row("surface element", "%.3f x %.3f x %.4f um"
                  % (el[0], el[1], el[2]))
        lo = min(el[0], el[1], el[2]) or 1e-12
        panel.row("aspect ratio", "%.1f : 1" % (max(el[0], el[1], el[2]) / lo))
    panel.gap()
    panel.row("depth of cut", "%.3f um" % ae)
    cost = plan.get("cost") or {}
    if cost:
        panel.row("stable increment", "%.2e s" % (cost.get("stable_dt_s") or 0))
        panel.row("increments", format(int(cost.get("increments") or 0), ","))
        panel.row("wall clock, 8 cores",
                  "%.1f h" % ((cost.get("est_hours") or {}).get("8", 0.0)))

    fig.suptitle("5 - The assembled model: %s"
                 % plan.get("title", "grinding deck"),
                 fontsize=13, fontweight="bold")
    return fig


# ---------------------------------------------------------------------------
# 5. the ductile / brittle physics
# ---------------------------------------------------------------------------

def chip_thickness(field, dc_mm: float, wp_length_mm: float,
                   *, title=None, figsize=(13, 5.2)):
    """h(u) along the scratch against dc -- the figure the model exists for.

    ``field`` is a :class:`semgrit.hybrid.ChipField`. The crossing drawn here is
    the same one written into the material card, so where this figure says the
    transition is, is where the subroutine will put it.
    """
    plt = _plt()
    hb = wp_length_mm / 2.0
    u = np.linspace(-hb, hb, 900)
    h = np.array([field.h_at(x) for x in u])
    u_um, h_nm, dc_nm = u * 1000.0, h * 1e6, dc_mm * 1e6

    fig, axes = plt.subplots(1, 2, figsize=figsize,
                             gridspec_kw={"width_ratios": [1.7, 1.0]})
    ax = axes[0]
    duct = h_nm < dc_nm
    ax.fill_between(u_um, 0, h_nm, where=duct, color=C_DUCTILE, alpha=0.30,
                    interpolate=True, label="ductile  h < dc")
    ax.fill_between(u_um, 0, h_nm, where=~duct, color=C_BRITTLE, alpha=0.30,
                    interpolate=True, label="brittle  h >= dc")
    ax.plot(u_um, h_nm, color="black", lw=1.8)
    ax.axhline(dc_nm, color=C_WARN, ls="--", lw=1.6)
    ax.annotate("dc = %.1f nm" % dc_nm, xy=(0.995, dc_nm),
                xycoords=("axes fraction", "data"), xytext=(0, 4),
                textcoords="offset points", color=C_WARN, fontsize=10,
                ha="right", va="bottom", fontweight="bold")

    # Only mark a transition the drawn curve actually has. transition_u_mm is
    # supplied by the caller, and a stale or mismatched one would put a
    # confident dotted line through a pass that never crosses dc at all -- the
    # figure would then contradict its own curve, which is worse than silence.
    tu = field.transition_u_mm
    crosses = bool(duct.any() and (~duct).any())
    if tu is not None and crosses and -hb <= tu <= hb:
        ax.axvline(tu * 1000.0, color="black", ls=":", lw=1.4)
        ax.annotate("transition\nu = %.2f um" % (tu * 1000.0),
                    xy=(tu * 1000.0, dc_nm), xytext=(14, -34),
                    textcoords="offset points", fontsize=9.5,
                    ha="left", va="top",
                    arrowprops=dict(arrowstyle="->", lw=1.0))
    elif not crosses:
        ax.annotate("h never crosses dc over this block:\n"
                    "the whole pass is %s"
                    % ("ductile" if duct.all() else "brittle"),
                    xy=(0.5, 0.06), xycoords="axes fraction", ha="center",
                    fontsize=9.5, color=C_GREY, style="italic")

    ax.set_xlabel("station along the scratch  u  [um]")
    ax.set_ylabel("undeformed chip thickness  h  [nm]")
    # Headroom so the dc line and its label clear the legend and the frame.
    ax.set_ylim(0, max(float(h_nm.max()), dc_nm) * 1.28)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9, frameon=True,
              facecolor="white")
    ax.grid(alpha=0.25)
    ax.set_title("h(u) = H0 + HG u - u$^2$/(2 R$_{tip}$)", pad=10)

    ax = axes[1]
    frac = 100.0 * float(duct.mean())
    panel = _TextPanel(ax, "THE SWITCH", line=0.078, val_x=0.52)
    panel.row("H0", "%.4f um" % (field.h0_mm * 1000))
    panel.row("HG", "%.3e" % field.hg)
    panel.row("R_tip", "%.4f mm" % field.rtip_mm)
    panel.row("h at entry", "%.1f nm" % (field.h_entry_mm * 1e6))
    panel.row("h at exit", "%.1f nm" % (field.h_exit_mm * 1e6))
    panel.row("dc", "%.1f nm" % dc_nm)
    panel.gap()
    panel.row("ductile", "%.0f %% of the pass" % frac, colour=C_DUCTILE)
    panel.row("brittle", "%.0f %% of the pass" % (100 - frac), colour=C_BRITTLE)
    panel.gap()
    # The sagitta is the reason the quadratic term is there at all, so say how
    # big it is rather than leaving it as an unexplained term in the caption.
    sag = (hb * hb) / (2.0 * field.rtip_mm) * 1e6 if field.rtip_mm > 0 else 0.0
    panel.note("path curvature over the block: %.1f nm -- the same order as dc, "
               "which is why the quadratic term is kept" % sag, C_GREY,
               size=8.6, wrap=42)

    fig.suptitle(title or "8 - Where removal stops being ductile",
                 fontsize=13, fontweight="bold")
    return fig


def dc_forms(materials_map, figsize=(12, 5.0)):
    """The three critical-depth formulas side by side, per material.

    They differ by (E/H)^1.5 -- about 17x on sandstone -- so lambda_c belongs to
    one form and is not transferable. Showing them together is the honest way to
    present a quantity the literature does not agree on.
    """
    plt = _plt()
    names = list(materials_map)
    fig, ax = plt.subplots(figsize=figsize)
    w = 0.36
    x = np.arange(len(names))
    for i, form in enumerate((1, 2)):
        vals = [materials_map[n].dc_nm(form) for n in names]
        b = ax.bar(x + (i - 0.5) * w, vals, w,
                   color=(C_DUCTILE if form == 1 else C_BRITTLE),
                   edgecolor="black", linewidth=0.6,
                   label=("form 1:  $\\lambda_c (H/E)^{1/2}(K_c/H)^2$" if form == 1
                          else "form 2 (Bifano):  $\\lambda_c (E/H)(K_c/H)^2$"))
        for r, v in zip(b, vals):
            ax.text(r.get_x() + r.get_width() / 2, v, " %.1f" % v, ha="center",
                    va="bottom", fontsize=9, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([materials_map[n].label for n in names])
    ax.set_ylabel("critical depth of cut  dc  [nm]")
    ax.set_yscale("log")
    ax.legend(fontsize=9.5)
    ax.grid(alpha=0.25, axis="y")
    ax.set_title("The two published forms differ by $(E/H)^{3/2}$ -- "
                 "$\\lambda_c$ is not transferable between them")
    fig.suptitle("9 - Critical depth of cut: which formula, and why it matters",
                 fontsize=13, fontweight="bold")
    return fig


__all__ = [
    "calibration", "segmentation_stages", "segmentation_overlay",
    "measurement_distributions", "outline_fidelity", "grain_gallery",
    "solid_verification", "assembly", "chip_thickness", "dc_forms",
]


def demo(image=None, outdir="_figures_demo"):
    """Render every figure from a real image and assert each one drew something.

    ``python -m semgrit.figures [image.tif]``

    Every figure here is called once per notebook run and nowhere else, so a
    typo in one of them surfaces in front of an audience rather than in a test.
    The check is deliberately blunt -- that each figure builds, has the axes it
    claims, and contains drawn artists -- because that is what actually breaks:
    a shadowed variable, a renamed plan key, a panel that silently comes out
    empty.
    """
    import glob
    import os

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .quick import measure_images

    if image is None:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        found = sorted(glob.glob(os.path.join(here, "*.tif")))
        if not found:
            raise SystemExit("no .tif found next to the repo root; pass one")
        image = found[0]

    os.makedirs(outdir, exist_ok=True)
    rec = measure_images([image], os.path.join(outdir, "meas"),
                         keep_stages=True, cache=False,
                         log=lambda *a: None)["per_image"][0]

    # (name, thunk, minimum number of axes that must carry artists)
    cases = [
        ("calibration", lambda: calibration(rec), 3),
        ("segmentation_stages", lambda: segmentation_stages(rec), 12),
        ("segmentation_overlay", lambda: segmentation_overlay(rec), 2),
        ("measurement_distributions",
         lambda: measurement_distributions(rec["grains"]), 6),
        ("outline_fidelity", lambda: outline_fidelity(rec), 3),
        ("grain_gallery", lambda: grain_gallery(rec["solids"]), 4),
        ("solid_verification", lambda: solid_verification(rec), 3),
    ]
    failed = 0
    for name, thunk, min_axes in cases:
        try:
            fig = thunk()
            drawn = sum(1 for ax in fig.axes
                        if ax.images or ax.lines or ax.patches or
                        ax.collections or ax.texts)
            assert drawn >= min_axes, ("%s drew %d populated axes, expected >= %d"
                                       % (name, drawn, min_axes))
            path = os.path.join(outdir, name + ".png")
            fig.savefig(path)
            plt.close(fig)
            print("  ok    %-26s %6.0f kB" % (name, os.path.getsize(path) / 1e3))
        except Exception as exc:                          # noqa: BLE001
            failed += 1
            print("  FAIL  %-26s %s" % (name, exc))

    # The physics figures need a deck, so they are checked against a synthetic
    # ChipField rather than a full build -- the drawing is what is under test.
    try:
        from .hybrid import ChipField, _transition_station
        dc = 8.775e-5                                  # sandstone, form 2, mm
        # A wedge that genuinely straddles dc across the block, so the figure's
        # transition branch is the one under test. The crossing is solved for
        # rather than asserted, exactly as chip_field does it.
        fld = ChipField(theta_c=0.0, h0_mm=dc * 1.25, hg=-2.6e-3, rtip_mm=25.0,
                        u_gov_mm=0.0, h_entry_mm=0.0, h_exit_mm=0.0)
        fld.h_entry_mm, fld.h_exit_mm = fld.h_at(0.024), fld.h_at(-0.024)
        fld.transition_u_mm = _transition_station(fld, -0.024, 0.024, dc)
        assert fld.transition_u_mm is not None, "demo field never crosses dc"
        fig = chip_thickness(fld, dc, 0.048)
        assert len(fig.axes) == 2
        # the crossing the figure marks must be the crossing the field has
        assert abs(fld.h_at(fld.transition_u_mm) - dc) < 1e-9 * dc
        fig.savefig(os.path.join(outdir, "chip_thickness.png"))
        plt.close(fig)
        print("  ok    %-26s transition at u = %+.2f um"
              % ("chip_thickness", fld.transition_u_mm * 1000))
    except Exception as exc:                              # noqa: BLE001
        failed += 1
        print("  FAIL  %-26s %s" % ("chip_thickness", exc))

    # assembly() needs a real plan, so it is checked against a small deck built
    # from the grains just measured -- cheap, because plan_deck writes nothing.
    try:
        from .analysis import AnalysisParams
        from .build_deck import DeckParams, plan_deck
        from .hybrid import HYBRID_DEPVAR
        from . import materials as _mat
        _hp = _mat.hybrid_params("sandstone", h_source=0, dc_form=2)
        _p = DeckParams(
            name="figures_demo", diameter_mm=50.0, sector_mode="arc",
            arc_length_mm=2.0, rim_depth_mm=0.012, width_mm=0.030,
            include_bond=False, grit_mode="single", single_grain_index=-1,
            single_grit_offset_mm=0.015, include_workpiece=True,
            wp_length_mm=0.048, wp_width_mm=0.020, wp_depth_mm=0.006,
            wp_element_size_length_mm=0.0003,
            wp_element_size_width_mm=0.0003,
            wp_element_size_depth_mm=0.0003,
            clearance_um=0.0, wp_position="centred",
            surface_speed_mm_s=30_000.0, cores=8,
            analysis=AnalysisParams(enabled=True, depth_of_cut_um=0.20,
                                    material_model="hybrid", hybrid=_hp,
                                    n_depvar=HYBRID_DEPVAR,
                                    element_deletion=True))
        _mat.apply(_p, "sandstone")
        fig = assembly(plan_deck(_p, rec["solids"]))
        assert len(fig.axes) == 3
        fig.savefig(os.path.join(outdir, "assembly.png"))
        plt.close(fig)
        print("  ok    %-26s" % "assembly")
    except Exception as exc:                              # noqa: BLE001
        failed += 1
        print("  FAIL  %-26s %s" % ("assembly", exc))

    try:
        from .materials import MATERIALS
        fig = dc_forms(MATERIALS)
        fig.savefig(os.path.join(outdir, "dc_forms.png"))
        plt.close(fig)
        print("  ok    %-26s" % "dc_forms")
    except Exception as exc:                              # noqa: BLE001
        failed += 1
        print("  FAIL  %-26s %s" % ("dc_forms", exc))

    print("%s  (%s)" % ("ALL FIGURES OK" if not failed else "%d FAILED" % failed,
                        outdir))
    return failed


if __name__ == "__main__":
    import sys
    raise SystemExit(demo(sys.argv[1] if len(sys.argv) > 1 else None))
