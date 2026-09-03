"""Figures for shape adaptive grinding.

``semgrit.figures`` covers the SEM measurement pipeline. These four cover the
SAG-specific physics, and each one exists to make a specific argument visible
rather than to decorate the notebook:

:func:`load_collapse`
    Why SAG works at all. The per-grain load against wheel compression, for
    all three pads, on log axes. A conventional wheel loads one grit with the
    whole force; a compliant one spreads it over a patch and the load per grain
    falls by four orders of magnitude. That collapse IS the process.

:func:`contact_patch`
    The patch, to scale on the wheel face, with the Hertzian pressure over it.
    Also shows the part the paper's own fits imply but cannot have: the
    ellipse is wider than the 10 mm face at every operating point, so the patch
    is clipped.

:func:`dc_comparison`
    The three critical-depth expressions on the same material. They differ by
    orders of magnitude, which is why a deck that does not say which one it
    used cannot be checked.

:func:`regime_map`
    Where this operating point sits: indentation against dc, with the ductile
    and brittle regions marked, and the paper's measured chip range for
    comparison.

Palette is the project's Okabe-Ito set, so a SAG figure and a rigid-wheel
figure of the same quantity agree on colour.
"""

from __future__ import annotations

import math
from typing import Optional

C_DUCTILE = "#0072B2"      # blue
C_BRITTLE = "#D55E00"      # vermillion
C_PAD = {6.0: "#0072B2", 15.0: "#009E73", 30.0: "#D55E00"}
C_GREY = "#7a828a"
C_FACE = "#E69F00"         # orange, for the wheel face limit


def _plt():
    import matplotlib.pyplot as plt
    return plt


def _pads(plan: dict, grains=(6.0, 15.0, 30.0), *, compressions=None):
    """Re-solve the contact for other pads/compressions off the same tool."""
    from . import sag as _sag

    from . import materials as _mat

    p = plan["params"]
    tool = p.tool()
    # The workpiece elastic constants live on the material's hybrid params,
    # not on SAGParams -- sagdeck.plan reads them from there and so must this.
    hp = _mat.get(p.material).hybrid_params()
    out = {}
    for dg in grains:
        t = _sag.Tool(diameter_mm=tool.diameter_mm, width_mm=tool.width_mm,
                      pad=_sag.Pad(dg), shore_a=tool.shore_a,
                      elastic_mpa=tool.elastic_mpa, poisson=tool.poisson,
                      layer_thickness_mm=tool.layer_thickness_mm)
        rows = []
        for T in (compressions if compressions is not None
                  else [0.05 * i for i in range(1, 17)]):
            # NOT wrapped in a bare except. Swallowing here is what hid a
            # wrong attribute name and produced four figures with no curves
            # on them -- the legend warning was the only symptom. A contact
            # that cannot be solved at a valid compression is a bug, so it
            # should surface.
            rows.append((T, _sag.solve_contact(
                t, compression_mm=T, speed_rpm=p.speed_rpm,
                work_modulus_mpa=hp.youngs_mpa, work_poisson=hp.poisson,
                bhn_kgf_mm2=p.bhn_kgf_mm2)))
        if not rows:
            raise RuntimeError("no contact solved for the %g um pad" % dg)
        out[dg] = rows
    return out


def load_collapse(plan: dict, figsize=(13.5, 5.0)):
    """Per-grain load and indentation against compression, all three pads."""
    plt = _plt()
    p = plan["params"]
    c = plan["contact"]
    data = _pads(plan)
    dc_nm = plan["material"]["dc_nm"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    for dg, rows in sorted(data.items()):
        if not rows:
            continue
        T = [r[0] for r in rows]
        ax1.plot(T, [r[1].load_per_grain_n for r in rows], "-",
                 color=C_PAD.get(dg, C_GREY), lw=2.0,
                 label="%g $\\mu$m pad" % dg)
        ax2.plot(T, [r[1].indentation_nm for r in rows], "-",
                 color=C_PAD.get(dg, C_GREY), lw=2.0,
                 label="%g $\\mu$m pad" % dg)

    # The paper's stated band for the per-grain normal force.
    ax1.axhspan(1e-5, 1e-4, color=C_GREY, alpha=0.15, zorder=0)
    ax1.annotate("measured range,\n30 $\\mu$m pad", xy=(0.06, 3e-5),
                 fontsize=8.5, color="#3a4148", va="center")
    ax1.set_yscale("log")
    ax1.set_xlabel("wheel compression $T$  (mm)")
    ax1.set_ylabel("normal load per grain $F_n$  (N)")
    ax1.set_title("The load each grain carries", fontsize=11)
    ax1.legend(fontsize=9, frameon=False)
    ax1.grid(alpha=0.25)

    ax2.axhline(dc_nm, color=C_BRITTLE, ls="--", lw=1.5)
    ax2.annotate("$d_c$ = %.0f nm" % dc_nm,
                 xy=(ax2.get_xlim()[1] * 0.98, dc_nm), ha="right",
                 va="bottom", fontsize=9, color=C_BRITTLE)
    ax2.set_yscale("log")
    ax2.set_xlabel("wheel compression $T$  (mm)")
    ax2.set_ylabel("indentation depth $d$  (nm)")
    ax2.set_title("...and how deep it cuts", fontsize=11)
    ax2.legend(fontsize=9, frameon=False)
    ax2.grid(alpha=0.25)

    for ax in (ax1, ax2):
        ax.axvline(p.compression_mm, color="#2b2b2b", lw=1.0, alpha=0.6)
    ax1.plot([p.compression_mm], [c.load_per_grain_n], "o", ms=7,
             color="#2b2b2b", zorder=5)
    ax2.plot([p.compression_mm], [c.indentation_nm], "o", ms=7,
             color="#2b2b2b", zorder=5)

    fig.suptitle("Why SAG can cut a brittle material ductilely: %s grains "
                 "share the load, so each takes %.3f nm"
                 % (format(int(c.active_grains), ","), c.indentation_nm),
                 fontsize=11.5, y=1.0)
    fig.tight_layout()
    return fig


def contact_patch(plan: dict, figsize=(13.5, 5.0)):
    """The patch on the wheel face, and the pressure across it."""
    plt = _plt()
    import numpy as np

    c = plan["contact"]
    p = plan["params"]
    a, b = c.semi_axis_a_mm, c.semi_axis_b_mm
    W = p.width_mm

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # --- the patch, to scale, against the face ------------------------
    th = np.linspace(0, 2 * math.pi, 400)
    b_free = c.spot_area_mm2 / (math.pi * a)
    ax1.fill(a * np.cos(th), b_free * np.sin(th), color=C_DUCTILE,
             alpha=0.18, lw=0)
    ax1.plot(a * np.cos(th), b_free * np.sin(th), color=C_DUCTILE, lw=1.6,
             label="ellipse from eqs. 6 and 7")
    ax1.axhline(+0.5 * W, color=C_FACE, lw=2.2)
    ax1.axhline(-0.5 * W, color=C_FACE, lw=2.2,
                label="wheel face, %.0f mm wide" % W)
    ax1.set_aspect("equal")
    ax1.set_xlabel("along the arc  (mm)")
    ax1.set_ylabel("across the face  (mm)")
    over = 2.0 * b_free / W
    ax1.set_title("The patch is %s the face"
                  % ("%.1f%% WIDER than" % (100.0 * (over - 1.0))
                     if over > 1.0 else "inside"), fontsize=11)
    ax1.legend(fontsize=8.5, frameon=False, loc="upper right")
    ax1.grid(alpha=0.2)
    if over > 1.0:
        ax1.annotate("clipped: %.1f%% of the\nnominal area is off the wheel"
                     % (100.0 * c.area_clipped_fraction),
                     xy=(0, 0.5 * (0.5 * W + b_free)), ha="center",
                     fontsize=8.5, color="#8a4b00")

    # --- the pressure -------------------------------------------------
    x = np.linspace(-a, a, 300)
    ax2.plot(x, c.max_pressure_mpa * np.sqrt(np.clip(1 - (x / a) ** 2, 0, 1)),
             color=C_DUCTILE, lw=2.0)
    ax2.fill_between(x, 0,
                     c.max_pressure_mpa
                     * np.sqrt(np.clip(1 - (x / a) ** 2, 0, 1)),
                     color=C_DUCTILE, alpha=0.15)
    ax2.axhline(c.mean_pressure_mpa, color=C_GREY, ls=":", lw=1.4)
    ax2.annotate("mean %.4f MPa" % c.mean_pressure_mpa,
                 xy=(-a * 0.95, c.mean_pressure_mpa), fontsize=8.5,
                 va="bottom", color="#3a4148")
    ax2.set_xlabel("along the arc  (mm)")
    ax2.set_ylabel("contact pressure  (MPa)")
    ax2.set_title("Hertzian pressure, peak %.4f MPa" % c.max_pressure_mpa,
                  fontsize=11)
    ax2.grid(alpha=0.25)

    fig.suptitle("$F_N$ = %.3f N over %.1f mm$^2$ -- a pressure three orders "
                 "below a conventional wheel's"
                 % (c.normal_load_n, c.spot_area_mm2), fontsize=11.5, y=1.0)
    fig.tight_layout()
    return fig


def dc_comparison(plan: dict, figsize=(12.0, 5.0)):
    """The three critical-depth expressions, on this material."""
    plt = _plt()

    forms = plan["dc_forms"]
    mat = plan["material"]
    # The energy form is not in dc_forms -- it is what the subroutine derives,
    # PSI*Kc^2/(E*H) with PSI pinned so it trips at the dc the card carries.
    # Computed here from the same numbers so all three sit on one axis.
    H, E = forms["hardness_mpa"], forms["youngs_mpa"]
    Kc = forms["kic_mpa_sqrt_mm"]
    energy_nm = forms["lambda_c"] * (H / E) * (Kc / H) ** 2 * 1e6
    names, vals, cols = [], [], []
    for val, label in (
            (forms["form1_nm"],
             "form 1\n$\\lambda_c(H/E)^{1/2}(K_c/H)^2$"),
            (forms["form2_nm"],
             "form 2, Bifano\n$\\lambda_c(E/H)(K_c/H)^2$"),
            (energy_nm,
             "energy\n$\\lambda_c(H/E)(K_c/H)^2$")):
        if val and val > 0:
            names.append(label)
            vals.append(val)
            cols.append(C_GREY)
    meas = None
    if forms.get("measured_lo_nm") and forms.get("measured_hi_nm"):
        meas = (forms["measured_lo_nm"], forms["measured_hi_nm"])

    fig, ax = plt.subplots(figsize=figsize)
    xs = range(len(vals))
    ax.bar(xs, vals, color=cols, width=0.55, zorder=3)
    for x, v in zip(xs, vals):
        ax.annotate("%.1f nm" % v, xy=(x, v), xytext=(0, 4),
                    textcoords="offset points", ha="center", fontsize=9.5)
    ax.set_xticks(list(xs))
    ax.set_xticklabels(names, fontsize=9)
    ax.set_yscale("log")
    ax.set_ylabel("critical depth of cut  (nm)")
    ax.grid(axis="y", alpha=0.25, zorder=0)

    used = mat["dc_nm"]
    ax.axhline(used, color=C_DUCTILE, lw=2.0, zorder=4)
    ax.annotate("this deck uses %.1f nm%s" % (used,
                                              "  (MEASURED)"
                                              if mat.get("dc_measured")
                                              else ""),
                xy=(len(vals) - 0.4, used), ha="right", va="bottom",
                fontsize=9.5, color=C_DUCTILE, zorder=5)
    if meas:
        ax.axhspan(meas[0], meas[1], color=C_DUCTILE, alpha=0.13, zorder=1)
        ax.annotate("measured, %g-%g nm" % (meas[0], meas[1]),
                    xy=(-0.42, meas[1]), fontsize=9, va="bottom",
                    color=C_DUCTILE)

    worst = max(vals) / used if used else 0.0
    ax.set_title("The three expressions disagree by %.0f$\\times$ on this "
                 "material -- which is why the deck records the one it used"
                 % worst, fontsize=11)
    fig.tight_layout()
    return fig


def regime_map(plan: dict, figsize=(12.0, 5.2)):
    """Indentation against dc across the compression range, per pad."""
    plt = _plt()

    dc = plan["material"]["dc_nm"]
    c = plan["contact"]
    p = plan["params"]
    data = _pads(plan)

    fig, ax = plt.subplots(figsize=figsize)
    for dg, rows in sorted(data.items()):
        if not rows:
            continue
        ax.plot([r[0] for r in rows],
                [r[1].indentation_nm / dc for r in rows], "-",
                color=C_PAD.get(dg, C_GREY), lw=2.0,
                label="%g $\\mu$m pad" % dg)

    ax.axhline(1.0, color=C_BRITTLE, lw=2.0)
    lo, hi = ax.get_xlim()
    ax.fill_between([lo, hi], 1.0, 1e4, color=C_BRITTLE, alpha=0.10, zorder=0)
    ax.fill_between([lo, hi], 1e-6, 1.0, color=C_DUCTILE, alpha=0.10, zorder=0)
    ax.set_xlim(lo, hi)
    ax.annotate("BRITTLE   $d > d_c$", xy=(hi * 0.98, 2.0), ha="right",
                fontsize=10, color=C_BRITTLE, weight="bold")
    ax.annotate("DUCTILE   $d < d_c$", xy=(hi * 0.98, 0.3), ha="right",
                fontsize=10, color=C_DUCTILE, weight="bold")

    ax.plot([p.compression_mm], [c.indentation_nm / dc], "o", ms=8,
            color="#2b2b2b", zorder=6)
    ax.annotate("this deck:\n$d/d_c$ = %.4f" % (c.indentation_nm / dc),
                xy=(p.compression_mm, c.indentation_nm / dc),
                xytext=(12, 14), textcoords="offset points", fontsize=9,
                arrowprops=dict(arrowstyle="-", color="#2b2b2b", lw=0.9))

    ax.set_yscale("log")
    ax.set_xlabel("wheel compression $T$  (mm)")
    ax.set_ylabel("$d\\,/\\,d_c$")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9, frameon=False, loc="lower right")
    ax.set_title("Every pad sits far inside the ductile region on indentation "
                 "alone -- the transition comes from ACCUMULATED work, not "
                 "depth", fontsize=10.5)
    fig.tight_layout()
    return fig


def demo(outdir: str = "_sagfig_demo") -> None:
    """Render all four on a real plan and check they came out."""
    import os

    import matplotlib
    matplotlib.use("Agg")
    plt = _plt()

    from . import sagdeck

    os.makedirs(outdir, exist_ok=True)
    plan = sagdeck.plan(sagdeck.SAGParams(grain_um=6.0, material="wc_co",
                                          micro_grains=1))
    made = []
    for fn in (load_collapse, contact_patch, dc_comparison, regime_map):
        fig = fn(plan)
        assert fig is not None, fn.__name__
        # every figure must actually have drawn something
        assert fig.get_axes(), fn.__name__
        for ax in fig.get_axes():
            assert (ax.lines or ax.patches or ax.collections), \
                "%s: an axis is empty" % fn.__name__
        # A legend with no entries means the labelled curves were never
        # drawn, which is exactly how the missing-attribute bug presented:
        # four plausible-looking figures with nothing plotted on them.
        for ax in fig.get_axes():
            h, _ = ax.get_legend_handles_labels()
            if ax.get_legend() is not None:
                assert h, "%s: a legend has no entries" % fn.__name__
        path = os.path.join(outdir, fn.__name__ + ".png")
        fig.savefig(path, dpi=100, bbox_inches="tight")
        plt.close(fig)
        assert os.path.getsize(path) > 12_000, path
        made.append((fn.__name__, os.path.getsize(path)))

    # the pad sweep must return real curves for every pad
    sweep = _pads(plan)
    assert set(sweep) == {6.0, 15.0, 30.0}, sorted(sweep)
    for dg, rows in sweep.items():
        assert len(rows) >= 10, (dg, len(rows))
        loads = [r[1].load_per_grain_n for r in rows]
        assert all(v > 0 for v in loads), dg
        # load must RISE with compression: more squash, more force
        assert loads[-1] > loads[0], dg
    # and a coarser pad must load each grain harder, since fewer share it
    assert (sweep[30.0][-1][1].load_per_grain_n
            > sweep[6.0][-1][1].load_per_grain_n), \
        "a coarse pad has fewer grains, so each must carry more"

    # the figures must agree with the plan they were drawn from
    c = plan["contact"]
    assert c.load_per_grain_n > 0 and c.indentation_nm > 0
    assert c.face_overrun > 1.0, \
        "the reference geometry clips; the patch figure says so"

    print("semgrit.sagfig: all checks passed")
    for n, b in made:
        print("  %-18s %6.1f KB" % (n, b / 1024))


if __name__ == "__main__":
    demo()
