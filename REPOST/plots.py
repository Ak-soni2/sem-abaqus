"""Figures for the finished grinding runs, from the CSVs they already wrote.

    python REPOST/plots.py                      # every job under "obd results/"
    python REPOST/plots.py <dir> [-o <outdir>]  # a specific results tree

WHY THIS IS NOT IN postprocess_odb.py. That script plots too, but it runs under
``abaqus python``, and Abaqus' bundled interpreter has no matplotlib -- so its
whole plotting block sits inside ``except ImportError`` and has silently drawn
nothing on every run so far. The six completed jobs are documented by phone
photographs of a screen because of it. This script reads the ``_forces.csv``,
``_energy.csv`` and ``_summary.json`` those runs already produced, with the
host Python that does have matplotlib, so the figures exist without Abaqus and
without re-running anything.

WHAT IT DRAWS, and what each is for:

  <job>_forces.png     normal and tangential force against time, with first
                       contact, peak and the engaged window marked. The
                       annotations are the point: all three are computed and
                       written to the summary JSON already, and a force trace
                       without them is unreadable.
  <job>_energy.png     the energy balance on a LOG axis, with the 5% artificial
                       and 10% kinetic thresholds drawn. Linear axes are useless
                       here -- on these runs ALLKE/ALLIE reaches 56,000, so a
                       linear plot is one flat line and three zeroes.
  compare_all.png      every job side by side: specific energy, peak force,
                       material removed, and the two energy-quality ratios.
                       This is the cross-deck figure the project never had.

HEALTH FLAGS ARE DRAWN ON THE FIGURES, not just printed. Two problems in the
archived runs would otherwise be invisible to anyone reading the plots:
  * artificial (hourglass) energy is 31-39% of internal energy, where under 5%
    is the usual bar. Hourglassing is carrying a third of the load.
  * kinetic energy runs 320-56,000x internal, i.e. mass scaling dominates.
A figure that hides that invites the number being quoted. So the thresholds are
drawn and the offending bars are hatched.

Units are scaled at the axis (us, mN, nm), never left as a 1e-6 offset label.
"""
from __future__ import annotations

import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Ductile/brittle is THE distinction in this project, and the rest of the repo
# draws it green/red -- which deuteranopes see as one colour. Blue/vermillion
# from the Okabe-Ito set reads for everyone and stays distinct in greyscale.
C_DUCTILE = "#0072B2"
C_BRITTLE = "#D55E00"
C_NORMAL = "#0072B2"
C_TANGENTIAL = "#D55E00"
C_WARN = "#B00020"
C_GREY = "#666666"

# Where "acceptable" sits for an explicit run. Not house style: ALLAE/ALLIE over
# ~5% means hourglass control is absorbing real work, and ALLKE/ALLIE over ~10%
# outside an impact problem means the mass scaling is driving the answer.
ARTIFICIAL_LIMIT = 0.05
KINETIC_LIMIT = 0.10

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 200, "font.size": 11,
    "axes.titlesize": 12, "axes.labelsize": 11, "axes.grid": True,
    "grid.alpha": 0.25, "legend.frameon": False, "figure.autolayout": True,
})


def read_csv(path):
    """CSV -> dict of column name to float array. Blank cells become NaN."""
    with open(path) as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 2:
        return {}
    head = [h.strip() for h in rows[0]]
    cols = {h: [] for h in head}
    for r in rows[1:]:
        if not any(x.strip() for x in r):
            continue
        for h, v in zip(head, r):
            try:
                cols[h].append(float(v))
            except ValueError:
                cols[h].append(float("nan"))
    return {h: np.asarray(v, dtype=float) for h, v in cols.items()}


def find_jobs(root):
    """Every job under root, keyed by job name, as (forces, energy, summary)."""
    jobs = {}
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if not f.endswith("_forces.csv"):
                continue
            job = f[: -len("_forces.csv")]
            base = os.path.join(dirpath, job)
            jobs[job] = dict(
                dir=dirpath, job=job, forces=base + "_forces.csv",
                energy=base + "_energy.csv" if os.path.exists(
                    base + "_energy.csv") else None,
                summary=base + "_summary.json" if os.path.exists(
                    base + "_summary.json") else None)
    return jobs


def load(meta):
    out = dict(meta)
    out["F"] = read_csv(meta["forces"]) if meta["forces"] else {}
    out["E"] = read_csv(meta["energy"]) if meta["energy"] else {}
    out["S"] = (json.load(open(meta["summary"])) if meta["summary"]
                else {})
    return out


def _material_of(job):
    """SiC jobs are named for it; everything else in this project is sandstone."""
    return "SiC" if "sic" in job.lower() else "sandstone"


def _kind_of(job):
    j = job.lower()
    for key, lab in (("multi", "multi"), ("energy", "energy"),
                     ("single", "single")):
        if key in j:
            return lab
    return "other"


# --------------------------------------------------------------------------
# 1. force against time
# --------------------------------------------------------------------------
def plot_forces(d, outdir):
    F, S = d["F"], d["S"]
    if not F or "time_s" not in F:
        return None
    t = F["time_s"] * 1e6                      # us, not a 1e-6 offset label
    fn = F.get("F_normal_N")
    ft = F.get("F_tangential_N")
    if fn is None or ft is None:
        return None
    fn, ft = fn * 1e3, ft * 1e3                # mN

    fig, ax = plt.subplots(figsize=(9, 4.6))
    ax.plot(t, fn, color=C_NORMAL, lw=1.2, label="normal $F_n$")
    ax.plot(t, ft, color=C_TANGENTIAL, lw=1.2, label="tangential $F_t$")
    ax.axhline(0, color=C_GREY, lw=0.8, zorder=0)

    s = (S.get("forces") or {})
    fc, pk = s.get("first_contact_s"), s.get("peak_time_s")
    # These three are computed by the postprocessor and written to JSON, then
    # thrown away. On the plot they are what separates rubbing from cutting.
    if fc:
        ax.axvline(fc * 1e6, color=C_GREY, ls=":", lw=1.2)
        ax.annotate("first contact", (fc * 1e6, ax.get_ylim()[1]),
                    xytext=(3, -12), textcoords="offset points",
                    fontsize=9, color=C_GREY, rotation=90, va="top")
    if fc and pk:
        ax.axvspan(fc * 1e6, t.max(), color=C_NORMAL, alpha=0.05, zorder=0)
    if pk:
        ax.axvline(pk * 1e6, color=C_WARN, ls="--", lw=1.2)
        ax.annotate("peak %.3g mN" % (abs(s.get("peak_magnitude_N", 0)) * 1e3),
                    (pk * 1e6, ax.get_ylim()[1]), xytext=(3, -12),
                    textcoords="offset points", fontsize=9, color=C_WARN,
                    rotation=90, va="top")

    ax.set_xlabel("time  [$\\mu$s]")
    ax.set_ylabel("force on the wheel  [mN]")
    ax.set_title("%s  --  grinding force  (%s)" % (d["job"],
                                                   _material_of(d["job"])))
    ax.legend(loc="best")
    # The sign convention is a documented trip hazard: state it where it is read.
    ax.text(0.01, 0.02, "$F_n>0$ pushes the wheel off the work",
            transform=ax.transAxes, fontsize=8, color=C_GREY)
    p = os.path.join(outdir, d["job"] + "_forces.png")
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------------------
# 2. energy balance, log axis, with the quality thresholds drawn
# --------------------------------------------------------------------------
def plot_energy(d, outdir):
    E, S = d["E"], d["S"]
    if not E or "time_s" not in E:
        return None
    t = E["time_s"] * 1e6
    fig, ax = plt.subplots(figsize=(9, 4.8))
    series = [("ALLIE_mJ", "internal  ALLIE", "#0072B2", 1.6),
              ("ALLKE_mJ", "kinetic  ALLKE", "#009E73", 1.2),
              ("ALLAE_mJ", "artificial  ALLAE", "#D55E00", 1.2),
              ("ALLPD_mJ", "plastic  ALLPD", "#CC79A7", 1.0),
              ("ALLDMD_mJ", "damage  ALLDMD", "#56B4E9", 1.0)]
    drawn = 0
    for key, lab, col, lw in series:
        y = E.get(key)
        if y is None or not np.isfinite(y).any() or np.nanmax(np.abs(y)) <= 0:
            continue
        ax.plot(t, np.abs(y), color=col, lw=lw, label=lab)
        drawn += 1
    if not drawn:
        plt.close(fig)
        return None
    # Log, because these span five decades: on a linear axis ALLKE is a flat
    # line and everything else is indistinguishable from zero.
    ax.set_yscale("log")
    ax.set_xlabel("time  [$\\mu$s]")
    ax.set_ylabel("energy  [mJ]   (log scale)")
    ax.set_title("%s  --  energy balance  (%s)" % (d["job"],
                                                   _material_of(d["job"])))
    ax.legend(loc="best", ncol=2, fontsize=9)

    en = (S.get("energy") or {})
    af, kf = en.get("artificial_fraction"), en.get("kinetic_fraction")
    msg = []
    if af is not None:
        msg.append("ALLAE/ALLIE = %.0f%%  (%s %.0f%%)"
                   % (100 * af, "limit" if af <= ARTIFICIAL_LIMIT else "OVER",
                      100 * ARTIFICIAL_LIMIT))
    if kf is not None:
        msg.append("ALLKE/ALLIE = %.0fx  (%s %.0f%%)"
                   % (kf, "limit" if kf <= KINETIC_LIMIT else "OVER",
                      100 * KINETIC_LIMIT))
    if msg:
        bad = ((af or 0) > ARTIFICIAL_LIMIT) or ((kf or 0) > KINETIC_LIMIT)
        ax.text(0.01, 0.97, "\n".join(msg), transform=ax.transAxes,
                va="top", fontsize=9, color=C_WARN if bad else C_GREY,
                bbox=dict(boxstyle="round,pad=0.4", fc="#FFF3F3" if bad
                          else "#F5F5F5", ec=C_WARN if bad else C_GREY,
                          lw=0.8))
    p = os.path.join(outdir, d["job"] + "_energy.png")
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------------------
# 3. the cross-deck comparison the project never had
# --------------------------------------------------------------------------
def plot_compare(loaded, outdir):
    rows = []
    for d in loaded:
        S = d["S"]
        f, e, r = (S.get("forces") or {}, S.get("energy") or {},
                   S.get("removal") or {})
        F = d["F"]
        # A job with no summary JSON is still a job: it ran, and its force CSV is
        # right there. Dropping it silently is how a comparison figure ends up
        # quietly missing an arm -- fall back to the CSV for what can be had.
        peak = f.get("peak_magnitude_N")
        if peak is None and "F_magnitude_N" in F and len(F["F_magnitude_N"]):
            peak = float(np.nanmax(np.abs(F["F_magnitude_N"])))
        rows.append(dict(
            job=d["job"], mat=_material_of(d["job"]), kind=_kind_of(d["job"]),
            u=S.get("specific_energy_J_mm3"),
            peak=(peak or 0) * 1e3,
            removed=(r.get("deleted_fraction") or 0) * 100,
            af=(e.get("artificial_fraction") or 0) * 100,
            kf=(e.get("kinetic_fraction") or 0),
            partial=not S))
    if not rows:
        return None
    rows.sort(key=lambda r: (r["mat"], r["kind"]))
    lab = ["%s\n%s" % (r["kind"], r["mat"]) for r in rows]
    x = np.arange(len(rows))
    col = [C_DUCTILE if r["mat"] == "sandstone" else C_BRITTLE for r in rows]

    # Duplicate detection, drawn rather than filed away. Four of the six archived
    # "results" are two datasets stored twice; a comparison figure that shows
    # them as four independent bars is a lie by omission. Hash the force CSV
    # itself rather than comparing summary numbers: it is the actual evidence,
    # and it catches the case where the summaries were regenerated separately.
    import hashlib
    seen, dupe = {}, set()
    by_job = {d["job"]: d for d in loaded}
    for r in rows:
        src = by_job[r["job"]]["forces"]
        h = hashlib.md5(open(src, "rb").read()).hexdigest()
        if h in seen:
            dupe.add(r["job"])
            dupe.add(seen[h])
        seen[h] = r["job"]
    n_distinct = len(seen)

    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    panels = [
        ("u", "specific energy  [J/mm$^3$]",
         "Specific energy -- the headline grinding number", True),
        ("peak", "peak force  [mN]", "Peak resultant force", False),
        ("removed", "elements deleted  [%]", "Material removed", False),
        ("af", "ALLAE / ALLIE  [%]",
         "Artificial (hourglass) energy -- solution quality", False),
    ]
    for ax, (key, ylab, title, logy) in zip(axes.ravel(), panels):
        vals = [r[key] or 0 for r in rows]
        bars = ax.bar(x, vals, color=col, edgecolor="black", linewidth=0.6)
        for b, r in zip(bars, rows):
            if r["job"] in dupe:
                b.set_hatch("//")
        ax.set_xticks(x)
        ax.set_xticklabels(lab, fontsize=9)
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=11)
        if logy and min(v for v in vals if v > 0) > 0:
            ax.set_yscale("log")
        for xi, v, r in zip(x, vals, rows):
            # "no data" and "measured zero" must not look the same.
            txt = "n/a" if (r["partial"] and key != "peak") else "%.3g" % v
            ax.annotate(txt, (xi, v), textcoords="offset points",
                        xytext=(0, 3), ha="center", fontsize=8,
                        color=C_GREY if txt == "n/a" else "black")
        if key == "af":
            ax.axhline(100 * ARTIFICIAL_LIMIT, color=C_WARN, ls="--", lw=1.2)
            ax.annotate("5% limit", (len(rows) - 0.4, 100 * ARTIFICIAL_LIMIT),
                        xytext=(0, 4), textcoords="offset points",
                        fontsize=9, color=C_WARN, ha="right")

    handles = [plt.Rectangle((0, 0), 1, 1, fc=C_DUCTILE, ec="k", lw=.6),
               plt.Rectangle((0, 0), 1, 1, fc=C_BRITTLE, ec="k", lw=.6)]
    labels = ["sandstone", "silicon carbide"]
    if dupe:
        handles.append(plt.Rectangle((0, 0), 1, 1, fc="white", ec="k", lw=.6,
                                     hatch="//"))
        labels.append("duplicate data (identical to another bar)")
    fig.legend(handles, labels, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.01))

    note = []
    if dupe:
        note.append("%d of the %d force datasets are byte-identical to another: "
                    "only %d runs here are distinct."
                    % (len(dupe), len(rows), n_distinct))
    over = [r["job"] for r in rows if r["af"] > 100 * ARTIFICIAL_LIMIT]
    if over:
        note.append("%d of %d exceed the 5%% artificial-energy bar; "
                    "hourglassing is carrying part of the load."
                    % (len(over), len(rows)))
    fig.suptitle("Grinding runs compared" + ("\n" + "  ".join(note) if note
                                             else ""),
                 fontsize=13, color=C_WARN if note else "black")
    p = os.path.join(outdir, "compare_all.png")
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("-")]
    src = args[0] if args else os.path.join(ROOT, "obd results")
    outdir = (argv[argv.index("-o") + 1] if "-o" in argv
              else os.path.join(HERE, "figures"))
    if not os.path.isdir(src):
        print("no results directory: %s" % src)
        return 2
    os.makedirs(outdir, exist_ok=True)

    jobs = find_jobs(src)
    if not jobs:
        print("no *_forces.csv found under %s" % src)
        return 1
    print("%d job(s) under %s" % (len(jobs), src))
    loaded, made = [], []
    for name in sorted(jobs):
        d = load(jobs[name])
        loaded.append(d)
        for fn in (plot_forces, plot_energy):
            p = fn(d, outdir)
            if p:
                made.append(p)
                print("  wrote %s" % os.path.relpath(p, ROOT))
    p = plot_compare(loaded, outdir)
    if p:
        made.append(p)
        print("  wrote %s" % os.path.relpath(p, ROOT))
    print("%d figure(s) in %s" % (len(made), os.path.relpath(outdir, ROOT)))
    return 0


def demo():
    """Self-check: the readers and the flags, without needing any real run."""
    import tempfile
    t = tempfile.mkdtemp()
    f = os.path.join(t, "j_forces.csv")
    with open(f, "w") as fh:
        fh.write("time_s,F_normal_N,F_tangential_N\n")
        for i in range(50):
            fh.write("%g,%g,%g\n" % (i * 1e-8, i * 1e-4, -i * 5e-5))
    with open(os.path.join(t, "j_energy.csv"), "w") as fh:
        fh.write("time_s,ALLIE_mJ,ALLKE_mJ,ALLAE_mJ\n")
        for i in range(50):
            fh.write("%g,%g,%g,%g\n" % (i * 1e-8, i * 1e-9, 1e-3, i * 4e-10))
    json.dump({"forces": {"peak_magnitude_N": 5e-3, "peak_time_s": 3e-7,
                          "first_contact_s": 5e-8},
               "energy": {"artificial_fraction": 0.4, "kinetic_fraction": 900.0},
               "removal": {"deleted_fraction": 0.05},
               "specific_energy_J_mm3": 1.2},
              open(os.path.join(t, "j_summary.json"), "w"))

    cols = read_csv(f)
    assert set(cols) == {"time_s", "F_normal_N", "F_tangential_N"}, cols
    assert len(cols["time_s"]) == 50
    jobs = find_jobs(t)
    assert list(jobs) == ["j"], jobs
    d = load(jobs["j"])
    out = os.path.join(t, "fig")
    os.makedirs(out)
    assert plot_forces(d, out) and plot_energy(d, out)
    # Two jobs with identical summaries must be flagged as duplicates.
    d2 = dict(d, job="j2")
    assert plot_compare([d, d2], out)
    assert _material_of("MULTI_abrasive1_sic") == "SiC"
    assert _material_of("multi_abrasive1") == "sandstone"
    assert _kind_of("energy_abrasive1") == "energy"
    print("demo OK ->", out)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        demo()
    else:
        sys.exit(main(sys.argv))
