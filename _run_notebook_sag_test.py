"""Execute SEM_TO_ABAQUS_SAG.ipynb's cells for real, in a clean directory.

    python _run_notebook_sag_test.py

AST-parsing a notebook proves it is valid Python. It does not prove the cells
RUN -- a wrong attribute name, a function that does not exist, a plan key that
was renamed all parse perfectly and fail at the first execution. This project
has been bitten by exactly that: ``sagdeck.summary`` and ``semgrit.sagfig``
were both referenced by the generator before either existed.

So this extracts the code cells, rewrites the interactive bits (uploads, the
huge MACRO deck, the paper build) into their cheapest equivalents, and executes
them in order in a temporary directory with the notebook's OWN embedded payload
-- not the working tree. If the payload is stale, this fails.

Matplotlib is forced to Agg and ``plt.show`` is stubbed, so figures are built
and thrown away rather than blocking.
"""
from __future__ import annotations

import io
import json
import os
import shutil
import sys
import tempfile
import traceback

NB = "SEM_TO_ABAQUS_SAG.ipynb"

# Cells rewritten for a headless, cheap run. Everything else executes verbatim.
REWRITE = {
    # cell 2: use the bundled images, never an upload dialog
    'SOURCE = "bundled"': 'SOURCE = "bundled"',
    # cell 6: MACRO is ~150 MB; the MICRO deck is what carries the physics
    "WRITE_MACRO = False": "WRITE_MACRO = False",
    # cell 10: building the paper's three decks takes many minutes
    "BUILD_PAPER = False": "BUILD_PAPER = False",
}

HEADER = """
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as _p
_p.show = lambda *a, **k: _p.close("all")


class _FakeDisplay(object):
    def __call__(self, *a, **k):
        return None


def display(*a, **k):
    return None


def HTML(x):
    return x
"""


def cells(path):
    nb = json.load(io.open(path, encoding="utf-8"))
    return [("".join(c["source"]), i)
            for i, c in enumerate(nb["cells"], 1)
            if c["cell_type"] == "code"]


def main():
    if not os.path.exists(NB):
        print("no %s -- run _make_notebook_sag.py first" % NB)
        return 1
    src = cells(NB)
    print("%s: %d code cells" % (NB, len(src)))

    tmp = tempfile.mkdtemp(prefix="_sagnb_")
    # The SEM images are inputs, not part of the payload, so they have to be
    # present for cell 2 -- same as a real Colab run with uploads.
    for f in sorted(os.listdir(".")):
        if f.startswith("B4C_1") and f.endswith(".tif"):
            shutil.copy2(f, tmp)
    here = os.getcwd()
    ns = {"__name__": "__main__"}
    bad = []
    try:
        os.chdir(tmp)
        exec(compile(HEADER, "<header>", "exec"), ns)
        for body, n in src:
            for a, b in REWRITE.items():
                if a in body:
                    body = body.replace(a, b)
            label = body.strip().split("\n", 1)[0][:66]
            try:
                exec(compile(body, "cell %d" % n, "exec"), ns)
                print("  [ok]   cell %-3d %s" % (n, label))
            except SystemExit as exc:
                print("  [stop] cell %-3d %s -- %s" % (n, label, exc))
                bad.append((n, "SystemExit: %s" % exc))
            except Exception as exc:                          # noqa: BLE001
                print("  [FAIL] cell %-3d %s" % (n, label))
                tb = traceback.format_exc().strip().split("\n")
                for ln in tb[-4:]:
                    print("         %s" % ln[:110])
                bad.append((n, "%s: %s" % (type(exc).__name__, exc)))
    finally:
        os.chdir(here)

    # What the run should have produced.
    print()
    checks = []
    for name, why in (("SOLIDS", "cell 3 measured the grains"),
                      ("PLAN", "cell 4 solved the contact"),
                      ("MICRO", "cell 6 wrote the MICRO deck")):
        checks.append((name in ns, "%s exists -- %s" % (name, why)))
    if "MICRO" in ns:
        mi = ns["MICRO"]
        # mi["path"] is relative to the temp dir the cells ran in, and we
        # have chdir'd back, so resolve it there rather than here.
        dp = mi["path"] if os.path.isabs(mi["path"])             else os.path.join(tmp, mi["path"])
        checks.append((os.path.exists(dp) and os.path.getsize(dp) > 1e6,
                       "the deck is on disk: %s (%.1f MB)"
                       % (os.path.basename(dp),
                          os.path.getsize(dp) / 1e6
                          if os.path.exists(dp) else 0.0)))
        checks.append((mi["elements"] > 1000,
                       "it has %s elements" % format(mi["elements"], ",")))
        checks.append((mi["swmode"] == 1,
                       "SWMODE = 1, the energy criterion"))
        checks.append((mi["n_passes"] >= 2,
                       "%d passes over one track" % mi["n_passes"]))
        checks.append((mi["dc_measured"],
                       "dc = %.1f nm is MEASURED" % mi["dc_nm"]))
    for ok, what in checks:
        print("  [%s] %s" % ("PASS" if ok else "FAIL", what))
        if not ok:
            bad.append((0, what))

    shutil.rmtree(tmp, ignore_errors=True)
    print()
    if bad:
        print("%d PROBLEM(S)" % len(bad))
        for n, why in bad:
            print("  cell %s: %s" % (n or "-", why))
        return 1
    print("ALL CELLS RAN")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
