"""Execute SEM_TO_ABAQUS_MULTI_ABRASIVE.ipynb headlessly, in a scratch dir.

Same idea as ``_run_notebook_test.py`` for the companion notebook: the notebook
is the deliverable, so it has to be run rather than eyeballed. The cells are
executed in order in a temporary directory with the dev copy of the package kept
off ``sys.path``, so what runs is the embedded payload and not the working tree.

Two substitutions, both unavoidable outside Colab: the image source is pointed at
a local ``.tif`` instead of an upload widget, and ``plt.show()`` writes a PNG.

    python _run_notebook2_test.py [image.tif]
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
NB = os.path.join(HERE, "SEM_TO_ABAQUS_MULTI_ABRASIVE.ipynb")

# Strings the run must produce. Each one is a claim that a whole stage worked,
# not a cosmetic check.
EXPECT = [
    "pipeline ready in",
    "modules   :",
    "grain library :",
    "SWEPT CHIP-THICHENSS" if False else "SWEPT CHIP-THICKNESS FIELD",
    "grits crossing the block",
    "elements the grits cut",
    "1/3  writing the deck",
    "2/3  sweeping the grit envelope",
    "3/3  injecting the field",
    "ALL",                       # a gate reporting all checks passed
    "to run it:",
]


def substitute(src: str, image: str, material: str = "") -> str:
    src = src.replace('SOURCE = "upload"', 'SOURCE = "already on disk"')
    src = src.replace('IMAGE_PATH = "/content/drive/MyDrive/sem/*.tif"',
                      'IMAGE_PATH = %r' % image)
    # No browser, so a figure is written instead of shown.
    src = src.replace("plt.show()",
                      "plt.gcf().savefig('fig%d.png' % _FIGN(), dpi=90)")
    src = src.replace("MA_DOWNLOAD = True", "MA_DOWNLOAD = False")
    if material:
        src = src.replace('MA_MATERIAL = "sandstone"',
                          'MA_MATERIAL = %r' % material)
    return src


MATERIAL = ""
"""Overrides the notebook's MA_MATERIAL dropdown, so the second material can be
run through the same pipeline: ``python _run_notebook2_test.py img.tif
silicon_carbide``."""


def main(argv) -> int:
    global MATERIAL
    image = argv[1] if len(argv) > 1 else os.path.join(HERE, "DIAMOND_11.tif")
    MATERIAL = argv[2] if len(argv) > 2 else ""
    if not os.path.exists(image):
        print("no such image: %s" % image)
        return 2
    if not os.path.exists(NB):
        print("no notebook: %s  (run _make_notebook2.py)" % NB)
        return 2

    cells = [c for c in json.load(open(NB, encoding="utf-8"))["cells"]
             if c["cell_type"] == "code"]
    scratch = tempfile.mkdtemp(prefix="nb2test_")
    print("#" * 78)
    print("# SEM_TO_ABAQUS_MULTI_ABRASIVE.ipynb   image=%s"
          % os.path.basename(image))
    print("# scratch dir: %s" % scratch)
    print("#" * 78)

    os.environ.setdefault("MPLBACKEND", "Agg")
    here_was = os.getcwd()
    # The dev tree must not be importable, or the test proves nothing about the
    # payload.
    saved_path = list(sys.path)
    sys.path[:] = [p for p in sys.path
                   if os.path.abspath(p or ".") != os.path.abspath(HERE)]
    for mod in [m for m in sys.modules
                if m.split(".")[0] in ("semgrit", "semgrit_multi")]:
        del sys.modules[mod]

    out_lines: list[str] = []
    fign = [0]

    class Tee:
        def write(self, s):
            out_lines.append(s)
            sys.__stdout__.write(s)

        def flush(self):
            sys.__stdout__.flush()

    g = {"__name__": "__main__", "_FIGN": lambda: fign.__setitem__(
        0, fign[0] + 1) or fign[0]}
    rc = 0
    try:
        os.chdir(scratch)
        sys.stdout = Tee()
        for n, cell in enumerate(cells, 1):
            src = substitute("".join(cell["source"]), image, MATERIAL)
            print("\n=============== CELL %d ===============" % n)
            exec(compile(src, "<cell %d>" % n, "exec"), g)
    except BaseException:
        sys.stdout = sys.__stdout__
        print("\nSTDERR:")
        traceback.print_exc()
        rc = 1
    finally:
        sys.stdout = sys.__stdout__
        os.chdir(here_was)
        sys.path[:] = saved_path

    text = "".join(out_lines)
    print("\n--- expectations ---")
    for want in EXPECT:
        ok = want in text
        print("  %-6s %s" % ("ok" if ok else "MISSING", want))
        if not ok:
            rc = 1
    print("\n--- artefacts ---")
    for root, _dirs, files in os.walk(scratch):
        for f in sorted(files):
            p = os.path.join(root, f)
            rel = os.path.relpath(p, scratch)
            if rel.startswith(("semgrit", "verify_", "vumat_", "_hybrid_test",
                               "_derive")):
                continue
            print("   %-58s %8.1f KB" % (rel, os.path.getsize(p) / 1024))

    print("=" * 78)
    print("MULTI-ABRASIVE NOTEBOOK: %s" % ("PASS" if rc == 0 else "FAILURE"))
    print("=" * 78)
    if rc == 0:
        shutil.rmtree(scratch, ignore_errors=True)
    else:
        print("scratch kept for inspection: %s" % scratch)
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv))
