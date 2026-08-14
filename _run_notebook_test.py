"""Execute the generated notebook headlessly in a clean directory, twice.

Simulates a fresh Colab session: a scratch directory containing nothing but the
notebook's own extracted payload, with the development copy of semgrit kept off
sys.path, so the run proves the notebook is genuinely self-contained.
"""
import json
import os
import re
import matplotlib
matplotlib.use('Agg')
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
NB = os.path.join(HERE, 'SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb')

RUNS = [
    # Run 1 is what a new user gets: SIMPLE MODE first (its defaults are on), then
    # the advanced cells over the top. Both paths in one pass.
    ("run1_simple_then_advanced", os.path.join(HERE, 'DIAMOND_11.tif'), {
        'SOURCE = "upload"': 'SOURCE = "already on disk"',
        'S_NAME = "wheel"': 'S_NAME = "run1_simple"',
        'S_SLICE_MM = 2.0': 'S_SLICE_MM = 1.0',
        'S_GRITS = "concentration"': 'S_GRITS = "a fixed number"',
        'S_GRIT_VALUE = 100.0': 'S_GRIT_VALUE = 40.0',
        'S_POSITION = "centred"': 'S_POSITION = "first grit at entry"',
        'MODEL_NAME = "wheel"': 'MODEL_NAME = "run1_wheel"',
        'WRITE_RUN_READY_INP = True': 'WRITE_RUN_READY_INP = False',  # CAE-only path
    }),
    ("run2_singlegrit_30deg_cad", os.path.join(HERE, 'B4C_15.tif'), {
        'SOURCE = "upload"': 'SOURCE = "already on disk"',
        'RUN_SIMPLE = True': 'RUN_SIMPLE = False',      # advanced path only
        'WORK_SPEED_MM_S = 0.0': 'WORK_SPEED_MM_S = 250.0',
        'CAD_MODE = "whole wheel"': 'CAD_MODE = "contact"',
        'MODEL_NAME = "wheel"': 'MODEL_NAME = "run2_wheel"',
        'SECTOR_MODE = "arc"': 'SECTOR_MODE = "angle"',
        'SECTOR_DEG = 30.0': 'SECTOR_DEG = 3.0',
        'GRIT_MODE = "concentration"': 'GRIT_MODE = "single"',
        'EDGE_RADIUS_UM = 0.0': 'EDGE_RADIUS_UM = 0.25',
        'WP_ELEMENT_SIZE_MM = 0.0003': 'WP_ELEMENT_SIZE_MM = 0.0006',
        # anisotropic mesh: coarse axially, coarser still through the depth
        'WP_ELEM_AXIAL_MM = 0.0': 'WP_ELEM_AXIAL_MM = 0.0015',
        'WP_ELEM_DEPTH_MM = 0.0': 'WP_ELEM_DEPTH_MM = 0.0012',
        # every artefact at once: both decks, wheel CAD, grit CAD
        'WRITE_WHEEL_STEP = False': 'WRITE_WHEEL_STEP = True',
        'WRITE_WHEEL_STL = False': 'WRITE_WHEEL_STL = True',
        'WRITE_GRAINS_STEP = False': 'WRITE_GRAINS_STEP = True',
        'WRITE_GRAIN_STLS = False': 'WRITE_GRAIN_STLS = True',
        'CORES = 8': 'CORES = 4',
        # run-ready path: a submittable deck, plus a shallower cut for this wheel
        'DEPTH_OF_CUT_UM = 3.0': 'DEPTH_OF_CUT_UM = 1.2',
        'MASS_SCALING = 10.0': 'MASS_SCALING = 20.0',
        'FIELD_FRAMES = 60': 'FIELD_FRAMES = 30',
        'VIEW_MODE = "contact"': 'VIEW_MODE = "wheel"',   # exercise both view modes
        # park the block clear of the grits, at the leading grain, and let the
        # automatic depth of cut close the standoff again
        'WP_POSITION = "centred"': 'WP_POSITION = "first grit at entry"',
        'CLEARANCE_UM = 0.0': 'CLEARANCE_UM = 0.6',
        'DEPTH_OF_CUT_UM = 1.2': 'DEPTH_OF_CUT_UM = 0.0',
        # A CAD-viewer edit, arriving the way a pasted Copy JSON does. It must reach the
        # deck: A13 rebuilt PARAMS from the widgets and dropped every edit on the floor.
        'PASTED_SETTINGS = ""': 'PASTED_SETTINGS = \'{"clearance_um": 0.45}\'',
        # Section B on the SECOND material. Run 1 covers sandstone, so
        # between the two runs both cards go through the writer and the
        # gate. Cheaper than a third full pass over the notebook, and the
        # plumbing it exercises -- dropdown -> registry -> JH-2 card ->
        # *Material name -> the 56 props -- is exactly where a per-material
        # bug would live.
        'SA_MATERIAL = "sandstone"': 'SA_MATERIAL = "silicon_carbide"',
    }),
]


def cells_source(path):
    nb = json.load(open(path, encoding='utf-8'))
    # nbformat concatenates the source list verbatim -- exactly what Jupyter does.
    return [c['source'] if isinstance(c['source'], str) else ''.join(c['source'])
            for c in nb['cells'] if c['cell_type'] == 'code']


def main():
    src = cells_source(NB)
    print('notebook has %d code cells' % len(src))
    overall = True

    for name, image, subs in RUNS:
        work = os.path.join(tempfile.gettempdir(), 'nbtest_' + name)
        shutil.rmtree(work, ignore_errors=True)
        os.makedirs(work)
        body = list(src)
        for i, cell in enumerate(body):
            for a, b in subs.items():
                cell = cell.replace(a, b)
            cell = cell.replace('IMAGE_PATH = "/content/drive/MyDrive/sem/*.tif"',
                                'IMAGE_PATH = r"%s"' % image)
            # Headless: there is no display, so pick a non-interactive backend and
            # save the preview to disk instead of calling plt.show().
            cell = cell.replace(
                'import matplotlib.pyplot as plt',
                'import matplotlib\nmatplotlib.use("Agg")\n'
                'import matplotlib.pyplot as plt')
            # headless: write the interactive view to HTML instead of opening it
            cell = cell.replace(
                '    FIG3D.show()',
                '    FIG3D.write_html("cad_view.html", include_plotlyjs="cdn")\n'
                '    print("3D CAD view saved: cad_view.html")')
            cell = cell.replace(
                '    plt.show()',
                '    fig.savefig("preview.png", dpi=100)\n'
                '    print("preview figure saved: preview.png")')
            body[i] = cell

        runner = os.path.join(work, 'runner.py')
        with open(runner, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('import sys, os\n')
            for i, cell in enumerate(body):
                fh.write('\nprint("\\n%s CELL %d %s")\n' % ('=' * 26, i + 1, '=' * 26))
                fh.write('exec(compile(%r, "<cell %d>", "exec"), globals())\n'
                         % (cell, i + 1))

        print('\n' + '#' * 78)
        print('# %s   image=%s' % (name, os.path.basename(image)))
        print('# scratch dir: %s' % work)
        print('#' * 78)
        env = dict(os.environ)
        env['PYTHONPATH'] = ''          # keep the dev copy of semgrit off the path
        r = subprocess.run([sys.executable, 'runner.py'], cwd=work, env=env,
                           capture_output=True, text=True)
        # Enough to hold the build report and both verifiers; a shorter tail hid the
        # cell that failed.
        print(r.stdout[-40000:])
        if r.stderr.strip():
            print('STDERR:\n' + r.stderr[-4000:])
        ok = r.returncode == 0 and 'DECK(S) GOOD' in r.stdout
        # Simple mode must genuinely have produced a deck, not just not crashed.
        if 'RUN_SIMPLE = False' not in str(subs):
            # simple mode must show before it builds, and build only in its own cell
            ok = ok and 'nothing written yet' in r.stdout
            ok = ok and 'run1_simple.inp' in r.stdout
            ok = ok and 'reusing' in r.stdout          # the grain cache is used
            ok = ok and 'GRAINS, counted off the deck' in r.stdout
            # No edits is the ordinary case and must not stop the run before A13.
            ok = ok and 'no edits found, so the build below' in r.stdout
            ok = ok and 'settings: the form widgets above' in r.stdout
        if 'PASTED_SETTINGS = ""' in subs:
            # An edit applied in A12b has to be the thing A13 writes, and the notebook has
            # to reach A13 at all -- A12b used to raise SystemExit on the ordinary "no
            # edits" case, which killed the run before the deck was ever built.
            for want in ('applied: clearance_um',
                         "the CAD viewer's edits from A12b",
                         'at a 0.450 um standoff'):
                if want not in r.stdout:
                    print('MISSING (viewer edit did not reach the deck): %r' % want)
                    ok = False
        produced = []
        for root, _, fs in os.walk(work):
            for f in fs:
                if f.endswith(('.inp', '.step', '.stl', '.json', '.csv', '.py', '.zip',
                               '.png', '.html', '.glb')) \
                        and 'semgrit' not in root and f != 'runner.py':
                    produced.append((os.path.relpath(os.path.join(root, f), work),
                                     os.path.getsize(os.path.join(root, f))))
        print('--- artefacts produced ---')
        for p, s in sorted(produced):
            print('   %-52s %9.1f KB' % (p, s / 1024))
        print('RESULT: %s' % ('PASS' if ok else 'FAIL'))
        overall = overall and ok

    print('\n' + '=' * 78)
    print('NOTEBOOK END-TO-END: %s' % ('ALL RUNS PASS' if overall else 'FAILURE'))
    print('=' * 78)
    return 0 if overall else 1


if __name__ == '__main__':
    sys.exit(main())
