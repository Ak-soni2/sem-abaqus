"""Full verification of the Colab notebook: structure, wiring, every feature, every mode.

The deck verifiers check a finished .inp. The pipeline verifiers check the measurement
half. Neither can catch the failures that live in the notebook itself:

  * the embedded payload going stale after the package is edited
  * a form widget wired to nothing, so turning it does nothing at all
  * a package feature with no widget, so it is unreachable from the notebook
  * a cell using a name an earlier cell never defines
  * a mode nobody exercised -- full wheel, wheel-only, single grit, uniform mesh
  * the preview disagreeing with the deck it claims to predict

  C1  structure    JSON, every cell parses, nbformat line endings, name ordering
  C2  payload      the blob matches the package on disk, file for file
  C3  wiring       every widget is read; every DeckParams/AnalysisParams field reachable
  C4  features     a matrix of modes, each built and put through both deck verifiers
  C5  sensitivity  changing a parameter actually changes the deck it should change
  C6  preview      plan_deck agrees with build_deck on every number it displays
  C7  determinism  same settings twice -> identical deck

Exits non-zero on any failure.
"""
import ast
import base64
import dataclasses
import glob
import gzip
import io
import json
import math
import os
import pickle
import re
import shutil
import subprocess
import sys
import tarfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

NB = 'SEM_TO_ABAQUS_GRINDING_WHEEL.ipynb'
OUT = '_colabcheck'
FAIL, WARN = [], []


def chk(name, ok, detail=''):
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', name,
                           (': ' + detail) if detail else ''))
    if not ok:
        FAIL.append(name)


def warn(name, ok, detail=''):
    print('  [%s] %s%s' % ('ok  ' if ok else 'WARN', name,
                           (': ' + detail) if detail else ''))
    if not ok:
        WARN.append(name)


def code_cells(nb):
    return [''.join(c['source']) if isinstance(c['source'], list) else c['source']
            for c in nb['cells'] if c['cell_type'] == 'code']


def main():
    print('=' * 78)
    print('COLAB NOTEBOOK  full verification')
    print('=' * 78)

    # ---------------- C1 structure ----------------
    print('C1  STRUCTURE')
    nb = json.load(open(NB, encoding='utf-8'))
    cells = code_cells(nb)
    chk('notebook is valid JSON, nbformat 4', nb.get('nbformat') == 4,
        'nbformat %s, %d cells' % (nb.get('nbformat'), len(nb['cells'])))
    bad = []
    for i, c in enumerate(cells, 1):
        try:
            ast.parse(c)
        except SyntaxError as e:
            bad.append((i, str(e)[:60]))
    chk('every code cell parses as Python', not bad, str(bad[:2]))
    # nbformat concatenates source lists verbatim; without trailing newlines every
    # cell collapses onto one line in real Jupyter even though our test runner copes.
    nonl = []
    for j, c in enumerate(nb['cells']):
        src = c['source']
        if isinstance(src, list) and len(src) > 1:
            if any(not l.endswith('\n') for l in src[:-1]):
                nonl.append(j)
    chk('source lines keep their trailing newlines', not nonl, str(nonl[:3]))

    # names must be defined before the cell that reads them
    defined, undef = set(), []
    SAFE = {'WORK', 'PAYLOAD', 'SOLIDS', 'IMAGES', 'PARAMS', 'PLAN', 'INFO',
            'ANALYSIS', 'OUT_DECK', 'OUT_MEAS', 'ALL_GRAINS', 'MODEL', 'SEARCH_ROOTS',
            'INP_NAME'}
    for i, c in enumerate(cells, 1):
        # Strip string and comment content first. The embedded base64 payload is one
        # long run of uppercase and digits, every fragment of which looks like an
        # undefined constant to a bare regex.
        clean = re.sub(r'"[^"\n]*"', '""', c)
        clean = re.sub(r"'[^'\n]*'", "''", clean)
        clean = re.sub(r'#[^\n]*', '', clean)
        used = set(re.findall(r'\b([A-Z][A-Z0-9_]{2,})\b', clean))
        # Names this cell binds itself, by assignment or by import. Collected before
        # the check, because a cell may legitimately import and use in one go.
        # Walk the AST rather than pattern-matching assignments: a regex misses tuple
        # unpacking (FIG3D, _drew = ...), for-targets and with-targets.
        local = set()
        try:
            for n in ast.walk(ast.parse(c)):
                if isinstance(n, (ast.Import, ast.ImportFrom)):
                    local |= {(a.asname or a.name).split('.')[0] for a in n.names}
                elif isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store):
                    local.add(n.id)
                elif isinstance(n, (ast.FunctionDef, ast.ClassDef)):
                    local.add(n.name)
        except SyntaxError:
            local |= set(re.findall(r'^\s*([A-Z][A-Z0-9_]*)\s*=', c, re.M))
        for u in sorted(used - defined - SAFE - local):
            if u in {'YES', 'NO', 'ALL', 'SPOS', 'ENCASTRE', 'STONE', 'JH2'}:
                continue
            undef.append((i, u))
        defined |= local
    chk('no cell reads a NAME an earlier cell never defined', not undef,
        str(undef[:4]))

    # ---------------- C2 payload freshness ----------------
    print()
    print('C2  EMBEDDED PAYLOAD')
    setup = next(c for c in cells if 'PAYLOAD = (' in c)
    blob = ''.join(re.findall(r'^\s*"([A-Za-z0-9+/=]+)"\s*$', setup, re.M))
    chk('payload blob found in the setup cell', len(blob) > 1000,
        '%.0f KB of base64' % (len(blob) / 1024))
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    with tarfile.open(fileobj=io.BytesIO(
            gzip.decompress(base64.b64decode(blob)))) as tf:
        tf.extractall(OUT)
        members = [m.name for m in tf.getmembers()]
    on_disk = sorted(glob.glob('semgrit/*.py')) + ['verify_rigid_deck.py',
                                                   'verify_rigid_deck2.py',
                                                   'verify_pipeline_A.py']
    chk('payload contains every package file',
        sorted(members) == sorted(f.replace(os.sep, '/') for f in on_disk),
        '%d files' % len(members))
    stale = [f for f in on_disk
             if open(f, 'rb').read() != open(os.path.join(OUT, f), 'rb').read()]
    # The single most likely way this notebook breaks: the package is edited and the
    # notebook is not regenerated, so Colab silently runs an older pipeline.
    chk('payload is up to date with the package on disk', not stale,
        'stale: %s' % (stale[:4] if stale else 'none'))

    # ---------------- C3 wiring ----------------
    print()
    print('C3  WIRING  (no dead widgets, no unreachable features)')
    widgets = []
    for c in cells:
        widgets += re.findall(r'^([A-Z][A-Z0-9_]*)\s*=.*#@param', c, re.M)
    body = '\n'.join(cells)
    dead = []
    for w in widgets:
        # count uses that are not the definition itself
        uses = len(re.findall(r'\b%s\b' % w, body)) - len(
            re.findall(r'^%s\s*=.*#@param' % w, body, re.M))
        if uses == 0:
            dead.append(w)
    chk('every form widget is read somewhere', not dead,
        '%d widgets, dead: %s' % (len(widgets), dead or 'none'))

    # Every third-party module the package imports must be in the notebook's install
    # list. mapbox_earcut is not preinstalled on Colab, and without it every grain
    # fails to triangulate -- which presents as an empty library and reads like a
    # segmentation problem, so it is worth checking structurally rather than hoping.
    STDLIB = set(sys.stdlib_module_names) | {'semgrit', '__future__'}
    imported = set()
    for f in glob.glob('semgrit/*.py'):
        tree = ast.parse(open(f, encoding='utf-8').read())
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                imported |= {a.name.split('.')[0] for a in n.names}
            elif isinstance(n, ast.ImportFrom) and n.level == 0 and n.module:
                imported.add(n.module.split('.')[0])
    third = sorted(m for m in imported - STDLIB
                   if m not in {'abaqus', 'abaqusConstants'})
    setup_cell = next(c for c in cells if 'pip' in c and 'missing' in c)
    ensured = set(re.findall(r'\("([A-Za-z_0-9]+)",\s*"[^"]+"\)', setup_cell))
    gap = [m for m in third if m not in ensured]
    chk('every third-party import is installed by the setup cell', not gap,
        'imports %s; not ensured: %s' % (third, gap or 'none'))

    from semgrit.analysis import AnalysisParams
    from semgrit.build_deck import DeckParams
    mp = body[body.index('def make_params'):]
    for cls, skip in ((DeckParams, {'analysis'}), (AnalysisParams, {'enabled'})):
        miss = [f.name for f in dataclasses.fields(cls)
                if f.name not in skip and not re.search(r'\b%s\s*=' % f.name, mp)]
        chk('every %s field is set from the notebook' % cls.__name__, not miss,
            'unreachable: %s' % (miss or 'none'))

    # ---------------- C4 feature matrix ----------------
    print()
    print('C4  FEATURE MATRIX  (each built, then both deck verifiers)')
    solids = pickle.load(open('WHEEL_FIXED/1_measurements/grain_library.pkl',
                              'rb'))['solids']
    from semgrit.build_deck import build_deck
    base = dict(diameter_mm=50.0, rim_depth_mm=2.0, width_mm=0.6,
                shell_circumferential_divisions=60, shell_axial_divisions=4,
                shell_radial_divisions=2, grit_arc_window_mm=0.4,
                grit_width_window_mm=0.3, seed=3,
                wp_length_mm=0.06, wp_width_mm=0.03, wp_depth_mm=0.03,
                wp_element_size_mm=0.004, travel_mm=0.02,
                surface_speed_mm_s=30000.0, cores=8)
    AN = dict(enabled=True, depth_of_cut_um=0.8, mass_scaling_factor=10.0,
              field_frames=20, restart_intervals=5)
    cases = {
        'arc + concentration':       dict(base, sector_mode='arc', arc_length_mm=1.0,
                                          grit_mode='concentration', concentration=60.0),
        'angle + areal density':     dict(base, sector_mode='angle', sector_deg=20.0,
                                          grit_mode='areal_density',
                                          areal_density_per_mm2=1500.0),
        'full wheel + count':        dict(base, sector_mode='full', grit_mode='count',
                                          grit_count=120, grit_arc_window_mm=0.0),
        'single grit':               dict(base, sector_mode='arc', arc_length_mm=1.0,
                                          grit_mode='single', single_grain_index=-1,
                                          single_grit_offset_mm=0.005,
                                          grit_arc_window_mm=0.0),
        'wheel only (no workpiece)': dict(base, sector_mode='arc', arc_length_mm=1.0,
                                          grit_mode='count', grit_count=60,
                                          include_workpiece=False),
        'graded depth mesh':         dict(base, sector_mode='arc', arc_length_mm=1.0,
                                          grit_mode='count', grit_count=60,
                                          wp_element_size_depth_mm=0.001,
                                          wp_surface_layer_mm=0.004,
                                          wp_depth_growth=1.4,
                                          wp_max_depth_element_mm=0.008),
        'anisotropic mesh':          dict(base, sector_mode='arc', arc_length_mm=1.0,
                                          grit_mode='count', grit_count=60,
                                          wp_element_size_width_mm=0.01,
                                          wp_element_size_depth_mm=0.006),
        'no inset band':             dict(base, sector_mode='arc', arc_length_mm=1.0,
                                          grit_mode='count', grit_count=60,
                                          inset_grit_band=False),
    }
    built = []
    for i, (label, kw) in enumerate(cases.items()):
        for ready in (False, True):
            if ready and not kw.get('include_workpiece', True):
                continue
            nm = 'c%d%s' % (i, 'r' if ready else 'g')
            try:
                info = build_deck(DeckParams(
                    name=nm, analysis=(AnalysisParams(**AN) if ready else None),
                    also_write_cae_deck=ready, **kw), solids, OUT)
            except Exception as e:
                chk('%s [%s]' % (label, 'run-ready' if ready else 'geometry'),
                    False, '%s: %s' % (type(e).__name__, e))
                continue
            built.append((label, ready, info))
    print('      built %d decks across %d feature combinations'
          % (len(built), len(cases)))
    bad = []
    total = 0
    for label, ready, info in built:
        decks = [info['path']] + ([info['cae_deck']] if info.get('cae_deck') else [])
        for d in decks:
            for v in ('verify_rigid_deck.py', 'verify_rigid_deck2.py'):
                r = subprocess.run([sys.executable, v, d], capture_output=True,
                                   text=True)
                total += r.stdout.count('[PASS]')
                if r.returncode != 0:
                    bad.append((label, os.path.basename(d), v,
                                [l.strip() for l in r.stdout.split('\n')
                                 if 'FAIL]' in l][:2]))
    chk('every deck in the matrix passes both verifiers', not bad,
        '%d checks across %d decks' % (total, sum(
            1 + (1 if i.get('cae_deck') else 0) for _, _, i in built)))
    for b in bad[:5]:
        print('      %s / %s / %s: %s' % b)

    # ---------------- C5 sensitivity ----------------
    print()
    print('C5  SENSITIVITY  (a setting that changes nothing is a broken setting)')
    sens = {
        'diameter_mm': (dict(diameter_mm=40.0), 'outer_radius_mm'),
        'rim_depth_mm': (dict(rim_depth_mm=3.0), 'rim_depth_mm'),
        'width_mm': (dict(width_mm=0.9), None),
        'arc_length_mm': (dict(arc_length_mm=1.5), 'arc_length_mm'),
        'grit_count': (dict(grit_count=90), 'n_grits'),
        'seed': (dict(seed=99), None),
        'wp_length_mm': (dict(wp_length_mm=0.09), None),
        'wp_element_size_mm': (dict(wp_element_size_mm=0.003),
                               'n_workpiece_elements'),
        'clearance_um': (dict(clearance_um=0.5), 'workpiece_ground_radius_mm'),
        'protrusion_mean': (dict(protrusion_mean=0.75), None),
        'wp_position': (dict(wp_position='first grit at entry'),
                        'theta_workpiece_deg'),
        'wp_position_deg': (dict(wp_position='custom angle', wp_position_deg=1.0),
                            'theta_workpiece_deg'),
        'max_tilt_deg': (dict(max_tilt_deg=10.0), None),
    }
    ref_kw = dict(base, sector_mode='arc', arc_length_mm=1.0, grit_mode='count',
                  grit_count=60)
    ref = build_deck(DeckParams(name='s_ref', **ref_kw), solids, OUT)
    ref_txt = open(ref['path'], encoding='ascii').read()
    for nm, (delta, key) in sens.items():
        kw = dict(ref_kw)
        kw.update(delta)
        got = build_deck(DeckParams(name='s_' + nm, **kw), solids, OUT)
        txt = open(got['path'], encoding='ascii').read()
        changed = txt != ref_txt
        keyed = (key is None) or (abs(float(got[key]) - float(ref[key])) > 1e-12)
        chk('%s changes the deck' % nm, changed and keyed,
            ('%s: %s -> %s' % (key, ref[key], got[key])) if key else 'file differs')

    an_sens = {'depth_of_cut_um': 0.5, 'mass_scaling_factor': 25.0,
               'friction': 0.5, 'field_frames': 11, 'restart_intervals': 7,
               'n_depvar': 14, 'contact_scope': 'all exterior',
               'hourglass': 'STIFFNESS', 'nlgeom': False}
    r0 = build_deck(DeckParams(name='a_ref', analysis=AnalysisParams(**AN),
                               **ref_kw), solids, OUT)
    t0 = open(r0['path'], encoding='ascii').read()
    for nm, val in an_sens.items():
        a = dict(AN)
        a[nm] = val
        g = build_deck(DeckParams(name='a_' + nm, analysis=AnalysisParams(**a),
                                  **ref_kw), solids, OUT)
        chk('analysis.%s changes the deck' % nm,
            open(g['path'], encoding='ascii').read() != t0)

    # ---------------- C6 preview fidelity ----------------
    print()
    print('C6  PREVIEW vs BUILD  (the preview must not promise a different model)')
    from semgrit.build_deck import plan_deck
    pp = DeckParams(name='pv', analysis=AnalysisParams(**AN), **ref_kw)
    plan = plan_deck(pp, solids)
    real = build_deck(pp, solids, OUT)
    for label, a, b, tol in (
            ('grit count', plan['n_grits'], real['n_grits'], 0),
            ('workpiece elements', plan['n_workpiece_elements'],
             real['n_workpiece_elements'], 0),
            ('theta of the block', plan['theta_workpiece_deg'],
             real['theta_workpiece_deg'], 1e-9),
            ('ground-face radius', plan['ground_radius_mm'],
             real['workpiece_ground_radius_mm'], 1e-9),
            ('bond clearance', plan['bond_clearance_um'],
             real['max_engaging_protrusion_um'], 1e-9),
            ('stable dt', plan['cost']['stable_dt_s'],
             real['cost']['stable_dt_s'], 1e-18),
            ('dressed band arc', plan['grit_band_arc_mm'],
             real['grit_band_arc_mm'], 1e-12)):
        chk('preview agrees on %s' % label, abs(float(a) - float(b)) <= tol,
            '%s vs %s' % (a, b))
    warn('preview file-size estimate is within 30%% of the deck written',
         abs(plan['estimated_mb'] - real['size_bytes'] / 1e6)
         <= 0.30 * real['size_bytes'] / 1e6,
         'estimated %.2f MB, actual %.2f MB'
         % (plan['estimated_mb'], real['size_bytes'] / 1e6))

    # ---------------- C6b the 3-D viewer ----------------
    print()
    print('C6b VIEWER vs DECK  (the 3-D view must be the deck, not a likeness)')
    from scipy.spatial import cKDTree
    from semgrit.viewer import view3d
    fig, drew = view3d(plan, mode='wheel', max_grits=10 ** 6)
    txt = open(real['path'], encoding='ascii').read()
    NL = chr(10)
    wb = re.search(r'\*Part, name=WHEEL\s*$(.*?)\*End Part', txt,
                   re.S | re.M).group(1)

    def _count(kind):
        blk = re.search(r'\*Element, type=%s[^%s]*%s(.*?)(?=\*)' % (kind, NL, NL),
                        wb, re.S).group(1)
        return len(blk.strip().split(NL))

    n3, n4 = _count('R3D3'), _count('R3D4')
    chk('viewer draws every grit facet in the deck',
        drew['grit_triangles'] == n3, '%d vs %d' % (drew['grit_triangles'], n3))
    chk('viewer draws every bond quad, split into triangles',
        drew['bond_triangles'] == 2 * n4,
        '%d vs 2 x %d' % (drew['bond_triangles'], n4))
    nb = re.search(r'\*Node\s*$(.*?)(?=\*Element)', wb, re.S | re.M).group(1)
    deck_nodes = np.array([[float(x) for x in l.split(',')[1:4]]
                           for l in nb.strip().split(NL)])
    gm = [t for t in fig.data if getattr(t, 'name', '') == 'abrasive grits'][0]
    V = np.column_stack([gm.x, gm.y, gm.z])
    d, _ = cKDTree(deck_nodes).query(V)
    chk('every grit vertex drawn is a node of the deck', float(d.max()) < 1e-12,
        '%d vertices, worst offset %.3e mm' % (len(V), d.max()))
    wpm = [t for t in fig.data if getattr(t, 'name', '') == 'workpiece']
    if wpm:
        box = np.column_stack([wpm[0].x, wpm[0].y, wpm[0].z])
        pb = re.search(r'\*Part, name=WORKPIECE\s*$\s*\*Node\s*$(.*?)(?=\*Element)',
                       txt, re.S | re.M).group(1)
        W = np.array([[float(x) for x in l.split(',')[1:4]]
                      for l in pb.strip().split(NL)])
        chk('workpiece box matches the meshed block',
            bool(np.abs((box.max(0) - box.min(0)) - (W.max(0) - W.min(0))).max() < 1e-9),
            'span %s' % np.round(box.max(0) - box.min(0), 6))
    fc, dc = view3d(plan, mode='contact', max_grits=50)
    chk('contact view clips and caps as asked',
        dc['grits_drawn'] <= 50 and dc['grits_drawn'] <= dc['grits_in_view']
        <= dc['grits_total'],
        'drawn %d, in view %d, total %d' % (dc['grits_drawn'], dc['grits_in_view'],
                                            dc['grits_total']))

    # ---------------- C6c the glTF export ----------------
    print()
    print('C6c glTF vs DECK  (the free CAD view must also be the deck)')
    from semgrit.glb import model_viewer_html, parts_from_plan, write_glb
    gp = os.path.join(OUT, 'check.glb')
    gi = write_glb(gp, parts_from_plan(plan, mode='wheel', max_grits=10 ** 6))
    chk('glb written with bond, grits and workpiece', gi['parts'] == 3,
        '%d parts, %.2f MB' % (gi['parts'], gi['bytes'] / 1e6))
    import trimesh
    sc = trimesh.load(gp)
    chk('an independent glTF reader parses it',
        set(sc.geometry) == {'bond rim', 'abrasive grits', 'workpiece'},
        str(sorted(sc.geometry)))
    g = sc.geometry['abrasive grits']
    chk('grit solids survive the round trip closed', bool(g.is_watertight),
        '%d faces' % len(g.faces))
    chk('glb grit triangles match the deck', len(g.faces) == n3,
        '%d vs %d' % (len(g.faces), n3))
    # undo the Y-up convention the writer applies
    GV = np.column_stack([g.vertices[:, 0], -g.vertices[:, 2], g.vertices[:, 1]])
    gd, _ = cKDTree(deck_nodes).query(GV)
    # glTF mandates float32 positions, so ~1.5e-6 mm at radius 25 is the floor
    chk('every glb vertex is a deck node, to float32 precision',
        float(gd.max()) < 5e-6, 'worst offset %.3e mm' % gd.max())
    html = model_viewer_html(gp, max_inline_mb=999)
    chk('model-viewer html embeds the model with no upload',
        'model-viewer' in html and 'data:model/gltf-binary' in html
        and 'http' not in html.split('src="data:')[1][:40])
    try:
        model_viewer_html(gp, max_inline_mb=1e-6)
        chk('oversize glb refused rather than inlined', False)
    except ValueError as exc:
        chk('oversize glb refused rather than inlined', 'too large to inline' in str(exc))

    # ---------------- C6d the three.js CAD viewer ----------------
    print()
    print('C6d CAD VIEWER  (the engineering tools must report the deck, not a copy)')
    from semgrit import cadviewer
    cvp = os.path.join(OUT, 'cadview.glb')
    cad_html, cmeta, cinfo = cadviewer.build(plan, cvp, mode='wheel',
                                             max_grits=10 ** 6, max_inline_mb=999)
    chk('the viewer builds bond, grits and workpiece', cinfo['parts'] == 3,
        '%d parts, %.2f MB, %s triangles'
        % (cinfo['parts'], cinfo['bytes'] / 1e6, format(cinfo['triangles'], ',')))

    # Picking maps a triangle index to a grain through these ranges. If they do not
    # tile the merged mesh exactly the panel shows one grain's numbers for another --
    # which looks perfectly plausible on screen, so nothing but this catches it.
    tri, gaps = 0, []
    for rec in cmeta['grains']:
        if rec['tri0'] != tri:
            gaps.append((rec['id'], rec['tri0'], tri))
        tri += rec['ntri']
    chk('grain triangle ranges tile the mesh with no gap or overlap',
        not gaps and tri == n3, '%d triangles in %d grains, %d gaps'
        % (tri, len(cmeta['grains']), len(gaps)))
    chk('every placed grit is clickable',
        len(cmeta['grains']) == real['n_grits'],
        '%d grains vs %d in the deck' % (len(cmeta['grains']), real['n_grits']))
    # The viewer draws nearest-the-centre first and caps the count, so the grain
    # records are not in placement order -- match them by id, not by index.
    by_id = {pl.placement_id: i for i, pl in enumerate(plan['_model'].placements)}
    chk('grain ids are the deck placement ids',
        {r['id'] for r in cmeta['grains']} == set(by_id),
        '%d ids, first drawn %s' % (len(cmeta['grains']),
                                    [r['id'] for r in cmeta['grains']][:5]))

    # The protrusion the panel shows must be the protrusion the deck geometry has.
    baked = plan['_place']['baked']
    worst = max(abs(r['protrusion_um']
                    - (float(np.hypot(baked[by_id[r['id']]][:, 0],
                                      baked[by_id[r['id']]][:, 1]).max())
                       - plan['outer_radius_mm']) * 1000.0)
                for r in cmeta['grains'])
    chk('the protrusion shown on click is measured off the deck geometry',
        worst < 1e-3, 'worst disagreement %.2e um' % worst)
    # ...and the grain the viewer names must be the grain whose triangles it lit up.
    tri0 = {r['id']: r['tri0'] for r in cmeta['grains']}
    off, mism = 0, []
    for r in cmeta['grains']:
        v = baked[by_id[r['id']]]
        if abs(tri0[r['id']] - off) > 0:
            mism.append(r['id'])
        off += r['ntri']
        if len(plan['_place']['faces'][by_id[r['id']]]) != r['ntri']:
            mism.append(r['id'])
    chk('each grain owns exactly its own triangles in the merged mesh', not mism,
        'mismatched %s' % mism[:5])

    fd = np.array(cmeta['face_dir'])
    er = plan['_place']['e_r']
    chk('the Face view looks along the contact normal',
        abs(np.linalg.norm(fd) - 1) < 1e-9
        and abs(fd[0] - er[0]) < 1e-12 and abs(fd[1] - er[2]) < 1e-12
        and abs(fd[2] + er[1]) < 1e-12, str([round(float(x), 5) for x in fd]))

    chk('the model is embedded, not fetched from anywhere',
        'data:model/gltf-binary;base64,' in cad_html
        and 'developer.api.autodesk.com' not in cad_html
        and cad_html.count('http') == cad_html.count('https://unpkg.com/three@'))
    for feature, needle in (('section plane', 'localClippingEnabled'),
                            ('feature edges', 'EdgesGeometry'),
                            ('parts tree', 'cadtree'),
                            ('click to inspect', 'intersectObjects'),
                            ('measurement', 'distanceTo'),
                            ('orthographic', 'OrthographicCamera'),
                            ('orientation triad', 'ArrowHelper'),
                            ('save png', 'toDataURL')):
        chk('the viewer really implements %s' % feature, needle in cad_html)
    try:
        cadviewer.viewer_html(cvp, cmeta, max_inline_mb=1e-6)
        chk('an oversize model is refused rather than inlined', False)
    except ValueError as exc:
        chk('an oversize model is refused rather than inlined',
            'too large to inline' in str(exc))

    # EdgesGeometry reads geometry.attributes.position.count in its constructor, so
    # handing it an empty BufferGeometry throws -- and GLTFLoader routes exceptions out
    # of onLoad into onError, so the whole model then reports as unloadable. That is
    # exactly what the edge budget did on the first deck whose grit mesh crossed the
    # threshold: 8,434 grains, 966k triangles, viewer dead.
    vsrc = cadviewer._TEMPLATE
    chk('EdgesGeometry is never built from an empty geometry',
        not re.search(r'EdgesGeometry\([^)]*new THREE\.BufferGeometry\(\)', vsrc),
        'the heavy path makes a bare LineSegments instead')
    chk('a heavy grit mesh still gets an edge geometry on demand',
        'userData.heavy' in vsrc and 'new THREE.EdgesGeometry(p.mesh.geometry' in vsrc)
    # The same size sweep the viewer will actually meet: build a model big enough to
    # trip the edge budget and confirm it still assembles.
    _over = dict(sector_mode='angle', sector_deg=4.0, grit_mode='count',
                 grit_count=1500, grit_width_window_mm=0.05)
    big = DeckParams(name='big', analysis=AnalysisParams(**AN),
                     **dict({k: v for k, v in ref_kw.items()
                             if k not in ('arc_length_mm',)}, **_over))
    bp = plan_deck(big, solids)
    bparts, bmeta = parts_from_plan(bp, mode='whole wheel', with_meta=True,
                                    budget_mb=24.0)
    btri = sum(len(np.asarray(x['faces'])) for x in bparts
               if 'grit' in x['name'])
    # Assert grain conservation, not a triangle threshold: how many grains a given
    # sector accepts depends on packing, so a hard-coded count made the check about the
    # fixture rather than the viewer. Whether the edge budget fires is settled by the
    # static guard above.
    chk('a dense grit mesh loses no grain to the budget',
        len(bmeta['grains']) + len(bmeta['grains_far']) == bp['n_grits'],
        '%s grit triangles, %d full + %d boxes = %d of %d grains'
        % (format(btri, ','), len(bmeta['grains']), len(bmeta['grains_far']),
           len(bmeta['grains']) + len(bmeta['grains_far']), bp['n_grits']))

    # And, if a browser is on this machine, drive the thing for real.
    probe = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_cadv_probe.py')
    if os.path.exists(probe):
        # Launching a headless browser while a previous one is still shutting down
        # sometimes fails outright, which says nothing about the viewer. One retry.
        # The probe drives three view modes, each a full browser session with
        # thousands of synthesised clicks, so it needs a generous ceiling -- and a
        # timeout has to be reported as "did not finish", never as a viewer fault.
        tail = []
        for attempt in (1, 2):
            try:
                r = subprocess.run([sys.executable, probe], capture_output=True,
                                   text=True, timeout=1800)
            except subprocess.TimeoutExpired:
                tail = []
                r = None
                break
            tail = [l for l in r.stdout.splitlines() if 'BROWSER PROBE' in l]
            if tail:
                break
            time.sleep(5)
        if r is None:
            warn('a real browser drives the viewer end to end', False,
                 'the probe did not finish in 30 min on this machine')
        else:
            warn('a real browser drives the viewer end to end', r.returncode == 0,
                 tail[-1] if tail else (r.stdout or r.stderr)[-160:].replace('\n', ' '))
    else:
        warn('a real browser drives the viewer end to end', False, 'probe not present')

    # ---------------- C6e workpiece placement and standoff ----------------
    print()
    print('C6e PLACEMENT + STANDOFF  (the block must go where it was asked, and the')
    print('    gap must be the gap that was asked for)')
    from semgrit.rigid_wheel import WP_POSITIONS

    pk = dict(ref_kw)
    pk.pop('clearance_um', None)
    seen_theta = {}
    for pos in WP_POSITIONS:
        kw = dict(pk, wp_position=pos, wp_position_deg=1.0)
        pl = plan_deck(DeckParams(name='p_' + pos.replace(' ', '_'),
                                  analysis=AnalysisParams(**AN), **kw), solids)
        bd = build_deck(DeckParams(name='b_' + pos.replace(' ', '_'),
                                   analysis=AnalysisParams(**AN), **kw), solids, OUT)
        seen_theta[pos] = pl['theta_workpiece_deg']
        chk("'%s' puts the block where the preview said" % pos,
            abs(pl['theta_workpiece_deg'] - bd['theta_workpiece_deg']) < 1e-12
            and pl['n_grits_under_block'] >= 1,
            'theta %.5f deg, %d grain(s) under it'
            % (pl['theta_workpiece_deg'], pl['n_grits_under_block']))
    # 'centred' and 'under the tallest grit' coincide whenever the centred footprint
    # holds no grit, because the fallback lands on exactly that grain. That is correct
    # behaviour, so only require that the modes are not all the same place.
    chk('the positions are actually different places',
        len({round(v, 6) for v in seen_theta.values()}) >= 3,
        str({k: round(v, 4) for k, v in seen_theta.items()}))
    chk("'custom angle' lands on the angle asked for",
        abs(seen_theta['custom angle'] - 1.0) < 1e-9,
        '%.9f deg' % seen_theta['custom angle'])

    # The entry edge is the high-theta side, because the surface travels toward
    # decreasing theta. 'first grit at entry' must put it at the leading grain.
    pe = plan_deck(DeckParams(name='p_entry', wp_position='first grit at entry',
                              analysis=AnalysisParams(**AN), **pk), solids)
    baked_e = pe['_place']['baked']
    hz_e = pe['workpiece']['width_mm'] / 2.0
    lead = max(float(np.arctan2(v[:, 1], v[:, 0]).max()) for v in baked_e
               if float(np.abs(v[:, 2]).min()) <= hz_e)
    chk('the entry edge sits at the leading grain',
        abs(pe['wp_entry_theta_deg'] - np.degrees(lead)) < 1e-9,
        'entry %.9f deg vs leading grain %.9f deg'
        % (pe['wp_entry_theta_deg'], np.degrees(lead)))
    # Only grains the block can reach across the face matter: one sitting outside the
    # block's width may well lie at a higher angle, and should.
    chk('no reachable grain is dressed past the entry edge',
        pe['wp_entry_theta_deg'] >= pe['grit_theta_reachable_deg'][1] - 1e-9,
        'entry %.6f vs furthest reachable grain %.6f'
        % (pe['wp_entry_theta_deg'], pe['grit_theta_reachable_deg'][1]))

    pt = plan_deck(DeckParams(name='p_tall', wp_position='under the tallest grit',
                              analysis=AnalysisParams(**AN), **pk), solids)
    bt = pt['_place']['baked']
    tall = max((i for i, v in enumerate(bt)
                if float(np.abs(v[:, 2]).min()) <= pt['workpiece']['width_mm'] / 2.0),
               key=lambda i: float(np.hypot(bt[i][:, 0], bt[i][:, 1]).max()))
    ct = bt[tall].mean(axis=0)
    chk("'under the tallest grit' centres on the tallest reachable grain",
        abs(pt['theta_workpiece_deg']
            - np.degrees(np.arctan2(ct[1], ct[0]))) < 1e-9,
        'grain %d' % pt['_model'].placements[tall].placement_id)

    # A standoff must lift the ground face by exactly itself, and nothing else.
    base = plan_deck(DeckParams(name='p_s0', clearance_um=0.0,
                                analysis=AnalysisParams(**AN), **pk), solids)
    for s in (0.25, 1.0, 3.0):
        ps = plan_deck(DeckParams(name='p_s%g' % s, clearance_um=s,
                                  analysis=AnalysisParams(**AN), **pk), solids)
        lift = (ps['ground_radius_mm'] - base['ground_radius_mm']) * 1000.0
        chk('a %.2f um standoff lifts the face by exactly that' % s,
            abs(lift - s) < 1e-9 and ps['n_grits'] == base['n_grits'],
            'lifted %.9f um' % lift)
        chk('and moves the whole depth window with it (%.2f um)' % s,
            abs((ps['first_contact_um'] - base['first_contact_um']) - s) < 1e-9
            and abs((ps['depth_ceiling_um'] - base['depth_ceiling_um']) - s) < 1e-9,
            'window %.3f..%.3f um' % (ps['first_contact_um'], ps['depth_ceiling_um']))

    # The heights the notebook reports have to be the heights of the geometry.
    st = base['protrusion_um']
    tip = build_deck(DeckParams(name='b_tip', clearance_um=0.0,
                                analysis=AnalysisParams(**AN), **pk),
                     solids, OUT)['tallest_tip_whole_arc_mm']
    chk('the tallest protrusion reported is the tallest grain in the deck',
        abs(st['max'] - (tip - base['outer_radius_mm']) * 1000.0) < 1e-9,
        'max %.4f um over %d grains' % (st['max'], st['n']))
    chk('the governing grain cannot out-stand the tallest grain',
        base['bond_clearance_um'] <= st['max'] + 1e-9,
        'governing %.4f um vs tallest %.4f um' % (base['bond_clearance_um'], st['max']))
    chk('grains under the block are a subset of all grains',
        base['protrusion_under_block_um']['n'] <= st['n']
        and base['protrusion_under_block_um']['max'] <= st['max'] + 1e-9)

    # And the preview's verdict on the depth of cut must be the build's verdict.
    from semgrit.build_deck import DeckError
    from semgrit.preview import summary_text as _sumtxt
    disagreed = []
    for s, ae in ((0.0, 0.0), (0.5, 0.0), (2.0, 1.0), (0.0, 99.0), (0.5, 1.2)):
        a2 = dict(AN)
        a2['depth_of_cut_um'] = ae
        kw = dict(pk, clearance_um=s)
        said = '***' in _sumtxt(plan_deck(
            DeckParams(name='v', analysis=AnalysisParams(**a2), **kw), solids))
        try:
            build_deck(DeckParams(name='v', analysis=AnalysisParams(**a2), **kw),
                       solids, OUT)
            refused = False
        except DeckError:
            refused = True
        if said != refused:
            disagreed.append((s, ae, said, refused))
    chk('the preview refuses exactly the depths the build refuses', not disagreed,
        str(disagreed))

    # ---------------- C6f whole wheel and every grit ----------------
    print()
    print('C6f WHOLE WHEEL + ALL GRITS  (context must be labelled, and nothing may')
    print('    be dropped without saying so)')
    from semgrit.glb import FAR_PART, GHOST_PART, MARK_PART, MODES

    names = {}
    for md in MODES:
        prts, mm = parts_from_plan(plan, mode=md, with_meta=True)
        names[md] = [x['name'] for x in prts]
        drawn = len(mm['grains']) + len(mm['grains_far'])
        chk("mode '%s' draws every grain in view" % md,
            drawn == mm['grits_in_view'] and drawn == mm['grits_drawn'],
            '%d drawn of %d in view, %d on the wheel'
            % (drawn, mm['grits_in_view'], mm['grits_total']))
    chk('wheel mode draws every grain on the wheel',
        len(parts_from_plan(plan, mode='wheel', with_meta=True)[1]['grains'])
        == real['n_grits'], '%d grains' % real['n_grits'])
    chk('only whole-wheel mode adds context parts',
        GHOST_PART not in names['wheel'] and GHOST_PART not in names['contact']
        and GHOST_PART in names['whole wheel'] and MARK_PART in names['whole wheel'],
        str(names['whole wheel']))

    # The context parts are the one thing drawn that is not in the deck. That is
    # allowed only because they are separately named and excluded here; if either
    # ever merged into a real part this check is what would notice.
    wparts, wmeta = parts_from_plan(plan, mode='whole wheel', with_meta=True)
    deck_named = [x for x in wparts if x['name'] not in (GHOST_PART, MARK_PART)]
    ctx = [x for x in wparts if x['name'] in (GHOST_PART, MARK_PART)]
    chk('context is separable from the deck by name alone',
        len(ctx) == 2 and len(deck_named) == 3, str([x['name'] for x in ctx]))
    gp2 = os.path.join(OUT, 'ww.glb')
    write_glb(gp2, deck_named)
    sc2 = trimesh.load(gp2)
    GV2 = np.column_stack([sc2.geometry['abrasive grits'].vertices[:, 0],
                           -sc2.geometry['abrasive grits'].vertices[:, 2],
                           sc2.geometry['abrasive grits'].vertices[:, 1]])
    gd2, _ = cKDTree(deck_nodes).query(GV2)
    chk('with the context removed, what is left is still exactly the deck',
        float(gd2.max()) < 5e-6, 'worst offset %.3e mm' % gd2.max())
    chk('the ghost is a wheel, not the rim ring',
        abs(max(abs(ctx[0]['vertices'][:, 0]).max(),
                abs(ctx[0]['vertices'][:, 1]).max())
            - plan['outer_radius_mm']) < 1e-9
        and float(np.hypot(ctx[0]['vertices'][:, 0],
                           ctx[0]['vertices'][:, 1]).min()) < 1e-9,
        'radius 0 to %.3f mm' % plan['outer_radius_mm'])

    # Degrade, never drop: squeezing the budget must convert grains to boxes and keep
    # the count whole.
    for mb in (24.0, 0.9, 0.3):
        _p, _m = parts_from_plan(plan, mode='wheel', with_meta=True, budget_mb=mb)
        tot = len(_m['grains']) + len(_m['grains_far'])
        _i = write_glb(os.path.join(OUT, 'b.glb'), _p)
        chk('a %.1f MB budget simplifies rather than drops' % mb,
            tot == real['n_grits'] and _i['bytes'] * 4 / 3 <= mb * 1e6,
            '%d full + %d boxes = %d grains, %.2f MB encoded'
            % (len(_m['grains']), len(_m['grains_far']), tot,
               _i['bytes'] * 4 / 3 / 1e6))
        if _m['grains_far']:
            chk('and says so', bool(_m['notes']), _m['notes'][0][:70])
            chk('a simplified grain still reports its real size',
                all(g['proxy'] for g in _m['grains_far'])
                and all(not g['proxy'] for g in _m['grains']))

    # ---------------- C8 simple mode ----------------
    print()
    print('C8  SIMPLE MODE  (must be the advanced path with defaults, not a second one)')
    from semgrit.quick import (GRIT_KINDS, WORKPIECE_SIZES, QuickError,
                               simple_params)

    for kind, val in (('a fixed number', 60), ('single grain', 0),
                      ('concentration', 100.0), ('grains per mm2', 3000.0)):
        sp = simple_params(diameter_mm=50.0, slice_mm=1.0, grit_kind=kind,
                           grit_value=val,
                           workpiece_mm=WORKPIECE_SIZES['small  48 x 15 x 6 um'],
                           wp_position='centred', standoff_um=0.0, run_ready=True,
                           cae_deck=False, cad=False, name='q_' + kind.split()[0])
        got = build_deck(sp, solids, OUT)
        vok = True
        for v in ('verify_rigid_deck.py', 'verify_rigid_deck2.py'):
            r = subprocess.run([sys.executable, v, got['path']],
                               capture_output=True, text=True)
            vok = vok and r.returncode == 0
        chk("simple mode '%s' builds a deck both verifiers accept" % kind, vok,
            '%d grits, %.1f MB' % (got['n_grits'], got['size_bytes'] / 1e6))
    chk('every grit kind the widget offers is handled',
        sorted(GRIT_KINDS) == sorted(
            ['a fixed number', 'single grain', 'concentration', 'grains per mm2']))

    # The anti-divergence check: the same intent through both paths, same bytes out.
    sp = simple_params(diameter_mm=50.0, slice_mm=1.0, grit_kind='a fixed number',
                       grit_value=60,
                       workpiece_mm=WORKPIECE_SIZES['small  48 x 15 x 6 um'],
                       wp_position='centred', standoff_um=0.0, run_ready=True,
                       cae_deck=False, cad=False, name='xsame')
    adv = dataclasses.replace(sp, name='xsame')
    a1 = build_deck(sp, solids, OUT)['path']
    a2 = build_deck(adv, solids, OUT)['path']
    chk('simple and advanced produce the same deck for the same intent',
        [l for l in open(a1, encoding='ascii') if not l.startswith('**')]
        == [l for l in open(a2, encoding='ascii') if not l.startswith('**')])
    for bad, why in ((dict(slice_mm=0.02), 'a slice shorter than the block'),
                     (dict(grit_kind='lots'), 'an unknown grit kind'),
                     (dict(wp_position='sideways'), 'an unknown position'),
                     (dict(workpiece_mm=(0.048, 0.0, 0.006)), 'a zero dimension')):
        kw = dict(diameter_mm=50.0, slice_mm=1.0, grit_kind='a fixed number',
                  grit_value=60, workpiece_mm=(0.048, 0.015, 0.006),
                  wp_position='centred', standoff_um=0.0, run_ready=True,
                  cae_deck=False, cad=False)
        kw.update(bad)
        try:
            simple_params(**kw)
            chk('simple mode refuses %s' % why, False)
        except QuickError:
            chk('simple mode refuses %s' % why, True)

    # ---------------- C9 grinding theory ----------------
    print()
    print('C9  GRINDING THEORY  (the measured column must be the deck, not an estimate)')
    from semgrit.grinding_theory import format_report
    from semgrit.grinding_theory import report as theory_report

    th = theory_report(plan, work_speed_mm_s=200.0)
    m = th['measured']
    chk('grains placed is the deck grit count', m['grains_placed'] == real['n_grits'],
        '%d' % m['grains_placed'])
    chk('grains under the block matches the plan',
        m['grains_under_the_block'] == plan['n_grits_under_block'])
    sw = np.asarray(plan['swept_clearances_um'])
    chk('active grains are the ones the infeed reaches',
        m['active_grains'] == int((sw <= plan['depth_of_cut_um']).sum()),
        '%d of %d in the swept band' % (m['active_grains'], sw.size))
    chk('static grain density matches the dressed band',
        abs(m['static_grain_density_per_mm2']
            - real['n_grits'] / (plan['grit_band_arc_mm']
                                 * plan['grit_band_width_mm'])) < 1e-6,
        '%.0f /mm2' % m['static_grain_density_per_mm2'])
    chk('wheel speed and rpm agree with the deck kinematics',
        abs(m['wheel_speed_m_s'] * 1000
            - plan['cost']['surface_speed_mm_s']) < 1e-9)
    chk('mesh resolution is read off the written elements',
        abs(th['mesh']['element_depth_um'] - plan['element_um'][2]) < 1e-12)
    # Without a work speed the classical rows must be withheld, not zeroed.
    t0 = theory_report(plan, work_speed_mm_s=0.0)
    chk('with no work speed the classical rows are withheld, not faked',
        not t0['theory']['applicable']
        and 'max_chip_thickness_um' not in t0['theory']
        and any('plunge' in n for n in t0['notes']))
    chk('the report renders', 'GRAINS, counted off the deck' in format_report(th))

    # ---------------- C10 the odb post-processor ----------------
    print()
    print('C10 POST-PROCESSING  (it must ask only for output the deck writes)')
    from semgrit.odbpost import history_outputs_required

    post = build_deck(DeckParams(name='pp', analysis=AnalysisParams(**AN), **ref_kw),
                      solids, OUT)
    script = post.get('postprocess_script')
    chk('a run-ready deck ships its post-processor',
        bool(script) and os.path.exists(script),
        os.path.basename(script) if script else 'missing')
    src = open(script, encoding='ascii').read()
    try:
        ast.parse(src)
        chk('the script parses', True, '%d lines' % src.count('\n'))
    except SyntaxError as exc:
        chk('the script parses', False, str(exc)[:70])
    py3only = []
    if re.search(r'\bf"', src) or re.search(r"\bf'", src):
        py3only.append('f-string')
    if ':=' in src:
        py3only.append('walrus')
    if re.search(r'\)\s*->\s*\w', src):
        py3only.append('return annotation')
    # Abaqus 2023 and earlier run Python 2.7, so anything 3-only makes the script
    # unusable on exactly the versions most likely to be installed.
    chk('nothing in it is Python-3 only', not py3only, str(py3only))
    chk('it does not need a CAE licence', 'odbAccess' in src
        and 'from abaqus import' not in src and 'session.' not in src)

    deck_txt = open(post['path'], encoding='ascii').read()
    nb = re.search(r'\*Node Output, nset=A_WHEEL_REF\n([^\n]*)', deck_txt)
    asked = set(x.strip() for x in nb.group(1).split(',')) if nb else set()
    need_rf = [n for n in history_outputs_required() if n[:2] in ('RF', 'RM')]
    chk('the deck requests the forces the script reads',
        all(n in asked for n in need_rf), 'asked for %s' % sorted(asked))
    chk('and the node set it names exists',
        bool(re.search(r'^\*Nset, nset=A_WHEEL_REF', deck_txt, re.M)))
    chk('the deck requests the energies the script reads',
        '*Output, history, variable=PRESELECT' in deck_txt)
    # grind11b died here: NUMBER INTERVAL is a field-output parameter and Abaqus
    # refuses it on history output, so the job never reached the solver.
    _hti = re.search(r'\*Output,\s*history,\s*time interval=\S+', deck_txt, re.I)
    chk('history output uses TIME INTERVAL, not NUMBER INTERVAL',
        not re.search(r'\*Output,\s*history[^\n]*number\s*interval', deck_txt, re.I)
        and _hti is not None,
        _hti.group(0) if _hti else 'absent')
    _eo = re.search(r'\*Element Output[^\n]*', deck_txt)
    chk('element output is scoped to the workpiece, not the rigid wheel',
        _eo is not None and 'elset=A_WP_ALL' in _eo.group(0),
        _eo.group(0) if _eo else 'no *Element Output')
    chk('and that elset is defined in the assembly',
        bool(re.search(r'^\*Elset, elset=A_WP_ALL, instance=WP-1$', deck_txt, re.M)))
    chk('the script reports rather than dies when an output is absent',
        'is None' in src and 'return None' in src)

    # ---------------- C12 boundary conditions ----------------
    print()
    print('C12 BOUNDARY CONDITIONS  (draw only what the deck writes, and draw it the')
    print('    way the deck writes it)')
    from semgrit.bcspec import build as bc_build
    from semgrit.bcspec import to_viewer

    bcp = DeckParams(name='bc', analysis=AnalysisParams(**AN), **ref_kw)
    bplan = plan_deck(bcp, solids)
    bdeck = build_deck(bcp, solids, OUT)
    btxt = open(bdeck['path'], encoding='ascii').read()
    spec = bc_build(bplan)
    chk('a run-ready deck reports boundary conditions', spec['has_analysis']
        and len(spec['items']) >= 4, '%d items' % len(spec['items']))

    # Nothing may be drawn that the deck does not contain.
    ghost_kw = [it for it in spec['items'] if it['keyword'].split(',')[0] not in btxt]
    chk('every glyph stands for a keyword the deck really writes', not ghost_kw,
        str([(i['kind'], i['keyword']) for i in ghost_kw]))
    unknown = []
    for it in spec['items']:
        nm = it.get('set') or ''
        if nm in ('ALL EXTERIOR', ''):
            continue
        if not re.search(r'^\*(Nset|Elset|Surface)[^\n]*(nset|elset|name)=%s\b'
                         % re.escape(nm), btxt, re.M | re.I):
            unknown.append(nm)
    chk('every set a glyph names is defined in the deck', not unknown, str(unknown))

    # The held faces must be exactly the ones the *Boundary block lists.
    mb = re.search(r'^\*Boundary\s*\n((?:\w+, ENCASTRE\s*\n)+)', btxt, re.M)
    deck_held = set(re.findall(r'(\w+), ENCASTRE', mb.group(1))) if mb else set()
    drawn_held = {it['set'] for it in spec['items'] if it['kind'] == 'encastre'}
    chk('the held faces drawn are exactly the held faces in the deck',
        deck_held == drawn_held, 'deck %s vs drawn %s'
        % (sorted(deck_held), sorted(drawn_held)))

    # The rotation glyph's sense must match the sign of VR3 the deck applies.
    vr3 = float(re.search(r'A_WHEEL_REF, 6, 6, (-?[\d.eE+]+)', btxt).group(1))
    rot = next(it for it in spec['items'] if it['kind'] == 'rotation')
    chk('the rotation glyph carries the sign the BC actually applies',
        (rot['sign'] < 0) == (vr3 < 0)
        and abs(abs(vr3) - rot['magnitude']) < 1e-6,
        'VR3 = %+g rad/s, glyph sign %+g, magnitude %g' % (vr3, rot['sign'],
                                                           rot['magnitude']))
    chk('and says which way the surface travels',
        'decreasing theta' in rot['detail'] if vr3 < 0
        else 'increasing theta' in rot['detail'], rot['detail'][:60])

    # The infeed arrow must point the way the velocity BC feeds -- the bug that made
    # every deck run without cutting, now checked in the picture as well as the deck.
    v1 = float(re.search(r'A_WHEEL_REF, 1, 1, (-?[\d.eE+]+)', btxt).group(1))
    v2 = float(re.search(r'A_WHEEL_REF, 2, 2, (-?[\d.eE+]+)', btxt).group(1))
    vel = next((it for it in spec['items'] if it['kind'] == 'velocity'), None)
    chk('an infeed arrow exists when the deck feeds', vel is not None)
    if vel:
        deck_dir = to_viewer((v1 / np.hypot(v1, v2), v2 / np.hypot(v1, v2), 0.0))
        dot = float(np.dot(deck_dir, vel['dir']))
        chk('the infeed arrow points the way the BC actually feeds', dot > 0.999,
            'arrow . BC = %+.6f' % dot)
        chk('and its magnitude is the commanded speed',
            abs(vel['magnitude'] - np.hypot(v1, v2)) < 1e-6,
            '%.3f vs %.3f mm/s' % (vel['magnitude'], np.hypot(v1, v2)))

    # The contact highlight must be the deck's own engaging set, not a lookalike.
    cnt = next((it for it in spec['items'] if it['kind'] == 'contact'), None)
    if cnt and cnt['set'] == 'A_GRITS_ENGAGE_SURF':
        eng_blk = re.search(r'^\*Elset, elset=ES_GRITS_ENGAGE[^\n]*\n((?:[^*]+\n)+)',
                            btxt, re.M)
        n_deck = len(re.findall(r'\d+', eng_blk.group(1))) if eng_blk else -1
        # ES_GRITS_ENGAGE lists facets; the meta flags grains. Compare grain counts via
        # the plan's own engage list, which is the shared helper the writer calls.
        chk('the contact glyph counts the deck engaging set',
            cnt['n_engaging'] == len(bplan['_engage']),
            '%d grains, %d facets in ES_GRITS_ENGAGE'
            % (cnt['n_engaging'], n_deck))
    _, bmeta = parts_from_plan(bplan, mode='wheel', with_meta=True)
    chk('the grains flagged engaging are the plan engaging set',
        sum(1 for g in bmeta['grains'] if g['engage']) == len(bplan['_engage']),
        '%d flagged of %d drawn' % (sum(1 for g in bmeta['grains'] if g['engage']),
                                    len(bmeta['grains'])))

    # A geometry-only deck has no BCs, and must say so rather than invent any.
    gspec = bc_build(plan_deck(DeckParams(name='g', **ref_kw), solids))
    chk('a geometry-only deck draws no boundary conditions',
        not gspec['has_analysis'] and not gspec['items'] and gspec['notes'],
        gspec['notes'][0][:64] if gspec['notes'] else 'no note')

    # The viewer must not have gained a glTF part: glyphs travel as numbers.
    chk('glyphs did not become geometry in the glTF',
        len([x for x in parts_from_plan(bplan, mode='wheel')]) == 3,
        '%d parts' % len(parts_from_plan(bplan, mode='wheel')))

    # ---------------- C13 editing ----------------
    print()
    print('C13 EDITING  (a number typed in the browser must reach the deck, by one path)')
    from semgrit.editable import FIELDS, LIVE, REBUILD
    from semgrit.editable import apply as edit_apply
    from semgrit.editable import load as edit_load
    from semgrit.editable import params_from_settings, settings_from_params

    ebase = DeckParams(name='ed', analysis=AnalysisParams(**AN), **ref_kw)
    s0 = settings_from_params(ebase)
    # Re-applying the settings unchanged must be a no-op. It was not: every dict carries
    # sector_deg, and applying it switched sector_mode from 'arc' to 'angle'.
    p_same = params_from_settings(s0, ebase)
    chk('re-applying unchanged settings changes nothing',
        settings_from_params(p_same) == s0
        and p_same.sector_mode == ebase.sector_mode,
        'sector_mode %s -> %s' % (ebase.sector_mode, p_same.sector_mode))

    # The edit the browser was driven through, verbatim from the probe run.
    browser = {'clearance_um': 0.45, 'wp_length_mm': 0.030, 'wp_position_deg': 0.9,
               'wp_position': 'custom angle', 'rotation_reversed': True,
               'grit_count': 25}
    got = edit_apply(dict(s0, **browser), ebase, solids)
    chk('every edited field is reported as changed',
        set(got['changed']) == set(browser), str(sorted(got['changed'])))
    chk('a grit-count change is classified as needing a rebuild',
        got['tier'] == 'rebuild', got['tier'])
    edeck = build_deck(got['params'], solids, OUT)
    etxt = open(edeck['path'], encoding='ascii').read()
    chk('the deck carries the edited standoff',
        abs(edeck['clearance_um'] - browser['clearance_um']) < 1e-12,
        '%.4f um' % edeck['clearance_um'])
    chk('the deck carries the edited grit count',
        edeck['n_grits'] <= browser['grit_count'],
        '%d grits asked, %d placed' % (browser['grit_count'], edeck['n_grits']))
    chk('the deck carries the edited block length',
        abs(edeck['params']['wp_length_mm'] - browser['wp_length_mm']) < 1e-12)
    evr3 = float(re.search(r'A_WHEEL_REF, 6, 6, (-?[\d.eE+]+)', etxt).group(1))
    chk('reversing the rotation in the browser flips VR3 in the deck', evr3 > 0,
        'VR3 = %+g rad/s' % evr3)
    chk('and the deck header follows the reversed sense',
        'increasing theta' in etxt and 'decreasing theta' not in etxt)
    chk('and the entry edge moves to the other end',
        got['plan']['wp_entry_theta_deg'] < got['plan']['wp_exit_theta_deg'],
        'entry %.4f vs exit %.4f deg' % (got['plan']['wp_entry_theta_deg'],
                                        got['plan']['wp_exit_theta_deg']))
    for v in ('verify_rigid_deck.py', 'verify_rigid_deck2.py'):
        r = subprocess.run([sys.executable, v, edeck['path']], capture_output=True,
                           text=True)
        chk('an edited deck still passes %s' % v, r.returncode == 0,
            str([l.strip() for l in r.stdout.splitlines() if '[FAIL]' in l][:1]))

    # The two commit paths must not diverge: what the file path produces and what the
    # kernel path produces have to be the same deck.
    sp = os.path.join(OUT, 'viewer_settings.json')
    with open(sp, 'w') as fh:
        json.dump({'settings': dict(s0, **browser)}, fh)
    g2 = edit_apply(edit_load(sp), ebase, solids)
    chk('the file path and the direct path give the same settings',
        g2['settings'] == got['settings'])
    d2 = build_deck(dataclasses.replace(g2['params'], name='ed2'), solids, OUT)
    a = [l for l in open(edeck['path'], encoding='ascii') if not l.startswith('**')]
    b = [l for l in open(d2['path'], encoding='ascii') if not l.startswith('**')]
    chk('...and therefore the same deck, line for line', a == b, '%d lines' % len(a))

    # Edits that would make an invalid model must be refused, with a reason.
    for bad, why in (({'wp_length_mm': 5.0}, 'a block longer than the slice'),
                     ({'wp_width_mm': 9.0}, 'a block wider than the face'),
                     ({'clearance_um': -1.0}, 'a negative standoff'),
                     ({'nope': 1.0}, 'a field that is not editable'),
                     ({'wp_position': 'sideways'}, 'an unknown placement')):
        try:
            params_from_settings(dict(s0, **bad), ebase)
            chk('editing refuses %s' % why, False)
        except DeckError:
            chk('editing refuses %s' % why, True)

    # Every field the panel offers must be one the deck actually reacts to.
    dead = []
    for f in FIELDS:
        if f.key == '_sector_mode':
            continue
        v = s0[f.key]
        nv = (not v) if isinstance(v, bool) else (
            f.choices[(list(f.choices).index(v) + 1) % len(f.choices)]
            if f.choices else (v * 1.5 + 1.0))
        try:
            q = edit_apply(dict(s0, **{f.key: nv}), ebase, solids)
            if f.key not in q['changed']:
                dead.append(f.key)
        except DeckError:
            pass          # refused is fine; silently ignored is not
    chk('no editable field is inert', not dead, str(dead))
    chk('the tiers cover every field', len(LIVE) + len(REBUILD) == len(FIELDS))

    # ---------------- C14 guardrails and dragging ----------------
    print()
    print('C14 GUARDRAILS  (the viewer must know the window it is editing inside, and')
    print('    must not misdescribe the wheel)')
    from semgrit.editable import FIELDS as EFIELDS
    from semgrit.editable import param_block

    gp = DeckParams(name='gd', analysis=AnalysisParams(**dict(AN, depth_of_cut_um=1.2)),
                    **ref_kw)
    gplan = plan_deck(gp, solids)
    _, gmeta = parts_from_plan(gplan, mode='whole wheel', with_meta=True)
    ge = gmeta['edit']
    chk('the viewer carries the deck depth-of-cut window',
        abs(ge['first_contact_um'] - gplan['first_contact_um']) < 1e-12
        and abs(ge['depth_ceiling_um'] - gplan['depth_ceiling_um']) < 1e-12,
        '%.4f .. %.4f um' % (ge['first_contact_um'], ge['depth_ceiling_um']))
    chk('and the deck engaging count',
        ge['engaging_now'] == len(gplan['_engage']), '%d grits' % ge['engaging_now'])

    # A depth inside the band must build; outside it must be refused, and the message
    # the viewer shows has to be the message Python gives.
    inside = 0.5 * (max(ge['first_contact_um'], 0.0) + ge['depth_ceiling_um'])
    ok_deck = build_deck(DeckParams(
        name='gd_in', analysis=AnalysisParams(**dict(AN, depth_of_cut_um=inside)),
        **ref_kw), solids, OUT)
    chk('a depth inside the band builds', bool(ok_deck['path']),
        '%.4f um' % inside)
    for ae, why in ((ge['depth_ceiling_um'] * 1.5, 'above the face-to-bond gap'),):
        try:
            build_deck(DeckParams(
                name='gd_out', analysis=AnalysisParams(**dict(AN, depth_of_cut_um=ae)),
                **ref_kw), solids, OUT)
            chk('a depth %s is refused' % why, False)
        except DeckError as exc:
            chk('a depth %s is refused' % why, True, str(exc)[:58])

    # The standoff shift the viewer applies must be the shift Python computes.
    delta = 0.6
    shifted = plan_deck(DeckParams(
        name='gd_s', clearance_um=delta,
        analysis=AnalysisParams(**dict(AN, depth_of_cut_um=1.2)),
        **{k: v for k, v in ref_kw.items() if k != 'clearance_um'}), solids)
    chk('shifting the standoff shifts both ends of the window by exactly that',
        abs((shifted['first_contact_um'] - gplan['first_contact_um']) - delta) < 1e-9
        and abs((shifted['depth_ceiling_um'] - gplan['depth_ceiling_um'])
                - delta) < 1e-9,
        'lo %+.6f, hi %+.6f um'
        % (shifted['first_contact_um'] - gplan['first_contact_um'],
           shifted['depth_ceiling_um'] - gplan['depth_ceiling_um']))

    # The sector field used to read 30 deg on a 2.29 deg arc wheel.
    chk('the viewer reports the resolved wheel extent, not the raw field',
        abs(ge['sector_resolved_deg'] - gp.resolved_sector_deg()) < 1e-12
        and abs(ge['arc_length_mm'] - gplan['arc_length_mm']) < 1e-12,
        '%.4f mm arc = %.4f deg (raw sector_deg field is %g)'
        % (ge['arc_length_mm'], ge['sector_resolved_deg'], gp.sector_deg))
    chk('and says which of arc/angle is in force',
        ge['settings']['_sector_mode'] == gp.sector_mode)

    # Every field the panel offers must map to a widget that exists, and the paste block
    # must reproduce the same params -- this is what catches the m/s unit trap.
    nb_body = '\n'.join(cells)
    missing = [f.widget for f in EFIELDS
               if f.widget and not re.search(r'^%s\s*=.*#@param' % f.widget,
                                             nb_body, re.M)]
    chk('every editable field maps to a widget the notebook really has', not missing,
        str(missing))
    edited = dict(settings_from_params(gp), surface_speed_mm_s=20_000.0,
                  width_mm=0.05, clearance_um=0.45)
    blk = param_block(edited, gp)
    ns = {}
    exec(blk, {}, ns)
    chk('the paste block round-trips through the widget units',
        abs(ns['SURFACE_SPEED_M_S'] * 1000.0 - 20_000.0) < 1e-9
        and abs(ns['WHEEL_WIDTH_MM'] - 0.05) < 1e-12,
        'SURFACE_SPEED_M_S = %g m/s, WHEEL_WIDTH_MM = %g' % (ns['SURFACE_SPEED_M_S'],
                                                             ns['WHEEL_WIDTH_MM']))
    chk('the paste block names only changed fields',
        'DIAMETER_MM' not in blk and 'CLEARANCE_UM' in blk, blk.replace('\n', ' · '))

    # A drag is only ever a parameter delta, so dragging and typing must agree exactly.
    r_g = gplan['ground_radius_mm']
    drag_um = 25.0
    typed = math.degrees(drag_um / 1000.0 / r_g)
    dkw = {k: v for k, v in ref_kw.items()
           if k not in ('wp_position', 'wp_position_deg')}
    dplan = plan_deck(DeckParams(
        name='gd_drag', wp_position='custom angle', wp_position_deg=typed,
        analysis=AnalysisParams(**dict(AN, depth_of_cut_um=1.2)), **dkw), solids)
    # The requested angle must arrive intact. The *seated* angle may differ, because a
    # footprint with no grit under it is relocated to the nearest reachable grain -- and
    # when that happens the deck says so rather than silently honouring the request.
    chk('a drag arrives as the angle that was asked for',
        abs(dplan['wp_requested_theta_deg'] - typed) < 1e-9,
        'requested %.6f deg for a %.1f um drag at r = %.4f mm'
        % (typed, drag_um, r_g))
    chk('and any relocation away from it is reported, not hidden',
        (abs(dplan['theta_workpiece_deg'] - typed) < 1e-9)
        or bool(dplan['wp_relocated']),
        'seated %.6f deg, relocated %s'
        % (dplan['theta_workpiece_deg'], dplan['wp_relocated']))
    d1 = build_deck(DeckParams(
        name='gd_drag', wp_position='custom angle', wp_position_deg=typed,
        analysis=AnalysisParams(**dict(AN, depth_of_cut_um=1.2)), **dkw), solids, OUT)
    chk('and the dragged deck builds', bool(d1['path']),
        '%d grits, %.1f MB' % (d1['n_grits'], d1['size_bytes'] / 1e6))

    # An edit that is applied and then discarded is worse than one that is refused: the
    # viewer says "applied", the deck says nothing, and the .inp is the unedited one. A13
    # rebuilt PARAMS from the form widgets, which did exactly that.
    def _cell(tag):
        m = [c for c in cells if tag in c]
        return m[0] if m else ''

    # Match the titles, not a mention: A7 now prints the words "re-run A12b", and a
    # substring match on 'A12b' picked that cell up instead of the apply cell.
    a7, a12b, a13 = (_cell('A7 · PREVIEW'), _cell('A12b · Rebuild'),
                     _cell('A13 · Build'))
    chk('the three cells the edit path runs through are all present',
        bool(a7) and bool(a12b) and bool(a13),
        'A7 %s, A12b %s, A13 %s' % (bool(a7), bool(a12b), bool(a13)))
    chk('the build cell honours the edits the apply cell made',
        bool(a13) and 'EDITED_SETTINGS' in a13
        and 'PARAMS = make_params()' not in a13,
        'A13 reads EDITED_SETTINGS: %s; rebuilds blindly: %s'
        % ('EDITED_SETTINGS' in a13, 'PARAMS = make_params()' in a13))
    chk('and the apply cell publishes them for it',
        bool(a12b) and 'EDITED_SETTINGS =' in a12b and 'EDITED_BASE' in a12b)
    chk('having no edits does not stop the notebook before the deck is built',
        'SystemExit' not in a12b, 'A12b raises SystemExit: %s'
        % ('SystemExit' in a12b))
    chk('and re-running the preview resets the edit state',
        bool(a7) and 'EDITED_SETTINGS = {}' in a7)

    # ...and nothing downstream of the apply cell may report a widget value either. The
    # build report printed the block size from WP_LENGTH_MM and the infeed from
    # DEPTH_OF_CUT_UM, so an edited deck was described by the numbers it had replaced.
    _widget_names = {f.widget for f in EFIELDS if f.widget}
    _after = cells[cells.index(a12b) + 1:] if a12b in cells else []
    leaks = []
    for c in _after:
        bare = '\n'.join(l.split('#')[0] for l in c.splitlines())
        for nm in sorted(_widget_names):
            if re.search(r'\b%s\b' % nm, bare):
                leaks.append('%s in %s' % (nm, c.split(chr(10))[0][:40]))
    chk('no cell after the apply cell describes the deck with a widget value',
        not leaks, '; '.join(leaks) or 'all read PARAMS/INFO')

    # The other branch of that cell: a widget changed after the edit was applied, so the
    # edited fields are re-applied on top of the new widget values and both survive.
    w0 = DeckParams(name='drift0', analysis=AnalysisParams(**AN), **ref_kw)
    w1 = dataclasses.replace(w0, name='drift1', rim_depth_mm=w0.rim_depth_mm * 1.5)
    merged = params_from_settings({'clearance_um': 0.45}, w1)
    chk('a widget changed after an edit keeps both the edit and the widget',
        abs(merged.clearance_um - 0.45) < 1e-12
        and abs(merged.rim_depth_mm - w1.rim_depth_mm) < 1e-12
        and merged.name == 'drift1',
        'standoff %.3f um (edited), rim %.4f mm (widget)'
        % (merged.clearance_um, merged.rim_depth_mm))

    # ---------------- C15 the Apply reply ----------------
    print()
    print('C15 APPLY REPLY  (Colab formats a callback return value before the page sees')
    print('    it, so the reply must survive that or the viewer learns nothing)')
    from semgrit.editable import CommitReply, commit_reply

    try:
        from IPython.core.formatters import DisplayFormatter
    except ImportError:
        DisplayFormatter = None
    if DisplayFormatter is None:
        warn('IPython is here, so the reply mimetypes can be checked', False,
             'install ipython to cover the live-kernel reply')
    else:
        fmt = DisplayFormatter()
        good, _ = fmt.format(commit_reply(False, error='depth of cut 9 um exceeds 3.8'))
        chk('a reply arrives under application/json, which is what the viewer reads',
            'application/json' in good
            and good['application/json']['error'].startswith('depth of cut'),
            ', '.join(sorted(good)))
        chk('and the same payload arrives as text, exactly, behind the sentinel',
            good['text/plain'].startswith(CommitReply.SENTINEL)
            and json.loads(good['text/plain'][len(CommitReply.SENTINEL):])
            == good['application/json'],
            good['text/plain'][:58])
        # The regression itself: this is what the callback used to return, and why the
        # viewer reported every Apply -- successful ones included -- as a refusal.
        plain, _ = fmt.format({'ok': True, 'message': 'applied'})
        chk('a plain dict really does not, which is the bug this replaces',
            'application/json' not in plain, ', '.join(sorted(plain)))
        ok_reply, _ = fmt.format(commit_reply(True, message='2 changed (live)'))
        chk('a success says so and carries its message',
            ok_reply['application/json'] == {'ok': True, 'message': '2 changed (live)',
                                             'error': ''},
            str(ok_reply['application/json'])[:58])
        silent, _ = fmt.format(commit_reply(False))
        chk('a refusal is never allowed to arrive without a reason',
            len(silent['application/json']['error']) > 20,
            silent['application/json']['error'][:58])
        # And it must be the real callback that does this, not just the helper: the
        # notebook cell is the only place with a live kernel and no test coverage.
        cell = [c for c in cells if 'register_callback' in c]
        chk('the notebook callback returns the reply object, not a dict',
            bool(cell) and 'commit_reply as _cad_reply' in cell[0]
            and 'return {"ok"' not in cell[0],
            'return _cad_reply(...) in the viewer cell' if cell else 'cell not found')

    # Having no viewer edits is the ordinary case -- most runs never open the panel. A12b
    # raising on it stops the notebook before A13 and no deck is ever written, which is a
    # whole run lost to a non-event. The end-to-end test catches it in eight minutes; this
    # catches it in milliseconds.
    a12b = [c for c in cells if 'APPLY_VIEWER_EDITS' in c and 'no edits found' in c]
    chk('the no-edits case is reported, not raised, so A13 still builds',
        bool(a12b) and 'SystemExit' not in a12b[0],
        ('A12b found, %d SystemExit' % a12b[0].count('SystemExit')) if a12b
        else 'A12b not found')

    # ---------------- C11 the measurement cache ----------------
    print()
    print('C11 MEASUREMENT CACHE  (re-running after a wheel change must be free,')
    print('    and must never serve the wrong grains)')
    from semgrit.quick import SIMPLE_MEASURE, measure_images

    cdir = os.path.join(OUT, 'cache')
    img = [os.path.abspath('DIAMOND_11.tif')]
    if os.path.exists(img[0]):
        t0 = time.time()
        m1 = measure_images(img, cdir, log=lambda *a: None, **SIMPLE_MEASURE)
        t_cold = time.time() - t0
        t0 = time.time()
        m2 = measure_images(img, cdir, log=lambda *a: None, **SIMPLE_MEASURE)
        t_warm = time.time() - t0
        chk('a repeat run is served from the cache',
            m2['cached'] and not m1['cached'] and t_warm < max(0.5, t_cold / 5),
            '%.2f s cold, %.2f s warm' % (t_cold, t_warm))
        chk('and returns the same grains',
            len(m1['solids']) == len(m2['solids'])
            and all(abs(a.height_um - b.height_um) < 1e-12
                    for a, b in zip(m1['solids'], m2['solids'])),
            '%d solids' % len(m2['solids']))
        m3 = measure_images(img, cdir, log=lambda *a: None,
                            **dict(SIMPLE_MEASURE, max_vertices=32))
        chk('changing a measurement setting invalidates it', not m3['cached'])
        m4 = measure_images(img, cdir, cache=False, log=lambda *a: None,
                            **SIMPLE_MEASURE)
        chk('cache=False always re-measures', not m4['cached'])
        # Simple mode and the Advanced defaults must agree, or the two paths silently
        # build different grain libraries and neither ever hits the other's cache.
        a1 = next(c for c in cells if 'A1 ' in c.splitlines()[0])
        wid = {}
        for mm in re.finditer(r'^([A-Z][A-Z0-9_]*)\s*=\s*(.+?)\s*#@param', a1, re.M):
            try:
                wid[mm.group(1)] = ast.literal_eval(mm.group(2))
            except Exception:
                pass
        chk('simple mode meshes grains exactly as A1 does by default',
            wid.get('SIMPLIFY_UM') == SIMPLE_MEASURE['simplify_um']
            and wid.get('MAX_VERTICES') == SIMPLE_MEASURE['max_vertices']
            and (not wid.get('INCLUDE_BORDER_GRAINS'))
            == SIMPLE_MEASURE['interior_only'],
            'A1 %s/%s vs simple %s/%s'
            % (wid.get('SIMPLIFY_UM'), wid.get('MAX_VERTICES'),
               SIMPLE_MEASURE['simplify_um'], SIMPLE_MEASURE['max_vertices']))
    else:
        warn('measurement cache exercised', False, 'DIAMOND_11.tif not present')

    # ---------------- C7 determinism ----------------
    print()
    print('C7  DETERMINISM')
    d1 = build_deck(DeckParams(name='d1', analysis=AnalysisParams(**AN), **ref_kw),
                    solids, OUT)
    d2 = build_deck(DeckParams(name='d2', analysis=AnalysisParams(**AN), **ref_kw),
                    solids, OUT)
    s1 = [l for l in open(d1['path'], encoding='ascii') if not l.startswith('**')]
    s2 = [l for l in open(d2['path'], encoding='ascii') if not l.startswith('**')]
    chk('same settings twice give an identical deck', s1 == s2,
        '%d lines' % len(s1))
    p1 = plan_deck(pp, solids)
    chk('plan_deck is deterministic too',
        p1['n_grits'] == plan['n_grits']
        and abs(p1['ground_radius_mm'] - plan['ground_radius_mm']) < 1e-15)

    print()
    print('=' * 78)
    print('COLAB TOTAL: %d failure(s), %d warning(s)%s'
          % (len(FAIL), len(WARN), '' if not FAIL else '  -> ' + str(FAIL)))
    print('=' * 78)
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
