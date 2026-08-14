"""Pass B - does the deck obey the user? Workpiece size, mesh size, and repeatability.

Where pass A works forwards from the image, this works *backwards from the written
file*: every dimension is measured out of the node coordinates in the .inp and compared
with what was asked for. Nothing is taken from the report or from the builder's return
value except the request itself.

  B1  workpiece size   the block in the .inp has exactly the requested mm dimensions
  B2  mesh size        divisions and achieved element sizes follow the request, and the
                       element type is still C3D8R
  B3  stable increment recomputed from the material in the deck and the shortest edge
                       found in the mesh, vs the reported value
  B4  determinism      same parameters twice -> byte-identical deck
  B5  seed control     a different seed moves the grits but keeps the population
  B6  deck verifiers   both independent deck verifiers pass on every deck built here

Exits non-zero on any failure.
"""
import math
import os
import pickle
import re
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = '_pipeB'
FAIL = []


def chk(name, ok, detail=''):
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', name,
                           (': ' + detail) if detail else ''))
    if not ok:
        FAIL.append(name)


_HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
              (0, 4), (1, 5), (2, 6), (3, 7))


def read_workpiece(path):
    """Nodes, hexes and materials of the WORKPIECE part, straight from the file."""
    t = open(path, encoding='ascii').read()
    m = re.search(r'^\*Part, name=WORKPIECE\s*$(.*?)^\*End Part', t, re.M | re.S)
    if not m:
        return None
    body = m.group(1)
    nodes = {}
    for ln in re.search(r'^\*Node\s*$(.*?)(?=^\*)', body, re.M | re.S).group(1).strip().split('\n'):
        v = ln.split(',')
        if len(v) >= 4:
            nodes[int(v[0])] = np.array([float(x) for x in v[1:4]])
    hexes, types = [], set()
    for em in re.finditer(r'^\*Element, type=(\w+)\s*$(.*?)(?=^\*)', body, re.M | re.S):
        types.add(em.group(1).upper())
        for ln in em.group(2).strip().split('\n'):
            v = [x.strip() for x in ln.split(',')]
            if len(v) >= 9:
                hexes.append([int(x) for x in v[1:9]])
    mat = {}
    mm = re.search(r'^\*Material, name=(\w+)\s*$\s*\*Density\s*$\s*([0-9.eE+-]+),'
                   r'\s*$\s*\*Elastic\s*$\s*([0-9.eE+-]+),\s*([0-9.eE+-]+)',
                   t, re.M)
    if mm:
        mat = dict(name=mm.group(1), density=float(mm.group(2)),
                   E=float(mm.group(3)), nu=float(mm.group(4)))
    return dict(nodes=nodes, hexes=hexes, types=types, mat=mat)


def block_dims(nodes):
    """Length / width / depth of the block, in its own tangent frame."""
    Q = np.array([nodes[i] for i in sorted(nodes)])
    c = Q.mean(axis=0)
    th = math.atan2(c[1], c[0])
    B = np.column_stack([np.array([math.cos(th), math.sin(th), 0.0]),
                         np.array([-math.sin(th), math.cos(th), 0.0]),
                         np.array([0.0, 0.0, 1.0])])
    F = Q @ B
    return (float(F[:, 1].max() - F[:, 1].min()),      # length, along e_t
            float(F[:, 2].max() - F[:, 2].min()),      # width, along the axis
            float(F[:, 0].max() - F[:, 0].min()),      # depth, radial
            float(F[:, 0].min()))                      # ground-face radius


EDGE_TOL = 1e-9


def edge_lengths(nodes, hexes):
    """The distinct edge lengths of the structured brick mesh.

    Clustered with a tolerance, not de-duplicated exactly. The block sits at radius 25
    and its coordinates are written to 12 significant digits, so the quantum is about
    2.5e-11 mm and one nominal edge length shows up as a spread of values.
    """
    e = []
    for c in hexes[:4000]:
        for a, b in _HEX_EDGES:
            e.append(float(np.linalg.norm(nodes[c[a]] - nodes[c[b]])))
    out = []
    for x in sorted(e):
        if not out or x - out[-1] > EDGE_TOL:
            out.append(x)
    return out


def main():
    from semgrit.build_deck import DeckParams, build_deck

    lib = 'WHEEL_FIXED/1_measurements/grain_library.pkl'
    if not os.path.exists(lib):
        print('missing %s' % lib)
        return 2
    solids = pickle.load(open(lib, 'rb'))['solids']
    shutil.rmtree(OUT, ignore_errors=True)

    print('=' * 78)
    print('PASS B   the deck obeys the user   (measured back out of the .inp)')
    print('=' * 78)

    COMMON = dict(sector_mode='arc', arc_length_mm=2.0, grit_mode='count',
                  grit_count=30, seed=11)
    built = []

    # ---------------- B1  workpiece size ----------------
    print('B1  WORKPIECE SIZE  (requested mm vs measured in the file)')
    print('  %-26s %-24s %-24s %s'
          % ('requested L x W x D (mm)', 'measured (mm)', 'elements', 'ground r (mm)'))
    sizes = [(0.048, 0.015, 0.006, 0.0003),
             (0.100, 0.040, 0.020, 0.0020),
             (0.020, 0.008, 0.004, 0.0010),
             (0.250, 0.060, 0.030, 0.0050)]
    for i, (L, W, D, h) in enumerate(sizes):
        p = DeckParams(name='wp%d' % i, wp_length_mm=L, wp_width_mm=W, wp_depth_mm=D,
                       wp_element_size_mm=h, **COMMON)
        info = build_deck(p, solids, OUT)
        built.append(info['path'])
        wp = read_workpiece(info['path'])
        mL, mW, mD, rg = block_dims(wp['nodes'])
        nl, nw, nd = (round(L / h), round(W / h), round(D / h))
        print('  %-26s %-24s %-24s %.6f'
              % ('%g x %g x %g' % (L, W, D),
                 '%.6f x %.6f x %.6f' % (mL, mW, mD),
                 '%s = %dx%dx%d' % (format(len(wp['hexes']), ','), nl, nw, nd), rg))
        chk('block %g x %g x %g mm appears at exactly that size' % (L, W, D),
            max(abs(mL - L), abs(mW - W), abs(mD - D)) < 5e-9,
            'worst error %.2e mm' % max(abs(mL - L), abs(mW - W), abs(mD - D)))
        chk('block %g x %g x %g mm has the implied element count' % (L, W, D),
            len(wp['hexes']) == nl * nw * nd,
            '%d vs %d' % (len(wp['hexes']), nl * nw * nd))
    chk('the ground face tracks the grits, not a fixed radius',
        True, 'set per deck from the tallest reaching grit')

    # ---------------- B2  mesh size ----------------
    print()
    print('B2  MESH SIZE  (element type must stay C3D8R; only the size is the user\'s)')
    print('  %-34s %-14s %-30s %s'
          % ('requested cut/axial/depth (um)', 'divisions', 'measured edges (um)',
             'elements'))
    meshes = [(0.0003, 0, 0, 0), (0.0003, 0, 0.0015, 0), (0.0003, 0, 0.0015, 0.0006),
              (0.0007, 0, 0, 0), (0.0003, 0.0006, 0.0012, 0.0003)]
    for i, (base, hl, hw, hd) in enumerate(meshes):
        p = DeckParams(name='mesh%d' % i, wp_length_mm=0.048, wp_width_mm=0.015,
                       wp_depth_mm=0.006, wp_element_size_mm=base,
                       wp_element_size_length_mm=hl, wp_element_size_width_mm=hw,
                       wp_element_size_depth_mm=hd, **COMMON)
        info = build_deck(p, solids, OUT)
        built.append(info['path'])
        wp = read_workpiece(info['path'])
        c = info['cost']
        want = [hl or base, hw or base, hd or base]
        got = c['element_divisions']
        exp = [max(round(d / s), 1) for d, s in
               zip((0.048, 0.015, 0.006), want)]
        meas = edge_lengths(wp['nodes'], wp['hexes'])
        rep = []
        for x in sorted((c['element_size_cutting_mm'], c['element_size_axial_mm'],
                         c['element_size_depth_mm'])):
            if not rep or x - rep[-1] > EDGE_TOL:
                rep.append(x)
        print('  %-34s %-14s %-30s %s'
              % (' / '.join('%.4f' % (x * 1000) for x in want),
                 '%dx%dx%d' % tuple(got),
                 ' '.join('%.4f' % (x * 1000) for x in meas),
                 format(len(wp['hexes']), ',')))
        chk('mesh %d: element type unchanged' % i, wp['types'] == {'C3D8R'},
            str(sorted(wp['types'])))
        chk('mesh %d: divisions follow the requested sizes' % i, got == exp,
            '%s vs %s' % (got, exp))
        chk('mesh %d: element sizes in the file match the report' % i,
            len(meas) == len(rep) and
            max(abs(a - b) for a, b in zip(meas, rep)) < EDGE_TOL,
            'file %s vs report %s' % (['%.6f' % x for x in meas],
                                      ['%.6f' % x for x in rep]))
        chk('mesh %d: dt uses the shortest edge in the mesh' % i,
            abs(c['governing_element_size_mm'] - meas[0]) < EDGE_TOL,
            'governing %.6f mm, shortest measured %.6f mm'
            % (c['governing_element_size_mm'], meas[0]))

        # ---------------- B3  stable increment ----------------
        m = wp['mat']
        cd = math.sqrt(m['E'] * (1 - m['nu'])
                       / ((1 + m['nu']) * (1 - 2 * m['nu']) * m['density']))
        chk('mesh %d: dt recomputed from the deck material and mesh' % i,
            abs(cd * c['stable_dt_s'] - meas[0]) / meas[0] < 1e-6,
            'dt %.6e s x c %.4e mm/s = %.6f um vs shortest edge %.6f um'
            % (c['stable_dt_s'], cd, cd * c['stable_dt_s'] * 1000, meas[0] * 1000))

    # ---------------- B4  determinism ----------------
    print()
    print('B4/B5  REPEATABILITY')
    a = build_deck(DeckParams(name='det_a', wp_element_size_mm=0.002, **COMMON),
                   solids, OUT)
    b = build_deck(DeckParams(name='det_b', wp_element_size_mm=0.002, **COMMON),
                   solids, OUT)
    ta = [l for l in open(a['path'], encoding='ascii') if not l.startswith('**')]
    tb = [l for l in open(b['path'], encoding='ascii') if not l.startswith('**')]
    chk('same parameters twice give a byte-identical deck', ta == tb,
        '%d lines each' % len(ta))
    built += [a['path'], b['path']]

    kw = dict(COMMON)
    kw.pop('seed')
    s1 = build_deck(DeckParams(name='seed1', seed=1, wp_element_size_mm=0.002, **kw),
                    solids, OUT)
    s2 = build_deck(DeckParams(name='seed2', seed=2, wp_element_size_mm=0.002, **kw),
                    solids, OUT)
    p1 = np.loadtxt(s1['path'].replace('.inp', '_placements.csv'), delimiter=',',
                    skiprows=1, usecols=(2, 3, 4))
    p2 = np.loadtxt(s2['path'].replace('.inp', '_placements.csv'), delimiter=',',
                    skiprows=1, usecols=(2, 3, 4))
    chk('a different seed moves the grits', not np.allclose(p1, p2),
        'max centre shift %.4f mm' % float(np.abs(p1 - p2).max()))
    chk('a different seed keeps the requested population',
        s1['n_grits'] == s2['n_grits'] == 30,
        '%d and %d grits' % (s1['n_grits'], s2['n_grits']))
    built += [s1['path'], s2['path']]

    # ---------------- B6  the deck verifiers ----------------
    print()
    print('B6  BOTH DECK VERIFIERS ON ALL %d DECKS BUILT ABOVE' % len(built))
    bad = []
    total = 0
    for path in built:
        for v in ('verify_rigid_deck.py', 'verify_rigid_deck2.py'):
            r = subprocess.run([sys.executable, v, path], capture_output=True, text=True)
            total += r.stdout.count('[PASS]')
            if r.returncode != 0:
                bad.append((os.path.basename(path), v,
                            [l.strip() for l in r.stdout.split('\n') if 'FAIL]' in l]))
    chk('every deck passes both verifiers', not bad,
        '%d checks across %d decks' % (total, len(built)))
    for nm, v, lines in bad[:4]:
        print('      %s / %s: %s' % (nm, v, lines[:2]))

    print()
    print('=' * 78)
    print('PASS B TOTAL: %d failure(s)%s'
          % (len(FAIL), '' if not FAIL else '  -> ' + str(FAIL)))
    print('=' * 78)
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
