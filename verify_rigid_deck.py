"""Independent verifier for an all-rigid-wheel semgrit deck.

    python verify_rigid_deck.py <deck.inp>

Reads the finished file with its own parser and re-derives every geometric claim from
the node coordinates, so a bug in the writer cannot hide behind a shared helper.

  P1  mesh        element types, ids, connectivity, degenerate facets, hex Jacobians
  P2  syntax      block balance, every set/surface/material reference resolves
  P3  rigid body  exactly one, covers all wheel elements, ref node, mass, inertia
  P4  geometry     true arc, outward normals, tangency, zero penetration, overlaps
  P5  regression  every failure mode previously hit in this project, by name

Exits non-zero on any failure.
"""
import json
import math
import os
import re
import sys
from collections import Counter

import numpy as np

PATH = sys.argv[1] if len(sys.argv) > 1 else 'FINAL_RIGID/wheel_rigid_2mm.inp'

raw = open(PATH, 'rb').read()
t = raw.decode('ascii', errors='replace')
lines = t.split('\n')
kw = [l.split(',')[0].strip().lower() for l in lines
      if l.startswith('*') and not l.startswith('**')]
FAIL = []
WARN = []

# The report JSON written beside the deck says what was *asked for*. Anything it
# claims is re-derived from the node coordinates below rather than trusted; it is used
# only to know which optional parts of the deck should exist.
REPORT = {}
_rj = re.sub(r'\.inp$', '', PATH) + '_report.json'
if os.path.exists(_rj):
    try:
        REPORT = json.load(open(_rj))
    except ValueError:
        pass


def chk(name, ok, detail=''):
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', name, (': ' + detail) if detail else ''))
    if not ok:
        FAIL.append(name)


def warn(name, ok, detail=''):
    """A design guideline, not a correctness property: reported, never fatal."""
    print('  [%s] %s%s' % ('ok  ' if ok else 'WARN', name,
                           (': ' + detail) if detail else ''))
    if not ok:
        WARN.append(name)


def parse_ints(block):
    out = []
    for ln in block.strip().split('\n'):
        ln = ln.strip()
        if not ln or ln.startswith('*'):
            continue
        for x in ln.split(','):
            x = x.strip()
            if x:
                out.append(int(x))
    return out


def parse_part(body):
    nodes = {}
    for m in re.finditer(r'^\*Node\s*$(.*?)(?=^\*)', body, re.M | re.S):
        for ln in m.group(1).strip().split('\n'):
            v = [x.strip() for x in ln.split(',')]
            if len(v) >= 4:
                nodes[int(v[0])] = np.array([float(x) for x in v[1:4]])
    els, el_type = {}, {}
    for m in re.finditer(r'^\*Element,\s*type=(\w+)(?:,\s*elset=([^\s,\n]+))?\s*$(.*?)(?=^\*)',
                         body, re.M | re.S):
        ty = m.group(1).upper()
        for ln in m.group(3).strip().split('\n'):
            v = [x.strip() for x in ln.split(',')]
            if len(v) >= 2:
                eid = int(v[0])
                els[eid] = [int(x) for x in v[1:] if x]
                el_type[eid] = ty
    elsets = {}
    for m in re.finditer(r'^\*Elset,\s*elset=([^\s,\n]+)(,\s*generate)?\s*$(.*?)(?=^\*)',
                         body, re.M | re.S):
        name, gen = m.group(1).upper(), bool(m.group(2))
        vals = parse_ints(m.group(3))
        if gen:
            ids = []
            for i in range(0, len(vals), 3):
                a, b, s = vals[i], vals[i + 1], vals[i + 2]
                ids.extend(range(a, b + 1, s))
            elsets[name] = set(ids)
        else:
            elsets[name] = set(vals)
    # elsets declared inline on *Element
    for m in re.finditer(r'^\*Element,\s*type=\w+,\s*elset=([^\s,\n]+)\s*$(.*?)(?=^\*)',
                         body, re.M | re.S):
        ids = set(int(l.split(',')[0]) for l in m.group(2).strip().split('\n') if l.strip())
        elsets.setdefault(m.group(1).upper(), set()).update(ids)
    nsets = {}
    for m in re.finditer(r'^\*Nset,\s*nset=([^\s,\n]+)\s*$(.*?)(?=^\*)', body, re.M | re.S):
        nsets[m.group(1).upper()] = set(parse_ints(m.group(2)))
    surfs = {}
    for m in re.finditer(r'^\*Surface,\s*type=ELEMENT,\s*name=([^\s,\n]+)\s*$(.*?)(?=^\*)',
                         body, re.M | re.S):
        rows = [[x.strip().upper() for x in l.split(',')]
                for l in m.group(2).strip().split('\n') if l.strip()]
        surfs[m.group(1).upper()] = rows
    return dict(nodes=nodes, els=els, ty=el_type, elsets=elsets, nsets=nsets, surfs=surfs,
                body=body)


parts = {}
for m in re.finditer(r'^\*Part,\s*name=([^\s,\n]+)\s*$(.*?)^\*End Part', t, re.M | re.S):
    parts[m.group(1).upper()] = parse_part(m.group(2))

print('=' * 78)
print('VERIFYING  %s   (%.2f MB)' % (PATH, len(raw) / 1e6))
print('=' * 78)
print('P1  MESH')

RIGID = {'R3D3', 'R3D4'}
SOLID = {'C3D8R'}
INERTIA = {'MASS', 'ROTARYI'}
types = Counter()
for p in parts.values():
    types.update(p['ty'].values())
chk('Explicit element types only', set(types) <= RIGID | SOLID | INERTIA,
    ' '.join('%s=%s' % (k, format(v, ',')) for k, v in sorted(types.items())))
chk('no Standard-only element types',
    not (set(types) & {'C3D8', 'C3D8I', 'R3D3S', 'C3D4H'}), 'none present')

dup = bad_conn = 0
for nm, p in parts.items():
    ids = list(p['els'])
    dup += len(ids) - len(set(ids))
    for e, c in p['els'].items():
        if any(n not in p['nodes'] for n in c):
            bad_conn += 1
chk('element ids unique within each part', dup == 0, '%d duplicates' % dup)
chk('every connectivity node exists', bad_conn == 0, '%d bad elements' % bad_conn)

worst = math.inf
nfac = 0
for nm, p in parts.items():
    for e, c in p['els'].items():
        if p['ty'][e] not in RIGID:
            continue
        q = [p['nodes'][i] for i in c]
        n = len(q)
        nfac += 1
        worst = min(worst, min(float(np.linalg.norm(q[i] - q[(i + 1) % n]))
                              for i in range(n)))
chk('no coincident-node rigid facets', worst >= 1e-5,
    '%s facets, min edge %.4e mm' % (format(nfac, ','), worst))

zero_area = 0
for nm, p in parts.items():
    for e, c in p['els'].items():
        if p['ty'][e] not in RIGID:
            continue
        q = np.array([p['nodes'][i] for i in c])
        a = 0.0
        for i in range(1, len(q) - 1):
            a += 0.5 * float(np.linalg.norm(np.cross(q[i] - q[0], q[i + 1] - q[0])))
        if a < 1e-12:
            zero_area += 1
chk('no zero-area rigid facets', zero_area == 0, '%d found' % zero_area)


def hexvol(nd, c):
    q = [nd[i] for i in c]
    v = 0.0
    for a, b, cc, dd in ((0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6),
                         (3, 4, 6, 7), (1, 4, 5, 6)):
        v += float(np.dot(np.cross(q[b] - q[a], q[cc] - q[a]), q[dd] - q[a])) / 6.0
    return v


WPN = 'WORKPIECE'
HAS_WP = WPN in parts
chk('workpiece presence agrees with the report',
    ('has_workpiece' not in REPORT) or bool(REPORT['has_workpiece']) == HAS_WP,
    'deck %s a workpiece' % ('has' if HAS_WP else 'is wheel-only'))
if HAS_WP:
    vv = np.array([hexvol(parts[WPN]['nodes'], c) for e, c in parts[WPN]['els'].items()
                   if parts[WPN]['ty'][e] == 'C3D8R'])
    chk('workpiece hex volumes all positive', bool((vv > 0).all()),
        '%s elems, min %.4e mm3' % (format(len(vv), ','), vv.min()))
chk('workpiece is the only deformable part',
    all(not (set(p['ty'].values()) & SOLID) for nm, p in parts.items() if nm != WPN),
    'wheel carries no solid elements')

print()
print('P2  SYNTAX / REFERENCES')
c = Counter(kw)
chk('blocks balanced',
    c['*part'] == c['*end part'] and c['*instance'] == c['*end instance']
    and c['*assembly'] == c['*end assembly'] == 1,
    'part %d/%d  instance %d/%d  assembly %d/%d'
    % (c['*part'], c['*end part'], c['*instance'], c['*end instance'],
       c['*assembly'], c['*end assembly']))
chk('pure ASCII, no en-dash, lines <= 256',
    all(ch < 128 for ch in raw) and max(len(l) for l in lines) <= 256,
    'longest line %d chars' % max(len(l) for l in lines))

bs = [(nm, s, r) for nm, p in parts.items() for s, rows in p['surfs'].items()
      for r in rows if r[0] not in p['elsets']]
chk('part *Surface -> *Elset in same part', not bs, str(bs[:3]))
bside = [(nm, s, r) for nm, p in parts.items() for s, rows in p['surfs'].items()
         for r in rows if len(r) < 2]
chk('part *Surface rows carry a side/face id', not bside, str(bside[:3]))

asm = t.split('*Assembly')[1].split('*End Assembly')[0]
inst = {m.group(1).upper(): m.group(2).upper() for m in
        re.finditer(r'\*Instance,\s*name=([^\s,]+),\s*part=([^\s,\n]+)', asm)}
chk('instances present (not 0)', len(inst) == (2 if HAS_WP else 1),
    '%d: %s' % (len(inst), sorted(inst)))

ba = []
for m in re.finditer(r'^\*Surface,\s*type=ELEMENT,\s*name=([^\s,\n]+)\s*$(.*?)(?=^\*)',
                     asm, re.M | re.S):
    for ln in m.group(2).strip().split('\n'):
        pr = [x.strip().upper() for x in ln.split(',') if x.strip()]
        if not pr:
            continue
        if '.' not in pr[0]:
            ba.append((m.group(1), 'not instance-qualified'))
        else:
            i, s = pr[0].split('.', 1)
            if i not in inst:
                ba.append((m.group(1), 'unknown instance ' + i))
            elif s not in parts[inst[i]]['elsets']:
                ba.append((m.group(1), 'elset %s missing in %s' % (s, inst[i])))
            elif len(pr) < 2:
                ba.append((m.group(1), 'no side/face id'))
chk('assembly *Surface -> instance.ELSET + side id', not ba, str(ba[:3]))

bn = []
for m in re.finditer(r'\*Nset,\s*nset=([^\s,\n]+),\s*instance=([^\s,\n]+)\s*\n([^\n*]+)',
                     asm, re.I):
    i = m.group(2).upper()
    ref = m.group(3).split(',')[0].strip().upper()
    if i not in inst or ref not in parts[inst[i]]['nsets']:
        bn.append(m.group(1))
chk('assembly *Nset,instance= resolves', not bn, str(bn[:3]))

mats = set(x.upper() for x in re.findall(r'^\*Material,\s*name=([^\s,\n]+)', t, re.M))
used = set(x.upper() for x in
           re.findall(r'\*Solid Section,[^\n]*material=([^\s,\n]+)', t, re.I))
chk('every section material is defined', not (used - mats),
    'sections use %s, defined %s' % (sorted(used), sorted(mats)))
secs = re.findall(r'^\*Solid Section,[^\n]*elset=([^\s,\n]+)', t, re.M | re.I)
chk('only the workpiece has a *Solid Section',
    [s.upper() for s in secs] == (['WP_ALL'] if HAS_WP else []),
    '%d section(s): %s' % (len(secs), secs))

print()
print('P3  RIGID BODY')
W = parts['WHEEL']
rb = re.findall(r'\*Rigid Body,\s*ref node=([^\s,\n]+),\s*elset=([^\s,\n]+)', t)
chk('exactly one *Rigid Body', len(rb) == 1, '%d found' % len(rb))
rn, re_ = (rb[0][0].upper(), rb[0][1].upper()) if rb else ('', '')
chk('*Rigid Body ref node nset exists', rn in W['nsets'], rn)
chk('*Rigid Body elset exists', re_ in W['elsets'], re_)

rigid_ids = set(e for e, ty in W['ty'].items() if ty in RIGID)
inert_ids = set(e for e, ty in W['ty'].items() if ty in INERTIA)
chk('rigid elset covers every wheel facet exactly',
    W['elsets'].get(re_) == rigid_ids,
    '%s facets in elset, %s R3D3/R3D4 in part'
    % (format(len(W['elsets'].get(re_, ())), ','), format(len(rigid_ids), ',')))
chk('mass/inertia elements excluded from the rigid elset',
    not (W['elsets'].get(re_, set()) & inert_ids), '%d inertia elements' % len(inert_ids))

refn = sorted(W['nsets'].get(rn, []))
chk('ref node is a single node', len(refn) == 1, str(refn[:3]))
rid = refn[0] if refn else -1
chk('ref node lies on the wheel axis at the origin',
    rid in W['nodes'] and float(np.linalg.norm(W['nodes'][rid])) < 1e-12,
    '%s' % (W['nodes'].get(rid)))
in_facet = any(rid in c for e, c in W['els'].items() if W['ty'][e] in RIGID)
chk('ref node belongs to no rigid facet', not in_facet)
chk('ref node carries the mass and rotary inertia',
    all(rid in c for e, c in W['els'].items() if W['ty'][e] in INERTIA)
    and len(inert_ids) == 2,
    '%d inertia elements on node %d' % (len(inert_ids), rid))

mval = re.search(r'\*Mass,\s*elset=([^\s,\n]+)\s*\n\s*([0-9.eE+-]+)', t)
chk('*Mass elset exists and value > 0',
    bool(mval) and mval.group(1).upper() in W['elsets'] and float(mval.group(2)) > 0,
    '%s = %s tonne' % (mval.group(1), mval.group(2)) if mval else 'missing')
riv = re.search(r'\*Rotary Inertia,\s*elset=([^\s,\n]+)\s*\n([^\n]+)', t)
if riv:
    v = [float(x) for x in riv.group(2).split(',') if x.strip()]
    I = np.array([[v[0], v[3], v[4]], [v[3], v[1], v[5]], [v[4], v[5], v[2]]])
    ev = np.linalg.eigvalsh(I)
    tri = (v[0] + v[1] >= v[2] - 1e-18 and v[1] + v[2] >= v[0] - 1e-18
           and v[0] + v[2] >= v[1] - 1e-18)
    chk('*Rotary Inertia elset exists', riv.group(1).upper() in W['elsets'], riv.group(1))
    chk('inertia tensor positive definite', bool(ev.min() > 0), 'eigenvalues %s' % ev)
    chk('inertia tensor satisfies the triangle inequalities', tri,
        'I11=%.3e I22=%.3e I33=%.3e' % (v[0], v[1], v[2]))
else:
    chk('*Rotary Inertia present', False, 'missing')

orphan = set(W['nodes']) - set(n for c in W['els'].values() for n in c)
chk('no orphan nodes in the wheel part', not orphan, '%d orphans' % len(orphan))

print()
print('P4  GEOMETRY  (re-derived from the node coordinates)')

# --- separate the shell from the grits by element type -----------------------
shell_e = {e: c for e, c in W['els'].items() if W['ty'][e] == 'R3D4'}
grit_e = {e: c for e, c in W['els'].items() if W['ty'][e] == 'R3D3'}
sh_nodes = sorted(set(n for c in shell_e.values() for n in c))

# An include_bond=False deck has no rim, so every check below that re-derives the
# wheel from bond-shell nodes has nothing to read. Skipping is stated out loud
# rather than passed silently: a gate that reports PASS on a file it never looked
# at is worse than one that says it did not look.
BONDLESS = not shell_e
if BONDLESS:
    print('      no R3D4 bond shell in this deck (include_bond=False):')
    print('      skipping the rim geometry, sector and shell-normal checks.')
    chk('bondless deck still carries grit facets', bool(grit_e),
        '%d R3D3 facets' % len(grit_e))
else:
    P = np.array([W['nodes'][i] for i in sh_nodes])
    rr = np.hypot(P[:, 0], P[:, 1])
    # Cluster with a tolerance rather than exact equality: the coordinates in the file are
    # quantised by the write format, so a node on a circle of radius 25 recovers that
    # radius only to within the quantum, not to the last bit.
    TOL_R = 1e-9
    srt = np.sort(rr)
    radii = [float(srt[0])]
    for x in srt[1:]:
        if x - radii[-1] > TOL_R:
            radii.append(float(x))
    # A thin shell has two radial levels; a deep wedge subdivided radially has one per
    # division, and its axial and cut faces legitimately carry the intermediate ones. What
    # must hold either way is that every node sits on one of a few *exact* circles, so the
    # outer face is a true arc and not a faceted polygon.
    R_out, R_in = max(radii), min(radii)
    chk('bond shell nodes lie on a few exact circles (true arc, not faceted)',
        2 <= len(radii) <= 64, '%d radial levels (tol %.0e): %.6f .. %.6f mm'
        % (len(radii), TOL_R, R_in, R_out))
    lev = np.array(radii)
    dev = float(np.abs(lev[np.abs(rr[:, None] - lev[None, :]).argmin(axis=1)] - rr).max())
    chk('shell node radii deviate from those exact circles by < 1 pm', dev < 1e-9,
        'max deviation %.3e mm across %d levels' % (dev, len(radii)))
    n_out = int((np.abs(rr - R_out) < TOL_R).sum())
    chk('the outer face is a full exact circle at the largest radius', n_out > 0,
        '%d nodes on r = %.6f mm' % (n_out, R_out))
    # A full wheel is identified from the mesh, not from the report: it is the only case
    # with no sector cut faces, because the last column of nodes is the first.
    FULL = 'ES_BOND_SECTOR_START' not in W['elsets']
    th = np.degrees(np.arctan2(P[:, 1], P[:, 0]))
    rim = R_out - R_in
    if FULL:
        span = 360.0
        arc = 2.0 * math.pi * R_out
        sag = None
        print('      R = %.6f mm, FULL WHEEL, circumference %.6f mm, rim %.6f mm'
              % (R_out, arc, rim))
        # Cluster before differencing: every angular station carries many nodes (one per
        # radius and axial station) whose recovered angles differ in the last bits, so a
        # raw diff is mostly ~1e-14 noise and its median is meaningless.
        uth = np.sort(th)
        keep = np.r_[True, np.diff(uth) > 1e-7]
        uth = uth[keep]
        # Include the wrap across +-180 deg: that gap is the seam, if there is one.
        gaps = np.r_[np.diff(uth), 360.0 - (uth[-1] - uth[0])]
        chk('full wheel closes without a seam',
            float(gaps.max()) <= 1.5 * float(np.median(gaps)) + 1e-9,
            '%d angular stations, largest gap %.6f deg vs median %.6f'
            % (len(uth), gaps.max(), np.median(gaps)))
    else:
        span = float(th.max() - th.min())
        arc = R_out * math.radians(span)
        sag = arc * arc / (8.0 * R_out)
        print('      R = %.6f mm, sector %.6f deg, arc %.6f mm, rim %.6f mm'
              % (R_out, span, arc, rim))
        print('      sagitta %.3f um = %.0f%% of the rim depth'
              % (sag * 1000, 100 * sag / rim))
    chk('sector agrees with the report', ('resolved_sector_deg' not in REPORT)
        or abs(span - float(REPORT['resolved_sector_deg'])) < 1e-6,
        '%.6f deg in the mesh' % span)
    chk('arc length agrees with the report', ('arc_length_mm' not in REPORT)
        or abs(arc - float(REPORT['arc_length_mm'])) < 1e-6, '%.6f mm' % arc)
    if sag is not None:
        # Not a correctness property -- a flat-looking sector is still a correct sector.
        # It is flagged because "it renders as a rectangle" is the usual complaint.
        #
        # Two independent ways for a sector to read as an arc, and either suffices:
        #  * a thin rim, where the outer face visibly bows across the band (sagitta > rim);
        #  * a wide angle, where the two radial cut faces obviously converge. Past about
        #    10 deg that convergence carries the shape on its own, and demanding
        #    sagitta > rim there would wrongly condemn a deep chunky wedge -- exactly the
        #    geometry that looks most like a slice of a real wheel.
        bows = sag > rim
        wide = span >= 10.0
        warn('sector reads as an arc rather than a rectangle', bows or wide,
             'sector %.3f deg, sagitta %.3f um vs rim %.3f um; widen the angle, lengthen '
             'the arc, shrink the radius or thin the rim'
             % (span, sag * 1000, rim * 1000))
        print('      reads as an arc because %s'
              % (' and '.join(filter(None, [
                  'the outer face bows %.2fx its rim depth' % (sag / rim) if bows else '',
                  'the cut faces converge at %.1f deg' % span if wide else '']))
                 or 'NOTHING -- it will look flat'))

    # --- every shell quad normal must point out of the sector -------------------
    zmax = float(P[:, 2].max())
    zmin = float(P[:, 2].min())
    th0, th1 = math.radians(float(th.min())), math.radians(float(th.max()))
    badn = []
    for e, cc in shell_e.items():
        q = np.array([W['nodes'][i] for i in cc])
        n = np.cross(q[1] - q[0], q[2] - q[0])
        n = n / max(float(np.linalg.norm(n)), 1e-30)
        m = q.mean(axis=0)
        rm = float(np.hypot(m[0], m[1]))
        er = np.array([m[0] / rm, m[1] / rm, 0.0])
        et = np.array([-m[1] / rm, m[0] / rm, 0.0])
        # Classify by the quad's four *nodes*, not its centroid. A quad spanning an arc
        # has its centroid on the chord, 5e-7 mm inside the circle here, so a
        # centroid-radius test with a quantum-sized tolerance rejects every face on the
        # cylinder. The nodes themselves lie on the circle exactly.
        nr = np.hypot(q[:, 0], q[:, 1])
        nth = np.arctan2(q[:, 1], q[:, 0])
        if np.all(np.abs(nr - R_out) < 1e-9):
            want = er
        elif np.all(np.abs(nr - R_in) < 1e-9):
            want = -er
        elif np.all(np.abs(q[:, 2] - zmin) < 1e-9):
            want = np.array([0.0, 0.0, -1.0])
        elif np.all(np.abs(q[:, 2] - zmax) < 1e-9):
            want = np.array([0.0, 0.0, 1.0])
        elif np.all(np.abs(nth - th0) < 1e-9):
            want = -et
        elif np.all(np.abs(nth - th1) < 1e-9):
            want = et
        else:
            badn.append((e, 'quad on no recognised face'))
            continue
        if float(np.dot(n, want)) < 0.99:
            badn.append((e, 'normal %s wanted %s' % (np.round(n, 4), np.round(want, 4))))
    chk('every bond shell quad normal points outward', not badn,
        '%d shell quads checked, %d bad %s' % (len(shell_e), len(badn), badn[:2]))

# --- grit facet orientation: outward, by the divergence theorem -------------
groups = {}
for e, cc in grit_e.items():
    groups.setdefault(min(cc) // 100000, []).append((e, cc))
by_node = {}
for e, cc in grit_e.items():
    by_node.setdefault(tuple(sorted(cc)), e)
# group facets into grits via connected components on shared nodes
adj = {}
for e, cc in grit_e.items():
    for n in cc:
        adj.setdefault(n, []).append(e)
seen, grits = set(), []
for e0 in grit_e:
    if e0 in seen:
        continue
    stack, comp = [e0], []
    seen.add(e0)
    while stack:
        e = stack.pop()
        comp.append(e)
        for n in grit_e[e]:
            for f in adj[n]:
                if f not in seen:
                    seen.add(f)
                    stack.append(f)
    grits.append(comp)
neg = 0
for comp in grits:
    v = 0.0
    for e in comp:
        a, b, cc2 = (W['nodes'][i] for i in grit_e[e])
        v += float(np.dot(np.cross(a, b), cc2)) / 6.0
    if v <= 0:
        neg += 1
chk('grits are closed with outward normals (positive enclosed volume)', neg == 0,
    '%d grits found by connectivity, %d mis-oriented' % (len(grits), neg))
chk('grit count matches the header', len(grits) == int(
    re.search(r'\*\* grits\s+:\s*(\d+)', t).group(1)),
    '%d by connectivity' % len(grits))

# --- workpiece box in the tangent frame -------------------------------------
if not HAS_WP:
    print('      wheel-only deck: no workpiece, so no tangency or penetration test')
    a_lo = best = None
else:
  WPp = parts[WPN]
  Q = np.array([WPp['nodes'][i] for i in sorted(WPp['nodes'])])
  # Recover the block's own tangent frame from its geometry instead of assuming
  # mid-arc: the writer relocates the block to the tallest grit when the nominal
  # angle has nothing under it, and a full wheel has no mid-arc at all.
  Qc = Q.mean(axis=0)
  theta_c = math.atan2(Qc[1], Qc[0])
  e_r = np.array([math.cos(theta_c), math.sin(theta_c), 0.0])
  e_t = np.array([-math.sin(theta_c), math.cos(theta_c), 0.0])
  e_z = np.array([0.0, 0.0, 1.0])
  B = np.column_stack([e_r, e_t, e_z])
  Qf = Q @ B
  a_lo, a_hi = float(Qf[:, 0].min()), float(Qf[:, 0].max())
  hb, hz = float(np.abs(Qf[:, 1]).max()), float(np.abs(Qf[:, 2]).max())
  chk('workpiece is an axis-aligned box in the tangent frame (centred on the arc)',
      abs(float(Qf[:, 1].min()) + hb) < 1e-9 and abs(float(Qf[:, 2].min()) + hz) < 1e-9,
      'b +-%.6f mm, z +-%.6f mm, a %.6f..%.6f mm' % (hb, hz, a_lo, a_hi))

  # --- tangency and penetration, re-derived by dense sampling ----------------
  NS = 12
  def clip_to_footprint(poly):
      """Sutherland-Hodgman against |t| <= hb and |z| <= hz."""
      for axis, lim, keep_low in ((1, hb, True), (1, -hb, False),
                                  (2, hz, True), (2, -hz, False)):
          if not len(poly):
              return poly
          out = []
          n = len(poly)
          for i in range(n):
              a, b = poly[i], poly[(i + 1) % n]
              av = a[axis] - lim if keep_low else lim - a[axis]
              bv = b[axis] - lim if keep_low else lim - b[axis]
              a_in, b_in = av <= 1e-12, bv <= 1e-12
              if a_in:
                  out.append(a)
              if a_in != b_in:
                  t = av / (av - bv)
                  out.append(a + t * (b - a))
          poly = np.array(out) if out else np.zeros((0, 3))
      return poly

  # The radial coordinate is a linear function over a planar facet and the footprint
  # is convex, so the maximum over the clipped facet is exactly at a vertex of the
  # clipped polygon. Sampling the triangle on a grid instead -- which this used to do
  # -- misses a narrow clipped corner by up to the sample spacing, which showed up as
  # a few nanometres of phantom gap on some layouts and nothing on others.
  best = -math.inf
  pen = 0
  worst_pen = 0.0
  for comp in grits:
      for e in comp:
          tri = np.array([W['nodes'][i] for i in grit_e[e]]) @ B
          # Trivial reject before clipping. A dressed wheel can carry millions of
          # facets and only a few thousand lie under the block; clipping every one of
          # them in Python turned a ten-minute check into hours. This is exact: a
          # triangle entirely on the far side of a footprint edge cannot contribute.
          t_, z_ = tri[:, 1], tri[:, 2]
          if t_.min() > hb or t_.max() < -hb or z_.min() > hz or z_.max() < -hz:
              continue
          poly = clip_to_footprint(tri)
          if not len(poly):
              continue
          m = float(poly[:, 0].max())
          if m > best:
              best = m
          d = m - a_lo
          if d > 1e-9:
              pen += 1
              worst_pen = max(worst_pen, d)
  chk('no grit penetrates the workpiece', pen == 0,
      '%d facets past the ground face, worst %.4e mm' % (pen, worst_pen))
  # Tangent, or exactly the standoff the build was asked for. The invariant is that
  # the gap is *chosen*, not that it is zero -- a deck built with a standoff is
  # correct precisely when the face sits that far off the tallest grit.
  _sd = float(REPORT.get('clearance_um') or 0.0) / 1000.0
  chk('tallest grit inside the footprint sits at the requested standoff',
      abs((a_lo - best) - _sd) < 2e-6,
      'gap %.4f nm, asked for %.4f nm  (exact, facets clipped to the footprint)'
      % ((a_lo - best) * 1e6, _sd * 1e6))
  print('      ground face at a = %.9f mm, closest grit material a = %.9f mm'
        % (a_lo, best))
  if not BONDLESS:
      print('      max engaging protrusion = %.4f um above r = %.4f'
            % ((best - R_out) * 1000, R_out))

# --- grits stay on the wheel -----------------------------------------------
gn = np.array([W['nodes'][i] for i in sorted(set(n for c in grit_e.values() for n in c))])
gth = np.degrees(np.arctan2(gn[:, 1], gn[:, 0]))
gr = np.hypot(gn[:, 0], gn[:, 1])
if BONDLESS:
    # "On the wheel" is a statement about the rim, and there is no rim to be on.
    print('      no bond rim: the arc, width and bore containment checks do not '
          'apply')
    print('      grits span %.4f..%.4f deg, r %.6f..%.6f mm'
          % (gth.min(), gth.max(), gr.min(), gr.max()))
elif FULL:
    # Every angle is on the wheel, so there is nothing to overhang.
    print('      full wheel: grits span %.4f..%.4f deg, no cut faces to overhang'
          % (gth.min(), gth.max()))
else:
    chk('grits lie within the bond arc', float(gth.min()) >= th.min() - 1e-9
        and float(gth.max()) <= th.max() + 1e-9,
        'grits %.4f..%.4f deg, bond %.4f..%.4f deg'
        % (gth.min(), gth.max(), th.min(), th.max()))
if not BONDLESS:
    chk('grits lie within the wheel width',
        float(np.abs(gn[:, 2]).max()) <= zmax + 1e-9,
        'grit |z| max %.6f mm, wheel half-width %.6f mm'
        % (np.abs(gn[:, 2]).max(), zmax))
    chk('no grit is buried below the bore', float(gr.min()) >= R_in - 1e-9,
        'min grit radius %.6f mm vs bore %.6f mm' % (gr.min(), R_in))

# --- grit-grit interpenetration, tested on the actual polyhedra -------------
# A bounding-sphere test is far too pessimistic for these shapes: the library's
# largest grain has a 6.6 um bounding radius but a 3.7 um mean footprint, so
# neighbours whose spheres overlap by microns need not touch at all. Test instead
# whether any vertex of one grit lies inside the closed surface of another, by
# counting ray crossings -- the grits are meshed at ~1 um, so a real
# interpenetration cannot avoid putting a vertex inside.
verts, tris, cent, brad = [], [], [], []
for comp in grits:
    ids = sorted(set(i for e in comp for i in grit_e[e]))
    idx = {n: k for k, n in enumerate(ids)}
    v = np.array([W['nodes'][n] for n in ids])
    verts.append(v)
    tris.append(np.array([[idx[n] for n in grit_e[e]] for e in comp]))
    c0 = v.mean(axis=0)
    cent.append(c0)
    brad.append(float(np.linalg.norm(v - c0, axis=1).max()))
cent = np.array(cent)
brad = np.array(brad)


def inside(pts, v, f):
    """Which of ``pts`` lie inside the closed triangle mesh (v, f). Ray along +X."""
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    e1, e2 = b - a, c - a
    d = np.array([1.0, 0.0, 0.0])
    h = np.cross(d, e2)
    det = np.einsum('ij,ij->i', e1, h)
    ok = np.abs(det) > 1e-18
    res = np.zeros(len(pts), dtype=bool)
    for k, p in enumerate(pts):
        s = p - a
        u = np.einsum('ij,ij->i', s, h) / np.where(ok, det, 1.0)
        q = np.cross(s, e1)
        vv = (q @ d) / np.where(ok, det, 1.0)
        tt = np.einsum('ij,ij->i', e2, q) / np.where(ok, det, 1.0)
        hit = ok & (u >= 0) & (u <= 1) & (vv >= 0) & (u + vv <= 1) & (tt > 1e-12)
        res[k] = bool(int(hit.sum()) % 2)
    return res


from scipy.spatial import cKDTree
pairs = cKDTree(cent).query_pairs(r=2 * brad.max(), output_type='ndarray')
sphere_ov = 0
if len(pairs):
    dd = np.linalg.norm(cent[pairs[:, 0]] - cent[pairs[:, 1]], axis=1)
    sphere_ov = int((dd < brad[pairs[:, 0]] + brad[pairs[:, 1]]).sum())
real_ov = 0
checked = 0
for i, j in pairs:
    lo_i, hi_i = verts[i].min(axis=0), verts[i].max(axis=0)
    lo_j, hi_j = verts[j].min(axis=0), verts[j].max(axis=0)
    if (hi_i < lo_j).any() or (hi_j < lo_i).any():
        continue          # axis-aligned boxes already disjoint
    checked += 1
    sel = np.all((verts[i] >= lo_j) & (verts[i] <= hi_j), axis=1)
    if sel.any() and inside(verts[i][sel], verts[j], tris[j]).any():
        real_ov += 1
        continue
    sel = np.all((verts[j] >= lo_i) & (verts[j] <= hi_i), axis=1)
    if sel.any() and inside(verts[j][sel], verts[i], tris[i]).any():
        real_ov += 1
print('      %d neighbour pairs, %d with overlapping bounding spheres, %d with '
      'overlapping boxes' % (len(pairs), sphere_ov, checked))
chk('no grit polyhedra actually interpenetrate', real_ov == 0,
    '%d box-overlapping pairs tested vertex-in-solid, %d interpenetrating'
    % (checked, real_ov))

print()
print('P5  REGRESSION  (every failure mode previously hit)')
RUN_READY = '*step' in kw
if not RUN_READY:
    forb = set(kw) & {'*step', '*boundary', '*contact', '*contact inclusions',
                      '*dynamic', '*static', '*shear failure', '*shearfailure',
                      '*plastic', '*drucker prager', '*surface interaction',
                      '*output', '*contact pair', '*restart'}
    chk('no step / BC / interaction / output keywords', not forb,
        str(sorted(forb)) if forb else 'geometry only, as asked')
else:
    chk('*Static / *Plastic / *Shear Failure never appear',
        not (set(kw) & {'*static', '*plastic', '*shear failure', '*shearfailure',
                        '*drucker prager'}), 'Explicit + user material only')
chk('*SHEARFAILURE outside *Plastic', '*shear failure' not in kw
    and 'SHEARFAILURE' not in t.upper())
# NUMBER INTERVAL is a field-output parameter. On *Output, history Abaqus rejects it
# ("THE PARAMETER HISTORY CANNOT BE USED WITH THE PARAMETER NUMBER INTERVAL") and the
# job dies in the pre-processor -- which is exactly how grind11b failed.
_hni = re.findall(r'^\*Output,\s*history[^\n]*number\s*interval', t, re.M | re.I)
chk('NUMBER INTERVAL on *Output, history', not _hni,
    'history output uses TIME INTERVAL' if not _hni else str(_hni[:2]))
# An unscoped *Element Output asks the rigid facets and the mass elements for stress
# and damage, which is twelve warnings per run.
_eo = re.findall(r'^\*Element Output([^\n]*)$', t, re.M | re.I)
_unscoped = [x for x in _eo if 'elset=' not in x.lower()]
chk('unscoped *Element Output over the rigid wheel', not _unscoped or not RUN_READY,
    'scoped to an elset' if not _unscoped else str(_unscoped[:2]))
for _m in re.finditer(r'^\*Element Output[^\n]*elset=([^\s,\n]+)', t, re.M | re.I):
    _nm = _m.group(1).upper()
    chk('the elset *Element Output names exists (%s)' % _nm,
        bool(re.search(r'^\*Elset,\s*elset=%s\b' % re.escape(_nm), t, re.M | re.I)))
chk('coincident facet nodes (min edge >= 1e-5 mm)', worst >= 1e-5,
    'min edge %.4e mm' % worst)
chk('assembly surface referencing a part *Surface', not ba, 'all resolve to elsets')
chk('assembly nset referencing a missing part nset', not bn, 'all resolve')
chk('Import -> Part giving 0 instances', len(inst) == (2 if HAS_WP else 1),
    'assembly holds %d instances; load with ModelFromInputFile' % len(inst))
chk('non-ASCII / en-dash in the deck', all(ch < 128 for ch in raw))
chk('bond left deformable (the standing instruction)',
    not (set(parts['WHEEL']['ty'].values()) & SOLID)
    and len(rb) == 1 and len(parts) == (2 if HAS_WP else 1),
    'wheel = 1 rigid body, %d parts total, %s' % (
        len(parts), ('only %s deformable' % WPN) if HAS_WP else 'nothing deformable'))

chk('per-grit rigid bodies needing hundreds of BCs',
    len(rb) == 1, 'one ref node drives all %d grits + bond' % len(grits))

if RUN_READY:
    print()
    print('P6  RUN-READY  (submittable from the terminal, no CAE)')
    c = Counter(kw)
    chk('exactly one step, opened and closed',
        c['*step'] == c['*end step'] == 1, '%d/%d' % (c['*step'], c['*end step']))
    dyn = re.search(r'^\*Dynamic, Explicit\s*\n\s*,\s*([0-9.eE+-]+)', t, re.M)
    chk('*Dynamic, Explicit with a positive step time',
        dyn is not None and float(dyn.group(1)) > 0,
        '%s s' % (dyn.group(1) if dyn else 'missing'))
    t_step = float(dyn.group(1)) if dyn else 0.0

    bc = re.search(r'^\*Boundary, type=VELOCITY\s*\n((?:A_WHEEL_REF[^\n]*\n)+)', t, re.M)
    rows = [r.strip() for r in bc.group(1).strip().split('\n')] if bc else []
    chk('wheel driven through its own rigid-body reference node',
        len(rows) == 6 and all(r.startswith('A_WHEEL_REF') for r in rows),
        '%d BC rows on A_WHEEL_REF' % len(rows))

    def dof(n):
        for r in rows:
            f = [x.strip() for x in r.split(',')]
            if len(f) >= 3 and int(f[1]) == n:
                return float(f[3]) if len(f) > 3 else 0.0
        return None

    v1, v2, vr3 = dof(1), dof(2), dof(6)
    chk('wheel actually rotates', vr3 is not None and abs(vr3) > 0,
        'VR3 = %s rad/s' % vr3)
    # The deck's own header describes the rotation sense. It said "+X-to-+Y travel" for
    # a negative VR3 for months, which is backwards, and nothing compared the sentence
    # with the number underneath it. A positive rotation about +Z carries +X toward +Y,
    # so VR3 < 0 means the surface travels toward DECREASING theta.
    if vr3:
        says_dec = 'decreasing theta' in t.lower()
        says_inc = 'increasing theta' in t.lower()
        chk('the header describes the rotation sense the BC actually applies',
            (vr3 < 0 and says_dec and not says_inc)
            or (vr3 > 0 and says_inc and not says_dec),
            'VR3 = %+g rad/s, header says %s' % (
                vr3, 'decreasing theta' if says_dec else
                     ('increasing theta' if says_inc else 'nothing about theta')))
    # The single most important one: without infeed the wheel spins on the spot,
    # one grit grazes at t=0 and nothing ever cuts. This is what a completed-but-
    # empty run looks like, and it cost a full job to discover.
    chk('wheel is fed INTO the work (there is a depth of cut)',
        v1 is not None and v2 is not None and math.hypot(v1, v2) > 0,
        'V1 = %s, V2 = %s mm/s' % (v1, v2))
    if v1 and v2:
        # Which way the block actually is, taken from its own node coordinates rather
        # than from the report or from an assumed sign. This check used to assert the
        # infeed pointed radially *inward* -- the same mistake the writer was making,
        # so the two agreed with each other and the deck ground nothing for months.
        # The block sits OUTSIDE the rim, so feeding in means moving out along +e_r.
        if HAS_WP:
            wc = np.array(list(parts[WPN]['nodes'].values())).mean(axis=0)[:2]
            er = wc / (np.linalg.norm(wc) or 1.0)
        else:
            th = math.radians(float(REPORT.get('theta_workpiece_deg', 0.0)))
            er = np.array([math.cos(th), math.sin(th)])
        vmag = math.hypot(v1, v2)
        along = (v1 * er[0] + v2 * er[1]) / vmag
        chk('the infeed closes the gap instead of opening it', along > 0.999,
            'infeed . e_r = %+.6f (e_r = %+.4f, %+.4f points from the axis to the '
            'block; +1 means straight at it)' % (along, er[0], er[1]))
        # And the amount it closes by must be the depth of cut that was asked for.
        ae_asked = float((REPORT.get('params', {}).get('analysis') or {})
                         .get('depth_of_cut_um') or 0.0)
        if ae_asked:
            chk('and closes it by exactly the depth of cut',
                abs(vmag * t_step * 1000.0 - ae_asked) < 1e-3,
                '%.4f um travelled vs %.4f um asked' % (vmag * t_step * 1000.0,
                                                        ae_asked))
        ae_um = math.hypot(v1, v2) * t_step * 1000.0
        clr = REPORT.get('max_engaging_protrusion_um')
        if clr:
            chk('depth of cut stays clear of the bond rim', ae_um < clr,
                'ae = %.3f um vs %.3f um clearance, margin %.3f um'
                % (ae_um, clr, clr - ae_um))
        band = REPORT.get('grit_band_arc_mm')
        wpl = (REPORT.get('params') or {}).get('wp_length_mm')
        if band and wpl and vr3:
            sweep = abs(vr3) * t_step * REPORT['outer_radius_mm']
            warn('the block stays on the dressed band for the whole pass',
                 sweep <= (band - wpl) / 2.0 + 1e-12,
                 'sweeps %.4f mm, %.4f mm available before it runs onto bare bond'
                 % (sweep, (band - wpl) / 2.0))

    sc = re.search(r'^\*Section Controls,[^\n]*', t, re.M)
    dv = re.search(r'^\*Depvar, delete=(\d+)\s*\n\s*(\d+)', t, re.M)
    chk('*Solid Section names the section controls',
        bool(sc) and 'controls=EC-1' in t, sc.group(0) if sc else 'missing')
    # All three must line up or damaged material simply never leaves the mesh.
    chk('element deletion armed consistently (section controls + delete= + SDV count)',
        bool(sc) and 'ELEMENT DELETION=YES' in sc.group(0).upper()
        and dv is not None and int(dv.group(1)) <= int(dv.group(2)),
        '*Depvar delete=%s of %s' % (dv.group(1), dv.group(2)) if dv else 'no *Depvar')
    chk('workpiece uses a user material (the JH-2 VUMAT hook)',
        '*User Material' in t,
        re.search(r'^\*User Material[^\n]*', t, re.M).group(0)
        if '*User Material' in t else 'none')

    asurf = set(re.findall(r'^\*Surface, type=ELEMENT, name=([^\s,\n]+)', asm, re.M))
    ci = re.search(r'^\*Contact Inclusions\s*\n([^\n*]+)', t, re.M)
    if ci:
        named = [x.strip() for x in ci.group(1).split(',') if x.strip()]
        chk('contact surfaces exist in the assembly',
            all(n in asurf for n in named), '%s vs %s' % (named, sorted(asurf)))
    else:
        chk('contact defined', 'ALL EXTERIOR' in t.upper(), 'ALL EXTERIOR')
    chk('the workpiece is held', 'A_WP_' in t and 'ENCASTRE' in t.upper(),
        re.search(r'^[A-Z_]*A_WP_[A-Z_]*, ENCASTRE', t, re.M).group(0)
        if 'ENCASTRE' in t.upper() else 'nothing fixed')

    rs = re.search(r'^\*Restart, write, number interval=(\d+)', t, re.M)
    # interval=1 writes the only restart state at the END of the step, so an
    # interrupted run cannot be resumed at all.
    warn('the run can be resumed after an interruption',
         rs is not None and int(rs.group(1)) > 1,
         'restart interval = %s' % (rs.group(1) if rs else 'none'))
    fo = re.search(r'^\*Output, field, number interval=(\d+)', t, re.M)
    chk('field output requested', fo is not None,
        '%s frames' % (fo.group(1) if fo else 'none'))
    chk('SDV requested, so JH-2 damage can be plotted', 'SDV' in t)

print()
print('=' * 78)
print('TOTAL: %d failure(s), %d warning(s)%s%s'
      % (len(FAIL), len(WARN), '' if not FAIL else '  -> ' + str(FAIL),
         '' if not WARN else '  warn -> ' + str(WARN)))
print('=' * 78)
sys.exit(1 if FAIL else 0)
