"""Standalone verifier for a semgrit wheel+workpiece .inp deck.

    python verify_inp_geometry.py <deck.inp> [outer_radius_mm]

Four passes, all reading the finished file with an independent parser:

  P1  mesh      element types, facet edge lengths, hex Jacobians, rigid ref nodes
  P2  syntax    block balance, line format, every set/surface reference resolves
  P3  geometry   grit/workpiece tangency, true-arc test, grit window vs bond arc
  P4  regression every failure mode previously hit, re-checked by name

Exits non-zero if anything fails.
"""
import math
import re
import sys
from collections import Counter

import numpy as np

PATH = sys.argv[1] if len(sys.argv) > 1 else 'FINAL_MODEL/wheel_and_workpiece.inp'
R_NOM = float(sys.argv[2]) if len(sys.argv) > 2 else None

t = open(PATH, encoding='ascii', errors='replace').read()
lines = t.split('\n')
names = [l.split(',')[0].strip().lower() for l in lines
         if l.startswith('*') and not l.startswith('**')]
FAIL = []


def chk(n, ok, d=''):
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', n, (': ' + d) if d else ''))
    if not ok:
        FAIL.append(n)


parts = {}
for pm in re.finditer(r'^\*Part,\s*name=([^\s,]+)(.*?)^\*End Part', t, re.M | re.S):
    nm, b = pm.group(1), pm.group(2)
    nd = {}
    for nq in re.finditer(r'^\*Node\s*$(.*?)(?=^\*)', b, re.M | re.S):
        for ln in nq.group(1).strip().split('\n'):
            v = [x.strip() for x in ln.split(',')]
            if len(v) >= 4:
                nd[int(v[0])] = np.array([float(x) for x in v[1:4]])
    els = {}
    for em in re.finditer(r'^\*Element,\s*type=(\w+)\s*$(.*?)(?=^\*)', b, re.M | re.S):
        for ln in em.group(2).strip().split('\n'):
            v = [x.strip() for x in ln.split(',')]
            if len(v) >= 4:
                els.setdefault(em.group(1), {})[int(v[0])] = [int(x) for x in v[1:]]
    es = set(x.upper() for x in re.findall(r'\*Elset,\s*elset=([^\s,\n]+)', b, re.I))
    ns = set(x.upper() for x in re.findall(r'\*Nset,\s*nset=([^\s,\n]+)', b, re.I))
    sf = {}
    for m in re.finditer(r'\*Surface,\s*type=ELEMENT,\s*name=([^\s,\n]+)\s*\n([^\n]+)',
                         b, re.I):
        sf[m.group(1).upper()] = [x.strip().upper() for x in m.group(2).split(',')]
    parts[nm] = dict(nd=nd, els=els, es=es, ns=ns, sf=sf)

print('=' * 74)
print('VERIFYING %s' % PATH)
print('=' * 74)
print('P1  MESH')

types = set()
for d in parts.values():
    types |= set(d['els'])
chk('Explicit element types only', types == {'C3D8R', 'R3D3'}, str(sorted(types)))

worst = 1e30
nfac = 0
for d in parts.values():
    for c in d['els'].get('R3D3', {}).values():
        p = [d['nd'][i] for i in c]
        n = len(p)
        nfac += 1
        worst = min(worst, min(np.linalg.norm(p[i] - p[(i + 1) % n]) for i in range(n)))
chk('no coincident-node facets', worst >= 1e-5,
    '%d facets, min edge %.4e mm' % (nfac, worst))


def hexvol(nd, c):
    p = [nd[i] for i in c]
    v = 0.0
    for a, b, cc, dd in ((0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6),
                         (3, 4, 6, 7), (1, 4, 5, 6)):
        v += np.dot(np.cross(p[b] - p[a], p[cc] - p[a]), p[dd] - p[a]) / 6.0
    return v


for nm in ('WHEEL_BOND', 'WORKPIECE'):
    if nm not in parts:
        continue
    vv = np.array([hexvol(parts[nm]['nd'], c)
                   for c in parts[nm]['els']['C3D8R'].values()])
    chk('%s hex volumes positive' % nm, bool((vv > 0).all()),
        '%d elems, min %.3e mm3' % (len(vv), vv.min()))

bad = 0
ngp = 0
for nm, d in parts.items():
    if not nm.startswith('GRAIN-'):
        continue
    ngp += 1
    used = set(n for c in d['els'].get('R3D3', {}).values() for n in c)
    if len(set(d['nd']) - used) != 1:
        bad += 1
chk('every grit part has a rigid-body ref node', bad == 0,
    '%d grit parts, %d bad' % (ngp, bad))

print()
print('P2  SYNTAX / REFERENCES')
c = Counter(names)
chk('blocks balanced',
    c['*part'] == c['*end part'] and c['*instance'] == c['*end instance'],
    'part %d/%d instance %d/%d' % (c['*part'], c['*end part'],
                                   c['*instance'], c['*end instance']))
chk('line length <=256 and pure ASCII',
    max(len(l) for l in lines) <= 256 and all(ord(ch) < 128 for l in lines for ch in l),
    '%d chars' % max(len(l) for l in lines))
bsf = [(nm, s, r[0]) for nm, d in parts.items() for s, r in d['sf'].items()
       if r[0] not in d['es']]
chk('part *Surface -> *Elset in same part', not bsf, str(bsf[:3]))

asm = t.split('*Assembly')[1].split('*End Assembly')[0]
inst = {}
for m in re.finditer(r'\*Instance,\s*name=([^\s,]+),\s*part=([^\s,]+)', asm):
    inst[m.group(1).upper()] = m.group(2)
ba = []
for sm in re.finditer(r'\*Surface,\s*type=ELEMENT,\s*name=([^\s,\n]+)\s*\n([^\n]+)',
                      asm, re.I):
    pr = [x.strip() for x in sm.group(2).split(',') if x.strip()]
    if '.' not in pr[0]:
        ba.append((sm.group(1), 'not instance-qualified'))
        continue
    i, s = pr[0].split('.', 1)
    if i.upper() not in inst:
        ba.append((sm.group(1), 'unknown instance'))
    elif s.upper() not in parts[inst[i.upper()]]['es']:
        ba.append((sm.group(1), 'elset missing'))
    elif len(pr) < 2:
        ba.append((sm.group(1), 'no face identifier'))
chk('assembly *Surface -> instance.ELSET + face id', not ba, str(ba[:3]))

bn = []
for m in re.finditer(r'\*Nset,\s*nset=([^\s,\n]+),\s*instance=([^\s,\n]+)\s*\n([^\n]+)',
                     asm, re.I):
    i = m.group(2).upper()
    if i not in inst or m.group(3).split(',')[0].strip().upper() not in parts[inst[i]]['ns']:
        bn.append(m.group(1))
chk('assembly *Nset,instance= resolves', not bn, str(bn[:3]))

mats = set(x.upper() for x in re.findall(r'^\*Material,\s*name=([^\s,\n]+)', t, re.M))
usedm = set(x.upper() for x in re.findall(r'\*Solid Section,[^\n]*material=([^\s,\n]+)', t))
chk('section materials defined', not (usedm - mats), str(sorted(usedm)))

print()
print('P3  GEOMETRY / ARC')


def rotm(a, d):
    n = np.linalg.norm(a)
    if n < 1e-15:
        return np.eye(3)
    x, y, z = a / n
    cs, sn = math.cos(math.radians(d)), math.sin(math.radians(d))
    return np.array([[cs + x * x * (1 - cs), x * y * (1 - cs) - z * sn, x * z * (1 - cs) + y * sn],
                     [y * x * (1 - cs) + z * sn, cs + y * y * (1 - cs), y * z * (1 - cs) - x * sn],
                     [z * x * (1 - cs) - y * sn, z * y * (1 - cs) + x * sn, cs + z * z * (1 - cs)]])


W = np.array([parts['WORKPIECE']['nd'][i] for i in sorted(parts['WORKPIECE']['nd'])])
r_ground = float(np.hypot(W[:, 0], W[:, 1]).min())
tips, gth = [], []
for m in re.finditer(r'\*Instance,\s*name=(G-\d+),\s*part=([^\s,]+)\s*\n([^*]*?)\*End Instance', t):
    nd = parts[m.group(2)]['nd']
    X = np.array([nd[i] for i in sorted(nd)])
    r = [[float(x) for x in ln.split(',') if x.strip()]
         for ln in m.group(3).strip().split('\n') if ln.strip()]
    Xw = X + np.array(r[0][:3])
    if len(r) > 1 and len(r[1]) >= 7:
        a = np.array(r[1][0:3]); b = np.array(r[1][3:6])
        Xw = (Xw - a) @ rotm(b - a, r[1][6]).T + a
    tips.append(float(np.hypot(Xw[:, 0], Xw[:, 1]).max()))
    gth.append(math.degrees(math.atan2(Xw[:, 1].mean(), Xw[:, 0].mean())))
tips = np.array(tips); gth = np.array(gth)

chk('tallest grit tangent to the ground face', abs(r_ground - tips.max()) < 1e-6,
    '%+.4f nm' % ((r_ground - tips.max()) * 1e6))
chk('no grit penetrates the workpiece', tips.max() <= r_ground + 1e-6,
    '%d past the surface' % int((tips > r_ground + 1e-9).sum()))

bnd = parts['WHEEL_BOND']['nd']
rb = np.array([float(np.hypot(bnd[i][0], bnd[i][1])) for i in sorted(bnd)])
nrad = len(set(round(x, 6) for x in rb))
R = R_NOM if R_NOM else float(rb.max())
chk('bond nodes lie on exact circles (true arc)', nrad <= 4,
    '%d distinct radii, r = %.4f .. %.4f mm' % (nrad, rb.min(), rb.max()))
th = np.array([math.degrees(math.atan2(bnd[i][1], bnd[i][0])) for i in sorted(bnd)])
span = th.max() - th.min()
L = R * math.radians(span)
rim = rb.max() - rb.min()
sag = L * L / (8 * R)
print('      bond: R=%.3f mm, %.4f deg, arc %.4f mm, rim %.4f mm' % (R, span, L, rim))
print('      sagitta %.2f um = %.0f%% of the rim depth  -> %s'
      % (sag * 1000, 100 * sag / rim, 'visibly curved' if sag / rim > 0.25 else 'renders flat'))
chk('grits lie inside the bond arc',
    th.min() - 0.02 <= gth.min() and gth.max() <= th.max() + 0.02,
    'grits %.4f..%.4f deg, bond %.4f..%.4f deg' % (gth.min(), gth.max(), th.min(), th.max()))
wth = math.degrees(math.atan2(W[:, 1].mean(), W[:, 0].mean()))
chk('workpiece centred on the grit window',
    abs(wth - 0.5 * (gth.min() + gth.max())) < 0.05,
    'workpiece %.4f deg vs grit-window centre %.4f deg'
    % (wth, 0.5 * (gth.min() + gth.max())))

print()
print('P4  REGRESSION (every failure mode previously hit)')
forb = set(names) & {'*step', '*boundary', '*contact', '*dynamic', '*static',
                     '*shear failure', '*plastic', '*surface interaction', '*output'}
chk('facet coincident nodes', worst >= 1e-5, 'min edge %.4e mm' % worst)
chk('assembly surface referencing a part surface', not ba, 'all resolve')
chk('*SHEARFAILURE / step / contact ordering', not forb,
    str(sorted(forb)) if forb else 'none present')
chk('Standard-library elements', types == {'C3D8R', 'R3D3'})
chk('instances present (not 0)', len(inst) > 0, '%d instances' % len(inst))

print()
print('=' * 74)
print('TOTAL: %d failure(s)%s' % (len(FAIL), '' if not FAIL else ' -> ' + str(FAIL)))
print('=' * 74)
sys.exit(1 if FAIL else 0)
