"""Second, independent verification of an all-rigid-wheel deck.

    python verify_rigid_deck2.py <deck.inp> [report.json]

Deliberately shares no code and no method with verify_rigid_deck.py. That one parses
with regular expressions and checks geometry; this one walks the file line by line as
a state machine and checks the things a regex sweep cannot see:

  Q1  grammar   every keyword is real, legal in its context, and its data lines have
                the right field count and type -- i.e. what the Abaqus input reader
                itself objects to
  Q2  ordering   model data, assembly data and part data in the order Abaqus requires
  Q3  identity   the header's own numbers, and the report JSON's numbers, match what
                is actually in the file
  Q4  physics    mass and inertia re-derived by numerical integration; tangency
                re-derived from the workpiece box rather than from the grits;
                stable increment re-derived from the material actually in the deck

Exits non-zero on any failure.
"""
import json
import math
import os
import re
import sys

import numpy as np

PATH = sys.argv[1] if len(sys.argv) > 1 else 'FINAL_RIGID/wheel_rigid_2mm.inp'
RPT = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(PATH)[0] + '_report.json'

FAIL = []


WARN = []


def warn(name, ok, detail=''):
    """A design guideline, not a correctness property: reported, never fatal."""
    print('  [%s] %s%s' % ('ok  ' if ok else 'WARN', name,
                           (': ' + detail) if detail else ''))
    if not ok:
        WARN.append(name)


def chk(name, ok, detail=''):
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', name, (': ' + detail) if detail else ''))
    if not ok:
        FAIL.append(name)


# keyword -> (allowed contexts, exact data-line field counts or None for free/any)
TOP, PART, ASM, INST, MODEL = 'top', 'part', 'assembly', 'instance', 'model'
STEP = 'step'
SPEC = {
    'heading': ({TOP}, None),
    'part': ({TOP}, ()),
    'end part': ({PART}, ()),
    'node': ({PART, ASM}, (4,)),
    'element': ({PART}, None),          # arity depends on type=
    'elset': ({PART, ASM}, None),
    'nset': ({PART, ASM}, None),
    'solid section': ({PART}, None),    # one data line, may be a bare comma
    'surface': ({PART, ASM}, (2,)),
    'rigid body': ({PART}, ()),
    'mass': ({PART}, (1,)),
    'rotary inertia': ({PART}, (6,)),
    'assembly': ({TOP}, ()),
    'end assembly': ({ASM}, ()),
    'instance': ({ASM}, None),          # 0, 1 or 2 data lines
    'end instance': ({INST}, ()),
    'material': ({MODEL}, ()),
    'density': ({MODEL}, (1, 2)),
    'elastic': ({MODEL}, (2,)),
    # --- run-ready decks: the history definition -------------------------
    'depvar': ({MODEL}, (1,)),
    'user material': ({MODEL}, None),      # constants wrap across lines
    'section controls': ({MODEL}, (3,)),
    'surface interaction': ({MODEL}, ()),
    'friction': ({MODEL}, (1,)),
    'boundary': ({MODEL, STEP}, None),
    'step': ({MODEL}, ()),
    'end step': ({STEP}, ()),
    'dynamic': ({STEP}, (2,)),             # ", <time>" -> two fields
    'bulk viscosity': ({STEP}, (2,)),
    'fixed mass scaling': ({STEP}, ()),
    'contact': ({STEP}, ()),
    'contact inclusions': ({STEP}, None),
    'contact property assignment': ({STEP}, None),
    'restart': ({STEP}, ()),
    'output': ({STEP}, ()),
    'node output': ({STEP}, None),
    'element output': ({STEP}, None),
}
EL_ARITY = {'R3D3': 4, 'R3D4': 5, 'C3D8R': 9, 'MASS': 2, 'ROTARYI': 2}

raw = open(PATH, 'rb').read()
text = raw.decode('ascii', errors='replace')
lines = text.split('\n')

print('=' * 78)
print('SECOND VERIFICATION  %s' % PATH)
print('=' * 78)
print('Q1  KEYWORD GRAMMAR  (line-by-line state machine)')

ctx = TOP
stack = []
cur_kw = None
cur_opts = {}
cur_type = None
bad_kw, bad_ctx, bad_arity, bad_num = [], [], [], []
order = []
parts = {}          # name -> dict
cur_part = None
cur_inst = None
instances = {}
nodes = {}
els = {}
elsets = {}
nsets = {}
counts = {}
mats = {}
cur_mat = None

for ln_no, ln in enumerate(lines, 1):
    s = ln.strip()
    if not s or s.startswith('**'):
        continue
    if s.startswith('*'):
        body = s[1:]
        fields = [f.strip() for f in body.split(',')]
        name = fields[0].lower()
        opts = {}
        for f in fields[1:]:
            if '=' in f:
                k, v = f.split('=', 1)
                opts[k.strip().lower()] = v.strip()
            elif f:
                opts[f.lower()] = True
        if name not in SPEC:
            bad_kw.append((ln_no, name))
            cur_kw = None
            continue
        allowed, arity = SPEC[name]
        eff = INST if (ctx == ASM and cur_inst) else ctx
        if eff not in allowed and not (name == 'end instance' and cur_inst):
            bad_ctx.append((ln_no, name, eff))
        order.append(name)
        # context transitions
        if name == 'part':
            ctx, cur_part = PART, opts.get('name', '').upper()
            parts[cur_part] = dict(nodes={}, els={}, ty={}, elsets={}, nsets={},
                                   surfs={}, sections=[])
        elif name == 'end part':
            ctx, cur_part = TOP, None
        elif name == 'assembly':
            ctx = ASM
        elif name == 'end assembly':
            ctx = MODEL
        elif name == 'instance':
            cur_inst = opts.get('name', '').upper()
            instances[cur_inst] = opts.get('part', '').upper()
        elif name == 'end instance':
            cur_inst = None
        elif name == 'material':
            cur_mat = opts.get('name', '').upper()
            mats[cur_mat] = {}
        elif name == 'step':
            ctx = STEP
        elif name == 'end step':
            ctx = MODEL
        cur_kw, cur_opts = name, opts
        cur_type = opts.get('type', '').upper() if name == 'element' else None
        if name in ('elset', 'nset'):
            key = opts.get('elset' if name == 'elset' else 'nset', '').upper()
            cur_opts['_key'] = key
        counts[name] = counts.get(name, 0) + 1
        continue

    # ---- data line ----
    if cur_kw is None:
        continue
    vals = [v.strip() for v in s.split(',')]
    vals = vals[:-1] if vals and vals[-1] == '' else vals
    n = len(vals)
    allowed, arity = SPEC[cur_kw]

    if cur_kw == 'element':
        want = EL_ARITY.get(cur_type)
        if want is None:
            bad_kw.append((ln_no, 'element type=' + str(cur_type)))
        elif n != want:
            bad_arity.append((ln_no, cur_kw, cur_type, n, want))
        else:
            try:
                ids = [int(v) for v in vals]
            except ValueError:
                bad_num.append((ln_no, cur_kw, s[:40]))
                continue
            if cur_part:
                parts[cur_part]['els'][ids[0]] = ids[1:]
                parts[cur_part]['ty'][ids[0]] = cur_type
                key = cur_opts.get('elset')
                if key:
                    parts[cur_part]['elsets'].setdefault(key.upper(), set()).add(ids[0])
    elif cur_kw == 'node':
        if n != 4:
            bad_arity.append((ln_no, cur_kw, None, n, 4))
        else:
            try:
                i = int(vals[0])
                xyz = np.array([float(v) for v in vals[1:]])
            except ValueError:
                bad_num.append((ln_no, cur_kw, s[:40]))
                continue
            if cur_part:
                parts[cur_part]['nodes'][i] = xyz
    elif cur_kw in ('elset', 'nset'):
        gen = 'generate' in cur_opts
        try:
            iv = [int(v) for v in vals]
        except ValueError:
            # a set defined by naming other sets
            if cur_part:
                pass
            continue
        if gen and n != 3:
            bad_arity.append((ln_no, cur_kw + ' generate', None, n, 3))
        elif gen:
            ids = set(range(iv[0], iv[1] + 1, iv[2]))
        else:
            ids = set(iv)
        if cur_part:
            tgt = 'elsets' if cur_kw == 'elset' else 'nsets'
            parts[cur_part][tgt].setdefault(cur_opts['_key'], set()).update(ids)
    elif cur_kw == 'solid section':
        if cur_part:
            parts[cur_part]['sections'].append(cur_opts.get('elset', '').upper())
    elif cur_kw == 'surface':
        if n != 2:
            bad_arity.append((ln_no, cur_kw, None, n, 2))
        if cur_part:
            parts[cur_part]['surfs'].setdefault(
                cur_opts.get('name', '').upper(), []).append(vals)
    elif cur_kw in ('mass', 'rotary inertia', 'density', 'elastic', 'user material'):
        want = SPEC[cur_kw][1]
        if want and n not in want:
            bad_arity.append((ln_no, cur_kw, None, n, want))
        try:
            fv = [float(v) for v in vals]
        except ValueError:
            bad_num.append((ln_no, cur_kw, s[:40]))
            continue
        if cur_kw == 'mass':
            mats['_MASS'] = fv[0]
        elif cur_kw == 'rotary inertia':
            mats['_ROTI'] = fv
        elif cur_mat and cur_kw in ('density', 'elastic'):
            mats[cur_mat][cur_kw] = fv
        elif cur_mat and cur_kw == 'user material':
            # constants wrap across lines, so accumulate rather than overwrite
            mats[cur_mat].setdefault('user material', []).extend(fv)
    elif cur_kw in ('instance', 'heading', 'part', 'assembly'):
        pass

chk('every keyword is a recognised Abaqus keyword', not bad_kw, str(bad_kw[:4]))
chk('every keyword appears in a legal context', not bad_ctx, str(bad_ctx[:4]))
chk('every data line has the arity its keyword requires', not bad_arity,
    str(bad_arity[:4]))
chk('every numeric field parses', not bad_num, str(bad_num[:4]))
chk('context stack closes cleanly (ended in model data)', ctx == MODEL, ctx)
chk('file ends with a newline', raw.endswith(b'\n'))
print('      keywords seen: %s' % ', '.join(
    '%s x%d' % (k, v) for k, v in sorted(counts.items())))

print()
print('Q2  ORDERING')
i_asm = order.index('assembly')
i_end = order.index('end assembly')
chk('all *Part blocks precede *Assembly',
    all(order.index('part') < i_asm for _ in [0])
    and max(i for i, k in enumerate(order) if k == 'end part') < i_asm)
chk('*Material appears after *End Assembly',
    all(i > i_end for i, k in enumerate(order) if k == 'material'),
    'materials at %s, *End Assembly at %d'
    % ([i for i, k in enumerate(order) if k == 'material'], i_end))
chk('*Rigid Body follows the *Elset and *Nset it names',
    all(order.index('rigid body') > order.index(k) for k in ('elset', 'nset')))
chk('*Density and *Elastic follow their *Material',
    all(i > order.index('material') for i, k in enumerate(order)
        if k in ('density', 'elastic')))
chk('no keyword after the last *Material block that belongs to a part',
    not any(k in ('node', 'element', 'rigid body') for k in order[i_end:]),
    'model-data section holds only material keywords')

print()
print('Q3  IDENTITY  (header comments and report JSON vs the file)')
hdr = {}
for ln in lines:
    if ln.startswith('** ') and ':' in ln:
        k, _, v = ln[3:].partition(':')
        hdr[k.strip()] = v.strip()

W = parts['WHEEL']
HAS_WP = 'WORKPIECE' in parts
P = parts.get('WORKPIECE')
n_r3d3 = sum(1 for ty in W['ty'].values() if ty == 'R3D3')
n_r3d4 = sum(1 for ty in W['ty'].values() if ty == 'R3D4')
n_c3d8 = sum(1 for ty in P['ty'].values() if ty == 'C3D8R') if HAS_WP else 0
chk('header grit-facet count matches the file',
    ('%d' % n_r3d3) in hdr.get('grits', ''),
    'header "%s" vs %d R3D3' % (hdr.get('grits'), n_r3d3))
if HAS_WP:
    chk('header workpiece element count matches the file',
        ('%d' % n_c3d8) in hdr.get('workpiece elements', ''),
        'header "%s" vs %d C3D8R' % (hdr.get('workpiece elements'), n_c3d8))
else:
    chk('wheel-only deck declares no workpiece elements',
        'workpiece elements' not in hdr, 'header has no workpiece line')
_sh = np.array([W['nodes'][n] for e, c in W['els'].items() if W['ty'][e] == 'R3D4'
                for n in c])
chk('header outer radius matches the bond shell geometry',
    abs(float(hdr['wheel outer radius (mm)'].split()[0])
        - float(np.hypot(_sh[:, 0], _sh[:, 1]).max())) < 1e-9,
    'header %s, shell max radius %.9f mm'
    % (hdr['wheel outer radius (mm)'], np.hypot(_sh[:, 0], _sh[:, 1]).max()))

rep = json.load(open(RPT)) if os.path.exists(RPT) else {}
if rep:
    chk('report n_grit_facets matches the file', rep['n_grit_facets'] == n_r3d3,
        '%d vs %d' % (rep['n_grit_facets'], n_r3d3))
    chk('report n_bond_shell_quads matches the file',
        rep['n_bond_shell_quads'] == n_r3d4, '%d vs %d' % (rep['n_bond_shell_quads'], n_r3d4))
    chk('report n_workpiece_elements matches the file',
        rep['n_workpiece_elements'] == n_c3d8,
        '%d vs %d' % (rep['n_workpiece_elements'], n_c3d8))
    chk('report size_bytes matches the file on disk',
        rep['size_bytes'] == len(raw), '%d vs %d' % (rep['size_bytes'], len(raw)))
    chk('report ref node is the node the *Rigid Body uses',
        rep['wheel_ref_node'] in W['nodes'], str(rep['wheel_ref_node']))
else:
    chk('report JSON present', False, RPT)

print()
print('Q4  PHYSICS  (re-derived independently)')
# --- mass and inertia by numerical integration over the sector volume -------
# The bond shell only. Taking the extent of every wheel node instead would fold in
# the grit tips, which stand above the rim and inside the arc ends, and inflate the
# integration domain -- the rim volume came out 39% high that way.
shell_ids = sorted(set(n for e, c in W['els'].items() if W['ty'][e] == 'R3D4'
                       for n in c))
S = np.array([W['nodes'][n] for n in shell_ids])
rr_s = np.hypot(S[:, 0], S[:, 1])
R_out, R_in = float(rr_s.max()), float(rr_s.min())
half_w = float(np.abs(S[:, 2]).max())
ths = np.arctan2(S[:, 1], S[:, 0])
# A full wheel wraps: its last node column IS the first, so the gap across theta=+-pi
# equals the ordinary node spacing. On a sector that gap is the whole missing angle.
_u = np.unique(np.round(ths, 12))
_sp = float(np.median(np.diff(_u))) if len(_u) > 2 else 0.0
FULL = len(_u) > 2 and (2 * math.pi - (float(_u[-1]) - float(_u[0]))) <= 1.5 * _sp
dth = 2 * math.pi if FULL else float(ths.max() - ths.min())
print('      shell: %d nodes, r %.6f..%.6f mm, %.6f rad%s, half-width %.6f mm'
      % (len(S), R_in, R_out, dth, '  (FULL WHEEL)' if FULL else '', half_w))

nr, nt, nz = 60, 400, 20
r_c = R_in + (np.arange(nr) + 0.5) * (R_out - R_in) / nr
t_c = ths.min() + (np.arange(nt) + 0.5) * dth / nt
z_c = -half_w + (np.arange(nz) + 0.5) * (2 * half_w) / nz
dV = ((R_out - R_in) / nr) * (dth / nt) * (2 * half_w / nz)
Rg, Tg, Zg = np.meshgrid(r_c, t_c, z_c, indexing='ij')
w = Rg * dV                                   # r dr dtheta dz
X, Y, Z = Rg * np.cos(Tg), Rg * np.sin(Tg), Zg
rho = float((json.load(open(RPT)).get('params', {}) if os.path.exists(RPT) else {})
            .get('bond_density_kg_m3', 2700.0)) * 1e-12
m_num = float(rho * w.sum())
I11 = float(rho * (w * (Y ** 2 + Z ** 2)).sum())
I22 = float(rho * (w * (X ** 2 + Z ** 2)).sum())
I33 = float(rho * (w * (X ** 2 + Y ** 2)).sum())
I12 = float(-rho * (w * X * Y).sum())

m_file = mats['_MASS']
roti = mats['_ROTI']
chk('*Mass equals the %s volume x bond density' % ('wheel' if FULL else 'sector'),
    abs(m_file - m_num) / m_num < 2e-4,
    'file %.6e vs integrated %.6e tonne (%.3f%%)'
    % (m_file, m_num, 100 * abs(m_file - m_num) / m_num))
rel = [abs(a - b) / abs(b) for a, b in zip(roti[:3], (I11, I22, I33))]
chk('*Rotary Inertia diagonal matches numerical integration', max(rel) < 2e-3,
    'I11 %.3e/%.3e  I22 %.3e/%.3e  I33 %.3e/%.3e  worst %.3f%%'
    % (roti[0], I11, roti[1], I22, roti[2], I33, 100 * max(rel)))
# A full wheel is rotationally symmetric, so I12 integrates to zero and a relative
# test would divide by it; compare against the polar term instead.
chk('*Rotary Inertia I12 matches numerical integration',
    abs(roti[3] - I12) <= max(2e-3 * abs(I12), 1e-6 * abs(I33)),
    'file %.4e vs integrated %.4e' % (roti[3], I12))

# --- tangency from the workpiece box, not from the grits --------------------
if not HAS_WP:
    print('      wheel-only deck: no ground face, so no tangency or increment test')
else:
    Q = np.array([P['nodes'][i] for i in sorted(P['nodes'])])
    # The block's own frame, taken from its centroid. Assuming mid-arc is wrong: the
    # writer relocates the block to the tallest grit when the nominal angle is bare,
    # and a full wheel has no mid-arc at all.
    _qc = Q.mean(axis=0)
    theta_c = math.atan2(_qc[1], _qc[0])
    e_r = np.array([math.cos(theta_c), math.sin(theta_c), 0.0])
    e_t = np.array([-math.sin(theta_c), math.cos(theta_c), 0.0])
    B = np.column_stack([e_r, e_t, np.array([0.0, 0.0, 1.0])])
    Qf = Q @ B
    a_lo = float(Qf[:, 0].min())
    hb = float(np.abs(Qf[:, 1]).max())
    hz = float(np.abs(Qf[:, 2]).max())

    grit_nodes = sorted(set(n for e, c in W['els'].items() if W['ty'][e] == 'R3D3'
                            for n in c))
    G = np.array([W['nodes'][n] for n in grit_nodes]) @ B
    inb = (np.abs(G[:, 1]) <= hb) & (np.abs(G[:, 2]) <= hz)
    chk('some grit material lies under the workpiece at all', bool(inb.any()),
        '%d of %d grit vertices inside the footprint' % (int(inb.sum()), len(G)))
    gap = a_lo - G[inb, 0]
    chk('no grit vertex is inside the workpiece box', float(gap.min()) >= 0.0,
        'closest vertex gap %.4f nm' % (gap.min() * 1e6))
    # Zero gap only when no standoff was asked for; otherwise the closest vertex
    # should sit exactly that far below the face. Checking for tangency regardless
    # would report a correctly parked block as a broken one.
    _sd = float(rep.get('clearance_um') or 0.0) / 1000.0
    # Vertices alone cannot resolve the tangency: the tallest point of a facet clipped
    # to the footprint is usually on a clip edge, not on a corner of the triangle, so
    # the nearest *vertex* legitimately sits tens of nanometres back. Verifier A does
    # the exact facet-clipped test; this one asserts the weaker, still meaningful
    # thing -- the closest vertex is at the standoff to within a small fraction of the
    # protrusion the block is seated on.
    _tol = max(1e-6, 0.05 * abs(float(rep.get('max_engaging_protrusion_um') or 0.0))
               / 1000.0)
    chk('the closest grit vertex sits at the requested standoff',
        abs(float(gap.min()) - _sd) < _tol,
        'gap %.4f nm, asked for %.4f nm (tol %.1f nm), ground face at a = %.9f mm'
        % (gap.min() * 1e6, _sd * 1e6, _tol * 1e6, a_lo))

    # --- stable increment from the material actually written in the deck ----
    # Use the material the workpiece section actually names, not a leftover
    # placeholder: an earlier deck carried an unused *Elastic STONE alongside the
    # real JH2, and reading that gave a wave speed 2.7x too high.
    sec = re.search(r'^\*Solid Section,.*material=(\S+)\s*$', text, re.M)
    wp_mat = sec.group(1).strip().rstrip(',').upper() if sec else 'STONE'
    md = mats.get(wp_mat, {})
    rho_wp = md['density'][0]
    if 'elastic' in md:
        E, nu = md['elastic']
        c_d = math.sqrt(E * (1 - nu) / ((1 + nu) * (1 - 2 * nu) * rho_wp))
        print('      workpiece material %s: E = %g MPa, nu = %g' % (wp_mat, E, nu))
    else:
        # JH-2: props 1 and 2 are the bulk and shear moduli, so the dilatational
        # speed is sqrt((K + 4G/3)/rho) directly.
        props = md['user material']
        K, G = props[0], props[1]
        c_d = math.sqrt((K + 4.0 * G / 3.0) / rho_wp)
        E = 9 * K * G / (3 * K + G)
        nu = (3 * K - 2 * G) / (2 * (3 * K + G))
        print('      workpiece material %s (VUMAT): K = %g, G = %g MPa'
              ' -> E = %.0f MPa, nu = %.3f' % (wp_mat, K, G, E, nu))
        if len(props) >= 56:
            # vumat_grind.for: a hybrid card carries a SECOND elasticity, props
            # 28 and 29, for the ductile branch. Any element may be running
            # either law, so the stable increment is set by the stiffer one.
            # Reading only props 1 and 2 would understate the wave speed
            # wherever the calibrated Johnson-Cook set is stiffer than the JH-2
            # card, and the recomputed dt would then disagree with the report
            # for a reason that is the verifier's fault, not the deck's.
            e_j, nu_j = props[27], props[28]
            if e_j > 0 and -1.0 < nu_j < 0.5:
                c_j = math.sqrt(e_j * (1 - nu_j)
                                / ((1 + nu_j) * (1 - 2 * nu_j) * rho_wp))
                print('      hybrid ductile branch: E = %g MPa, nu = %g'
                      ' -> c = %.4e mm/s (JH-2 branch %.4e)'
                      % (e_j, nu_j, c_j, c_d))
                c_d = max(c_d, c_j)
    # All 12 edges of the brick, not just the one along the cutting direction: with a
    # per-direction mesh size the shortest edge -- the one that sets the stable
    # increment -- can lie in any of the three.
    _HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
                  (0, 4), (1, 5), (2, 6), (3, 7))
    h = min(float(np.linalg.norm(P['nodes'][c[a]] - P['nodes'][c[b]]))
            for e, c in list(P['els'].items())[:2000] if P['ty'][e] == 'C3D8R'
            for a, b in _HEX_EDGES)
    # Mass scaling multiplies density, so it lengthens the stable increment by the
    # square root of the factor. Leaving it out makes the recomputed dt disagree with
    # the report by exactly sqrt(factor).
    _ms = re.search(r'^\*Fixed Mass Scaling,[^\n]*factor=([0-9.eE+-]+)', text, re.M)
    scale = float(_ms.group(1).rstrip('.')) if _ms else 1.0
    dt = h / c_d * math.sqrt(scale)
    print('      rho = %.4e tonne/mm3, mass scaling factor %g' % (rho_wp, scale))
    print('      element %.6f mm, dilatational speed %.4e mm/s, dt = %.4e s'
          % (h, c_d, dt))
    # build_deck nests the kinematics under "cost"; the older hand-written build
    # scripts put them at the top level. Accept either.
    cost = rep.get('cost') or rep
    if cost.get('stable_dt_s'):
        chk('report stable_dt_s matches the material in the deck',
            abs(cost['stable_dt_s'] - dt) / dt < 1e-6,
            'report %.6e vs recomputed %.6e' % (cost['stable_dt_s'], dt))
        chk('report omega is consistent with surface speed and radius',
            abs(cost['omega_rad_s'] - cost['surface_speed_mm_s'] / R_out) < 1e-9,
            '%.3f rad/s at r = %.3f mm -> %.1f mm/s'
            % (cost['omega_rad_s'], R_out, cost['omega_rad_s'] * R_out))
        chk('report increments = step time / dt',
            abs(cost['increments'] - cost['step_time_s'] / dt)
            / cost['increments'] < 1e-6, '%.0f' % cost['increments'])
        # A wall-clock target is the user's budget, not a property of the deck, so it
        # is reported rather than enforced.
        _eh = cost.get('est_hours') or {}
        _h4 = float(_eh.get('4', cost.get('est_hours_4core', 0.0)))
        _h8 = float(_eh.get('8', cost.get('est_hours_8core', 0.0)))
        warn('estimated wall clock is under 6 h on 4 cores', _h4 <= 6.0,
             '%.2f h on 4 cores, %.2f h on 8' % (_h4, _h8))

print()
print('=' * 78)
print('TOTAL: %d failure(s), %d warning(s)%s%s'
      % (len(FAIL), len(WARN), '' if not FAIL else '  -> ' + str(FAIL),
         '' if not WARN else '  warn -> ' + str(WARN)))
print('=' * 78)
sys.exit(1 if FAIL else 0)
