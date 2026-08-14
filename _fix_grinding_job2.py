"""Conservative fix: same three faults, but using deletions only.

The first attempt also *added* point-mass and rotary-inertia elements so that ALLKE
would report the wheel. That is the only part of the patch that introduces new
keywords, so it is the only part that can plausibly upset the input file processor --
and it is not needed for the wheel to turn. With all six degrees of freedom prescribed
on the reference node, a rigid body needs no mass at all.

So this version does the minimum that makes the wheel move, and nothing else:

  * delete the kinematic coupling that was silently discarded and left the wheel undriven
  * delete the now-unused reference point, its two node sets, and the orphaned node
    surface the coupling used -- so nothing dangles and no massless free node is left
  * point the existing velocity BC at A_WHEEL_REF, the rigid body's own reference node
  * shorten the step from 0.002 s to 6.0e-6 s, the intended 0.18 mm pass

Every change is a deletion or a substitution of a name/number. No new keyword is
introduced anywhere.

Consequence to be aware of: ALLKE will still read ~0 for the wheel, because a massless
rigid body has no kinetic energy however fast it turns. That is expected, not a fault.
Confirm the wheel is turning with UR3 at the reference node instead -- it should reach
0.00721 rad (0.413 deg) by the end of the step.
"""
import math
import re

SRC = 'grinding_actual_33.inp'
DST = 'grinding_actual_33_FIXED2.inp'
STEP_TIME = 6.0e-6

t = open(SRC, encoding='ascii', errors='replace').read()
orig = len(t)
log = []


def cut(pattern, label, flags=re.M):
    global t
    t, n = re.subn(pattern, '', t, flags=flags)
    log.append('%-52s %d' % (label, n))
    return n


# 1. the coupling that fought the rigid body, and the node surface it used
cut(r'\*\* Constraint: Constraint-2\s*\n\*Coupling[^\n]*\n\*Kinematic\s*\n',
    'deleted *Coupling + *Kinematic (Constraint-2)')
cut(r'\*Surface, type=NODE, name=s_Set-9_CNS_[^\n]*\n[^\n]*\n',
    'deleted the orphaned node surface s_Set-9_CNS_')

# 2. the reference point itself: the assembly node and the two sets naming it
cut(r'^\*Node\s*\n\s*1,\s*0\.,\s*0\.,\s*0\.\s*\n', 'deleted assembly node 1 (the old RP)')
cut(r'^\*Nset, nset=Set-13\s*\n\s*1,\s*\n', 'deleted *Nset Set-13')
cut(r'^\*Nset, nset=m_Set-12\s*\n\s*1,\s*\n', 'deleted *Nset m_Set-12')
cut(r'^\*Nset, nset=s_Set-9,[^\n]*\n[^\n]*\n', 'deleted *Nset s_Set-9 (coupling slaves)')

# 3. drive the rigid body's own reference node
m = re.search(r'(\*Boundary, type=VELOCITY\s*\n)((?:Set-13[^\n]*\n)+)', t)
if not m:
    raise SystemExit('velocity BC not found')
t = t[:m.start(2)] + m.group(2).replace('Set-13', 'A_WHEEL_REF') + t[m.end(2):]
log.append('%-52s %d' % ('BC retargeted Set-13 -> A_WHEEL_REF',
                         m.group(2).count('\n')))

# 4. a step length that matches the intended pass
t, n = re.subn(r'(\*Dynamic, Explicit\s*\n)\s*,\s*0\.002\s*\n',
               r'\g<1>, %g\n' % STEP_TIME, t)
log.append('%-52s %d' % ('step time 0.002 -> %g s' % STEP_TIME, n))

open(DST, 'w', encoding='ascii', newline='\n').write(t)
print('wrote %s  (%.1f MB, was %.1f MB)' % (DST, len(t) / 1e6, orig / 1e6))
for l in log:
    print('  ' + l)

# ---- checks -------------------------------------------------------------
print()
fails = []


def chk(name, ok, d=''):
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', name, (': ' + d) if d else ''))
    if not ok:
        fails.append(name)


chk('no *Coupling / *Kinematic left', '*Coupling' not in t and '*Kinematic' not in t)
chk('no reference to Set-13 anywhere', 'Set-13' not in t)
chk('no reference to m_Set-12 anywhere', 'm_Set-12' not in t)
chk('no reference to s_Set-9 anywhere', 's_Set-9' not in t)
rows = re.search(r'\*Boundary, type=VELOCITY\s*\n((?:[^\n*][^\n]*\n)+)', t).group(1)
rows = [r.strip() for r in rows.strip().split('\n')]
chk('velocity BC on A_WHEEL_REF, 6 DOF', len(rows) == 6
    and all(r.startswith('A_WHEEL_REF') for r in rows), rows[-1])
chk('A_WHEEL_REF is the rigid body ref node',
    re.search(r'\*Nset, nset=A_WHEEL_REF, instance=WHEEL-1\s*\n\s*286923', t) is not None
    and 'ref node=WHEEL-1.WHEEL_REF' in t)
chk('workpiece still encastred', 'A_WP_BACK_FACE, ENCASTRE' in t)
chk('JH-2 deletion still armed',
    '*Depvar, delete=12' in t and 'ELEMENT DELETION=YES' in t)
chk('step time is the 0.18 mm pass', ', 6e-06' in t)
chk('no new keyword introduced (no MASS/ROTARYI added)',
    'type=MASS' not in t and 'type=ROTARYI' not in t)
chk('assembly has no leftover free node',
    re.search(r'\*Assembly[^\n]*\n(?:\*\*[^\n]*\n)*\*Node', t) is None)
print()
print('  %d failure(s)' % len(fails))
print()
print('  wheel rotation over the step : %.4f rad = %.3f deg'
      % (1200 * STEP_TIME, math.degrees(1200 * STEP_TIME)))
print('  check UR3 at node 286923, NOT ALLKE - the wheel is massless by design')
