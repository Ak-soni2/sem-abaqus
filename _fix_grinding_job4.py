"""FIXED4: make the wheel actually cut.

The previous run completed with no contact because the wheel had no radial infeed.
V1 = V2 = V3 = 0 means the wheel centre never approaches the work: it spins in place,
one grit grazes at t = 0 and rotates away, and every other grit sits 1-2 um below the
surface for ever. No amount of extra time or rpm changes that -- they only sweep more
grits past at the same radius.

Four changes:

  depth of cut   3.0 um of radial infeed, applied as a constant approach velocity along
                 -(cos 15, sin 15). This is the change that creates contact. 3.0 um
                 leaves 0.56 um before the bond rim itself would touch (clearance is
                 3.557 um), and the depth grows through the step so engagement builds
                 from one grit to dozens -- which also animates better than a step jump.

  wheel speed    1200 -> 2000 rad/s (30 -> 50 m/s), as asked. This is free: the wheel is
                 rigid, so the stable increment comes from the workpiece mesh alone and
                 does not care how fast the wheel turns.

  step time      6 -> 10 us, as asked. Sweep = omega * T * R = 0.500 mm, which is just
                 inside the 0.543 mm the dressed band allows before the workpiece would
                 run off it into bare bond.

  mass scaling   30 -> 10. The budget is 10-12 h and the sweep is capped by the dressed
                 band, so extra simulated time buys nothing. Spending the budget on less
                 mass scaling instead is a real gain: factor 30 inflates density 30x and
                 badly distorts inertia at 50 m/s.
"""
import math
import re

SRC = 'grinding_actual_33_FIXED3.inp'
DST = 'grinding_actual_33_FIXED4.inp'

OMEGA = 2000.0          # rad/s  = 50 m/s at r = 25 mm
STEP = 1.0e-5           # s
AE = 3.0e-3             # mm of radial infeed over the step
THETA = math.radians(15.0)
SCALE = 10.0            # mass scaling factor

R, CLEAR = 25.0, 3.557e-3
v_r = AE / STEP
V1, V2 = -math.cos(THETA) * v_r, -math.sin(THETA) * v_r

t = open(SRC, encoding='ascii', errors='replace').read()
log = []

# 1. step time
t, n = re.subn(r'(\*Dynamic, Explicit\s*\n)\s*,\s*6e-06\s*\n',
               r'\g<1>, %g\n' % STEP, t)
log.append('step time 6e-06 -> %g s' % STEP)
assert n == 1

# 2. rotation + radial infeed, in one BC block
old = re.search(r'\*Boundary, type=VELOCITY\s*\n(?:A_WHEEL_REF[^\n]*\n){6}', t)
new = ('*Boundary, type=VELOCITY\n'
       'A_WHEEL_REF, 1, 1, %.6f\n'
       'A_WHEEL_REF, 2, 2, %.6f\n'
       'A_WHEEL_REF, 3, 3\n'
       'A_WHEEL_REF, 4, 4\n'
       'A_WHEEL_REF, 5, 5\n'
       'A_WHEEL_REF, 6, 6, %.1f\n' % (V1, V2, -OMEGA))
t = t[:old.start()] + new + t[old.end():]
log.append('infeed V1=%.3f V2=%.3f mm/s  (%.1f um over the step)' % (V1, V2, AE * 1000))
log.append('rotation VR3 -1200 -> %.0f rad/s (%.0f m/s)' % (-OMEGA, OMEGA * R / 1000))

# 3. mass scaling
t, n = re.subn(r'\*Fixed Mass Scaling, factor=30\.', '*Fixed Mass Scaling, factor=%g.' % SCALE, t)
log.append('mass scaling 30 -> %g' % SCALE)
assert n == 1

open(DST, 'w', encoding='ascii', newline='\n').write(t)
print('wrote %s  (%.1f MB)' % (DST, len(t) / 1e6))
for l in log:
    print('  ' + l)

# ---- checks ----------------------------------------------------------------
F = []
def chk(n, ok, d=''):
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', n, (': ' + d) if d else ''))
    if not ok:
        F.append(n)

print()
rows = re.search(r'\*Boundary, type=VELOCITY\s*\n((?:A_WHEEL_REF[^\n]*\n)+)', t).group(1)
rows = [r.strip() for r in rows.strip().split('\n')]
chk('BC still on the rigid-body ref node, all 6 DOF', len(rows) == 6
    and all(r.startswith('A_WHEEL_REF') for r in rows))
chk('radial infeed is now non-zero', rows[0].count(',') == 3 and rows[1].count(',') == 3,
    '%s | %s' % (rows[0], rows[1]))
chk('infeed points radially INWARD', V1 < 0 and V2 < 0,
    'unit (%.4f, %.4f) vs block at theta=15 deg' % (V1 / v_r, V2 / v_r))
chk('infeed magnitude matches the radial speed',
    abs(math.hypot(V1, V2) - v_r) < 1e-9, '%.3f mm/s' % math.hypot(V1, V2))
chk('depth of cut stays clear of the bond rim', AE < CLEAR,
    '%.2f um cut vs %.2f um clearance, margin %.2f um'
    % (AE * 1000, CLEAR * 1000, (CLEAR - AE) * 1000))
sweep = OMEGA * STEP * R
chk('sweep stays on the dressed band', sweep <= 0.543,
    '%.3f mm swept, %.3f mm available' % (sweep, 0.543))
chk('rotation increased as asked', OMEGA > 1200, '%.0f rad/s = %.0f m/s' % (OMEGA, OMEGA*R/1000))
chk('step time increased as asked', STEP > 6e-6, '%g s' % STEP)
chk('mass scaling reduced (better fidelity)', SCALE < 30, 'factor %g' % SCALE)
chk('JH-2 deletion still armed', '*Depvar, delete=12' in t and 'ELEMENT DELETION=YES' in t)
chk('SDV still requested', 'SDV' in t)
chk('wheel mass still present', 'type=MASS' in t and 'type=ROTARYI' in t)
print()
print('  %d failure(s)' % len(F))

c = math.sqrt(50000 * 0.75 / (1.25 * 0.5 * 2.65e-9))
dt = 0.001 / c * math.sqrt(SCALE)
inc = STEP / dt
print()
print('  dt (mass scaling %g)   : %.3e s' % (SCALE, dt))
print('  increments            : %s' % format(int(inc), ','))
print('  estimate              : %.1f h on 8 cores, %.1f h on 4' %
      (1e6 * inc / 1.05e6 / 3600, 1e6 * inc / 5.25e5 / 3600))
print('  wheel sweep           : %.3f mm of arc (%.2f deg)' % (sweep, math.degrees(OMEGA*STEP)))
print('  grits crossing the 0.4 mm block : about %d' % int(sweep / 0.4 * 267))
print('  depth of cut at end   : %.2f um   (grows from 0)' % (AE * 1000))
