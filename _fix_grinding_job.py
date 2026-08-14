"""Patch grinding_actual_33.inp so the wheel actually turns.

Three faults, found by reading the submitted deck:

1. THE VELOCITY BC DRIVES NOTHING.  CAE created a reference point (assembly node 1),
   kinematically coupled *every* wheel node (s_Set-9 = WHEEL-1 nodes 1..286923) to it,
   and put the -1200 rad/s BC on that point. But all of those nodes already belong to
   the rigid body `WHEEL-1.ES_WHEEL_ALL`, governed by `WHEEL-1.WHEEL_REF` (node
   286923). A node cannot be both a rigid-body member and a kinematic-coupling slave,
   so Abaqus discards the coupling -- the reference point spins and the wheel does not.
   Meanwhile the rigid body's own reference node carries no boundary condition at all.

2. THE WHEEL HAS NO MASS.  The point mass and rotary inertia written on the reference
   node did not survive the round trip through CAE. A massless rigid body has zero
   kinetic energy however fast it turns, so ALLKE would read zero even after fault 1
   is fixed. This is the symptom that was actually reported.

3. THE STEP IS 300x TOO LONG.  0.002 s at 1200 rad/s is 137.5 deg of rotation on a
   30 deg sector, sweeping the grits 60 mm past a 0.4 mm workpiece, and with mass
   scaling 30 it is 1.74 million increments -- about 19 days on 8 cores. The intended
   0.18 mm pass is 6.0e-6 s.

Everything else in the deck is left exactly as submitted: the JH-2 material, DEPVAR
with delete=12, ELEMENT DELETION=YES, general contact, mass scaling, output requests
and the encastre on the workpiece back face are all correct and untouched.
"""
import math
import os
import re
import sys

sys.path.insert(0, '.')

from semgrit.rigid_wheel import rim_mass_properties
from semgrit.wheel import WheelSpec

SRC = 'grinding_actual_33.inp'
DST = 'grinding_actual_33_FIXED.inp'
REF_NODE = 286923
STEP_TIME = 6.0e-6          # 0.18 mm pass at 30 m/s

t = open(SRC, encoding='ascii', errors='replace').read()
orig = len(t)
changes = []

# ---- 1. delete the kinematic coupling that fights the rigid body -------------
pat = re.compile(r'\*\* Constraint: Constraint-2\s*\n\*Coupling[^\n]*\n\*Kinematic\s*\n')
t, n = pat.subn('** Constraint-2 (kinematic coupling) removed: it slaved every wheel\n'
                '** node to a reference point, but those nodes are already in the rigid\n'
                '** body, so Abaqus discarded it and the wheel was never driven.\n', t)
changes.append('removed the kinematic coupling            : %d' % n)

# ---- 2. drive the rigid body's own reference node ----------------------------
m = re.search(r'(\*Boundary, type=VELOCITY\s*\n)((?:Set-13[^\n]*\n)+)', t)
if not m:
    raise SystemExit('could not find the velocity BC')
block = m.group(2).replace('Set-13', 'A_WHEEL_REF')
t = t[:m.start(2)] + block + t[m.end(2):]
changes.append('BC moved Set-13 -> A_WHEEL_REF (node %d) : %d lines'
               % (REF_NODE, block.count('\n')))

# The old reference point is now unused. Fully fix it so it cannot be reported as a
# node with no mass.
t = t.replace('** \n** INTERACTIONS\n',
              '** \n** The old reference point is no longer driving anything; fix it so it\n'
              '** cannot be flagged as a massless free node.\n'
              '*Boundary\nSet-13, ENCASTRE\n** \n** INTERACTIONS\n', 1)
changes.append('old reference point Set-13 encastred      : 1')

# ---- 3. put the mass and rotary inertia back on the reference node -----------
spec = WheelSpec(diameter_mm=50.0, width_mm=3.0, sector_deg=30.0, rim_depth_mm=15.0)
mass, inertia = rim_mass_properties(spec, 2700.0)

wheel = re.search(r'(\*Part, name=WHEEL\s*\n)(.*?)(\*End Part)', t, re.S)
body = wheel.group(2)
eid = max(int(l.split(',')[0]) for blk in
          re.finditer(r'\*Element,[^\n]*\n((?:[ \t]*\d[^\n]*\n)+)', body)
          for l in blk.group(1).strip().split('\n'))
add = ('*Element, type=MASS, elset=ES_WHEEL_MASS\n%d, %d\n'
       '*Mass, elset=ES_WHEEL_MASS\n%.9e,\n'
       '*Element, type=ROTARYI, elset=ES_WHEEL_ROTI\n%d, %d\n'
       '*Rotary Inertia, elset=ES_WHEEL_ROTI\n'
       '%.9e, %.9e, %.9e, %.9e, %.9e, %.9e\n'
       % (eid + 1, REF_NODE, mass, eid + 2, REF_NODE,
          inertia[0, 0], inertia[1, 1], inertia[2, 2],
          inertia[0, 1], inertia[0, 2], inertia[1, 2]))
t = t[:wheel.start(3)] + add + t[wheel.start(3):]
changes.append('mass %.4e t + rotary inertia restored : elements %d,%d'
               % (mass, eid + 1, eid + 2))

# ---- 4. a step length that matches the intended pass -------------------------
t, n = re.subn(r'(\*Dynamic, Explicit\s*\n)\s*,\s*0\.002\s*\n',
               r'\g<1>, %g\n' % STEP_TIME, t)
changes.append('step time 0.002 -> %g s                : %d' % (STEP_TIME, n))

open(DST, 'w', encoding='ascii', newline='\n').write(t)

print('wrote %s  (%.1f MB, was %.1f MB)' % (DST, len(t) / 1e6, orig / 1e6))
for c in changes:
    print('  ' + c)

omega = 1200.0
print()
print('after the fix:')
print('  wheel rotation over the step : %.3f deg' % math.degrees(omega * STEP_TIME))
print('  grit sweep at r = 25 mm      : %.3f mm' % (25 * omega * STEP_TIME))
print('  wheel kinetic energy         : %.4e mJ  (0.5*I33*w^2, was 0 with no mass)'
      % (0.5 * inertia[2, 2] * omega ** 2))
c_d = math.sqrt(50000 * 0.75 / (1.25 * 0.5 * 2.65e-9))
dt = 0.001 / c_d * math.sqrt(30.0)
print('  increments (mass scaling 30) : %s' % format(int(STEP_TIME / dt), ','))
print('  estimate                     : %.1f h on 8 cores, %.1f h on 4'
      % (1e6 * STEP_TIME / dt / 1.05e6 / 3600, 1e6 * STEP_TIME / dt / 5.25e5 / 3600))
