# THEORY — the physics and mathematics of this project

This document is self-contained. It assumes mechanical-engineering background
but nothing about this codebase, and it derives or states **every** equation the
model uses, with the provenance of each. Read it and you know what the model
computes and why.

`README.md` is how to run it. `context.md` is how the code is laid out.
This is what it *means*.

---

## Table of contents

1. [The physical question](#1-the-physical-question)
2. [Why a size effect exists at all](#2-why-a-size-effect-exists-at-all)
3. [The critical depth of cut, three ways](#3-the-critical-depth-of-cut-three-ways)
4. [Chip thickness: the kinematics](#4-chip-thickness-the-kinematics)
5. [The brittle branch: Johnson–Holmquist II](#5-the-brittle-branch-johnsonholmquist-ii)
6. [The ductile branch: Johnson–Cook + strain gradients](#6-the-ductile-branch-johnsoncook--strain-gradients)
7. [The switch](#7-the-switch)
8. [From SEM micrograph to grain solid](#8-from-sem-micrograph-to-grain-solid)
9. [Wheel assembly and geometry](#9-wheel-assembly-and-geometry)
10. [Against classical grinding theory](#10-against-classical-grinding-theory)
11. [Discretisation: mesh, time step, cost](#11-discretisation-mesh-time-step-cost)
12. [Material cards](#12-material-cards)
13. [What is measured, what is modelled, what is assumed](#13-what-is-measured-what-is-modelled-what-is-assumed)
14. [Units](#14-units)
15. [Symbols](#15-symbols)
16. [References](#16-references)

---

## 1. The physical question

Grind a brittle solid — silicon carbide, silicon, glass, a quartz-bonded rock —
and the material does not always break. Below a critical depth of cut it flows
plastically and leaves a smooth groove; above it, it fractures and leaves a
chipped one. That is the **ductile-to-brittle transition**, and it is the whole
subject of this project.

It matters because ductile-regime grinding produces optical-quality surfaces in
materials that are nominally unmachinable. The engineering question is:
**where along a single grit's path does the transition happen, and what sets
that location?**

The model answers it by simulating one abrasive grain — a grain whose shape was
**measured from an SEM micrograph**, not idealised as a sphere or a cone —
scratching a workpiece, with a constitutive law that follows plasticity where
the cut is shallow and fracture where it is deep.

The output that answers the question is a single field: **SDV13**, the branch
each material point took. Plotted over the groove, it is a map of where the
transition sits.

---

## 2. Why a size effect exists at all

The transition is not a threshold someone inserted into the physics. It falls
out of a competition between two ways of getting rid of material.

**Plastic removal** costs work proportional to **volume**. Push a grit a depth
*d* over an area *A*: the plastic work is roughly

    W_plastic  ~  H · A · d

where *H* is the hardness (a flow stress in indentation).

**Fracture removal** costs energy proportional to **new surface area**. Creating
crack area *A* in a material of toughness *K_c* and modulus *E* costs

    W_fracture  ~  G_c · A  =  (K_c² / E) · A

using *G_c = K_c²/E*, Griffith's relation between toughness and critical energy
release rate.

Material takes the cheaper route. Setting the two equal:

    H · A · d_c  =  (K_c² / E) · A

           ┌────────────────────────┐
           │   d_c  ~  K_c² / (E H) │
           └────────────────────────┘

The area cancels. What remains is a **length**, built only from material
properties. Below it, plastic flow is cheaper because the volume is small.
Above it, fracture is cheaper because volume grows as *d* while surface does
not.

This is the essential insight: **the transition is a length scale, and it exists
because two energy costs scale differently with depth.** For SiC that length is
tens of nanometres — which is why ductile-regime grinding needs such fine
control, and why the mesh in this project is sized in nanometres.

The size effect has a second, independent origin in the plastic branch itself
(strain-gradient hardening, §6), which makes shallow cuts *harder* rather than
merely cheaper. Both are in the model.

---

## 3. The critical depth of cut, three ways

The scaling argument gives `d_c ~ K_c²/(EH)` but not the constant, and the
literature contains **several inequivalent expressions**. This project
implements three and does not pretend they agree.

**Form 1** — the (H/E)^½ variant:

    d_c  =  λ_c · (H/E)^0.5 · (K_c/H)²

**Form 2** — Bifano, Dow & Scattergood (1991), with a calibrated λ_c = 0.15
across a wide range of brittle solids:

    d_c  =  λ_c · (E/H) · (K_c/H)²

**Form 3** — the local energy criterion, written pointwise on plastic work
rather than geometry (§7.2):

    d_c  =  Ψ · K_c²/(E H)  =  Ψ · (H/E)^1 · (K_c/H)²

### These are not the same formula

All three share the `(K_c/H)²` factor — the part carrying the dimension of
length, since `K_c/H` has units of √length. They differ in the **exponent on
H/E**: +½, −1, and +1 respectively.

The consequence is large. Forms 1 and 2 differ by `(E/H)^1.5`:

| material | H (MPa) | E (MPa) | E/H | d_c form 1 | d_c form 2 | ratio |
|---|---|---|---|---|---|---|
| sandstone | 1 000 | 6 500 | 6.5 | 5.295 nm | 87.750 nm | 16.6× |
| silicon carbide | 25 000 | 450 000 | 18 | 0.693 nm | 52.920 nm | 76.4× |

**λ_c is therefore not transferable between forms.** A λ_c calibrated for form 2
inserted into form 1 gives an answer wrong by more than an order of magnitude.
The card carries a selector (`PROPS(51)`, `IDCF`) and every deck records which
form produced its d_c.

Bifano's form 2 with λ_c = 0.15 is used throughout this project's shipped decks,
because that is the combination with published calibration behind it.
**SiC: d_c = 52.920 nm. Sandstone: d_c = 87.750 nm.**

> The three forms are *not* reconciled into one authority, deliberately. They
> come from the same indentation-fracture family but reference crack initiation
> differently, and choosing one silently would hide a real modelling decision.

---

## 4. Chip thickness: the kinematics

The constitutive law needs to know, at every material point, **how deep the cut
is where that point sits**. Call it *h*, the undeformed chip thickness. A VUMAT
sees no kinematics — it gets strain increments, not the wheel — so *h* must be
handed to it in closed form.

### 4.1 Derivation for a single grit

With one grit the trajectory is exact. The tip rides a circle of radius
`r_tip` about an axis translating outward at infeed speed `v_r`, and sweeps the
workpiece at surface speed `v_s`. Let *u* be the tangential station of a
material point, measured along the workpiece surface.

The grit arrives at station *u* at time

    t(u)  =  (u₀ − u) / v_s

Its reach along the radial direction at that moment is the circle plus the
infeed:

    r(u)  =  r_tip · cos(u / r_tip)  +  v_r · t(u)

Chip thickness is that reach minus the already-ground radius. Expanding the
cosine for `u ≪ r_tip` — true here, microns against 25 mm —

    r_tip · cos(u/r_tip)  ≈  r_tip − u²/(2 r_tip)

and substituting *t(u)* gives the form the card carries:

           ┌──────────────────────────────────────┐
           │  h(u) = H₀ + H_G · u − u²/(2·r_tip)  │
           └──────────────────────────────────────┘

with

    H_G  =  dh/du  =  −v_r / v_s
    H₀   =  (r_tip − r_ground)  +  v_r · u₀ / v_s

**Three terms, three meanings:**

- `H₀` — the depth at the block centre.
- `H_G · u` — the **wedge**. This is the rubbing → ploughing → shearing ramp in
  every textbook figure of a grit trajectory. Here it comes from radial infeed
  rather than table feed.
- `−u²/(2 r_tip)` — the **sagitta** of the grit's circular path. Tens of
  nanometres over a 50 µm block on a 50 mm wheel. Negligible in most grinding
  analyses — **but d_c is of that same order**, so here it is not.

### 4.2 Consistency with traverse grinding

The classical kinematic chip thickness for traverse grinding,

    h(θ)  =  L_g · (v_w/v_s) · sin θ

linearises to `H₀ + H_G·u` over a block much shorter than the contact arc. So a
traverse case enters through the same three constants; the form is general.

### 4.3 The sagitta relation, used for geometry

The same quadratic gives the chord needed to reach a given depth on a wheel of
radius *R*. For sagitta *D*:

    D = L²/(8R)        ⟹        L = √(8 R D)

This sets the length of the trajectory jobs (§9.3) and is the reason a 26 µm arc
on a 25 mm wheel is 2.28 mm long, not 17 mm.

### 4.4 Where the transition sits

The transition stations are the roots of `h(u) = d_c`:

    H₀ + H_G·u − u²/(2 r_tip)  =  d_c

a quadratic, so there can be **zero, one, or two** crossings inside the block.
A monotone ramp gives one. **A prescribed arc that rises above d_c and comes
back down gives two** — ductile, then brittle, then ductile again.

> This mattered. The original implementation bisected between the block
> endpoints assuming *h* was monotone. On an arc, where both endpoints are at
> zero depth, it reported *"entirely ductile"* for a deck that was 58 % brittle.
> The code now scans for every root and the deck header lists them all.

### 4.5 Why H₀ is pinned, not computed

`H₀` is **not** evaluated from `r_tip` and `r_ground` independently. It is
pinned to the deck's own tangency: at *t* = 0 the governing grit vertex sits
exactly the standoff clear of the ground face, so *h* at that vertex's station
is known to the picometre and `H₀` follows.

Deriving it any other way lets a sub-micron disagreement between "tallest tip"
and "tallest tip inside the footprint" leak into a quantity compared against a
53 nm threshold.

---

## 5. The brittle branch: Johnson–Holmquist II

JH-2 is a pressure-dependent, damage-coupled strength model for brittle solids
under high pressure and high rate. Implementation follows JH94, with Gazonas'
algorithmic ordering and Cronin's clarifications.

### 5.1 Normalisation

All strengths are normalised by the equivalent stress at the Hugoniot elastic
limit; pressures by the pressure at the HEL:

    σ* = σ / σ_HEL        P* = P / P_HEL        T* = T / P_HEL

with *P* positive in compression, so `P = −⅓·tr(σ)` under Abaqus'
tension-positive convention.

### 5.2 The two strength surfaces

**Intact:**

    σ_i*  =  A · (P* + T*)^N · (1 + C·ln ε̇*)

**Fractured:**

    σ_f*  =  B · (P*)^M · (1 + C·ln ε̇*)        capped at σ_f*max

**Current**, interpolated by damage `D ∈ [0,1]`:

    σ*  =  σ_i*  −  D · (σ_i* − σ_f*)

The yield surface in Mises stress is `q_lim = σ* · σ_HEL`, enforced by radial
return in J2 deviatoric space.

### 5.3 Damage

Damage accumulates as plastic strain against a pressure-dependent failure
strain:

    D  =  Σ  Δε_p / ε_f(P)          ε_f(P)  =  D₁ · (P* + T*)^D₂

More pressure, more ductility before failure — the observed behaviour of
confined brittle solids.

### 5.4 Bulking (the dilatancy that makes JH-2 distinctive)

When damage lowers the strength surface, the deviatoric elastic energy released
is converted to pressure. Fractured material occupies more volume than intact;
under confinement that generates pressure. The energy released is

    ΔU  =  (q_pre² − q_post²) / (6G)

and the bulking pressure follows from the energy balance

    ΔP_new  =  −K₁μ  +  √[ (K₁μ + ΔP_old)²  +  2·β·K₁·ΔU ]

with β ∈ [0,1] the bulking fraction and `μ = ρ/ρ₀ − 1`. Generated in
compression only.

*Verified: this reproduces the JH94 bulking benchmark — 0.559 / 0.711 /
0.646 GPa — to 1.3 %.*

### 5.5 Equation of state and tensile cutoff

    P  =  K₁μ + K₂μ² + K₃μ³  +  ΔP        (μ > 0)
    P  =  K₁μ                             (μ ≤ 0)

Tension is capped at `P ≥ −T(1−D)`: a fully damaged point carries no hydrostatic
tension.

---

## 6. The ductile branch: Johnson–Cook + strain gradients

### 6.1 Johnson–Cook flow stress

    σ_JC  =  [ A + B·(ε̄_p)^n ] · [ 1 + C·ln ε̇* ] · [ 1 − (T*)^m ]

with `ε̇* = ε̄̇_p / ε̇₀` and the homologous temperature
`T* = (T − T₀)/(T_melt − T₀)`. Three multiplicative factors: strain hardening,
rate hardening, thermal softening.

Adiabatic heating from plastic work:

    ΔT  =  β_TQ · σ · Δε̄_p / (ρ · c_p)

with β_TQ the Taylor–Quinney fraction (0.9 here).

### 6.2 Strain-gradient enhancement — the second size effect

Plain Johnson–Cook has **no length scale**, so it cannot represent the "smaller
is stronger" behaviour that dominates at cut depths of tens of nanometres. The
strain-gradient term supplies one.

The mechanism is **geometrically necessary dislocations**. A plastic strain
gradient across a small volume requires extra dislocations purely to keep the
lattice compatible. Their density is

    η  =  4 · ε̄_p / ℓ

where ℓ is the characteristic length — **here the local uncut chip thickness
h**, which is what ties this term to the grinding kinematics. η is a *density
per unit length*; multiplied by the Burgers vector *b* it becomes a dislocation
density.

By Taylor hardening, the GND contribution to flow stress adds **in quadrature**
with the statistical (Johnson–Cook) part:

    σ_eff²  =  σ_JC²  +  r' · η · b · (M·α·G)²

Written as an amplification factor, with Λ the SGE exponent:

    ┌────────────────────────────────────────────────────┐
    │  σ_eff = σ_JC · √[ 1 + (Σ·η / σ_JC²)^Λ ]           │
    │  where  Σ = r' · b · (M·α·G)²                      │
    └────────────────────────────────────────────────────┘

At Λ = 1 this is exactly the quadrature sum above. *M* is the Taylor factor, α
the Taylor hardening constant, *G* the shear modulus, *r'* the GND coefficient.

**The amplification is reported as SDV19.** A value of 1 means no size effect;
larger means the shallow cut is being hardened. Plot it to see how hard the size
effect is working.

ℓ is floored at one Burgers vector — as `h → 0` in the rubbing zone, `η → ∞`
and the stress would diverge. This is the regularisation that keeps the rubbing
zone finite.

### 6.3 Johnson–Cook damage

    ε_f  =  [ D₁ + D₂·exp(D₃·σ*) ] · [ 1 + D₄·ln ε̇* ] · [ 1 + D₅·T* ]

with `σ* = −P/q` the stress triaxiality. Damage accumulates as
`D = Σ Δε̄_p / ε_f`, and the carried stress is `(1−D)·σ_eff`. The element is
deleted at `D ≥ D_crit`.

*Provenance: the SGE formulation follows Yadav et al. across three papers
(micro-milling 2022, peening 2024, blanking 2026), which write the same
Taylor/GND hardening with different characteristic lengths. Here ℓ = h, the
micro-milling form.*

---

## 7. The switch

### 7.1 The geometric criterion

    h(u)  <  d_c    →    ductile   (Johnson–Cook + SGE)
    h(u)  ≥  d_c    →    brittle   (Johnson–Holmquist II)

*h* is evaluated **once**, at the first call, from *undeformed* coordinates, and
latched in SDV14. Recomputing it later would track the deformed position and
drift as the material moves.

### 7.2 The energy criterion (multi-abrasive)

The geometric criterion needs a known trajectory. With hundreds of grits, a
traverse, or a second pass over the same groove, there isn't one.
`vumat_grind2.for` adds a purely **local** criterion needing no geometry at all:

           ┌───────────────────────────────────────┐
           │  W_p · L_c  ≥  Ψ · K_c²/E  →  brittle │
           └───────────────────────────────────────┘

`W_p` is accumulated plastic work per unit volume; `L_c` the element's own
characteristic length, which Abaqus supplies.

Ψ is **not** defaulted to λ_c — the exponents on H/E do not line up (§3).
Instead it is defaulted to the value making the local criterion trip at exactly
the d_c the deck already chose:

    Ψ  =  d_c · E · H / K_c²        ⟹        W_p · L_c  ≥  H · d_c

which reads plainly: **brittle once the plastic work per unit area exceeds the
cost of plastically removing a layer of thickness d_c at flow stress H.**
Whichever d_c form the card carries, the two criteria then agree by
construction.

**This criterion triggers on history, not position.** A point starts ductile and
turns brittle as the cut deepens under it — which is the physical transition,
rather than being told in advance where it happens. Once triggered it latches.

On a flip, the JH-2 branch **inherits** the Johnson–Cook damage. The two damages
are different mechanisms with the same meaning — the fraction of the way to
failure — and carrying it over is the only continuation that neither forgives
nor double-counts the damage already accumulated.

`SWMODE` (`PROPS(57)`): 0 geometric only · 1 energy only · 2 either.

### 7.3 Making the two branches meet

The Johnson–Cook yield *A* is **not** invented. It is set to the quasi-static
uniaxial compressive strength of the *same JH-2 card* — the intersection of the
elastic path with the intact strength surface. Starting the ductile branch's
yield there makes the two laws **agree at the transition** in uniaxial
compression, instead of stepping discontinuously across it.

### 7.4 Verification modes

`IHMODE` (`PROPS(56)`) forces a branch for testing: 0 from coordinates · 1 from
field variable · 2 all ductile · 3 all brittle. **At `IHMODE = 3` every stress
component matches standalone `vumat_jh2.for` to 0 ulp** — the hybrid is provably
a strict superset of the brittle model.

---

## 8. From SEM micrograph to grain solid

The grain shapes are **measured**, not idealised. This is what distinguishes the
model from a sphere-on-a-wheel analysis.

### 8.1 Calibration

Pixel size comes from the SEM databar, parsed and cross-checked. Every length
downstream inherits it, so this is checked first and reported per image.

### 8.2 Segmentation — fourteen recorded stages

A marker-controlled watershed, with every intermediate retained for inspection
(`STAGE_KEYS`):

1. raw → 2. denoised → 3. threshold (Otsu) → 4. morphology (open/close) →
5. foreground → 6. distance transform → 7. seeds (local maxima) → 8. gradient →
9. elevation → 10. watershed → 11. boundary evidence → 12. merge → 13. labels →
14. border classification

The distance transform,

    dist(x)  =  min ‖x − y‖   over all y outside the grain

is what separates *touching* grains: its local maxima are grain centres, and the
watershed floods from them, so a cluster splits at the necks rather than being
counted as one grain.

Boundaries carry an **evidence test** — a boundary with no gradient support is
merged away, and each verdict is recorded, so over-segmentation is visible
rather than silent.

Grains touching the image border are **classified and excluded** from
statistics: a truncated grain has a real area but not its own, and including it
biases the size distribution downward.

### 8.3 Descriptors

25 per grain, in microns: area, perimeter, equivalent diameter, major/minor
axis, aspect ratio, orientation, solidity, convexity, circularity `4πA/P²`,
eccentricity, extent, Feret diameters, and moment invariants.

### 8.4 Reconstruction to a 3D solid

The measured 2D outline is lofted to a closed polyhedron, preserving the
**concave** features — the actual cutting edges. Height follows from the
in-plane size via the measured aspect statistics.

Each solid is validated: watertight (every edge in exactly two faces),
consistently wound, positive volume, and its projected outline recovers the
measured outline to tolerance. **Negative controls exist** — an inward-wound
solid must be *rejected* — so the checks cannot pass vacuously.

### 8.5 Protrusion

A grit cuts only as deep as it stands proud of the bond. With mean protrusion
fraction *f* and grain height *H_g*, usable depth is `f · H_g`.

**This is a hard physical limit, and it bites.** The measured B4C grains are
~7 µm tall, so a 26 µm cut from one is not a meshing problem — it is impossible.
`verify_rigid_deck.py` refuses such a deck outright (the bond rim would be
driven 22 µm into the workpiece). The remedy is geometric scaling (§9.3), which
preserves shape exactly.

---

## 9. Wheel assembly and geometry

### 9.1 Placement

Grains are placed on the rim of a wheel of diameter *D* and width *W*, either at
a specified **areal density** ρ_A (grains/mm²) or as a **single grit**. For a
sector of angle θ and rim depth *t*, the rim area is

    A_rim  =  (θ/360°) · π · D · W        N_grits  =  ρ_A · A_rim

Only the outer annulus of depth *t* is modelled. Grinding is confined to a
shallow surface layer, and this is what makes high grain counts tractable.

### 9.2 Rigid wheel

The wheel is **discrete rigid** (R3D3/R3D4 facets). It is not the object of
study, it is three orders of magnitude stiffer than the cut, and making it
deformable would cost the whole element budget for no physics. All wheel motion
is applied at one reference node, `A_WHEEL_REF`, where `RF`/`RM` give the
grinding force and moment directly.

### 9.3 The trajectory jobs — a worked example

Three prescribed-trajectory jobs demonstrate the transition. They also show how
the theory constrains geometry.

A supplied figure showed a scallop groove: 0 → 26 µm depth over 0 → 17 mm.
Taken literally, §4.3 gives the required radius:

    R  =  L²/(8D)  =  (17 mm)² / (8 × 0.026 mm)  =  1389 mm

— a **2.8 m wheel**. The figure is vertically exaggerated about **650:1**. On
the actual 25 mm radius, 17 mm of arc would cut 1445 µm deep, 56× the figure.

So the **shape** is reproduced exactly and the **length** follows from the real
wheel via `L = √(8RD)`: 2.28 mm at 26 µm, 126.5 µm at 80 nm. Both are honest
arcs of the actual wheel; neither invents a radius.

Imposing the trajectory uses the card's existing three terms:

- an **arc** is `H_G = 0`: no infeed, the parabola alone carries *h* up and back
  down. The depth is the wheel's own curvature.
- a **ramp** is the linear term, with `H₀` and `H_G` solved so *h* = 0 at entry
  and peak at exit **with the parabola present** rather than pretending it is
  not.

These constants are **prescribed**, not derived from seating and infeed, because
here the trajectory is the specification rather than the consequence. Every deck
header says which, so an imposed number is never mistaken for a derived one.

| job | trajectory | peak h | h/d_c | ductile | crossings |
|---|---|---|---|---|---|
| `ARC_26UM` | arc | 26 µm | 491 | 0.2 % | 2 (at the ends) |
| `ARC_80NM` | arc | 80 nm | 1.51 | 41.9 % | 2, at u = ±36.8 µm |
| `RAMP_80NM` | linear | 80 nm | 1.51 | 63.8 % | 1, at 63.8 % along |

Job 1 is brittle essentially everywhere by design — it verifies path-following.
Job 2 is the picture the model exists to produce: **ductile → brittle →
ductile**. Job 3 is the control, where depth is linear in position so the
transition station can be read off the plot and checked against d_c by hand.

Job 1 scales the grain library by **×7.716** for the reason in §8.5. Every
length scales together, so the grain keeps its measured shape exactly — same
aspect ratio, corner angles, concave features — and only its size changes. Jobs
2 and 3 cut 80 nm, which any measured grain clears by two orders of magnitude,
and use the grains **exactly as measured**.

---

## 10. Against classical grinding theory

A deck can be geometrically perfect and still describe a process nobody would
call grinding — a depth of cut no grain reaches, a mesh too coarse to carry a
chip, a contact zone shorter than the block. Those are the first questions a
reviewer asks, and verifying the `.inp` answers none of them.
`semgrit/grinding_theory.py` reports them in two columns, deliberately
separated: **measured** rows counted off the geometry the deck actually
contains, and **theory** rows from the textbook expressions.

### 10.1 The classical expressions

Geometric contact length, for depth of cut a_e on equivalent diameter d_e:

    l_c  =  √(a_e · d_e)

Equivalent chip thickness — the thickness of a continuous layer removed at the
same rate:

    h_eq  =  a_e · v_w / v_s

Maximum chip thickness (Malkin), for active grain density C and grain shape
factor r:

    h_max  =  √[ (4 v_w)/(v_s · C · r) · √(a_e/d_e) ]

Specific material removal rate:

    Q'_w  =  a_e · v_w

### 10.2 Why the theory column needs a work speed

Every chip-thickness formula above assumes a **traverse** grind at work speed
v_w. This model is a rotating wheel with radial infeed against a fixed block —
a **plunge** configuration, where v_w does not exist.

So v_w must be supplied to say *which traverse case* the numbers correspond to.
With `work_speed_mm_s = 0` those rows are reported as **not applicable** rather
than quietly computed from a speed of zero. This is the difference between a
comparison and a fabrication.

### 10.3 A second transition at almost the same depth

This one is a genuine confound and worth stating plainly.

A blunt edge stops cutting and starts **ploughing** below a minimum chip
thickness set by the edge radius and the friction — not by the material's
fracture behaviour at all. From the stagnation point on a round edge
(Son/Lim's form):

    h_min / r_e  =  1 − cos(π/4 − β/2),        β = arctan μ

On the shipped decks, r_e = 0.35 µm and μ = 0.2 give

    h_min  =  0.07933 µm      against      d_c = 0.08775 µm

**9.6 % apart.** Two mechanisms, two thresholds, one force trace — and they
cannot be told apart in it. Only one of them (d_c) is this project's subject.

The implications are worth being explicit about:

- A force inflection near 80 nm on sandstone is **not** evidence of the
  ductile-brittle transition on its own. It could be the ploughing threshold.
- What *does* distinguish them is **SDV13**, which reports the branch the
  constitutive law actually took, independent of the force. This is a second
  reason the branch map is the quotable result and the force is not.
- The two can be separated experimentally by changing r_e (dressing the wheel)
  without changing the material: d_c stays put, h_min moves.

Nothing in the repo computed h_min before this was added, so the confound was
invisible. It is now reported alongside d_c on every deck.

## 11. Discretisation: mesh, time step, cost

### 10.1 Resolving d_c

The transition happens across a length of d_c. A mesh that cannot resolve d_c
cannot resolve the transition, no matter how good the constitutive law:

    Δx  =  d_c / N          N = 5   (ELEMENTS_PER_DC)

For SiC that is 52.92/5 = 10.58 nm through the surface layer. **This single
requirement sets the entire computational cost of the project.**

Away from the surface, elements grade coarser geometrically
(`wp_surface_layer_mm`, `wp_depth_growth`, `wp_width_band_mm`,
`wp_width_growth`), because resolution is only needed where the transition is.

### 10.2 Explicit stability

Abaqus/Explicit is conditionally stable. The step is bounded by the time for a
dilatational wave to cross the smallest element:

    Δt  ≤  L_min / c_d

    c_d  =  √[  E(1−ν)  /  ( ρ (1+ν)(1−2ν) )  ]

**This is why SiC is expensive.** Its wave speed is ~6.7× sandstone's, so on the
same mesh the stable increment is 6.7× smaller — 6.7× the CPU time for identical
geometry. That is material physics, not a mesh defect.

Total cost scales as

    cost  ~  N_elements · t_sim · c_d / L_min

Halving Δx costs 8× in elements (3D) **and** 2× in steps: **16× overall.** Hence
the graded mesh, and hence quoting wall-clock in the deck reports.

### 10.3 Precision

**`double=both` is mandatory on every run.** *h* is compared against a 53 nm
threshold on a 25 mm radius — a ratio of 2×10⁻⁶. Single precision has ~7 decimal
digits and does not have the digits to resolve that. The failure mode is silent:
the branch flag comes out wrong, the job does not crash.

### 10.4 Energy diagnostics

Explicit results are not trustworthy without checking:

- **Artificial (hourglass) energy** should be < 5 % of internal energy. Reduced
  integration (C3D8R) has zero-energy modes; if hourglass control is doing much
  work, the deformation is not physical.
- **Kinetic energy** should be small relative to internal for a
  quasi-static-in-spirit process, or mass scaling is dominating.

> Both are **currently out of bound** on the archived runs — hourglass 31–39 %,
> KE 320–56 000× internal. These must be resolved before any force or
> specific-energy figure from those runs is quoted. See `context.md`.

---

## 12. Material cards

56 constants, written **8 per line** — 4 per line is silently rejected by
Abaqus.

| slots | contents |
|---|---|
| 1–21 | JH-2. Identical to standalone `vumat_jh2.for`, so an existing card is a prefix of this one. |
| 22–40 | Johnson–Cook flow + strain-gradient enhancement |
| 41–46 | Johnson–Cook damage D₁..D₅, D_crit |
| 47–51 | the switch: d_c, λ_c, H, K_c, form selector |
| 52–55 | the chip field: θ_c, H₀, H_G, r_tip |
| 56 | `IHMODE` — how *h* is obtained |
| 57 | `SWMODE` (grind2 only) — which criterion applies |

**State variables** (20, `*Depvar delete=12`). The ones to plot:

| SDV | meaning |
|---|---|
| **13** | **branch: 1 ductile, 2 brittle — the result** |
| 14 | *h* this point was given |
| 15 | d_c actually used |
| **19** | SGE amplification σ_eff/σ_JC |
| 1–5 | D, ε̄_p, P, q, ε̇ |
| 12 | deletion flag |

### A unit trap worth naming

K_c is a stress intensity, so in a mm–MPa deck it is MPa·√mm — **31.623× the
MPa·√m** toughness is normally quoted in. Getting it wrong scales d_c by
**1000**: the difference between a transition at 5 nm and one at 5 µm, i.e.
between a model that is mostly brittle and one that is mostly ductile. Use
`kic_from_mpa_sqrt_m()`.

---

## 13. What is measured, what is modelled, what is assumed

Honesty about provenance is what makes the model usable in a paper.

**Measured** — from the SEM micrographs: grain outlines, all 25 descriptors,
size distributions, aspect statistics, pixel calibration.

**Modelled** — from measured inputs plus stated theory: 3D grain solids, wheel
placement, h(u), d_c, mesh grading, the constitutive response.

**Assumed / placeholder** — and this is the important list:

- **Johnson–Cook constants B, n, C, m and D₁..D₅ are placeholders** for both
  sandstone and SiC. Only *A* is derived (§7.3). They are defensible orders of
  magnitude and nothing more.
- λ_c = 0.15 is Bifano's calibration, not this material's.
- Grain height from in-plane size via aspect statistics, not measured directly.
- Bond is elastic; wheel is rigid.
- No coolant, no wheel wear, no thermal field beyond adiabatic JC heating.

### What may and may not be quoted

> **The branch map (SDV13) is the result.** It follows from *h* vs d_c, both of
> which are derived from measured geometry and published material properties.
>
> **Force magnitudes are not**, until B, n, C, m, D₁..D₅ are calibrated against
> real nanoindentation or scratch data. The *shape* of a force curve is
> informative; its absolute value is not.

An open question sits in the same category: an archived SiC run reached ~40 GPa
peak Mises, 2.8× SiC's HEL. Whether that is legitimate uncapped intact-surface
JH-2 behaviour or an SGE artefact from `h → 0` in the rubbing zone is unresolved
— SDV19 is the diagnostic.

---

## 14. Units

**mm – MPa – tonne – s**, consistently. This is Abaqus' standard mm system.

| quantity | unit | note |
|---|---|---|
| length | mm | d_c, h, H₀, b all in mm |
| stress | MPa | |
| mass | tonne | = 10³ kg |
| density | tonne/mm³ | ρ_SI × 10⁻¹² |
| force | N | |
| toughness | MPa·√mm | **×31.623** from MPa·√m |
| specific heat | given J/(kg·K) | converted on the way into the card |

Two conversions cause most unit bugs here: toughness (§12) and density.

---

## 15. Symbols

| symbol | meaning | unit |
|---|---|---|
| h(u) | undeformed chip thickness at station u | mm |
| d_c | critical depth of cut | mm |
| u | tangential station along the scratch | mm |
| H₀, H_G | chip-field offset and wedge slope | mm, – |
| r_tip | grit tip radius | mm |
| θ_c | workpiece centre angle | rad |
| H | hardness (indentation flow stress) | MPa |
| E, G, ν | Young's modulus, shear modulus, Poisson | MPa, MPa, – |
| K_c | fracture toughness | MPa·√mm |
| λ_c | d_c prefactor | – |
| Ψ | energy-criterion prefactor | – |
| D | damage, 0..1 | – |
| ε̄_p | equivalent plastic strain | – |
| P, q | pressure (+compression), Mises stress | MPa |
| μ | volumetric compression ρ/ρ₀ − 1 | – |
| β | JH-2 bulking factor | – |
| η | GND density 4ε̄_p/ℓ | 1/mm |
| b | Burgers vector | mm |
| M, α, r', Λ | Taylor factor, hardening const., GND coeff., SGE exponent | – |
| β_TQ | Taylor–Quinney fraction | – |
| c_d | dilatational wave speed | mm/s |
| v_s, v_r, v_w | surface, radial infeed, workpiece speed | mm/s |
| W_p | plastic work per unit volume | MPa |
| L_c | element characteristic length | mm |

---

## 16. References

**Brittle constitutive model**

- Johnson, G.R. & Holmquist, T.J. (1994). *An improved computational
  constitutive model for brittle materials.* AIP Conf. Proc. **309**, 981–984.
- Gazonas, G.A. (2002). *Implementation of the Johnson–Holmquist II (JH-2)
  constitutive model in LS-DYNA.* ARL-TR-2699.
- Cronin, D.S. et al. (2003). *Implementation and validation of the
  Johnson–Holmquist ceramic material model in LS-DYNA.* 4th European LS-DYNA
  Conference.

**Ductile constitutive model**

- Johnson, G.R. & Cook, W.H. (1983). *A constitutive model and data for metals
  subjected to large strains, high strain rates and high temperatures.* 7th Int.
  Symp. Ballistics.
- Johnson, G.R. & Cook, W.H. (1985). *Fracture characteristics of three metals
  subjected to various strains, strain rates, temperatures and pressures.*
  Eng. Fract. Mech. **21**(1), 31–48.

**Strain-gradient / size effect**

- Yadav, Chakladar & Paul (2022). Int. J. Mech. Sci. **231**, 107582, eqs. 8–10.
  *(micro-milling; ℓ = uncut chip thickness — the form used here)*
- Yadav, Das Chakladar & Paul (2024). Int. J. Mach. Tools Manuf. **194**,
  104100, eqs. 24–26. *(peening; r' = 2.0)*
- Yadav, Jewell, Jones & Ghadbeigi (2026). Int. J. Mech. Sci. **314**, 111375,
  eqs. 4–7. *(blanking; Λ = 1.0)*

**Ductile-to-brittle transition**

- Bifano, T.G., Dow, T.A. & Scattergood, R.O. (1991). *Ductile-regime grinding:
  a new technology for machining brittle materials.* J. Eng. Ind. **113**,
  285–308. *(form 2, λ_c = 0.15)*
- Griffith, A.A. (1921). *The phenomena of rupture and flow in solids.*
  Phil. Trans. R. Soc. A **221**, 163–198. *(G_c = K_c²/E)*

---

## Where the equations live in the code

| theory | implementation | independent check |
|---|---|---|
| d_c, all three forms | `semgrit/hybrid.py::critical_depth_mm` | `verify_vumat_grind.py` |
| h(u) and its roots | `semgrit/hybrid.py::chip_field`, `transition_stations` | `verify_hybrid_deck.py` |
| JH-2 | `vumat_grind.for` | `verify_vumat_grind.py` (JH94 benchmark, 0-ulp vs `vumat_jh2.for`) |
| JC + SGE | `vumat_grind.for::grsge` | `verify_vumat_grind.py` (closed form) |
| energy criterion | `vumat_grind2.for` | `verify_vumat_grind2.py` |
| segmentation | `semgrit/segment.py` | `verify_all.py` |
| grain solids | `semgrit/grain3d.py` | `verify_all.py` (+ negative controls) |
| wheel geometry | `semgrit/rigid_wheel.py` | `verify_rigid_deck.py`, `verify_rigid_deck2.py` |
| material cards | `semgrit/materials.py` | `_check_presets.py` |
| classical comparison, h_min | `semgrit/grinding_theory.py` | reported per deck |
| figures | `semgrit/figures.py` | `python -m semgrit.figures` (`demo()`) |

Every verifier is deliberately written **without sharing code** with what it
checks, so a bug in the pipeline cannot be baked into its own verifier.
