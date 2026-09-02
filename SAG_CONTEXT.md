# SAG — Shape Adaptive Grinding: context

Everything gathered before writing a line of SAG code, from the four sources
in the brief: the reference paper, the presentation deck, the reference `.inp`,
and the wider literature. Findings that change what we should build are marked
**FINDING**.

Companion to [`THEORY.md`](THEORY.md) (the physics of the rigid-wheel model).
This file is SAG-specific and states where SAG **departs** from that model.

---

## 1. What SAG is

**Shape Adaptive Grinding** is compliant abrasive finishing: the tool deforms
under load to match the workpiece, instead of the workpiece being forced to
match the tool. It sits between grinding (removes material) and polishing
(produces finish).

The tool is **three layers**:

| layer | material | role |
|---|---|---|
| ① base wheel | PMMA (paper) / aluminium (deck) | rigid, transmits spindle torque |
| ② **elastic layer** | **polyurethane** | **the whole point** — compresses, lets the tool conform |
| ③ abrasive layer | resin-bonded diamond pad | contacts and cuts |

### The mechanism, in one causal chain

    elastic layer compresses
      -> line contact becomes an ELLIPTICAL AREA contact
      -> the same total load is shared by many more grains
      -> force per grain collapses (1e-4 - 1e-5 N)
      -> penetration depth per grain collapses (nanometres)
      -> h < d_c at every grain
      -> DUCTILE removal, even on a brittle material

That last step is the same `h` vs `d_c` competition the rigid-wheel model
already implements. **SAG does not need a new transition criterion — it needs a
new way of computing `h`.** In the rigid model `h` comes from wheel kinematics;
in SAG it comes from *contact mechanics*. This is the single most important
structural insight for the port.

### Why it matters

Conventional rigid grinding on brittle materials: concentrated force → fracture,
chipping, subsurface cracks, high residual stress, and no conformity to freeform
surfaces. SAG: Ra ~20 nm, no detectable subsurface damage, residual stress cut
to 28 % of the ground value.

---

## 2. The reference paper

> Ghosh, G., Sidpara, A. & Bandyopadhyay, P.P. (2021). *Brittle-ductile
> transition in compliant finishing of HVOF sprayed hard WC-Co coating.*
> Int. J. Refractory Metals and Hard Materials **99**, 105610.
> https://doi.org/10.1016/j.ijrmhm.2021.105610

IIT Kharagpur. 14 pages, read in full.

### 2.1 Experimental setup

- Wheel **D = 125 mm, W = 10 mm**, PMMA base + PU foam + diamond pad.
- Attached to a CNC spindle; raster path, feed 15 mm/min; load cell 0.1 N.
- Workpiece: **HVOF-sprayed WC-12Co coating**, 300–400 µm thick, 1.5 % porosity.
- Mean WC carbide size **D_WC = 1.36 ± 0.29 µm**.

**Process parameters (Table 3)**

| parameter | levels |
|---|---|
| wheel speed N (rpm) | 300, 550, 800, 1050 |
| wheel compression T (mm) | 0.2, 0.4, 0.6 |
| abrasive size d_g (µm) | 6, 15, 30 |

**WC-Co coating properties (Table 2)**

| property | value |
|---|---|
| microhardness | 11.02 ± 1.2 GPa |
| elastic modulus | 200 ± 21 GPa |
| fracture toughness | 7.78 ± 0.9 MPa·√m |
| BHN | 581 kg/mm² |

### 2.2 The MRR model — every equation

Hertzian pressure over the elliptical patch, semi-axes a and b:

    (1)   p = p0 · [ 1 − (x/a)² − (y/b)² ]^(1/2)          p_av = (2/3) p0

Total normal load from wheel compression T on wheel radius R:

    (2)   F_N = 1.44 · E_eq · R^(1/2) · T^(3/2)

Equivalent modulus of the contacting pair:

    (3)   E_eq = [ (1−ν_w²)/E_w  +  (1−ν_t²)/E_t ]^(−1)

Backing-pad modulus from shore hardness S:

    (4)   E_t = 0.0981 · (56 + 7.62 S) / (0.1375 · (254 − 2.54 S))

Tangential load:

    (5)   F_T = μ · F_N

Empirical fits to the measured finishing spot (these are **fits, not theory**):

    (6)   A_s = 138.22 · ( T^0.151 · N^0.009 )      mm²,  R² = 0.991
    (7)   L_s = 17.69  · ( T^0.232 · N^0.012 )      mm,   R² = 0.981

Active abrasives over the spot (active density is 80 % of areal density):

    (8)   N_abr = C_a · A_s          C_a = 0.8 · C_o

Load per grain — the quantity the whole process turns on:

    (9)   F_n = F_N / N_abr
   (10)   F_t = F_T / N_abr

Penetration depth by Brinell indentation of a sphere d_g:

   (11)   BHN = 2 F_n / [ π d_g ( d_g − sqrt(d_g² − d_i²) ) ]
   (12)   d   = d_g/2 − (1/2) sqrt(d_g² − d_i²)

Groove cross-section cut by one spherical grain:

   (13)   A' = (d_g²/4)·asin( 2 sqrt(d(d_g−d)) / d_g )
              − sqrt(d(d_g−d)) · (d_g/2 − d)

   (14)   V_a = L_s · A'
   (15)   n_t = (π D W) · C_a          grains engaged per revolution
   (16)   MRR = V_a · n_t · N

Validated against experiment: **R² = 0.991**.

### 2.3 Abrasive pad densities (measured from SEM)

| d_g | areal density C_o | active density C_a |
|---|---|---|
| 6 µm | 1750 mm⁻² | 1400 mm⁻² |
| 15 µm | 720 mm⁻² | 576 mm⁻² |
| 30 µm | 312 mm⁻² | 250 mm⁻² |

Consistent with `C_a = 0.8 C_o`.

### 2.4 The brittle–ductile result — the core of the paper

Chip size measured directly from SEM of the collected chips:

| d_g | chip size | observed mode |
|---|---|---|
| 30 µm | 240–350 nm | plastic flow **+ brittle fracture** |
| 15 µm | 160–230 nm | plastic flow + **fewer** fractures |
| **6 µm** | **60–100 nm** | **pure ductile, no fracture** |

> **d_c (WC-Co) = 60–100 nm**, measured.

The energy-competition argument (their Eq. 17) is the same one `THEORY.md` §2
derives, written for indentation — plastic work goes as d³, fracture as d²,
so the ratio goes as d:

   (17)   G_p/G_f  ~  f(d³, d_g, F_n) / f(d², d_g, F_n)  ~  f(d, d_g, F_n)

A second, purely empirical criterion — the carbide-size ratio:

    k = d_g / D_WC  ;   k = 22, 11, 4 for the three pads
    pure ductile requires  k < 5

### 2.5 FINDING — Bifano's model fails on WC-Co by 17×

The paper applies the same Bifano form our codebase uses as its default:

    d_c = 0.15 · (E/H) · (K_c/H)²

With E = 200 GPa, H = 11.02 GPa, K_c = 7.78 MPa·√m it gives **1.37 µm** —
against a **measured 60–100 nm**. I reproduced this through our own
`semgrit.hybrid.critical_depth_mm(form=2)`: **1.3569 µm, agreeing with the
paper to 0.96 %.** So our implementation is right and the *model* is wrong for
this material, by a factor of **17** against the 80 nm midpoint.

The paper's stated reasons:

1. Thermally sprayed coating ≠ sintered bulk — mechanical bonding between
   splats, porosity, decarburisation.
2. WC-Co is **multiphase**: hard brittle WC in soft ductile Co. A
   single-phase (E, H, K_c) triple cannot represent it.
3. WC grains are anisotropic — prismatic and basal planes have different
   nanohardness.

Note 1.37 µm ≈ D_WC = 1.36 µm. Bifano is returning something close to the
*carbide size*, not the transition depth. Suggestive, and the paper does not
claim more than that.

**Consequences for what we build:**

- We must **not** ship a SAG deck whose `d_c` comes from Bifano on WC-Co.
  Default to the **measured 60–100 nm** and make the form selectable, exactly
  as `THEORY.md` §3 already argues for the three-form disagreement.
- This is a **publishable result in its own right**: an independent
  reimplementation reproducing the failure to 0.96 % is stronger evidence than
  the original single-source claim.
- It strengthens the existing `dc_form` selector rather than contradicting it.

### 2.6 Other results

- Best finish: T = 0.4 mm, N = 1050 rpm.
- Multistep 30 → 15 → 6 µm: Sa 50 → 32 → **21 nm**.
- Residual stress: ground 280 ± 24 MPa → finished **−78 ± 8 MPa** (28 %).
- Max surface temperature **45 °C** (room 21 °C) → thermal damage negligible.
  **This matters for modelling: SAG is essentially isothermal.**
- Lower T and higher N both push toward ductile.

---

## 3. The presentation deck

9 slides, consistent with the paper and adding no new physics. Useful framing:

- *"SAG combines grinding precision with polishing flexibility — the tool
  adapts to the workpiece, not the other way around."*
- *"Compliance is a feature, not a limitation."*
- Contact evolution: **line contact → elastic compression → elliptical area
  contact**, "like Hertzian contact, but compliant".
- Parameter effects: ↑N → more interactions → ↑MRR; ↑T → ↑contact area →
  ↑force, but too much risks brittle damage; ↓d_g → ↓penetration → better Ra.
- Applications: aerospace TBCs, optics, AM post-processing, freeform moulds.

---

## 4. The reference `.inp`

`SAG(input file).inp` — 135 MB, 2 271 691 lines, Abaqus/CAE 2024, job `FINAL_35`.

### 4.1 Architecture

Three parts, all **deformable C3D8R** — note there is **no rigid wheel body**
and **no discrete abrasive**:

| part | nodes | elements | bbox (mm) | material |
|---|---|---|---|---|
| `inner wheel` | 249 | 108 C3D8R | Ø200 × 10 | `ALUMINUM` |
| `outer_ring` | 41 437 | 31 670 C3D8R | r = 100→105, W = 10 | `polyurathane` |
| `workpeice` | 1 110 525 | 880 000 C3D8R + 176 000 C3D8RT | 50 × 12 × 50 | `ALUMINUM` |

Total **1 087 778 elements**.

- Wheel **Ø 210 mm** overall; PU ring is a **5 mm** wall from r = 100 to 105.
- `*Rigid Body, ref node=_PickedSet40, elset=b_Set-24` — the inner wheel is
  *made* rigid by constraint, not by element type.
- `*Tie, name=tie_bt_inner_outer_wheel, adjust=yes` — PU ring tied to hub.
- Workpiece mesh: X 0.25 mm uniform, Y graded 0.30–1.50 mm,
  Z graded **0.00181–0.227 mm** (the fine surface layer).

### 4.2 Materials, verbatim

```
*Material, name=polyurathane
*Density
 2.2e-13,
*Hyperelastic, neo hooke, moduli=LONG TERM
 0.0575, 1.
*Viscoelastic, time=PRONY
 0.11, 0.05, 0.01
```

Neo-Hookean C10 = 0.0575 MPa, D1 = 1 → μ = 2C10 = 0.115 MPa, E ≈ 6C10 =
0.345 MPa. Prony: g1 = 0.11, k1 = 0.05, τ1 = 0.01 s → 89 % of the modulus
retained long-term.

```
*Material, name=ALUMINUM
*Elastic          69000., 0.34, 20.  (to 60000. at 300 C)
*Plastic, hardening=JOHNSON COOK    22., 68., 0.33, 1., 933., 293.
*Rate Dependent, type=JOHNSON COOK  0.015, 2000.
*Density          2.7e-09
*Damage Initiation, criterion=SHEAR / *Damage Evolution, type=DISPLACEMENT 1e-09
*Conductivity, *Specific Heat, *Inelastic Heat Fraction 0.9
```

### 4.3 Steps, contact, BCs

```
*Step, name=for_wheel, nlgeom=YES     *Dynamic, Explicit  , 0.0005
  *Fixed Mass Scaling, factor=50.
  *Boundary, amplitude=Amp-1   Set-27, 2,2, -51600.   (press-in)
*Step, name=grinding, nlgeom=YES      *Dynamic, Explicit  , 0.02
  *Boundary, amplitude=Amp-1   Set-27, 2,2 ; Set-27, 6,6, -500.   (rotate)
```

- **Two steps** — press the wheel in, then rotate. This is exactly the physical
  sequence: establish compression T, then grind.
- `*Contact, op=NEW` + `*Contact Inclusions, ALL EXTERIOR` — **general
  contact**, not contact pairs.
- `*Friction 2.` with `*Surface Behavior, pressure-overclosure=HARD`.
- `*Section Controls EC-1: DISTORTION CONTROL=YES, ELEMENT DELETION=YES,
  hourglass=ENHANCED`.
- `*Boundary Set-7, ENCASTRE` fixes the workpiece.
- ω = 500 rad/s = **4775 rpm**.

### 4.4 FINDINGS on the reference deck

These are defects in the reference, not in our work. Each needs a decision
before it propagates into a generated deck.

**(a) Polyurethane density is wrong by ~5000×.** `2.2e-13` tonne/mm³ =
**0.22 kg/m³**. Real polyurethane is ~1100 kg/m³ (`1.1e-9` tonne/mm³); even
open-cell PU *foam* is 30–100 kg/m³. At 0.22 kg/m³ the ring is lighter than
air. In Abaqus/Explicit density sets the stable increment, so this inflates
Δt for the PU ring and removes its inertia. Almost certainly a deliberate
speed hack, but it is not defensible in a paper.

**(b) The material is aluminium, labelled titanium, and the subject is
WC-Co.** The card is named `ALUMINUM`; the section comments say
`** Section: titanium alloy`; and the paper it references is about WC-Co
coating. The constants (E = 69 GPa, ν = 0.34, ρ = 2.7e-9, T_melt = 933 K,
JC A = 22 MPa) are unambiguously **aluminium** — A = 22 MPa is annealed pure
Al. Ti-6Al-4V would be E = 113.8 GPa, T_melt = 1878 K, A = 1098 MPa. So all
three names disagree and the numbers win.

**(c) There is no abrasive.** No grain geometry, no abrasive layer, no
`*User Material`. The workpiece is cut by *the polyurethane ring itself*
through general contact. So the deck models **wheel compliance only** — it
cannot produce a per-grain chip thickness, and therefore cannot show a
ductile–brittle transition at all. **This is the single biggest gap, and it is
precisely what our existing pipeline supplies.**

**(d) ~10⁷ increments.** c_d(Al) ≈ 6.27e6 mm/s, L_min = 1.81 µm →
Δt = 2.9e-10 s; with mass scaling 50 → 2.0e-9 s over 0.0205 s ≈ **10 million
increments** on 1.09 M elements. That is a very long run. Our
`ELEMENTS_PER_DC` grading exists to spend elements only where they are needed.

**(e) Output requests ask for everything.** `*Element Output` lists ~70
variables including many irrelevant ones (`BURNF`, `IWCONWEP`, `EFABRIC`), and
`*Output, history, variable=ALL`. On 1.09 M elements this produces an enormous
`.odb`.

**(f) Duplicate surfaces.** `m_Surf-1/3/5` are element-for-element identical, as
are `s_Surf-1/3/5`. Harmless, but CAE-generated clutter.

---

## 5. Wider literature

Consistent with the paper and the deck; nothing contradicts them.

- SAG originated at **Cranfield University** (Beaucamp, Namba and co-workers,
  ~2013–2015) for freeform optics and hard ceramics. Elastic tool + bonded
  abrasive, deliberately positioned between grinding and polishing.
- Applied to WC, SiC, CVD-SiC, sapphire, Zerodur, Ti alloys, AM surfaces.
- Ductile-regime removal on nominally brittle materials is the recurring claim,
  and the mechanism is always the same: compliance spreads load, load per grain
  falls below the fracture threshold.
- **Preston's equation**, `MRR ∝ p · v`, is the standard first-order model and
  the paper invokes it to explain MRR rising with N at constant force.
- Hertzian elliptical contact for non-conforming bodies is textbook (Johnson,
  *Contact Mechanics*, 1985).

---

## 6. What this means for the build

### 6.1 What we already have that maps directly

| SAG needs | we already have |
|---|---|
| `d_c` for the workpiece | `semgrit.hybrid.critical_depth_mm`, 3 selectable forms |
| ductile/brittle switch | `vumat_grind.for` — `h < d_c` per material point |
| measured grain shapes | the whole SEM pipeline (this is the differentiator) |
| grain size distributions | `semgrit.measure` — 25 descriptors |
| graded mesh resolving `d_c` | `ELEMENTS_PER_DC = 5` |
| deck writing + verification | `build_deck.py` + 11 `verify_*.py` gates |
| per-stage figures | `semgrit/figures.py` |
| CAD viewer, notebooks | `semgrit/cadviewer.py`, `_make_notebook*.py` |

### 6.2 What is genuinely new

1. **A compliant layer.** Hyperelastic + viscoelastic PU between hub and
   abrasive. New material model path — Abaqus built-in, no VUMAT needed.
2. **`h` from contact, not kinematics.** The rigid model gets `h(u)` from wheel
   geometry and infeed. SAG must get it from Hertzian contact: T → F_N →
   F_N/N_abr → indentation depth d. The chain is paper Eqs. (2)–(12).
3. **Two-step loading.** Press-in to establish T, then rotate.
4. **General contact** rather than the current explicit pairing.
5. **A contact-mechanics preview** — spot area, pressure distribution, active
   grain count, force per grain — the SAG analogue of the chip-thickness
   preview.

### 6.3 The key structural insight

The rigid model computes `h` from **kinematics** (`h(u) = H₀ + H_G u − u²/2r`)
and hands three constants to the VUMAT. SAG computes `h` from **contact
mechanics** (Eqs. 2–12), which yields a *single* `h` per grain, roughly uniform
over the contact patch.

**The VUMAT does not need to change.** It already accepts `h` from a field
variable (`IHMODE = 1`, `PROPS(56)`). So SAG can reuse the hybrid constitutive
law unmodified, supplying `h` per element through a field variable computed
from the contact model. That is a small, well-tested path — not a rewrite.

### 6.4 Open decisions for the user

1. **Workpiece material.** The paper is WC-Co; the reference deck is aluminium
   mislabelled titanium. WC-Co needs a new material card (E = 200 GPa,
   H = 11.02 GPa, K_c = 7.78 MPa·√m, and a **measured** `d_c` of 60–100 nm
   rather than Bifano's 1.37 µm).
2. **`d_c` source.** Given §2.5, default to measured. Keep Bifano available and
   report both, so the 17× discrepancy is visible rather than hidden.
3. **Explicit abrasives, or a smeared abrasive layer?** Explicit measured grains
   are our differentiator but expensive at 250–1750 grains/mm². A hybrid —
   explicit grains in the contact patch, smeared elsewhere — is likely right.
4. **Fix the PU density?** Recommend yes (1.1e-9 tonne/mm³) with mass scaling
   used openly instead of a disguised density hack.

---

## 7. Sources

1. Ghosh, Sidpara & Bandyopadhyay (2021), IJRMHM **99**, 105610. Read in full.
2. `Shape-Adaptive-Grinding-SAG (1).pptx` — 9 slides, read in full.
3. `SAG(input file).inp` — 2 271 691 lines, parsed in full.
4. Literature on SAG (Cranfield lineage), Hertzian contact, Preston's equation.
