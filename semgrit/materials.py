"""Workpiece materials a hybrid deck can be built on.

One entry per material, carrying the three things a hybrid deck needs and
nothing else:

  * the 17-constant JH-2 card for the brittle branch,
  * the Johnson-Cook + SGE constants for the ductile branch,
  * hardness and toughness, which is what dc is computed from.

    from semgrit.materials import MATERIALS, hybrid_params
    hp = hybrid_params("silicon_carbide", h_source=0, dc_form=2)
    jh2 = MATERIALS["silicon_carbide"].jh2

UNITS, because this is the file where a unit error is silent and fatal:
the whole project is **mm - tonne - s - MPa - N**. So

    GPa      -> MPa            x 1000
    kg/m^3   -> tonne/mm^3     x 1e-12   (done on the way into the card)
    MPa*m^.5 -> MPa*mm^.5      x sqrt(1000) = 31.623
    J/(kg K) -> mJ/(tonne K)   x 1e6      (done on the way into the card)

Everything stored here is already in MPa and mm except ``density_kg_m3``,
``specific_heat_j_kgk`` and ``kic_mpa_sqrt_m``, which keep their familiar
units and are converted where they are used. ``kic`` inside
:class:`~semgrit.hybrid.HybridParams` is MPa*sqrt(mm), and getting that wrong
scales dc by 1000.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

from .analysis import JH2_SANDSTONE
from .hybrid import HybridParams, kic_from_mpa_sqrt_m

# JH-2 card order, all stresses MPa:
#   1 K1  2 G  3 HEL  4 PHEL  5 T  6 A  7 B  8 C  9 N  10 M
#  11 beta  12 D1  13 D2  14 K2  15 K3  16 SFMAX  17 SIGHEL
#
# SIGHEL = 1.5*(HEL - PHEL) in both cards, which is the JH-2 identity, not a
# coincidence -- it is asserted in _check_jh2_card below.

JH2_SILICON_CARBIDE = (
    204785.0,   # K1      204.785 GPa
    183000.0,   # G       183 GPa
    14457.0,    # HEL     14.457 GPa
    5900.0,     # PHEL    5.9 GPa
    370.0,      # T       0.37 GPa
    0.96,       # A
    0.35,       # B
    0.0,        # C       this card has NO strain-rate term
    0.65,       # N
    1.0,        # M
    1.0,        # beta
    0.48,       # D1
    0.48,       # D2
    0.0,        # K2
    0.0,        # K3
    0.8,        # SFMAX   DIMENSIONLESS, see the note on the entry below
    12835.5,    # SIGHEL  1.5*(14457 - 5900)
)


def quasi_static_ucs_mpa(jh2: Sequence[float]) -> float:
    """Uniaxial compressive strength implied by a JH-2 card, MPa.

    Where the uniaxial-stress elastic path meets the intact surface: at axial
    stress s the pressure is s/3 and the equivalent stress is s, so

        s = A * ((s/3 + T)/PHEL)^N * SIGHEL

    which is monotone on each side, hence the bisection. This is the number
    the ductile branch's Johnson-Cook ``A`` is set to, so the two laws agree in
    uniaxial compression at the transition instead of stepping across it.
    ``verify_vumat_grind.py`` re-derives the same 90.0 MPa for sandstone from
    the subroutine's own surface.
    """
    _, _, hel, phel, t, a, _, _, n, _, _, _, _, _, _, _, sighel = jh2[:17]
    if sighel <= 0:
        sighel = 1.5 * (hel - phel)
    tstar = t / phel
    lo, hi = 1e-9, 1.0e6
    for _ in range(300):
        s = 0.5 * (lo + hi)
        q = a * max(s / 3.0 / phel + tstar, 1e-30) ** n * sighel
        lo, hi = (s, hi) if q > s else (lo, s)
    return 0.5 * (lo + hi)


@dataclass(frozen=True)
class Workpiece:
    """One selectable workpiece material."""

    key: str
    label: str
    inp_material: str
    """Name written on the deck's ``*Material`` line, so an .inp says which
    material it is without opening the report."""
    jh2: tuple
    density_kg_m3: float
    jc: dict
    """Johnson-Cook + SGE constants. ``a_mpa`` is omitted deliberately: it is
    derived from the JH-2 card so the two branches meet."""
    dc: dict
    """``hardness_mpa``, ``kic_mpa_sqrt_m``, and the default ``dc_form``."""
    notes: tuple = ()
    placeholder: bool = True

    def hybrid_params(self, **over) -> HybridParams:
        kw = dict(self.jc)
        kw.setdefault("a_mpa", quasi_static_ucs_mpa(self.jh2))
        kw.update(
            enabled=True,
            density_kg_m3=self.density_kg_m3,
            hardness_mpa=self.dc["hardness_mpa"],
            kic=kic_from_mpa_sqrt_m(self.dc["kic_mpa_sqrt_m"]),
            dc_form=self.dc.get("dc_form", 2),
            lambda_c=self.dc.get("lambda_c", 0.15),
            placeholder=self.placeholder,
        )
        kw.update(over)
        p = HybridParams(**kw)
        p.validate()
        return p

    def dc_nm(self, form: int = 0) -> float:
        p = self.hybrid_params()
        if form:
            p = HybridParams(**{**p.__dict__, "dc_form": form})
        return p.critical_depth_mm() * 1e6


MATERIALS: dict[str, Workpiece] = {}


def _add(w: Workpiece) -> Workpiece:
    _check_jh2_card(w)
    MATERIALS[w.key] = w
    return w


def _check_jh2_card(w: Workpiece) -> None:
    """Refuse a card that is internally inconsistent.

    Cheap, and it catches the two ways a hand-typed card goes wrong: a unit
    slip on one stress (SIGHEL stops matching HEL and PHEL) and SFMAX entered
    as a stress instead of a normalised strength.
    """
    j = w.jh2
    if len(j) != 17:
        raise ValueError("%s: JH-2 card is 17 constants, got %d"
                         % (w.key, len(j)))
    hel, phel, sighel = j[2], j[3], j[16]
    if not math.isclose(sighel, 1.5 * (hel - phel), rel_tol=1e-6):
        raise ValueError(
            "%s: SIGHEL %.6g does not equal 1.5*(HEL-PHEL) = %.6g. One of the "
            "three stresses is in the wrong unit." % (w.key, sighel,
                                                      1.5 * (hel - phel)))
    if not 0.0 < j[15] <= 1.5:
        raise ValueError(
            "%s: SFMAX = %.6g. SFMAX is the maximum NORMALISED fractured "
            "strength, dimensionless and of order 0.2-1.0. A value in MPa "
            "here would make the fractured branch effectively uncapped."
            % (w.key, j[15]))
    if j[4] >= j[3]:
        raise ValueError(
            "%s: T = %.6g MPa is not below PHEL = %.6g MPa, so T* >= 1 and "
            "the intact surface is nonsense." % (w.key, j[4], j[3]))


# ---------------------------------------------------------------------------
# sandstone -- the material every existing deck in this project uses
# ---------------------------------------------------------------------------
_add(Workpiece(
    key="sandstone",
    label="Sandstone (JH-2, the project's original card)",
    inp_material="STONE",
    jh2=JH2_SANDSTONE,
    density_kg_m3=2350.0,
    jc=dict(
        b_mpa=50.0, n=0.50, c=0.020, m=1.0, edot0=1.0,
        youngs_mpa=6500.0, poisson=0.21,
        specific_heat_j_kgk=800.0, taylor_quinney=0.9,
        t0_k=293.15, tmelt_k=1473.15,
        burgers_mm=5.0e-7, taylor_factor=3.0, alpha=0.3, sge_exponent=1.0,
        r_prime=2.0, sge_shear_mpa=0.0,
        d1=0.0, d2=0.15, d3=-1.5, d4=0.0, d5=0.0, dcrit=1.0,
    ),
    dc=dict(hardness_mpa=1000.0, kic_mpa_sqrt_m=0.3, dc_form=2),
    notes=(
        "E = 6500 MPa and nu = 0.21 are exactly this card's K1 and G, so both",
        "branches share one elasticity and the stable increment is",
        "unambiguous. A = 90 MPa is the card's own quasi-static uniaxial",
        "compressive strength. B, n, C, m and D1..D5 are PLACEHOLDERS.",
    ),
))

# ---------------------------------------------------------------------------
# silicon carbide
#
# NAMING. This card was supplied as "monocrystalline silicon", but the numbers
# are the standard Holmquist & Johnson SiC-N set: rho 3163 kg/m^3, G 183 GPa,
# HEL ~14.5 GPa, K1 204.785 GPa, A 0.96, B 0.35, N 0.65, T 0.37 GPa,
# D1 = D2 = 0.48. Monocrystalline silicon is a different material -- rho 2329,
# E ~ 130-188 GPa, H ~ 11 GPa, Kc ~ 0.9 MPa*m^0.5 -- so the entry is keyed and
# labelled as silicon carbide. Say "silicon carbide" in the paper, or replace
# the card.
#
# "Fracture strength 0.8 GPa" in the source table is SFMAX, which is
# DIMENSIONLESS in JH-2. 0.8 is the published SiC value; read as a stress it
# would be 0.8/SIGHEL = 0.062 and the fractured branch would be 13x weaker.
# The GPa in that column is a units-column error, not a different quantity.
#
# E = 450 GPa and nu = 0.16 are the table's. The K1/G pair implies 423 GPa and
# 0.156 instead, a 6% disagreement typical of these tables. The table's values
# are kept -- they are what was supplied -- and the stable increment is taken
# from the stiffer of the two branches, which build_deck.cost_model already
# does.
#
# HARDNESS AND TOUGHNESS are not in the supplied table and are literature
# values for sintered SiC: HV ~ 25 GPa, K_Ic ~ 3.5 MPa*m^0.5. dc follows from
# them, so override them if you have measured your own.
# ---------------------------------------------------------------------------
_add(Workpiece(
    key="silicon_carbide",
    label="Silicon carbide, SiC-N (supplied as 'monocrystalline silicon')",
    inp_material="SIC",
    jh2=JH2_SILICON_CARBIDE,
    density_kg_m3=3163.0,
    jc=dict(
        # A is derived: 8004 MPa, this card's own quasi-static uniaxial
        # compressive strength. Worth noticing that H/3 = 8333 MPa, the
        # indentation-plasticity estimate of the flow stress under a grit, is
        # within 4% of it -- two independent routes to the same yield.
        b_mpa=800.0, n=0.35, c=0.0, m=1.0, edot0=1.0,
        # C = 0 mirrors the JH-2 card, which also has no rate term for SiC.
        youngs_mpa=450000.0, poisson=0.16,
        specific_heat_j_kgk=690.0, taylor_quinney=0.9,
        t0_k=293.15, tmelt_k=3103.15,   # SiC sublimes ~2830 C
        burgers_mm=3.08e-7,             # 3C-SiC a/sqrt(2), a = 0.436 nm
        taylor_factor=3.0, alpha=0.3, sge_exponent=1.0,
        r_prime=2.0, sge_shear_mpa=0.0,
        d1=0.0, d2=0.15, d3=-1.5, d4=0.0, d5=0.0, dcrit=1.0,
    ),
    dc=dict(hardness_mpa=25000.0, kic_mpa_sqrt_m=3.5, dc_form=2),
    notes=(
        "The JH-2 card is the supplied one, unchanged. The Johnson-Cook",
        "constants other than A are PLACEHOLDERS: B, n, m and D1..D5 are",
        "order-of-magnitude values for a covalent ceramic in the ductile",
        "regime. H = 25 GPa and Kc = 3.5 MPa*m^0.5 are literature values for",
        "sintered SiC and are what dc is computed from -- override them if",
        "you have measured your own.",
    ),
))


def get(key: str) -> Workpiece:
    try:
        return MATERIALS[key]
    except KeyError:
        raise KeyError("unknown workpiece material %r; have %s"
                       % (key, ", ".join(sorted(MATERIALS)))) from None


def hybrid_params(key: str, **over) -> HybridParams:
    """:class:`~semgrit.hybrid.HybridParams` for a named material."""
    return get(key).hybrid_params(**over)


def apply(params, key: str, *, check_hybrid: bool = True) -> None:
    """Point a built :class:`~semgrit.build_deck.DeckParams` at a material.

    Mutates in place, because that is how ``depth_of_cut_um`` is already set at
    every call site, and because the JH-2 card, its density and the
    ``*Material`` name have to move together. Set two of the three and the deck
    runs, silently mixing two materials.

    An existing ``analysis.hybrid`` is NOT overwritten -- it may carry hardness
    or toughness you measured -- but it is checked against this material, so a
    ductile branch left on the previous card is an error rather than a result.
    Pass ``check_hybrid=False`` when the ductile constants are deliberately
    hand-entered and are meant to disagree with the registry.
    """
    w = get(key)
    an = params.analysis
    an.jh2_constants = w.jh2
    an.jh2_density_kg_m3 = w.density_kg_m3
    params.wp_material = w.inp_material
    h = getattr(an, "hybrid", None)
    if h is None or not check_hybrid:
        return
    ref = w.hybrid_params()
    for name in ("youngs_mpa", "poisson", "density_kg_m3", "a_mpa"):
        if not math.isclose(getattr(h, name), getattr(ref, name),
                            rel_tol=1e-9):
            raise ValueError(
                "analysis.hybrid.%s is %.6g but %s wants %.6g -- the ductile "
                "branch belongs to a different material than the JH-2 card. "
                "Build it with materials.hybrid_params(%r, ...)."
                % (name, getattr(h, name), key, getattr(ref, name), key))


def summary_text(key: str) -> str:
    w = get(key)
    p = w.hybrid_params()
    L = ["%s  [%s]" % (w.label, w.key), ""]
    L.append("  JH-2   K1 %.6g  G %.6g  HEL %.6g  PHEL %.6g  T %.6g MPa"
             % (w.jh2[0], w.jh2[1], w.jh2[2], w.jh2[3], w.jh2[4]))
    L.append("         A %.4g  B %.4g  C %.4g  N %.4g  M %.4g  SFMAX %.4g"
             % (w.jh2[5], w.jh2[6], w.jh2[7], w.jh2[8], w.jh2[9], w.jh2[15]))
    L.append("         D1 %.4g  D2 %.4g   rho %.6g kg/m3"
             % (w.jh2[11], w.jh2[12], w.density_kg_m3))
    L.append("  JC     A %.6g MPa (derived from the JH-2 card)  B %.6g  n %.4g"
             % (p.a_mpa, p.b_mpa, p.n))
    L.append("         E %.6g MPa  nu %.4g  b %.4g mm"
             % (p.youngs_mpa, p.poisson, p.burgers_mm))
    L.append("  dc     H %.6g MPa  Kc %.6g MPa*sqrt(m)  ->  %.4f nm (form 2)"
             % (p.hardness_mpa, w.dc["kic_mpa_sqrt_m"], w.dc_nm(2)))
    L.append("         %.4f nm with form 1" % w.dc_nm(1))
    if w.notes:
        L.append("")
        L += ["  " + n for n in w.notes]
    return "\n".join(L)


def demo() -> None:
    """Self-check: run ``python -m semgrit.materials``."""
    sand = get("sandstone")
    assert math.isclose(quasi_static_ucs_mpa(sand.jh2), 90.0, rel_tol=2e-4), \
        "sandstone UCS should reproduce the 90.0 MPa in vumat_jh2.for"

    sic = get("silicon_carbide")
    p = sic.hybrid_params()
    # A within 5% of H/3, the indentation-plasticity flow stress. Independent
    # of the bisection, so it is a real check and not a restatement.
    assert abs(p.a_mpa - p.hardness_mpa / 3.0) / (p.hardness_mpa / 3.0) < 0.05

    # dc has to sit where a mesh can resolve it. Both materials land in tens of
    # nanometres on form 2, which is why the same graded mesh serves both.
    for k in MATERIALS:
        assert 10.0 < get(k).dc_nm(2) < 500.0, k

    # Unit tripwires. Each of these is a mistake that produces a deck that
    # runs and is wrong.
    assert math.isclose(p.kic, 3.5 * math.sqrt(1000.0)), "Kc must be MPa*mm^.5"
    props = __import__("semgrit.hybrid", fromlist=["x"]).hybrid_props(
        sic.jh2, p, None)
    from .hybrid import N_HYBRID_PROPS
    assert len(props) == N_HYBRID_PROPS
    assert math.isclose(props[29], 3163.0 * 1e-12), "density -> tonne/mm^3"
    assert math.isclose(props[30], 690.0 * 1e6), "cp -> mJ/(tonne K)"

    # NOTE, deliberately not codified here: T is the one JH-2 constant the
    # source paper (papers/s41598-023-49668-z, Table 2 footnote) flags as
    # mesh-dependent, and the shipped cards carry the 1.0 mm value on 0.03 um
    # elements. See the memory note jh2-mesh-dependent-parameters for the
    # crack-band argument and the numbers. It is NOT applied anywhere, and it
    # should not be until the saturation limit is measured rather than
    # extrapolated -- see the deck headers, which state the element size.

    for k in sorted(MATERIALS):
        print(summary_text(k) + "\n")
    print("materials: ok")


if __name__ == "__main__":
    demo()
