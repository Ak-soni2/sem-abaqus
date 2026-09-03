C=======================================================================
C  VUMAT_GRIND2.FOR
C
C  Hybrid ductile / brittle constitutive model for grinding with ANY
C  number of abrasives, Abaqus/Explicit VUMAT.
C
C  DERIVED FROM vumat_grind.for BY _derive_grind2.py. Do not edit the
C  shared body here: edit vumat_grind.for and re-run the derivation, or
C  the two will drift. verify_vumat_grind2.py proves this file is still
C  bit-identical to vumat_grind.for wherever SWMODE = 0.
C
C  It adds a SECOND, purely local criterion that needs no geometry at
C  all, so the switch works for one abrasive or seven hundred, for a
C  plunge or a traverse, and for a second pass over the same groove:
C
C        W_p * L_c  >=  PSI * Kc^2 / E     ->  brittle
C
C  W_p is the accumulated plastic work per unit volume and L_c the
C  element's own characteristic length, which Abaqus hands the VUMAT.
C
C  This is the statement Bifano's critical depth came from, written
C  pointwise -- but the exponents do NOT line up and it would be wrong to
C  claim they do. Plastic work per volume of order H over a depth d,
C  against the fracture energy Kc^2/E, gives
C
C        d_c(energy) = PSI Kc^2/(E H) = PSI (H/E)^1 (Kc/H)^2
C
C  an exponent of +1 on H/E, where the two published geometric forms use
C  +0.5 and -1. The local criterion is a THIRD member of the family.
C
C  So PSI is not defaulted to lambda_c. It is defaulted to the value that
C  makes the local criterion trip at exactly the dc the deck already
C  chose:
C
C        PSI = dc E H / Kc^2      giving      W_p L_c >= H dc
C
C  which reads as plainly as it should: brittle once the plastic work per
C  unit area exceeds the cost of plastically removing a layer of
C  thickness dc at flow stress H. Whichever dc form the card carries, the
C  two criteria then agree by construction. Override PSI to calibrate it
C  against scratch or nanoindentation data instead.
C
C  Unlike the geometric switch, this one triggers on HISTORY, so a point
C  starts ductile and turns brittle as the cut deepens under it -- which
C  is the physical transition, rather than being told in advance where it
C  happens. Once triggered it latches.
C
C  SWMODE, PROPS(57):
C     0  geometric only: h vs dc. Identical to vumat_grind.for.
C     1  energy only: every point starts ductile and flips when the work
C        criterion is met. Needs no chip thickness at all.
C     2  both: brittle if either criterion says so.
C
C  On a flip the JH-2 branch inherits the Johnson-Cook damage as its own.
C  The two damages are different mechanisms with the same meaning -- the
C  fraction of the way to failure -- and carrying it is the only
C  continuation that does not either forgive or double-count the damage
C  the point already has.
C
C  ONE POINT OF PHYSICS TO WATCH. The energy criterion is regularised by
C  L_c, so it is mesh-dependent by construction, exactly as every
C  energy-based failure criterion is. Halving the element size halves the
C  work density needed to trigger. That is the correct behaviour for a
C  fracture-energy criterion and it means PSI is calibrated FOR A MESH.
C  State the element size alongside PSI.
C
C-----------------------------------------------------------------------
C  ORIGINAL HEADER OF vumat_grind.for FOLLOWS
C=======================================================================
C  VUMAT_GRIND.FOR
C
C  Hybrid ductile / brittle constitutive model for SINGLE-GRIT grinding,
C  Abaqus/Explicit VUMAT.
C
C  One material point follows ONE of two laws, chosen once, from the
C  undeformed chip thickness h that the grit takes at that point's
C  station along the scratch:
C
C        h  <  dc   ->  Johnson-Cook with strain-gradient enhancement
C                       (ductile regime: plastic flow, no fracture)
C        h  >= dc   ->  Johnson-Holmquist II
C                       (brittle regime: damage, bulking, chipping)
C
C  dc is the critical depth of cut for the ductile-brittle transition.
C
C-----------------------------------------------------------------------
C  PROVENANCE OF EACH PIECE
C
C  JH-2 branch      verbatim from vumat_jh2.for (this project), which is
C                   itself after
C                     Johnson & Holmquist, AIP Conf. Proc. 309 (1994) [JH94]
C                     Gazonas, ARL-TR-2699 (2002)                    [GAZ]
C                     Cronin et al., 4th Eur. LS-DYNA Conf. (2003)   [CRO]
C                   and reproduces the JH94 bulking benchmark to 1.3 %.
C
C  JC+SGE branch    Johnson & Cook (1983) flow and (1985) damage, with the
C                   strain-gradient (geometrically-necessary-dislocation)
C                   term of
C                     Yadav, Chakladar & Paul, Int. J. Mech. Sci. 231
C                       (2022) 107582, eqs. 8-10          [micro-milling]
C                     Yadav, Das Chakladar & Paul, Int. J. Mach. Tools
C                       Manuf. 194 (2024) 104100, eqs. 24-26  [peening]
C                     Yadav, Jewell, Jones & Ghadbeigi, Int. J. Mech.
C                       Sci. 314 (2026) 111375, eqs. 4-7      [blanking]
C                   All three write the same Taylor/GND hardening; only
C                   the characteristic length changes (uncut chip
C                   thickness / indent diameter / non-local length). Here
C                   it is the local uncut chip thickness h, which is the
C                   micro-milling form, eq. 10, with f_i -> h.
C
C  dc               Bifano, Dow & Scattergood, J. Eng. Ind. 113 (1991)
C                   285-308, and the (H/E)^1/2 variant of the same
C                   indentation-fracture family. BOTH are provided,
C                   selected by PROPS(51), because they differ by
C                   (E/H)^3/2 and therefore carry different lambda_c.
C
C-----------------------------------------------------------------------
C  HOW THE MATERIAL POINT LEARNS h
C
C  A VUMAT sees no kinematics, so h has to be handed to it. With ONE
C  grit the trajectory is known in closed form, so h is a function of
C  the point's position along the scratch and nothing else:
C
C        u  = X . e_t                         (tangential station)
C        h(u) = H0 + HG*u - u^2 / (2*RTIP)
C
C  e_t = (-sin THC, cos THC, 0) is the tangential unit vector at the
C  workpiece centre angle THC, H0 is the chip thickness at the block
C  centre, HG = dh/du is the wedge slope, and the quadratic term is the
C  curvature of the grit's circular path (sagitta), which is tens of
C  nanometres over a 50 um block on a 50 mm wheel and therefore of the
C  same order as dc itself.
C
C  For the plunge deck this project writes, the wedge comes from the
C  radial infeed rather than a table feed:
C
C        HG = -v_r / v_s        H0 = (r_tip - r_ground) + v_r*u0/v_s
C
C  and for a traverse grind the classical kinematic form
C  h(theta) = L_g (v_w/v_s) sin(theta) linearises to the same H0 + HG*u
C  over a block much shorter than the contact arc. semgrit/hybrid.py
C  computes H0, HG and RTIP and writes them into the card, so the
C  Fortran carries no process knowledge of its own.
C
C  h is evaluated ONCE, at the first call, from the undeformed
C  coordinates, and latched in SDV14. Recomputing it later would track
C  the deformed position and drift.
C
C  PROPS(56) selects the source:
C     0  h from the coordinates as above          (single grit, default)
C     1  h from field variable 1                  (per-element map)
C     2  force ductile everywhere                 (verification)
C     3  force brittle everywhere                 (verification, this
C        reproduces vumat_jh2.for exactly)
C
C-----------------------------------------------------------------------
C  SIGN AND ORDERING CONVENTIONS
C    Abaqus stress is positive in tension; pressure P is positive in
C    compression, so P = -mean(stress).
C    Component order, 3D solids:  1=11 2=22 3=33 4=12 5=23 6=13.
C    VUMAT shear strain increments are tensor shear strains.
C    Units: whatever the deck uses, consistently. semgrit writes
C    mm - MPa - tonne - s, so lengths (b, dc, h, H0) are in mm.
C
C-----------------------------------------------------------------------
C  PROPS
C   1..21  JH-2, IDENTICAL to vumat_jh2.for so an existing card is a
C          prefix of this one:
C     1  K1      bulk modulus / EOS coefficient 1            [stress]
C     2  G       shear modulus                               [stress]
C     3  HEL     Hugoniot elastic limit                      [stress]
C     4  PHEL    pressure at the HEL                         [stress]
C     5  T       maximum hydrostatic tensile pressure        [stress]
C     6  A       intact strength coefficient
C     7  B       fractured strength coefficient
C     8  C       strain-rate coefficient
C     9  N       intact strength exponent
C    10  M       fractured strength exponent
C    11  beta    bulking factor, 0..1
C    12  D1      damage coefficient
C    13  D2      damage exponent
C    14  K2      EOS coefficient 2                           [stress]
C    15  K3      EOS coefficient 3                           [stress]
C    16  SFMAX   maximum normalised fractured strength
C    17  SIGHEL  equivalent stress at the HEL                [stress]
C                (<=0 or absent -> 1.5*(HEL-PHEL))
C    18  EDOT0   reference strain rate, default 1.0          [1/s]
C    19  EDMIN   floor on edot*, default 1.0. Use 1.0e-6 for
C                quasi-static work; see vumat_jh2.for header.
C    20  ITCUT   tensile cutoff: 1 (default) P >= -T*(1-D), 0 P >= -T
C    21  FSMAX   JH-2 deletion: >0 delete when EPBAR > FSMAX,
C                <0 delete when D reaches 1 (chip separation),
C                 0 no deletion
C
C  22..40  Johnson-Cook flow and strain-gradient enhancement:
C    22  A_JC    yield stress                                [stress]
C    23  B_JC    hardening modulus                           [stress]
C    24  n_JC    hardening exponent
C    25  C_JC    strain-rate coefficient
C    26  m_JC    thermal-softening exponent
C    27  EDOT0J  reference plastic strain rate               [1/s]
C    28  E       Young's modulus                             [stress]
C    29  NU      Poisson's ratio
C    30  RHO     density                                     [mass/vol]
C    31  CP      specific heat                    [energy/(mass*temp)]
C    32  BETAQ   Taylor-Quinney fraction (0..1)
C    33  T0      reference temperature                       [temp]
C    34  TMELT   melting temperature                         [temp]
C    35  BVEC    Burgers vector b                            [length]
C    36  MTAY    Taylor factor M
C    37  ALPHT   Taylor hardening constant alpha
C    38  LAMSGE  SGE exponent (1.0 reproduces the blanking paper)
C    39  RPRIME  GND coefficient r' (2.0 reproduces the peening paper)
C    40  GSGE    shear modulus used in the SGE term; <=0 -> from E,NU
C
C  41..46  Johnson-Cook damage:
C    41  JD1  42  JD2  43  JD3  44  JD4  45  JD5
C    46  DCRITJ  damage at which the element is deleted (default 1.0)
C
C  47..56  the switch and the chip-thickness field:
C    47  DCUT    critical depth of cut                       [length]
C                <=0 -> computed from 48..51
C    48  LAMC    lambda_c
C    49  HARDN   hardness H                                  [stress]
C    50  KIC     fracture toughness            [stress*sqrt(length)]
C    51  IDCF    1 -> dc = LAMC*(H/E)^0.5*(KIC/H)^2   (slide form)
C                2 -> dc = LAMC*(E/H)    *(KIC/H)^2   (Bifano 1991,
C                     whose calibrated LAMC is 0.15)
C    52  THC     workpiece centre angle theta_c              [rad]
C    53  H0      chip thickness at u = 0                     [length]
C    54  HG      dh/du, the wedge slope                      [-]
C    55  RTIP    grit tip radius; <=0 drops the curvature term [length]
C    56  IHMODE  0 coords, 1 field variable 1, 2 all ductile,
C                3 all brittle
C
C  57..58  the local energy criterion (vumat_grind2 only):
C    57  SWMODE  0 geometric only (this file then equals vumat_grind.for),
C                1 energy only, 2 both
C    58  PSI     calibration constant in W_p L_c >= PSI Kc^2/E.
C                <=0 -> dc*E*H/Kc^2, which makes the local criterion trip
C                at the same dc the geometric one uses; if H or Kc is
C                absent it falls back to LAMC, PROPS(48).
C
C-----------------------------------------------------------------------
C  STATEV  (22; *Depvar, delete=12)
C    1  D       scalar damage, 0..1                       both branches
C    2  EPBAR   accumulated equivalent plastic strain      both
C    3  P       pressure, positive in compression          both
C    4  Q       von Mises equivalent stress                both
C    5  EDOT    equivalent strain rate                     both
C    6  SIGI    normalised intact strength                 JH-2
C    7  SIGF    normalised fractured strength              JH-2
C    8  SIGD    normalised current strength                JH-2
C    9  DEP     equivalent plastic strain increment        both
C   10  MU      volumetric compression rho/rho0 - 1        JH-2
C   11  DELTAP  bulking pressure                           JH-2
C   12  STATUS  1 = active, 0 = deleted                    both
C   13  MODE    1 = ductile JC+SGE, 2 = brittle JH-2
C   14  HLOC    local undeformed chip thickness            [length]
C   15  DCLOC   critical depth of cut actually used        [length]
C   16  TEMP    temperature                                JC
C   17  SJC     Johnson-Cook flow stress, before SGE and   JC
C               before damage
C   18  SEFF    flow stress after SGE, before damage. The  JC
C               stress actually carried is (1-D)*SEFF, and
C               that is what SDV4 holds.
C   19  FSGE    SGE amplification, SEFF/SJC (1 = no size effect)
C   20  INIT    1 once the point has been initialised
C   21  WPLAS   accumulated plastic work per unit volume  [stress]
C   22  ERATIO  W_p L_c E / Kc^2, the energy criterion's own ratio.
C               Reaching PSI is what flips the point. Plot it to see the
C               transition coming.
C
C  Plot SDV13 to see where the transition sits along the scratch, and
C  SDV19 to see how hard the size effect is working.
C
C-----------------------------------------------------------------------
C  VERIFICATION STATUS
C
C  Exercised outside Abaqus by verify_vumat_grind.py, which compiles
C  this file with a single-material-point driver and checks, among
C  others:
C   * with IHMODE=3 every stress component matches vumat_jh2.for to
C     0 ulp over uniaxial, triaxial and load-unload histories;
C   * the JH94 bulking benchmark, 0.559 / 0.711 / 0.646 GPa;
C   * JC uniaxial yield, hardening, rate and thermal terms against the
C     closed-form Johnson-Cook expression;
C   * the SGE factor against sqrt(1 + (r' eta b (M alpha G)^2/s^2)^L);
C   * dc against both published expressions;
C   * h(u) against the closed-form wedge, and mode latching.
C=======================================================================
      subroutine vumat(
C Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     2  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      integer nblock, ndir, nshr, nstatev, nfieldv, nprops
      integer lanneal, ncomp, km, i, it
      integer itcut, idcf, ihmode, imode, iswm
      integer jInfoArray(*)
C
      parameter (i_info_AnnealFlag = 1)
      parameter (z0=0.d0, z1=1.d0, z2=2.d0, z3=3.d0)
      parameter (z4=4.d0, z6=6.d0, half=0.5d0, third=1.d0/3.d0)
      parameter (tiny=1.d-16)
C
      character*80 cmname
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), dtArray(*),
     2  strainInc(nblock,ndir+nshr), relSpinInc(nblock,nshr),
     3  tempOld(nblock), stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     9  defgradNew(nblock,ndir+nshr+nshr),
     1  fieldNew(nblock,nfieldv), stressNew(nblock,ndir+nshr),
     2  stateNew(nblock,nstatev), enerInternNew(nblock),
     3  enerInelasNew(nblock)
C
      dimension sold(6), strl(6), snew(6), de(6)
C
      ncomp = ndir + nshr
      dt = dtArray(1)
      lanneal = jInfoArray(i_info_AnnealFlag)
C
C-----------------------------------------------------------------------
C     JH-2 constants. Table 2 sandstone values are the fallback so that
C     an incomplete card still runs; supply all 17 for any other
C     material.
C-----------------------------------------------------------------------
      rk1   = 3735.6d0
      rg    = 2686.d0
      hel   = 1982.d0
      phel  = 1374.d0
      rt    = 8.d0
      ra    = 0.71d0
      rb    = 0.30d0
      rc    = 0.022d0
      rn    = 0.55d0
      rm    = 0.40d0
      beta  = 1.d0
      rd1   = 0.002d0
      rd2   = 1.20d0
      rk2   = 9000.d0
      rk3   = 22000.d0
      sfmax = 0.25d0
C
      if (nprops .ge. 16) then
        rk1   = props(1)
        rg    = props(2)
        hel   = props(3)
        phel  = props(4)
        rt    = props(5)
        ra    = props(6)
        rb    = props(7)
        rc    = props(8)
        rn    = props(9)
        rm    = props(10)
        beta  = props(11)
        rd1   = props(12)
        rd2   = props(13)
        rk2   = props(14)
        rk3   = props(15)
        sfmax = props(16)
      endif
C
C     sigma_HEL from HEL = PHEL + (2/3)*sigma_HEL   (JH94 eq. 16)
      sighel = 1.5d0*(hel - phel)
      if (nprops .ge. 17) then
        if (props(17) .gt. z0) sighel = props(17)
      endif
C
      edot0 = z1
      if (nprops .ge. 18) then
        if (props(18) .gt. z0) edot0 = props(18)
      endif
C
      edmin = z1
      if (nprops .ge. 19) then
        if (props(19) .gt. z0) edmin = props(19)
      endif
C
      itcut = 1
      if (nprops .ge. 20) then
        if (props(20) .lt. half) itcut = 0
      endif
C
      fsmax = z0
      if (nprops .ge. 21) fsmax = props(21)
C
      tstar = rt/phel
C
C-----------------------------------------------------------------------
C     Johnson-Cook flow, strain-gradient enhancement and JC damage.
C     Defaults are inert: with AJC = 0 the ductile branch yields at
C     zero stress, which is obvious in the output rather than silently
C     plausible, so an incomplete card cannot be mistaken for a
C     calibrated one.
C-----------------------------------------------------------------------
      ajc   = z0
      bjc   = z0
      rnjc  = z1
      cjc   = z0
      rmjc  = z1
      ed0jc = z1
      ejc   = z0
      rnu   = 0.25d0
      rho   = z1
      cp    = z1
      btq   = 0.9d0
      t0k   = 293.15d0
      tmk   = 1873.15d0
      bvec  = z0
      rmtay = z3
      alpht = 0.3d0
      rlam  = z1
      rprim = z2
      gsge  = z0
C
      if (nprops .ge. 22) ajc   = props(22)
      if (nprops .ge. 23) bjc   = props(23)
      if (nprops .ge. 24) rnjc  = props(24)
      if (nprops .ge. 25) cjc   = props(25)
      if (nprops .ge. 26) rmjc  = props(26)
      if (nprops .ge. 27) then
        if (props(27) .gt. z0) ed0jc = props(27)
      endif
      if (nprops .ge. 28) ejc   = props(28)
      if (nprops .ge. 29) rnu   = props(29)
      if (nprops .ge. 30) then
        if (props(30) .gt. z0) rho = props(30)
      endif
      if (nprops .ge. 31) then
        if (props(31) .gt. z0) cp = props(31)
      endif
      if (nprops .ge. 32) btq   = props(32)
      if (nprops .ge. 33) t0k   = props(33)
      if (nprops .ge. 34) tmk   = props(34)
      if (nprops .ge. 35) bvec  = props(35)
      if (nprops .ge. 36) then
        if (props(36) .gt. z0) rmtay = props(36)
      endif
      if (nprops .ge. 37) alpht = props(37)
      if (nprops .ge. 38) then
        if (props(38) .gt. z0) rlam = props(38)
      endif
      if (nprops .ge. 39) rprim = props(39)
      if (nprops .ge. 40) gsge  = props(40)
C
      dj1 = z0
      dj2 = z0
      dj3 = z0
      dj4 = z0
      dj5 = z0
      dcrtj = z1
      if (nprops .ge. 41) dj1 = props(41)
      if (nprops .ge. 42) dj2 = props(42)
      if (nprops .ge. 43) dj3 = props(43)
      if (nprops .ge. 44) dj4 = props(44)
      if (nprops .ge. 45) dj5 = props(45)
      if (nprops .ge. 46) then
        if (props(46) .gt. z0) dcrtj = props(46)
      endif
C
C     Elastic constants of the ductile branch. The JH-2 branch keeps its
C     own K1 and G: the two laws are different descriptions of the same
C     solid and must not be forced to share an elasticity that suits
C     neither.
      if (ejc .gt. z0) then
        gd = ejc/(z2*(z1 + rnu))
        rkd = ejc/(z3*(z1 - z2*rnu))
      else
        gd = rg
        rkd = rk1
      endif
      if (gsge .le. z0) gsge = gd
C
C     Taylor/GND group, constant over the block:
C        sigma_e^2 = sigma_jc^2 * (1 + (SGEC*eta/sigma_jc^2)^LAM)
C     with SGEC = r' * b * (M alpha G)^2      [stress^2 * length]
      sgec = rprim*bvec*(rmtay*alpht*gsge)**2
C
C-----------------------------------------------------------------------
C     The switch.
C-----------------------------------------------------------------------
      dcut  = z0
      rlamc = 0.15d0
      hardn = z0
      rkic  = z0
      idcf  = 1
      thc   = z0
      h0    = z0
      hg    = z0
      rtip  = z0
      ihmode = 0
      if (nprops .ge. 47) dcut  = props(47)
      if (nprops .ge. 48) rlamc = props(48)
      if (nprops .ge. 49) hardn = props(49)
      if (nprops .ge. 50) rkic  = props(50)
      if (nprops .ge. 51) then
        if (props(51) .gt. 1.5d0) idcf = 2
      endif
      if (nprops .ge. 52) thc   = props(52)
      if (nprops .ge. 53) h0    = props(53)
      if (nprops .ge. 54) hg    = props(54)
      if (nprops .ge. 55) rtip  = props(55)
      if (nprops .ge. 56) ihmode = int(props(56) + half)
C
C     The local energy criterion. Off unless asked for, so a card written
C     for vumat_grind.for behaves identically here.
      iswm = 0
      if (nprops .ge. 57) iswm = int(props(57) + half)
      if (iswm .lt. 0 .or. iswm .gt. 2) iswm = 0
      psi = z0
      if (nprops .ge. 58) psi = props(58)
      if (psi .le. z0) then
C       Default: the value that makes the local criterion trip at the same
C       critical depth the geometric one uses, so the two agree instead of
C       being two different thresholds with one name.
        if (dcut .gt. z0 .and. hardn .gt. z0 .and. rkic .gt. z0
     1      .and. ejc .gt. z0) then
          psi = dcut*ejc*hardn/(rkic*rkic)
        else
          psi = rlamc
        endif
      endif
C     Fracture energy per unit area, Kc^2/E, times PSI: the work per area
C     a point has to do before it is allowed to fracture. With the default
C     PSI this is exactly H*dc.
      gcrit = z0
      if (rkic .gt. z0 .and. ejc .gt. z0) gcrit = psi*rkic*rkic/ejc
C
C     dc from hardness and toughness when it was not given directly.
C     Both published forms are offered because they differ by (E/H)^3/2
C     -- roughly 17x on this project's sandstone -- so lambda_c is NOT
C     transferable between them. 0.15 is Bifano's calibration and
C     belongs to form 2.
      if (dcut .le. z0) then
        if (hardn .gt. z0 .and. rkic .gt. z0 .and. ejc .gt. z0) then
          if (idcf .eq. 2) then
            dcut = rlamc*(ejc/hardn)*(rkic/hardn)**2
          else
            dcut = rlamc*sqrt(hardn/ejc)*(rkic/hardn)**2
          endif
        endif
      endif
C
      et1 = -sin(thc)
      et2 =  cos(thc)
C
C-----------------------------------------------------------------------
C     Annealing: reset stress and state.
C-----------------------------------------------------------------------
      if (lanneal .ne. 0) then
        do 20 km = 1, nblock
          do 10 i = 1, ncomp
            stressNew(km,i) = z0
   10     continue
          do 15 i = 1, nstatev
            stateNew(km,i) = z0
   15     continue
          if (nstatev .ge. 12) stateNew(km,12) = z1
          enerInternNew(km) = enerInternOld(km)
          enerInelasNew(km) = enerInelasOld(km)
   20   continue
        return
      endif
C
C-----------------------------------------------------------------------
C     Supported element types: three direct components (3D solids,
C     axisymmetric, plane strain). Anything else passes through.
C-----------------------------------------------------------------------
      if (ndir .ne. 3) then
        do 40 km = 1, nblock
          do 30 i = 1, ncomp
            stressNew(km,i) = stressOld(km,i)
   30     continue
          do 35 i = 1, nstatev
            stateNew(km,i) = stateOld(km,i)
   35     continue
          enerInternNew(km) = enerInternOld(km)
          enerInelasNew(km) = enerInelasOld(km)
   40   continue
        return
      endif
C
C=======================================================================
      do 900 km = 1, nblock
C
C-----------------------------------------------------------------------
C       Carry the old state forward, then overwrite what changes.
C-----------------------------------------------------------------------
        do 45 i = 1, nstatev
          stateNew(km,i) = stateOld(km,i)
   45   continue
C
        dold   = z0
        epold  = z0
        rmuold = z0
        dpold  = z0
        if (nstatev .ge. 1)  dold   = stateOld(km,1)
        if (nstatev .ge. 2)  epold  = stateOld(km,2)
        if (nstatev .ge. 10) rmuold = stateOld(km,10)
        if (nstatev .ge. 11) dpold  = stateOld(km,11)
        if (dold .lt. z0) dold = z0
        if (dold .gt. z1) dold = z1
        if (dpold .lt. z0) dpold = z0
C
C-----------------------------------------------------------------------
C       Decide the branch, once, from the undeformed position.
C
C       Latched in SDV20 rather than recomputed: coordMp follows the
C       deformation, so a point that has been dragged 2 um by the grit
C       would re-read a chip thickness belonging to somewhere else and
C       could change constitutive law mid-scratch.
C-----------------------------------------------------------------------
        iinit = 0
        if (nstatev .ge. 20) iinit = int(stateOld(km,20) + half)
        if (iinit .eq. 0 .or. stepTime .le. z0) then
          if (ihmode .eq. 1) then
            hloc = z0
            if (nfieldv .ge. 1) hloc = fieldOld(km,1)
          else
            u = coordMp(km,1)*et1 + coordMp(km,2)*et2
            hloc = h0 + hg*u
            if (rtip .gt. z0) hloc = hloc - u*u/(z2*rtip)
          endif
          if (hloc .lt. z0) hloc = z0
          imode = 1
          if (hloc .ge. dcut) imode = 2
C         Energy only: start every point ductile and let the work
C         criterion decide. The chip thickness is then never consulted,
C         which is the whole point of the mode.
          if (iswm .eq. 1) imode = 1
          if (ihmode .eq. 2) imode = 1
          if (ihmode .eq. 3) imode = 2
          if (nstatev .ge. 13) stateNew(km,13) = dble(imode)
          if (nstatev .ge. 14) stateNew(km,14) = hloc
          if (nstatev .ge. 15) stateNew(km,15) = dcut
          if (nstatev .ge. 16) stateNew(km,16) = t0k
          if (nstatev .ge. 20) stateNew(km,20) = z1
          if (nstatev .ge. 21) stateNew(km,21) = z0
          if (nstatev .ge. 22) stateNew(km,22) = z0
        else
          imode = 2
          if (nstatev .ge. 13) imode = int(stateOld(km,13) + half)
          hloc = z0
          if (nstatev .ge. 14) hloc = stateOld(km,14)
        endif
        if (imode .ne. 1) imode = 2
C
        tloc = t0k
        if (nstatev .ge. 16) then
          if (stateOld(km,16) .gt. z0) tloc = stateOld(km,16)
        endif
        if (iinit .eq. 0 .or. stepTime .le. z0) tloc = t0k
C
C-----------------------------------------------------------------------
C       Kinematics shared by both branches.
C-----------------------------------------------------------------------
        dvol = strainInc(km,1) + strainInc(km,2) + strainInc(km,3)
        pold = -(stressOld(km,1) + stressOld(km,2)
     1         + stressOld(km,3))*third
        do 50 i = 1, 3
          sold(i) = stressOld(km,i) + pold
          de(i)   = strainInc(km,i) - dvol*third
   50   continue
        do 55 i = 4, ncomp
          sold(i) = stressOld(km,i)
          de(i)   = strainInc(km,i)
   55   continue
C
C       Equivalent deviatoric strain rate, used by both laws.
        een = z0
        do 60 i = 1, 3
          een = een + de(i)*de(i)
   60   continue
        do 65 i = 4, ncomp
          een = een + z2*de(i)*de(i)
   65   continue
        if (dt .gt. z0) then
          edot = sqrt((z2/z3)*een)/dt
        else
          edot = edot0
        endif
C
        if (imode .eq. 2) goto 500
C
C=======================================================================
C       DUCTILE BRANCH: Johnson-Cook + strain-gradient enhancement
C=======================================================================
C
C       Thermal softening, evaluated at the start-of-increment
C       temperature. The heating within one explicit increment is of
C       order 1e-4 K, so resolving it inside the return map would buy
C       nothing and would make the map temperature-implicit.
        thom = (tloc - t0k)/max(tmk - t0k, tiny)
        if (thom .lt. z0) thom = z0
        if (thom .gt. z1) thom = z1
        fthm = z1 - thom**rmjc
        if (fthm .lt. tiny) fthm = tiny
C
C       Gradient length. It cannot fall below a Burgers vector: eta is a
C       dislocation density divided by b, so a shorter length would
C       describe a lattice curvature no lattice can hold. Without this
C       floor the rubbing zone, where h -> 0, returns an infinite flow
C       stress.
        hlen = hloc
        if (hlen .lt. bvec) hlen = bvec
        if (hlen .le. z0) hlen = z1
C
C       Elastic predictor.
        pnew = pold - rkd*dvol
        do 70 i = 1, ncomp
          strl(i) = sold(i) + z2*gd*de(i)
   70   continue
        ss = z0
        do 75 i = 1, 3
          ss = ss + strl(i)*strl(i)
   75   continue
        do 80 i = 4, ncomp
          ss = ss + z2*strl(i)*strl(i)
   80   continue
        qtrl = sqrt(1.5d0*ss)
C
C       Radial return onto the surface carried by the OLD damage.
C
C       Damage degrades the SURFACE, not the stress. Scaling the stress
C       tensor by (1-D) at the end of every increment looks equivalent
C       and is not: the scaled stress is what comes back as stressOld,
C       so after k increments the point carries (1-D)^k, and the
C       per-increment loss soon cancels the elastic increment exactly.
C       The point then parks just under yield, stops accumulating
C       plastic strain, and never fails. Measured on the pure-shear
C       driver case: plastic strain froze at 3.1e-4 and D at 0.0017 for
C       the remaining 39,000 increments. That is the same construction
C       vumat_jc_damage.for uses.
C
C       Degrading the surface also makes the two branches consistent:
C       JH-2 has always moved its strength surface with damage rather
C       than touching the stress.
C
C       sigma_e is monotonically increasing in dep (isotropic hardening,
C       rate hardening and the SGE term all grow with plastic strain),
C       so f(dep) = qtrl - 3G dep - (1-D) sigma_e(dep) decreases
C       monotonically and has a single root in [0, qtrl/3G]. Newton is
C       safeguarded onto that bracket rather than trusted: an
C       unbracketed Newton on a stiff hardening curve can step negative
C       and hand back an elastic answer for a plastic point.
        fdold = z1 - dold
        if (fdold .lt. z0) fdold = z0
        dep = z0
        call grsge(epold, z0, dt, hlen, fthm, ajc, bjc, rnjc, cjc,
     1             ed0jc, sgec, rlam, sy0, sjc0, fsg0, eta0)
        sy   = sy0
        sjc  = sjc0
        fsge = fsg0
        eta  = eta0
        if (qtrl .gt. fdold*sy0) then
          dlo = z0
          dhi = qtrl/(z3*gd)
          dep = min(dhi, max(z0, (qtrl - fdold*sy0)/(z3*gd)))
          rtolq = 1.d-12*max(qtrl, z1)
          do 100 it = 1, 30
            call grsge(epold, dep, dt, hlen, fthm, ajc, bjc, rnjc, cjc,
     1                 ed0jc, sgec, rlam, sy, sjc, fsge, eta)
            f = qtrl - z3*gd*dep - fdold*sy
            if (f .gt. z0) then
              dlo = dep
            else
              dhi = dep
            endif
            if (abs(f) .le. rtolq) goto 110
            ddp = 1.d-7*max(dep, 1.d-8)
            call grsge(epold, dep+ddp, dt, hlen, fthm, ajc, bjc, rnjc,
     1                 cjc, ed0jc, sgec, rlam, sy2, s2, f2, e2)
            dfdd = -z3*gd - fdold*(sy2 - sy)/ddp
            if (dfdd .ge. -tiny) then
              dtry = half*(dlo + dhi)
            else
              dtry = dep - f/dfdd
              if (dtry .le. dlo .or. dtry .ge. dhi)
     1          dtry = half*(dlo + dhi)
            endif
            if (abs(dtry - dep) .le. 1.d-16*max(dep, 1.d-16)) goto 110
            dep = dtry
  100     continue
  110     continue
          call grsge(epold, dep, dt, hlen, fthm, ajc, bjc, rnjc, cjc,
     1               ed0jc, sgec, rlam, sy, sjc, fsge, eta)
        endif
C
        qnew = qtrl - z3*gd*dep
        if (qnew .lt. z0) qnew = z0
        if (qtrl .gt. tiny) then
          scal = qnew/qtrl
        else
          scal = z1
        endif
        do 120 i = 1, ncomp
          snew(i) = strl(i)*scal
  120   continue
        qpl = qnew
C
        epnew = epold + dep
C
C       Adiabatic heating, from the stress that did the plastic work.
C       Explicit in temperature, one increment behind the return map; at
C       dt ~ 1e-10 s that is a 1e-4 K lag.
        tnew = tloc
        if (dep .gt. z0) then
          tnew = tloc + btq*qpl*dep/(max(rho, tiny)*max(cp, tiny))
        endif
C
C       Johnson-Cook damage. The triaxiality is sigma_m/sigma_eq, and
C       sigma_m = -P, so a compressed point has a negative triaxiality
C       and a long failure strain -- which is exactly why the ductile
C       regime survives under a grit and the free surface does not.
        dnew = dold
        if (dep .gt. z0) then
          if (qnew .gt. tiny) then
            trix = -pnew/qnew
          else
            trix = z0
          endif
          if (trix .gt. 1.5d0)  trix = 1.5d0
          if (trix .lt. -1.5d0) trix = -1.5d0
          edr = (dep/max(dt, tiny))/ed0jc
          if (edr .lt. z1) edr = z1
          epf = (dj1 + dj2*exp(dj3*trix))*(z1 + dj4*log(edr))
     1          *(z1 + dj5*thom)
          if (epf .gt. tiny) then
            dnew = dold + dep/epf
          else
            dnew = z1
          endif
          if (dnew .gt. z1) dnew = z1
        endif
        if (dnew .lt. z0) dnew = z0
C
C       Reduce the deviator onto the surface for the NEW damage. This is
C       the softening: the surface falls as D grows and the stress
C       follows it down to zero at D = 1. It mirrors GAZ Table 1 step 14
C       in the brittle branch, and unlike a per-increment (1-D) scaling
C       of the stress it cannot compound, because the target is
C       recomputed from sigma_e each time rather than applied to
C       whatever the point happened to be carrying.
C
C       The pressure is left alone. Johnson-Cook damage is a ductile
C       void mechanism carried by the deviator, and a point that has
C       failed in shear must still hold the grit up: the element is
C       removed by STATUS on this same increment anyway.
        qlim2 = (z1 - dnew)*sy
        if (qlim2 .lt. z0) qlim2 = z0
        if (qnew .gt. qlim2 .and. qnew .gt. tiny) then
          scal = qlim2/qnew
          do 130 i = 1, ncomp
            snew(i) = snew(i)*scal
  130     continue
          qnew = qlim2
        endif
C
        do 140 i = 1, 3
          stressNew(km,i) = snew(i) - pnew
  140   continue
        do 145 i = 4, ncomp
          stressNew(km,i) = snew(i)
  145   continue
C
        if (nstatev .ge. 1)  stateNew(km,1)  = dnew
        if (nstatev .ge. 2)  stateNew(km,2)  = epnew
        if (nstatev .ge. 3)  stateNew(km,3)  = pnew
        if (nstatev .ge. 4)  stateNew(km,4)  = qnew
        if (nstatev .ge. 5)  stateNew(km,5)  = edot
        if (nstatev .ge. 9)  stateNew(km,9)  = dep
        if (nstatev .ge. 16) stateNew(km,16) = tnew
        if (nstatev .ge. 17) stateNew(km,17) = sjc
        if (nstatev .ge. 18) stateNew(km,18) = sy
        if (nstatev .ge. 19) stateNew(km,19) = fsge
C
C       The local energy criterion. Plastic work per unit volume times the
C       element's characteristic length is work per unit area; once it
C       reaches the fracture energy the point is allowed to crack, and
C       from the next increment it follows JH-2.
C
C       Deliberately evaluated AFTER this increment's stress: the point
C       has already done this increment's work under the ductile law, and
C       re-solving it under the brittle one would need the increment
C       repeated. One increment of lag at dt ~ 1e-10 s is nothing.
        wpold = z0
        if (nstatev .ge. 21) wpold = stateOld(km,21)
        if (wpold .lt. z0) wpold = z0
        wpnew = wpold + qpl*dep
        clen = charLength(km)
        if (clen .le. z0) clen = z1
        eratio = z0
        if (gcrit .gt. z0) eratio = wpnew*clen/gcrit
        if (nstatev .ge. 21) stateNew(km,21) = wpnew
        if (nstatev .ge. 22) stateNew(km,22) = eratio
        if (iswm .ne. 0 .and. eratio .ge. z1 .and. nstatev .ge. 13) then
          stateNew(km,13) = z2
        endif
C
C       Deletion. Driven by D alone, which is monotonically
C       non-decreasing, so the flag latches without a sticky read of
C       stateOld -- which would misfire on the first increment, before
C       Abaqus has initialised the SDVs to 1, and delete the model.
        if (nstatev .ge. 12) then
          stateNew(km,12) = z1
          if (dnew .ge. dcrtj) stateNew(km,12) = z0
        endif
C
        qavg = qpl
        goto 800
C
C=======================================================================
C       BRITTLE BRANCH: Johnson-Holmquist II
C       Algorithm and comments as in vumat_jh2.for; the numbering of
C       the steps follows GAZ Table 1.
C=======================================================================
  500   continue
C
C       1. Volumetric compression mu = rho/rho0 - 1 = 1/J - 1.
        u11 = stretchNew(km,1)
        u22 = stretchNew(km,2)
        u33 = stretchNew(km,3)
        if (nshr .eq. 3) then
          u12 = stretchNew(km,4)
          u23 = stretchNew(km,5)
          u13 = stretchNew(km,6)
        else
          u12 = stretchNew(km,4)
          u23 = z0
          u13 = z0
        endif
        detu = u11*(u22*u33 - u23*u23)
     1       - u12*(u12*u33 - u23*u13)
     2       + u13*(u12*u23 - u22*u13)
        if (detu .gt. tiny) then
          rmunew = z1/detu - z1
        else
          rmunew = rmuold - dvol
        endif
C
C       2. Polynomial EOS (JH94 eqs. 7 and 8). K2 and K3 are dropped in
C          tension. The accumulated bulking pressure is always included.
        if (rmunew .ge. z0) then
          pbase = rk1*rmunew + rk2*rmunew*rmunew
     1          + rk3*rmunew*rmunew*rmunew
        else
          pbase = rk1*rmunew
        endif
        pnew = pbase + dpold
C
C       3. Tensile cutoff. GAZ: T* approaches zero as D approaches 1.
        if (itcut .eq. 1) then
          tcut = rt*(z1 - dold)
        else
          tcut = rt
        endif
        if (pnew .lt. -tcut) pnew = -tcut
C
C       4. Elastic deviatoric predictor.
        do 510 i = 1, ncomp
          strl(i) = sold(i) + z2*rg*de(i)
  510   continue
        ss = z0
        do 515 i = 1, 3
          ss = ss + strl(i)*strl(i)
  515   continue
        do 520 i = 4, ncomp
          ss = ss + z2*strl(i)*strl(i)
  520   continue
        qtrl = sqrt(1.5d0*ss)
C
C       5. JH-2 rate factor, on the shared equivalent strain rate.
        eds = edot/edot0
        if (eds .lt. edmin) eds = edmin
        rfac = z1 + rc*log(eds)
        if (rfac .lt. z0) rfac = z0
C
C       6. Strength surfaces at the current pressure, OLD damage.
        pstar = pnew/phel
        pit   = pstar + tstar
        if (pit .lt. tiny) pit = tiny
        pf = pstar
        if (pf .lt. z0) pf = z0
C
C       Optional per-element strength heterogeneity through field
C       variable 1. Not available when field variable 1 is being used to
C       carry the chip thickness.
        het = z1
        if (nfieldv .ge. 1 .and. ihmode .ne. 1) then
          if (fieldOld(km,1) .gt. z0) het = fieldOld(km,1)
        endif
C
        sigis = ra*pit**rn*rfac
        sigfs = rb*pf**rm*rfac
        if (sigfs .gt. sfmax) sigfs = sfmax
        if (sigfs .lt. z0) sigfs = z0
        sigis = sigis*het
        sigfs = sigfs*het
        sigds = sigis - dold*(sigis - sigfs)
        if (sigds .lt. z0) sigds = z0
C
        qlim = sigds*sighel
C
C       7. Radial return in J2 deviatoric space.
        dep = z0
        if (qtrl .gt. qlim .and. qtrl .gt. tiny) then
          scal = qlim/qtrl
          do 525 i = 1, ncomp
            snew(i) = strl(i)*scal
  525     continue
          dep = (qtrl - qlim)/(z3*rg)
          qpre = qlim
        else
          do 530 i = 1, ncomp
            snew(i) = strl(i)
  530     continue
          qpre = qtrl
        endif
C
C       8. Damage accumulation (JH94 eqs. 5 and 6).
        dnew = dold
        if (dep .gt. z0) then
          epf = rd1*pit**rd2
          if (epf .lt. tiny) epf = tiny
          dnew = dold + dep/epf
          if (dnew .gt. z1) dnew = z1
        endif
        if (dnew .lt. z0) dnew = z0
C
C       9. Reduce the deviator onto the surface for the NEW damage
C          (GAZ Table 1 step 14).
        sigd2 = sigis - dnew*(sigis - sigfs)
        if (sigd2 .lt. z0) sigd2 = z0
        qlim2 = sigd2*sighel
        qnew  = qpre
        if (qnew .gt. qlim2 .and. qnew .gt. tiny) then
          scal = qlim2/qnew
          do 535 i = 1, ncomp
            snew(i) = snew(i)*scal
  535     continue
          qnew = qlim2
        endif
C
C       10-11. Bulking, from the deviatoric elastic energy released by
C          the strength drop, U = sigma^2/(6G) (JH94 eqs. 9, 10, 12).
C          Generated in compression only.
        du = (qpre*qpre - qnew*qnew)/(z6*rg)
        if (du .lt. z0) du = z0
        dpnew = dpold
        if (rmunew .gt. z0 .and. du .gt. z0 .and. beta .gt. z0) then
          arg = (rk1*rmunew + dpold)*(rk1*rmunew + dpold)
     1        + z2*beta*rk1*du
          if (arg .lt. z0) arg = z0
          dpnew = -rk1*rmunew + sqrt(arg)
          if (dpnew .lt. dpold) dpnew = dpold
        endif
C
C       12. Final pressure and stress.
        pnew = pbase + dpnew
        if (itcut .eq. 1) then
          tcut = rt*(z1 - dnew)
        else
          tcut = rt
        endif
        if (pnew .lt. -tcut) pnew = -tcut
C
        do 540 i = 1, 3
          stressNew(km,i) = snew(i) - pnew
  540   continue
        do 545 i = 4, ncomp
          stressNew(km,i) = snew(i)
  545   continue
C
C       Reported strengths at the final pressure.
        pstar = pnew/phel
        pit   = pstar + tstar
        if (pit .lt. tiny) pit = tiny
        pf = pstar
        if (pf .lt. z0) pf = z0
        sigis = ra*pit**rn*rfac
        sigfs = rb*pf**rm*rfac
        if (sigfs .gt. sfmax) sigfs = sfmax
        if (sigfs .lt. z0) sigfs = z0
        sigis = sigis*het
        sigfs = sigfs*het
        sigds = sigis - dnew*(sigis - sigfs)
        if (sigds .lt. z0) sigds = z0
C
        epnew = epold + dep
C
        if (nstatev .ge. 1)  stateNew(km,1)  = dnew
        if (nstatev .ge. 2)  stateNew(km,2)  = epnew
        if (nstatev .ge. 3)  stateNew(km,3)  = pnew
        if (nstatev .ge. 4)  stateNew(km,4)  = qnew
        if (nstatev .ge. 5)  stateNew(km,5)  = edot
        if (nstatev .ge. 6)  stateNew(km,6)  = sigis
        if (nstatev .ge. 7)  stateNew(km,7)  = sigfs
        if (nstatev .ge. 8)  stateNew(km,8)  = sigds
        if (nstatev .ge. 9)  stateNew(km,9)  = dep
        if (nstatev .ge. 10) stateNew(km,10) = rmunew
        if (nstatev .ge. 11) stateNew(km,11) = dpnew
C
C       Deletion. Both criteria are driven by monotonically
C       non-decreasing quantities, so the flag latches by itself.
        if (nstatev .ge. 12) then
          stateNew(km,12) = z1
          if (fsmax .gt. z0) then
            if (epnew .gt. fsmax) stateNew(km,12) = z0
          elseif (fsmax .lt. z0) then
            if (dnew .ge. z1) stateNew(km,12) = z0
          endif
        endif
C
        qavg = qlim
C
C=======================================================================
C       Specific internal and inelastic energy, shared.
C=======================================================================
  800   continue
        spow = z0
        do 810 i = 1, 3
          spow = spow + (stressOld(km,i) + stressNew(km,i))
     1                 *strainInc(km,i)
  810   continue
        do 815 i = 4, ncomp
          spow = spow + z2*(stressOld(km,i) + stressNew(km,i))
     1                 *strainInc(km,i)
  815   continue
        spow = half*spow
C
        if (density(km) .gt. tiny) then
          enerInternNew(km) = enerInternOld(km) + spow/density(km)
          enerInelasNew(km) = enerInelasOld(km)
     1                      + qavg*dep/density(km)
        else
          enerInternNew(km) = enerInternOld(km)
          enerInelasNew(km) = enerInelasOld(km)
        endif
C
  900 continue
C
      return
      end
C
C=======================================================================
C  GRSGE -- Johnson-Cook flow stress with the strain-gradient term.
C
C     sigma_jc = [A + B ep^n] [1 + C ln(max(edot/edot0, 1))] * FTHM
C     eta      = 4 ep / hlen
C     sigma_e  = sigma_jc * sqrt(1 + (SGEC*eta/sigma_jc^2)^LAM)
C
C  with SGEC = r' b (M alpha G)^2 assembled by the caller. At LAM = 1
C  this is exactly sqrt(sigma_jc^2 + r' eta b (M alpha G)^2), the
C  blanking paper's eq. 7, and with r' = 2 it is the peening paper's
C  eq. 25 and the micro-milling paper's eq. 8 with Lc = 2b(MaG)^2/s^2.
C
C  The rate ratio is floored at 1 so the logarithm can never soften the
C  material below its quasi-static strength. Johnson-Cook has the same
C  pathology JH-2 has -- a point that stops straining sees edot -> 0 and
C  its yield surface collapses below the stress it is already carrying.
C
C  ep is the plastic strain at the START of the increment and dep the
C  increment being solved for, so the caller can evaluate the same
C  function at trial values without disturbing the state.
C=======================================================================
      subroutine grsge(ep, dep, dt, hlen, fthm, ajc, bjc, rnjc, cjc,
     1                 ed0jc, sgec, rlam, sy, sjc, fsge, eta)
      include 'vaba_param.inc'
      parameter (z0=0.d0, z1=1.d0, z4=4.d0, tiny=1.d-16)
C
      epn = ep + dep
      if (epn .lt. z0) epn = z0
C
      if (epn .gt. z0) then
        hard = ajc + bjc*epn**rnjc
      else
        hard = ajc
      endif
C
      if (dt .gt. z0) then
        edr = (dep/dt)/ed0jc
      else
        edr = z1
      endif
      if (edr .lt. z1) edr = z1
C
      sjc = hard*(z1 + cjc*log(edr))*fthm
      if (sjc .lt. tiny) sjc = tiny
C
      eta = z4*epn/max(hlen, tiny)
      arg = sgec*eta/(sjc*sjc)
      if (arg .lt. z0) arg = z0
      if (arg .gt. z0) then
        fsge = sqrt(z1 + arg**rlam)
      else
        fsge = z1
      endif
      sy = sjc*fsge
      return
      end
