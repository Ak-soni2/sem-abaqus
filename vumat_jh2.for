C=======================================================================
C  VUMAT_JH2.FOR
C
C  Johnson-Holmquist II (JH-2) constitutive model for brittle materials,
C  Abaqus/Explicit VUMAT.  Parameters below are the sandstone set of
C
C    Baranowski, Kucewicz & Janiszewski, "JH-2 constitutive model of
C    sandstone for dynamic problems", Sci. Rep. 14, 3339 (2024).
C
C  Model equations follow the primary sources:
C    Johnson & Holmquist, AIP Conf. Proc. 309, 981-984 (1994)   [JH94]
C    Gazonas, ARL-TR-2699 (2002)                                [GAZ]
C    Cronin, Bui, Kaufmann, McIntosh, Berstad, 4th European
C      LS-DYNA Users Conf. (2003), D-I-47..60                   [CRO]
C      (this is the LS-DYNA *MAT_110 implementation note, i.e. the
C       code that produced the reference paper's results)
C
C-----------------------------------------------------------------------
C  VERIFICATION STATUS
C
C  This routine was compiled and exercised outside Abaqus through a
C  single material point driver.  Results:
C
C  1. Johnson & Holmquist (1994) 1 m cube, uniaxial strain, load to 5%
C     then unload.  Published bulking pressures vs. this routine:
C         Case A   published 0.56 GPa   computed 0.559 GPa   (-0.11 %)
C         Case B   published 0.72 GPa   computed 0.711 GPa   (-1.30 %)
C         Case C   published 0.65 GPa   computed 0.646 GPa   (-0.61 %)
C     (published values are quoted to two significant figures)
C
C  2. Closed-form intersection of the elastic path with the intact
C     strength surface, sandstone parameters, rate factor 1:
C         uniaxial compression      exact  90.000   VUMAT  89.968 MPa
C         uniaxial tension          exact  17.932   VUMAT  17.931 MPa
C         triaxial, 10 MPa confine  exact 119.886   VUMAT 119.876 MPa
C         triaxial, 17 MPa confine  exact 138.380   VUMAT 138.353 MPa
C         triaxial, 25 MPa confine  exact 158.049   VUMAT 158.017 MPa
C
C  3. Fully fractured (D=1) residual, triaxial: 57.6 / 59.8 / 61.5 % of
C     peak at 10 / 17 / 25 MPa confinement.  The reference paper reports
C     60 % for the single-element triaxial test.
C
C-----------------------------------------------------------------------
C  IMPORTANT - THE STRAIN RATE FLOOR CONTROLS QUASI-STATIC STRENGTH
C
C  JH-2 scales both strength surfaces by (1 + C*ln(edot*)), with
C  edot* = edot/edot0 and edot0 = 1.0 1/s.  Below the reference rate this
C  factor is LESS than one, so the model softens in quasi-static loading.
C  Neither JH94, GAZ nor CRO imposes any floor on edot*.
C
C  With Table 2 sandstone constants the single-element uniaxial
C  compressive strength is
C        edot* floored at 1.0   ->  90.0 MPa
C        edot ~ 0.04 1/s        ->  80.4 MPa  (measured UCS is 79.1 MPa)
C
C  So the paper's quasi-static results are only reproducible with the
C  floor OFF.  But leaving it off is dangerous in a dynamic model: a
C  material point that is loaded and then goes still sees edot -> 0, the
C  rate factor collapses, and the yield surface falls below the stress
C  the point is already carrying.  Measured with this routine, a point
C  held at 70.6 MPa and then given 2000 zero-strain increments decays to
C  62.8 MPa and reaches D = 1.0 without any real deformation.
C
C  DEFAULT IS THEREFORE THE SAFE CLAMP, PROPS(19) = 1.0.
C
C  Choose per analysis:
C    quasi-static UC / UT / TXC   PROPS(19) = 1.0e-6   (reproduces the
C                                 paper; every point strains at the same
C                                 rate, so nothing goes spuriously quiet)
C    SHPB, drop weight, blast,    PROPS(19) = 1.0      (leave the default;
C    projectile impact            edot* >> 1 throughout, so the clamp is
C                                 inert and results are unaffected)
C
C-----------------------------------------------------------------------
C  SIGN AND ORDERING CONVENTIONS
C    Abaqus stress is positive in tension; JH-2 pressure P is positive
C    in compression, so P = -mean(stress).
C    Component order, 3D solids:  1=11 2=22 3=33 4=12 5=23 6=13.
C    Axisymmetric / plane strain: 1=11 2=22 3=33 4=12  (nshr = 1).
C    VUMAT shear strain increments are tensor shear strains.
C
C-----------------------------------------------------------------------
C  MATERIAL CARD  (mm - MPa - tonne - s)
C
C    *Density
C     2.35e-9
C    *Depvar
C     11
C    *User Material, constants=17
C     3735.6, 2686.0, 1982.0, 1374.0, 8.0, 0.71, 0.30, 0.022
C     0.55, 0.40, 1.0, 0.002, 1.20, 9000.0, 22000.0, 0.25
C     912.0
C
C  To also enable element deletion, use 21 constants and *Depvar 12
C  with the delete attribute:
C
C    *Depvar, delete=12
C     12
C    *User Material, constants=21
C     3735.6, 2686.0, 1982.0, 1374.0, 8.0, 0.71, 0.30, 0.022
C     0.55, 0.40, 1.0, 0.002, 1.20, 9000.0, 22000.0, 0.25
C     912.0, 1.0, 1.0, 1.0, 0.5
C
C  For the quasi-static UC / UT / TXC tests replace the 19th constant
C  (1.0) with 1.0e-6.
C
C-----------------------------------------------------------------------
C  PROPS
C    1  K1      bulk modulus / EOS pressure coefficient 1   [stress]
C    2  G       shear modulus                               [stress]
C    3  HEL     Hugoniot elastic limit                      [stress]
C    4  PHEL    pressure at the HEL                         [stress]
C    5  T       maximum hydrostatic tensile pressure        [stress]
C    6  A       intact strength coefficient
C    7  B       fractured strength coefficient
C    8  C       strain-rate coefficient
C    9  N       intact strength exponent
C    10 M       fractured strength exponent
C    11 beta    bulking factor, 0..1
C    12 D1      damage coefficient
C    13 D2      damage exponent
C    14 K2      EOS pressure coefficient 2                  [stress]
C    15 K3      EOS pressure coefficient 3                  [stress]
C    16 SFMAX   maximum normalized fractured strength
C    17 SIGHEL  equivalent stress at the HEL                [stress]
C               (<=0 or absent -> 1.5*(HEL-PHEL))
C  optional:
C    18 EDOT0   reference strain rate, default 1.0          [1/s]
C    19 EDMIN   floor on edot*, DEFAULT 1.0 (rate factor clamped at >= 1).
C               Set 1.0e-6 for quasi-static problems (see note above).
C    20 ITCUT   tensile cutoff mode:
C                 1 (default) P >= -T*(1-D)   per GAZ / JH94 Fig. 1
C                 0           P >= -T         constant cutoff
C    21 FSMAX   element deletion.  Requires *Depvar,delete=12 and 12 SDVs.
C                 > 0  delete when EPBAR exceeds this value
C                 < 0  delete when D reaches 1  (CHIP SEPARATION - use
C                      this for machining / orthogonal cutting)
C                 = 0  or absent: no deletion
C               Erosion is NOT part of the constitutive law (CRO).
C
C  FIELD VARIABLE 1 (optional) - strength heterogeneity multiplier
C    Multiplies the intact and fractured strength surfaces.  Absent, or
C    <= 0, means 1.0, i.e. the homogeneous material of the paper.
C
C    Why you may want it: a homogeneous specimen under FRICTIONLESS
C    compression has an exactly uniform stress field, so every point
C    reaches the yield surface on the same increment.  At zero
C    confinement the fractured strength B(P*)^M is zero, so softening
C    dissipates ~0.0002 N-mm/mm^3 (Gf ~ 0.2 J/m^2 at 1 mm, against
C    10-50 J/m^2 for real sandstone).  Zero fracture energy means no
C    length scale, so the crack pattern is set by round-off and mesh
C    topology.  A few percent of Weibull scatter gives localisation a
C    physical seed and makes the failure mode statistically meaningful.
C
C    Assign with *INITIAL CONDITIONS, TYPE=FIELD, VARIABLE=1 and
C    *USER DEFINED FIELD, or generate an element-wise table.  Note this
C    is a deliberate departure from the reference paper, which used a
C    homogeneous specimen.
C
C  STATEV
C    1  D       scalar damage, 0..1
C    2  EPBAR   accumulated equivalent plastic strain
C    3  P       pressure, positive in compression
C    4  Q       von Mises equivalent stress
C    5  EDOT    equivalent strain rate
C    6  SIGI    normalized intact strength
C    7  SIGF    normalized fractured strength
C    8  SIGD    normalized current (damaged) strength
C    9  DEP     equivalent plastic strain increment
C    10 MU      volumetric compression, rho/rho0 - 1
C    11 DELTAP  bulking pressure   (LS-DYNA history variable 1)
C    12 STATUS  1 = active, 0 = deleted (only if FSMAX > 0)
C
C-----------------------------------------------------------------------
C  ALGORITHM, one increment, after CRO Table/narrative and GAZ Table 1
C    1  mu from det(U)
C    2  polynomial EOS + accumulated bulking pressure
C    3  tensile cutoff
C    4  elastic deviatoric predictor
C    5  equivalent strain rate
C    6  intact / fractured / current strength using the OLD damage
C    7  radial return -> equivalent plastic strain increment
C    8  damage update
C    9  reduce the deviator onto the NEW damage surface  (GAZ step 14)
C    10 energy released by that reduction, dU = (q_old^2-q_new^2)/(6G)
C    11 bulking pressure, JH94 eq. 12
C    12 final pressure, tensile cutoff, assemble stress
C
C  Note on ordering: GAZ (Table 1 steps 12-13) applies the new bulking
C  pressure within the increment; CRO (eq. 2a) lags it one cycle.  The
C  sources disagree; this routine follows GAZ, which is the more accurate
C  of the two.  With explicit time increments the difference is one
C  cycle and is numerically irrelevant - both reproduce the JH94
C  benchmark above.
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
      integer lanneal, ncomp, km, i
      integer jInfoArray(*)
C
      parameter (i_info_AnnealFlag = 1)
      parameter (z0=0.d0, z1=1.d0, z2=2.d0, z3=3.d0)
      parameter (z6=6.d0, half=0.5d0, third=1.d0/3.d0)
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
C     Material constants.  Table 2 sandstone values are the fallback so
C     that an incomplete card still runs; supply all 17 for any other
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
C     Floor on the dimensionless strain rate.  Default clamps the rate
C     factor at >= 1; lower it for quasi-static work (see header).
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
C     axisymmetric, plane strain).  Anything else passes through.
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
      do 200 km = 1, nblock
C
C-----------------------------------------------------------------------
C       Old state.
C-----------------------------------------------------------------------
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
C       The deletion flag is derived from EPBAR alone (below).  EPBAR is
C       monotonically non-decreasing, so the flag latches without needing
C       a sticky read of stateOld, which would misfire on the first
C       increment when Abaqus has not yet initialised the SDVs to 1.
C
C-----------------------------------------------------------------------
C       1. Volumetric compression mu = rho/rho0 - 1 = 1/J - 1, with J
C          from the right stretch tensor.  Fall back to the incremental
C          volume change if U is unavailable.
C-----------------------------------------------------------------------
        dvol = strainInc(km,1) + strainInc(km,2) + strainInc(km,3)
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
C-----------------------------------------------------------------------
C       2. Polynomial EOS (JH94 eqs. 7 and 8).  K2 and K3 are dropped in
C          tension.  The accumulated bulking pressure is always included.
C-----------------------------------------------------------------------
        if (rmunew .ge. z0) then
          pbase = rk1*rmunew + rk2*rmunew*rmunew
     1          + rk3*rmunew*rmunew*rmunew
        else
          pbase = rk1*rmunew
        endif
        pnew = pbase + dpold
C
C       3. Tensile cutoff.  GAZ: "T* approaches zero as D approaches 1".
        if (itcut .eq. 1) then
          tcut = rt*(z1 - dold)
        else
          tcut = rt
        endif
        if (pnew .lt. -tcut) pnew = -tcut
C
C-----------------------------------------------------------------------
C       4. Elastic deviatoric predictor.
C-----------------------------------------------------------------------
        pold = -(stressOld(km,1) + stressOld(km,2)
     1         + stressOld(km,3))*third
        do 60 i = 1, 3
          sold(i) = stressOld(km,i) + pold
          de(i)   = strainInc(km,i) - dvol*third
   60   continue
        do 65 i = 4, ncomp
          sold(i) = stressOld(km,i)
          de(i)   = strainInc(km,i)
   65   continue
C
        do 70 i = 1, ncomp
          strl(i) = sold(i) + z2*rg*de(i)
   70   continue
C
        ss = z0
        do 75 i = 1, 3
          ss = ss + strl(i)*strl(i)
   75   continue
        do 80 i = 4, ncomp
          ss = ss + z2*strl(i)*strl(i)
   80   continue
        qtrl = sqrt(1.5d0*ss)
C
C-----------------------------------------------------------------------
C       5. Equivalent (deviatoric) strain rate and the JH-2 rate factor.
C-----------------------------------------------------------------------
        een = z0
        do 85 i = 1, 3
          een = een + de(i)*de(i)
   85   continue
        do 90 i = 4, ncomp
          een = een + z2*de(i)*de(i)
   90   continue
        if (dt .gt. z0) then
          edot = sqrt((z2/z3)*een)/dt
        else
          edot = edot0
        endif
C
        eds = edot/edot0
        if (eds .lt. edmin) eds = edmin
        rfac = z1 + rc*log(eds)
        if (rfac .lt. z0) rfac = z0
C
C-----------------------------------------------------------------------
C       6. Strength surfaces at the current pressure, OLD damage
C          (JH94 eqs. 1, 3, 4).
C-----------------------------------------------------------------------
        pstar = pnew/phel
        pit   = pstar + tstar
        if (pit .lt. tiny) pit = tiny
        pf = pstar
        if (pf .lt. z0) pf = z0
C
C       Optional per-element strength heterogeneity.  Field variable 1
C       multiplies both strength surfaces.  Assign a Weibull-distributed
C       field with *INITIAL CONDITIONS, TYPE=FIELD, VARIABLE=1.
C       Absent or non-positive -> 1.0, i.e. the homogeneous material of
C       the reference paper, bit-for-bit.
C
C       A homogeneous specimen under frictionless compression has an
C       EXACTLY uniform stress field, so every point reaches the yield
C       surface on the same increment and the crack pattern is decided
C       by round-off rather than mechanics.  Scatter of a few percent
C       gives localisation a physical seed.
        het = z1
        if (nfieldv .ge. 1) then
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
C-----------------------------------------------------------------------
C       7. Radial return in J2 deviatoric space.
C-----------------------------------------------------------------------
        dep = z0
        if (qtrl .gt. qlim .and. qtrl .gt. tiny) then
          scal = qlim/qtrl
          do 95 i = 1, ncomp
            snew(i) = strl(i)*scal
   95     continue
          dep = (qtrl - qlim)/(z3*rg)
          qpre = qlim
        else
          do 100 i = 1, ncomp
            snew(i) = strl(i)
  100     continue
          qpre = qtrl
        endif
C
C-----------------------------------------------------------------------
C       8. Damage accumulation (JH94 eqs. 5 and 6).  epf -> 0 as
C          P* -> -T*, where no plastic strain is admissible, so damage
C          saturates immediately there.
C-----------------------------------------------------------------------
        dnew = dold
        if (dep .gt. z0) then
          epf = rd1*pit**rd2
          if (epf .lt. tiny) epf = tiny
          dnew = dold + dep/epf
          if (dnew .gt. z1) dnew = z1
        endif
        if (dnew .lt. z0) dnew = z0
C
C-----------------------------------------------------------------------
C       9. Reduce the deviator onto the surface for the new damage
C          (GAZ Table 1 step 14).
C-----------------------------------------------------------------------
        sigd2 = sigis - dnew*(sigis - sigfs)
        if (sigd2 .lt. z0) sigd2 = z0
        qlim2 = sigd2*sighel
        qnew  = qpre
        if (qnew .gt. qlim2 .and. qnew .gt. tiny) then
          scal = qlim2/qnew
          do 105 i = 1, ncomp
            snew(i) = snew(i)*scal
  105     continue
          qnew = qlim2
        endif
C
C-----------------------------------------------------------------------
C       10-11. Bulking.  dU is the deviatoric elastic energy released by
C          the strength drop, U = sigma^2/(6G)  (JH94 eqs. 9, 10), and
C          the pressure increment follows from energy conservation
C          (JH94 eq. 12):
C            dP_new = -K1*mu + sqrt((K1*mu + dP_old)^2 + 2*beta*K1*dU)
C          Bulking is generated only in compression (JH94, JHB).
C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
C       12. Final pressure and stress.
C-----------------------------------------------------------------------
        pnew = pbase + dpnew
        if (itcut .eq. 1) then
          tcut = rt*(z1 - dnew)
        else
          tcut = rt
        endif
        if (pnew .lt. -tcut) pnew = -tcut
C
        do 110 i = 1, 3
          stressNew(km,i) = snew(i) - pnew
  110   continue
        do 115 i = 4, ncomp
          stressNew(km,i) = snew(i)
  115   continue
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
C-----------------------------------------------------------------------
C       State variables.
C-----------------------------------------------------------------------
        do 120 i = 1, nstatev
          stateNew(km,i) = stateOld(km,i)
  120   continue
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
C       Deletion flag.  Both criteria are driven by monotonically
C       non-decreasing quantities (EPBAR and D only ever grow), so the
C       flag latches by itself without reading stateOld back - which
C       would misfire on the first increment, before Abaqus has
C       initialised the SDVs to 1, and delete the whole model.
C         FSMAX > 0 : delete when EPBAR exceeds FSMAX
C         FSMAX < 0 : delete when D reaches 1  (chip separation; this is
C                     the criterion needed for machining / cutting, and
C                     the one used in vumat_jh2_3)
C         FSMAX = 0 : no deletion
        if (nstatev .ge. 12) then
          stateNew(km,12) = z1
          if (fsmax .gt. z0) then
            if (epnew .gt. fsmax) stateNew(km,12) = z0
          elseif (fsmax .lt. z0) then
            if (dnew .ge. z1) stateNew(km,12) = z0
          endif
        endif
C
C-----------------------------------------------------------------------
C       Specific internal and inelastic energy.
C-----------------------------------------------------------------------
        spow = z0
        do 125 i = 1, 3
          spow = spow + (stressOld(km,i) + stressNew(km,i))
     1                 *strainInc(km,i)
  125   continue
        do 130 i = 4, ncomp
          spow = spow + z2*(stressOld(km,i) + stressNew(km,i))
     1                 *strainInc(km,i)
  130   continue
        spow = half*spow
C
        if (density(km) .gt. tiny) then
          enerInternNew(km) = enerInternOld(km) + spow/density(km)
          enerInelasNew(km) = enerInelasOld(km)
     1                      + qlim*dep/density(km)
        else
          enerInternNew(km) = enerInternOld(km)
          enerInelasNew(km) = enerInelasOld(km)
        endif
C
  200 continue
C
      return
      end
