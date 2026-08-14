c=======================================================================
c  Single material point driver for the grinding VUMATs.
c
c  Reads one problem from stdin, calls VUMAT once per increment with
c  nblock = 1, and writes the stress and every state variable to stdout.
c  Nothing here knows anything about the constitutive law, so the same
c  driver exercises vumat_grind.for and the original vumat_jh2.for and
c  the two can be compared component by component.
c
c  stdin:
c     nprops
c     props(1..nprops)
c     nstatev nfieldv
c     field(1..nfieldv)                       (omitted when nfieldv = 0)
c     x y z                                   material point coordinates
c     charlen density
c     dt nseg nout
c     nstep de1 de2 de3 de4 de5 de6           x nseg
c
c  Strain increments are TENSOR components in Abaqus order
c  (11 22 33 12 23 13), matching VUMAT's strainInc.
c
c  stretchNew is built as exp(total logarithmic strain) on the diagonal
c  with the accumulated tensor shears off it, which is the exact right
c  stretch tensor for the coaxial histories these tests use and is what
c  the JH-2 branch reads to get J.
c=======================================================================
      program mpdrv
      implicit real*8(a-h,o-z)
      parameter (mxp=200, mxs=40, mxseg=64)
      dimension props(mxp)
      dimension stateOld(1,mxs), stateNew(1,mxs)
      dimension stressOld(1,6), stressNew(1,6)
      dimension strainInc(1,6), stretchOld(1,6), stretchNew(1,6)
      dimension defgradOld(1,9), defgradNew(1,9)
      dimension fieldOld(1,8), fieldNew(1,8)
      dimension coordMp(1,3), charLength(1), density(1)
      dimension tempOld(1), tempNew(1), relSpinInc(1,3)
      dimension enerInternOld(1), enerInternNew(1)
      dimension enerInelasOld(1), enerInelasNew(1)
      dimension dtArray(2), jInfoArray(8)
      dimension eps(6), dseg(mxseg,6), nsseg(mxseg)
      character*80 cmname
c
      read(*,*) nprops
      read(*,*) (props(i), i = 1, nprops)
      read(*,*) nstatev, nfieldv
      if (nfieldv .gt. 0) read(*,*) (fieldOld(1,i), i = 1, nfieldv)
      read(*,*) coordMp(1,1), coordMp(1,2), coordMp(1,3)
      read(*,*) charLength(1), density(1)
      read(*,*) dt, nseg, nout
      do 10 k = 1, nseg
        read(*,*) nsseg(k), (dseg(k,i), i = 1, 6)
   10 continue
c
      cmname = 'HYBRID'
      ndir  = 3
      nshr  = 3
      ncomp = 6
      dtArray(1) = dt
      dtArray(2) = dt
      do 20 i = 1, 8
        jInfoArray(i) = 0
   20 continue
      do 30 i = 1, 6
        eps(i) = 0.d0
        stressOld(1,i) = 0.d0
        stressNew(1,i) = 0.d0
        strainInc(1,i) = 0.d0
        stretchOld(1,i) = 0.d0
        stretchNew(1,i) = 0.d0
   30 continue
      stretchOld(1,1) = 1.d0
      stretchOld(1,2) = 1.d0
      stretchOld(1,3) = 1.d0
      do 35 i = 1, mxs
        stateOld(1,i) = 0.d0
        stateNew(1,i) = 0.d0
   35 continue
      do 40 i = 1, 9
        defgradOld(1,i) = 0.d0
        defgradNew(1,i) = 0.d0
   40 continue
      defgradOld(1,1) = 1.d0
      defgradOld(1,2) = 1.d0
      defgradOld(1,3) = 1.d0
      defgradNew(1,1) = 1.d0
      defgradNew(1,2) = 1.d0
      defgradNew(1,3) = 1.d0
      do 45 i = 1, 3
        relSpinInc(1,i) = 0.d0
   45 continue
      do 46 i = 1, 8
        fieldNew(1,i) = fieldOld(1,i)
   46 continue
      tempOld(1) = 0.d0
      tempNew(1) = 0.d0
      enerInternOld(1) = 0.d0
      enerInelasOld(1) = 0.d0
      enerInternNew(1) = 0.d0
      enerInelasNew(1) = 0.d0
c
      istep = 0
      totalTime = 0.d0
      do 200 k = 1, nseg
        do 150 is = 1, nsseg(k)
          istep = istep + 1
          stepTime = dble(istep - 1)*dt
          do 100 i = 1, 6
            strainInc(1,i) = dseg(k,i)
            eps(i) = eps(i) + dseg(k,i)
  100     continue
          stretchNew(1,1) = exp(eps(1))
          stretchNew(1,2) = exp(eps(2))
          stretchNew(1,3) = exp(eps(3))
          stretchNew(1,4) = eps(4)
          stretchNew(1,5) = eps(5)
          stretchNew(1,6) = eps(6)
c
          call vumat(
     1      1, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     2      stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     3      props, density, strainInc, relSpinInc,
     4      tempOld, stretchOld, defgradOld, fieldOld,
     5      stressOld, stateOld, enerInternOld, enerInelasOld,
     6      tempNew, stretchNew, defgradNew, fieldNew,
     7      stressNew, stateNew, enerInternNew, enerInelasNew )
c
          totalTime = totalTime + dt
          do 110 i = 1, 6
            stressOld(1,i) = stressNew(1,i)
            stretchOld(1,i) = stretchNew(1,i)
  110     continue
          do 120 i = 1, nstatev
            stateOld(1,i) = stateNew(1,i)
  120     continue
          enerInternOld(1) = enerInternNew(1)
          enerInelasOld(1) = enerInelasNew(1)
c
          if (mod(istep, nout) .eq. 0) then
            write(*,900) istep, totalTime,
     1        (stressNew(1,i), i = 1, 6),
     2        (stateNew(1,i), i = 1, nstatev),
     3        enerInternNew(1), enerInelasNew(1)
          endif
  150   continue
  200 continue
  900 format(i9, 1x, 60(1pe24.16, 1x))
      end
