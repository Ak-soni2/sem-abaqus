C===============================================================================
C VUMAT for Johnson-Cook material model with damage, temperature coupling,
C and element deletion (Fixed-form Fortran version for Abaqus/Explicit)
C===============================================================================

      SUBROUTINE VUMAT(
     * STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,
     * DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,
     * TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,
     * NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     * CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,
     * JSTEP,KINC)

      INCLUDE 'VABA_PARAM.INC'

C----- standard Abaqus arguments -----
      CHARACTER*80 CMNAME
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS
      INTEGER NOEL,NPT,LAYER,KSPT,JSTEP,KINC
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV)
      DOUBLE PRECISION DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT
      DOUBLE PRECISION STRAN(NTENS),DSTRAN(NTENS)
      DOUBLE PRECISION TIME(2),DTIME,TEMP,DTEMP
      DOUBLE PRECISION PREDEF(*),DPRED(*)
      DOUBLE PRECISION PROPS(NPROPS)
      DOUBLE PRECISION COORDS(3),DROT(3,3)
      DOUBLE PRECISION PNEWDT,CELENT
      DOUBLE PRECISION DFGRD0(3,3),DFGRD1(3,3)

C----- locals -----
      INTEGER I,J,IERR
      DOUBLE PRECISION A,B,N,CRATE,M,EPS_DOT0,RHO,T0_C,TM_C,BETA
      DOUBLE PRECISION T0K,TMK,TLOC
      DOUBLE PRECISION EPS_P,ACCUM_PW
      DOUBLE PRECISION E,G,KBULK,NU
      DOUBLE PRECISION DEPS(6),SIG_OLD(6),SIG_TRIAL(6)
      DOUBLE PRECISION DEV_TRIAL(6),DEV_NEW(6)
      DOUBLE PRECISION MEAN_STRESS,SEQ_TRIAL,SEQ_NEW,SIGMA_BAR
      DOUBLE PRECISION YIELD_STRESS,RATE_FACTOR,TEMP_FACTOR
      DOUBLE PRECISION DP_GAMMA,H_TANGENT,DSYDEPS
      DOUBLE PRECISION DEPS_P_EQ,EPS_DOT_EQ
      DOUBLE PRECISION CP,K_COND
      DOUBLE PRECISION TPTS(4),EPTS(4),CPPTS(4),KPTS(4)
      INTEGER NTPTS
      PARAMETER (NTPTS=4)
      DOUBLE PRECISION TINY,ONE,TWO,THREE
      DOUBLE PRECISION DTEMP_LOCAL
      INTEGER DEBUG_FLAG
      DOUBLE PRECISION D1,D2,D3,D4,D5,DISP_FAIL,DCRIT
      DOUBLE PRECISION D,DELTA_ACC,DELTA_U_EFF,CELENT_LOC
      INTEGER DEL_FLAG
      DOUBLE PRECISION EPS_F
      DOUBLE PRECISION INTERP_LINEAR
      EXTERNAL INTERP_LINEAR

C----- constants -----
      TINY = 1.0D-16
      ONE = 1.0D0
      TWO = 2.0D0
      THREE = 3.0D0

C----- defaults -----
      A=0.0D0
      B=0.0D0
      N=0.0D0
      CRATE=0.0D0
      M=0.0D0
      EPS_DOT0=1.0D0
      RHO=1.0D0
      T0_C=20.0D0
      TM_C=1600.0D0
      BETA=0.9D0
      DEBUG_FLAG=0
      D1=0.0D0
      D2=0.0D0
      D3=0.0D0
      D4=0.0D0
      D5=0.0D0
      DISP_FAIL=0.0D0
      DCRIT=0.98D0

C----- read PROPS -----
      IF (NPROPS.GE.1)  A=PROPS(1)
      IF (NPROPS.GE.2)  B=PROPS(2)
      IF (NPROPS.GE.3)  N=PROPS(3)
      IF (NPROPS.GE.4)  CRATE=PROPS(4)
      IF (NPROPS.GE.5)  M=PROPS(5)
      IF (NPROPS.GE.6)  EPS_DOT0=PROPS(6)
      IF (NPROPS.GE.9)  RHO=PROPS(9)
      IF (NPROPS.GE.10) T0_C=PROPS(10)
      IF (NPROPS.GE.11) TM_C=PROPS(11)
      IF (NPROPS.GE.12) BETA=PROPS(12)
      IF (NPROPS.GE.14) DEBUG_FLAG=INT(PROPS(14))
      IF (NPROPS.GE.15) D1=PROPS(15)
      IF (NPROPS.GE.16) D2=PROPS(16)
      IF (NPROPS.GE.17) D3=PROPS(17)
      IF (NPROPS.GE.18) D4=PROPS(18)
      IF (NPROPS.GE.19) D5=PROPS(19)
      IF (NPROPS.GE.20) DISP_FAIL=PROPS(20)
      IF (NPROPS.GE.21) DCRIT=PROPS(21)

C----- convert temps -----
      T0K=T0_C+273.15D0
      TMK=TM_C+273.15D0

C----- temp tables -----
      TPTS(1)=293.15D0
      TPTS(2)=373.15D0
      TPTS(3)=573.15D0
      TPTS(4)=773.15D0

      EPTS(1)=109.0D9
      EPTS(2)=91.0D9
      EPTS(3)=75.0D9
      EPTS(4)=75.0D9

      CPPTS(1)=611.0D0
      CPPTS(2)=624.0D0
      CPPTS(3)=674.0D0
      CPPTS(4)=703.0D0

      KPTS(1)=6.8D0
      KPTS(2)=7.4D0
      KPTS(3)=9.8D0
      KPTS(4)=11.8D0

C----- read states -----
      EPS_P=0.0D0
      ACCUM_PW=0.0D0
      TLOC=TEMP
      D=0.0D0
      DELTA_ACC=0.0D0
      DEL_FLAG=0

      IF (NSTATV.GE.1) EPS_P=STATEV(1)
      IF (NSTATV.GE.2) ACCUM_PW=STATEV(2)
      IF (NSTATV.GE.3) THEN
         IF (STATEV(3).GT.0.0D0) TLOC=STATEV(3)
      ENDIF
      IF (NSTATV.GE.4) D=STATEV(4)
      IF (NSTATV.GE.5) DELTA_ACC=STATEV(5)
      IF (NSTATV.GE.6) DEL_FLAG=INT(STATEV(6))
      IF (TLOC.LE.0.0D0) TLOC=T0K

C----- interpolate -----
      E=INTERP_LINEAR(TLOC,TPTS,EPTS,NTPTS)
      CP=INTERP_LINEAR(TLOC,TPTS,CPPTS,NTPTS)
      K_COND=INTERP_LINEAR(TLOC,TPTS,KPTS,NTPTS)

C----- elastic constants -----
      NU=0.34D0
      G=E/(2.0D0*(1.0D0+NU))
      KBULK=E/(3.0D0*(1.0D0-2.0D0*NU))

C----- strain increment -----
      DO I=1,NTENS
         DEPS(I)=DSTRAN(I)
         SIG_OLD(I)=STRESS(I)
         SIG_TRIAL(I)=STRESS(I)
      END DO

C----- elastic predictor -----
      MEAN_STRESS=(DEPS(1)+DEPS(2)+DEPS(3))/3.0D0
      DEV_TRIAL(1)=DEPS(1)-MEAN_STRESS
      DEV_TRIAL(2)=DEPS(2)-MEAN_STRESS
      DEV_TRIAL(3)=DEPS(3)-MEAN_STRESS
      DEV_TRIAL(4)=DEPS(4)
      DEV_TRIAL(5)=DEPS(5)
      DEV_TRIAL(6)=DEPS(6)

      DO I=1,3
         SIG_TRIAL(I)=SIG_TRIAL(I)+2.0D0*G*DEV_TRIAL(I)
     1                 +KBULK*(DEPS(1)+DEPS(2)+DEPS(3))
      END DO
      DO I=4,6
         SIG_TRIAL(I)=SIG_TRIAL(I)+2.0D0*G*DEV_TRIAL(I)
      END DO

C----- deviatoric part -----
      MEAN_STRESS=(SIG_TRIAL(1)+SIG_TRIAL(2)+SIG_TRIAL(3))/3.0D0
      DO I=1,3
         DEV_TRIAL(I)=SIG_TRIAL(I)-MEAN_STRESS
      END DO
      DO I=4,6
         DEV_TRIAL(I)=SIG_TRIAL(I)
      END DO

      SIGMA_BAR=DEV_TRIAL(1)**2+DEV_TRIAL(2)**2+DEV_TRIAL(3)**2+
     1  2.0D0*(DEV_TRIAL(4)**2+DEV_TRIAL(5)**2+DEV_TRIAL(6)**2)
      SEQ_TRIAL=SQRT(1.5D0*SIGMA_BAR+TINY)

C----- temperature factor -----
      TEMP_FACTOR=(TLOC-T0K)/MAX((TMK-T0K),TINY)
      IF (TEMP_FACTOR.LT.0.0D0) TEMP_FACTOR=0.0D0
      IF (TEMP_FACTOR.GT.1.0D0) TEMP_FACTOR=1.0D0

C----- yield stress -----
      RATE_FACTOR=1.0D0
      YIELD_STRESS=(A+B*(MAX(EPS_P,1.0D-20)**N))*(1.0D0-TEMP_FACTOR)**M

C----- elastic/plastic -----
      IF (SEQ_TRIAL.LE.YIELD_STRESS+1.0D-12) THEN
         DO I=1,NTENS
            STRESS(I)=SIG_TRIAL(I)
         END DO
         DTEMP=0.0D0
      ELSE
         IF (DTIME.GT.0.0D0) THEN
            EPS_DOT_EQ=MAX(1.0D-12,EPS_P/DTIME)
         ELSE
            EPS_DOT_EQ=EPS_DOT0
         ENDIF

         RATE_FACTOR=1.0D0+CRATE*LOG(MAX(EPS_DOT_EQ/EPS_DOT0,1.0D-12))
         IF (EPS_P.GT.0.0D0) THEN
            DSYDEPS=B*N*(EPS_P**(N-1.0D0))*RATE_FACTOR*(1.0D0-TEMP_FACTOR)**M
         ELSE
            DSYDEPS=B*N*(1.0D-12**(N-1.0D0))*RATE_FACTOR*(1.0D0-TEMP_FACTOR)**M
         ENDIF
         H_TANGENT=DSYDEPS

         DP_GAMMA=(SEQ_TRIAL-YIELD_STRESS)/(THREE*G+H_TANGENT+TINY)
         IF (DP_GAMMA.LT.0.0D0) DP_GAMMA=0.0D0

         DO I=1,6
            DEV_NEW(I)=DEV_TRIAL(I)*(1.0D0-(THREE*G*DP_GAMMA)/(SEQ_TRIAL+TINY))
         END DO

         SIGMA_BAR=DEV_NEW(1)**2+DEV_NEW(2)**2+DEV_NEW(3)**2+
     1     2.0D0*(DEV_NEW(4)**2+DEV_NEW(5)**2+DEV_NEW(6)**2)
         SEQ_NEW=SQRT(1.5D0*SIGMA_BAR+TINY)

         MEAN_STRESS=(SIG_TRIAL(1)+SIG_TRIAL(2)+SIG_TRIAL(3))/3.0D0
         STRESS(1)=MEAN_STRESS+DEV_NEW(1)
         STRESS(2)=MEAN_STRESS+DEV_NEW(2)
         STRESS(3)=MEAN_STRESS+DEV_NEW(3)
         STRESS(4)=DEV_NEW(4)
         STRESS(5)=DEV_NEW(5)
         STRESS(6)=DEV_NEW(6)

         DEPS_P_EQ=SQRT(3.0D0/2.0D0)*DP_GAMMA
         EPS_P=EPS_P+DEPS_P_EQ
         ACCUM_PW=ACCUM_PW+SEQ_NEW*DEPS_P_EQ/MAX(RHO,TINY)
         DTEMP_LOCAL=BETA*SEQ_NEW*DEPS_P_EQ/(MAX(RHO,TINY)*MAX(CP,TINY))
         TLOC=TLOC+DTEMP_LOCAL
         DTEMP=DTEMP_LOCAL

C----- Damage -----
         CELENT_LOC=MAX(CELENT,1.0D-12)
         DELTA_U_EFF=DEPS_P_EQ*CELENT_LOC
         DELTA_ACC=DELTA_ACC+DELTA_U_EFF

         IF (DISP_FAIL.GT.0.0D0) THEN
            D=D+(DELTA_U_EFF/DISP_FAIL)
         ELSE
            IF (SEQ_NEW.GT.TINY) THEN
               TEMP_FACTOR=(TLOC-T0K)/MAX((TMK-T0K),TINY)
               IF (TEMP_FACTOR.LT.0.0D0) TEMP_FACTOR=0.0D0
               IF (TEMP_FACTOR.GT.1.0D0) TEMP_FACTOR=1.0D0
               IF (ABS(SEQ_NEW).GT.TINY) THEN
                  SIGMA_BAR=MEAN_STRESS/SEQ_NEW
               ELSE
                  SIGMA_BAR=0.0D0
               ENDIF
               IF (EPS_DOT_EQ.LE.0.0D0) EPS_DOT_EQ=EPS_DOT0
               EPS_F=(D1+D2*EXP(D3*SIGMA_BAR))*
     1               (1.0D0+D4*LOG(MAX(EPS_DOT_EQ/EPS_DOT0,1.0D-12)))*
     2               (1.0D0+D5*TEMP_FACTOR)
               IF (EPS_F.LE.0.0D0) EPS_F=1.0D-12
               D=D+DEPS_P_EQ/EPS_F
            ENDIF
         ENDIF

         IF (D.LT.0.0D0) D=0.0D0
         IF (D.GT.1.0D0) D=1.0D0

         DO I=1,NTENS
            STRESS(I)=STRESS(I)*(1.0D0-D)
         END DO

         IF (D.GE.DCRIT) THEN
            DO I=1,NTENS
               STRESS(I)=1.0D-12
            END DO
            DEL_FLAG=1
         ENDIF

         IF (NSTATV.GE.1) STATEV(1)=EPS_P
         IF (NSTATV.GE.2) STATEV(2)=ACCUM_PW
         IF (NSTATV.GE.3) STATEV(3)=TLOC
         IF (NSTATV.GE.4) STATEV(4)=D
         IF (NSTATV.GE.5) STATEV(5)=DELTA_ACC
         IF (NSTATV.GE.6) STATEV(6)=DBLE(DEL_FLAG)
      ENDIF

C----- zero DDSDDE -----
      DO I=1,NTENS
         DO J=1,NTENS
            DDSDDE(I,J)=0.0D0
         END DO
      END DO

C----- debug output for single-element tests (disable in parallel) -----
      IF (DEBUG_FLAG .EQ. 1) THEN
         OPEN(UNIT=77, FILE='vumat_debug.log', STATUS='UNKNOWN', IOSTAT=IERR)
         IF (IERR .EQ. 0) THEN
            WRITE(77,*) 'NOEL=', NOEL, ' NPT=', NPT, 'STEP=', JSTEP, 'KINC=', KINC
            WRITE(77,'(A,6E15.6)') 'STRESS=', STRESS(1),STRESS(2),STRESS(3),
     &                              STRESS(4),STRESS(5),STRESS(6)
            WRITE(77,'(A,3E12.6)') 'eps_p, D, Tloc =', EPS_P, D, TLOC
            CLOSE(77)
         END IF
      END IF

      RETURN
      END SUBROUTINE VUMAT


C-------------------------------------------------------------------------------
C Linear interpolation helper (returns double precision)
C-------------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION INTERP_LINEAR(T, TARRAY, YARRAY, NPTS)
      DOUBLE PRECISION T
      DOUBLE PRECISION TARRAY(*), YARRAY(*)
      INTEGER NPTS
      INTEGER II
      IF (T .LE. TARRAY(1)) THEN
         INTERP_LINEAR = YARRAY(1)
         RETURN
      ELSEIF (T .GE. TARRAY(NPTS)) THEN
         INTERP_LINEAR = YARRAY(NPTS)
         RETURN
      END IF
      DO II = 1, NPTS - 1
         IF (T .GE. TARRAY(II) .AND. T .LT. TARRAY(II+1)) THEN
            INTERP_LINEAR = YARRAY(II) + (YARRAY(II+1)-YARRAY(II))
     &         * ((T - TARRAY(II)) / (TARRAY(II+1) - TARRAY(II)))
            RETURN
         END IF
      END DO
      INTERP_LINEAR = YARRAY(NPTS)
      RETURN
      END
