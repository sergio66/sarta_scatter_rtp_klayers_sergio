C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    YCALT5 (set5=FWO sun bfsw) version with trace gases (no SO2 or HNO3)
C
!F77====================================================================


!ROUTINE NAME:
C    YCALT5


!ABSTRACT:
C    Calculate the transmittance for set5 using the predictors and the
C    fast transmittance coefficients.


!CALL PROTOCOL:
C    YCALT5( INDCHN, NLAY, NCHN5, CLIST5, COEF5,
C       FIXMUL, CONPD5, FPRED5, WPRED5, OPRED5, TRCPRD, INDCO2, COFCO2,
C       CO2MLT, INDN2O, COFN2O, N2OMLT, TAU, TAUZ )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INT arr   INDCHN  channel indices             none
C    INTEGER   NLAY    number of layers to bottom  none
C    INTEGER   NCHN5   set5 number of channels     none
C    INT arr   CLIST5  set5 channel list           none
C    REAL arr  COEF5   set5 fast trans coefs       various
C    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
C    REAL arr  CONPD5  set5 H2O continuum preds    various
C    REAL arr  FPRED5  set5 fixed gases preds      various
C    REAL arr  WPRED5  set5 water predictors       various
C    REAL arr  OPRED5  set5 ozone predictors       various
C    REAL arr  TRCPRD  trace gas pert predictors   various
C    INT arr   INDCO2  CO2 pert chan indices       none
C    REAL arr  COFCO2  CO2 pert coefs              various
C    REAL arr  CO2MLT  CO2 pert multiplier         none
C    INT arr   INDN2O  N2O pert chan indices       none
C    REAL arr  COFN2O  N2O pert coefs              various
C    REAL arr  N2OMLT  N2O pert multiplier         none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  TAU     effective layer opt depth   none
C    REAL arr  TAUZ    layer-to-space opt depth    none


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    USEFAST


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    August 2000 version of the 100 layer AIRS Fast Transmittance
C    Code by L.L.Strow/S.Hannon.
C
C    The fast trans coefficents and predictors are multiplied
C    together and summed to calculate the effective layer
C    transmittances. Fixed, water, and ozone transmittances are each
C    checked individually to be sure they give 0 < trans < 1.
C
C    ===================================================================
C    Loops downward over all the layers for each of the NCHN5 channels
C    to compute the layer transmittances TAU.
C
C    The water continuum absorption coefficient is:
C       k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
C
C    The layer effective fixed gas absorption coefficient is:
C       k_fixed = the sum i=1 to 11 of { COEF(5+i)*FPRED(i) }
C
C    The layer effective water lines absorption coefficient is:
C       k_water = the sum i=1 to 3 of { COEF(5+11+i)*WPRED(i) }
C
C    The layer effective ozone absorption coefficient is:
C       k_ozone = COEF(5+11+3+1)*OPRED(1)
C
C    where
C      "COEF" are the fast transmittance coefficients COEF5
C      "CONPRD" are the water continuum predictors CONPRD
C      "FPRED" are the fixed gases predictors FPRED5
C      "WPRED" are the water lines predictors WPRED5
C      "OPRED" are the ozone predictors OPRED5
C
C    The total layer effective optical depth TAU is:
C       TAU = [ k_con + k_fixed + k_water + k_ozone ]
C
C    ===================================================================


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C Date        Programmer     Comments
C ----------- -------------- -------------------------------------------
C 26 Jun 1997 Scott Hannon   Created for set5
C 03 Sep 1997 Scott Hannon   Added TAUZ and BLMULT
C 30 Sep 1997 Scott Hannon   Added variable CO2
C 27 Feb 1998 Scott Hannon   Added LTAU
C 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
C 12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
C 03 Jan 2003 Scott Hannon   Add XZ
C 12 Oct 2004 Scott Hannon   Change CO2MLT from scaler to vector
C 28 Jun 2005 Scott Hannon   "trace" version for CO2,N2O
C 28 Mar 2006 Scott Hannon   Change TAU from trans to optical depth
C 22 Dec 2006 Scott Hannon   Change TAUZ from trans to optical depth
C                            and from (1 x n) to (m x n) array;
C                            delete func QIKEXP & arguments BLMULT
C                            & LTAU & XZ.
C 02 Sep 2008 Scott Hannon   Add 5th CO2 predictor

!END====================================================================

C      =================================================================
       SUBROUTINE YCALT5 ( INDCHN, NLAY, NCHN5, CLIST5,
     $    COEF5, FIXMUL, CONPD5, FPRED5, WPRED5, OPRED5, TRCPRD, INDCO2,
     $    COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, TAU, TAUZ, IY, 
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
C      =================================================================

C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
       include 'incFTC.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input
       INTEGER IY
       INTEGER INDCHN(MXCHAN)
       INTEGER   NLAY
       INTEGER  NCHN5
       INTEGER CLIST5(MXCHN5)
       REAL  COEF5(N5COEF,MAXLAY,MXCHN5)
       REAL FIXMUL(MAXLAY)
       REAL CONPD5( N5CON,MAXLAY)
       REAL FPRED5( N5FIX,MAXLAY)
       REAL WPRED5( N5H2O,MAXLAY)
       REAL OPRED5(  N5O3,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2(  NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT(MAXLAY)
       INTEGER INDN2O(MXCHAN)
       REAL COFN2O(  NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
C
C      Output
       REAL    TAU(MAXLAY,MXCHAN)
       REAL   TAUZ(MAXLAY,MXCHAN)
       REAL DTAU_DTZ(MAXLAY,MXCHAN),DTAU_DG1(MAXLAY,MXCHAN),DTAU_DG2(MAXLAY,MXCHAN),DTAU_DG3(MAXLAY,MXCHAN)
       REAL DTAU_DG4(MAXLAY,MXCHAN),DTAU_DG5(MAXLAY,MXCHAN),DTAU_DG6(MAXLAY,MXCHAN),DTAU_DG9(MAXLAY,MXCHAN)
       REAL DTAU_DG12(MAXLAY,MXCHAN)

c input
       LOGICAL DOJAC
       INTEGER NWANTJ          ! number of wanted jacs (default 0=none)
       INTEGER  LISTJ(MAXPRO)  ! list of wanted channels
       REAL CO2JACMLT(MAXLAY)
       REAL SO2JACMLT(MAXLAY)
       REAL HNOJACMLT(MAXLAY)
       REAL N2OJACMLT(MAXLAY)
       REAL NH3JACMLT(MAXLAY)
       REAL HDOJACMLT(MAXLAY)
       !!! first index is the d/dT   second deriv is the d/dQ
       REAL CONJACPRD(MAXJAC, N1CON,MAXLAY)
       REAL FJACPRED1(MAXJAC, N1FIX,MAXLAY)
       REAL FJACPRED2(MAXJAC, N2FIX,MAXLAY)
       REAL FJACPRED3(MAXJAC, N3FIX,MAXLAY)
       REAL FJACPRED4(MAXJAC, N4FIX,MAXLAY)
       REAL FJACPRED5(MAXJAC, N5FIX,MAXLAY)
       REAL FJACPRED6(MAXJAC, N6FIX,MAXLAY)
       REAL FJACPRED7(MAXJAC, N7FIX,MAXLAY)
       REAL WJACPRED1(MAXJAC, N1H2O,MAXLAY)
       REAL WJACPRED2(MAXJAC, N2H2O,MAXLAY)
       REAL WJACPRED3(MAXJAC, N3H2O,MAXLAY)
       REAL WJACPRED4(MAXJAC, N4H2O,MAXLAY)
       REAL WJACPRED5(MAXJAC, N5H2O,MAXLAY)
       REAL WJACPRED6(MAXJAC, N6H2O,MAXLAY)
       REAL WJACPRED7(MAXJAC, N7H2O,MAXLAY)
       REAL OJACPRED1(MAXJAC,  N1O3,MAXLAY)
       REAL OJACPRED2(MAXJAC,  N2O3,MAXLAY)
       REAL OJACPRED4(MAXJAC,  N4O3,MAXLAY)
       REAL OJACPRED5(MAXJAC,  N5O3,MAXLAY)
       REAL OJACPRED6(MAXJAC,  N6O3,MAXLAY)
       REAL OJACPRED7(MAXJAC,  N7O3,MAXLAY)
       REAL  DJACPRED(MAXJAC,  NHDO,MAXLAY)
       REAL MJACPRED3(MAXJAC, N3CH4,MAXLAY)
       REAL CJACPRED4(MAXJAC,  N4CO,MAXLAY)
       REAL TRCJACPRD(MAXJAC,NTRACE,MAXLAY)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      IWHICHJAC,ITRYJAC,INTERSECT
       INTEGER      I
       INTEGER   ICO2
       INTEGER   ILAY
       INTEGER   IN2O
       INTEGER      J
       REAL     DK
       REAL  DKCO2, QDKCO2, RAQDKCO2(MAXLAY)
       REAL  DKN2O, QDKN2O, RAQDKN2O(MAXLAY)
c       REAL  DKCO2
c       REAL  DKN2O
       REAL   KCON
       REAL   KFIX
       REAL KLAYER
       REAL   KOZO
       REAL   KWAT
       REAL     KZ
       LOGICAL   LCO2
       LOGICAL   LN2O


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
C      ---------------------------
C      Loop on channel (frequency)
C      ---------------------------
C       DO I=1,NCHN5
C        DO I=IY,IY
         I = IY
C
C         Index for TAU
          J=INDCHN( CLIST5(I) )
C
C         Determine whether or not to do variable CO2
          ICO2=INDCO2( CLIST5(I) )
          IF (ICO2 .GT. 0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable CO2
          IN2O=INDN2O( CLIST5(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
C
C         Initialize the layer-to-space optical depth
          KZ=0.0E+0
C
C         ------------------------------
C         Loop on layers (top to ground)
C         ------------------------------
          DO ILAY=1,NLAY
C
C            ---------------------------
C            Compute the water continuum
C            ---------------------------
             KCON=( COEF5(1,ILAY,I)*CONPD5(1,ILAY) ) +
     $            ( COEF5(2,ILAY,I)*CONPD5(2,ILAY) ) +
     $            ( COEF5(3,ILAY,I)*CONPD5(3,ILAY) ) +
     $            ( COEF5(4,ILAY,I)*CONPD5(4,ILAY) ) +
     $            ( COEF5(5,ILAY,I)*CONPD5(5,ILAY) ) +
     $            ( COEF5(6,ILAY,I)*CONPD5(6,ILAY) ) +
     $            ( COEF5(7,ILAY,I)*CONPD5(7,ILAY) )
C
             IF (KCON .LT. 0.0E+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 1.0E+1) THEN
                KCON=1.0E+1
             ENDIF
C
C            -----------------------------
C            Calc the fixed gases abs coef
C            -----------------------------
             KFIX=( COEF5( 8,ILAY,I)*FPRED5( 1,ILAY) ) +
     $            ( COEF5( 9,ILAY,I)*FPRED5( 2,ILAY) ) +
     $            ( COEF5(10,ILAY,I)*FPRED5( 3,ILAY) ) +
     $            ( COEF5(11,ILAY,I)*FPRED5( 4,ILAY) ) +
     $            ( COEF5(12,ILAY,I)*FPRED5( 5,ILAY) ) +
     $            ( COEF5(13,ILAY,I)*FPRED5( 6,ILAY) ) +
     $            ( COEF5(14,ILAY,I)*FPRED5( 7,ILAY) ) +
     $            ( COEF5(15,ILAY,I)*FPRED5( 8,ILAY) ) +
     $            ( COEF5(16,ILAY,I)*FPRED5( 9,ILAY) ) +
     $            ( COEF5(17,ILAY,I)*FPRED5(10,ILAY) ) +
     $            ( COEF5(18,ILAY,I)*FPRED5(11,ILAY) )
C
             KFIX=KFIX*FIXMUL(ILAY)
C
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
C
C
C            --------------------------
C            Compute the water abs coef
C            --------------------------
             KWAT=( COEF5(19,ILAY,I)*WPRED5( 1,ILAY) ) +
     $            ( COEF5(20,ILAY,I)*WPRED5( 2,ILAY) ) +
     $            ( COEF5(21,ILAY,I)*WPRED5( 3,ILAY) )
C
             IF (KWAT .LT. 0.0E+0) THEN
                KWAT=0.0E+0
             ELSEIF( KWAT .GT. 1.0E+1) THEN
                KWAT=1.0E+1
             ENDIF
C
C
C            --------------------------
C            Compute the ozone abs coef
C            --------------------------
             KOZO=( COEF5(22,ILAY,I)*OPRED5(1,ILAY) )
C
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+1) THEN
                KOZO=1.0E+1
             ENDIF
C
C
C            ----------------------------------
C            Calc the total layer transmittance
C            ----------------------------------
c
ccccc
c This block is usually commented out and is only uncommented for
c testing purposes.
c
c           kcon=0.0E+0
c           kfix=0.0E+0
c           kwat=0.0E+0
c           kozo=0.0E+0
ccccc
C
C            ----------------------------
C            Calc change in total optical
C            depth due to variable CO2
C            ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0.0) THEN
                RAQDKCO2(ILAY)=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) +
     $                ( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) +
     $                ( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) +
     $                ( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) +
     $                ( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
                DKCO2=RAQDKCO2(ILAY)*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
C
C            ----------------------------
C            Calc change in total optical
C            depth due to variable N2O
C            ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0.0) THEN
                RAQDKN2O(ILAY)=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) +
     $                ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) +
     $                ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) +
     $                ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) +
     $                ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) +
     $                ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) +
     $                ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=RAQDKN2O(ILAY)*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
C
ccc
c this block for testing
c      DKCO2=0.0
c      DKN2O=0.0
ccc
C            Limit -DK so it can never totally totally cancel KFIX
             DK = DKCO2 + DKN2O
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF

C            Calc total layer optical depth
             KLAYER = KCON + KFIX + KWAT + KOZO + DK
             TAU(ILAY,J)=KLAYER
C
C            Calc layer-to-space optical depth
             KZ=KZ + KLAYER
             TAUZ(ILAY,J)=KZ
C
          ENDDO
C         End loop on levels
C
C       ENDDO  ! DO I = IY,IY
C      End loops on channel number (frequency)
C

c************************************************************************

          IF (DOJAC) THEN
            DO IWHICHJAC = 1,9
              IF (IWHICHJAC .EQ. 1) THEN 
                ITRYJAC = 100  !! TZ
              ELSEIF (IWHICHJAC .EQ. 2) THEN 
                ITRYJAC = 1    !! GID 1
              ELSEIF (IWHICHJAC .EQ. 3) THEN 
                ITRYJAC = 3    !! GID 3
              ELSEIF (IWHICHJAC .EQ. 4) THEN 
                ITRYJAC = 2    !! GID 2
              ELSEIF (IWHICHJAC .EQ. 5) THEN 
                ITRYJAC = 4    !! GID 4
              ELSEIF (IWHICHJAC .EQ. 6) THEN 
                ITRYJAC = 5    !! GID 5
              ELSEIF (IWHICHJAC .EQ. 7) THEN 
                ITRYJAC = 6    !! GID 6
              ELSEIF (IWHICHJAC .EQ. 8) THEN 
                ITRYJAC = 9    !! GID 9
              ELSEIF (IWHICHJAC .EQ. 9) THEN 
                ITRYJAC = 12   !! GID 12
              END IF
              IF (INTERSECT(ITRYJAC,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
cbaba
C         Initialize the layer-to-space optical depth
          KZ=0.0E+0
C
C         ------------------------------
C         Loop on layers (top to ground)
C         ------------------------------
          DO ILAY=1,NLAY
C
C            ---------------------------
C            Compute the water continuum
C            ---------------------------
             KCON=( COEF5(1,ILAY,I)*CONJACPRD(IWHICHJAC,1,ILAY) ) +
     $            ( COEF5(2,ILAY,I)*CONJACPRD(IWHICHJAC,2,ILAY) ) +
     $            ( COEF5(3,ILAY,I)*CONJACPRD(IWHICHJAC,3,ILAY) ) +
     $            ( COEF5(4,ILAY,I)*CONJACPRD(IWHICHJAC,4,ILAY) ) +
     $            ( COEF5(5,ILAY,I)*CONJACPRD(IWHICHJAC,5,ILAY) ) +
     $            ( COEF5(6,ILAY,I)*CONJACPRD(IWHICHJAC,6,ILAY) ) +
     $            ( COEF5(7,ILAY,I)*CONJACPRD(IWHICHJAC,7,ILAY) )
C
c             IF (KCON .LT. 0.0E+0) THEN
c                KCON=0.0E+0
c             ELSEIF (KCON .GT. 1.0E+1) THEN
c                KCON=1.0E+1
c             ENDIF
C
C            -----------------------------
C            Calc the fixed gases abs coef
C            -----------------------------
             KFIX=( COEF5( 8,ILAY,I)*FJACPRED5(IWHICHJAC, 1,ILAY) ) +
     $            ( COEF5( 9,ILAY,I)*FJACPRED5(IWHICHJAC, 2,ILAY) ) +
     $            ( COEF5(10,ILAY,I)*FJACPRED5(IWHICHJAC, 3,ILAY) ) +
     $            ( COEF5(11,ILAY,I)*FJACPRED5(IWHICHJAC, 4,ILAY) ) +
     $            ( COEF5(12,ILAY,I)*FJACPRED5(IWHICHJAC, 5,ILAY) ) +
     $            ( COEF5(13,ILAY,I)*FJACPRED5(IWHICHJAC, 6,ILAY) ) +
     $            ( COEF5(14,ILAY,I)*FJACPRED5(IWHICHJAC, 7,ILAY) ) +
     $            ( COEF5(15,ILAY,I)*FJACPRED5(IWHICHJAC, 8,ILAY) ) +
     $            ( COEF5(16,ILAY,I)*FJACPRED5(IWHICHJAC, 9,ILAY) ) +
     $            ( COEF5(17,ILAY,I)*FJACPRED5(IWHICHJAC,10,ILAY) ) +
     $            ( COEF5(18,ILAY,I)*FJACPRED5(IWHICHJAC,11,ILAY) )
C
             KFIX=KFIX*FIXMUL(ILAY)
C
c             IF (KFIX .LT. 0.0E+0) THEN
c                KFIX=0.0E+0
c             ELSEIF (KFIX .GT. 1.0E+1) THEN
c                KFIX=1.0E+1
c             ENDIF
C
C
C            --------------------------
C            Compute the water abs coef
C            --------------------------
             KWAT=( COEF5(19,ILAY,I)*WJACPRED5(IWHICHJAC, 1,ILAY) ) +
     $            ( COEF5(20,ILAY,I)*WJACPRED5(IWHICHJAC, 2,ILAY) ) +
     $            ( COEF5(21,ILAY,I)*WJACPRED5(IWHICHJAC, 3,ILAY) )
C
c             IF (KWAT .LT. 0.0E+0) THEN
c                KWAT=0.0E+0
c             ELSEIF( KWAT .GT. 1.0E+1) THEN
c                KWAT=1.0E+1
c             ENDIF
C
C
C            --------------------------
C            Compute the ozone abs coef
C            --------------------------
             KOZO=( COEF5(22,ILAY,I)*OJACPRED5(IWHICHJAC, 1,ILAY) )
C
c             IF (KOZO .LT. 0.0E+0) THEN
c                KOZO=0.0E+0
c             ELSEIF (KOZO .GT. 1.0E+1) THEN
c                KOZO=1.0E+1
c             ENDIF
C
C
C            ----------------------------------
C            Calc the total layer transmittance
C            ----------------------------------
c
ccccc
c This block is usually commented out and is only uncommented for
c testing purposes.
c
c           kcon=0.0E+0
c           kfix=0.0E+0
c           kwat=0.0E+0
c           kozo=0.0E+0
ccccc
C
C            ----------------------------
C            Calc change in total optical
C            depth due to variable CO2
C            ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0.0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,1,ILAY) ) +
     $                ( COFCO2(2,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,2,ILAY) ) +
     $                ( COFCO2(3,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,3,ILAY) ) +
     $                ( COFCO2(4,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,4,ILAY) ) +
     $                ( COFCO2(5,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF

             IF (LCO2) THEN
               QDKCO2 = RAQDKCO2(ILAY)*CO2JACMLT(ILAY)
             END IF
C
C            ----------------------------
C            Calc change in total optical
C            depth due to variable N2O
C            ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0.0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCJACPRD(IWHICHJAC,1,ILAY) ) +
     $                ( COFN2O(2,ILAY,IN2O)*TRCJACPRD(IWHICHJAC,2,ILAY) ) +
     $                ( COFN2O(3,ILAY,IN2O)*TRCJACPRD(IWHICHJAC,3,ILAY) ) +
     $                ( COFN2O(4,ILAY,IN2O)*TRCJACPRD(IWHICHJAC,4,ILAY) ) +
     $                ( COFN2O(5,ILAY,IN2O)*TRCJACPRD(IWHICHJAC,5,ILAY) ) +
     $                ( COFN2O(6,ILAY,IN2O)*TRCJACPRD(IWHICHJAC,6,ILAY) ) +
     $                ( COFN2O(7,ILAY,IN2O)*TRCJACPRD(IWHICHJAC,7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF

             IF (LN2O) THEN
               QDKN2O = RAQDKN2O(ILAY)*N2OJACMLT(ILAY)
             END IF
C
ccc
c this block for testing
c      DKCO2=0.0
c      DKN2O=0.0
ccc
C            Limit -DK so it can never totally totally cancel KFIX
c             DK = DKCO2 + DKN2O
c             IF (-DK .GE. KFIX) THEN
c                DK = -0.999*KFIX
c             ENDIF
              DK = 0

C            Calc total layer optical depth
             KLAYER = KCON + KFIX + KWAT + KOZO + DK

                   IF (IWHICHJAC .EQ. 1) THEN 
                     DTAU_DTZ(ILAY,J)=KLAYER
                   ELSEIF (IWHICHJAC .EQ. 2) THEN 
                     DTAU_DG1(ILAY,J)=KLAYER
                   ELSEIF (IWHICHJAC .EQ. 3) THEN 
                     DTAU_DG3(ILAY,J)=KLAYER
                   ELSEIF (IWHICHJAC .EQ. 4) THEN 
                     DTAU_DG2(ILAY,J)=QDKCO2
                   ELSEIF (IWHICHJAC .EQ. 5) THEN 
                     DTAU_DG4(ILAY,J)=QDKN2O
                   ELSEIF (IWHICHJAC .EQ. 6) THEN 
                     DTAU_DG5(ILAY,J)=0
                   ELSEIF (IWHICHJAC .EQ. 7) THEN 
                     DTAU_DG6(ILAY,J)=0
                   ELSEIF (IWHICHJAC .EQ. 8) THEN 
                     DTAU_DG9(ILAY,J)=0
                   ELSEIF (IWHICHJAC .EQ. 9) THEN 
                     DTAU_DG12(ILAY,J)=0
                   END IF

c             TAU(ILAY,J)=KLAYER
C
C            Calc layer-to-space optical depth
c             KZ=KZ + KLAYER
c             TAUZ(ILAY,J)=KZ
C
          ENDDO
C         End loop on levels
cbaba
              END IF   !!!! IF (INTERSECT(ITRYJAC,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) : is this a jac to work on?????
            END DO     !!! DO IWHICHJAC = 1,3
          END IF  !! DOJAC

c************************************************************************

       RETURN
       END
