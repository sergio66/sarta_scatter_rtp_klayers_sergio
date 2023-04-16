C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALT6 (for set6 = FWO sun mfmw)
C
!F77====================================================================


!ROUTINE NAME:
C    CALT6


!ABSTRACT:
C    Calculate the transmittance for set6 using the predictors and the
C    fast transmittance coefficients.


!CALL PROTOCOL:
C    CALT6( LTAU, INDCHN, NLAY, BMULT, NCHN6, CLIST6, COEF6,
C       FIXMUL, CONPD6, FPRED6, WPRED6, OPRED6, CO2PRD, INDCO2, COFCO2,
C       CO2MLT, XZ, TAU, TAUZ )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    LOGICAL   LTAU    Calc all layer trans?       none
C    INT arr   INDCHN  channel indices             none
C    INTEGER   NLAY    number of layers to bottom  none
C    REAL      BLMULT  bottom layer opt depth mult none
C    INTEGER   NCHN6   set6 number of channels     none
C    INT arr   CLIST6  set6 channel list           none
C    REAL arr  COEF6   set6 fast trans coefs       various
C    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
C    REAL arr  CONPD6  set6 H2O continuum preds    various
C    REAL arr  FPRED6  set6 fixed gases preds      various
C    REAL arr  WPRED6  set6 water predictors       various
C    REAL arr  OPRED6  set6 ozone predictors       various
C    REAL arr  CO2PRD  CO2 pert predictors         various
C    INT arr   INDCO2  CO2 pert chan indices       none
C    REAL arr  COFCO2  CO2 pert coefs              various
C    REAL      CO2MLT  CO2 pert multiplier         none
C    REAL      XZ      optical depth mult for TAUZ none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  TAU     effective layer trans       none
C    REAL arr  TAUZ    layer-to-space trans        none


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
C    Loops downward over all the layers for each of the NCHN6 channels
C    to compute the layer transmittances TAU.
C
C    The water continuum absorption coefficient is:
C       k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
C
C    The layer effective fixed gas absorption coefficient is:
C       k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
C
C    The layer effective water lines absorption coefficient is:
C       k_water = the sum i=1 to 7 of { COEF(5+8+i)*WPRED(i) }
C
C    The layer effective ozone absorption coefficient is:
C       k_ozone = COEF(5+8+7+1)*OPRED(1)
C
C    where
C      "COEF" are the fast transmittance coefficients COEF5
C      "CONPRD" are the water continuum predictors CONPRD
C      "FPRED" are the fixed gases predictors FPRED5
C      "WPRED" are the water lines predictors WPRED5
C      "OPRED" are the ozone predictors OPRED5
C
C    The total layer effective transmittance TAU is:
C       TAU = exp( -[ k_con + k_fixed + k_water + k_ozone ])
C    TAU is only calc'ed if LTAU is TRUE.
C
C    To help speed up the exponential calculations, we use our own
C    "EXP" replacement function called QIKEXP which uses just the
C    first few series expansion terms for exp(x) if x is suitably small.
C    ===================================================================


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date        Programmer     Comments
C    ----------- -------------- ----------------------------------------
C     3 Jul 1997 Scott Hannon   Created for set6
C     3 Sep 1997 Scott Hannon   Added TAUZ and BLMULT
C    30 Sep 1997 Scott Hannon   Added variable CO2
C    27 Feb 1998 Scott Hannon   Added LTAU
C    11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
C    12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
C     3 Jan 2003 Scott Hannon   Add XZ
C     6 Feb 2003 Scott Hannon   Bug fix - ozone use coef 23, not 21


!END====================================================================

C      =================================================================
       SUBROUTINE CALT6 ( LTAU, INDCHN, NLAY, BLMULT, NCHN6, CLIST6,
     $    COEF6, FIXMUL, CONPD6, FPRED6, WPRED6, OPRED6, CO2PRD, INDCO2,
     $    COFCO2, CO2MLT, XZ, TAU, TAUZ )
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
C      QIKEXP


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input
       LOGICAL   LTAU
       INTEGER INDCHN(MXCHAN)
       INTEGER   NLAY
       REAL BLMULT
       INTEGER  NCHN6
       INTEGER CLIST6(MXCHN6)
       REAL  COEF6(N6COEF,MAXLAY,MXCHN6)
       REAL FIXMUL(MAXLAY)
       REAL CONPD6( N6CON,MAXLAY)
       REAL FPRED6( N6FIX,MAXLAY)
       REAL WPRED6( N6H2O,MAXLAY)
       REAL OPRED6(  N6O3,MAXLAY)
       REAL CO2PRD(  NCO2,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2(  NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT
       REAL     XZ
C
C      Output
       REAL    TAU(MAXLAY,MXCHAN)
       REAL   TAUZ(MXCHAN)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      I
       INTEGER   ICO2
       INTEGER   ILAY
       INTEGER      J
       REAL  DKCO2
       REAL   KCON
       REAL   KFIX
       REAL KLAYER
       REAL   KOZO
       REAL   KWAT
       REAL     KZ
       LOGICAL   LCO2
C
C      for function QIKEXP
       REAL QIKEXP


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
       DO I=1,NCHN6
C
C         Index for TAU
          J=INDCHN( CLIST6(I) )
C
C         Determine whether or not to do variable CO2
          ICO2=INDCO2( CLIST6(I) )
          IF (ICO2 .GT. 0 .AND. CO2MLT .NE. 0.0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
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
             KCON=( COEF6(1,ILAY,I)*CONPD6(1,ILAY) ) +
     $            ( COEF6(2,ILAY,I)*CONPD6(2,ILAY) ) +
     $            ( COEF6(3,ILAY,I)*CONPD6(3,ILAY) ) +
     $            ( COEF6(4,ILAY,I)*CONPD6(4,ILAY) ) +
     $            ( COEF6(5,ILAY,I)*CONPD6(5,ILAY) ) +
     $            ( COEF6(6,ILAY,I)*CONPD6(6,ILAY) ) +
     $            ( COEF6(7,ILAY,I)*CONPD6(7,ILAY) )
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
             KFIX=( COEF6( 8,ILAY,I)*FPRED6( 1,ILAY) ) +
     $            ( COEF6( 9,ILAY,I)*FPRED6( 2,ILAY) ) +
     $            ( COEF6(10,ILAY,I)*FPRED6( 3,ILAY) ) +
     $            ( COEF6(11,ILAY,I)*FPRED6( 4,ILAY) ) +
     $            ( COEF6(12,ILAY,I)*FPRED6( 5,ILAY) ) +
     $            ( COEF6(13,ILAY,I)*FPRED6( 6,ILAY) ) +
     $            ( COEF6(14,ILAY,I)*FPRED6( 7,ILAY) ) +
     $            ( COEF6(15,ILAY,I)*FPRED6( 8,ILAY) )
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
             KWAT=( COEF6(16,ILAY,I)*WPRED6( 1,ILAY) ) +
     $            ( COEF6(17,ILAY,I)*WPRED6( 2,ILAY) ) +
     $            ( COEF6(18,ILAY,I)*WPRED6( 3,ILAY) ) +
     $            ( COEF6(19,ILAY,I)*WPRED6( 4,ILAY) ) +
     $            ( COEF6(20,ILAY,I)*WPRED6( 5,ILAY) ) +
     $            ( COEF6(21,ILAY,I)*WPRED6( 6,ILAY) ) +
     $            ( COEF6(22,ILAY,I)*WPRED6( 7,ILAY) )
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
             KOZO=( COEF6(23,ILAY,I)*OPRED6(1,ILAY) )
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
             IF (LCO2) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*CO2PRD(1,ILAY) ) +
     $                ( COFCO2(2,ILAY,ICO2)*CO2PRD(2,ILAY) ) +
     $                ( COFCO2(3,ILAY,ICO2)*CO2PRD(3,ILAY) ) +
     $                ( COFCO2(4,ILAY,ICO2)*CO2PRD(4,ILAY) )
                DKCO2=DKCO2*CO2MLT*FIXMUL(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
C
C            ------------------------------------------
C            Calc total optical depth and transmittance
C            ------------------------------------------
C            Calc total layer optical depth
             KLAYER=KCON + KFIX + KWAT + KOZO + DKCO2
C
C            Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
C
C            Calc layer-to-space optical depth
             KZ=KZ + KLAYER
C
C            Calc effective layer transmittance
             IF (LTAU) TAU(ILAY,J)=QIKEXP(-KLAYER)
C
          ENDDO
C         End loop on levels
C
C         Convert KZ to TAUZ
          TAUZ(J)=QIKEXP(-KZ*XZ)
C
       ENDDO
C      End loops on channel number (frequency)
C
       RETURN
       END
