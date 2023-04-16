C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALT1 (for set1 = FWO)
C
!F77====================================================================


!ROUTINE NAME:
C    CALT1


!ABSTRACT:
C    Calculate the transmittance for set1 using the predictors
C    and the fast transmittance coefficients.


!CALL PROTOCOL:
C    CALT1 ( INDCHN, NLAY, BLMULT, NCHN1, CLIST1, COEF1,
C      FIXMUL, CONPD1, FPRED1, WPRED1, OPRED1, CO2PRD, INDCO2,
C      COFCO2, CO2MLT, INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX,
C      LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      BLMULT  bottom layer opt depth mult none
C    INT arr   CLIST1  set channel list            none
C    REAL arr  COEF1   set1 fast trans coefs       various
C    REAL arr  CONPD1  set1 H2O continuum preds    various
C    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
C    REAL arr  FPRED1  set1 fixed gases preds      various
C    INT arr   INDCHN  channel indices             none
C    INTEGER   NLAY    Number of layers to bottom  none
C    INTEGER   NCHN1   set1 number of channels     none
C    REAL arr  OPRED1  set1 ozone predictors       various
C    REAL arr  WPRED1  set1 water predictors       various
C    REAL arr  CO2PRD  CO2 pert predictors         various
C    INT arr   INDCO2  CO2 pert chan indices       none
C    REAL arr  COFCO2  CO2 pert coefs              various
C    REAL      CO2MLT  CO2 pert multiplier         none


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
C    The FTC coefficents and profile FTC variables are multiplied
C    together and summed to calculate the effective layer
C    transmittances. Fixed, water, and ozone transmittances are each
C    checked individually to be sure they give 0 < trans < 1.
C
C    ===================================================================
C    The routine loops over the selected channels of set1.  For each
C    channel, it first decides if needs to do a calculation for
C    variable CO2, and also if it needs to do an OPTRAN water calc (if
C    so, it does so immediately).  The program then loops downward over
C    all the layers and computes the layer transmittances TAU.
C
C    The water continuum absorption coefficient is:
C       k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
C
C    The layer effective fixed gas absorption coefficient is:
C       k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
C
C    The layer effective water lines absorption coefficient is:
C       k_water = the sum i=1 to 11 of { COEF(5+8+i)*WPRED(i) }
C
C    The layer effective ozone absorption coefficient is:
C       k_ozone = the sum i=1 to 5 of { COEF(5+8+11+i)*OPRED(i) }
C
C    where
C      "COEF" are the fast transmittance coefficients COEF1
C      "CONPRD" are the water continuum predictors CONPRD
C      "FPRED" are the fixed gases predictors FPRED1
C      "WPRED" are the water lines predictors WPRED1
C      "OPRED" are the ozone predictors OPRED1
C
C    The total layer effective transmittance TAU is:
C       TAU = exp( -[ k_con + k_fixed + k_water + k_ozone ])
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
C    Dec  1 1994 Scott Hannon   Created
C     3 Feb 1997 Scott Hannon   Re-wrote (from CALTAU) for FWO
C     3 Sep 1997 Scott Hannon   Added TAUZ and BLMULT
C    30 Sep 1997 Scott Hannon   Added variable CO2
C    27 Feb 1998 Scott Hannon   Added OPTRAN H2O
C     6 Mar 1998 Scott Hannon   Deleted water preds 12 & 13 and shifted
C                               ozone coef indices to 24-28 (was 26-30)
C     4 May 1998 Scott Hannon   Fix error: INDH2O(MXCHAN) not (MXCHNW)
C    26 Aug 1998 Scott Hannon   Add NLAY to call to CALOKW
C    11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
C    12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
ca

!END====================================================================

C      =================================================================
       SUBROUTINE CALT1 ( INDCHN, NLAY, BLMULT, NCHN1, CLIST1, COEF1,
     $    FIXMUL, CONPD1, FPRED1, WPRED1, OPRED1,
     $    CO2PRD, INDCO2, COFCO2, CO2MLT,
     $    INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX,
     $    LOPLOW, LOPUSE, WAOP, DAOP, WAANG,
     $    TAU, TAUZ)
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
       INTEGER INDCHN(MXCHAN)
       INTEGER   NLAY
       REAL BLMULT
       INTEGER  NCHN1
       INTEGER CLIST1(MXCHN1)
       REAL  COEF1(N1COEF,MAXLAY,MXCHN1)
       REAL FIXMUL(MAXLAY)
       REAL CONPD1( N1CON,MAXLAY)
       REAL FPRED1( N1FIX,MAXLAY)
       REAL WPRED1( N1H2O,MAXLAY)
       REAL OPRED1(  N1O3,MAXLAY)
       REAL CO2PRD(  NCO2,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2(  NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT
       INTEGER INDH2O(MXCHAN)
       REAL H2OPRD(  NH2O,MXOWLY)
       REAL COFH2O(  NH2O,MXOWLY,MXCHNW)
       INTEGER LOPMIN
       INTEGER LOPMAX
       INTEGER LOPLOW(MAXLAY)
       LOGICAL LOPUSE(MXOWLY)
       REAL   WAOP(MXOWLY)
       REAL   DAOP(MAXLAY)
       REAL  WAANG(MAXLAY)
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
       REAL     KZ
       REAL   KZFW
       LOGICAL   LCO2
       LOGICAL   LH2O
C
C      for function QIKEXP
       REAL QIKEXP
C
C      for CALOKW
       INTEGER   IH2O
       REAL     KW(MAXLAY)


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
       DO I=1,NCHN1
C
C         Array index of channel in TAU
          J=INDCHN( CLIST1(I) )
C
C         Determine whether or not to do variable CO2 calc
          ICO2=INDCO2( CLIST1(I) )
          IF (ICO2 .GT. 0 .AND. CO2MLT .NE. 0.0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
C
C         -------------------------
C         Do OPTRAN water if needed
C         -------------------------
          IH2O=INDH2O( CLIST1(I) )
          IF (IH2O .GT. 0) THEN
             LH2O=.FALSE.
C            Calc OPTRAN water
C
             CALL CALOKW( NLAY, IH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE,
     $          H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW )
C
          ELSE
             LH2O=.TRUE.
          ENDIF
C
C         Initialize the layer-to-space optical depth
          KZ=0.0E+0
          KZFW=0.0E+0
C
C         ------------------------------
C         Loop on layers (top to bottom)
C         ------------------------------
          DO ILAY=1,NLAY
C
C            ---------------------------
C            Compute the water continuum
C            ---------------------------
             KCON=( COEF1(1,ILAY,I)*CONPD1(1,ILAY) ) +
     $            ( COEF1(2,ILAY,I)*CONPD1(2,ILAY) ) +
     $            ( COEF1(3,ILAY,I)*CONPD1(3,ILAY) ) +
     $            ( COEF1(4,ILAY,I)*CONPD1(4,ILAY) ) +
     $            ( COEF1(5,ILAY,I)*CONPD1(5,ILAY) ) +
     $            ( COEF1(6,ILAY,I)*CONPD1(6,ILAY) ) +
     $            ( COEF1(7,ILAY,I)*CONPD1(7,ILAY) )
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
             KFIX=( COEF1( 8,ILAY,I)*FPRED1(1,ILAY) ) +
     $            ( COEF1( 9,ILAY,I)*FPRED1(2,ILAY) ) +
     $            ( COEF1(10,ILAY,I)*FPRED1(3,ILAY) ) +
     $            ( COEF1(11,ILAY,I)*FPRED1(4,ILAY) ) +
     $            ( COEF1(12,ILAY,I)*FPRED1(5,ILAY) ) +
     $            ( COEF1(13,ILAY,I)*FPRED1(6,ILAY) ) +
     $            ( COEF1(14,ILAY,I)*FPRED1(7,ILAY) ) +
     $            ( COEF1(15,ILAY,I)*FPRED1(8,ILAY) )
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
             IF (LH2O) THEN
C               Not an OPTRAN water channel
                KW(ILAY)=
     $               ( COEF1(16,ILAY,I)*WPRED1( 1,ILAY) ) +
     $               ( COEF1(17,ILAY,I)*WPRED1( 2,ILAY) ) +
     $               ( COEF1(18,ILAY,I)*WPRED1( 3,ILAY) ) +
     $               ( COEF1(19,ILAY,I)*WPRED1( 4,ILAY) ) +
     $               ( COEF1(20,ILAY,I)*WPRED1( 5,ILAY) ) +
     $               ( COEF1(21,ILAY,I)*WPRED1( 6,ILAY) ) +
     $               ( COEF1(22,ILAY,I)*WPRED1( 7,ILAY) ) +
     $               ( COEF1(23,ILAY,I)*WPRED1( 8,ILAY) ) +
     $               ( COEF1(24,ILAY,I)*WPRED1( 9,ILAY) ) +
     $               ( COEF1(25,ILAY,I)*WPRED1(10,ILAY) ) +
     $               ( COEF1(26,ILAY,I)*WPRED1(11,ILAY) )
C
                IF (KW(ILAY) .LT. 0.0E+0) KW(ILAY)=0.0E+0
             ENDIF
C
C            --------------------------
C            Compute the ozone abs coef
C            --------------------------
             KOZO=( COEF1(27,ILAY,I)*OPRED1(1,ILAY) ) +
     $            ( COEF1(28,ILAY,I)*OPRED1(2,ILAY) ) +
     $            ( COEF1(29,ILAY,I)*OPRED1(3,ILAY) ) +
     $            ( COEF1(30,ILAY,I)*OPRED1(4,ILAY) ) +
     $            ( COEF1(31,ILAY,I)*OPRED1(5,ILAY) )
C
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+1) THEN
                KOZO=1.0E+1
             ENDIF
C
C            Update KZFW
             KZFW=KZFW + KFIX + KW(ILAY)
C
C            ----------------------------------
C            Calc the total layer transmittance
C            ----------------------------------
c
ccccc
c This block is usually commented out and is only uncommented for
c testing purposes.
cc
c           kcon=0.0E+0
c           kfix=0.0E+0
c           kw(ilay)=0.0E+0
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
             KLAYER=KCON + KFIX + KW(ILAY) + KOZO + DKCO2
C
C            Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
C
C            Calc layer-to-space optical depth
             KZ=KZ + KLAYER
C
C            Calc effective layer transmittance
             TAU(ILAY,J)=QIKEXP(-KLAYER)
C
          ENDDO
C         End loop on levels
C
C         Calc surface-to-space transmittance
          TAUZ(J)=QIKEXP(-KZ)
C
       ENDDO
C      End loops on channel number (frequency)
C
       RETURN
       END
