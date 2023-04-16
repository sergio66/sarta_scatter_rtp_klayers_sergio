C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALT3 (for set3 = FMW)
C
!F77====================================================================


!ROUTINE NAME:
C    CALT3


!ABSTRACT:
C    Calculate the transmittance for set3 using the prdictor and the
C    fast transmittance coefficients.


!CALL PROTOCOL:
C    CALT3 ( INDCHN, NLAY, BLMULT, NCHN3, CLIST3, COEF3,
C       FIXMUL, CONPD3, FPRED3, MPRED3, WPRED3,
C       INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX,
C       LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INT arr   INDCHN  channel indices             none
C    INTEGER   NLAY    number of layers to bottom  none
C    REAL      BLMULT  bottom layer opt depth mult none
C    INTEGER   NCHN3   set3 number of channels     none
C    INT arr   CLIST3  set3 channel list           none
C    REAL arr  COEF3   set3 fast trans coefs       various
C    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
C    REAL arr  CONPD3  set3 H2O continuum preds    various
C    REAL arr  FPRED3  set3 fixed gases preds      various
C    REAL arr  MPRED3  set3 methane predictors     various
C    REAL arr  WPRED3  set3 water predictors       various


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
C    transmittances. Fixed, methane, and water transmittances are each
C    checked individually to be sure they give 0 < trans < 1.
C
C    ===================================================================
C    Loops downward over all the layers for each of the NCHN3 channels
C    to compute the layer transmittances TAU.
C
C    The water continuum absorption coefficient is:
C       k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
C
C    The layer effective fixed gas absorption coefficient is:
C       k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
C
C    The layer effective methane absorption coefficient is:
C       k_methane = the sum i=1 to 9 of { COEF(5+8+i)*OPRED(i) }
C
C    The layer effective water lines absorption coefficient is:
C       k_water = the sum i=1 to 11 of { COEF(5+8+9+i)*WPRED(i) }
C
C    where
C      "COEF" are the fast transmittance coefficients COEF3
C      "CONPRD" are the water continuum predictors CONPRD
C      "FPRED" are the fixed gases predictors FPRED3
C      "MPRED" are the methane predictors OPRED3
C      "WPRED" are the water lines predictors WPRED3
C
C    The total layer effective transmittance is:
C       TAU = exp( -[ k_con + k_fixed + k_methane + k_water])
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
C     3 Feb 1997 Scott Hannon   Re-wrote (from CALTAU) for FMW
C     3 Sep 1997 Scott Hannon   Added TAUZ and BLMULT
C     5 Mar 1998 Scott Hannon   Added OPTRAN water and deleted water
C                               preds 12 & 13
C     4 May 1998 Scott Hannon   Fix error: INDH2O(MXCHAN) not (MXCHNW)
C    26 Aug 1998 Scott Hannon   Fix mistake: loop on NLAY not MAXLAY;
C                               Add NLAY to call to CALOKW
C    11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
C    12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con


!END====================================================================

C      =================================================================
       SUBROUTINE CALT3 ( INDCHN, NLAY, BLMULT, NCHN3, CLIST3, COEF3,
     $    FIXMUL, CONPD3, FPRED3, MPRED3, WPRED3,
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
       INTEGER  NCHN3
       INTEGER CLIST3(MXCHN3)
       REAL  COEF3(N3COEF,MAXLAY,MXCHN3)
       REAL FIXMUL(MAXLAY)
       REAL CONPD3( N3CON,MAXLAY)
       REAL FPRED3( N3FIX,MAXLAY)
       REAL MPRED3( N3CH4,MAXLAY)
       REAL WPRED3( N3H2O,MAXLAY) 
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
       INTEGER      J
       INTEGER   ILAY
       REAL   KCON
       REAL   KFIX
       REAL   KMET
       REAL     KZ
       REAL  KZFMW
       REAL KLAYER
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
       DO I=1,NCHN3
C
C         Index for TAU
          J=INDCHN( CLIST3(I) )
C
C         -------------------------
C         Do OPTRAN water if needed
C         -------------------------
          IH2O=INDH2O( CLIST3(I) )
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
          KZFMW=0.0E+0
C
C         ------------------------------
C         Loop on layers (top to ground)
C         ------------------------------
          DO ILAY=1,NLAY
C
C            ---------------------------
C            Compute the water continuum
C            ---------------------------
             KCON=( COEF3(1,ILAY,I)*CONPD3(1,ILAY) ) +
     $            ( COEF3(2,ILAY,I)*CONPD3(2,ILAY) ) +
     $            ( COEF3(3,ILAY,I)*CONPD3(3,ILAY) ) +
     $            ( COEF3(4,ILAY,I)*CONPD3(4,ILAY) ) +
     $            ( COEF3(5,ILAY,I)*CONPD3(5,ILAY) ) +
     $            ( COEF3(6,ILAY,I)*CONPD3(6,ILAY) ) +
     $            ( COEF3(7,ILAY,I)*CONPD3(7,ILAY) )
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
             KFIX=( COEF3( 8,ILAY,I)*FPRED3(1,ILAY) ) +
     $            ( COEF3( 9,ILAY,I)*FPRED3(2,ILAY) ) +
     $            ( COEF3(10,ILAY,I)*FPRED3(3,ILAY) ) +
     $            ( COEF3(11,ILAY,I)*FPRED3(4,ILAY) ) +
     $            ( COEF3(12,ILAY,I)*FPRED3(5,ILAY) ) +
     $            ( COEF3(13,ILAY,I)*FPRED3(6,ILAY) ) +
     $            ( COEF3(14,ILAY,I)*FPRED3(7,ILAY) ) +
     $            ( COEF3(15,ILAY,I)*FPRED3(8,ILAY) )
C
             KFIX=KFIX*FIXMUL(ILAY)
C
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
C
C            ----------------------------
C            Compute the methane abs coef
C            ----------------------------
             KMET=( COEF3(16,ILAY,I)*MPRED3(1,ILAY) ) +
     $            ( COEF3(17,ILAY,I)*MPRED3(2,ILAY) ) +
     $            ( COEF3(18,ILAY,I)*MPRED3(3,ILAY) ) +
     $            ( COEF3(19,ILAY,I)*MPRED3(4,ILAY) ) +
     $            ( COEF3(20,ILAY,I)*MPRED3(5,ILAY) ) +
     $            ( COEF3(21,ILAY,I)*MPRED3(6,ILAY) ) +
     $            ( COEF3(22,ILAY,I)*MPRED3(7,ILAY) ) +
     $            ( COEF3(23,ILAY,I)*MPRED3(8,ILAY) ) +
     $            ( COEF3(24,ILAY,I)*MPRED3(9,ILAY) )
C
             IF (KMET .LT. 0.0E+0) THEN
                KMET=0.0E+0
             ELSEIF (KMET .GT. 1.0E+1) THEN
                KMET=1.0E+1
             ENDIF
C
C            --------------------------
C            Compute the water abs coef
C            --------------------------
             IF (LH2O) THEN
C               Not an OPTRAN water channel
                KW(ILAY)=
     $               ( COEF3(25,ILAY,I)*WPRED3( 1,ILAY) ) +
     $               ( COEF3(26,ILAY,I)*WPRED3( 2,ILAY) ) +
     $               ( COEF3(27,ILAY,I)*WPRED3( 3,ILAY) ) +
     $               ( COEF3(28,ILAY,I)*WPRED3( 4,ILAY) ) +
     $               ( COEF3(29,ILAY,I)*WPRED3( 5,ILAY) ) +
     $               ( COEF3(30,ILAY,I)*WPRED3( 6,ILAY) ) +
     $               ( COEF3(31,ILAY,I)*WPRED3( 7,ILAY) ) +
     $               ( COEF3(32,ILAY,I)*WPRED3( 8,ILAY) ) +
     $               ( COEF3(33,ILAY,I)*WPRED3( 9,ILAY) ) +
     $               ( COEF3(34,ILAY,I)*WPRED3(10,ILAY) ) +
     $               ( COEF3(35,ILAY,I)*WPRED3(11,ILAY) )
C
                IF (KW(ILAY) .LT. 0.0E+0) KW(ILAY)=0.0E+0
             ENDIF
C
C            Update KZFMW
             KZFMW=KZFMW + KFIX + KMET + KW(ILAY)
C
C            ----------------------------------
C            Calc the total layer transmittance
C            ----------------------------------
c
ccccc
c This block is usually commented out and is only uncommented for
c testing purposes.
c
c           kcon=0.0
c           kfix=0.0
c           kmet=0.0
c           kw(ilay)=0.0
ccccc
C
C            Calc the total layer effective optical depth
             KLAYER=KCON + KFIX + KMET + KW(ILAY)
C
C            Adjust the optical depth of the bottom layer
             IF (ILAY .EQ. NLAY) KLAYER=BLMULT*KLAYER
C
             KZ=KZ + KLAYER
             TAU(ILAY,J)=QIKEXP(-KLAYER)
C
          ENDDO
C         End loop on levels
C
C         Convert KZ to TAUZ
          TAUZ(J)=QIKEXP(-KZ)
C
       ENDDO
C      End loops on channel number (frequency)
C
       RETURN
       END
