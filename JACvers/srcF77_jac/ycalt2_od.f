C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    YCALT2 (for set2 = FOW) version with trace gases
C
!F77====================================================================


!ROUTINE NAME:
C    YCALT2


!ABSTRACT:
C    Calculate the transmittance for set2 using the predictors and the
C    fast transmittance coefficients.


!CALL PROTOCOL:
C    YCALT2 ( INDCHN, NLAY, NCHN2, CLIST2, COEF2, FIXMUL,
C       CONPD2, FPRED2, OPRED2, WPRED2, DPRED,  TRCPRD,
C       INDCO2, COFCO2, CO2MLT, INDSO2, SOFCO2, SO2MLT,
C       INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT,
C       INDNH3, COFNH3, NH3MLT, INDHDO, COFHDO, HDOMLT,
C       TAU, TAUZ )

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INT arr   INDCHN  channel indices             none
C    INTEGER   NLAY    number of layers to bottom  none
C    INTEGER   NCHN2   set2 number of channels     none
C    INT arr   CLIST2  set2 channel list           none
C    REAL arr  COEF2   set2 fast trans coefs       various
C    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
C    REAL arr  CONPD2  set2 H2O continuum preds    various
C    REAL arr  FPRED2  set2 fixed gases preds      various
C    REAL arr  OPRED2  set2 ozone predictors       various
C    REAL arr  WPRED2  set2 water predictors       various
C    REAL arr  DPRED   HDO predictors              various
C    REAL arr  TRCPRD  Trace gas pert predictors   various
C    INT arr   INDCO2  CO2 pert chan indices       none
C    REAL arr  COFCO2  CO2 pert coefs              various
C    REAL arr  CO2MLT  CO2 pert multiplier         none
C    INT arr   INDSO2  SO2 pert chan indices       none
C    REAL arr  COFSO2  SO2 pert coefs              various
C    REAL arr  SO2MLT  SO2 pert multiplier         none
C    INT arr   INDHNO  HNO3 pert chan indices      none
C    REAL arr  COFHNO  HNO3 pert coefs             various
C    REAL arr  HNOMLT  HNO3 pert multiplier        none
C    INT arr   INDN2O  N2O pert chan indices       none
C    REAL arr  COFN2O  N2O pert coefs              various
C    REAL arr  N2OMLT  N2O pert multiplier         none
C    INT arr   INDNH3  NH3 pert chan indices       none
C    REAL arr  COFNH3  NH3 pert coefs              various
C    REAL arr  NH3MLT  NH3 pert multiplier         none
C    INT arr   INDHDO  HDO pert chan indices       none
C    REAL arr  COFHDO  HDO pert coefs              various
C    REAL arr  HDOMLT  HDO pert multiplier         none

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
C    transmittances. Fixed, ozone, and water transmittances are each
C    checked individually to be sure they give 0 < trans < 1.
C
C    ===================================================================
C    Loops downward over all the layers for each of the NCHN2 channels
C    to compute the layer transmittances TAU.
C
C    The water continuum absorption coefficient is:
C       k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
C
C    The layer effective fixed gas absorption coefficient is:
C       k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
C
C    The layer effective ozone absorption coefficient is:
C       k_ozone = the sum i=1 to 10 of { COEF(5+8+i)*OPRED(i) }
C
C    The layer effective water lines absorption coefficient is:
C       k_water = the sum i=1 to 11 of { COEF(5+8+10+i)*WPRED(i) }
C
C    where
C      "COEF" are the fast transmittance coefficients COEF2
C      "CONPRD" are the water continuum predictors CONPRD
C      "FPRED" are the fixed gases predictors FPRED2
C      "OPRED" are the ozone predictors OPRED2
C      "WPRED" are the water lines predictors WPRED2
C
C    The total layer effective optical depth TAU is:
C       TAU = [ k_con + k_fixed + k_ozone + k_water ]
C
C    ===================================================================


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date        Programmer     Comments
C    ----------- -------------- ----------------------------------------
C    Dec  1 1994 Scott Hannon   Created
C     3 Feb 1997 Scott Hannon   Re-wrote (from YCALTAU) for FOW
C     3 Sep 1997 Scott Hannon   Added TAUZ and BLMULT
C    30 Sep 1997 Scott Hannon   Added variable CO2
C    11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
C    12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
C    18 May 2005 Scott Hannon   Add HNO3 based on SO2 code
C    28 Jun 2005 Scott Hannon   "trace" version for CO2,SO2,HNO3,N2O.
C    28 Mar 2006 Scott Hannon   Change TAU from trans to optical depth
C    22 Dec 2006 Scott Hannon   Change TAUZ from trans to optical depth
C                               and from (1 x n) to (m x n) array;
C                               delete func QIKEXP & argument BLMULT.
C    14 Sep 2010 Scott Hannon   Add 5th CO2 coef
C    10 May 2018 C Hepplewhite  Add NH3
C    1  Feb 2019 C Hepplewhite  Add HDO

!END====================================================================

C      =================================================================
       SUBROUTINE YCALT2 ( INDCHN, NLAY, NCHN2, CLIST2, COEF2,
     $    FIXMUL, CONPD2, FPRED2, OPRED2, WPRED2, DPRED, TRCPRD,
     $    INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $    INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT,
     $    INDNH3, COFNH3, NH3MLT, INDHDO, COFHDO, HDOMLT, TAU, TAUZ, IY, 
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,DTAU_DTZ,DTAU_DG1,DTAU_DG3)

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
       INTEGER  NCHN2
       INTEGER CLIST2(MXCHN2)
       REAL  COEF2(N2COEF,MAXLAY,MXCHN2)
       REAL FIXMUL(MAXLAY)
       REAL CONPD2( N2CON,MAXLAY)
       REAL FPRED2( N2FIX,MAXLAY)
       REAL OPRED2(  N2O3,MAXLAY)
       REAL WPRED2( N2H2O,MAXLAY)
       REAL DPRED(   NHDO,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2(  NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT(MAXLAY)
       INTEGER INDSO2(MXCHAN)
       REAL COFSO2(  NSO2,MAXLAY,MXCHNS)
       REAL SO2MLT(MAXLAY)
       INTEGER INDHNO(MXCHAN)
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH)
       REAL HNOMLT(MAXLAY)
       INTEGER INDN2O(MXCHAN)
       REAL COFN2O(  NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
       INTEGER INDNH3(MXCHAN)
       REAL COFNH3(  NNH3,MAXLAY,MXCHNA)
       REAL NH3MLT(MAXLAY)
       INTEGER INDHDO(MXCHAN)
       REAL COFHDO(  NHDO,MAXLAY,MXCHND)
       REAL HDOMLT(MAXLAY)
C
C      Output
       REAL    TAU(MAXLAY,MXCHAN)
       REAL   TAUZ(MAXLAY,MXCHAN)
       REAL DTAU_DTZ(MAXLAY,MXCHAN),DTAU_DG1(MAXLAY,MXCHAN),DTAU_DG3(MAXLAY,MXCHAN)

c input
       LOGICAL DOJAC
       INTEGER NWANTJ          ! number of wanted jacs (default 0=none)
       INTEGER  LISTJ(MAXPRO)  ! list of wanted channels
       !!! first index is the d/dT   second deriv is the d/dQ
       REAL CONJACPRD(3, N1CON,MAXLAY)
       REAL FJACPRED1(3, N1FIX,MAXLAY)
       REAL FJACPRED2(3, N2FIX,MAXLAY)
       REAL FJACPRED3(3, N3FIX,MAXLAY)
       REAL FJACPRED4(3, N4FIX,MAXLAY)
       REAL FJACPRED5(3, N5FIX,MAXLAY)
       REAL FJACPRED6(3, N6FIX,MAXLAY)
       REAL FJACPRED7(3, N7FIX,MAXLAY)
       REAL WJACPRED1(3, N1H2O,MAXLAY)
       REAL WJACPRED2(3, N2H2O,MAXLAY)
       REAL WJACPRED3(3, N3H2O,MAXLAY)
       REAL WJACPRED4(3, N4H2O,MAXLAY)
       REAL WJACPRED5(3, N5H2O,MAXLAY)
       REAL WJACPRED6(3, N6H2O,MAXLAY)
       REAL WJACPRED7(3, N7H2O,MAXLAY)
       REAL OJACPRED1(3,  N1O3,MAXLAY)
       REAL OJACPRED2(3,  N2O3,MAXLAY)
       REAL OJACPRED4(3,  N4O3,MAXLAY)
       REAL OJACPRED5(3,  N5O3,MAXLAY)
       REAL OJACPRED6(3,  N6O3,MAXLAY)
       REAL OJACPRED7(3,  N7O3,MAXLAY)
       REAL  DJACPRED(3,  NHDO,MAXLAY)
       REAL MJACPRED3(3, N3CH4,MAXLAY)
       REAL CJACPRED4(3,  N4CO,MAXLAY)
       REAL TRCJACPRD(3,NTRACE,MAXLAY)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      IWHICHJAC,ITRYJAC,INTERSECT
       INTEGER      I
       INTEGER   ICO2
       INTEGER  IHNO3
       INTEGER   ILAY
       INTEGER   IN2O
       INTEGER   INH3
       INTEGER   ISO2
       INTEGER   IHDO
       INTEGER      J
       REAL     DK
       REAL  DKCO2
       REAL DKHNO3
       REAL  DKN2O
       REAL  DKNH3
       REAL  DKSO2
       REAL  DKHDO
       REAL   KHDO
       REAL   KCON
       REAL   KFIX
       REAL KLAYER
       REAL   KOZO
       REAL   KWAT
       REAL     KZ
       LOGICAL   LCO2
       LOGICAL  LHNO3
       LOGICAL   LN2O
       LOGICAL   LNH3
       LOGICAL   LSO2
       LOGICAL   LHDO

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
C       DO I=1,NCHN2
C        DO I=IY,IY
         I = IY
C
C         Index for TAU
          J=INDCHN( CLIST2(I) )
C
C         Determine whether or not to do variable CO2
          ICO2=INDCO2( CLIST2(I) )
          IF (ICO2 .GT. 0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable CO2
          ISO2=INDSO2( CLIST2(I) )
          IF (ISO2 .GT. 0) THEN
             LSO2=.TRUE.
          ELSE
             LSO2=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable HNO3
          IHNO3=INDHNO( CLIST2(I) )
          IF (IHNO3 .GT. 0) THEN
             LHNO3=.TRUE.
          ELSE
             LHNO3=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable N2O
          IN2O=INDN2O( CLIST2(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable NH3
          INH3=INDNH3( CLIST2(I) )
          IF (INH3 .GT. 0) THEN
             LNH3=.TRUE.
          ELSE
             LNH3=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable HDO calc
          IHDO=INDHDO( CLIST2(I) )
          IF (IHDO .GT. 0) THEN
             LHDO=.TRUE.
          ELSE
             LHDO=.FALSE.
          ENDIF
C
c************************************************************************
C         Initialize the layer-to-space optical depth
          KZ=0.0E+0
C
C         ------------------------------
C         Loop on layers (top to ground)
C         ------------------------------
          DO ILAY=1,NLAY

C            ---------------------------
C            Compute the water continuum
C            ---------------------------
             KCON=( COEF2(1,ILAY,I)*CONPD2(1,ILAY) ) +
     $            ( COEF2(2,ILAY,I)*CONPD2(2,ILAY) ) +
     $            ( COEF2(3,ILAY,I)*CONPD2(3,ILAY) ) +
     $            ( COEF2(4,ILAY,I)*CONPD2(4,ILAY) ) +
     $            ( COEF2(5,ILAY,I)*CONPD2(5,ILAY) ) +
     $            ( COEF2(6,ILAY,I)*CONPD2(6,ILAY) ) +
     $            ( COEF2(7,ILAY,I)*CONPD2(7,ILAY) )
C
             IF (KCON .LT. 0.0+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 0.6E+0 .AND. ILAY .EQ. 1) THEN
c%%%%%                KCON=1.0E+1
                KCON=1.0E-10
             ENDIF
C

C            -----------------------------
C            Calc the fixed gases abs coef
C            -----------------------------
             KFIX=( COEF2( 8,ILAY,I)*FPRED2(1,ILAY) ) +
     $            ( COEF2( 9,ILAY,I)*FPRED2(2,ILAY) ) +
     $            ( COEF2(10,ILAY,I)*FPRED2(3,ILAY) ) +
     $            ( COEF2(11,ILAY,I)*FPRED2(4,ILAY) ) +
     $            ( COEF2(12,ILAY,I)*FPRED2(5,ILAY) ) +
     $            ( COEF2(13,ILAY,I)*FPRED2(6,ILAY) ) +
     $            ( COEF2(14,ILAY,I)*FPRED2(7,ILAY) ) +
     $            ( COEF2(15,ILAY,I)*FPRED2(8,ILAY) )
C
             KFIX=KFIX*FIXMUL(ILAY)
C
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 0.6E+0 .AND. ILAY .EQ. 1) THEN
c%%%%                KFIX=1.0E+1
                KFIX=1.0E-10
             ENDIF
C

C            --------------------------
C            Compute the ozone abs coef
C            --------------------------
             KOZO=( COEF2(16,ILAY,I)*OPRED2( 1,ILAY) ) +
     $            ( COEF2(17,ILAY,I)*OPRED2( 2,ILAY) ) +
     $            ( COEF2(18,ILAY,I)*OPRED2( 3,ILAY) ) +
     $            ( COEF2(19,ILAY,I)*OPRED2( 4,ILAY) ) +
     $            ( COEF2(20,ILAY,I)*OPRED2( 5,ILAY) ) +
     $            ( COEF2(21,ILAY,I)*OPRED2( 6,ILAY) ) +
     $            ( COEF2(22,ILAY,I)*OPRED2( 7,ILAY) ) +
     $            ( COEF2(23,ILAY,I)*OPRED2( 8,ILAY) ) +
     $            ( COEF2(24,ILAY,I)*OPRED2( 9,ILAY) ) +
     $            ( COEF2(25,ILAY,I)*OPRED2(10,ILAY) )
C
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+0) THEN
c%%%%%%                KOZO=1.0E+1
                KOZO=1.0E-10
             ENDIF
C

C            --------------------------
C            Compute the water abs coef
C            --------------------------
             KWAT=( COEF2(26,ILAY,I)*WPRED2( 1,ILAY) ) +
     $            ( COEF2(27,ILAY,I)*WPRED2( 2,ILAY) ) +
     $            ( COEF2(28,ILAY,I)*WPRED2( 3,ILAY) ) +
     $            ( COEF2(29,ILAY,I)*WPRED2( 4,ILAY) ) +
     $            ( COEF2(30,ILAY,I)*WPRED2( 5,ILAY) ) +
     $            ( COEF2(31,ILAY,I)*WPRED2( 6,ILAY) ) +
     $            ( COEF2(32,ILAY,I)*WPRED2( 7,ILAY) ) +
     $            ( COEF2(33,ILAY,I)*WPRED2( 8,ILAY) ) +
     $            ( COEF2(34,ILAY,I)*WPRED2( 9,ILAY) ) +
     $            ( COEF2(35,ILAY,I)*WPRED2(10,ILAY) ) +
     $            ( COEF2(36,ILAY,I)*WPRED2(11,ILAY) )
C
             IF (KWAT .LT. 0.0E+0) THEN
                KWAT=0.0E+0
             ELSEIF( KWAT .GT. 1.0E+0) THEN
c%%%%%%%%                KWAT=1.0E+1
                KWAT=1.0E-10
             ENDIF
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
c           kozo=0.0E+0
c           kwat=0.0E+0
ccccc
C
C            ----------------------------
C            Calc change in total optical
C            depth due to variable CO2
C            ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) +
     $                ( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) +
     $                ( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) +
     $                ( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) +
     $                ( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
C

C            ----------------------------
C            Calc change in total optical
C            depth due to variable SO2
C            ----------------------------
             IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0) THEN
                DKSO2=( COFSO2(1,ILAY,ISO2)*TRCPRD(1,ILAY) ) +
     $                ( COFSO2(2,ILAY,ISO2)*TRCPRD(2,ILAY) ) +
     $                ( COFSO2(3,ILAY,ISO2)*TRCPRD(3,ILAY) ) +
     $                ( COFSO2(4,ILAY,ISO2)*TRCPRD(4,ILAY) )
                DKSO2=DKSO2*SO2MLT(ILAY)
             ELSE
                DKSO2=0.0
             ENDIF
C

C            ----------------------------
C            Calc change in total optical
C            depth due to variable HNO3
C            ----------------------------
             IF (LHNO3 .AND. HNOMLT(ILAY) .NE. 0) THEN
                DKHNO3=( COFHNO(1,ILAY,IHNO3)*TRCPRD(1,ILAY) ) +
     $                 ( COFHNO(2,ILAY,IHNO3)*TRCPRD(2,ILAY) ) +
     $                 ( COFHNO(3,ILAY,IHNO3)*TRCPRD(3,ILAY) ) +
     $                 ( COFHNO(4,ILAY,IHNO3)*TRCPRD(4,ILAY) )
                DKHNO3=DKHNO3*HNOMLT(ILAY)
             ELSE
                DKHNO3=0.0
             ENDIF
C

C            ----------------------------
C            Calc change in total optical
C            depth due to variable N2O
C            ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) +
     $                ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) +
     $                ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) +
     $                ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) +
     $                ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) +
     $                ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) +
     $                ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
C
C            ----------------------------
C            Calc change in total optical
C            depth due to variable NH3
C            ----------------------------
             IF (LNH3 .AND. NH3MLT(ILAY) .NE. 0) THEN
                DKNH3=( COFNH3(1,ILAY,INH3)*TRCPRD(1,ILAY) ) +
     $                ( COFNH3(2,ILAY,INH3)*TRCPRD(2,ILAY) ) +
     $                ( COFNH3(3,ILAY,INH3)*TRCPRD(3,ILAY) ) +
     $                ( COFNH3(4,ILAY,INH3)*TRCPRD(4,ILAY) )
                DKNH3=DKNH3*NH3MLT(ILAY)
             ELSE
                DKNH3=0.0
             ENDIF
C
C            ------------------------------------------
C            Calc total optical depth and transmittance
C            ------------------------------------------
C            Calc total layer optical depth
ccc
c this block for testing
c      DKCO2=0.0
c      DKSO2=0.0
c      DKHNO3=0.0
c      DKN2O=0.0
C       DKNH3=0.0
C       DKHDO=0.0
       KHDO=0.0
ccc
C            Limit -DK so it can never totally totally cancel KFIX
             DK = DKCO2 + DKSO2 + DKHNO3 + DKN2O + DKNH3
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF

C            Calc effective layer optical depth
             KLAYER = KCON + KFIX + KOZO + KWAT + DK
             TAU(ILAY,J)=KLAYER
C
C            Calc layer-to-space optical depth
             KZ=KZ + KLAYER
             TAUZ(ILAY,J)=KZ
C
          ENDDO
C         End loop on levels
C
C       ENDDO    !DO I = IY,IY
C      End loops on channel number (frequency)
C

c************************************************************************
          IF (DOJAC) THEN
            DO IWHICHJAC = 1,3
              IF (IWHICHJAC .EQ. 1) THEN 
                ITRYJAC = 100  !! TZ
              ELSEIF (IWHICHJAC .EQ. 2) THEN 
                ITRYJAC = 1    !! GID 1
              ELSEIF (IWHICHJAC .EQ. 3) THEN 
                ITRYJAC = 3    !! GID 3
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

C            ---------------------------
C            Compute the water continuum
C            ---------------------------
             KCON=( COEF2(1,ILAY,I)*CONJACPRD(IWHICHJAC,1,ILAY) ) +
     $            ( COEF2(2,ILAY,I)*CONJACPRD(IWHICHJAC,2,ILAY) ) +
     $            ( COEF2(3,ILAY,I)*CONJACPRD(IWHICHJAC,3,ILAY) ) +
     $            ( COEF2(4,ILAY,I)*CONJACPRD(IWHICHJAC,4,ILAY) ) +
     $            ( COEF2(5,ILAY,I)*CONJACPRD(IWHICHJAC,5,ILAY) ) +
     $            ( COEF2(6,ILAY,I)*CONJACPRD(IWHICHJAC,6,ILAY) ) +
     $            ( COEF2(7,ILAY,I)*CONJACPRD(IWHICHJAC,7,ILAY) )
C
c             IF (KCON .LT. 0.0+0) THEN
c                KCON=0.0E+0
c             ELSEIF (KCON .GT. 0.6E+0 .AND. ILAY .EQ. 1) THEN
cc%%%%%                KCON=1.0E+1
c                KCON=1.0E-10
c             ENDIF
C

C            -----------------------------
C            Calc the fixed gases abs coef
C            -----------------------------
             KFIX=( COEF2( 8,ILAY,I)*FJACPRED2(IWHICHJAC,1,ILAY) ) +
     $            ( COEF2( 9,ILAY,I)*FJACPRED2(IWHICHJAC,2,ILAY) ) +
     $            ( COEF2(10,ILAY,I)*FJACPRED2(IWHICHJAC,3,ILAY) ) +
     $            ( COEF2(11,ILAY,I)*FJACPRED2(IWHICHJAC,4,ILAY) ) +
     $            ( COEF2(12,ILAY,I)*FJACPRED2(IWHICHJAC,5,ILAY) ) +
     $            ( COEF2(13,ILAY,I)*FJACPRED2(IWHICHJAC,6,ILAY) ) +
     $            ( COEF2(14,ILAY,I)*FJACPRED2(IWHICHJAC,7,ILAY) ) +
     $            ( COEF2(15,ILAY,I)*FJACPRED2(IWHICHJAC,8,ILAY) )
C
             KFIX=KFIX*FIXMUL(ILAY)
C
c             IF (KFIX .LT. 0.0E+0) THEN
c                KFIX=0.0E+0
c             ELSEIF (KFIX .GT. 0.6E+0 .AND. ILAY .EQ. 1) THEN
cc%%%%                KFIX=1.0E+1
c                KFIX=1.0E-10
c             ENDIF
C

C            --------------------------
C            Compute the ozone abs coef
C            --------------------------
             KOZO=( COEF2(16,ILAY,I)*OJACPRED2(IWHICHJAC, 1,ILAY) ) +
     $            ( COEF2(17,ILAY,I)*OJACPRED2(IWHICHJAC, 2,ILAY) ) +
     $            ( COEF2(18,ILAY,I)*OJACPRED2(IWHICHJAC, 3,ILAY) ) +
     $            ( COEF2(19,ILAY,I)*OJACPRED2(IWHICHJAC, 4,ILAY) ) +
     $            ( COEF2(20,ILAY,I)*OJACPRED2(IWHICHJAC, 5,ILAY) ) +
     $            ( COEF2(21,ILAY,I)*OJACPRED2(IWHICHJAC, 6,ILAY) ) +
     $            ( COEF2(22,ILAY,I)*OJACPRED2(IWHICHJAC, 7,ILAY) ) +
     $            ( COEF2(23,ILAY,I)*OJACPRED2(IWHICHJAC, 8,ILAY) ) +
     $            ( COEF2(24,ILAY,I)*OJACPRED2(IWHICHJAC, 9,ILAY) ) +
     $            ( COEF2(25,ILAY,I)*OJACPRED2(IWHICHJAC,10,ILAY) )
C              IF ((ILAY .EQ. 1) .AND. (IWHICHJAC .EQ. 1)) write(*,'(A,10(ES12.5))') 'OJACPRED2',OJACPRED2(IWHICHJAC,1:10,ILAY)
C
c             IF (KOZO .LT. 0.0E+0) THEN
c                KOZO=0.0E+0
c             ELSEIF (KOZO .GT. 1.0E+0) THEN
cc%%%%%%                KOZO=1.0E+1
c                KOZO=1.0E-10
c             ENDIF
C

C            --------------------------
C            Compute the water abs coef
C            --------------------------
             KWAT=( COEF2(26,ILAY,I)*WJACPRED2(IWHICHJAC, 1,ILAY) ) +
     $            ( COEF2(27,ILAY,I)*WJACPRED2(IWHICHJAC, 2,ILAY) ) +
     $            ( COEF2(28,ILAY,I)*WJACPRED2(IWHICHJAC, 3,ILAY) ) +
     $            ( COEF2(29,ILAY,I)*WJACPRED2(IWHICHJAC, 4,ILAY) ) +
     $            ( COEF2(30,ILAY,I)*WJACPRED2(IWHICHJAC, 5,ILAY) ) +
     $            ( COEF2(31,ILAY,I)*WJACPRED2(IWHICHJAC, 6,ILAY) ) +
     $            ( COEF2(32,ILAY,I)*WJACPRED2(IWHICHJAC, 7,ILAY) ) +
     $            ( COEF2(33,ILAY,I)*WJACPRED2(IWHICHJAC, 8,ILAY) ) +
     $            ( COEF2(34,ILAY,I)*WJACPRED2(IWHICHJAC, 9,ILAY) ) +
     $            ( COEF2(35,ILAY,I)*WJACPRED2(IWHICHJAC,10,ILAY) ) +
     $            ( COEF2(36,ILAY,I)*WJACPRED2(IWHICHJAC,11,ILAY) )
C
c             IF (KWAT .LT. 0.0E+0) THEN
c                KWAT=0.0E+0
c             ELSEIF( KWAT .GT. 1.0E+0) THEN
cc%%%%%%%%                KWAT=1.0E+1
c                KWAT=1.0E-10
c             ENDIF
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
c           kozo=0.0E+0
c           kwat=0.0E+0
ccccc
C
C            ----------------------------
C            Calc change in total optical
C            depth due to variable CO2
C            ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,1,ILAY) ) +
     $                ( COFCO2(2,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,2,ILAY) ) +
     $                ( COFCO2(3,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,3,ILAY) ) +
     $                ( COFCO2(4,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,4,ILAY) ) +
     $                ( COFCO2(5,ILAY,ICO2)*TRCJACPRD(IWHICHJAC,5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
C

C            ----------------------------
C            Calc change in total optical
C            depth due to variable SO2
C            ----------------------------
             IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0) THEN
                DKSO2=( COFSO2(1,ILAY,ISO2)*TRCJACPRD(IWHICHJAC,1,ILAY) ) +
     $                ( COFSO2(2,ILAY,ISO2)*TRCJACPRD(IWHICHJAC,2,ILAY) ) +
     $                ( COFSO2(3,ILAY,ISO2)*TRCJACPRD(IWHICHJAC,3,ILAY) ) +
     $                ( COFSO2(4,ILAY,ISO2)*TRCJACPRD(IWHICHJAC,4,ILAY) )
                DKSO2=DKSO2*SO2MLT(ILAY)
             ELSE
                DKSO2=0.0
             ENDIF
C

C            ----------------------------
C            Calc change in total optical
C            depth due to variable HNO3
C            ----------------------------
             IF (LHNO3 .AND. HNOMLT(ILAY) .NE. 0) THEN
                DKHNO3=( COFHNO(1,ILAY,IHNO3)*TRCJACPRD(IWHICHJAC,1,ILAY) ) +
     $                 ( COFHNO(2,ILAY,IHNO3)*TRCJACPRD(IWHICHJAC,2,ILAY) ) +
     $                 ( COFHNO(3,ILAY,IHNO3)*TRCJACPRD(IWHICHJAC,3,ILAY) ) +
     $                 ( COFHNO(4,ILAY,IHNO3)*TRCJACPRD(IWHICHJAC,4,ILAY) )
                DKHNO3=DKHNO3*HNOMLT(ILAY)
             ELSE
                DKHNO3=0.0
             ENDIF
C

C            ----------------------------
C            Calc change in total optical
C            depth due to variable N2O
C            ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0) THEN
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
C
C            ----------------------------
C            Calc change in total optical
C            depth due to variable NH3
C            ----------------------------
             IF (LNH3 .AND. NH3MLT(ILAY) .NE. 0) THEN
                DKNH3=( COFNH3(1,ILAY,INH3)*TRCJACPRD(IWHICHJAC,1,ILAY) ) +
     $                ( COFNH3(2,ILAY,INH3)*TRCJACPRD(IWHICHJAC,2,ILAY) ) +
     $                ( COFNH3(3,ILAY,INH3)*TRCJACPRD(IWHICHJAC,3,ILAY) ) +
     $                ( COFNH3(4,ILAY,INH3)*TRCJACPRD(IWHICHJAC,4,ILAY) )
                DKNH3=DKNH3*NH3MLT(ILAY)
             ELSE
                DKNH3=0.0
             ENDIF
C
C            ------------------------------------------
C            Calc total optical depth and transmittance
C            ------------------------------------------
C            Calc total layer optical depth
ccc
c this block for testing
c      DKCO2=0.0
c      DKSO2=0.0
c      DKHNO3=0.0
c      DKN2O=0.0
C       DKNH3=0.0
C       DKHDO=0.0
       KHDO=0.0
ccc
C            Limit -DK so it can never totally totally cancel KFIX
c             DK = DKCO2 + DKSO2 + DKHNO3 + DKN2O + DKNH3
c             IF (-DK .GE. KFIX) THEN
c                DK = -0.999*KFIX
c             ENDIF

C            Calc effective layer optical depth
             KLAYER = KCON + KFIX + KOZO + KWAT + DK

                   IF (IWHICHJAC .EQ. 1) THEN 
                     DTAU_DTZ(ILAY,J)=KLAYER
C                     IF (IY .EQ. 18) THEN
C                       write(*,'(A,2(I4),7(ES12.5))') 'moo',ILAY,IY,FIXMUL(ILAY),KCON,KFIX,KOZO,KWAT,DK,KLAYER
C                     END IF
                   ELSEIF (IWHICHJAC .EQ. 2) THEN 
                     DTAU_DG1(ILAY,J)=KLAYER
                   ELSEIF (IWHICHJAC .EQ. 3) THEN 
                     DTAU_DG3(ILAY,J)=KLAYER
                   END IF

c             TAU(ILAY,J)=KLAYER
C
C            Calc layer-to-space optical depth
c             KZ=KZ + KLAYER
c             TAUZ(ILAY,J)=KZ
C
          ENDDO
C         End loop on levels
C
C       ENDDO    !DO I = IY,IY
C      End loops on channel number (frequency)
C
cbaba
              END IF   !!!! IF (INTERSECT(ITRYJAC,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) : is this a jac to work on?????
            END DO     !!! DO IWHICHJAC = 1,3
          END IF  !! DOJAC

c************************************************************************

       RETURN
       END
