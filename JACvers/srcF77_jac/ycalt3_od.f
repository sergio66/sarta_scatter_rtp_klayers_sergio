C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    YCALT3 (for set3 = FMW) version with trace gases (no CO2)
C
!F77====================================================================


!ROUTINE NAME:
C    YCALT3


!ABSTRACT:
C    Calculate the transmittance for set3 using the prdictor and the
C    fast transmittance coefficients.


!CALL PROTOCOL:
C    YCALT3 ( INDCHN, NLAY, NCHN3, CLIST3, COEF3,
C       FIXMUL, CONPD3, FPRED3, MPRED3, WPRED3, DPRED, TRCPRD,
C       INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT,
C       INDN2O, COFN2O, N2OMLT, 
C       INDNH3, COFNH3, NH3MLT, INDHDO, COFHDO, HDOMLT, 
C       INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX,
C       LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INT arr   INDCHN  channel indices             none
C    INTEGER   NLAY    number of layers to bottom  none
C    INTEGER   NCHN3   set3 number of channels     none
C    INT arr   CLIST3  set3 channel list           none
C    REAL arr  COEF3   set3 fast trans coefs       various
C    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
C    REAL arr  CONPD3  set3 H2O continuum preds    various
C    REAL arr  FPRED3  set3 fixed gases preds      various
C    REAL arr  MPRED3  set3 methane predictors     various
C    REAL arr  WPRED3  set3 water predictors       various
C    REAL arr  DPRED   HDO predictors              various
C    REAL arr  TRCPRD  trace gas pert predictors   various
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
C    INT arr   INDH2O  OPTRAN H2O chan indices     none
C    REAL arr  H2OPRD  OPTRAN H2O predictors       various
C    REAL arr  COFH2O  OPTRAN H2O coefs            various
C    INTEGER   LOPMAX  OPTRAN max level            none
C    INTEGER   LOPLOW  OPTRAN low bracketing level none
C    LOG arr   LOPUSE  OPTRAN level needed?        none
C    REAL arr  WAOP    OPTRAN layer water amounts  kilomoles/cm^2
C    REAL arr  DAOP    OPTRAN-to-AIRS interp fac   none
C    REAL arr  WAANG   AIRS layer water amounts    kilomoles/cm^2


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
C    The total layer effective optical depth is:
C       TAU = [ k_con + k_fixed + k_methane + k_water ]
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
C     3 Feb 1997 Scott Hannon   Re-wrote (from YCALTAU) for FMW
C     3 Sep 1997 Scott Hannon   Added TAUZ and BLMULT
C     5 Mar 1998 Scott Hannon   Added OPTRAN water and deleted water
C                               preds 12 & 13
C     4 May 1998 Scott Hannon   Fix error: INDH2O(MXCHAN) not (MXCHNW)
C    26 Aug 1998 Scott Hannon   Fix mistake: loop on NLAY not MAXLAY;
C                               Add NLAY to call to CALOKW
C    11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
C    12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
C    25 Apr 2003 Scott Hannon   Add HNO3 based on SO2 code
C    28 Jun 2005 Scott Hannon   "trace" version for SO2,HNO3,N2O
C    28 Mar 2006 Scott Hannon   Change TAU from trans to optical depth
C    22 Dec 2006 Scott Hannon   Change TAUZ from trans to optical depth
C                               and from (1 x n) to (m x n) array;
C                               delete func QIKEXP & argument BLMULT.
C    10 May 2018 C Hepplewhite  Add NH3
C    1  Feb 2019 C Hepplewhite  Add HDO

!END====================================================================

C      =================================================================
       SUBROUTINE YCALT3 ( INDCHN, NLAY, NCHN3, CLIST3, COEF3,
     $    FIXMUL, CONPD3, FPRED3, MPRED3, WPRED3, DPRED, TRCPRD,
     $    INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT,
     $    INDN2O, COFN2O, N2OMLT, INDNH3, COFNH3, NH3MLT,
     $    INDHDO, COFHDO, HDOMLT, INDH2O, H2OPRD, COFH2O,
     $    LOPMIN, LOPMAX, LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ, 
     $    IY,
     $    DOJAC,LISTJ,NWANTJ,SECANG,CONJACPRD,DJACPRED,H2OJACPRD,DAOPJAC,
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
       INTEGER  NCHN3
       INTEGER CLIST3(MXCHN3)
       REAL SECANG(MAXLAY)
       REAL  COEF3(N3COEF,MAXLAY,MXCHN3)
       REAL FIXMUL(MAXLAY)
       REAL CONPD3( N3CON,MAXLAY)
       REAL FPRED3( N3FIX,MAXLAY)
       REAL MPRED3( N3CH4,MAXLAY)
       REAL WPRED3( N3H2O,MAXLAY) 
       REAL DPRED(   NHDO,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
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
c optran
       REAL   H2OJACPRD(OPTRANJAC,NH2O,MXOWLY)
       REAL  DAOPJAC(OPTRANJAC, MAXLAY)           !! OPTRAN, used by ycalt1_od, ycalt3_od
       REAL  KW_T(MAXLAY)
       REAL  KW_1(MAXLAY)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      IWHICHJAC,ITRYJAC,INTERSECT
       INTEGER      I
       INTEGER  IHNO3
       INTEGER   IN2O
       INTEGER   ILAY
       INTEGER   ISO2
       INTEGER   INH3
       INTEGER   IHDO
       INTEGER      J
       REAL     DK
C       REAL  DKCO2, QDKCO2
       REAL DKHNO3, QDKHNO3
       REAL  DKN2O, QDKN2O
       REAL  DKSO2, QDKSO2
       REAL  DKNH3, QDKNH3
       REAL  DKHDO, QDKHDO
       REAL   KHDO
       REAL   KCON
       REAL   KFIX
       REAL   KMET
       REAL     KZ
       REAL  KZFMW
       REAL KLAYER
       LOGICAL   LH2O
       LOGICAL  LHNO3
       LOGICAL   LN2O
       LOGICAL   LSO2
       LOGICAL   LNH3
       LOGICAL   LHDO
C
C      for CALOKW
       INTEGER   IH2O
       REAL     KW(MAXLAY)
      
       REAL f1,f2,f3,f4,f5

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
C       DO I=1,NCHN3
C        DO I=IY,IY
         I = IY
C
C         Index for TAU
          J=INDCHN( CLIST3(I) )
C
C
C         Determine whether or not to do variable SO2
          ISO2=INDSO2( CLIST3(I) )
          IF (ISO2 .GT. 0) THEN
             LSO2=.TRUE.
          ELSE
             LSO2=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable HNO3
          IHNO3=INDHNO( CLIST3(I) )
          IF (IHNO3 .GT. 0) THEN
             LHNO3=.TRUE.
          ELSE
             LHNO3=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable N2O
          IN2O=INDN2O( CLIST3(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable NH3
          INH3=INDNH3( CLIST3(I) )
          IF (INH3 .GT. 0) THEN
             LNH3=.TRUE.
          ELSE
             LNH3=.FALSE.
          ENDIF
C
C         Determine whether or not to do variable HDO calc
          IHDO=INDHDO( CLIST3(I) )
          IF (IHDO .GT. 0) THEN
             LHDO=.TRUE.
          ELSE
             LHDO=.FALSE.
          ENDIF
C
C         -------------------------
C         Do OPTRAN water if needed
C         -------------------------
          IH2O=INDH2O( CLIST3(I) )
          IF (IH2O .GT. 0) THEN
             LH2O=.FALSE.
C            Calc OPTRAN water
C
             CALL YCALOKW( NLAY, IH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE,
     $          H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW, 
     $          DOJAC, SECANG, H2OJACPRD, DAOPJAC, KW_T, KW_1 )
C
          ELSE
             LH2O=.TRUE.
          ENDIF
C
c************************************************************************
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
C            --------------------------
C            Compute the HDO abs coef
C            --------------------------
             IF (LHDO) THEN
                KHDO=( COFHDO(1,ILAY,IHDO)*DPRED( 1,ILAY) ) +
     $               ( COFHDO(2,ILAY,IHDO)*DPRED( 2,ILAY) ) +
     $               ( COFHDO(3,ILAY,IHDO)*DPRED( 3,ILAY) ) +
     $               ( COFHDO(4,ILAY,IHDO)*DPRED( 4,ILAY) ) +
     $               ( COFHDO(5,ILAY,IHDO)*DPRED( 5,ILAY) ) +
     $               ( COFHDO(6,ILAY,IHDO)*DPRED( 6,ILAY) ) +
     $               ( COFHDO(7,ILAY,IHDO)*DPRED( 7,ILAY) ) +
     $               ( COFHDO(8,ILAY,IHDO)*DPRED( 8,ILAY) )
C     $               ( COFHDO(9,ILAY,IHDO)*DPRED( 9,ILAY) ) +
C     $               ( COFHDO(10,ILAY,IHDO)*DPRED(10,ILAY) ) +
C     $               ( COFHDO(11,ILAY,IHDO)*DPRED(11,ILAY) )
C
C                IF (KHDO .LT. 0.0E+0) KHDO=0.0E+0
                KHDO=KHDO*HDOMLT(ILAY)
             ELSE
                KHDO=0.0
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
c           kcon=0.0
c           kfix=0.0
c           kmet=0.0
c           kw(ilay)=0.0
ccccc
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
ccc
c this block for testing
c      DKSO2=0.0
c      DKHNO3=0.0
c      DKN2O=0.0
C       DKNH3=0.0
C       DKHDO=0.0
      KHDO=0.0
ccc
C            Limit -DK so it can never totally totally cancel KFIX
             DK = DKSO2 + DKHNO3 + DKN2O + DKNH3
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF

C            Calc total layer optical depth
             KLAYER = KCON + KFIX + KMET + KW(ILAY) + DK
cDEBUG             IF (LH2O) THEN
cDEBUG               !! line
cDEBUG               print *,1,ILAY,KW(ILAY),KLAYER      !! WATER DEBUG, see test_wv_6um_jacs.m
cDEBUG             ELSE
cDEBUG               !! optran
cDEBUG               print *,2,ILAY,KW(ILAY),KLAYER      !! WATER DEBUG, see test_wv_6um_jacs.m
cDEBUG             END IF
             TAU(ILAY,J)=KLAYER
C
C            Calc layer-to-space optical depth
             KZ=KZ + KLAYER
             TAUZ(ILAY,J)=KZ
C
          ENDDO
C         End loop on levels
C
C       ENDDO     ! DO I = IY,IY
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
             KCON=( COEF3(1,ILAY,I)*CONJACPRD(IWHICHJAC,1,ILAY) ) +
     $            ( COEF3(2,ILAY,I)*CONJACPRD(IWHICHJAC,2,ILAY) ) +
     $            ( COEF3(3,ILAY,I)*CONJACPRD(IWHICHJAC,3,ILAY) ) +
     $            ( COEF3(4,ILAY,I)*CONJACPRD(IWHICHJAC,4,ILAY) ) +
     $            ( COEF3(5,ILAY,I)*CONJACPRD(IWHICHJAC,5,ILAY) ) +
     $            ( COEF3(6,ILAY,I)*CONJACPRD(IWHICHJAC,6,ILAY) ) +
     $            ( COEF3(7,ILAY,I)*CONJACPRD(IWHICHJAC,7,ILAY) )
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
             KFIX=( COEF3( 8,ILAY,I)*FJACPRED3(IWHICHJAC,1,ILAY) ) +
     $            ( COEF3( 9,ILAY,I)*FJACPRED3(IWHICHJAC,2,ILAY) ) +
     $            ( COEF3(10,ILAY,I)*FJACPRED3(IWHICHJAC,3,ILAY) ) +
     $            ( COEF3(11,ILAY,I)*FJACPRED3(IWHICHJAC,4,ILAY) ) +
     $            ( COEF3(12,ILAY,I)*FJACPRED3(IWHICHJAC,5,ILAY) ) +
     $            ( COEF3(13,ILAY,I)*FJACPRED3(IWHICHJAC,6,ILAY) ) +
     $            ( COEF3(14,ILAY,I)*FJACPRED3(IWHICHJAC,7,ILAY) ) +
     $            ( COEF3(15,ILAY,I)*FJACPRED3(IWHICHJAC,8,ILAY) )
C
             KFIX=KFIX*FIXMUL(ILAY)
C
c             IF (KFIX .LT. 0.0E+0) THEN
c                KFIX=0.0E+0
c             ELSEIF (KFIX .GT. 1.0E+1) THEN
c                KFIX=1.0E+1
c             ENDIF
C
C            ----------------------------
C            Compute the methane abs coef
C            ----------------------------
             KMET=( COEF3(16,ILAY,I)*MJACPRED3(IWHICHJAC,1,ILAY) ) +
     $            ( COEF3(17,ILAY,I)*MJACPRED3(IWHICHJAC,2,ILAY) ) +
     $            ( COEF3(18,ILAY,I)*MJACPRED3(IWHICHJAC,3,ILAY) ) +
     $            ( COEF3(19,ILAY,I)*MJACPRED3(IWHICHJAC,4,ILAY) ) +
     $            ( COEF3(20,ILAY,I)*MJACPRED3(IWHICHJAC,5,ILAY) ) +
     $            ( COEF3(21,ILAY,I)*MJACPRED3(IWHICHJAC,6,ILAY) ) +
     $            ( COEF3(22,ILAY,I)*MJACPRED3(IWHICHJAC,7,ILAY) ) +
     $            ( COEF3(23,ILAY,I)*MJACPRED3(IWHICHJAC,8,ILAY) ) +
     $            ( COEF3(24,ILAY,I)*MJACPRED3(IWHICHJAC,9,ILAY) )

c              IF ((ILAY .EQ. 1) .AND. (IWHICHJAC .EQ. 7)) print *,MJACPRED3(IWHICHJAC,1:9,ILAY)
C
c             IF (KMET .LT. 0.0E+0) THEN
c                KMET=0.0E+0
c             ELSEIF (KMET .GT. 1.0E+1) THEN
c                KMET=1.0E+1
c             ENDIF
C
C            --------------------------
C            Compute the water abs coef
C            --------------------------
             IF (LH2O) THEN
C               Not an OPTRAN water channel
                KW(ILAY)=
     $               ( COEF3(25,ILAY,I)*WJACPRED3(IWHICHJAC, 1,ILAY) ) +
     $               ( COEF3(26,ILAY,I)*WJACPRED3(IWHICHJAC, 2,ILAY) ) +
     $               ( COEF3(27,ILAY,I)*WJACPRED3(IWHICHJAC, 3,ILAY) ) +
     $               ( COEF3(28,ILAY,I)*WJACPRED3(IWHICHJAC, 4,ILAY) ) +
     $               ( COEF3(29,ILAY,I)*WJACPRED3(IWHICHJAC, 5,ILAY) ) +
     $               ( COEF3(30,ILAY,I)*WJACPRED3(IWHICHJAC, 6,ILAY) ) +
     $               ( COEF3(31,ILAY,I)*WJACPRED3(IWHICHJAC, 7,ILAY) ) +
     $               ( COEF3(32,ILAY,I)*WJACPRED3(IWHICHJAC, 8,ILAY) ) +
     $               ( COEF3(33,ILAY,I)*WJACPRED3(IWHICHJAC, 9,ILAY) ) +
     $               ( COEF3(34,ILAY,I)*WJACPRED3(IWHICHJAC,10,ILAY) ) +
     $               ( COEF3(35,ILAY,I)*WJACPRED3(IWHICHJAC,11,ILAY) )
c                IF (KW(ILAY) .LT. 0.0E+0) KW(ILAY)=0.0E+0
ccDEBUG                print *,101,ILAY,KW(ILAY),0         !!!! WATER DEBUG, see test_wv_6um_jacs.m
ccDEBUG              ELSE
ccDEBUG                print *,102,ILAY,KW_1(ILAY),0         !!!! WATER DEBUG, see test_wv_6um_jacs.m
             ENDIF
C
C            Update KZFMW
             KZFMW=KZFMW + KFIX + KMET + KW(ILAY)
C
C            --------------------------
C            Compute the HDO abs coef
C            --------------------------
             IF (LHDO) THEN
                KHDO=( COFHDO(1,ILAY,IHDO)*DJACPRED(IWHICHJAC, 1,ILAY) ) +
     $               ( COFHDO(2,ILAY,IHDO)*DJACPRED(IWHICHJAC, 2,ILAY) ) +
     $               ( COFHDO(3,ILAY,IHDO)*DJACPRED(IWHICHJAC, 3,ILAY) ) +
     $               ( COFHDO(4,ILAY,IHDO)*DJACPRED(IWHICHJAC, 4,ILAY) ) +
     $               ( COFHDO(5,ILAY,IHDO)*DJACPRED(IWHICHJAC, 5,ILAY) ) +
     $               ( COFHDO(6,ILAY,IHDO)*DJACPRED(IWHICHJAC, 6,ILAY) ) +
     $               ( COFHDO(7,ILAY,IHDO)*DJACPRED(IWHICHJAC, 7,ILAY) ) +
     $               ( COFHDO(8,ILAY,IHDO)*DJACPRED(IWHICHJAC, 8,ILAY) )
C     $               ( COFHDO(9,ILAY,IHDO)*DJACPRED(IWHICHJAC, 9,ILAY) ) +
C     $               ( COFHDO(10,ILAY,IHDO)*DJACPRED(IWHICHJAC,10,ILAY) ) +
C     $               ( COFHDO(11,ILAY,IHDO)*DJACPRED(IWHICHJAC,11,ILAY) )
C
C                IF (KHDO .LT. 0.0E+0) KHDO=0.0E+0
                KHDO=KHDO*HDOMLT(ILAY)
c             ELSE
c                KHDO=0.0
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
c           kcon=0.0
c           kfix=0.0
c           kmet=0.0
c           kw(ilay)=0.0
ccccc
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

             IF (LSO2) THEN
                QDKSO2=( COFSO2(1,ILAY,ISO2)*TRCPRD(1,ILAY) ) +
     $                ( COFSO2(2,ILAY,ISO2)*TRCPRD(2,ILAY) ) +
     $                ( COFSO2(3,ILAY,ISO2)*TRCPRD(3,ILAY) ) +
     $                ( COFSO2(4,ILAY,ISO2)*TRCPRD(4,ILAY) )
                QDKSO2=QDKSO2*SO2JACMLT(ILAY)
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

             IF (LHNO3) THEN
                QDKHNO3=( COFHNO(1,ILAY,IHNO3)*TRCPRD(1,ILAY) ) +
     $                 ( COFHNO(2,ILAY,IHNO3)*TRCPRD(2,ILAY) ) +
     $                 ( COFHNO(3,ILAY,IHNO3)*TRCPRD(3,ILAY) ) +
     $                 ( COFHNO(4,ILAY,IHNO3)*TRCPRD(4,ILAY) )
                QDKHNO3=QDKHNO3*HNOJACMLT(ILAY)
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

             IF (LN2O) THEN
                QDKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) +
     $                ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) +
     $                ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) +
     $                ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) +
     $                ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) +
     $                ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) +
     $                ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                QDKN2O=QDKN2O*N2OJACMLT(ILAY)
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
             IF (LNH3) THEN
                QDKNH3=( COFNH3(1,ILAY,INH3)*TRCPRD(1,ILAY) ) +
     $                ( COFNH3(2,ILAY,INH3)*TRCPRD(2,ILAY) ) +
     $                ( COFNH3(3,ILAY,INH3)*TRCPRD(3,ILAY) ) +
     $                ( COFNH3(4,ILAY,INH3)*TRCPRD(4,ILAY) )
                QDKNH3=QDKNH3*NH3JACMLT(ILAY)
             ENDIF
C
ccc
c this block for testing
c      DKSO2=0.0
c      DKHNO3=0.0
c      DKN2O=0.0
C       DKNH3=0.0
C       DKHDO=0.0

      KHDO=0.0
ccc
C            Limit -DK so it can never totally totally cancel KFIX
c             DK = DKSO2 + DKHNO3 + DKN2O + DKNH3
c             IF (-DK .GE. KFIX) THEN
c                DK = -0.999*KFIX
c             ENDIF
             DK = 0

C            Calc total layer optical depth
                 IF (LH2O .EQ. .FALSE.) THEN
                   IF (IWHICHJAC .EQ. 1) THEN 
                     KW(ILAY) = KW_T(ILAY)
                   ELSEIF (IWHICHJAC .EQ. 2) THEN 
                     KW(ILAY) = KW_1(ILAY)
                   ELSEIF (IWHICHJAC .GE. 3) THEN 
                     KW(ILAY) = 0.0
                   END IF
                 END IF

c             f1 = 1; f2 = 1; f3 = 1; f4 = 1; f5 = 1;
c             KLAYER = KCON + KFIX + KMET + KW(ILAY) + DK

c             f1 = 1; f2 = 0; f3 = 0; f4 = 0; f5 = 0;
c             f1 = 1; f2 = 1; f3 = 1; f4 = 1; f5 = 0;
c             f1 = 0; f2 = 1; f3 = 1; f4 = 1; f5 = 0;  !!! hardly any change
c             f1 = 0; f2 = 1; f3 = 1; f4 = 0; f5 = 0;  !!! big change
c             KLAYER = f1*KCON + f2*KFIX + f3*KMET + f4*KW(ILAY) + f5*DK

             KLAYER = KCON + KFIX + KMET + KW(ILAY) + DK


                   IF (IWHICHJAC .EQ. 1) THEN 
                     DTAU_DTZ(ILAY,J)=KLAYER
                   ELSEIF (IWHICHJAC .EQ. 2) THEN 
                     DTAU_DG1(ILAY,J)=KLAYER
                   ELSEIF (IWHICHJAC .EQ. 3) THEN 
                     DTAU_DG3(ILAY,J)=KLAYER
                   ELSEIF (IWHICHJAC .EQ. 4) THEN 
                     DTAU_DG2(ILAY,J)=0
                   ELSEIF (IWHICHJAC .EQ. 5) THEN 
                     DTAU_DG4(ILAY,J)=QDKN2O
                   ELSEIF (IWHICHJAC .EQ. 6) THEN 
                     DTAU_DG5(ILAY,J)=0
                   ELSEIF (IWHICHJAC .EQ. 7) THEN 
                     DTAU_DG6(ILAY,J)=KMET
                   ELSEIF (IWHICHJAC .EQ. 8) THEN 
                     DTAU_DG9(ILAY,J)=QDKSO2
                   ELSEIF (IWHICHJAC .EQ. 9) THEN 
                     DTAU_DG12(ILAY,J)=QDKHNO3
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
C       ENDDO     ! DO I = IY,IY
C      End loops on channel number (frequency)
C
cbaba
              END IF   !!!! IF (INTERSECT(ITRYJAC,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) : is this a jac to work on?????
            END DO     !!! DO IWHICHJAC = 1,3
          END IF  !! DOJAC

c************************************************************************

       RETURN
       END
