!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS
!
!    CALT6 (set6=FWO sun mfmw) version for trace gases (no HNO3)
!
!F77====================================================================


!ROUTINE NAME:
!    CALT6


!ABSTRACT:
!    Calculate the transmittance for set6 using the predictors and the
!    fast transmittance coefficients.


!CALL PROTOCOL:
!    CALT6( INDCHN, NLAY, NCHN6, CLIST6, COEF6,
!       FIXMUL, CONPD6, FPRED6, WPRED6, OPRED6, TRCPRD,
!       INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
!       INDN2O, COFN2O, N2OMLT, TAU, TAUZ )


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INT arr   INDCHN  channel indices             none
!    INTEGER   NLAY    number of layers to bottom  none
!    INTEGER   NCHN6   set6 number of channels     none
!    INT arr   CLIST6  set6 channel list           none
!    REAL arr  COEF6   set6 fast trans coefs       various
!    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
!    REAL arr  CONPD6  set6 H2O continuum preds    various
!    REAL arr  FPRED6  set6 fixed gases preds      various
!    REAL arr  WPRED6  set6 water predictors       various
!    REAL arr  OPRED6  set6 ozone predictors       various
!    REAL arr  TRCPRD  trace gas pert predictors   various
!    INT arr   INDCO2  CO2 pert chan indices       none
!    REAL arr  COFCO2  CO2 pert coefs              various
!    REAL arr  CO2MLT  CO2 pert multiplier         none
!    REAL arr  SO2PRD  SO2 pert predictors         various
!    INT arr   INDSO2  SO2 pert chan indices       none
!    REAL arr  COFSO2  SO2 pert coefs              various
!    REAL arr  SO2MLT  SO2 pert multiplier         none
!    INT arr   INDN2O  N2O pert chan indices       none
!    REAL arr  COFN2O  N2O pert coefs              various
!    REAL arr  N2OMLT  N2O pert multiplier         none


!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL arr  TAU     effective layer opt depth   none
!    REAL arr  TAUZ    layer-to-space opt depth    none


!INPUT/OUTPUT PARAMETERS:
!    none


!RETURN VALUES:
!    none


!PARENT(S):
!    USEFAST


!ROUTINES CALLED:
!    none


!FILES ACCESSED:
!    incFTC.f : include file of parameter statements accessed during
!       compilation only.


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    August 2000 version of the 100 layer AIRS Fast Transmittance
!    Code by L.L.Strow/S.Hannon.
!
!    The fast trans coefficents and predictors are multiplied
!    together and summed to calculate the effective layer
!    transmittances. Fixed, water, and ozone transmittances are each
!    checked individually to be sure they give 0 < trans < 1.
!
!    ===================================================================
!    Loops downward over all the layers for each of the NCHN6 channels
!    to compute the layer transmittances TAU.
!
!    The water continuum absorption coefficient is:
!       k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
!
!    The layer effective fixed gas absorption coefficient is:
!       k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
!    The layer effective water lines absorption coefficient is:
!       k_water = the sum i=1 to 7 of { COEF(5+8+i)*WPRED(i) }
!
!    The layer effective ozone absorption coefficient is:
!       k_ozone = COEF(5+8+7+1)*OPRED(1)
!
!    where
!      "COEF" are the fast transmittance coefficients COEF5
!      "CONPRD" are the water continuum predictors CONPRD
!      "FPRED" are the fixed gases predictors FPRED5
!      "WPRED" are the water lines predictors WPRED5
!      "OPRED" are the ozone predictors OPRED5
!
!    The total layer effective optical depth TAU is:
!       TAU = [ k_con + k_fixed + k_water + k_ozone ]
!
!    ===================================================================


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
! Date        Programmer     Comments
! ----------- -------------- -------------------------------------------
! 03 Jul 1997 Scott Hannon   Created for set6
! 03 Sep 1997 Scott Hannon   Added TAUZ and BLMULT
! 30 Sep 1997 Scott Hannon   Added variable CO2
! 27 Feb 1998 Scott Hannon   Added LTAU
! 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
! 12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
! 03 Jan 2003 Scott Hannon   Add XZ
! 06 Feb 2003 Scott Hannon   Bug fix - ozone use coef 23, not 21
! 25 Apr 2003 Scott Hannon   Add SO2
! 28 Jun 2005 Scott Hannon   "trace" version with CO2,SO2,N2O
! 28 Mar 2006 Scott Hannon   Change TAU from trans to optical depth
! 22 Dec 2006 Scott Hannon   Change TAUZ from trans to optical depth
!                            and from (1 x n) to (m x n) array;
!                            delete func QIKEXP & arguments BLMULT
!                            & LTAU & XZ.
! 02 Sep 2008 Scott Hannon   Add 5th CO2 predictor
! 05 Sep 2017 Bill Irion     Conversion from F77 to F90
!END====================================================================

!      =================================================================
       SUBROUTINE CALT6 ( INDCHN, NLAY, NCHN6, CLIST6, &
         COEF6, FIXMUL, CONPD6, FPRED6, WPRED6, OPRED6, TRCPRD, &
         INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
         INDN2O, COFN2O, N2OMLT, TAU, TAUZ)
!      =================================================================

!-----------------------------------------------------------------------
!      IMPLICIT NONE
!-----------------------------------------------------------------------
       USE INCFTC
       IMPLICIT NONE

!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
!      Input
       INTEGER INDCHN(MXCHAN)
       INTEGER   NLAY
       INTEGER  NCHN6
       INTEGER CLIST6(MXCHN6)
       REAL  COEF6(N6COEF,MAXLAY,MXCHN6)
       REAL FIXMUL(MAXLAY)
       REAL CONPD6( N6CON,MAXLAY)
       REAL FPRED6( N6FIX,MAXLAY)
       REAL WPRED6( N6H2O,MAXLAY)
       REAL OPRED6(  N6O3,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)
       INTEGER INDCO2(MXCHAN)
       REAL COFCO2(  NCO2,MAXLAY,MXCHNC)
       REAL CO2MLT(MAXLAY)
       INTEGER INDSO2(MXCHAN)
       REAL COFSO2(  NSO2,MAXLAY,MXCHNS)
       REAL SO2MLT(MAXLAY)
       INTEGER INDN2O(MXCHAN)
       REAL COFN2O(  NN2O,MAXLAY,MXCHNN)
       REAL N2OMLT(MAXLAY)
!
!      Output
       REAL    TAU(MAXLAY,MXCHAN)
       REAL   TAUZ(MAXLAY,MXCHAN)


!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER      I
       INTEGER   ICO2
       INTEGER   ILAY
       INTEGER   IN2O
       INTEGER   ISO2
       INTEGER      J
       REAL     DK
       REAL  DKCO2
       REAL  DKN2O
       REAL  DKSO2
       REAL   KCON
       REAL   KFIX
       REAL KLAYER
       REAL   KOZO
       REAL   KWAT
       REAL     KZ
       LOGICAL   LCO2
       LOGICAL   LN2O
       LOGICAL   LSO2


!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none


!***********************************************************************
!***********************************************************************
!                    EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
!      ---------------------------
!      Loop on channel (frequency)
!      ---------------------------
       DO I=1,NCHN6
!
!         Index for TAU
          J=INDCHN( CLIST6(I) )
!
!         Determine whether or not to do variable CO2
          ICO2=INDCO2( CLIST6(I) )
          IF (ICO2 .GT. 0) THEN
             LCO2=.TRUE.
          ELSE
             LCO2=.FALSE.
          ENDIF
!
!         Determine whether or not to do variable SO2
          ISO2=INDSO2( CLIST6(I) )
          IF (ISO2 .GT. 0) THEN
             LSO2=.TRUE.
          ELSE
             LSO2=.FALSE.
          ENDIF
!
!         Determine whether or not to do variable N2O
          IN2O=INDN2O( CLIST6(I) )
          IF (IN2O .GT. 0) THEN
             LN2O=.TRUE.
          ELSE
             LN2O=.FALSE.
          ENDIF
!
!         Initialize the layer-to-space optical depth
          KZ=0.0E+0
!
!         ------------------------------
!         Loop on layers (top to ground)
!         ------------------------------
          DO ILAY=1,NLAY
!
!            ---------------------------
!            Compute the water continuum
!            ---------------------------
             KCON=( COEF6(1,ILAY,I)*CONPD6(1,ILAY) ) + &
                 ( COEF6(2,ILAY,I)*CONPD6(2,ILAY) ) + &
                 ( COEF6(3,ILAY,I)*CONPD6(3,ILAY) ) + &
                 ( COEF6(4,ILAY,I)*CONPD6(4,ILAY) ) + &
                 ( COEF6(5,ILAY,I)*CONPD6(5,ILAY) ) + &
                 ( COEF6(6,ILAY,I)*CONPD6(6,ILAY) ) + &
                 ( COEF6(7,ILAY,I)*CONPD6(7,ILAY) )
!
             IF (KCON .LT. 0.0E+0) THEN
                KCON=0.0E+0
             ELSEIF (KCON .GT. 1.0E+1) THEN
                KCON=1.0E+1
             ENDIF
!
!            -----------------------------
!            Calc the fixed gases abs coef
!            -----------------------------
             KFIX=( COEF6( 8,ILAY,I)*FPRED6( 1,ILAY) ) + &
                 ( COEF6( 9,ILAY,I)*FPRED6( 2,ILAY) ) + &
                 ( COEF6(10,ILAY,I)*FPRED6( 3,ILAY) ) + &
                 ( COEF6(11,ILAY,I)*FPRED6( 4,ILAY) ) + &
                 ( COEF6(12,ILAY,I)*FPRED6( 5,ILAY) ) + &
                 ( COEF6(13,ILAY,I)*FPRED6( 6,ILAY) ) + &
                 ( COEF6(14,ILAY,I)*FPRED6( 7,ILAY) ) + &
                 ( COEF6(15,ILAY,I)*FPRED6( 8,ILAY) )
!
             KFIX=KFIX*FIXMUL(ILAY)
!
             IF (KFIX .LT. 0.0E+0) THEN
                KFIX=0.0E+0
             ELSEIF (KFIX .GT. 1.0E+1) THEN
                KFIX=1.0E+1
             ENDIF
!
!
!            --------------------------
!            Compute the water abs coef
!            --------------------------
             KWAT=( COEF6(16,ILAY,I)*WPRED6( 1,ILAY) ) + &
                 ( COEF6(17,ILAY,I)*WPRED6( 2,ILAY) ) + &
                 ( COEF6(18,ILAY,I)*WPRED6( 3,ILAY) ) + &
                 ( COEF6(19,ILAY,I)*WPRED6( 4,ILAY) ) + &
                 ( COEF6(20,ILAY,I)*WPRED6( 5,ILAY) ) + &
                 ( COEF6(21,ILAY,I)*WPRED6( 6,ILAY) ) + &
                 ( COEF6(22,ILAY,I)*WPRED6( 7,ILAY) )
!
             IF (KWAT .LT. 0.0E+0) THEN
                KWAT=0.0E+0
             ELSEIF( KWAT .GT. 1.0E+1) THEN
                KWAT=1.0E+1
             ENDIF
!
!
!            --------------------------
!            Compute the ozone abs coef
!            --------------------------
             KOZO=( COEF6(23,ILAY,I)*OPRED6(1,ILAY) )
!
             IF (KOZO .LT. 0.0E+0) THEN
                KOZO=0.0E+0
             ELSEIF (KOZO .GT. 1.0E+1) THEN
                KOZO=1.0E+1
             ENDIF
!
!
!            ----------------------------------
!            Calc the total layer transmittance
!            ----------------------------------
!
!cccc
! This block is usually commented out and is only uncommented for
! testing purposes.
!
!           kcon=0.0E+0
!           kfix=0.0E+0
!           kwat=0.0E+0
!           kozo=0.0E+0
!cccc
!
!            ----------------------------
!            Calc change in total optical
!            depth due to variable CO2
!            ----------------------------
             IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0.0) THEN
                DKCO2=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) + &
                     ( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) + &
                     ( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) + &
                     ( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) + &
                     ( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
                DKCO2=DKCO2*CO2MLT(ILAY)
             ELSE
                DKCO2=0.0
             ENDIF
!
!            ----------------------------
!            Calc change in total optical
!            depth due to variable SO2
!            ----------------------------
             IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0.0) THEN
                DKSO2=( COFSO2(1,ILAY,ISO2)*TRCPRD(1,ILAY) ) + &
                     ( COFSO2(2,ILAY,ISO2)*TRCPRD(2,ILAY) ) + &
                     ( COFSO2(3,ILAY,ISO2)*TRCPRD(3,ILAY) ) + &
                     ( COFSO2(4,ILAY,ISO2)*TRCPRD(4,ILAY) )
                DKSO2=DKSO2*SO2MLT(ILAY)
             ELSE
                DKSO2=0.0
             ENDIF
!
!            ----------------------------
!            Calc change in total optical
!            depth due to variable N2O
!            ----------------------------
             IF (LN2O .AND. N2OMLT(ILAY) .NE. 0.0) THEN
                DKN2O=( COFN2O(1,ILAY,IN2O)*TRCPRD(1,ILAY) ) + &
                     ( COFN2O(2,ILAY,IN2O)*TRCPRD(2,ILAY) ) + &
                     ( COFN2O(3,ILAY,IN2O)*TRCPRD(3,ILAY) ) + &
                     ( COFN2O(4,ILAY,IN2O)*TRCPRD(4,ILAY) ) + &
                     ( COFN2O(5,ILAY,IN2O)*TRCPRD(5,ILAY) ) + &
                     ( COFN2O(6,ILAY,IN2O)*TRCPRD(6,ILAY) ) + &
                     ( COFN2O(7,ILAY,IN2O)*TRCPRD(7,ILAY) )
                DKN2O=DKN2O*N2OMLT(ILAY)
             ELSE
                DKN2O=0.0
             ENDIF
!
!cc
! this block for testing
!      DKCO2=0.0
!      DKSO2=0.0
!      DKN2O=0.0
!cc
!            Limit -DK so it can never totally totally cancel KFIX
             DK = DKCO2 + DKSO2 + DKN2O
             IF (-DK .GE. KFIX) THEN
                DK = -0.999*KFIX
             ENDIF

!            Calc total layer optical depth
             KLAYER = KCON + KFIX + KWAT + KOZO + DK
             TAU(ILAY,J)=KLAYER
!
!            Calc layer-to-space optical depth
             KZ=KZ + KLAYER
             TAUZ(ILAY,J)=KZ
!
          ENDDO
!         End loop on levels
!
       ENDDO
!      End loops on channel number (frequency)
!
       RETURN
       END
