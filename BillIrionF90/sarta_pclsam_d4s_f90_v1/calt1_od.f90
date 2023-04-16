!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS
!
!    CALT1 (for set1 = FWO)  version with trace gases
!
!F77====================================================================


!ROUTINE NAME:
!    CALT1


!ABSTRACT:
!    Calculate the transmittance for set1 using the predictors
!    and the fast transmittance coefficients.


!CALL PROTOCOL:
!    CALT1 ( INDCHN, NLAY, NCHN1, CLIST1, COEF1,
!      FIXMUL, CONPD1, FPRED1, WPRED1, OPRED1, TRCPRD,
!      INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
!      INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT,
!      INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX,
!      LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ)


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INT arr   CLIST1  set channel list            none
!    REAL arr  COEF1   set1 fast trans coefs       various
!    REAL arr  CONPD1  set1 H2O continuum preds    various
!    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
!    REAL arr  FPRED1  set1 fixed gases preds      various
!    INT arr   INDCHN  channel indices             none
!    INTEGER   NLAY    Number of layers to bottom  none
!    INTEGER   NCHN1   set1 number of channels     none
!    REAL arr  OPRED1  set1 ozone predictors       various
!    REAL arr  WPRED1  set1 water predictors       various
!    REAL arr  TRCPRD  trace gas pert predictors   various
!    INT arr   INDCO2  CO2 pert chan indices       none
!    REAL arr  COFCO2  CO2 pert coefs              various
!    REAL arr  CO2MLT  CO2 pert multiplier         none
!    INT arr   INDSO2  SO2 pert chan indices       none
!    REAL arr  COFSO2  SO2 pert coefs              various
!    REAL      SO2MLT  SO2 pert multiplier         none
!    INT arr   INDHNO  HNO3 pert chan indices      none
!    REAL arr  COFHNO  HNO3 pert coefs             various
!    REAL arr  HNOMLT  HNO3 pert multiplier        none
!    INT arr   INDN2O  N2O pert chan indices       none
!    REAL arr  COFN2O  N2O pert coefs              various
!    REAL arr  N2OMLT  N2O pert multiplier         none
!    INT arr   INDH2O  OPTRAN H2O chan indices     none
!    REAL arr  H2OPRD  OPTRAN H2O predictors       various
!    REAL arr  COFH2O  OPTRAN H2O coefs            various
!    INTEGER   LOPMAX  OPTRAN max level            none
!    INTEGER   LOPLOW  OPTRAN low bracketing level none
!    LOG arr   LOPUSE  OPTRAN level needed?        none
!    REAL arr  WAOP    OPTRAN layer water amounts  kilomoles/cm^2
!    REAL arr  DAOP    OPTRAN-to-AIRS interp fac   none
!    REAL arr  WAANG   AIRS layer water amounts    kilomoles/cm^2


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
!    The FTC coefficents and profile FTC variables are multiplied
!    together and summed to calculate the effective layer
!    transmittances. Fixed, water, and ozone transmittances are each
!    checked individually to be sure they give 0 < trans < 1.
!
!    ===================================================================
!    The routine loops over the selected channels of set1.  For each
!    channel, it first decides if needs to do a calculation for
!    variable CO2, and also if it needs to do an OPTRAN water calc (if
!    so, it does so immediately).  The program then loops downward over
!    all the layers and computes the layer transmittances TAU.
!
!    The water continuum absorption coefficient is:
!       k_con = the sum i=1 to 5 of { COEF(i)*CONPRD(i) }
!
!    The layer effective fixed gas absorption coefficient is:
!       k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
!    The layer effective water lines absorption coefficient is:
!       k_water = the sum i=1 to 11 of { COEF(5+8+i)*WPRED(i) }
!
!    The layer effective ozone absorption coefficient is:
!       k_ozone = the sum i=1 to 5 of { COEF(5+8+11+i)*OPRED(i) }
!
!    where
!      "COEF" are the fast transmittance coefficients COEF1
!      "CONPRD" are the water continuum predictors CONPRD
!      "FPRED" are the fixed gases predictors FPRED1
!      "WPRED" are the water lines predictors WPRED1
!      "OPRED" are the ozone predictors OPRED1
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
!    Date        Programmer     Comments
!    ----------- -------------- ----------------------------------------
!    Dec  1 1994 Scott Hannon   Created
!     3 Feb 1997 Scott Hannon   Re-wrote (from CALTAU) for FWO
!     3 Sep 1997 Scott Hannon   Added TAUZ and BLMULT
!    30 Sep 1997 Scott Hannon   Added variable CO2
!    27 Feb 1998 Scott Hannon   Added OPTRAN H2O
!     6 Mar 1998 Scott Hannon   Deleted water preds 12 & 13 and shifted
!                               ozone coef indices to 24-28 (was 26-30)
!     4 May 1998 Scott Hannon   Fix error: INDH2O(MXCHAN) not (MXCHNW)
!    26 Aug 1998 Scott Hannon   Add NLAY to call to CALOKW
!    11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
!    12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
!    18 May 2005 Scott Hannon   Add HNO3 based on SO2 code
!    28 Jun 2005 Scott Hannon   "trace" version for CO2,SO2,HNO3,N2O
!    28 Mar 2006 Scott Hannon   Change TAU from trans to optical depth
!    22 Dec 2006 Scott Hannon   Change TAUZ from trans to optical depth
!                               and from (1 x n) to (m x n) array;
!                               delete func QIKEXP & argument BLMULT.
!    14 Sep 2010 Scott Hannon   Add 5th CO2 coef


!END====================================================================

SUBROUTINE CALT1 (  &
	INDCHN, NLAY, NCHN1, CLIST1, COEF1, &
	FIXMUL, CONPD1, FPRED1, WPRED1, OPRED1, TRCPRD, &
	INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
	INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, &
	INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE, &
	  WAOP,   DAOP,  WAANG,    TAU,   TAUZ) 

	USE INCFTC
	IMPLICIT NONE
 

	!-----------------------------------------------------------------------
	! ARGUMENTS
	!-----------------------------------------------------------------------
	! Input
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDCHN
	INTEGER, INTENT(IN) :: NLAY
	INTEGER, INTENT(IN) :: NCHN1
	INTEGER, INTENT(IN), DIMENSION(MXCHN1) :: CLIST1
	REAL,    INTENT(IN), DIMENSION(N1COEF, MAXLAY, MXCHN1) :: COEF1
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: FIXMUL
	REAL, INTENT(IN), DIMENSION(N1CON, MAXLAY) :: CONPD1
	REAL, INTENT(IN), DIMENSION(N1FIX, MAXLAY) :: FPRED1
	REAL, INTENT(IN), DIMENSION(N1H2O, MAXLAY) :: WPRED1
	REAL, INTENT(IN), DIMENSION(N1O3, MAXLAY) :: OPRED1
	REAL, INTENT(IN), DIMENSION(NTRACE, MAXLAY) :: TRCPRD
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDCO2
	REAL, INTENT(IN), DIMENSION(NCO2, MAXLAY, MXCHNC) :: COFCO2
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: CO2MLT
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDSO2
	REAL, INTENT(IN), DIMENSION(NSO2, MAXLAY, MXCHNS) :: COFSO2
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: SO2MLT
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDHNO
	REAL, INTENT(IN), DIMENSION(NHNO3, MAXLAY, MXCHNH) :: COFHNO
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: HNOMLT
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDN2O
	REAL, INTENT(IN), DIMENSION(NN2O, MAXLAY, MXCHNN) :: COFN2O
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: N2OMLT
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDH2O
	REAL, INTENT(IN), DIMENSION(NH2O, MXOWLY) :: H2OPRD
	REAL, INTENT(IN), DIMENSION(NH2O, MXOWLY, MXCHNW) :: COFH2O
	INTEGER, INTENT(IN) :: LOPMIN
	INTEGER, INTENT(IN) :: LOPMAX
	INTEGER, INTENT(IN), DIMENSION(MAXLAY) :: LOPLOW
	LOGICAL, INTENT(IN), DIMENSION(MXOWLY) :: LOPUSE
	REAL, INTENT(IN), DIMENSION(MXOWLY) :: WAOP
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: DAOP
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: WAANG

	! Output
	REAL, INTENT(OUT), DIMENSION(MAXLAY, MXCHAN) :: TAU
	REAL, INTENT(OUT), DIMENSION(MAXLAY, MXCHAN) :: TAUZ


	! LOCAL VARIABLES
	INTEGER      I
	INTEGER   ICO2
	INTEGER  IHNO3
	INTEGER   ILAY
	INTEGER   IN2O
	INTEGER   ISO2
	INTEGER      J
	REAL     DK
	REAL  DKCO2
	REAL DKHNO3
	REAL  DKN2O
	REAL  DKSO2
	REAL   KCON
	REAL   KFIX
	REAL KLAYER
	REAL   KOZO
	REAL     KZ
	REAL   KZFW
	LOGICAL   LCO2
	LOGICAL   LH2O
	LOGICAL  LHNO3
	LOGICAL   LN2O
	LOGICAL   LSO2
	!
	! for CALOKW
	INTEGER   IH2O
	REAL     KW(MAXLAY)


	!-----------------------------------------------------------------------
	!      SAVE STATEMENTS
	!-----------------------------------------------------------------------
	!      none


	!***********************************************************************
	!***********************************************************************
	!                    EXECUTABLE CODE
	!***********************************************************************
	!***********************************************************************

	!---------------------------
	! Loop on channel (frequency)
	!---------------------------
	DO I=1,NCHN1

		! Array index of channel in TAU
		J=INDCHN( CLIST1(I) )

		!PRINT *, "CALT1_OD NCHN1, I CLIST1(I) INDCHN( CLIST1(I) )  J ", NCHN1, I, CLIST1(I), INDCHN( CLIST1(I) ),  J

		! Determine whether or not to do variable CO2 calc
		ICO2=INDCO2( CLIST1(I) )
		!PRINT *, "CALT1_OD I CLIST1(I) INDCO2(CLIST1(I)) ", I, CLIST1(I), INDCO2(CLIST1(I)) 
		IF (ICO2 .GT. 0) THEN
			LCO2=.TRUE.
		ELSE
			LCO2=.FALSE.
		ENDIF

		! Determine whether or not to do variable SO2 calc
		ISO2=INDSO2( CLIST1(I) )
		IF (ISO2 .GT. 0) THEN
			LSO2=.TRUE.
		ELSE
			LSO2=.FALSE.
		ENDIF

		! Determine whether or not to do variable HNO3 calc
		IHNO3=INDHNO( CLIST1(I) )
		IF (IHNO3 .GT. 0) THEN
			LHNO3=.TRUE.
		ELSE
			LHNO3=.FALSE.
		ENDIF

		! Determine whether or not to do variable N2O calc
		IN2O=INDN2O( CLIST1(I) )
		IF (IN2O .GT. 0) THEN
			LN2O=.TRUE.
		ELSE
			LN2O=.FALSE.
		ENDIF

		!-------------------------
		! Do OPTRAN water if needed
		!-------------------------
		IH2O=INDH2O( CLIST1(I) )
		IF (IH2O .GT. 0) THEN
			LH2O=.FALSE.
			! Calc OPTRAN water

			CALL CALOKW( NLAY, IH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE, &
				H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW )

		ELSE
			LH2O=.TRUE.
		ENDIF

		! Initialize the layer-to-space optical depth
		KZ=0.0E+0
		KZFW=0.0E+0

		!------------------------------
		! Loop on layers (top to bottom)
		!------------------------------
		DO ILAY=1,NLAY

			!---------------------------
			! Compute the water continuum
			!---------------------------
			KCON=( COEF1(1,ILAY,I)*CONPD1(1,ILAY) ) + &
				( COEF1(2,ILAY,I)*CONPD1(2,ILAY) ) + &
				( COEF1(3,ILAY,I)*CONPD1(3,ILAY) ) + &
				( COEF1(4,ILAY,I)*CONPD1(4,ILAY) ) + &
				( COEF1(5,ILAY,I)*CONPD1(5,ILAY) ) + &
				( COEF1(6,ILAY,I)*CONPD1(6,ILAY) ) + &
				( COEF1(7,ILAY,I)*CONPD1(7,ILAY) ) 

			IF (KCON .LT. 0.0E+0) THEN
				KCON=0.0E+0
			ELSEIF (KCON .GT. 1.0E+1) THEN
				KCON=1.0E+1
			ENDIF


			!-----------------------------
			! Calc the fixed gases abs coef
			!-----------------------------
			KFIX=( COEF1( 8,ILAY,I)*FPRED1(1,ILAY) ) + &
				( COEF1( 9,ILAY,I)*FPRED1(2,ILAY) ) + &
				( COEF1(10,ILAY,I)*FPRED1(3,ILAY) ) + &
				( COEF1(11,ILAY,I)*FPRED1(4,ILAY) ) + &
				( COEF1(12,ILAY,I)*FPRED1(5,ILAY) ) + &
				( COEF1(13,ILAY,I)*FPRED1(6,ILAY) ) + &
				( COEF1(14,ILAY,I)*FPRED1(7,ILAY) ) + &
				( COEF1(15,ILAY,I)*FPRED1(8,ILAY) ) 

			KFIX=KFIX*FIXMUL(ILAY)

			IF (KFIX .LT. 0.0E+0) THEN
				KFIX=0.0E+0
			ELSEIF (KFIX .GT. 1.0E+1) THEN
				KFIX=1.0E+1
			ENDIF


			!            --------------------------
			!            Compute the water abs coef
			!            --------------------------
			IF (LH2O) THEN
				! Not an OPTRAN water channel
				KW(ILAY)= &
				( COEF1(16,ILAY,I)*WPRED1( 1,ILAY) ) + &
				( COEF1(17,ILAY,I)*WPRED1( 2,ILAY) ) + &
				( COEF1(18,ILAY,I)*WPRED1( 3,ILAY) ) + &
				( COEF1(19,ILAY,I)*WPRED1( 4,ILAY) ) + &
				( COEF1(20,ILAY,I)*WPRED1( 5,ILAY) ) + &
				( COEF1(21,ILAY,I)*WPRED1( 6,ILAY) ) + &
				( COEF1(22,ILAY,I)*WPRED1( 7,ILAY) ) + &
				( COEF1(23,ILAY,I)*WPRED1( 8,ILAY) ) + &
				( COEF1(24,ILAY,I)*WPRED1( 9,ILAY) ) + &
				( COEF1(25,ILAY,I)*WPRED1(10,ILAY) ) + &
				( COEF1(26,ILAY,I)*WPRED1(11,ILAY) ) 

				IF (KW(ILAY) .LT. 0.0E+0) KW(ILAY)=0.0E+0
			ENDIF


			!--------------------------
			! Compute the ozone abs coef
			!--------------------------
			KOZO=( COEF1(27,ILAY,I)*OPRED1(1,ILAY) ) + &
				( COEF1(28,ILAY,I)*OPRED1(2,ILAY) ) + &
				( COEF1(29,ILAY,I)*OPRED1(3,ILAY) ) + &
				( COEF1(30,ILAY,I)*OPRED1(4,ILAY) ) + &
				( COEF1(31,ILAY,I)*OPRED1(5,ILAY) ) 

			IF (KOZO .LT. 0.0E+0) THEN
				KOZO=0.0E+0
			ELSEIF (KOZO .GT. 1.0E+1) THEN
				KOZO=1.0E+1
			ENDIF

			! Update KZFW
			KZFW=KZFW + KFIX + KW(ILAY)

			!----------------------------------
			! Calc the total layer transmittance
			!----------------------------------

			!cccc
			! This block is usually commented out and is only uncommented for
			! testing purposes.
			!c
			!           kcon=0.0E+0
			!           kfix=0.0E+0
			!           kw(ilay)=0.0E+0
			!           kozo=0.0E+0
			!cccc
			!

			!----------------------------
			! Calc change in total optical
			! depth due to variable CO2
			!----------------------------
			IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0) THEN
				DKCO2=( COFCO2(1,ILAY,ICO2)*TRCPRD(1,ILAY) ) + &
				( COFCO2(2,ILAY,ICO2)*TRCPRD(2,ILAY) ) + &
				( COFCO2(3,ILAY,ICO2)*TRCPRD(3,ILAY) ) + &
				( COFCO2(4,ILAY,ICO2)*TRCPRD(4,ILAY) ) + &
				( COFCO2(5,ILAY,ICO2)*TRCPRD(5,ILAY) )
				DKCO2=DKCO2*CO2MLT(ILAY)
			ELSE
				DKCO2=0.0
			ENDIF


			!----------------------------
			! Calc change in total optical
			! depth due to variable SO2
			!----------------------------
			IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0) THEN
				DKSO2=( COFSO2(1,ILAY,ISO2)*TRCPRD(1,ILAY) ) + &
				( COFSO2(2,ILAY,ISO2)*TRCPRD(2,ILAY) ) + &
				( COFSO2(3,ILAY,ISO2)*TRCPRD(3,ILAY) ) + &
				( COFSO2(4,ILAY,ISO2)*TRCPRD(4,ILAY) ) 
				DKSO2=DKSO2*SO2MLT(ILAY)
			ELSE
				DKSO2=0.0
			ENDIF


			!----------------------------
			! Calc change in total optical
			! depth due to variable HNO3
			!----------------------------
			IF (LHNO3 .AND. HNOMLT(ILAY) .NE. 0) THEN
				DKHNO3=( COFHNO(1,ILAY,IHNO3)*TRCPRD(1,ILAY) ) + &
				( COFHNO(2,ILAY,IHNO3)*TRCPRD(2,ILAY) ) + &
				( COFHNO(3,ILAY,IHNO3)*TRCPRD(3,ILAY) ) + &
				( COFHNO(4,ILAY,IHNO3)*TRCPRD(4,ILAY) )
				DKHNO3=DKHNO3*HNOMLT(ILAY)
			ELSE
				DKHNO3=0.0
			ENDIF


			!----------------------------
			! Calc change in total optical
			! depth due to variable N2O
			!----------------------------
			IF (LN2O .AND. N2OMLT(ILAY) .NE. 0) THEN
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


			!------------------------------------------
			! Calc total optical depth and transmittance
			!------------------------------------------
			! Calc total layer optical depth
			!cc
			! this block for testing
			!      DKHNO3=0.0
			!      DKSO2=0.0
			!      DKCO2=0.0
			!      DKN2O=0.0
			!cc
			! Limit -DK so it can never totally totally cancel KFIX
			DK = DKCO2 + DKSO2 + DKHNO3 + DKN2O
			IF (-DK .GE. KFIX) THEN
				DK = -0.999*KFIX
			ENDIF

			! Calc effective layer optical
			KLAYER = KCON + KFIX + KW(ILAY) + KOZO + DK
			TAU(ILAY,J)=KLAYER

			! Calc layer-to-space optical depth
			KZ=KZ + KLAYER
			TAUZ(ILAY,J)=KZ

		ENDDO
		! End loop on levels

	ENDDO
	!End loops on channel number (frequency)

	RETURN

END SUBROUTINE CALT1
