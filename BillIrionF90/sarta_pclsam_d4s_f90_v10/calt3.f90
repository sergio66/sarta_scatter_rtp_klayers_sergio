!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS
!
!    CALT3 (for set3 = FMW) version with trace gases (no CO2)
!
!F77====================================================================


!ROUTINE NAME:
!    CALT3


!ABSTRACT:
!    Calculate the transmittance for set3 using the prdictor and the
!    fast transmittance coefficients.


!CALL PROTOCOL:
!    CALT3 ( INDCHN, NLAY, NCHN3, CLIST3, COEF3,
!       FIXMUL, CONPD3, FIXED_PRED3, CH4_PRED3, H2O_PRED3, TRACEGAS_PRED,
!       INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT,
!       INDN2O, COFN2O, N2OMLT, INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX,
!       LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ)


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INT arr   INDCHN  channel indices             none
!    INTEGER   NLAY    number of layers to bottom  none
!    INTEGER   NCHN3   set3 number of channels     none
!    INT arr   CLIST3  set3 channel list           none
!    REAL arr  COEF3   set3 fast trans coefs       various
!    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
!    REAL arr  CONPD3  set3 H2O continuum preds    various
!    REAL arr  FIXED_PRED3  set3 fixed gases preds      various
!    REAL arr  CH4_PRED3  set3 methane predictors     various
!    REAL arr  H2O_PRED3  set3 water predictors       various
!    REAL arr  TRACEGAS_PRED  trace gas pert predictors   various
!    INT arr   INDSO2  SO2 pert chan indices       none
!    REAL arr  COFSO2  SO2 pert coefs              various
!    REAL arr  SO2MLT  SO2 pert multiplier         none
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
!    The fast trans coefficents and predictors are multiplied
!    together and summed to calculate the effective layer
!    transmittances. Fixed, methane, and water transmittances are each
!    checked individually to be sure they give 0 < trans < 1.
!
!    ===================================================================
!    Loops downward over all the layers for each of the NCHN3 channels
!    to compute the layer transmittances TAU.
!
!    The water continuum absorption coefficient is:
!       k_con = the sum i=1 to 5 of { COEF(i)*H2O_CONTINUUM_PRED(i) }
!
!    The layer effective fixed gas absorption coefficient is:
!       k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
!    The layer effective methane absorption coefficient is:
!       k_methane = the sum i=1 to 9 of { COEF(5+8+i)*OPRED(i) }
!
!    The layer effective water lines absorption coefficient is:
!       k_water = the sum i=1 to 11 of { COEF(5+8+9+i)*WPRED(i) }
!
!    where
!      "COEF" are the fast transmittance coefficients COEF3
!      "H2O_CONTINUUM_PRED" are the water continuum predictors H2O_CONTINUUM_PRED
!      "FPRED" are the fixed gases predictors FIXED_PRED3
!      "MPRED" are the methane predictors O3_PRED3
!      "WPRED" are the water lines predictors H2O_PRED3
!
!    The total layer effective optical depth is:
!       TAU = [ k_con + k_fixed + k_methane + k_water ]
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
!     3 Feb 1997 Scott Hannon   Re-wrote (from CALTAU) for FMW
!     3 Sep 1997 Scott Hannon   Added TAUZ and BLMULT
!     5 Mar 1998 Scott Hannon   Added OPTRAN water and deleted water
!                               preds 12 & 13
!     4 May 1998 Scott Hannon   Fix error: INDH2O(MXCHAN) not (MXCHNW)
!    26 Aug 1998 Scott Hannon   Fix mistake: loop on NLAY not MAXLAY;
!                               Add NLAY to call to CALOKW
!    11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
!    12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
!    25 Apr 2003 Scott Hannon   Add HNO3 based on SO2 code
!    28 Jun 2005 Scott Hannon   "trace" version for SO2,HNO3,N2O
!    28 Mar 2006 Scott Hannon   Change TAU from trans to optical depth
!    22 Dec 2006 Scott Hannon   Change TAUZ from trans to optical depth
!                               and from (1 x n) to (m x n) array;
!                               delete func QIKEXP & argument BLMULT.
!	  5 Jan 2017 Bill Irion     Conversion to F90 and made into a module
!	                            Added in code to save OD's etc so that it can
!                               be started in a layer greater than 1.


!END====================================================================

MODULE CALT3

	USE INCFTC

 	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: TAU_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: TAUZ_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KZFMW_SAVED

	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KCON_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KFIX_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KMET_SAVED

	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: DKSO2_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: DKHNO3_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: DKN2O_SAVED

	LOGICAL :: CALT3_ALLOCATED

CONTAINS

	SUBROUTINE CALT3__ALLOCATE(NCHN3)

		INTEGER, INTENT(IN) :: NCHN3

		IF (ALLOCATED(TAU_SAVED)) DEALLOCATE(TAU_SAVED)
		IF (ALLOCATED(TAUZ_SAVED)) DEALLOCATE(TAUZ_SAVED)
		IF (ALLOCATED(KZFMW_SAVED)) DEALLOCATE(KZFMW_SAVED)

		IF (ALLOCATED(KCON_SAVED)) 	DEALLOCATE(KCON_SAVED)
		IF (ALLOCATED(KFIX_SAVED)) 	DEALLOCATE(KFIX_SAVED)
		IF (ALLOCATED(KMET_SAVED)) 	DEALLOCATE(KMET_SAVED)

		IF (ALLOCATED(DKSO2_SAVED)) DEALLOCATE(DKSO2_SAVED)
		IF (ALLOCATED(DKHNO3_SAVED)) DEALLOCATE(DKHNO3_SAVED)
		IF (ALLOCATED(DKN2O_SAVED)) DEALLOCATE(DKN2O_SAVED)

		ALLOCATE(TAU_SAVED(MAXLAY, NCHN3))
		ALLOCATE(TAUZ_SAVED(MAXLAY, NCHN3))
		ALLOCATE(KZFMW_SAVED(MAXLAY, NCHN3))

		ALLOCATE(KCON_SAVED(MAXLAY, NCHN3))
		ALLOCATE(KFIX_SAVED(MAXLAY, NCHN3))
		ALLOCATE(KMET_SAVED(MAXLAY, NCHN3))

		ALLOCATE(DKSO2_SAVED(MAXLAY, NCHN3))
		ALLOCATE(DKHNO3_SAVED(MAXLAY, NCHN3))
		ALLOCATE(DKN2O_SAVED(MAXLAY, NCHN3))

		CALT3_ALLOCATED = .TRUE.

		PRINT *, "CALT3 allocated"

	END SUBROUTINE CALT3__ALLOCATE


	SUBROUTINE CALT3__CALCULATE_OD ( INDCHN, NLAY, NCHN3, CLIST3, COEF3, &
		FIXMUL, CONPD3, FIXED_PRED3, CH4_PRED3, H2O_PRED3, TRACEGAS_PRED, &
		INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT, &
		INDN2O, COFN2O, N2OMLT, INDH2O, H2OPRD, COFH2O, &
		LOPMIN, LOPMAX, LOPLOW, LOPUSE, WAOP, DAOP, WAANG, TAU, TAUZ, &
		START_LAYER, DO_SPECIES)

		USE SELECT_SPECIES
		IMPLICIT NONE


		!
		! ARGUMENTS
		!
		! Input
		INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDCHN(MXCHAN)
		INTEGER, INTENT(IN) :: NLAY
		INTEGER, INTENT(IN) ::  NCHN3
		INTEGER, INTENT(IN), DIMENSION(MXCHN3) :: CLIST3
		REAL, INTENT(IN), DIMENSION(N3COEF,MAXLAY,MXCHN3) :: COEF3
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: FIXMUL
		REAL, INTENT(IN), DIMENSION( N3CON,MAXLAY) :: CONPD3
		REAL, INTENT(IN), DIMENSION( N3FIX,MAXLAY) :: FIXED_PRED3
		REAL, INTENT(IN), DIMENSION( N3CH4,MAXLAY) :: CH4_PRED3
		REAL, INTENT(IN), DIMENSION ( N3H2O,MAXLAY) :: H2O_PRED3
		REAL, INTENT(IN), DIMENSION (NTRACE,MAXLAY) :: TRACEGAS_PRED
		INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDSO2
		REAL, INTENT(IN), DIMENSION(  NSO2,MAXLAY,MXCHNS) :: COFSO2
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: SO2MLT
		INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDHNO
		REAL, INTENT(IN), DIMENSION( NHNO3,MAXLAY,MXCHNH) :: COFHNO
		REAL, INTENT(IN), DIMENSION(MAXLAY) ::HNOMLT
		INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDN2O
		REAL, INTENT(IN), DIMENSION(  NN2O,MAXLAY,MXCHNN) :: COFN2O
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: N2OMLT
		INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDH2O
		REAL, INTENT(IN), DIMENSION(  NH2O,MXOWLY) :: H2OPRD
		REAL, INTENT(IN), DIMENSION(  NH2O,MXOWLY,MXCHNW) :: COFH2O
		INTEGER, INTENT(IN) :: LOPMIN
		INTEGER, INTENT(IN) :: LOPMAX
		INTEGER, INTENT(IN), DIMENSION(MAXLAY) :: LOPLOW
		LOGICAL, INTENT(IN), DIMENSION(MXOWLY) :: LOPUSE
		REAL, INTENT(IN), DIMENSION(MXOWLY) ::   WAOP
		REAL, INTENT(IN), DIMENSION(MAXLAY) ::   DAOP
		REAL, INTENT(IN), DIMENSION(MAXLAY) ::  WAANG

		INTEGER, INTENT(IN) :: START_LAYER
		TYPE (DO_SPECIES_LOGICAL_ARRAY_TYPE), INTENT(IN) :: DO_SPECIES

		! Output
		REAL    TAU(MAXLAY,MXCHAN)
		REAL   TAUZ(MAXLAY,MXCHAN)

		!
		! LOCAL VARIABLES
		!
		INTEGER      I
		INTEGER  IHNO3
		INTEGER   IN2O
		INTEGER   ILAY
		INTEGER   ISO2
		INTEGER      J
		REAL     DK
		REAL DKHNO3
		REAL  DKN2O
		REAL  DKSO2
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

		! for CALOKW
		INTEGER   IH2O
		REAL     KW(MAXLAY)

		INTEGER :: INIT_LAYER ! (= START_LAYER - 1)

		!
		!
		! EXECUTABLE CODE
		!
		!

		!
		! Loop on channel (frequency)
		!
		DO I = 1, NCHN3

			! Index for TAU
			J = INDCHN( CLIST3(I) )


			! Determine whether or not to do variable SO2
			ISO2=INDSO2( CLIST3(I) )
			IF (ISO2 .GT. 0) THEN
				LSO2=.TRUE.
			ELSE
				LSO2=.FALSE.
			ENDIF

			! Determine whether or not to do variable HNO3
			IHNO3=INDHNO( CLIST3(I) )
			IF (IHNO3 .GT. 0) THEN
				LHNO3=.TRUE.
			ELSE
				LHNO3=.FALSE.
			ENDIF

			! Determine whether or not to do variable N2O
			IN2O=INDN2O( CLIST3(I) )
			IF (IN2O .GT. 0) THEN
				LN2O=.TRUE.
			ELSE
				LN2O=.FALSE.
			ENDIF

			! 
			! Do OPTRAN water if needed
			!
			IH2O=INDH2O( CLIST3(I) )
			IF (IH2O .GT. 0) THEN
				LH2O=.FALSE.
				! Calc OPTRAN water
				CALL CALOKW( NLAY, IH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE, &
					H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW )			
			ELSE
				LH2O=.TRUE.
			ENDIF

			IF (START_LAYER .EQ. 1) THEN
				! Initialize the layer-to-space optical depth
				KZ=0.0E+0
				KZFMW=0.0E+0
			ELSE
				! Restore TAU and TAUZ above START_LAYER
				INIT_LAYER = START_LAYER - 1
				TAU(1:INIT_LAYER, J) = TAU_SAVED(1:INIT_LAYER, I)
				TAUZ(1:INIT_LAYER, J) = TAUZ_SAVED(1:INIT_LAYER, I)
				KZ = TAUZ_SAVED(INIT_LAYER, I)
				KZFMW = KZFMW_SAVED(INIT_LAYER, I)
			ENDIF

			!
			! Loop on layers (top to ground)
			!
			DO ILAY = START_LAYER, NLAY

				IF (DO_SPECIES % DO_H2O) THEN
					!
					! Compute the water continuum
					!
					KCON = &
						( COEF3(1,ILAY,I)*CONPD3(1,ILAY) ) + &
						( COEF3(2,ILAY,I)*CONPD3(2,ILAY) ) + &
						( COEF3(3,ILAY,I)*CONPD3(3,ILAY) ) + &
						( COEF3(4,ILAY,I)*CONPD3(4,ILAY) ) + &
						( COEF3(5,ILAY,I)*CONPD3(5,ILAY) ) + &
						( COEF3(6,ILAY,I)*CONPD3(6,ILAY) ) + &
						( COEF3(7,ILAY,I)*CONPD3(7,ILAY) )

					IF (KCON .LT. 0.0E+0) THEN
						KCON=0.0E+0
					ELSEIF (KCON .GT. 1.0E+1) THEN
						KCON=1.0E+1
					ENDIF
				ELSE
					KCON = KCON_SAVED(ILAY, I)
				ENDIF

				IF (DO_SPECIES % DO_TEMPERATURE) THEN
					!
					! Calc the fixed gases abs coef
					!
					KFIX= &
						( COEF3( 8,ILAY,I)*FIXED_PRED3(1,ILAY) ) + &
						( COEF3( 9,ILAY,I)*FIXED_PRED3(2,ILAY) ) + &
						( COEF3(10,ILAY,I)*FIXED_PRED3(3,ILAY) ) + &
						( COEF3(11,ILAY,I)*FIXED_PRED3(4,ILAY) ) + &
						( COEF3(12,ILAY,I)*FIXED_PRED3(5,ILAY) ) + &
						( COEF3(13,ILAY,I)*FIXED_PRED3(6,ILAY) ) + &
						( COEF3(14,ILAY,I)*FIXED_PRED3(7,ILAY) ) + &
						( COEF3(15,ILAY,I)*FIXED_PRED3(8,ILAY) )

					KFIX=KFIX*FIXMUL(ILAY)

					IF (KFIX .LT. 0.0E+0) THEN
						KFIX=0.0E+0
					ELSEIF (KFIX .GT. 1.0E+1) THEN
						KFIX=1.0E+1
					ENDIF
				ELSE
					KFIX = KFIX_SAVED(ILAY, I)
				ENDIF

				IF (DO_SPECIES % DO_CH4) THEN
					!
					! Compute the methane abs coef
					!
					KMET = &
						( COEF3(16,ILAY,I)*CH4_PRED3(1,ILAY) ) + &
						( COEF3(17,ILAY,I)*CH4_PRED3(2,ILAY) ) + &
						( COEF3(18,ILAY,I)*CH4_PRED3(3,ILAY) ) + &
						( COEF3(19,ILAY,I)*CH4_PRED3(4,ILAY) ) + &
						( COEF3(20,ILAY,I)*CH4_PRED3(5,ILAY) ) + &
						( COEF3(21,ILAY,I)*CH4_PRED3(6,ILAY) ) + &
						( COEF3(22,ILAY,I)*CH4_PRED3(7,ILAY) ) + &
						( COEF3(23,ILAY,I)*CH4_PRED3(8,ILAY) ) + &
						( COEF3(24,ILAY,I)*CH4_PRED3(9,ILAY) )

					IF (KMET .LT. 0.0E+0) THEN
						KMET=0.0E+0
					ELSEIF (KMET .GT. 1.0E+1) THEN
						KMET=1.0E+1
					ENDIF
				ELSE
					KMET = KMET_SAVED(ILAY, I)
				ENDIF
		
				!
				! Compute the water abs coef
				!
				IF (LH2O) THEN
					! Not an OPTRAN water channel
					KW(ILAY) = &
						( COEF3(25,ILAY,I)*H2O_PRED3( 1,ILAY) ) + &
						( COEF3(26,ILAY,I)*H2O_PRED3( 2,ILAY) ) + &
						( COEF3(27,ILAY,I)*H2O_PRED3( 3,ILAY) ) + &
						( COEF3(28,ILAY,I)*H2O_PRED3( 4,ILAY) ) + &
						( COEF3(29,ILAY,I)*H2O_PRED3( 5,ILAY) ) + &
						( COEF3(30,ILAY,I)*H2O_PRED3( 6,ILAY) ) + &
						( COEF3(31,ILAY,I)*H2O_PRED3( 7,ILAY) ) + &
						( COEF3(32,ILAY,I)*H2O_PRED3( 8,ILAY) ) + &
						( COEF3(33,ILAY,I)*H2O_PRED3( 9,ILAY) ) + &
						( COEF3(34,ILAY,I)*H2O_PRED3(10,ILAY) ) + &
						( COEF3(35,ILAY,I)*H2O_PRED3(11,ILAY) )
	
					IF (KW(ILAY) .LT. 0.0E+0) KW(ILAY)=0.0E+0
				ENDIF
				!
				!            Update KZFMW
				KZFMW=KZFMW + KFIX + KMET + KW(ILAY)

				IF (DO_SPECIES % DO_ALL .AND. START_LAYER .EQ. 1) KZFMW_SAVED(ILAY, I) = KZFMW

				!
				! Calc the total layer transmittance
				!
				!
				!cccc
				! This block is usually commented out and is only uncommented for
				! testing purposes.
				!
				!           kcon=0.0
				!           kfix=0.0
				!           kmet=0.0
				!           kw(ilay)=0.0
				!cccc

				IF (DO_SPECIES % DO_SO2) THEN
					!
					! Calc change in total optical
					! depth due to variable SO2
					!
					IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0) THEN
						DKSO2 = &
							( COFSO2(1,ILAY,ISO2)*TRACEGAS_PRED(1,ILAY) ) + &
							( COFSO2(2,ILAY,ISO2)*TRACEGAS_PRED(2,ILAY) ) + &
							( COFSO2(3,ILAY,ISO2)*TRACEGAS_PRED(3,ILAY) ) + &
							( COFSO2(4,ILAY,ISO2)*TRACEGAS_PRED(4,ILAY) )
						DKSO2=DKSO2*SO2MLT(ILAY)
					ELSE
						DKSO2=0.0
					ENDIF
				ELSE
					DKSO2 = DKSO2_SAVED(ILAY, I)
				ENDIF

				IF (DO_SPECIES % DO_HNO3) THEN
					!
					! Calc change in total optical
					! depth due to variable HNO3
					!
					IF (LHNO3 .AND. HNOMLT(ILAY) .NE. 0) THEN
						DKHNO3 = &
							( COFHNO(1,ILAY,IHNO3)*TRACEGAS_PRED(1,ILAY) ) + &
							( COFHNO(2,ILAY,IHNO3)*TRACEGAS_PRED(2,ILAY) ) + &
							( COFHNO(3,ILAY,IHNO3)*TRACEGAS_PRED(3,ILAY) ) + &
							( COFHNO(4,ILAY,IHNO3)*TRACEGAS_PRED(4,ILAY) )
						DKHNO3 = DKHNO3*HNOMLT(ILAY)
					ELSE
						DKHNO3=0.0
					ENDIF
				ELSE
					DKHNO3 = DKHNO3_SAVED(ILAY, I)
				ENDIF
		
				IF (DO_SPECIES % DO_N2O) THEN
					!
					! Calc change in total optical
					! depth due to variable N2O
					!
					IF (LN2O .AND. N2OMLT(ILAY) .NE. 0) THEN
						DKN2O = &
							( COFN2O(1,ILAY,IN2O)*TRACEGAS_PRED(1,ILAY) ) + &
							( COFN2O(2,ILAY,IN2O)*TRACEGAS_PRED(2,ILAY) ) + &
							( COFN2O(3,ILAY,IN2O)*TRACEGAS_PRED(3,ILAY) ) + &
							( COFN2O(4,ILAY,IN2O)*TRACEGAS_PRED(4,ILAY) ) + &
							( COFN2O(5,ILAY,IN2O)*TRACEGAS_PRED(5,ILAY) ) + &
							( COFN2O(6,ILAY,IN2O)*TRACEGAS_PRED(6,ILAY) ) + &
							( COFN2O(7,ILAY,IN2O)*TRACEGAS_PRED(7,ILAY) )
						DKN2O=DKN2O*N2OMLT(ILAY)
					ELSE
						DKN2O = 0.0
					ENDIF
				ELSE
					DKN2O = DKN2O_SAVED(ILAY, I)
				ENDIF

				!cc
				! this block for testing
				!      DKSO2=0.0
				!      DKHNO3=0.0
				!      DKN2O=0.0
				!cc

				! Limit -DK so it can never totally totally cancel KFIX
				DK = DKSO2 + DKHNO3 + DKN2O
				IF (-DK .GE. KFIX) THEN
					DK = -0.999*KFIX
				ENDIF

				! Calc total layer optical depth
				KLAYER = KCON + KFIX + KMET + KW(ILAY) + DK
				TAU(ILAY,J)=KLAYER

				! Calc layer-to-space optical depth
				KZ=KZ + KLAYER
				TAUZ(ILAY,J)=KZ

				! Save KCON if DO_SPECIES % DO_ALL  and the loop started from layer 1
				IF (DO_SPECIES % DO_ALL .AND. START_LAYER .EQ. 1) THEN
					KCON_SAVED(ILAY, I) = KCON
					KFIX_SAVED(ILAY, I) = KFIX
					KMET_SAVED(ILAY, I) = KMET

					DKSO2_SAVED(ILAY, I) = DKSO2
					DKHNO3_SAVED(ILAY, I) = DKHNO3
					DKN2O_SAVED(ILAY, I) = DKN2O
				ENDIF

			ENDDO ! End loop on levels

			! Save TAU and TAUZ if loop started from the layer 1
			IF (DO_SPECIES % DO_ALL .AND. START_LAYER .EQ. 1) THEN
				TAU_SAVED(:, I) = TAU(:, J)
				TAUZ_SAVED(:, I) = TAUZ(:, J)
			ENDIF

		ENDDO ! End loops on channel number (frequency)

		RETURN

	END SUBROUTINE CALT3__CALCULATE_OD

END MODULE CALT3
