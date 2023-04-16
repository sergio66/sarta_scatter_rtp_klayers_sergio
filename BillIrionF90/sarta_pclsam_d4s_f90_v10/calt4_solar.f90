!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS
!
!    CALT4_SOLAR (for set4 = FCOW) version with trace gases (no SO2 or HNO3)
!
!F77====================================================================


!ROUTINE NAME:
!    CALT4_SOLAR


!ABSTRACT:
!    Calculate the transmittance for set4 using the predictors and the
!    fast transmittance coefficients.


!CALL PROTOCOL:
!    CALT4_SOLAR( INDCHN, NLAY, NCHN4, CLIST4, COEF4,
!       FIXMUL, CONPD4, FIXED_PRED4, CO_PRED4, O3_PRED4, H2O_PRED4, TRACEGAS_PRED, INDCO2,
!       COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, TAU, TAUZ )


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INT arr   INDCHN  channel indices             none
!    INTEGER   NLAY    number of layers to bottom  none
!    INTEGER   NCHN4   set4 number of channels     none
!    INT arr   CLIST4  set4 channel list           none
!    REAL arr  COEF4   set4 fast trans coefs       various
!    REAL arr  FIXMUL  fixed amount mult (~1.0)    none
!    REAL arr  CONPD4  set4 H2O continuum preds    various
!    REAL arr  FIXED_PRED4  set4 fixed gases preds      various
!    REAL arr  CO_PRED4  set4 carbon monoxide preds  various
!    REAL arr  O3_PRED4  set4 ozone predictors       various
!    REAL arr  H2O_PRED4  set4 water predictors       various
!    REAL arr  TRACEGAS_PRED  trace gases pert predictors various
!    INT arr   INDCO2  CO2 pert chan indices       none
!    REAL arr  COFCO2  CO2 pert coefs              various
!    REAL arr  CO2MLT  CO2 pert multiplier         none
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
!    transmittances. Fixed, CO, ozone, and water transmittances are
!    each checked individually to be sure they give 0 < trans < 1.
!
!    ===================================================================
!    Loops downward over all the layers for each of the NCHN4 channels
!    to compute the layer transmittances TAU.
!
!    The water continuum absorption coefficient is:
!       k_con = the sum i=1 to 5 of { COEF(i)*H2O_CONTINUUM_PRED(i) }
!
!    The layer effective fixed gas absorption coefficient is:
!       k_fixed = the sum i=1 to 8 of { COEF(5+i)*FPRED(i) }
!
!    The layer effective CO absorption coefficient is:
!       k_co = the sum i=1 to 9 of { COEF(5+8+i)*OPRED(i) }
!
!    The layer effective ozone absorption coefficient is:
!       k_ozone = the sum i=1 to 3 of { COEF(5+8+9+i)*OPRED(i) }
!
!    The layer effective water lines absorption coefficient is:
!       k_water = the sum i=1 to 11 of { COEF(5+8+9+3+i)*WPRED(i) }
!
!    where
!      "COEF" are the fast transmittance coefficients COEF4
!      "H2O_CONTINUUM_PRED" are the water continuum predictors H2O_CONTINUUM_PRED
!      "FPRED" are the fixed gases predictors FIXED_PRED4
!      "CPRED" are the carbon monoxide predictors CO_PRED4
!      "OPRED" are the ozone predictors O3_PRED4
!      "WPRED" are the water lines predictors H2O_PRED4
!
!    The total layer effective optical depth is:
!       TAU = [ k_con + k_fixed + k_co + k_ozone + k_water ]
!
!    ===================================================================


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
! Date        Programmer     Comments
! ----------- -------------- -------------------------------------------
! 01 Dec 1994 Scott Hannon   Created
! 03 Feb 1997 Scott Hannon   Re-wrote (from CALTAU) for FCOW
! 03 Sep 1997 Scott Hannon   Re-wrote for sun and added TAUZ & BLMULT
! 30 Sep 1997 Scott Hannon   Added variable CO2
! 27 Feb 1998 Scott Hannon   Added LTAU
! 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
! 12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
! 03 Jan 2003 Scott Hannon   Add XZ
! 12 Oct 2004 Scott Hannon   Change CO2MLT from scaler to vector
! 28 Jun 2005 Scott Hannon   "trace" version for CO2,N2O
! 28 Mar 2006 Scott Hannon   Change TAU from trans to optical depth
! 22 Dec 2006 Scott Hannon   Change TAUZ from trans to optical depth
!                               and from (1 x n) to (m x n) array;
!                               delete func QIKEXP & arguments BLMULT
!                               & LTAU & XZ.
! 02 Sep 2008 Scott Hannon   Add 5th CO2 predictor
!  5 Jan 2017 Bill Irion     Conversion to F90 and made into a module
!	                            Added in code to save OD's etc so that it can
!                               be started in a layer greater than 1.

!END====================================================================

MODULE CALT4_SOLAR

	USE INCFTC

 	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: TAU_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: TAUZ_SAVED


	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KCON_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KFIX_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KCO_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KOZO_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KWAT_SAVED

	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: DKCO2_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: DKN2O_SAVED


	LOGICAL :: CALT4_SOLAR_ALLOCATED

CONTAINS

	SUBROUTINE CALT4_SOLAR__ALLOCATE(NCHN4)

		INTEGER, INTENT(IN) :: NCHN4

		IF (ALLOCATED(TAU_SAVED)) DEALLOCATE(TAU_SAVED)
		IF (ALLOCATED(TAUZ_SAVED)) DEALLOCATE(TAUZ_SAVED)

		IF (ALLOCATED(KCON_SAVED)) 	DEALLOCATE(KCON_SAVED)
		IF (ALLOCATED(KFIX_SAVED)) 	DEALLOCATE(KFIX_SAVED)
		IF (ALLOCATED(KCO_SAVED)) 	DEALLOCATE(KCO_SAVED)
		IF (ALLOCATED(KOZO_SAVED)) 	DEALLOCATE(KOZO_SAVED)
		IF (ALLOCATED(KWAT_SAVED)) 	DEALLOCATE(KWAT_SAVED)

		IF (ALLOCATED(DKCO2_SAVED)) DEALLOCATE(DKCO2_SAVED)
		IF (ALLOCATED(DKN2O_SAVED)) DEALLOCATE(DKN2O_SAVED)

		ALLOCATE(TAU_SAVED(MAXLAY, NCHN4))
		ALLOCATE(TAUZ_SAVED(MAXLAY, NCHN4))

		ALLOCATE(KCON_SAVED(MAXLAY, NCHN4))
		ALLOCATE(KFIX_SAVED(MAXLAY, NCHN4))
		ALLOCATE(KCO_SAVED(MAXLAY, NCHN4))
		ALLOCATE(KOZO_SAVED(MAXLAY, NCHN4))
		ALLOCATE(KWAT_SAVED(MAXLAY, NCHN4))

		ALLOCATE(DKCO2_SAVED(MAXLAY, NCHN4))
		ALLOCATE(DKN2O_SAVED(MAXLAY, NCHN4))

		CALT4_SOLAR_ALLOCATED = .TRUE.

		PRINT *, "CALT4_SOLAR allocated"	


	END SUBROUTINE CALT4_SOLAR__ALLOCATE


	SUBROUTINE CALT4_SOLAR__CALCULATE_OD ( INDCHN, NLAY, NCHN4, CLIST4, &
		COEF4, FIXMUL, CONPD4, FIXED_PRED4, CO_PRED4, O3_PRED4, H2O_PRED4, TRACEGAS_PRED, &
		INDCO2, COFCO2, CO2MLT,INDN2O, COFN2O, N2OMLT, TAU, TAUZ, START_LAYER, DO_SPECIES)

		USE SELECT_SPECIES
		IMPLICIT NONE

		!
		! ARGUMENTS
		!
		! Input
		INTEGER INDCHN(MXCHAN)
		INTEGER   NLAY
		INTEGER  NCHN4
		INTEGER CLIST4(MXCHN4)
		REAL  COEF4(N4COEF,MAXLAY,MXCHN4)
		REAL FIXMUL(MAXLAY)
		REAL CONPD4( N4CON,MAXLAY)
		REAL FIXED_PRED4( N4FIX,MAXLAY)
		REAL CO_PRED4(  N4CO,MAXLAY)
		REAL O3_PRED4(  N4O3,MAXLAY)
		REAL H2O_PRED4( N4H2O,MAXLAY)
		REAL TRACEGAS_PRED(NTRACE,MAXLAY)
		INTEGER INDCO2(MXCHAN)
		REAL COFCO2(  NCO2,MAXLAY,MXCHNC)
		REAL CO2MLT(MAXLAY)
		INTEGER INDN2O(MXCHAN)
		REAL COFN2O(  NN2O,MAXLAY,MXCHNN)
		REAL N2OMLT(MAXLAY)

		INTEGER, INTENT(IN) :: START_LAYER
		TYPE (DO_SPECIES_LOGICAL_ARRAY_TYPE), INTENT(IN) :: DO_SPECIES

		! Output
		REAL    TAU(MAXLAY,MXCHAN)
		REAL   TAUZ(MAXLAY,MXCHAN)


!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
		INTEGER  ::  I
		INTEGER  :: ICO2
		INTEGER  :: ILAY
		INTEGER  :: IN2O
		INTEGER  ::  J
		REAL     DK
		REAL  DKCO2
		REAL  DKN2O
		REAL    KCO
		REAL   KCON
		REAL   KFIX
		REAL KLAYER
		REAL   KOZO
		REAL   KWAT
		REAL     KZ
		LOGICAL   LCO2
		LOGICAL   LN2O

		INTEGER :: INIT_LAYER ! (= START_LAYER - 1)


		!
		!
		! EXECUTABLE CODE
		!
		!

		!---------------------------
		! Loop on channel (frequency)
		!---------------------------
       DO I = 1 , NCHN4

			! Index for TAU
			J=INDCHN( CLIST4(I) )

			! Determine whether or not to do variable CO2
			ICO2=INDCO2( CLIST4(I) )
			IF (ICO2 .GT. 0) THEN
				LCO2=.TRUE.
			ELSE
				LCO2=.FALSE.
			ENDIF

			! Determine whether or not to do variable N2O
			IN2O=INDN2O( CLIST4(I) )
			IF (IN2O .GT. 0) THEN
				LN2O=.TRUE.
			ELSE
				LN2O=.FALSE.
			ENDIF

			IF (START_LAYER .EQ. 1) THEN
				! Initialize the layer-to-space optical depth
				KZ=0.0E+0
			ELSE
				! Restore TAU, TAUZ and KZ above START_LAYER
				INIT_LAYER = START_LAYER - 1
				TAU(1:INIT_LAYER, J) = TAU_SAVED(1:INIT_LAYER, I)
				TAUZ(1:INIT_LAYER, J) = TAUZ_SAVED(1:INIT_LAYER, I)
				KZ = TAUZ_SAVED(INIT_LAYER, I)
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
						( COEF4(1,ILAY,I)*CONPD4(1,ILAY) ) + &
						( COEF4(2,ILAY,I)*CONPD4(2,ILAY) ) + &
						( COEF4(3,ILAY,I)*CONPD4(3,ILAY) ) + &
						( COEF4(4,ILAY,I)*CONPD4(4,ILAY) ) + &
						( COEF4(5,ILAY,I)*CONPD4(5,ILAY) ) + &
						( COEF4(6,ILAY,I)*CONPD4(6,ILAY) ) + &
						( COEF4(7,ILAY,I)*CONPD4(7,ILAY) )
					
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
					KFIX = &
						( COEF4( 8,ILAY,I)*FIXED_PRED4( 1,ILAY) ) + &
						( COEF4( 9,ILAY,I)*FIXED_PRED4( 2,ILAY) ) + &
						( COEF4(10,ILAY,I)*FIXED_PRED4( 3,ILAY) ) + &
						( COEF4(11,ILAY,I)*FIXED_PRED4( 4,ILAY) ) + &
						( COEF4(12,ILAY,I)*FIXED_PRED4( 5,ILAY) ) + &
						( COEF4(13,ILAY,I)*FIXED_PRED4( 6,ILAY) ) + &
						( COEF4(14,ILAY,I)*FIXED_PRED4( 7,ILAY) ) + &
						( COEF4(15,ILAY,I)*FIXED_PRED4( 8,ILAY) ) + &
						( COEF4(16,ILAY,I)*FIXED_PRED4( 9,ILAY) ) + &
						( COEF4(17,ILAY,I)*FIXED_PRED4(10,ILAY) ) + &
						( COEF4(18,ILAY,I)*FIXED_PRED4(11,ILAY) )

					KFIX = KFIX * FIXMUL(ILAY)

					 IF (KFIX .LT. 0.0E+0) THEN
						KFIX=0.0E+0
					 ELSEIF (KFIX .GT. 1.0E+1) THEN
						KFIX=1.0E+1
					 ENDIF
				ELSE
					KFIX = KFIX_SAVED(ILAY, I)
				ENDIF

				IF (DO_SPECIES % DO_CO) THEN
					!
					! Compute the CO abs coef
					!
					KCO = &
						( COEF4(19,ILAY,I)*CO_PRED4( 1,ILAY) ) + &
						( COEF4(20,ILAY,I)*CO_PRED4( 2,ILAY) ) + &
						( COEF4(21,ILAY,I)*CO_PRED4( 3,ILAY) ) + &
						( COEF4(22,ILAY,I)*CO_PRED4( 4,ILAY) ) + &
						( COEF4(23,ILAY,I)*CO_PRED4( 5,ILAY) ) + &
						( COEF4(24,ILAY,I)*CO_PRED4( 6,ILAY) ) + &
						( COEF4(25,ILAY,I)*CO_PRED4( 7,ILAY) ) + &
						( COEF4(26,ILAY,I)*CO_PRED4( 8,ILAY) ) + &
						( COEF4(27,ILAY,I)*CO_PRED4( 9,ILAY) ) + &
						( COEF4(28,ILAY,I)*CO_PRED4(10,ILAY) ) + &
						( COEF4(29,ILAY,I)*CO_PRED4(11,ILAY) )

					IF (KCO .LT. 0.0E+0) THEN
						KCO=0.0E+0
					ELSEIF (KCO .GT. 1.0E+1) THEN
						KCO=1.0E+1
					ENDIF
				ELSE
					KCO = KCO_SAVED(ILAY, I)
				ENDIF

				IF (DO_SPECIES % DO_O3) THEN
					!
					! Compute the ozone abs coef
					!
					KOZO = &
						( COEF4(30,ILAY,I)*O3_PRED4(1,ILAY) ) + &
						( COEF4(31,ILAY,I)*O3_PRED4(2,ILAY) ) + &
						( COEF4(32,ILAY,I)*O3_PRED4(3,ILAY) )

					IF (KOZO .LT. 0.0E+0) THEN
						KOZO=0.0E+0
					ELSEIF (KOZO .GT. 1.0E+1) THEN
						KOZO=1.0E+1
					ENDIF
				ELSE
					KOZO = KOZO_SAVED(ILAY, I)
				ENDIF

				IF (DO_SPECIES % DO_H2O) THEN
					!
					! Compute the water abs coef
					!
					KWAT = &
						( COEF4(33,ILAY,I)*H2O_PRED4( 1,ILAY) ) + &
						( COEF4(34,ILAY,I)*H2O_PRED4( 2,ILAY) ) + &
						( COEF4(35,ILAY,I)*H2O_PRED4( 3,ILAY) ) + &
						( COEF4(36,ILAY,I)*H2O_PRED4( 4,ILAY) ) + &
						( COEF4(37,ILAY,I)*H2O_PRED4( 5,ILAY) ) + &
						( COEF4(38,ILAY,I)*H2O_PRED4( 6,ILAY) ) + &
						( COEF4(39,ILAY,I)*H2O_PRED4( 7,ILAY) ) + &
						( COEF4(40,ILAY,I)*H2O_PRED4( 8,ILAY) ) + &
						( COEF4(41,ILAY,I)*H2O_PRED4( 9,ILAY) ) + &
						( COEF4(42,ILAY,I)*H2O_PRED4(10,ILAY) ) + &
						( COEF4(43,ILAY,I)*H2O_PRED4(11,ILAY) ) + &
						( COEF4(44,ILAY,I)*H2O_PRED4(12,ILAY) ) + &
						( COEF4(45,ILAY,I)*H2O_PRED4(13,ILAY) )

					IF (KWAT .LT. 0.0E+0) THEN
						KWAT=0.0E+0
					ELSEIF( KWAT .GT. 1.0E+1) THEN
						KWAT=1.0E+1
					ENDIF
				ELSE
					KWAT = KWAT_SAVED(ILAY, I)
				ENDIF

				!
				!
				! Calc the total layer transmittance
				! 
				!
				!cccc
				! This block is usually commented out and is only uncommented for
				! testing purposes.
				!
				!           kcon=0.0E+0
				!           kfix=0.0E+0
				!           kco =0.0E+0
				!           kozo=0.0E+0
				!           kwat=0.0E+0
				!cccc
				!
				!
				! Calc change in total optical
				! depth due to variable CO2
				!
				IF (DO_SPECIES % DO_CO2) THEN
					IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0) THEN
						DKCO2 = &
							( COFCO2(1,ILAY,ICO2)*TRACEGAS_PRED(1,ILAY) ) + &
							( COFCO2(2,ILAY,ICO2)*TRACEGAS_PRED(2,ILAY) ) + &
							( COFCO2(3,ILAY,ICO2)*TRACEGAS_PRED(3,ILAY) ) + &
							( COFCO2(4,ILAY,ICO2)*TRACEGAS_PRED(4,ILAY) ) + &
							( COFCO2(5,ILAY,ICO2)*TRACEGAS_PRED(5,ILAY) )
						DKCO2=DKCO2*CO2MLT(ILAY)
					ELSE
						DKCO2=0.0
					ENDIF
				ELSE
					DKCO2 = DKCO2_SAVED(ILAY, I)
				ENDIF


				IF (DO_SPECIES % DO_N2O) THEN
					!
					! Calc change in total optical
					! depth due to variable N2O
					!
					IF (LN2O .AND. N2OMLT(ILAY) .NE. 0) THEN
						DKN2O=( COFN2O(1,ILAY,IN2O)*TRACEGAS_PRED(1,ILAY) ) + &
						 ( COFN2O(2,ILAY,IN2O)*TRACEGAS_PRED(2,ILAY) ) + &
						 ( COFN2O(3,ILAY,IN2O)*TRACEGAS_PRED(3,ILAY) ) + &
						 ( COFN2O(4,ILAY,IN2O)*TRACEGAS_PRED(4,ILAY) ) + &
						 ( COFN2O(5,ILAY,IN2O)*TRACEGAS_PRED(5,ILAY) ) + &
						 ( COFN2O(6,ILAY,IN2O)*TRACEGAS_PRED(6,ILAY) ) + &
						 ( COFN2O(7,ILAY,IN2O)*TRACEGAS_PRED(7,ILAY) )
						DKN2O=DKN2O*N2OMLT(ILAY)
					ELSE
						DKN2O=0.0
					ENDIF
				ELSE
					DKN2O = DKN2O_SAVED(ILAY, I)
				ENDIF

				!cc
				! this block for testing
				!      DKCO2=0.0
				!      DKN2O=0.0
				!cc

				! Limit -DK so it can never totally totally cancel KFIX
				DK = DKCO2 + DKN2O
				IF (-DK .GE. KFIX) THEN
					DK = -0.999*KFIX
				ENDIF

				! Calc total layer optical depth
				KLAYER = KCON + KFIX + KCO + KOZO + KWAT + DK
				TAU(ILAY,J)=KLAYER
				!
				! Calc layer-to-space optical depth
				KZ=KZ + KLAYER
				TAUZ(ILAY,J)=KZ

				! Save KCON if DO_SPECIES % DO_ALL  and the loop started from layer 1
				IF (DO_SPECIES % DO_ALL .AND. START_LAYER .EQ. 1) THEN
					KCON_SAVED(ILAY, I) = KCON
					KFIX_SAVED(ILAY, I) = KFIX
					KCO_SAVED(ILAY, I) = KCO
					KOZO_SAVED(ILAY, I) = KOZO
					KWAT_SAVED(ILAY, I) = KWAT

					DKCO2_SAVED(ILAY, I) = DKCO2
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

	END SUBROUTINE CALT4_SOLAR__CALCULATE_OD

END MODULE CALT4_SOLAR
