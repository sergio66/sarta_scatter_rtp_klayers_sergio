!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS
!
!    CALT6_SOLAR (set6=FWO sun mfmw) version for trace gases (no HNO3)
!
!F77====================================================================


!ROUTINE NAME:
!    CALT6_SOLAR


!ABSTRACT:
!    Calculate the transmittance for set6 using the predictors and the
!    fast transmittance coefficients.


!CALL PROTOCOL:
!    CALT6_SOLAR( INDCHN, NLAY, NCHN6, CLIST6, COEF6,
!       FIXMUL, CONPD6, FIXED_PRED6, H2O_PRED6, O3_PRED6, TRACEGAS_PRED,
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
!    REAL arr  FIXED_PRED6  set6 fixed gases preds      various
!    REAL arr  H2O_PRED6  set6 water predictors       various
!    REAL arr  O3_PRED6  set6 ozone predictors       various
!    REAL arr  TRACEGAS_PRED  trace gas pert predictors   various
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
!       k_con = the sum i=1 to 5 of { COEF(i)*H2O_CONTINUUM_PRED(i) }
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
!      "H2O_CONTINUUM_PRED" are the water continuum predictors H2O_CONTINUUM_PRED
!      "FPRED" are the fixed gases predictors FIXED_PRED5
!      "WPRED" are the water lines predictors H2O_PRED5
!      "OPRED" are the ozone predictors O3_PRED5
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

MODULE CALT6_SOLAR

	USE INCFTC

 	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: TAU_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: TAUZ_SAVED

	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KCON_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KFIX_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KOZO_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: KWAT_SAVED

	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: DKCO2_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: DKSO2_SAVED
	REAL, PRIVATE, DIMENSION(:, :), ALLOCATABLE :: DKN2O_SAVED


 	LOGICAL :: CALT6_SOLAR_ALLOCATED

CONTAINS

	SUBROUTINE CALT6_SOLAR__ALLOCATE(NCHN6)

		INTEGER, INTENT(IN) :: NCHN6

		IF (ALLOCATED(TAU_SAVED)) DEALLOCATE(TAU_SAVED)
		IF (ALLOCATED(TAUZ_SAVED)) DEALLOCATE(TAUZ_SAVED)

		IF (ALLOCATED(KCON_SAVED)) 	DEALLOCATE(KCON_SAVED)
		IF (ALLOCATED(KFIX_SAVED)) 	DEALLOCATE(KFIX_SAVED)
		IF (ALLOCATED(KOZO_SAVED)) 	DEALLOCATE(KOZO_SAVED)
		IF (ALLOCATED(KWAT_SAVED)) 	DEALLOCATE(KWAT_SAVED)

		IF (ALLOCATED(DKCO2_SAVED)) DEALLOCATE(DKCO2_SAVED)
		IF (ALLOCATED(DKSO2_SAVED)) DEALLOCATE(DKSO2_SAVED)
		IF (ALLOCATED(DKN2O_SAVED)) DEALLOCATE(DKN2O_SAVED)

		ALLOCATE(TAU_SAVED(MAXLAY, NCHN6))
		ALLOCATE(TAUZ_SAVED(MAXLAY, NCHN6))

		ALLOCATE(KCON_SAVED(MAXLAY, NCHN6))
		ALLOCATE(KFIX_SAVED(MAXLAY, NCHN6))
		ALLOCATE(KOZO_SAVED(MAXLAY, NCHN6))
		ALLOCATE(KWAT_SAVED(MAXLAY, NCHN6))

		ALLOCATE(DKCO2_SAVED(MAXLAY, NCHN6))
		ALLOCATE(DKSO2_SAVED(MAXLAY, NCHN6))
		ALLOCATE(DKN2O_SAVED(MAXLAY, NCHN6))

		CALT6_SOLAR_ALLOCATED = .TRUE.

		PRINT *, "CALT6_SOLAR allocated"	

	END SUBROUTINE CALT6_SOLAR__ALLOCATE


	SUBROUTINE CALT6_SOLAR__CALCULATE_OD ( INDCHN, NLAY, NCHN6, CLIST6, &
		COEF6, FIXMUL, CONPD6, FIXED_PRED6, H2O_PRED6, O3_PRED6, TRACEGAS_PRED, &
		INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
		INDN2O, COFN2O, N2OMLT, TAU, TAUZ, START_LAYER, DO_SPECIES)

		USE SELECT_SPECIES
		IMPLICIT NONE

		!-----------------------------------------------------------------------
		! ARGUMENTS
		!-----------------------------------------------------------------------
		! Input
		INTEGER INDCHN(MXCHAN)
		INTEGER   NLAY
		INTEGER  NCHN6
		INTEGER CLIST6(MXCHN6)
		REAL  COEF6(N6COEF,MAXLAY,MXCHN6)
		REAL FIXMUL(MAXLAY)
		REAL CONPD6( N6CON,MAXLAY)
		REAL FIXED_PRED6( N6FIX,MAXLAY)
		REAL H2O_PRED6( N6H2O,MAXLAY)
		REAL O3_PRED6(  N6O3,MAXLAY)
		REAL TRACEGAS_PRED(NTRACE,MAXLAY)
		INTEGER INDCO2(MXCHAN)
		REAL COFCO2(  NCO2,MAXLAY,MXCHNC)
		REAL CO2MLT(MAXLAY)
		INTEGER INDSO2(MXCHAN)
		REAL COFSO2(  NSO2,MAXLAY,MXCHNS)
		REAL SO2MLT(MAXLAY)
		INTEGER INDN2O(MXCHAN)
		REAL COFN2O(  NN2O,MAXLAY,MXCHNN)
		REAL N2OMLT(MAXLAY)


		INTEGER, INTENT(IN) :: START_LAYER
		TYPE (DO_SPECIES_LOGICAL_ARRAY_TYPE), INTENT(IN) :: DO_SPECIES

		! Output
		REAL    TAU(MAXLAY,MXCHAN)
		REAL   TAUZ(MAXLAY,MXCHAN)


		!
		! LOCAL VARIABLES
		!
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

		INTEGER :: INIT_LAYER ! (= START_LAYER - 1)

		!
		!
		! EXECUTABLE CODE
		!
		!

		!---------------------------
		! Loop on channel (frequency)
		!---------------------------
		DO I=1,NCHN6

			! Index for TAU
			J=INDCHN( CLIST6(I) )

			! Determine whether or not to do variable CO2
			ICO2=INDCO2( CLIST6(I) )
			IF (ICO2 .GT. 0) THEN
				LCO2=.TRUE.
			ELSE
				LCO2=.FALSE.
			ENDIF

			! Determine whether or not to do variable SO2
			ISO2=INDSO2( CLIST6(I) )
			IF (ISO2 .GT. 0) THEN
				LSO2=.TRUE.
			ELSE
				LSO2=.FALSE.
			ENDIF

			! Determine whether or not to do variable N2O
			IN2O=INDN2O( CLIST6(I) )
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
			!-
			DO ILAY= START_LAYER, NLAY

				IF (DO_SPECIES % DO_H2O) THEN
					!
					! Compute the water continuum
					!
					KCON = &
						( COEF6(1,ILAY,I)*CONPD6(1,ILAY) ) + &
						( COEF6(2,ILAY,I)*CONPD6(2,ILAY) ) + &
						( COEF6(3,ILAY,I)*CONPD6(3,ILAY) ) + &
						( COEF6(4,ILAY,I)*CONPD6(4,ILAY) ) + &
						( COEF6(5,ILAY,I)*CONPD6(5,ILAY) ) + &
						( COEF6(6,ILAY,I)*CONPD6(6,ILAY) ) + &
						( COEF6(7,ILAY,I)*CONPD6(7,ILAY) )

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
						( COEF6( 8,ILAY,I)*FIXED_PRED6( 1,ILAY) ) + &
						( COEF6( 9,ILAY,I)*FIXED_PRED6( 2,ILAY) ) + &
						( COEF6(10,ILAY,I)*FIXED_PRED6( 3,ILAY) ) + &
						( COEF6(11,ILAY,I)*FIXED_PRED6( 4,ILAY) ) + &
						( COEF6(12,ILAY,I)*FIXED_PRED6( 5,ILAY) ) + &
						( COEF6(13,ILAY,I)*FIXED_PRED6( 6,ILAY) ) + &
						( COEF6(14,ILAY,I)*FIXED_PRED6( 7,ILAY) ) + &
						( COEF6(15,ILAY,I)*FIXED_PRED6( 8,ILAY) )

					KFIX=KFIX*FIXMUL(ILAY)

					IF (KFIX .LT. 0.0E+0) THEN
						KFIX=0.0E+0
					ELSEIF (KFIX .GT. 1.0E+1) THEN
						KFIX=1.0E+1
					ENDIF
				ELSE
					KFIX = KFIX_SAVED(ILAY, I)
				ENDIF

				IF (DO_SPECIES % DO_H2O) THEN
					! 
					! Compute the water abs coef
					!
					KWAT=( COEF6(16,ILAY,I)*H2O_PRED6( 1,ILAY) ) + &
						( COEF6(17,ILAY,I)*H2O_PRED6( 2,ILAY) ) + &
						( COEF6(18,ILAY,I)*H2O_PRED6( 3,ILAY) ) + &
						( COEF6(19,ILAY,I)*H2O_PRED6( 4,ILAY) ) + &
						( COEF6(20,ILAY,I)*H2O_PRED6( 5,ILAY) ) + &
						( COEF6(21,ILAY,I)*H2O_PRED6( 6,ILAY) ) + &
						( COEF6(22,ILAY,I)*H2O_PRED6( 7,ILAY) )
				
					IF (KWAT .LT. 0.0E+0) THEN
						KWAT=0.0E+0
					ELSEIF (KWAT .GT. 1.0E+1) THEN
						KWAT=1.0E+1
					ENDIF

				ELSE
					KWAT = KWAT_SAVED(ILAY, I)
				ENDIF

				IF (DO_SPECIES % DO_H2O) THEN
					!
					! Compute the ozone abs coef
					!
					KOZO=( COEF6(23,ILAY,I)*O3_PRED6(1,ILAY) )

					IF (KOZO .LT. 0.0E+0) THEN
						KOZO=0.0E+0
					ELSEIF (KOZO .GT. 1.0E+1) THEN
						KOZO=1.0E+1
					ENDIF
				ELSE
					KOZO = KOZO_SAVED(ILAY, I)
				ENDIF

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
				!           kwat=0.0E+0
				!           kozo=0.0E+0
				!cccc
				!
				!
				! Calc change in total optical
				! depth due to variable CO2
				!
				IF (DO_SPECIES % DO_CO2) THEN
					IF (LCO2 .AND. CO2MLT(ILAY) .NE. 0.0) THEN
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

				IF (DO_SPECIES % DO_SO2) THEN
					!
					! Calc change in total optical
					! depth due to variable SO2
					!
					IF (LSO2 .AND. SO2MLT(ILAY) .NE. 0.0) THEN
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

				IF (DO_SPECIES % DO_N2O) THEN
					!
					! Calc change in total optical
					! depth due to variable N2O
					!
					IF (LN2O .AND. N2OMLT(ILAY) .NE. 0.0) THEN
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
						DKN2O=0.0
					ENDIF
				ELSE
					DKN2O = DKN2O_SAVED(ILAY, I)
				ENDIF

				!cc
				! this block for testing
				!      DKCO2=0.0
				!      DKSO2=0.0
				!      DKN2O=0.0
				!cc

				! Limit -DK so it can never totally totally cancel KFIX
				DK = DKCO2 + DKSO2 + DKN2O
				IF (-DK .GE. KFIX) THEN
					DK = -0.999*KFIX
				ENDIF

				! Calc total layer optical depth
				KLAYER = KCON + KFIX + KWAT + KOZO + DK
				TAU(ILAY,J)=KLAYER

				! Calc layer-to-space optical depth
				KZ=KZ + KLAYER
				TAUZ(ILAY,J)=KZ

				! Save KCON if DO_SPECIES % DO_ALL  and the loop started from layer 1
				IF (DO_SPECIES % DO_ALL .AND. START_LAYER .EQ. 1) THEN
					KCON_SAVED(ILAY, I) = KCON
					KFIX_SAVED(ILAY, I) = KFIX
					KOZO_SAVED(ILAY, I) = KOZO
					KWAT_SAVED(ILAY, I) = KWAT

					DKCO2_SAVED(ILAY, I) = DKCO2
					DKSO2_SAVED(ILAY, I) = DKSO2
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

	END SUBROUTINE CALT6_SOLAR__CALCULATE_OD

END MODULE CALT6_SOLAR
