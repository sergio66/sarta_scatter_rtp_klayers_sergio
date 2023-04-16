!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS
!
!    SUNPAR version with trace gases
!
!F77====================================================================


!ROUTINE NAME:
!    SUNPAR


!ABSTRACT:
!    Calculate the fast transmittance code temperature/amount/angle
!    dependent predictors for a profile at the effective sun angle.


!CALL PROTOCOL:
!    SUNPAR ( LBOT, 
!  $                TEMPERATURE_PROFILE, H2O_PROFILE, O3_PROFILE, CO_PROFILE,
!  $          REF_LAYER_PRESSURE_PROFILE,   SECANG, H2O_CONTINUUM_PRED,
!  $          FIXED_PRED4, FIXED_PRED5, FIXED_PRED6, FIXED_PRED7,
!  $          H2O_PRED4, H2O_PRED5, H2O_PRED6, H2O_PRED7,
!  $          O3_PRED4, O3_PRED5, O3_PRED6, O3_PRED7,
!  $          CO_PRED4, TRACEGAS_PRED )


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INTEGER   LBOT    bottom layer number         none
!    REAL arr  TEMPERATURE_PROFILE   profile temperature         Kelvin
!    REAL arr  CO_PROFILE  prof carbon monoxide amnt   kiloMoles/cm^2
!    REAL arr  O3_PROFILE  profile ozone amount        kiloMoles/cm^2
!    REAL arr  REF_LAYER_PRESSURE_PROFILE    layer pressures             atmospheres
!    REAL arr  H2O_PROFILE  profile water amount        kiloMoles/cm^2
!    REAL arr  REF_TEMPERATURE_PROFILE   reference temperature       Kelvin
!    REAL arr  REF_CO_AMNT  ref carbon monoxide amount  kiloMoles/cm^2
!    REAL arr  REF_O3_AMNT  reference ozone amount      kiloMoles/cm^2
!    REAL arr  REF_H2O_AMNT  reference water amount      kiloMoles/cm^2
!    REAL arr  SECANG  secant of path angle        none


!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL arr  CO_PRED4  carbon monoxide pred set4   various
!    REAL arr  FIXED_PRED4  fixed predictors set4       various
!    REAL arr  FIXED_PRED5  fixed predictors set5       various
!    REAL arr  FIXED_PRED6  fixed predictors set6       various
!    REAL arr  FIXED_PRED7  fixed predictors set7       various
!    REAL arr  O3_PRED4  ozone predictors set4       various
!    REAL arr  O3_PRED5  ozone predictors set5       various
!    REAL arr  O3_PRED6  ozone predictors set6       variou
!    REAL arr  O3_PRED7  ozone predictors set7       various
!    REAL arr  TRACEGAS_PRED  trace gas pert predictors   various
!    REAL arr  H2O_PRED4  water predictors set4       various
!    REAL arr  H2O_PRED5  water predictors set5       various
!    REAL arr  H2O_PRED6  water predictors set6       various
!    REAL arr  H2O_PRED7  water predictors set7       various


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
!    Rapid transmittace algorithm predictors consisting of various gas
!    amount and temperature ratios and offsets relative to a reference
!    profile are calculated.  Only sets 4 - 7 are calculated, as these
!    are the only sets fit for extreme sun angles (up to secant 9).
!
!    ===================================================================
!    The FTC profile variables computed for each layer are:
!
!    ---------------------------------
!    H2O_CONTINUUM_PRED: water continuum predictors (7 terms)
!       1) a*W/Tr^2    2) a*(W/Tr^2)^2   3) a*W/Tr  4) a*W^2/Tr
!       5) a*(W/Tr)^2  6) a*W/Tr^4       7) a*Wr
!
!    -------------------------------
!    Fixed predictors
!
!    FIXED_PRED4: FCOW (11 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a^2*Trz
!       9) a^2*Tr  10) a^3     11) sqrt(a)
!
!    FIXED_PRED5: FWO (11 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
!       9) a^2*Tr  10) sqrt(a) 11) Trz
!
!    FIXED_PRED6: FWO (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) sqrt(a)
!
!    FIXED_PRED7: FWO (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) sqrt(a)
!
!    ---------------------------------
!    Water predictors
!
!    H2O_PRED4: FCOW (13 terms):
!       1) W*a             2) W             3) sqrt(W*a)
!       4) W*a*dT          5) (W*a)^2       6) sqrt(W*a)*dT
!       7) root^4(W*a)     8) W*a*W/Wz      9) W*a^2
!      10) (W*a)^3        11) W*a*Cz*a     12) sqrt(W*a)*W/Wz
!      13) W*a^2*dT
!
!    H2O_PRED5: FWO bfsw (3 terms):
!       1) W*a           2) (W*a)^3/2       3) W*a*dT
!
!    H2O_PRED6: FWO mfmw (7 terms):
!       1) W*a           2) (W*a)^3/2       3) W*a*dT
!       4) (W*a)^2       5) (W*a)^3/2*dT    6) (W*a)^3
!       7) W*a^2
!
!    H2O_PRED7: FWO mfbw (13 terms):
!       1) W*a           2) (W*a)^3/2       3) W*a*dT
!       4) (W*a)^2       5) (W*a)^3/2*dT    6) (W*a)^3
!       7) W*a^2         8) W*a*W/Wz        9) (W*a)^3/2*W/Wz
!      10) (W*a)^5/4    11) (W*a)^2*W/Wz   12) W^2*a
!      13) (W*a)^7/4
!
!    ---------------------------
!    Ozone predictors
!
!    O3_PRED4: FCOW (3 terms):
!       1) O*a         2) sqrt(O*a)     3) O*a*dT
!
!    O3_PRED5: FWO bfsw (1 term):
!       1) O*a
!
!    O3_PRED6: FWO mfmw (1 term):
!       1) O*a
!
!    O3_PRED7: FWO mfbw (1 term):
!       1) O*a
!
!    ---------------------------
!    CO_PRED4: carbon monoxide predictors (11 terms):
!       1) C*a           2) sqrt(C*a)       3) C*a*dT
!       4) (C*a)^2       5) C*a*C/Cz        6) sqrt(C*a)*dT
!       7) root^4(C*a)   8) sqrt(C*a)*C/Cz  9) C
!
!    ---------------------------
!    CO2PRD: CO2 perturbation coefs (4 terms):
!       1) a        2) Tr      3) a*Tr    4) a*Tr^2
!
!    -----
!    where:
!    "a" is the secant of the viewing angle SECANG
!    "Tr" is the temperature ratio TEMPERATURE_PROFILE/REF_TEMPERATURE_PROFILE
!    "Trz" is the pressure weighted temperature ratio above, i.e.
!      the sum i=2 to i=L of { P(i) * ( P(i) -  P(i-1) )* Tr(i-1) }
!      where "P" is the pressure REF_LAYER_PRESSURE_PROFILE and "L" is the layer number, and
!      Trz(L=1)=0
!    "W" is the water amount ratio H2O_PROFILE/REF_H2O_AMNT
!    "dT" is the temperature offset TEMPERATURE_PROFILE-REF_TEMPERATURE_PROFILE
!    "Wz" is the pressure weighted water amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * H2O_PROFILE(i) },
!      divided by the same sum except using REF_H2O_AMNT instead of H2O_PROFILE.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!    "O" is the ozone amount ratio O3_PROFILE/REF_O3_AMNT
!    "C" is the carbon monoxide amount ratio O3_PROFILE/REF_O3_AMNT
!    "Cz" is the pressure weighted CO amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * CO_PROFILE(i) },
!      divided by the same sum except using REF_CO_AMNT instead of CO_PROFILE.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!
!    ===================================================================


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    Assumes vaguely realistic profile amounts and temperatures, else
!    there might be divide by zero problems, etc.


!ROUTINE HISTORY:
!    Date        Programmer     Comments
!    ----------- -------------- ----------------------------------------
!    27 Aug 1997 Scott Hannon   Created from calpar
!    30 Sep 1997 Scott Hannon   Added variable CO2
!    26 Aug 1998 Scott Hannon   Add LBOT to call; loop on LBOT instead
!                               of MAXLAY
!    24 Aug 2000 Scott Hannon   Remove FIXMUL (calc'ed in CALPAR)
!    18 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
!    25 Apr 2003 Scott Hannon   Add SO2
!    23 Jun 2005 Scott Hannon   "trace" version for CO2, SO2, & HNO3,
!                               with all using the same predictors.
!    13 Oct 2005 S.Hannon/C.Barnet bug fix: assign TRACEGAS_PRED 1-7 (was 1-4)
!    21 Dec 2017 Bill Irion	    Conversion to F90. Added IF statements
!                               to skip channel sets that are not needed.
!    02 Jan 2018 Bill Irion     Made INITIALIZE subroutine to calculate and store
!                               those array elements that don't change by retrieval


!END====================================================================

MODULE SUNPAR


	USE INCFTC
	USE REF_PROFILES

	IMPLICIT NONE

	! Elements that remain the same retrieval-to-retrieval
	REAL, PRIVATE, DIMENSION(MAXLAY) :: PDP_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: PNORM_ARRAY

	REAL, PRIVATE, DIMENSION(MAXLAY) :: WZREF_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: CZREF_ARRAY

	LOGICAL :: SUNPAR_INITIALIZED

	! Elements that are saved ONLY when START_LAYER == 1 and SELECT_PRED == "ALL"
	! These are substituted in when this routine is needed to calculate predictors
	! only below a layer > 1 (as would happen in calculating a Jacobian.)

	REAL, PRIVATE, DIMENSION(MAXLAY) :: TR_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: DT_SAVED_ARRAY

	!REAL, PRIVATE, DIMENSION(MAXLAY) :: A_O_SAVED_ARRAY

	REAL, PRIVATE, DIMENSION(MAXLAY) :: TZ_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: WZ_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: CZ_SAVED_ARRAY

	REAL, TARGET, DIMENSION(N1CON, MAXLAY) :: H2O_CONTINUUM_PRED_SAVED
	REAL, TARGET, DIMENSION(N4FIX, MAXLAY) :: FIXED_PRED4_SAVED
	REAL, TARGET, DIMENSION(N5FIX, MAXLAY) :: FIXED_PRED5_SAVED
	REAL, TARGET, DIMENSION(N6FIX, MAXLAY) :: FIXED_PRED6_SAVED
	REAL, TARGET, DIMENSION(N7FIX, MAXLAY) :: FIXED_PRED7_SAVED
	REAL, TARGET, DIMENSION(N4H2O, MAXLAY) :: H2O_PRED4_SAVED
	REAL, TARGET, DIMENSION(N5H2O, MAXLAY) :: H2O_PRED5_SAVED
	REAL, TARGET, DIMENSION(N6H2O, MAXLAY) :: H2O_PRED6_SAVED
	REAL, TARGET, DIMENSION(N7H2O, MAXLAY) :: H2O_PRED7_SAVED
	REAL, TARGET, DIMENSION( N4O3, MAXLAY) :: O3_PRED4_SAVED
	REAL, TARGET, DIMENSION( N5O3, MAXLAY) :: O3_PRED5_SAVED
	REAL, TARGET, DIMENSION( N6O3, MAXLAY) :: O3_PRED6_SAVED
	REAL, TARGET, DIMENSION( N7O3, MAXLAY) :: O3_PRED7_SAVED
	REAL, TARGET, DIMENSION( N4CO, MAXLAY) :: CO_PRED4_SAVED
	REAL, TARGET, DIMENSION(NTRACE,MAXLAY) :: TRACEGAS_PRED_SAVED

	REAL, TARGET, DIMENSION(N1CON, MAXLAY) :: H2O_CONTINUUM_PRED
	REAL, TARGET, DIMENSION(N4FIX, MAXLAY) :: FIXED_PRED4
	REAL, TARGET, DIMENSION(N5FIX, MAXLAY) :: FIXED_PRED5
	REAL, TARGET, DIMENSION(N6FIX, MAXLAY) :: FIXED_PRED6
	REAL, TARGET, DIMENSION(N7FIX, MAXLAY) :: FIXED_PRED7
	REAL, TARGET, DIMENSION(N4H2O, MAXLAY) :: H2O_PRED4
	REAL, TARGET, DIMENSION(N5H2O, MAXLAY) :: H2O_PRED5
	REAL, TARGET, DIMENSION(N6H2O, MAXLAY) :: H2O_PRED6
	REAL, TARGET, DIMENSION(N7H2O, MAXLAY) :: H2O_PRED7
	REAL, TARGET, DIMENSION( N4O3, MAXLAY) :: O3_PRED4
	REAL, TARGET, DIMENSION( N5O3, MAXLAY) :: O3_PRED5
	REAL, TARGET, DIMENSION( N6O3, MAXLAY) :: O3_PRED6
	REAL, TARGET, DIMENSION( N7O3, MAXLAY) :: O3_PRED7
	REAL, TARGET, DIMENSION( N4CO, MAXLAY) :: CO_PRED4
	REAL, TARGET, DIMENSION(NTRACE,MAXLAY) :: TRACEGAS_PRED

CONTAINS

	! Initialize those elements that do not change from retrieval to retrieval
	SUBROUTINE SUNPAR__INITIALIZE

		INTEGER :: L

		PDP_ARRAY(1) = REF_LAYER_PRESSURE_PROFILE(1) * (REF_LAYER_PRESSURE_PROFILE(2) - REF_LAYER_PRESSURE_PROFILE(1))
		PNORM_ARRAY(1) = 0.
		WZREF_ARRAY(1) = PDP_ARRAY(1) + REF_H2O_PROFILE(1)
		CZREF_ARRAY(1) = PDP_ARRAY(1) * REF_CO_PROFILE(1)
		DO L = 2, MAXLAY
			PDP_ARRAY(L) = REF_LAYER_PRESSURE_PROFILE(L)*( REF_LAYER_PRESSURE_PROFILE(L) - REF_LAYER_PRESSURE_PROFILE(L-1) )
			PNORM_ARRAY(L) = PNORM_ARRAY(L-1) + PDP_ARRAY(L)
			WZREF_ARRAY(L) = WZREF_ARRAY(L-1) + PDP_ARRAY(L) * REF_H2O_PROFILE(L)
			CZREF_ARRAY(L) = CZREF_ARRAY(L-1) + PDP_ARRAY(L) * REF_CO_PROFILE(L)
		ENDDO

		SUNPAR_INITIALIZED = .TRUE.
		
		RETURN

	END SUBROUTINE SUNPAR__INITIALIZE

	SUBROUTINE SUNPAR__CALCULATE_PREDICTORS ( LBOT, &
		TEMPERATURE_PROFILE,  H2O_PROFILE, O3_PROFILE, CO_PROFILE, &
		SECANG, H2O_CONTINUUM_PRED_PTR, &
		FIXED_PRED4_PTR, FIXED_PRED5_PTR, FIXED_PRED6_PTR, FIXED_PRED7_PTR, &
		H2O_PRED4_PTR,   H2O_PRED5_PTR,   H2O_PRED6_PTR,   H2O_PRED7_PTR, &
		O3_PRED4_PTR,    O3_PRED5_PTR,    O3_PRED6_PTR,    O3_PRED7_PTR, &
		CO_PRED4_PTR,    TRACEGAS_PRED_PTR, &
		START_LAYER, DO_SPECIES)


		USE AIRS_SARTA_VARIABLES, ONLY: NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7

		USE SELECT_SPECIES

		!
		! ARGUMENTS
		!
		! Input
		INTEGER, INTENT(IN) :: LBOT
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: TEMPERATURE_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) ::  H2O_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) ::  O3_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) ::  CO_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) ::  SECANG

		INTEGER, INTENT(IN) :: START_LAYER
		TYPE (DO_SPECIES_LOGICAL_ARRAY_TYPE), INTENT(IN) :: DO_SPECIES

		! Output
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: H2O_CONTINUUM_PRED_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: FIXED_PRED4_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: FIXED_PRED5_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: FIXED_PRED6_PTR 
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: FIXED_PRED7_PTR 
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: H2O_PRED4_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: H2O_PRED5_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: H2O_PRED6_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: H2O_PRED7_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: O3_PRED4_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: O3_PRED5_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: O3_PRED6_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: O3_PRED7_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: CO_PRED4_PTR
		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: TRACEGAS_PRED_PTR


		!
		! LOCAL VARIABLES
		!
		INTEGER      L
		REAL    PDP
		REAL  PNORM
		REAL     DT
		REAL     TR
		REAL     TZ
		REAL    TRZ
		!       REAL    A_F ! unused so removed 14 Feb 2001
		REAL    A_W
		REAL  WZREF
		REAL     WZ
		REAL   AZ_W
		REAL    A_O
		REAL    A_C
		REAL     CZ
		REAL  CZREF
		REAL   AZ_C
		REAL TJUNKS
		REAL WJUNKA
		REAL WJUNKR
		REAL WJUNKS
		REAL WJUNKZ
		REAL WJUNK4
		REAL OJUNKA
		REAL CJUNKA
		REAL CJUNKR
		REAL CJUNKS
		REAL CJUNKZ


		INTEGER :: INIT_LAYER

		!
		!
		! EXECUTABLE CODE
		!
		!


		! If you're not selecting for any of temperature, H2O, O3 or CO, 
		! return all saved predictors and multipliers without making any new calculations

		IF ( &
			(.NOT. DO_SPECIES % DO_TEMPERATURE) &
			.AND. &
			(.NOT. DO_SPECIES % DO_H2O) &
			.AND. &
			(.NOT. DO_SPECIES % DO_O3) &
			.AND. &
			(.NOT. DO_SPECIES % DO_CO) ) THEN
				H2O_CONTINUUM_PRED_PTR => H2O_CONTINUUM_PRED_SAVED
				TRACEGAS_PRED_PTR => TRACEGAS_PRED_SAVED		
				IF (NCHN4 .GT. 0) THEN
					FIXED_PRED4_PTR => FIXED_PRED4_SAVED
					H2O_PRED4_PTR => H2O_PRED4_SAVED
					O3_PRED4_PTR => O3_PRED4_SAVED
					CO_PRED4_PTR => CO_PRED4_SAVED
				ENDIF
				IF (NCHN5 .GT. 0) THEN
					FIXED_PRED5_PTR => FIXED_PRED5_SAVED
					H2O_PRED5_PTR => H2O_PRED5_SAVED
					O3_PRED5_PTR => O3_PRED5_SAVED
				ENDIF
				IF (NCHN6 .GT. 0) THEN
					FIXED_PRED6_PTR => FIXED_PRED6_SAVED
					H2O_PRED6_PTR => H2O_PRED6_SAVED
					O3_PRED6_PTR => O3_PRED6_SAVED
				ENDIF
				IF (NCHN7 .GT. 0) THEN
					FIXED_PRED7_PTR => FIXED_PRED7_SAVED
					H2O_PRED7_PTR => H2O_PRED7_SAVED
					O3_PRED7_PTR => O3_PRED7_SAVED
				ENDIF
				RETURN
		ENDIF

		IF (START_LAYER .EQ. 1) THEN
			! Initialize the sum terms to zero
			!PNORM=0.0E+0
			TZ=0.0E+0
			WZ=0.0E+0
			CZ=0.0E+0
		ELSE
			INIT_LAYER = START_LAYER - 1
			! Pull the values from saved arrays
			TZ = TZ_SAVED_ARRAY(INIT_LAYER)
			WZ = WZ_SAVED_ARRAY(INIT_LAYER)
			CZ = CZ_SAVED_ARRAY(INIT_LAYER)

			TR = TR_SAVED_ARRAY(INIT_LAYER)
			DT = DT_SAVED_ARRAY(INIT_LAYER)
			!A_O = A_O_SAVED_ARRAY(INIT_LAYER)

			! Restore saved predictors above start layer
			! For now, restore all layers. Later on, we can be more selective
			! depending on what SELECT_PRED is

			IF (DO_SPECIES % DO_H2O) H2O_CONTINUUM_PRED(:, 1:INIT_LAYER) = H2O_CONTINUUM_PRED_SAVED(:, 1:INIT_LAYER)

			IF (DO_SPECIES % DO_TEMPERATURE) TRACEGAS_PRED(:, 1:INIT_LAYER) = TRACEGAS_PRED_SAVED(:, 1:INIT_LAYER)
			
			IF (NCHN4 .GT. 0) THEN
				IF (DO_SPECIES % DO_TEMPERATURE) FIXED_PRED4(:, 1:INIT_LAYER) = FIXED_PRED4_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_H2O) H2O_PRED4(:, 1:INIT_LAYER) = H2O_PRED4_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_O3) O3_PRED4(:, 1:INIT_LAYER) = O3_PRED4_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_CO) CO_PRED4(:, 1:INIT_LAYER) = CO_PRED4_SAVED(:, 1:INIT_LAYER)
			ENDIF

			IF (NCHN5 .GT. 0) THEN
				IF (DO_SPECIES % DO_TEMPERATURE) FIXED_PRED5(:, 1:INIT_LAYER) = FIXED_PRED5_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_H2O) H2O_PRED5(:, 1:INIT_LAYER) = H2O_PRED5_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_O3) O3_PRED5(:, 1:INIT_LAYER) = O3_PRED5_SAVED(:, 1:INIT_LAYER)
			ENDIF

			IF (NCHN6 .GT. 0) THEN
				IF (DO_SPECIES % DO_TEMPERATURE) FIXED_PRED6(:, 1:INIT_LAYER) = FIXED_PRED6_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_H2O) H2O_PRED6(:, 1:INIT_LAYER) = H2O_PRED6_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_O3) O3_PRED6(:, 1:INIT_LAYER) = O3_PRED6_SAVED(:, 1:INIT_LAYER)
			ENDIF

			IF (NCHN7 .GT. 0) THEN
				IF (DO_SPECIES % DO_TEMPERATURE) FIXED_PRED7(:, 1:INIT_LAYER) = FIXED_PRED7_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_H2O) H2O_PRED7(:, 1:INIT_LAYER) = H2O_PRED7_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_O3) O3_PRED7(:, 1:INIT_LAYER) = O3_PRED7_SAVED(:, 1:INIT_LAYER)
			ENDIF


		ENDIF


		! 
		! Loop over the layers
		! 
		DO L = START_LAYER, LBOT

			! 
			! Calculate the basic profile
			! dependent predictors.
			! 
			IF (L .EQ. 1) THEN
				!PDP=REF_LAYER_PRESSURE_PROFILE(1)*( REF_LAYER_PRESSURE_PROFILE(2) - REF_LAYER_PRESSURE_PROFILE(1))
				PDP=PDP_ARRAY(1)
				TRZ=0.0E+0
			ELSE
				!PDP=REF_LAYER_PRESSURE_PROFILE(L)*( REF_LAYER_PRESSURE_PROFILE(L) - REF_LAYER_PRESSURE_PROFILE(L-1) )
				!PNORM=PNORM + PDP
				PDP = PDP_ARRAY(L)
				PNORM = PNORM_ARRAY(L)

				! Note: TRZ use layer-above terms
				TZ = TZ + PDP*TR
				TRZ = TZ / PNORM
			ENDIF

			! Temperature terms
			DT=TEMPERATURE_PROFILE(L) - REF_TEMPERATURE_PROFILE(L)
			TR=TEMPERATURE_PROFILE(L)/REF_TEMPERATURE_PROFILE(L)

			! Water terms
			A_W=H2O_PROFILE(L)/REF_H2O_PROFILE(L)
			!WZREF=WZREF + PDP*REF_H2O_PROFILE(L)
			WZ=WZ + PDP*H2O_PROFILE(L)
			!AZ_W=WZ/WZREF
			AZ_W=WZ/WZREF_ARRAY(L)

			! Ozone terms
			A_O=O3_PROFILE(L)/REF_O3_PROFILE(L)

			! Carbon monoxide terms
			A_C=CO_PROFILE(L)/REF_CO_PROFILE(L)
			!CZREF=CZREF + PDP*REF_CO_PROFILE(L)
			CZ=CZ + PDP*CO_PROFILE(L)
			!AZ_C=CZ/CZREF
			AZ_C=CZ/CZREF_ARRAY(L)


			IF ((START_LAYER .EQ. 1) .AND. DO_SPECIES % DO_ALL) THEN
				! Save above layer terms in array
				! Save  WZ, XZ, OZ, CZ, MZ, TZ, TOZ and TMZ(L) 
				TR_SAVED_ARRAY(L) = TR
				DT_SAVED_ARRAY(L) = DT

				!A_O_SAVED_ARRAY(L) = A_O

				TZ_SAVED_ARRAY(L) = TZ
				WZ_SAVED_ARRAY(L) = WZ
				CZ_SAVED_ARRAY(L) = CZ
			ENDIF

			!
			! 
			! Load up the predictors
			!
			!
			! 
			! Fixed
			! 
			IF (DO_SPECIES % DO_TEMPERATURE) THEN
				TJUNKS=TR*TR

				IF (NCHN4 .GT. 0) THEN
					FIXED_PRED4(1,L)=SECANG(L)
					FIXED_PRED4(2,L)=SECANG(L)*SECANG(L)
					FIXED_PRED4(3,L)=SECANG(L)*TR
					FIXED_PRED4(4,L)=SECANG(L)*TJUNKS
					FIXED_PRED4(5,L)=TR
					FIXED_PRED4(6,L)=TJUNKS
					FIXED_PRED4(7,L)=SECANG(L)*TRZ
					FIXED_PRED4(8,L)=SECANG(L)*SECANG(L)*TRZ
					FIXED_PRED4(9,L)=SECANG(L)*SECANG(L)*TR
					FIXED_PRED4(10,L)=SECANG(L)*SECANG(L)*SECANG(L)
					FIXED_PRED4(11,L)=SQRT(SECANG(L))
				ENDIF

				! Fixed predictors for FWO sun bfsw = set5
				IF (NCHN5 .GT. 0) THEN
					FIXED_PRED5(1,L)=SECANG(L)
					FIXED_PRED5(2,L)=SECANG(L)*SECANG(L)
					FIXED_PRED5(3,L)=SECANG(L)*TR
					FIXED_PRED5(4,L)=SECANG(L)*TJUNKS
					FIXED_PRED5(5,L)=TR
					FIXED_PRED5(6,L)=TJUNKS
					FIXED_PRED5(7,L)=SECANG(L)*TRZ
					FIXED_PRED5(8,L)=SECANG(L)*TRZ/TR
					FIXED_PRED5(9,L)=SECANG(L)*SECANG(L)*TR
					FIXED_PRED5(10,L)=SQRT(SECANG(L))
					FIXED_PRED5(11,L)=TRZ
				ENDIF

				! Fixed predictors for FWO sun mfmw = set6
				IF (NCHN6 .GT. 0) THEN
					FIXED_PRED6(1,L)=SECANG(L)
					FIXED_PRED6(2,L)=SECANG(L)*SECANG(L)
					FIXED_PRED6(3,L)=SECANG(L)*TR
					FIXED_PRED6(4,L)=SECANG(L)*TJUNKS
					FIXED_PRED6(5,L)=TR
					FIXED_PRED6(6,L)=TJUNKS
					FIXED_PRED6(7,L)=SECANG(L)*TRZ
					FIXED_PRED6(8,L)=SQRT(SECANG(L))
				ENDIF

				! Fixed predictors for FWO sun mfbw = set7
				IF (NCHN7 .GT. 0) THEN
					FIXED_PRED7(1,L)=SECANG(L)
					FIXED_PRED7(2,L)=SECANG(L)*SECANG(L)
					FIXED_PRED7(3,L)=SECANG(L)*TR
					FIXED_PRED7(4,L)=SECANG(L)*TJUNKS
					FIXED_PRED7(5,L)=TR
					FIXED_PRED7(6,L)=TJUNKS
					FIXED_PRED7(7,L)=SECANG(L)*TRZ
					FIXED_PRED7(8,L)=SQRT(SECANG(L))
				ENDIF
			ENDIF

			!
			!         Ozone
			! 
			IF (DO_SPECIES % DO_O3) THEN
				OJUNKA=SECANG(L)*A_O

				! ozone predictors for FCOW = set4
				IF (NCHN4 .GT. 0) THEN
					O3_PRED4(1,L)=OJUNKA
					O3_PRED4(2,L)=SQRT( OJUNKA )
					O3_PRED4(3,L)=OJUNKA*DT
				ENDIF

				! ozone predictors for FWO sun bfsw = set5
				O3_PRED5(1,L)=OJUNKA

				! ozone predictors for FWO sun mfmw = set6
				O3_PRED6(1,L)=OJUNKA

				! ozone predictors for FWO sun mfbw = set7
				O3_PRED7(1,L)=OJUNKA
			ENDIF
		
			! 
			! Water
			!
			IF (DO_SPECIES % DO_H2O) THEN
				WJUNKA=SECANG(L)*A_W
				WJUNKR=SQRT( WJUNKA )
				WJUNKS=WJUNKA*WJUNKA
				WJUNKZ=WJUNKA*A_W/AZ_W
				WJUNK4=SQRT( WJUNKR )

				! water predictors for FCOW = set4
				IF (NCHN4 .GT. 0) THEN
					H2O_PRED4( 1,L)=WJUNKA
					H2O_PRED4( 2,L)=A_W
					H2O_PRED4( 3,L)=WJUNKR
					H2O_PRED4( 4,L)=WJUNKA*DT
					H2O_PRED4( 5,L)=WJUNKS
					H2O_PRED4( 6,L)=WJUNKR*DT
					H2O_PRED4( 7,L)=WJUNK4
					H2O_PRED4( 8,L)=WJUNKZ
					H2O_PRED4( 9,L)=WJUNKA*SECANG(L)
					H2O_PRED4(10,L)=WJUNKS*WJUNKA
					H2O_PRED4(11,L)=WJUNKA*AZ_C*SECANG(L)
					H2O_PRED4(12,L)=WJUNKZ/WJUNKR
					H2O_PRED4(13,L)=WJUNKA*DT*SECANG(L)
				ENDIF

				! Water predictors for FWO sun bfsw = set5
				IF (NCHN5 .GT. 0) THEN
					H2O_PRED5( 1,L)=WJUNKA
					H2O_PRED5( 2,L)=WJUNKA*WJUNKR
					H2O_PRED5( 3,L)=WJUNKA*DT
				ENDIF

				! Water predictors for FWO sun mfmw = set6
				IF (NCHN6 .GT. 0) THEN
					H2O_PRED6( 1,L)=WJUNKA
					H2O_PRED6( 2,L)=WJUNKA*WJUNKR
					H2O_PRED6( 3,L)=WJUNKA*DT
					H2O_PRED6( 4,L)=WJUNKS
					H2O_PRED6( 5,L)=WJUNKA*WJUNKR*DT
					H2O_PRED6( 6,L)=WJUNKA*WJUNKS
					H2O_PRED6( 7,L)=WJUNKA*SECANG(L)
				ENDIF

				! Water predictors for FWO sun mfbw = set7
				IF (NCHN7 .GT. 0) THEN
					H2O_PRED7( 1,L)=WJUNKA
					H2O_PRED7( 2,L)=WJUNKA*WJUNKR
					H2O_PRED7( 3,L)=WJUNKA*DT
					H2O_PRED7( 4,L)=WJUNKS
					H2O_PRED7( 5,L)=WJUNKA*WJUNKR*DT
					H2O_PRED7( 6,L)=WJUNKA*WJUNKS
					H2O_PRED7( 7,L)=WJUNKA*SECANG(L)
					H2O_PRED7( 8,L)=WJUNKZ
					H2O_PRED7( 9,L)=WJUNKZ*WJUNKR
					H2O_PRED7(10,L)=WJUNKA*WJUNK4
					H2O_PRED7(11,L)=WJUNKA*WJUNKZ
					H2O_PRED7(12,L)=WJUNKA*A_W
					H2O_PRED7(13,L)=WJUNKS/WJUNK4
				ENDIF

				!
				!         Water continuum (for FWO, FOW, FMW, FCOW)
				!
				TJUNKS = TR * TR
				H2O_CONTINUUM_PRED(1,L)=WJUNKA/TJUNKS
				H2O_CONTINUUM_PRED(2,L)=H2O_CONTINUUM_PRED(1,L)*A_W/TJUNKS
				H2O_CONTINUUM_PRED(3,L)=WJUNKA/TR
				H2O_CONTINUUM_PRED(4,L)=H2O_CONTINUUM_PRED(3,L)*A_W
				H2O_CONTINUUM_PRED(5,L)=H2O_CONTINUUM_PRED(1,L)*A_W
				H2O_CONTINUUM_PRED(6,L)=H2O_CONTINUUM_PRED(1,L)/TJUNKS
				H2O_CONTINUUM_PRED(7,L)=WJUNKA
			ENDIF

			!         ---------------
			!         Carbon monoxide for FCOW = set4
			!         ---------------
			IF (DO_SPECIES % DO_CO) THEN
				IF (NCHN4 .GT. 0) THEN
					CJUNKA=SECANG(L)*A_C
					CJUNKR=SQRT( CJUNKA )
					CJUNKS=CJUNKA*CJUNKA
					CJUNKZ=CJUNKA*A_C/AZ_C
				!
					CO_PRED4(1,L)=CJUNKA
					CO_PRED4(2,L)=CJUNKR
					CO_PRED4(3,L)=CJUNKA*DT
					CO_PRED4(4,L)=CJUNKS
					CO_PRED4(5,L)=CJUNKZ
					CO_PRED4(6,L)=CJUNKR*DT
					CO_PRED4(7,L)=SQRT( CJUNKR )
					CO_PRED4(8,L)=CJUNKZ/CJUNKR
					CO_PRED4(9,L)=A_C
					CO_PRED4(10,L)=CJUNKA*SECANG(L)
					CO_PRED4(11,L)=CJUNKR*SECANG(L)
				ENDIF
			ENDIF

			!
			!         trace gas perturbation coefs
			! 
			!         The first 4 trace predictors are used by all trace gases
			IF (DO_SPECIES % DO_TEMPERATURE) THEN
				TRACEGAS_PRED(1,L)=SECANG(L)
				TRACEGAS_PRED(2,L)=TR
				TRACEGAS_PRED(3,L)=SECANG(L)*TR
				TRACEGAS_PRED(4,L)=SECANG(L)*TJUNKS
				!         The last 3 trace predictors are only used by N2O
				TRACEGAS_PRED(5,L)=SECANG(L)*SECANG(L)
				TRACEGAS_PRED(6,L)=1.0
				TRACEGAS_PRED(7,L)=SQRT( SECANG(L) )
			ENDIF

		ENDDO
		! End loop over layers

		! Save the predictors IF calculated for all layers and all constituents. 
		IF ((START_LAYER .EQ. 1) .AND. DO_SPECIES % DO_ALL) THEN

			H2O_CONTINUUM_PRED_SAVED(:, :) = H2O_CONTINUUM_PRED(:, :)

			IF (NCHN4 .GT. 0) THEN
				FIXED_PRED4_SAVED(:, :) = FIXED_PRED4(:, :)
				H2O_PRED4_SAVED(:, :) = H2O_PRED4(:, :)
				O3_PRED4_SAVED(:, :) = O3_PRED4(:, :)
				CO_PRED4_SAVED(:, :) = CO_PRED4(:, :)
			ENDIF

			IF (NCHN5 .GT. 0) THEN
				FIXED_PRED5_SAVED(:, :) = FIXED_PRED5(:, :)
				H2O_PRED5_SAVED(:, :) = H2O_PRED5(:, :)
				O3_PRED5_SAVED(:, :) = O3_PRED5(:, :)
			ENDIF

			IF (NCHN6 .GT. 0) THEN
				FIXED_PRED6_SAVED(:, :) = FIXED_PRED6(:, :)
				H2O_PRED6_SAVED(:, :) = H2O_PRED6(:, :)
				O3_PRED6_SAVED(:, :) = O3_PRED6(:, :)
			ENDIF
 
			IF (NCHN7 .GT. 0) THEN
				FIXED_PRED7_SAVED(:, :) = FIXED_PRED7(:, :)
				H2O_PRED7_SAVED(:, :) = H2O_PRED7(:, :)
				O3_PRED7_SAVED(:, :) = O3_PRED7(:, :)
			ENDIF

			TRACEGAS_PRED_SAVED(:, :) = TRACEGAS_PRED(:, :)
		ENDIF

		! Now decide what array to assign the pointers to

		IF (DO_SPECIES % DO_TEMPERATURE) THEN
			FIXED_PRED4_PTR => FIXED_PRED4
			FIXED_PRED5_PTR => FIXED_PRED5
			FIXED_PRED6_PTR => FIXED_PRED6
			FIXED_PRED7_PTR => FIXED_PRED7
			TRACEGAS_PRED_PTR => TRACEGAS_PRED
		ELSE
			FIXED_PRED4_PTR => FIXED_PRED4_SAVED
			FIXED_PRED5_PTR => FIXED_PRED5_SAVED
			FIXED_PRED6_PTR => FIXED_PRED6_SAVED
			FIXED_PRED7_PTR => FIXED_PRED7_SAVED
			TRACEGAS_PRED_PTR => TRACEGAS_PRED_SAVED
		ENDIF

		IF (DO_SPECIES % DO_H2O) THEN
			H2O_CONTINUUM_PRED_PTR => H2O_CONTINUUM_PRED
			H2O_PRED4_PTR => H2O_PRED4
			H2O_PRED5_PTR => H2O_PRED5
			H2O_PRED6_PTR => H2O_PRED6
			H2O_PRED7_PTR => H2O_PRED7
		ELSE
			H2O_CONTINUUM_PRED_PTR => H2O_CONTINUUM_PRED_SAVED
			H2O_PRED4_PTR => H2O_PRED4_SAVED
			H2O_PRED5_PTR => H2O_PRED5_SAVED
			H2O_PRED6_PTR => H2O_PRED6_SAVED
			H2O_PRED7_PTR => H2O_PRED7_SAVED
		ENDIF

		IF (DO_SPECIES % DO_O3) THEN
			O3_PRED4_PTR => O3_PRED4
			O3_PRED5_PTR => O3_PRED5
			O3_PRED6_PTR => O3_PRED6
			O3_PRED7_PTR => O3_PRED7
		ELSE
			O3_PRED4_PTR => O3_PRED4_SAVED
			O3_PRED5_PTR => O3_PRED5_SAVED
			O3_PRED6_PTR => O3_PRED6_SAVED
			O3_PRED7_PTR => O3_PRED7_SAVED
		ENDIF

		IF (DO_SPECIES % DO_CO) THEN
			CO_PRED4_PTR => CO_PRED4
		ELSE
			CO_PRED4_PTR => CO_PRED4_SAVED
		ENDIF

		RETURN

	END SUBROUTINE SUNPAR__CALCULATE_PREDICTORS

END MODULE SUNPAR
