!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS
!
!    CALPAR
!
!F90====================================================================


!ROUTINE NAME:
!    CALPAR


!ABSTRACT:
!    Calculate the fast transmittance code temperature/amount/angle
!    dependent variables for a profile.


!CALL PROTOCOL:
!    CALPAR (LBOT, REF_TEMPERATURE_PROFILE,REF_CO2_PROFILE,REF_H2O_PROFILE,REF_O3_PROFILE,REF_CO_PROFILE,REF_CH4_PROFILE,REF_N2O_PROFILE,
!  $               TEMPERATURE_PROFILE,CO2_PROFILE,H2O_PROFILE,O3_PROFILE,CO_PROFILE,CH4_PROFILE,N2O_PROFILE,
!  $          RES,SECANG,  ALAT,    FX, REF_LAYER_THICKNESS_PROFILE,
!  $         LCO2,  LN2O,  LSO2, LHNO3,LCO2PM,FIXMUL,H2O_CONTINUUM_PRED,
!  $       FIXED_PRED1,FIXED_PRED2,FIXED_PRED3,FIXED_PRED4,FIXED_PRED5,FIXED_PRED6,FIXED_PRED7,
!  $       H2O_PRED1,H2O_PRED2,H2O_PRED3,H2O_PRED4,H2O_PRED5,H2O_PRED6,H2O_PRED7,
!  $       O3_PRED1,O3_PRED2,       O3_PRED4,O3_PRED5,O3_PRED6,O3_PRED7,
!  $       CH4_PRED3,CO_PRED4,TRACEGAS_PRED,CO2MLT,SO2MLT,HNOMLT,N2OMLT )


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INTEGER   LBOT    bottom layer number         none
!    LOGICAL   LCO2    CO2 profile switch          none
!    LOGICAL   LN2O    N2O profile switch          none
!    LOGICAL   LSO2    SO2 profile switch          none
!    LOGICAL   LHNO3   HNO3 profile switch         none
!    LOGICAL   LCO2PM  CO2 ppmv profile switch     none
!    REAL      ALAT    profile latitude            degrees (-90 to +90)
!    REAL arr  REF_LAYER_THICKNESS_PROFILE   ref prof layer thickness    meters
!    REAL arr  FX      fixed gases adjustment      none
!    REAL arr  TEMPERATURE_PROFILE   profile temperature         K
!    REAL arr  CO_PROFILE  prof carbon monoxide amnt   kiloMoles/cm^2
!    REAL arr  CO2_PROFILE  profile CO2 gas amount      kiloMoles/cm^2
!    REAL arr  HNO3_PROFILE  profile HNO3 gas amount     kiloMoles/cm^2
!    REAL arr  CH4_PROFILE  profile methane amount      kiloMoles/cm^2
!    REAL arr  N2O_PROFILE  profile N2O amount          kiloMoles/cm^2
!    REAL arr  O3_PROFILE  profile ozone amount        kiloMoles/cm^2
!    REAL arr  REF_LAYER_PRESSURE_PROFILE    layer pressures             atm
!    REAL arr  SO2_PROFILE  profile SO2 amount          kiloMoles/cm^2
!    REAL arr  H2O_PROFILE  profile water amount        kiloMoles/cm^2
!    REAL arr  REF_TEMPERATURE_PROFILE   reference temperature       K
!    REAL arr  REF_CO_PROFILE  ref carbon monoxide amount  kiloMoles/cm^2
!    REAL arr  REF_CO2_PROFILE  reference CO2 amount        kiloMoles/cm^2
!    REAL arr  REF_HNO3_PROFILE  reference HNO3 amount       kiloMoles/cm^2
!    REAL arr  REF_CH4_PROFILE  reference methane amount    kiloMoles/cm^2
!    REAL arr  REF_N2O_PROFILE  reference N2O amount        kiloMoles/cm^2
!    REAL arr  REF_O3_PROFILE  reference ozone amount      kiloMoles/cm^2
!    REAL arr  REF_SO2_PROFILE  reference SO2 amount        kiloMoles/cm^2
!    REAL arr  REF_H2O_PROFILE  reference water amount      kiloMoles/cm^2
!    REAL arr  SECANG  secant of path angle        none


!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL arr  CO2MLT  CO2 multiplier              none
!    REAL arr  CO_PRED4  carbon monoxide pred set4   various
!    REAL arr  FIXMUL  fixed amount multiplier     none
!    REAL arr  FIXED_PRED1  fixed predictors set1       various
!    REAL arr  FIXED_PRED2  fixed predictors set2       various
!    REAL arr  FIXED_PRED3  fixed predictors set3       various
!    REAL arr  FIXED_PRED4  fixed predictors set4       various
!    REAL arr  FIXED_PRED5  fixed predictors set5       various
!    REAL arr  FIXED_PRED6  fixed predictors set6       various
!    REAL arr  FIXED_PRED7  fixed predictors set7       various
!    REAL arr  HNOMLT  HNO3 multiplier             none
!    REAL arr  CH4_PRED3  methane predictors set3     various
!    REAL arr  N2OMLT  N2O multiplier              none
!    REAL arr  O3_PRED1  ozone predictors set1       various
!    REAL arr  O3_PRED2  ozone predictors set2       various
!    REAL arr  O3_PRED4  ozone predictors set4       various
!    REAL arr  O3_PRED5  ozone predictors set5       various
!    REAL arr  O3_PRED6  ozone predictors set6       variou
!    REAL arr  O3_PRED7  ozone predictors set7       various
!    REAL arr  SO2MLT  SO2 multiplier              none
!    REAL arr  TRACEGAS_PRED  trace gas pert predictors   various
!    REAL arr  H2O_PRED1  water predictors set1       various
!    REAL arr  H2O_PRED2  water predictors set2       various
!    REAL arr  H2O_PRED3  water predictors set3       various
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
!    May 2008 version of the 100 layer AIRS Fast Transmittance
!    Code by L.L.Strow/S.Hannon.
!
!    Rapid transmittace algorithm predictors consisting of various gas
!    amount and temperature ratios and offsets relative to a reference
!    profile are calculated.
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
!    FIXED_PRED1: FWO (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
!
!    FIXED_PRED2: FOW (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
!
!    FIXED_PRED3: FMW (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
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
!    H2O_PRED1: FWO (11 terms):
!       1) W*a           2) sqrt(W*a)       3) W*a*W/Wz
!       4) W*a*dT        5) (W*a)^2         6) sqrt(W*a)*dT
!       7) root^4(W*a)   8) sqrt(W*a)*W/Wz  9) (W*a)^3
!      10) W            11) W*a*dT*|dT|
!
!    H2O_PRED2: FOW (11 terms):
!       1) W*a             2) sqrt(W*a)     3) W*a*dT
!       4) W*a*Ox*a        5) (W*a)^2       6) root^4(W*a)
!       7) sqrt(W*a)*dT    8) W*a*W/Wz      9) (W*a)^3
!      10) W*a*(Ox*a)^2   11) sqrt(W*a)*W/Wz
!
!    H2O_PRED3: FMW (11 terms):
!       1) W*a             2) sqrt(W*a)     3) W*a*W/Wz
!       4) W*a*dT          5) (W*a)^2       6) sqrt(W*a)*dT
!       7) root^4(W*a)     8) (W*a)^3       9) W
!      10) sqrt(W*a)*W/Wz 11) sqrt(W*a)*Mz*a
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
!    O3_PRED1: FWO (5 terms):
!       1) O*a             2) sqrt(O*a)     3) O*a*dT
!       4) (O*a)^2         5) sqrt(O*a)*dT
!
!    O3_PRED2: FOW (10 terms):
!       1) O*a             2) sqrt(O*a)     3) O*a*dT
!       4) (O*a)^2         5) sqrt(O*a)*dT  6) O*a*O/Ox
!       7) sqrt(O*a)*O/Ox  8) O*a*Oz/Ox     9) O*a*sqrt(Ox*a)
!      10) O*a*TOz*a
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
!    CH4_PRED3: methane predictors (9 terms):
!       1) M*a           2) sqrt(M*a)     3) M*a*dT
!       4) (M*a)^2       5) M*a^2         6) Mz*a
!       7) M*dT          8) TMz*a         9) sqrt(Mz*a)
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
!    "W" is the water amount ratio H2O_PROFILE/REF_H2O_PROFILE
!    "dT" is the temperature offset TEMPERATURE_PROFILE-REF_TEMPERATURE_PROFILE
!    "Wz" is the pressure weighted water amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * H2O_PROFILE(i) },
!      divided by the same sum except using REF_H2O_PROFILE instead of H2O_PROFILE.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!    "O" is the ozone amount ratio O3_PROFILE/REF_O3_PROFILE
!    "Oz" is the pressure weighted ozone amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * O3_PROFILE(i) },
!      divided by the same sum except using REF_O3_PROFILE instead of O3_PROFILE.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2)
!    "Ox" is the unweighted ozone amount above ratio, the
!      sum i=1 to i=L of { O3_PROFILE(i) },
!      divided by the same sum except using REF_O3_PROFILE instead of O3_PROFILE.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!    "TOz" is the pressure and ozone weighted temperature ratio above,
!      sum i=2 to i=L of { P(i) * ( P(i)-P(i-1) )* dT(i-1) * O(i-1) }
!      and TOz(L=1)=0
!    "C" is the carbon monoxide amount ratio O3_PROFILE/REF_O3_PROFILE
!    "Cz" is the pressure weighted CO amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * CO_PROFILE(i) },
!      divided by the same sum except using REF_CO_PROFILE instead of CO_PROFILE.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!    "M" is the methane amount ratio CH4_PROFILE/REF_CH4_PROFILE
!    "Mz" is the pressure weighted methane amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * CH4_PROFILE(i) },
!      divided by the same sum except using REF_CH4_PROFILE instead of CH4_PROFILE.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!    "TMz" is the pressure and methane weighted temperature ratio above,
!      sum i=2 to i=L of { P(i) * ( P(i)-P(i-1) )* Tr(i-1) * M(i-1) }
!      and TMz(L=1)=0
!
!    -----------------------------------------------------
!    FIXMUL: the not-quite-fixed "fixed" amount multiplier.  The
!       value should be close to (within a few percent of) unity.
!       This term adjusts for the effects of water vapor displacement
!       and latitude dependent gravity.  The equations used below are
!       a combination of analytic adjustments for water and gravity,
!       as well as trial-and-error fudge factors to make it all work
!       accurately for any realistic surface pressure and altitude.
!    ===================================================================


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    Assumes the user has supplied vaguely realistic profile amounts
!    and temperatures.


!ROUTINE HISTORY:
! Date        Programmer     Comments
! ----------- -------------- --------------------------------------
!  1 Dec 1994 Scott Hannon   Created
! 10 Apr 1995 Scott Hannon   New header comments; redefined WZ;
!                               changed SECANG to array
!  6 Sep 1995 Scott Hannon   Correct WZ for top layer
!  3 Feb 1997 Scott Hannon   Re-wrote it for FWO+FOW+FMW+FCOW
!  7 Jul 1997 Scott Hannon   Re-wrote it for sets 1 thru 7
! 30 Sep 1997 Scott Hannon   Added CO2PRD
!  5 Mar 1998 Scott Hannon   Deleted water pred 12 & 13 of set 1 & 3
! 26 Aug 1998 Scott Hannon   Add LBOT to call; loop on LBOT instead
!                               of MAXLAY
! 31 Mar 2000 Scott Hannon   Change FIXMUL equation; add ALAT input
!                               var; add FX and REF_LAYER_THICKNESS_PROFILE data; add
!                               PWATER, PMULT, and GSCAL local vars.
! 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
! 17 Aug 2000 Scott Hannon   Add FX & REF_LAYER_THICKNESS_PROFILE input vars
! 12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
! 28 Jun 2005 Scott Hannon   "trace" version for CO2,SO2,HNO3,N2O
! 23 Jan 2008 Scott Hannon   Add LCO2,LN2O,LSO2,LHNO3 switches for
!                               perturbation multiplier calcs; add
!                               LCO2PM to allow CO2 ppmv profile
! 14 May 2008 Scott Hannon   Add no prof CO2MLT calc; add CO2TOP and
!                               CO2PPM to call; add CO2TOP calc
! 30 Nov 2017 Bill Irion     Modified to be F90 compliant
! 22 Dec 2017 Bill Irion	 Put in IF statements to calculate only those
!                            predictors for needed molecules and channel sets.
!                            Made INITIALIZE subroutine to calculate and store
!                            those array elements that don't change by retrieval
!END====================================================================

MODULE CALPAR

	USE INCFTC      
	USE REF_PROFILES
	USE SELECT_SPECIES
	

	IMPLICIT NONE

	! Elements that remain the same retrieval-to-retrieval
	REAL, PRIVATE, DIMENSION(MAXLAY) :: PDP_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: PNORM_ARRAY

	REAL, PRIVATE, DIMENSION(MAXLAY) :: WZREF_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: XZREF_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: OZREF_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: CZREF_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: MZREF_ARRAY

	REAL :: GSCAL

	LOGICAL :: CALPAR_INITIALIZED

	! Elements that are saved ONLY when START_LAYER == 1 and SELECT_PRED == "ALL"
	! These are substituted in when this routine is needed to calculate predictors
	! only below a layer > 1 (as would happen in calculating a Jacobian.)

	REAL, PRIVATE, DIMENSION(MAXLAY) :: TR_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: DT_SAVED_ARRAY

	REAL, PRIVATE, DIMENSION(MAXLAY) :: A_O_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: A_M_SAVED_ARRAY

	REAL, PRIVATE, DIMENSION(MAXLAY) :: CO2TOP_SAVED_ARRAY


	REAL, PRIVATE, DIMENSION(MAXLAY) :: TZ_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: TOZ_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: TMZ_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: WZ_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: XZ_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: OZ_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: CZ_SAVED_ARRAY
	REAL, PRIVATE, DIMENSION(MAXLAY) :: MZ_SAVED_ARRAY

	REAL, TARGET, DIMENSION(N1CON, MAXLAY) :: H2O_CONTINUUM_PRED_SAVED
	REAL, TARGET, DIMENSION(N1FIX, MAXLAY) :: FIXED_PRED1_SAVED
	REAL, TARGET, DIMENSION(N2FIX, MAXLAY) :: FIXED_PRED2_SAVED
	REAL, TARGET, DIMENSION(N3FIX, MAXLAY) :: FIXED_PRED3_SAVED
	REAL, TARGET, DIMENSION(N4FIX, MAXLAY) :: FIXED_PRED4_SAVED
	REAL, TARGET, DIMENSION(N5FIX, MAXLAY) :: FIXED_PRED5_SAVED
	REAL, TARGET, DIMENSION(N6FIX, MAXLAY) :: FIXED_PRED6_SAVED
	REAL, TARGET, DIMENSION(N7FIX, MAXLAY) :: FIXED_PRED7_SAVED
	REAL, TARGET, DIMENSION(N1H2O, MAXLAY) :: H2O_PRED1_SAVED
	REAL, TARGET, DIMENSION(N2H2O, MAXLAY) :: H2O_PRED2_SAVED
	REAL, TARGET, DIMENSION(N3H2O, MAXLAY) :: H2O_PRED3_SAVED
	REAL, TARGET, DIMENSION(N4H2O, MAXLAY) :: H2O_PRED4_SAVED
	REAL, TARGET, DIMENSION(N5H2O, MAXLAY) :: H2O_PRED5_SAVED
	REAL, TARGET, DIMENSION(N6H2O, MAXLAY) :: H2O_PRED6_SAVED
	REAL, TARGET, DIMENSION(N7H2O, MAXLAY) :: H2O_PRED7_SAVED
	REAL, TARGET, DIMENSION( N1O3, MAXLAY) :: O3_PRED1_SAVED
	REAL, TARGET, DIMENSION( N2O3, MAXLAY) :: O3_PRED2_SAVED
	!                        Note: There is no O3_PRED3_SAVED
	REAL, TARGET, DIMENSION( N4O3, MAXLAY) :: O3_PRED4_SAVED
	REAL, TARGET, DIMENSION( N5O3, MAXLAY) :: O3_PRED5_SAVED
	REAL, TARGET, DIMENSION( N6O3, MAXLAY) :: O3_PRED6_SAVED
	REAL, TARGET, DIMENSION( N7O3, MAXLAY) :: O3_PRED7_SAVED
	REAL, TARGET, DIMENSION(N3CH4, MAXLAY) :: CH4_PRED3_SAVED
	REAL, TARGET, DIMENSION( N4CO, MAXLAY) :: CO_PRED4_SAVED
	REAL, TARGET, DIMENSION(NTRACE,MAXLAY) :: TRACEGAS_PRED_SAVED

	! Saved multipliers ONLY when START_LAYER == 1 and SELECT_PRED == "ALL"

	REAL, PRIVATE, DIMENSION(MAXLAY) :: CO2MLT_SAVED
	REAL, PRIVATE, DIMENSION(MAXLAY) :: SO2MLT_SAVED
	REAL, PRIVATE, DIMENSION(MAXLAY) :: HNOMLT_SAVED
	REAL, PRIVATE, DIMENSION(MAXLAY) :: N2OMLT_SAVED
	REAL, PRIVATE, DIMENSION(MAXLAY) :: FIXMUL_SAVED

	REAL, PRIVATE :: CO2TOP_MEAN_SAVED

	! These predictors can change from calls to calculate Jacobians
	REAL, TARGET, DIMENSION(N1CON, MAXLAY) :: H2O_CONTINUUM_PRED
	REAL, TARGET, DIMENSION(N1FIX, MAXLAY) :: FIXED_PRED1
	REAL, TARGET, DIMENSION(N2FIX, MAXLAY) :: FIXED_PRED2
	REAL, TARGET, DIMENSION(N3FIX, MAXLAY) :: FIXED_PRED3
	REAL, TARGET, DIMENSION(N4FIX, MAXLAY) :: FIXED_PRED4
	REAL, TARGET, DIMENSION(N5FIX, MAXLAY) :: FIXED_PRED5
	REAL, TARGET, DIMENSION(N6FIX, MAXLAY) :: FIXED_PRED6
	REAL, TARGET, DIMENSION(N7FIX, MAXLAY) :: FIXED_PRED7
	REAL, TARGET, DIMENSION(N1H2O, MAXLAY) :: H2O_PRED1
	REAL, TARGET, DIMENSION(N2H2O, MAXLAY) :: H2O_PRED2
	REAL, TARGET, DIMENSION(N3H2O, MAXLAY) :: H2O_PRED3
	REAL, TARGET, DIMENSION(N4H2O, MAXLAY) :: H2O_PRED4
	REAL, TARGET, DIMENSION(N5H2O, MAXLAY) :: H2O_PRED5
	REAL, TARGET, DIMENSION(N6H2O, MAXLAY) :: H2O_PRED6
	REAL, TARGET, DIMENSION(N7H2O, MAXLAY) :: H2O_PRED7
	REAL, TARGET, DIMENSION( N1O3, MAXLAY) :: O3_PRED1
	REAL, TARGET, DIMENSION( N2O3, MAXLAY) :: O3_PRED2
	!                        Note: There is no O3_PRED3
	REAL, TARGET, DIMENSION( N4O3, MAXLAY) :: O3_PRED4
	REAL, TARGET, DIMENSION( N5O3, MAXLAY) :: O3_PRED5
	REAL, TARGET, DIMENSION( N6O3, MAXLAY) :: O3_PRED6
	REAL, TARGET, DIMENSION( N7O3, MAXLAY) :: O3_PRED7
	REAL, TARGET, DIMENSION(N3CH4, MAXLAY) :: CH4_PRED3
	REAL, TARGET, DIMENSION( N4CO, MAXLAY) :: CO_PRED4
	REAL, TARGET, DIMENSION(NTRACE,MAXLAY) :: TRACEGAS_PRED

CONTAINS

	! Initialize those elements that do not change from retrieval to retrieval
	SUBROUTINE CALPAR__INITIALIZE

		INTEGER :: L

		PDP_ARRAY(1) = REF_LAYER_PRESSURE_PROFILE(1) * (REF_LAYER_PRESSURE_PROFILE(2) - REF_LAYER_PRESSURE_PROFILE(1))
		PNORM_ARRAY(1) = 0.
		WZREF_ARRAY(1) = PDP_ARRAY(1) + REF_H2O_PROFILE(1)
		XZREF_ARRAY(1) = REF_O3_PROFILE(1)
		OZREF_ARRAY(1) = PDP_ARRAY(1) * REF_O3_PROFILE(1)
		CZREF_ARRAY(1) = PDP_ARRAY(1) * REF_CO_PROFILE(1)
		MZREF_ARRAY(1) = PDP_ARRAY(1) * REF_CH4_PROFILE(1)
		DO L = 2, MAXLAY
			PDP_ARRAY(L) = REF_LAYER_PRESSURE_PROFILE(L)*( REF_LAYER_PRESSURE_PROFILE(L) - REF_LAYER_PRESSURE_PROFILE(L-1) )
			PNORM_ARRAY(L) = PNORM_ARRAY(L-1) + PDP_ARRAY(L)
			WZREF_ARRAY(L) = WZREF_ARRAY(L-1) + PDP_ARRAY(L) * REF_H2O_PROFILE(L)
			XZREF_ARRAY(L) = XZREF_ARRAY(L-1) + REF_O3_PROFILE(L)
			OZREF_ARRAY(L) = OZREF_ARRAY(L-1) + PDP_ARRAY(L) * REF_O3_PROFILE(L)
			CZREF_ARRAY(L) = CZREF_ARRAY(L-1) + PDP_ARRAY(L) * REF_CO_PROFILE(L)
			MZREF_ARRAY(L) = MZREF_ARRAY(L-1) + PDP_ARRAY(L) * REF_CH4_PROFILE(L)
		ENDDO

		CALPAR_INITIALIZED = .TRUE.
		
		RETURN

	END SUBROUTINE CALPAR__INITIALIZE

	SUBROUTINE CALPAR__CALCULATE_PREDICTORS ( LBOT, &
		TEMPERATURE_PROFILE, CO2_PROFILE, H2O_PROFILE, O3_PROFILE, &
		CO_PROFILE, CH4_PROFILE, SO2_PROFILE, HNO3_PROFILE, N2O_PROFILE, &
		SECANG,  ALAT,    FX, &
		LCO2,  LN2O,  LSO2, LHNO3, LCO2PM, CO2PPM, CO2TOP, &
		FIXMUL, H2O_CONTINUUM_PRED_PTR, &
		FIXED_PRED1_PTR, FIXED_PRED2_PTR, FIXED_PRED3_PTR, FIXED_PRED4_PTR, &
		FIXED_PRED5_PTR, FIXED_PRED6_PTR, FIXED_PRED7_PTR,  &
		H2O_PRED1_PTR,  H2O_PRED2_PTR,  H2O_PRED3_PTR,  H2O_PRED4_PTR, &
		H2O_PRED5_PTR,  H2O_PRED6_PTR,  H2O_PRED7_PTR,  &
		O3_PRED1_PTR,   O3_PRED2_PTR,                   O3_PRED4_PTR, &
		O3_PRED5_PTR,   O3_PRED6_PTR,   O3_PRED7_PTR,  &
		CH4_PRED3_PTR,  CO_PRED4_PTR,   TRACEGAS_PRED_PTR, &
		CO2MLT, SO2MLT, HNOMLT, N2OMLT, &
		START_LAYER, DO_SPECIES)

		USE AIRS_SARTA_VARIABLES, ONLY: NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7

		!
		! ARGUMENTS
		!
		! Input
		INTEGER, INTENT(IN) :: LBOT
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: TEMPERATURE_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: CO2_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: H2O_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: O3_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: CO_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: CH4_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: SO2_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: HNO3_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: N2O_PROFILE
		!REAL, INTENT(IN), DIMENSION(MAXLAY) ::   REF_LAYER_PRESSURE_PROFILE
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: SECANG
		REAL, INTENT(IN) :: ALAT
		REAL, INTENT(IN), DIMENSION(MAXLAY) ::     FX
		!REAL, INTENT(IN), DIMENSION(MAXLAY) ::  REF_LAYER_THICKNESS_PROFILE
		LOGICAL, INTENT(IN) :: LCO2
		LOGICAL, INTENT(IN) :: LN2O
		LOGICAL, INTENT(IN) :: LSO2
		LOGICAL, INTENT(IN) :: LHNO3
		LOGICAL, INTENT(IN) :: LCO2PM
		REAL, INTENT(IN) :: CO2PPM

		INTEGER, INTENT(IN) :: START_LAYER
		TYPE (DO_SPECIES_LOGICAL_ARRAY_TYPE), INTENT(IN) :: DO_SPECIES

		! Output

		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: H2O_CONTINUUM_PRED_PTR ! water continuum predictors
		REAL, DIMENSION(:, :), POINTER :: FIXED_PRED1_PTR ! set1 "fixed" predictors
		REAL, DIMENSION(:, :), POINTER :: FIXED_PRED2_PTR ! set2 "fixed" predictors
		REAL, DIMENSION(:, :), POINTER :: FIXED_PRED3_PTR ! set3 "fixed" predictors
		REAL, DIMENSION(:, :), POINTER :: FIXED_PRED4_PTR ! set4 "fixed" predictors
		REAL, DIMENSION(:, :), POINTER :: FIXED_PRED5_PTR ! set5 "fixed" predictors
		REAL, DIMENSION(:, :), POINTER :: FIXED_PRED6_PTR ! set6 "fixed" predictors
		REAL, DIMENSION(:, :), POINTER :: FIXED_PRED7_PTR ! set7 "fixed" predictors
		REAL, DIMENSION(:, :), POINTER :: H2O_PRED1_PTR ! set1 water predictors
		REAL, DIMENSION(:, :), POINTER :: H2O_PRED2_PTR ! set2 water predictors
		REAL, DIMENSION(:, :), POINTER :: H2O_PRED3_PTR ! set3 water predictors
		REAL, DIMENSION(:, :), POINTER :: H2O_PRED4_PTR ! set4 water predictors
		REAL, DIMENSION(:, :), POINTER :: H2O_PRED5_PTR ! set5 water predictors
		REAL, DIMENSION(:, :), POINTER :: H2O_PRED6_PTR ! set6 water predictors
		REAL, DIMENSION(:, :), POINTER :: H2O_PRED7_PTR ! set7 water predictors
		REAL, DIMENSION(:, :), POINTER :: O3_PRED1_PTR ! set1 ozone predictors
		REAL, DIMENSION(:, :), POINTER :: O3_PRED2_PTR ! set2 ozone predictors
		REAL, DIMENSION(:, :), POINTER :: O3_PRED4_PTR ! set4 ozone predictors
		REAL, DIMENSION(:, :), POINTER :: O3_PRED5_PTR ! set5 ozone predictors
		REAL, DIMENSION(:, :), POINTER :: O3_PRED6_PTR ! set6 ozone predictors
		REAL, DIMENSION(:, :), POINTER :: O3_PRED7_PTR ! set7 ozone predictors
		REAL, DIMENSION(:, :), POINTER :: CH4_PRED3_PTR ! set3 methane predictors
		REAL, DIMENSION(:, :), POINTER :: CO_PRED4_PTR ! set4 carbon monoxide predictors
		REAL, DIMENSION(:, :), POINTER :: TRACEGAS_PRED_PTR ! trace gas pert perdictors

		REAL, INTENT(OUT) :: CO2TOP
		REAL, INTENT(OUT), DIMENSION(MAXLAY) ::  FIXMUL
!		REAL, INTENT(OUT), DIMENSION(N1CON, MAXLAY) :: H2O_CONTINUUM_PRED
!		REAL, INTENT(OUT), DIMENSION(:, :), POINTER :: FIXED_PRED1_PTR
!		REAL, INTENT(OUT), DIMENSION(N2FIX, MAXLAY) :: FIXED_PRED2
!		REAL, INTENT(OUT), DIMENSION(N3FIX, MAXLAY) :: FIXED_PRED3
!		REAL, INTENT(OUT), DIMENSION(N4FIX, MAXLAY) :: FIXED_PRED4
!		REAL, INTENT(OUT), DIMENSION(N5FIX, MAXLAY) :: FIXED_PRED5
!		REAL, INTENT(OUT), DIMENSION(N6FIX, MAXLAY) :: FIXED_PRED6
!		REAL, INTENT(OUT), DIMENSION(N7FIX, MAXLAY) :: FIXED_PRED7
!		REAL, INTENT(OUT), DIMENSION(N1H2O, MAXLAY) :: H2O_PRED1
!		REAL, INTENT(OUT), DIMENSION(N2H2O, MAXLAY) :: H2O_PRED2
!		REAL, INTENT(OUT), DIMENSION(N3H2O, MAXLAY) :: H2O_PRED3
!		REAL, INTENT(OUT), DIMENSION(N4H2O, MAXLAY) :: H2O_PRED4
!		REAL, INTENT(OUT), DIMENSION(N5H2O, MAXLAY) :: H2O_PRED5
!		REAL, INTENT(OUT), DIMENSION(N6H2O, MAXLAY) :: H2O_PRED6
!		REAL, INTENT(OUT), DIMENSION(N7H2O, MAXLAY) :: H2O_PRED7
!		REAL, INTENT(OUT), DIMENSION( N1O3, MAXLAY) :: O3_PRED1
!		REAL, INTENT(OUT), DIMENSION( N2O3, MAXLAY) :: O3_PRED2
!		!                            Note: There is no O3_PRED3
!		REAL, INTENT(OUT), DIMENSION( N4O3, MAXLAY) :: O3_PRED4
!		REAL, INTENT(OUT), DIMENSION( N5O3, MAXLAY) :: O3_PRED5
!		REAL, INTENT(OUT), DIMENSION( N6O3, MAXLAY) :: O3_PRED6
!		REAL, INTENT(OUT), DIMENSION( N7O3, MAXLAY) :: O3_PRED7
!		REAL, INTENT(OUT), DIMENSION(N3CH4, MAXLAY) :: CH4_PRED3
!		REAL, INTENT(OUT), DIMENSION( N4CO, MAXLAY) :: CO_PRED4
!		REAL, INTENT(OUT), DIMENSION(NTRACE,MAXLAY) :: TRACEGAS_PRED
		REAL, INTENT(OUT), DIMENSION(MAXLAY) :: CO2MLT
		REAL, INTENT(OUT), DIMENSION(MAXLAY) :: SO2MLT
		REAL, INTENT(OUT), DIMENSION(MAXLAY) :: HNOMLT
		REAL, INTENT(OUT), DIMENSION(MAXLAY) :: N2OMLT


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
		REAL    A_F
		REAL    A_W
		REAL  WZREF
		REAL     WZ
		REAL   AZ_W
		REAL    A_O
		REAL  XZREF
		REAL     XZ
		REAL   XZ_O
		REAL  OZREF
		REAL     OZ
		REAL   AZ_O
		REAL    TOZ
		REAL  TAZ_O
		REAL    A_C
		REAL     CZ
		REAL  CZREF
		REAL   AZ_C
		REAL    A_M
		REAL  MZREF
		REAL     MZ
		REAL   AZ_M
		REAL    TMZ
		REAL  TAZ_M
		REAL TR_SQUARED
		REAL WJUNKA
		REAL WJUNKR
		REAL WJUNKS
		REAL WJUNKZ
		REAL WJUNK4
		REAL OJUNKA
		REAL OJUNKR
		REAL OJUNKZ
		REAL OJUNKX
		REAL CO_JUNK_A
		REAL CO_JUNK_R
		REAL CO_JUNK_S
		REAL CO_JUNK_Z
		REAL MJUNKA
		REAL MJUNKR
		REAL MJUNKZ

		! Variables for fixed gases adjustment
		REAL PWATER
		!REAL  GSCAL

		! Variables with fixed assignments
		! Changed from DATA statements by FWI 11/30/2017
		REAL, PARAMETER :: PMULT = 0.58 		! fudge factor * (0.622=M_H2O/M_AIR)	
		REAL, PARAMETER :: STDDEN = 2.6867E+19	! Loschmidt aka standard density
		REAL, PARAMETER :: STDTMP = 273.15		! Standard Temperature
		REAL, PARAMETER :: KMOLE = 6.022045E+26	! 1000 * Avagadro's Number

		INTEGER :: INIT_LAYER


		!
		!
		! EXECUTABLE CODE
		!
		!


		IF (DO_SPECIES % DO_NONE) THEN ! Return all saved predictors and multipliers without making any new calculations

			H2O_CONTINUUM_PRED_PTR => H2O_CONTINUUM_PRED_SAVED
			TRACEGAS_PRED_PTR => TRACEGAS_PRED_SAVED		

			IF (NCHN1 .GT. 0) THEN
				FIXED_PRED1_PTR => FIXED_PRED1_SAVED
				H2O_PRED1_PTR => H2O_PRED1_SAVED
				O3_PRED1_PTR => O3_PRED1_SAVED
			ENDIF
			IF (NCHN2 .GT. 0) THEN
				FIXED_PRED2_PTR => FIXED_PRED2_SAVED
				H2O_PRED2_PTR => H2O_PRED2_SAVED
				O3_PRED2_PTR => O3_PRED2_SAVED
			ENDIF
			IF (NCHN3 .GT. 0) THEN
				FIXED_PRED3_PTR => FIXED_PRED3_SAVED
				H2O_PRED3_PTR => H2O_PRED3_SAVED
				CH4_PRED3_PTR => CH4_PRED3_SAVED
			ENDIF
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

			CO2MLT(:) = CO2MLT_SAVED(:)
			SO2MLT(:) = SO2MLT_SAVED(:)
			HNOMLT(:) = HNOMLT_SAVED(:)
			N2OMLT(:) = N2OMLT_SAVED(:)			
			FIXMUL(:) = FIXMUL_SAVED(:)			

			CO2TOP = CO2TOP_MEAN_SAVED

			RETURN
		ENDIF

		
		! Calc the fixed gases gravity correction factor
		IF (DO_SPECIES % DO_ALL) THEN
			GSCAL=( 9.78050518 + 0.0518017*( COS( (ALAT - 90.0)*PI/180) )**2 )/9.80683613 
		ENDIF

		IF (START_LAYER .EQ. 1) THEN
			! Initialize the sum terms to zero
			!PNORM = 0.E+0
			TZ = 0.0E+0
			!WZREF = 0.0E+0
			WZ = 0.0E+0
			!XZREF = 0.0E+0
			XZ = 0.0E+0
			!OZREF = 0.0E+0
			OZ = 0.0E+0
			TOZ = 0.0E+0
			!CZREF = 0.0E+0
			CZ = 0.0E+0
			!MZREF = 0.0E+0
			MZ = 0.0E+0
			TMZ = 0.0E+0
			CO2TOP = 0.0E+0
		ELSE
			INIT_LAYER = START_LAYER - 1
			! Pull the values from saved arrays
			TZ = TZ_SAVED_ARRAY(INIT_LAYER)
			WZ = WZ_SAVED_ARRAY(INIT_LAYER)
			XZ = XZ_SAVED_ARRAY(INIT_LAYER)
			OZ = OZ_SAVED_ARRAY(INIT_LAYER)
			TOZ = TOZ_SAVED_ARRAY(INIT_LAYER)
			CZ = CZ_SAVED_ARRAY(INIT_LAYER)
			MZ = MZ_SAVED_ARRAY(INIT_LAYER)
			TMZ = TMZ_SAVED_ARRAY(INIT_LAYER)
			CO2TOP = CO2TOP_SAVED_ARRAY(INIT_LAYER)

			TR = TR_SAVED_ARRAY(INIT_LAYER)
			DT = DT_SAVED_ARRAY(INIT_LAYER)
			A_O = A_O_SAVED_ARRAY(INIT_LAYER)
			A_M = A_M_SAVED_ARRAY(INIT_LAYER)

			! Restore saved predictors above start layer
			! For now, restore all layers. Later on, we can be more selective
			! depending on what SELECT_PRED is

			FIXMUL(1:INIT_LAYER) = FIXMUL_SAVED(1:INIT_LAYER)

			IF (DO_SPECIES % DO_H2O) H2O_CONTINUUM_PRED(:, :) = H2O_CONTINUUM_PRED_SAVED(:, :)

			IF (DO_SPECIES % DO_TEMPERATURE) TRACEGAS_PRED(:, :) = TRACEGAS_PRED_SAVED(:, :)
			
			IF (NCHN1 .GT. 0) THEN
				IF (DO_SPECIES % DO_TEMPERATURE) FIXED_PRED1(:, 1:INIT_LAYER) = FIXED_PRED1_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_H2O) H2O_PRED1(:, 1:INIT_LAYER) = H2O_PRED1_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_O3) O3_PRED1(:, 1:INIT_LAYER) = O3_PRED1_SAVED(:, 1:INIT_LAYER)
			ENDIF

			IF (NCHN2 .GT. 0) THEN
				IF (DO_SPECIES % DO_TEMPERATURE) FIXED_PRED2(:, 1:INIT_LAYER) = FIXED_PRED2_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_H2O) H2O_PRED2(:, 1:INIT_LAYER) = H2O_PRED2_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_O3) O3_PRED2(:, 1:INIT_LAYER) = O3_PRED2_SAVED(:, 1:INIT_LAYER)
			ENDIF

			IF (NCHN3 .GT. 0) THEN
				IF (DO_SPECIES % DO_TEMPERATURE) FIXED_PRED3(:, 1:INIT_LAYER) = FIXED_PRED3_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_H2O) H2O_PRED3(:, 1:INIT_LAYER) = H2O_PRED3_SAVED(:, 1:INIT_LAYER)
				IF (DO_SPECIES % DO_CH4) CH4_PRED3(:, 1:INIT_LAYER) = CH4_PRED3_SAVED(:, 1:INIT_LAYER)
			ENDIF

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

		!--------------------
		! Loop over the layers
		!--------------------
		DO L = START_LAYER, LBOT

			!
			! Calculate the basic profile
			! dependent predictors.
			!
			IF (L .EQ. 1) THEN
				!PDP = REF_LAYER_PRESSURE_PROFILE(1) * (REF_LAYER_PRESSURE_PROFILE(2) - REF_LAYER_PRESSURE_PROFILE(1))
				PDP = PDP_ARRAY(1)
				TRZ = 0.0E+0
				TAZ_O = 0.0E+0
				TAZ_M = 0.0E+0
			ELSE
				!PDP = REF_LAYER_PRESSURE_PROFILE(L)*( REF_LAYER_PRESSURE_PROFILE(L) - REF_LAYER_PRESSURE_PROFILE(L-1) )
				!PNORM = PNORM + PDP
				PDP = PDP_ARRAY(L)
				PNORM = PNORM_ARRAY(L)
			
				! Note: TRZ, TOZ, and TMZ use layer-above terms
				TZ = TZ + PDP*TR  ! Eq 20 in Strow et al., 2003
				TRZ = TZ/PNORM
			
				TOZ = TOZ + PDP*DT*A_O
				TAZ_O = TOZ/PNORM
			
				TMZ = TMZ + PDP*TR*A_M
				TAZ_M = TMZ/PNORM
			ENDIF

			!IF (START_LAYER .NE. 1)
				! Pull WZ, XZ, OZ, CZ, MZ, TZ, TOZ and TMZ(L-1) from array
			!ENDIF

			! Temperature terms
			DT=TEMPERATURE_PROFILE(L) - REF_TEMPERATURE_PROFILE(L)
			TR=TEMPERATURE_PROFILE(L)/REF_TEMPERATURE_PROFILE(L)

			! Calc the fixed gases correction term for this layer
			PWATER=KMOLE*H2O_PROFILE(L)*TEMPERATURE_PROFILE(L)/(STDDEN*STDTMP*100*REF_LAYER_THICKNESS_PROFILE(L))
			A_F=( 1 - PMULT*PWATER/REF_LAYER_PRESSURE_PROFILE(L) )/( FX(L)*GSCAL )
			FIXMUL(L)=A_F  ! Moved from outside IF(DO_SPECIES % DO_TEMPERTURE) block  FWI 1/8/17

			! for testing
			! A_F=1.0

			! Water terms
			A_W=H2O_PROFILE(L)/REF_H2O_PROFILE(L)
			!WZREF=WZREF + PDP*REF_H2O_PROFILE(L)
			WZ=WZ + PDP*H2O_PROFILE(L)
			!AZ_W=WZ/WZREF
			AZ_W=WZ/WZREF_ARRAY(L)

			! Ozone terms
			A_O=O3_PROFILE(L)/REF_O3_PROFILE(L)
			!XZREF=XZREF + REF_O3_PROFILE(L)
			XZ=XZ + O3_PROFILE(L)
			!XZ_O=XZ/XZREF
			XZ_O=XZ/XZREF_ARRAY(L)
			!OZREF=OZREF + PDP*REF_O3_PROFILE(L)
			OZ=OZ + PDP*O3_PROFILE(L)
			!AZ_O=OZ/OZREF
			AZ_O=OZ/OZREF_ARRAY(L)

			! Carbon monoxide terms
			A_C=CO_PROFILE(L)/REF_CO_PROFILE(L)
			!CZREF=CZREF + PDP*REF_CO_PROFILE(L)
			CZ=CZ + PDP*CO_PROFILE(L)
			!AZ_C=CZ/CZREF
			AZ_C=CZ/CZREF_ARRAY(L)

			! Methane terms
			A_M=CH4_PROFILE(L)/REF_CH4_PROFILE(L)
			!MZREF=MZREF + PDP*REF_CH4_PROFILE(L)
			MZ=MZ + PDP*CH4_PROFILE(L)
			!AZ_M=MZ/MZREF
			AZ_M=MZ/MZREF_ARRAY(L)

			IF ((START_LAYER .EQ. 1) .AND. DO_SPECIES % DO_ALL) THEN
				! Save above layer terms in array
				! Save  WZ, XZ, OZ, CZ, MZ, TZ, TOZ and TMZ(L) 
				TR_SAVED_ARRAY(L) = TR
				DT_SAVED_ARRAY(L) = DT

				A_O_SAVED_ARRAY(L) = A_O
				A_M_SAVED_ARRAY(L) = A_M

				TZ_SAVED_ARRAY(L) = TZ
				TOZ_SAVED_ARRAY(L) = TOZ
				TMZ_SAVED_ARRAY(L) = TMZ
				WZ_SAVED_ARRAY(L) = WZ
				XZ_SAVED_ARRAY(L) = XZ
				OZ_SAVED_ARRAY(L) = OZ
				CZ_SAVED_ARRAY(L) = CZ
				MZ_SAVED_ARRAY(L) = MZ
			ENDIF

			!----------------------
			! Load up the predictors
			!----------------------
			IF (DO_SPECIES % DO_TEMPERATURE) THEN
				!-----
				! Fixed (for FWO, FOW, FMW, & FCOW)
				!-----
				TR_SQUARED = TR * TR
				!FIXMUL(L)=A_F  ! Moved outside of block   FWI 1/8/17

				! Calculate Fixed predictors for set 1 even if they're not needed
				! as they may be copied for sets 2 or 3
				IF ( (NCHN1 .GT. 0) .OR. (NCHN2 .GT. 0) .OR. (NCHN3 .GT. 0)) THEN
					FIXED_PRED1(1,L)=SECANG(L)
					FIXED_PRED1(2,L)=SECANG(L)*SECANG(L)
					FIXED_PRED1(3,L)=SECANG(L)*TR
					FIXED_PRED1(4,L)=SECANG(L)*TR_SQUARED
					FIXED_PRED1(5,L)=TR
					FIXED_PRED1(6,L)=TR_SQUARED
					FIXED_PRED1(7,L)=SECANG(L)*TRZ
					FIXED_PRED1(8,L)=SECANG(L)*TRZ/TR
				ENDIF

				IF (NCHN2 .GT. 0) THEN
					FIXED_PRED2(:,L)=FIXED_PRED1(:,L)
				ENDIF

				IF (NCHN3 .GT. 0) THEN
					FIXED_PRED3(:,L)=FIXED_PRED1(:,L)
				ENDIF

				IF (NCHN4 .GT. 0) THEN
					FIXED_PRED4(1,L)=SECANG(L)
					FIXED_PRED4(2,L)=SECANG(L)*SECANG(L)
					FIXED_PRED4(3,L)=SECANG(L)*TR
					FIXED_PRED4(4,L)=SECANG(L)*TR_SQUARED
					FIXED_PRED4(5,L)=TR
					FIXED_PRED4(6,L)=TR_SQUARED
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
					FIXED_PRED5(4,L)=SECANG(L)*TR_SQUARED
					FIXED_PRED5(5,L)=TR
					FIXED_PRED5(6,L)=TR_SQUARED
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
					FIXED_PRED6(4,L)=SECANG(L)*TR_SQUARED
					FIXED_PRED6(5,L)=TR
					FIXED_PRED6(6,L)=TR_SQUARED
					FIXED_PRED6(7,L)=SECANG(L)*TRZ
					FIXED_PRED6(8,L)=SQRT(SECANG(L))
				ENDIF

				! Fixed predictors for FWO sun mfbw = set7
				IF (NCHN7 .GT. 0) THEN
					FIXED_PRED7(1,L)=SECANG(L)
					FIXED_PRED7(2,L)=SECANG(L)*SECANG(L)
					FIXED_PRED7(3,L)=SECANG(L)*TR
					FIXED_PRED7(4,L)=SECANG(L)*TR_SQUARED
					FIXED_PRED7(5,L)=TR
					FIXED_PRED7(6,L)=TR_SQUARED
					FIXED_PRED7(7,L)=SECANG(L)*TRZ
					FIXED_PRED7(8,L)=SQRT(SECANG(L))
				ENDIF

			ENDIF ! DO_SPECIES % DO_TEMPERATURE
		
			!
			! Ozone
			!
			IF (DO_SPECIES % DO_O3) THEN
				OJUNKA=SECANG(L)*A_O
				OJUNKR=SQRT( OJUNKA )
				OJUNKZ=OJUNKA/XZ_O
				OJUNKX=SECANG(L)*XZ_O

				! Ozone predictors for FWO = set1
				IF (NCHN1 .GT. 0) THEN
					O3_PRED1(1,L)=OJUNKA
					O3_PRED1(2,L)=OJUNKR
					O3_PRED1(3,L)=OJUNKA*DT
					O3_PRED1(4,L)=OJUNKA*OJUNKA
					O3_PRED1(5,L)=OJUNKR*DT
				ENDIF

				! Ozone predictors for FOW = set2
				IF (NCHN2 .GT. 0) THEN
					O3_PRED2( 1,L)=OJUNKA
					O3_PRED2( 2,L)=OJUNKR
					O3_PRED2( 3,L)=OJUNKA*DT
					O3_PRED2( 4,L)=OJUNKA*OJUNKA
					O3_PRED2( 5,L)=OJUNKR*DT
					O3_PRED2( 6,L)=OJUNKZ*A_O
					O3_PRED2( 7,L)=OJUNKR*A_O/XZ_O
					O3_PRED2( 8,L)=OJUNKZ*AZ_O
					O3_PRED2( 9,L)=OJUNKA*SQRT( OJUNKX )
					O3_PRED2(10,L)=OJUNKA*TAZ_O*SECANG(L)
				ENDIF

				! There are no ozone predictors for set3 = FMW (the ozone
				! absorption in the region covered by FMW is negligible).

				! Ozone predictors for FCOW = set4
				IF (NCHN4 .GT. 0) THEN
					O3_PRED4(1,L)=OJUNKA
					O3_PRED4(2,L)=OJUNKR
					O3_PRED4(3,L)=OJUNKA*DT
				ENDIF

				! Ozone predictors for FWO sun bfsw = set5
				O3_PRED5(1,L)=OJUNKA

				! Ozone predictors for FWO sun mfmw = set6
				O3_PRED6(1,L)=OJUNKA

				! Ozone predictors for FWO sun mfbw = set7
				O3_PRED7(1,L)=OJUNKA
			ENDIF ! DO_SPECIES % DO_O3 .OR. DO_SPECIES % DO_TEMPERATURE

			!
			! Methane for FMW = set3
			!
			IF (DO_SPECIES % DO_CH4) THEN
				IF (NCHN3 .GT. 0) THEN
					MJUNKA=SECANG(L)*A_M
					MJUNKR=SQRT(MJUNKA)
					MJUNKZ=SECANG(L)*AZ_M
					CH4_PRED3(1,L)=MJUNKA
					CH4_PRED3(2,L)=MJUNKR
					CH4_PRED3(3,L)=MJUNKA*DT
					CH4_PRED3(4,L)=MJUNKA*MJUNKA
					CH4_PRED3(5,L)=MJUNKA*SECANG(L)
					CH4_PRED3(6,L)=MJUNKZ
					CH4_PRED3(7,L)=A_M*DT
					CH4_PRED3(8,L)=TAZ_M*SECANG(L)
					CH4_PRED3(9,L)=SQRT( MJUNKZ )
				ENDIF
			ENDIF

			!
			! Water
			! 
			! Depends on temperature, water
			IF (DO_SPECIES % DO_H2O) THEN
				WJUNKA=SECANG(L)*A_W
				WJUNKR=SQRT( WJUNKA )
				WJUNKS=WJUNKA*WJUNKA
				WJUNKZ=WJUNKA*A_W/AZ_W
				WJUNK4=SQRT( WJUNKR )

				! Water predictors for FWO = set1
				IF (NCHN1 .GT. 0) THEN
					H2O_PRED1( 1,L)=WJUNKA
					H2O_PRED1( 2,L)=WJUNKR
					H2O_PRED1( 3,L)=WJUNKZ
					H2O_PRED1( 4,L)=WJUNKA*DT
					H2O_PRED1( 5,L)=WJUNKS
					H2O_PRED1( 6,L)=WJUNKR*DT
					H2O_PRED1( 7,L)=WJUNK4
					H2O_PRED1( 8,L)=WJUNKZ/WJUNKR
					H2O_PRED1( 9,L)=WJUNKS*WJUNKA
					H2O_PRED1(10,L)=A_W
					H2O_PRED1(11,L)=WJUNKA*DT*ABS( DT )
				ENDIF

				! Water predictors for FOW = set2
				IF (NCHN2 .GT. 0) THEN
					H2O_PRED2( 1,L)=WJUNKA
					H2O_PRED2( 2,L)=WJUNKR
					H2O_PRED2( 3,L)=WJUNKA*DT
					H2O_PRED2( 4,L)=WJUNKA*OJUNKX
					H2O_PRED2( 5,L)=WJUNKS
					H2O_PRED2( 6,L)=WJUNK4
					H2O_PRED2( 7,L)=WJUNKR*DT
					H2O_PRED2( 8,L)=WJUNKZ
					H2O_PRED2( 9,L)=WJUNKA*WJUNKS
					H2O_PRED2(10,L)=WJUNKA*OJUNKX*OJUNKX
					H2O_PRED2(11,L)=WJUNKZ/WJUNKR
				ENDIF

				! Water predictors for FMW = set3
				IF (NCHN3 .GT. 0) THEN
					H2O_PRED3( 1,L)=WJUNKA
					H2O_PRED3( 2,L)=WJUNKR
					H2O_PRED3( 3,L)=WJUNKZ
					H2O_PRED3( 4,L)=WJUNKA*DT
					H2O_PRED3( 5,L)=WJUNKS
					H2O_PRED3( 6,L)=WJUNKR*DT
					H2O_PRED3( 7,L)=WJUNK4
					H2O_PRED3( 8,L)=WJUNKS*WJUNKA
					H2O_PRED3( 9,L)=A_W
					H2O_PRED3(10,L)=WJUNKZ/WJUNKR
					H2O_PRED3(11,L)=WJUNKR*MJUNKZ
				ENDIF

				! Water predictors for FCOW = set4
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
				! Water continuum (for FWO, FOW, FMW, FCOW)
				!
				TR_SQUARED = TR * TR
				H2O_CONTINUUM_PRED(1,L)=WJUNKA/TR_SQUARED
				H2O_CONTINUUM_PRED(2,L)=H2O_CONTINUUM_PRED(1,L)*A_W/TR_SQUARED
				H2O_CONTINUUM_PRED(3,L)=WJUNKA/TR
				H2O_CONTINUUM_PRED(4,L)=H2O_CONTINUUM_PRED(3,L)*A_W
				H2O_CONTINUUM_PRED(5,L)=H2O_CONTINUUM_PRED(1,L)*A_W
				H2O_CONTINUUM_PRED(6,L)=H2O_CONTINUUM_PRED(1,L)/TR_SQUARED
				H2O_CONTINUUM_PRED(7,L)=WJUNKA
			ENDIF

			! 
			! Carbon monoxide for FCOW = set4
			! 
			! Depends on temperature
			IF (DO_SPECIES % DO_CO) THEN
				IF (NCHN4 .GT. 0) THEN
					CO_JUNK_A=SECANG(L)*A_C
					CO_JUNK_R=SQRT( CO_JUNK_A )
					CO_JUNK_S=CO_JUNK_A*CO_JUNK_A
					CO_JUNK_Z=CO_JUNK_A*A_C/AZ_C
					CO_PRED4(1,L)=CO_JUNK_A
					CO_PRED4(2,L)=CO_JUNK_R
					CO_PRED4(3,L)=CO_JUNK_A*DT
					CO_PRED4(4,L)=CO_JUNK_S
					CO_PRED4(5,L)=CO_JUNK_Z
					CO_PRED4(6,L)=CO_JUNK_R*DT
					CO_PRED4(7,L)=SQRT( CO_JUNK_R )
					CO_PRED4(8,L)=CO_JUNK_Z/CO_JUNK_R
					CO_PRED4(9,L)=A_C
					CO_PRED4(10,L)=CO_JUNK_A*SECANG(L)
					CO_PRED4(11,L)=CO_JUNK_R*SECANG(L)
				ENDIF
			ENDIF

			!
			! Trace gas perturbation predictors
			! 
			! The first 4 trace predictors are used by all trace gases
			IF (DO_SPECIES % DO_TEMPERATURE) THEN
				TRACEGAS_PRED(1,L)=SECANG(L)
				TRACEGAS_PRED(2,L)=TR
				TRACEGAS_PRED(3,L)=SECANG(L)*TR
				TRACEGAS_PRED(4,L)=SECANG(L)*TR_SQUARED
				! The last 3 trace predictors are only used by N2O
				TRACEGAS_PRED(5,L)=SECANG(L)*SECANG(L)
				TRACEGAS_PRED(6,L)=1.0
				TRACEGAS_PRED(7,L)=SQRT( SECANG(L) )
			ENDIF

			IF (DO_SPECIES % DO_CO2) THEN
				IF (LCO2) THEN
					IF (LCO2PM) THEN
						CO2MLT(L)=100.0*(CO2_PROFILE(L) - CO2STD)/(3.0*CO2STD)
					ELSE
						! CO2 mult=1 when prof amount = 1.03 * ref amount
						CO2MLT(L)=33.3333*( CO2_PROFILE(L) - FIXMUL(L)*REF_CO2_PROFILE(L) ) /  REF_CO2_PROFILE(L)
						! Ignore changes in CO2 of less than ~0.03%
						IF (ABS(CO2MLT(L)) .LT. 1E-2) CO2MLT(L)=0.0
					ENDIF
				ELSE
					CO2MLT(L)=100.0*(CO2PPM - CO2STD)/(3.0*CO2STD)
				ENDIF
				IF (L .LE. NTEBOT) THEN
					CO2TOP=CO2TOP + CO2STD*(1.0 + CO2MLT(L)*3.0E-2)
				ENDIF
			ELSE
				CO2MLT(L) = CO2MLT_SAVED(L)
				CO2TOP = CO2TOP_SAVED_ARRAY(L)
			ENDIF

			IF (DO_SPECIES % DO_N2O) THEN
				IF (LN2O) THEN
					! N2O mult=-1 when prof amount = 0.75 * ref amount
					N2OMLT(L)=4.0*( N2O_PROFILE(L) - FIXMUL(L)*REF_N2O_PROFILE(L) )/ REF_N2O_PROFILE(L)
					! Ignore changes in N2O less than ~0.3%
					IF (ABS(N2OMLT(L)) .LT. 1E-2) N2OMLT(L)=0.0
				ELSE
					N2OMLT(L)=0.0
				ENDIF
			ELSE
				N2OMLT(L) = N2OMLT_SAVED(L)
			ENDIF

			IF (DO_SPECIES % DO_SO2) THEN
				IF (LSO2) THEN
					! SO2 mult=1 when prof amount = 1000 * ref amount
					SO2MLT(L)=1.0010E-3*( SO2_PROFILE(L) - FIXMUL(L)*REF_SO2_PROFILE(L) )/ REF_SO2_PROFILE(L)
					! Ignore changes in SO2 of less than ~10%
					IF (ABS(SO2MLT(L)) .LT. 1E-4) SO2MLT(L)=0.0
				ELSE
					SO2MLT(L)=0.0
				ENDIF
			ELSE
				SO2MLT(L) = SO2MLT_SAVED(L)
			ENDIF

			IF (DO_SPECIES % DO_HNO3) THEN
				IF (LHNO3) THEN
					! HNO3 mult=1 when prof amount = 2 * ref amount
					HNOMLT(L)=( HNO3_PROFILE(L) - FIXMUL(L)*REF_HNO3_PROFILE(L) )/ REF_HNO3_PROFILE(L)
					! Ignore changes in HNO3 less than ~1%
					IF (ABS(HNOMLT(L)) .LT. 1E-2) HNOMLT(L)=0.0
				ELSE
					HNOMLT(L)=0.0
				ENDIF
			ELSE
				HNOMLT(L) = HNOMLT_SAVED(L)
			ENDIF
		
			!cc this block for testing
			!      N2OMLT(L)=0.0          
			!      SO2MLT(L)=0.0
			!      HNOMLT(L)=0.0

			! Save CO2TOP in array if intital run
			IF ((START_LAYER .EQ. 1) .AND. DO_SPECIES % DO_ALL) THEN
				CO2TOP_SAVED_ARRAY(L) = CO2TOP
			ENDIF

		ENDDO ! End loop over layers

		! Convert CO2TOP from sum to mean
		CO2TOP=CO2TOP/AMIN0(NTEBOT, LBOT)

		! Save the predictors IF calculated for all layers and all constituents. 
		IF ((START_LAYER .EQ. 1) .AND. DO_SPECIES % DO_ALL) THEN

			H2O_CONTINUUM_PRED_SAVED(:, :) = H2O_CONTINUUM_PRED(:, :)

			IF (NCHN1 .GT. 0) THEN
				FIXED_PRED1_SAVED(:, :) = FIXED_PRED1(:, :)
				H2O_PRED1_SAVED(:, :) = H2O_PRED1(:, :)
				O3_PRED1_SAVED(:, :) = O3_PRED1(:, :)
			ENDIF

			IF (NCHN2 .GT. 0) THEN
				FIXED_PRED2_SAVED(:, :) = FIXED_PRED2(:, :)
				H2O_PRED2_SAVED(:, :) = H2O_PRED2(:, :)
				O3_PRED2_SAVED(:, :) = O3_PRED2(:, :)
			ENDIF

			IF (NCHN3 .GT. 0) THEN
				FIXED_PRED3_SAVED(:, :) = FIXED_PRED3(:, :)
				H2O_PRED3_SAVED(:, :) = H2O_PRED3(:, :)
				CH4_PRED3_SAVED(:, :) = CH4_PRED3(:, :)
			ENDIF

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

			CO2TOP_MEAN_SAVED = CO2TOP

			CO2MLT_SAVED(:) = CO2MLT(:)
			SO2MLT_SAVED(:) = SO2MLT(:)
			HNOMLT_SAVED(:) = HNOMLT(:)
			N2OMLT_SAVED(:) = N2OMLT(:)
			FIXMUL_SAVED(:) = FIXMUL(:)

		ENDIF
		
		! Now decide what array to assign the pointers to

		IF (DO_SPECIES % DO_TEMPERATURE) THEN
			FIXED_PRED1_PTR => FIXED_PRED1
			FIXED_PRED2_PTR => FIXED_PRED2
			FIXED_PRED3_PTR => FIXED_PRED3
			FIXED_PRED4_PTR => FIXED_PRED4
			FIXED_PRED5_PTR => FIXED_PRED5
			FIXED_PRED6_PTR => FIXED_PRED6
			FIXED_PRED7_PTR => FIXED_PRED7
			TRACEGAS_PRED_PTR => TRACEGAS_PRED
		ELSE
			FIXED_PRED1_PTR => FIXED_PRED1_SAVED
			FIXED_PRED2_PTR => FIXED_PRED2_SAVED
			FIXED_PRED3_PTR => FIXED_PRED3_SAVED
			FIXED_PRED4_PTR => FIXED_PRED4_SAVED
			FIXED_PRED5_PTR => FIXED_PRED5_SAVED
			FIXED_PRED6_PTR => FIXED_PRED6_SAVED
			FIXED_PRED7_PTR => FIXED_PRED7_SAVED
			TRACEGAS_PRED_PTR => TRACEGAS_PRED_SAVED
		ENDIF

		IF (DO_SPECIES % DO_H2O) THEN
			H2O_CONTINUUM_PRED_PTR => H2O_CONTINUUM_PRED
			H2O_PRED1_PTR => H2O_PRED1
			H2O_PRED2_PTR => H2O_PRED2
			H2O_PRED3_PTR => H2O_PRED3
			H2O_PRED4_PTR => H2O_PRED4
			H2O_PRED5_PTR => H2O_PRED5
			H2O_PRED6_PTR => H2O_PRED6
			H2O_PRED7_PTR => H2O_PRED7
		ELSE
			H2O_CONTINUUM_PRED_PTR => H2O_CONTINUUM_PRED_SAVED
			H2O_PRED1_PTR => H2O_PRED1_SAVED
			H2O_PRED2_PTR => H2O_PRED2_SAVED
			H2O_PRED3_PTR => H2O_PRED3_SAVED
			H2O_PRED4_PTR => H2O_PRED4_SAVED
			H2O_PRED5_PTR => H2O_PRED5_SAVED
			H2O_PRED6_PTR => H2O_PRED6_SAVED
			H2O_PRED7_PTR => H2O_PRED7_SAVED
		ENDIF

		IF (DO_SPECIES % DO_O3) THEN
			O3_PRED1_PTR => O3_PRED1
			O3_PRED2_PTR => O3_PRED2
			O3_PRED4_PTR => O3_PRED4
			O3_PRED5_PTR => O3_PRED5
			O3_PRED6_PTR => O3_PRED6
			O3_PRED7_PTR => O3_PRED7
		ELSE
			O3_PRED1_PTR => O3_PRED1_SAVED
			O3_PRED2_PTR => O3_PRED2_SAVED
			O3_PRED4_PTR => O3_PRED4_SAVED
			O3_PRED5_PTR => O3_PRED5_SAVED
			O3_PRED6_PTR => O3_PRED6_SAVED
			O3_PRED7_PTR => O3_PRED7_SAVED
		ENDIF

		IF (DO_SPECIES % DO_CH4) THEN
			CH4_PRED3_PTR => CH4_PRED3
		ELSE
			CH4_PRED3_PTR => CH4_PRED3_SAVED
		ENDIF

		IF (DO_SPECIES % DO_CO) THEN
			CO_PRED4_PTR => CO_PRED4
		ELSE
			CO_PRED4_PTR => CO_PRED4_SAVED
		ENDIF

	END SUBROUTINE CALPAR__CALCULATE_PREDICTORS 

END MODULE CALPAR