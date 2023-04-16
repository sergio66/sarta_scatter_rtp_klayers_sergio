!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS
!
!    CALPAR
!
!F77====================================================================


!ROUTINE NAME:
!    CALPAR


!ABSTRACT:
!    Calculate the fast transmittance code temperature/amount/angle
!    dependent variables for a profile.


!CALL PROTOCOL:
!    CALPAR (LBOT, RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RNAMNT,
!  $               PTEMP,PFAMNT,PWAMNT,POAMNT,PCAMNT,PMAMNT,PNAMNT,
!  $          RES,SECANG,  ALAT,    FX, DZREF,
!  $         LCO2,  LN2O,  LSO2, LHNO3,LCO2PM,FIXMUL,CONPRD,
!  $       FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7,
!  $       WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7,
!  $       OPRED1,OPRED2,       OPRED4,OPRED5,OPRED6,OPRED7,
!  $       MPRED3,CPRED4,TRCPRD,CO2MLT,SO2MLT,HNOMLT,N2OMLT )


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
!    REAL arr  DZREF   ref prof layer thickness    meters
!    REAL arr  FX      fixed gases adjustment      none
!    REAL arr  PTEMP   profile temperature         K
!    REAL arr  PCAMNT  prof carbon monoxide amnt   kiloMoles/cm^2
!    REAL arr  PFAMNT  profile CO2 gas amount      kiloMoles/cm^2
!    REAL arr  PHAMNT  profile HNO3 gas amount     kiloMoles/cm^2
!    REAL arr  PMAMNT  profile methane amount      kiloMoles/cm^2
!    REAL arr  PNAMNT  profile N2O amount          kiloMoles/cm^2
!    REAL arr  POAMNT  profile ozone amount        kiloMoles/cm^2
!    REAL arr  PRES    layer pressures             atm
!    REAL arr  PSAMNT  profile SO2 amount          kiloMoles/cm^2
!    REAL arr  PWAMNT  profile water amount        kiloMoles/cm^2
!    REAL arr  RTEMP   reference temperature       K
!    REAL arr  RCAMNT  ref carbon monoxide amount  kiloMoles/cm^2
!    REAL arr  RFAMNT  reference CO2 amount        kiloMoles/cm^2
!    REAL arr  RHAMNT  reference HNO3 amount       kiloMoles/cm^2
!    REAL arr  RMAMNT  reference methane amount    kiloMoles/cm^2
!    REAL arr  RNAMNT  reference N2O amount        kiloMoles/cm^2
!    REAL arr  ROAMNT  reference ozone amount      kiloMoles/cm^2
!    REAL arr  RSAMNT  reference SO2 amount        kiloMoles/cm^2
!    REAL arr  RWAMNT  reference water amount      kiloMoles/cm^2
!    REAL arr  SECANG  secant of path angle        none


!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL arr  CO2MLT  CO2 multiplier              none
!    REAL arr  CPRED4  carbon monoxide pred set4   various
!    REAL arr  FIXMUL  fixed amount multiplier     none
!    REAL arr  FPRED1  fixed predictors set1       various
!    REAL arr  FPRED2  fixed predictors set2       various
!    REAL arr  FPRED3  fixed predictors set3       various
!    REAL arr  FPRED4  fixed predictors set4       various
!    REAL arr  FPRED5  fixed predictors set5       various
!    REAL arr  FPRED6  fixed predictors set6       various
!    REAL arr  FPRED7  fixed predictors set7       various
!    REAL arr  HNOMLT  HNO3 multiplier             none
!    REAL arr  MPRED3  methane predictors set3     various
!    REAL arr  N2OMLT  N2O multiplier              none
!    REAL arr  OPRED1  ozone predictors set1       various
!    REAL arr  OPRED2  ozone predictors set2       various
!    REAL arr  OPRED4  ozone predictors set4       various
!    REAL arr  OPRED5  ozone predictors set5       various
!    REAL arr  OPRED6  ozone predictors set6       variou
!    REAL arr  OPRED7  ozone predictors set7       various
!    REAL arr  SO2MLT  SO2 multiplier              none
!    REAL arr  TRCPRD  trace gas pert predictors   various
!    REAL arr  WPRED1  water predictors set1       various
!    REAL arr  WPRED2  water predictors set2       various
!    REAL arr  WPRED3  water predictors set3       various
!    REAL arr  WPRED4  water predictors set4       various
!    REAL arr  WPRED5  water predictors set5       various
!    REAL arr  WPRED6  water predictors set6       various
!    REAL arr  WPRED7  water predictors set7       various


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
!    CONPRD: water continuum predictors (7 terms)
!       1) a*W/Tr^2    2) a*(W/Tr^2)^2   3) a*W/Tr  4) a*W^2/Tr
!       5) a*(W/Tr)^2  6) a*W/Tr^4       7) a*Wr
!
!    -------------------------------
!    Fixed predictors
!
!    FPRED1: FWO (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
!
!    FPRED2: FOW (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
!
!    FPRED3: FMW (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
!
!    FPRED4: FCOW (11 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a^2*Trz
!       9) a^2*Tr  10) a^3     11) sqrt(a)
!
!    FPRED5: FWO (11 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
!       9) a^2*Tr  10) sqrt(a) 11) Trz
!
!    FPRED6: FWO (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) sqrt(a)
!
!    FPRED7: FWO (8 terms):
!       1) a        2) a^2      3) a*Tr    4) a*Tr^2
!       5) Tr       6) Tr^2     7) a*Trz   8) sqrt(a)
!
!    ---------------------------------
!    Water predictors
!
!    WPRED1: FWO (11 terms):
!       1) W*a           2) sqrt(W*a)       3) W*a*W/Wz
!       4) W*a*dT        5) (W*a)^2         6) sqrt(W*a)*dT
!       7) root^4(W*a)   8) sqrt(W*a)*W/Wz  9) (W*a)^3
!      10) W            11) W*a*dT*|dT|
!
!    WPRED2: FOW (11 terms):
!       1) W*a             2) sqrt(W*a)     3) W*a*dT
!       4) W*a*Ox*a        5) (W*a)^2       6) root^4(W*a)
!       7) sqrt(W*a)*dT    8) W*a*W/Wz      9) (W*a)^3
!      10) W*a*(Ox*a)^2   11) sqrt(W*a)*W/Wz
!
!    WPRED3: FMW (11 terms):
!       1) W*a             2) sqrt(W*a)     3) W*a*W/Wz
!       4) W*a*dT          5) (W*a)^2       6) sqrt(W*a)*dT
!       7) root^4(W*a)     8) (W*a)^3       9) W
!      10) sqrt(W*a)*W/Wz 11) sqrt(W*a)*Mz*a
!
!    WPRED4: FCOW (13 terms):
!       1) W*a             2) W             3) sqrt(W*a)
!       4) W*a*dT          5) (W*a)^2       6) sqrt(W*a)*dT
!       7) root^4(W*a)     8) W*a*W/Wz      9) W*a^2
!      10) (W*a)^3        11) W*a*Cz*a     12) sqrt(W*a)*W/Wz
!      13) W*a^2*dT
!
!    WPRED5: FWO bfsw (3 terms):
!       1) W*a           2) (W*a)^3/2       3) W*a*dT
!
!    WPRED6: FWO mfmw (7 terms):
!       1) W*a           2) (W*a)^3/2       3) W*a*dT
!       4) (W*a)^2       5) (W*a)^3/2*dT    6) (W*a)^3
!       7) W*a^2
!
!    WPRED7: FWO mfbw (13 terms):
!       1) W*a           2) (W*a)^3/2       3) W*a*dT
!       4) (W*a)^2       5) (W*a)^3/2*dT    6) (W*a)^3
!       7) W*a^2         8) W*a*W/Wz        9) (W*a)^3/2*W/Wz
!      10) (W*a)^5/4    11) (W*a)^2*W/Wz   12) W^2*a
!      13) (W*a)^7/4
!
!    ---------------------------
!    Ozone predictors
!
!    OPRED1: FWO (5 terms):
!       1) O*a             2) sqrt(O*a)     3) O*a*dT
!       4) (O*a)^2         5) sqrt(O*a)*dT
!
!    OPRED2: FOW (10 terms):
!       1) O*a             2) sqrt(O*a)     3) O*a*dT
!       4) (O*a)^2         5) sqrt(O*a)*dT  6) O*a*O/Ox
!       7) sqrt(O*a)*O/Ox  8) O*a*Oz/Ox     9) O*a*sqrt(Ox*a)
!      10) O*a*TOz*a
!
!    OPRED4: FCOW (3 terms):
!       1) O*a         2) sqrt(O*a)     3) O*a*dT
!
!    OPRED5: FWO bfsw (1 term):
!       1) O*a
!
!    OPRED6: FWO mfmw (1 term):
!       1) O*a
!
!    OPRED7: FWO mfbw (1 term):
!       1) O*a
!
!    ---------------------------
!    CPRED4: carbon monoxide predictors (11 terms):
!       1) C*a           2) sqrt(C*a)       3) C*a*dT
!       4) (C*a)^2       5) C*a*C/Cz        6) sqrt(C*a)*dT
!       7) root^4(C*a)   8) sqrt(C*a)*C/Cz  9) C
!
!    ---------------------------
!    MPRED3: methane predictors (9 terms):
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
!    "Tr" is the temperature ratio PTEMP/RTEMP
!    "Trz" is the pressure weighted temperature ratio above, i.e.
!      the sum i=2 to i=L of { P(i) * ( P(i) -  P(i-1) )* Tr(i-1) }
!      where "P" is the pressure PRES and "L" is the layer number, and
!      Trz(L=1)=0
!    "W" is the water amount ratio PWAMNT/RWAMNT
!    "dT" is the temperature offset PTEMP-RTEMP
!    "Wz" is the pressure weighted water amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PWAMNT(i) },
!      divided by the same sum except using RWAMNT instead of PWAMNT.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!    "O" is the ozone amount ratio POAMNT/ROAMNT
!    "Oz" is the pressure weighted ozone amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * POAMNT(i) },
!      divided by the same sum except using ROAMNT instead of POAMNT.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2)
!    "Ox" is the unweighted ozone amount above ratio, the
!      sum i=1 to i=L of { POAMNT(i) },
!      divided by the same sum except using ROAMNT instead of POAMNT.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!    "TOz" is the pressure and ozone weighted temperature ratio above,
!      sum i=2 to i=L of { P(i) * ( P(i)-P(i-1) )* dT(i-1) * O(i-1) }
!      and TOz(L=1)=0
!    "C" is the carbon monoxide amount ratio POAMNT/ROAMNT
!    "Cz" is the pressure weighted CO amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PCAMNT(i) },
!      divided by the same sum except using RCAMNT instead of PCAMNT.
!      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
!    "M" is the methane amount ratio PMAMNT/RMAMNT
!    "Mz" is the pressure weighted methane amount above ratio, the
!      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PMAMNT(i) },
!      divided by the same sum except using RMAMNT instead of PMAMNT.
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
!                               var; add FX and DZREF data; add
!                               PWATER, PMULT, and GSCAL local vars.
! 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
! 17 Aug 2000 Scott Hannon   Add FX & DZREF input vars
! 12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
! 28 Jun 2005 Scott Hannon   "trace" version for CO2,SO2,HNO3,N2O
! 23 Jan 2008 Scott Hannon   Add LCO2,LN2O,LSO2,LHNO3 switches for
!                               perturbation multiplier calcs; add
!                               LCO2PM to allow CO2 ppmv profile
! 14 May 2008 Scott Hannon   Add no prof CO2MLT calc; add CO2TOP and
!                               CO2PPM to call; add CO2TOP calc
! 30 Nov 2017 Bill Irion     Modified to be F90 compliant

!END====================================================================

!=================================================================
SUBROUTINE CALPAR ( LBOT, &
	 RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RSAMNT,RHAMNT,RNAMNT, &
	 PTEMP,PFAMNT,PWAMNT,POAMNT,PCAMNT,PMAMNT,PSAMNT,PHAMNT,PNAMNT, &
	  PRES,SECANG,  ALAT,    FX, DZREF, &
	  LCO2,  LN2O,  LSO2, LHNO3,LCO2PM,CO2PPM,CO2TOP, &
	FIXMUL,CONPRD, &
	FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7, &
	WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7, &
	OPRED1,OPRED2,       OPRED4,OPRED5,OPRED6,OPRED7, &
	MPRED3,CPRED4,TRCPRD,CO2MLT,SO2MLT,HNOMLT,N2OMLT )
!=================================================================


	USE INCFTC      
	IMPLICIT NONE

	!
	! ARGUMENTS
	!
	! Input
	INTEGER, INTENT(IN) :: LBOT
	REAL, INTENT(IN), DIMENSION(MAXLAY) ::  RTEMP
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: RFAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: RWAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: ROAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: RCAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: RMAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: RSAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: RHAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: RNAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) ::  PTEMP
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: PFAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: PWAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: POAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: PCAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: PMAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: PSAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: PHAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: PNAMNT
	REAL, INTENT(IN), DIMENSION(MAXLAY) ::   PRES
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: SECANG
	REAL, INTENT(IN) :: ALAT
	REAL, INTENT(IN), DIMENSION(MAXLAY) ::     FX
	REAL, INTENT(IN), DIMENSION(MAXLAY) ::  DZREF
	LOGICAL, INTENT(IN) :: LCO2
	LOGICAL, INTENT(IN) :: LN2O
	LOGICAL, INTENT(IN) :: LSO2
	LOGICAL, INTENT(IN) :: LHNO3
	LOGICAL, INTENT(IN) :: LCO2PM
	REAL, INTENT(IN) :: CO2PPM

	! Output
	REAL, INTENT(OUT) :: CO2TOP
	REAL, INTENT(OUT), DIMENSION(MAXLAY) ::  FIXMUL
	REAL, INTENT(OUT), DIMENSION(N1CON, MAXLAY) :: CONPRD
	REAL, INTENT(OUT), DIMENSION(N1FIX, MAXLAY) :: FPRED1
	REAL, INTENT(OUT), DIMENSION(N2FIX, MAXLAY) :: FPRED2
	REAL, INTENT(OUT), DIMENSION(N3FIX, MAXLAY) :: FPRED3
	REAL, INTENT(OUT), DIMENSION(N4FIX, MAXLAY) :: FPRED4
	REAL, INTENT(OUT), DIMENSION(N5FIX, MAXLAY) :: FPRED5
	REAL, INTENT(OUT), DIMENSION(N6FIX, MAXLAY) :: FPRED6
	REAL, INTENT(OUT), DIMENSION(N7FIX, MAXLAY) :: FPRED7
	REAL, INTENT(OUT), DIMENSION(N1H2O, MAXLAY) :: WPRED1
	REAL, INTENT(OUT), DIMENSION(N2H2O, MAXLAY) :: WPRED2
	REAL, INTENT(OUT), DIMENSION(N3H2O, MAXLAY) :: WPRED3
	REAL, INTENT(OUT), DIMENSION(N4H2O, MAXLAY) :: WPRED4
	REAL, INTENT(OUT), DIMENSION(N5H2O, MAXLAY) :: WPRED5
	REAL, INTENT(OUT), DIMENSION(N6H2O, MAXLAY) :: WPRED6
	REAL, INTENT(OUT), DIMENSION(N7H2O, MAXLAY) :: WPRED7
	REAL, INTENT(OUT), DIMENSION( N1O3, MAXLAY) :: OPRED1
	REAL, INTENT(OUT), DIMENSION( N2O3, MAXLAY) :: OPRED2
	!                            Note: There is no OPRED3
	REAL, INTENT(OUT), DIMENSION( N4O3, MAXLAY) :: OPRED4
	REAL, INTENT(OUT), DIMENSION( N5O3, MAXLAY) :: OPRED5
	REAL, INTENT(OUT), DIMENSION( N6O3, MAXLAY) :: OPRED6
	REAL, INTENT(OUT), DIMENSION( N7O3, MAXLAY) :: OPRED7
	REAL, INTENT(OUT), DIMENSION(N3CH4, MAXLAY) :: MPRED3
	REAL, INTENT(OUT), DIMENSION( N4CO, MAXLAY) :: CPRED4
	REAL, INTENT(OUT), DIMENSION(NTRACE,MAXLAY) :: TRCPRD
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
	REAL TJUNKS
	REAL WJUNKA
	REAL WJUNKR
	REAL WJUNKS
	REAL WJUNKZ
	REAL WJUNK4
	REAL OJUNKA
	REAL OJUNKR
	REAL OJUNKZ
	REAL OJUNKX
	REAL CJUNKA
	REAL CJUNKR
	REAL CJUNKS
	REAL CJUNKZ
	REAL MJUNKA
	REAL MJUNKR
	REAL MJUNKZ

	! Variables for fixed gases adjustment
	REAL PWATER
	REAL  GSCAL

	! Variables with fixed assignments
	! Changed from DATA statements by FWI 11/30/2017
	REAL, PARAMETER :: PMULT = 0.58 		! fudge factor * (0.622=M_H2O/M_AIR)	
	REAL, PARAMETER :: STDDEN = 2.6867E+19	! Loschmidt aka standard density
	REAL, PARAMETER :: STDTMP = 273.15		! Standard Temperature
	REAL, PARAMETER :: KMOLE = 6.022045E+26	! 1000 * Avagadro's Number



	!
	!
	! EXECUTABLE CODE
	!
	!

	! Calc the fixed gases gravity correction factor
	GSCAL=( 9.78050518 + 0.0518017*( COS( (ALAT - 90.0)*PI/180) )**2 )/9.80683613

	! Initialize the sum terms to zero
	PNORM = 0.E+0
	TZ = 0.0E+0
	WZREF = 0.0E+0
	WZ = 0.0E+0
	XZREF = 0.0E+0
	XZ = 0.0E+0
	OZREF = 0.0E+0
	OZ = 0.0E+0
	TOZ = 0.0E+0
	CZREF = 0.0E+0
	CZ = 0.0E+0
	MZREF = 0.0E+0
	MZ = 0.0E+0
	TMZ = 0.0E+0
	CO2TOP = 0.0E+0

	!--------------------
	! Loop over the layers
	!--------------------
	DO L = 1, LBOT

		!
		! Calculate the basic profile
		! dependent predictors.
		!
		IF (L .EQ. 1) THEN
			PDP = PRES(1) * (PRES(2) - PRES(1))
			TRZ = 0.0E+0
			TAZ_O = 0.0E+0
			TAZ_M = 0.0E+0
		ELSE
			PDP = PRES(L)*( PRES(L) - PRES(L-1) )
			PNORM = PNORM + PDP
			
			! Note: TRZ, TOZ, and TMZ use layer-above terms
			TZ = TZ + PDP*TR
			TRZ = TZ/PNORM
			
			TOZ = TOZ + PDP*DT*A_O
			TAZ_O = TOZ/PNORM
			
			TMZ = TMZ + PDP*TR*A_M
			TAZ_M = TMZ/PNORM
		ENDIF

		! Temperature terms
		DT=PTEMP(L) - RTEMP(L)
		TR=PTEMP(L)/RTEMP(L)

		! Calc the fixed gases correction term for this layer
		PWATER=KMOLE*PWAMNT(L)*PTEMP(L)/(STDDEN*STDTMP*100*DZREF(L))
		A_F=( 1 - PMULT*PWATER/PRES(L) )/( FX(L)*GSCAL )

		! for testing
		! A_F=1.0

		! Water terms
		A_W=PWAMNT(L)/RWAMNT(L)
		WZREF=WZREF + PDP*RWAMNT(L)
		WZ=WZ + PDP*PWAMNT(L)
		AZ_W=WZ/WZREF

		! Ozone terms
		A_O=POAMNT(L)/ROAMNT(L)
		XZREF=XZREF + ROAMNT(L)
		XZ=XZ + POAMNT(L)
		XZ_O=XZ/XZREF
		OZREF=OZREF + PDP*ROAMNT(L)
		OZ=OZ + PDP*POAMNT(L)
		AZ_O=OZ/OZREF

		! Carbon monoxide terms
		A_C=PCAMNT(L)/RCAMNT(L)
		CZREF=CZREF + PDP*RCAMNT(L)
		CZ=CZ + PDP*PCAMNT(L)
		AZ_C=CZ/CZREF

		! Methane terms
		A_M=PMAMNT(L)/RMAMNT(L)
		MZREF=MZREF + PDP*RMAMNT(L)
		MZ=MZ + PDP*PMAMNT(L)
		AZ_M=MZ/MZREF

		!----------------------
		! Load up the predictors
		!----------------------

		!-----
		! Fixed (for FWO, FOW, FMW, & FCOW)
		!-----
		TJUNKS=TR*TR
		FIXMUL(L)=A_F

		FPRED1(1,L)=SECANG(L)
		FPRED1(2,L)=SECANG(L)*SECANG(L)
		FPRED1(3,L)=SECANG(L)*TR
		FPRED1(4,L)=SECANG(L)*TJUNKS
		FPRED1(5,L)=TR
		FPRED1(6,L)=TJUNKS
		FPRED1(7,L)=SECANG(L)*TRZ
		FPRED1(8,L)=SECANG(L)*TRZ/TR

		FPRED2(1,L)=FPRED1(1,L)
		FPRED2(2,L)=FPRED1(2,L)
		FPRED2(3,L)=FPRED1(3,L)
		FPRED2(4,L)=FPRED1(4,L)
		FPRED2(5,L)=FPRED1(5,L)
		FPRED2(6,L)=FPRED1(6,L)
		FPRED2(7,L)=FPRED1(7,L)
		FPRED2(8,L)=FPRED1(8,L)

		FPRED3(1,L)=FPRED1(1,L)
		FPRED3(2,L)=FPRED1(2,L)
		FPRED3(3,L)=FPRED1(3,L)
		FPRED3(4,L)=FPRED1(4,L)
		FPRED3(5,L)=FPRED1(5,L)
		FPRED3(6,L)=FPRED1(6,L)
		FPRED3(7,L)=FPRED1(7,L)
		FPRED3(8,L)=FPRED1(8,L)

		FPRED4(1,L)=SECANG(L)
		FPRED4(2,L)=SECANG(L)*SECANG(L)
		FPRED4(3,L)=SECANG(L)*TR
		FPRED4(4,L)=SECANG(L)*TJUNKS
		FPRED4(5,L)=TR
		FPRED4(6,L)=TJUNKS
		FPRED4(7,L)=SECANG(L)*TRZ
		FPRED4(8,L)=SECANG(L)*SECANG(L)*TRZ
		FPRED4(9,L)=SECANG(L)*SECANG(L)*TR
		FPRED4(10,L)=SECANG(L)*SECANG(L)*SECANG(L)
		FPRED4(11,L)=SQRT(SECANG(L))

		! Fixed predictors for FWO sun bfsw = set5
		FPRED5(1,L)=SECANG(L)
		FPRED5(2,L)=SECANG(L)*SECANG(L)
		FPRED5(3,L)=SECANG(L)*TR
		FPRED5(4,L)=SECANG(L)*TJUNKS
		FPRED5(5,L)=TR
		FPRED5(6,L)=TJUNKS
		FPRED5(7,L)=SECANG(L)*TRZ
		FPRED5(8,L)=SECANG(L)*TRZ/TR
		FPRED5(9,L)=SECANG(L)*SECANG(L)*TR
		FPRED5(10,L)=SQRT(SECANG(L))
		FPRED5(11,L)=TRZ

		! Fixed predictors for FWO sun mfmw = set6
		FPRED6(1,L)=SECANG(L)
		FPRED6(2,L)=SECANG(L)*SECANG(L)
		FPRED6(3,L)=SECANG(L)*TR
		FPRED6(4,L)=SECANG(L)*TJUNKS
		FPRED6(5,L)=TR
		FPRED6(6,L)=TJUNKS
		FPRED6(7,L)=SECANG(L)*TRZ
		FPRED6(8,L)=SQRT(SECANG(L))

		! Fixed predictors for FWO sun mfbw = set7
		FPRED7(1,L)=SECANG(L)
		FPRED7(2,L)=SECANG(L)*SECANG(L)
		FPRED7(3,L)=SECANG(L)*TR
		FPRED7(4,L)=SECANG(L)*TJUNKS
		FPRED7(5,L)=TR
		FPRED7(6,L)=TJUNKS
		FPRED7(7,L)=SECANG(L)*TRZ
		FPRED7(8,L)=SQRT(SECANG(L))


		!
		! Ozone
		!
		OJUNKA=SECANG(L)*A_O
		OJUNKR=SQRT( OJUNKA )
		OJUNKZ=OJUNKA/XZ_O
		OJUNKX=SECANG(L)*XZ_O

		! Ozone predictors for FWO = set1
		OPRED1(1,L)=OJUNKA
		OPRED1(2,L)=OJUNKR
		OPRED1(3,L)=OJUNKA*DT
		OPRED1(4,L)=OJUNKA*OJUNKA
		OPRED1(5,L)=OJUNKR*DT

		! Ozone predictors for FOW = set2
		OPRED2( 1,L)=OJUNKA
		OPRED2( 2,L)=OJUNKR
		OPRED2( 3,L)=OJUNKA*DT
		OPRED2( 4,L)=OJUNKA*OJUNKA
		OPRED2( 5,L)=OJUNKR*DT
		OPRED2( 6,L)=OJUNKZ*A_O
		OPRED2( 7,L)=OJUNKR*A_O/XZ_O
		OPRED2( 8,L)=OJUNKZ*AZ_O
		OPRED2( 9,L)=OJUNKA*SQRT( OJUNKX )
		OPRED2(10,L)=OJUNKA*TAZ_O*SECANG(L)

		! There are no ozone predictors for set3 = FMW (the ozone
		! absorption in the region covered by FMW is negligible).

		! Ozone predictors for FCOW = set4
		OPRED4(1,L)=OJUNKA
		OPRED4(2,L)=OJUNKR
		OPRED4(3,L)=OJUNKA*DT

		! Ozone predictors for FWO sun bfsw = set5
		OPRED5(1,L)=OJUNKA

		! Ozone predictors for FWO sun mfmw = set6
		OPRED6(1,L)=OJUNKA

		! Ozone predictors for FWO sun mfbw = set7
		OPRED7(1,L)=OJUNKA

		!
		! Methane for FMW = set3
		!
		MJUNKA=SECANG(L)*A_M
		MJUNKR=SQRT(MJUNKA)
		MJUNKZ=SECANG(L)*AZ_M
		MPRED3(1,L)=MJUNKA
		MPRED3(2,L)=MJUNKR
		MPRED3(3,L)=MJUNKA*DT
		MPRED3(4,L)=MJUNKA*MJUNKA
		MPRED3(5,L)=MJUNKA*SECANG(L)
		MPRED3(6,L)=MJUNKZ
		MPRED3(7,L)=A_M*DT
		MPRED3(8,L)=TAZ_M*SECANG(L)
		MPRED3(9,L)=SQRT( MJUNKZ )

		!
		! Water
		! 
		WJUNKA=SECANG(L)*A_W
		WJUNKR=SQRT( WJUNKA )
		WJUNKS=WJUNKA*WJUNKA
		WJUNKZ=WJUNKA*A_W/AZ_W
		WJUNK4=SQRT( WJUNKR )

		! Water predictors for FWO = set1
		WPRED1( 1,L)=WJUNKA
		WPRED1( 2,L)=WJUNKR
		WPRED1( 3,L)=WJUNKZ
		WPRED1( 4,L)=WJUNKA*DT
		WPRED1( 5,L)=WJUNKS
		WPRED1( 6,L)=WJUNKR*DT
		WPRED1( 7,L)=WJUNK4
		WPRED1( 8,L)=WJUNKZ/WJUNKR
		WPRED1( 9,L)=WJUNKS*WJUNKA
		WPRED1(10,L)=A_W
		WPRED1(11,L)=WJUNKA*DT*ABS( DT )

		! Water predictors for FOW = set2
		WPRED2( 1,L)=WJUNKA
		WPRED2( 2,L)=WJUNKR
		WPRED2( 3,L)=WJUNKA*DT
		WPRED2( 4,L)=WJUNKA*OJUNKX
		WPRED2( 5,L)=WJUNKS
		WPRED2( 6,L)=WJUNK4
		WPRED2( 7,L)=WJUNKR*DT
		WPRED2( 8,L)=WJUNKZ
		WPRED2( 9,L)=WJUNKA*WJUNKS
		WPRED2(10,L)=WJUNKA*OJUNKX*OJUNKX
		WPRED2(11,L)=WJUNKZ/WJUNKR

		! Water predictors for FMW = set3
		WPRED3( 1,L)=WJUNKA
		WPRED3( 2,L)=WJUNKR
		WPRED3( 3,L)=WJUNKZ
		WPRED3( 4,L)=WJUNKA*DT
		WPRED3( 5,L)=WJUNKS
		WPRED3( 6,L)=WJUNKR*DT
		WPRED3( 7,L)=WJUNK4
		WPRED3( 8,L)=WJUNKS*WJUNKA
		WPRED3( 9,L)=A_W
		WPRED3(10,L)=WJUNKZ/WJUNKR
		WPRED3(11,L)=WJUNKR*MJUNKZ

		! Water predictors for FCOW = set4
		WPRED4( 1,L)=WJUNKA
		WPRED4( 2,L)=A_W
		WPRED4( 3,L)=WJUNKR
		WPRED4( 4,L)=WJUNKA*DT
		WPRED4( 5,L)=WJUNKS
		WPRED4( 6,L)=WJUNKR*DT
		WPRED4( 7,L)=WJUNK4
		WPRED4( 8,L)=WJUNKZ
		WPRED4( 9,L)=WJUNKA*SECANG(L)
		WPRED4(10,L)=WJUNKS*WJUNKA
		WPRED4(11,L)=WJUNKA*AZ_C*SECANG(L)
		WPRED4(12,L)=WJUNKZ/WJUNKR
		WPRED4(13,L)=WJUNKA*DT*SECANG(L)

		! Water predictors for FWO sun bfsw = set5
		WPRED5( 1,L)=WJUNKA
		WPRED5( 2,L)=WJUNKA*WJUNKR
		WPRED5( 3,L)=WJUNKA*DT

		! Water predictors for FWO sun mfmw = set6
		WPRED6( 1,L)=WJUNKA
		WPRED6( 2,L)=WJUNKA*WJUNKR
		WPRED6( 3,L)=WJUNKA*DT
		WPRED6( 4,L)=WJUNKS
		WPRED6( 5,L)=WJUNKA*WJUNKR*DT
		WPRED6( 6,L)=WJUNKA*WJUNKS
		WPRED6( 7,L)=WJUNKA*SECANG(L)

		! Water predictors for FWO sun mfbw = set7
		WPRED7( 1,L)=WJUNKA
		WPRED7( 2,L)=WJUNKA*WJUNKR
		WPRED7( 3,L)=WJUNKA*DT
		WPRED7( 4,L)=WJUNKS
		WPRED7( 5,L)=WJUNKA*WJUNKR*DT
		WPRED7( 6,L)=WJUNKA*WJUNKS
		WPRED7( 7,L)=WJUNKA*SECANG(L)
		WPRED7( 8,L)=WJUNKZ
		WPRED7( 9,L)=WJUNKZ*WJUNKR
		WPRED7(10,L)=WJUNKA*WJUNK4
		WPRED7(11,L)=WJUNKA*WJUNKZ
		WPRED7(12,L)=WJUNKA*A_W
		WPRED7(13,L)=WJUNKS/WJUNK4

		!
		! Water continuum (for FWO, FOW, FMW, FCOW)
		!
		CONPRD(1,L)=WJUNKA/TJUNKS
		CONPRD(2,L)=CONPRD(1,L)*A_W/TJUNKS
		CONPRD(3,L)=WJUNKA/TR
		CONPRD(4,L)=CONPRD(3,L)*A_W
		CONPRD(5,L)=CONPRD(1,L)*A_W
		CONPRD(6,L)=CONPRD(1,L)/TJUNKS
		CONPRD(7,L)=WJUNKA


		! 
		! Carbon monoxide for FCOW = set4
		! 
		CJUNKA=SECANG(L)*A_C
		CJUNKR=SQRT( CJUNKA )
		CJUNKS=CJUNKA*CJUNKA
		CJUNKZ=CJUNKA*A_C/AZ_C
		CPRED4(1,L)=CJUNKA
		CPRED4(2,L)=CJUNKR
		CPRED4(3,L)=CJUNKA*DT
		CPRED4(4,L)=CJUNKS
		CPRED4(5,L)=CJUNKZ
		CPRED4(6,L)=CJUNKR*DT
		CPRED4(7,L)=SQRT( CJUNKR )
		CPRED4(8,L)=CJUNKZ/CJUNKR
		CPRED4(9,L)=A_C
		CPRED4(10,L)=CJUNKA*SECANG(L)
		CPRED4(11,L)=CJUNKR*SECANG(L)
		
		!
		! Trace gas perturbation predictors
		! 
		! The first 4 trace predictors are used by all trace gases
		TRCPRD(1,L)=SECANG(L)
		TRCPRD(2,L)=TR
		TRCPRD(3,L)=SECANG(L)*TR
		TRCPRD(4,L)=SECANG(L)*TJUNKS
		! The last 3 trace predictors are only used by N2O
		TRCPRD(5,L)=SECANG(L)*SECANG(L)
		TRCPRD(6,L)=1.0
		TRCPRD(7,L)=SQRT( SECANG(L) )


		IF (LCO2) THEN
			IF (LCO2PM) THEN
				CO2MLT(L)=100.0*(PFAMNT(L) - CO2STD)/(3.0*CO2STD)
			ELSE
				! CO2 mult=1 when prof amount = 1.03 * ref amount
				CO2MLT(L)=33.3333*( PFAMNT(L) - FIXMUL(L)*RFAMNT(L) ) /  RFAMNT(L)
				! Ignore changes in CO2 of less than ~0.03%
				IF (ABS(CO2MLT(L)) .LT. 1E-2) CO2MLT(L)=0.0
			ENDIF
		ELSE
			CO2MLT(L)=100.0*(CO2PPM - CO2STD)/(3.0*CO2STD)
		ENDIF
		IF (L .LE. NTEBOT) THEN
			CO2TOP=CO2TOP + CO2STD*(1.0 + CO2MLT(L)*3.0E-2)
		ENDIF

		IF (LN2O) THEN
			! N2O mult=-1 when prof amount = 0.75 * ref amount
			N2OMLT(L)=4.0*( PNAMNT(L) - FIXMUL(L)*RNAMNT(L) )/ RNAMNT(L)
			! Ignore changes in N2O less than ~0.3%
			IF (ABS(N2OMLT(L)) .LT. 1E-2) N2OMLT(L)=0.0
		ELSE
			N2OMLT(L)=0.0
		ENDIF

		IF (LSO2) THEN
			! SO2 mult=1 when prof amount = 1000 * ref amount
			SO2MLT(L)=1.0010E-3*( PSAMNT(L) - FIXMUL(L)*RSAMNT(L) )/ RSAMNT(L)
			! Ignore changes in SO2 of less than ~10%
			IF (ABS(SO2MLT(L)) .LT. 1E-4) SO2MLT(L)=0.0
		ELSE
			SO2MLT(L)=0.0
		ENDIF

		IF (LHNO3) THEN
			! HNO3 mult=1 when prof amount = 2 * ref amount
			HNOMLT(L)=( PHAMNT(L) - FIXMUL(L)*RHAMNT(L) )/ RHAMNT(L)
			! Ignore changes in HNO3 less than ~1%
			IF (ABS(HNOMLT(L)) .LT. 1E-2) HNOMLT(L)=0.0
		ELSE
			HNOMLT(L)=0.0
		ENDIF
		
		!cc this block for testing
		!      N2OMLT(L)=0.0          
		!      SO2MLT(L)=0.0
		!      HNOMLT(L)=0.0

	ENDDO ! End loop over layers

	! Convert CO2TOP from sum to mean
	CO2TOP=CO2TOP/AMIN0(NTEBOT, LBOT)

	RETURN

END SUBROUTINE CALPAR 
