C=======================================================================
C ifort -c -O2 -convert big_endian -extend-source 132 -I/home/sergio/RTPV201/rtpV201_140levs/include ycalpar_jac_main.f
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALPAR version with NH3. HDO
C
!F77====================================================================


!ROUTINE NAME:
C    CALPAR


!ABSTRACT:
C    Calculate the fast transmittance code temperature/amount/angle
C    dependent variables for a profile.


!CALL PROTOCOL:
C    CALPAR (LBOT, RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RSAMNT,
C  $    RHAMNT,RNAMNT,RAAMNT, PTEMP,PFAMNT,PWAMNT,POAMNT,PCAMNT,
C  $    PMAMNT,PSAMNT,PHAMNT,PNAMNT,PAAMNT,
C  $     PRES,  SECANG, ALAT,  FX,   DZREF,
C  $     LCO2,  LN2O,  LSO2,  LNH3, LHDO, LHNO3,LCO2PM,FIXMUL,CONPRD,
C  $     FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7,
C  $     WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7,
C  $     OPRED1,OPRED2,       OPRED4,OPRED5,OPRED6,OPRED7,
C  $     MPRED3,CPRED4,TRCPRD,
C  $     CO2MLT,SO2MLT,HNOMLT,N2OMLT,NH3MLT,HDOMLT, 
C  $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
C  $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
C  $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
C  $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
C  $    MJACPRED3,CJACPRED4,TRCJACPRD,
C  $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   LBOT    bottom layer number         none
C    LOGICAL   LCO2    CO2 profile switch          none
C    LOGICAL   LN2O    N2O profile switch          none
C    LOGICAL   LSO2    SO2 profile switch          none
C    LOGICAL   LHNO3   HNO3 profile switch         none
C    LOGICAL   LNH3    NH3 profile switch          none
C    LOGICAL   LHDO    HDO profile switch          none
C    LOGICAL   LCO2PM  CO2 ppmv profile switch     none
C    REAL      ALAT    profile latitude            degrees (-90 to +90)
C    REAL arr  DZREF   ref prof layer thickness    meters
C    REAL arr  FX      fixed gases adjustment      none
C    REAL arr  PTEMP   profile temperature         K
C    REAL arr  PCAMNT  prof carbon monoxide amnt   kiloMoles/cm^2
C    REAL arr  PFAMNT  profile CO2 gas amount      kiloMoles/cm^2
C    REAL arr  PHAMNT  profile HNO3 gas amount     kiloMoles/cm^2
C    REAL arr  PMAMNT  profile methane amount      kiloMoles/cm^2
C    REAL arr  PNAMNT  profile N2O amount          kiloMoles/cm^2
C    REAL arr  POAMNT  profile ozone amount        kiloMoles/cm^2
C    REAL arr  PRES    layer pressures             atm
C    REAL arr  PSAMNT  profile SO2 amount          kiloMoles/cm^2
C    REAL arr  PAAMNT  prof ammonia (NH3) amnt     kiloMoles/cm^2
C    REAL arr  PWAMNT  profile water amount        kiloMoles/cm^2
C    REAL arr  RTEMP   reference temperature       K
C    REAL arr  RCAMNT  ref carbon monoxide amount  kiloMoles/cm^2
C    REAL arr  RFAMNT  reference CO2 amount        kiloMoles/cm^2
C    REAL arr  RHAMNT  reference HNO3 amount       kiloMoles/cm^2
C    REAL arr  RMAMNT  reference methane amount    kiloMoles/cm^2
C    REAL arr  RNAMNT  reference N2O amount        kiloMoles/cm^2
C    REAL arr  ROAMNT  reference ozone amount      kiloMoles/cm^2
C    REAL arr  RSAMNT  reference SO2 amount        kiloMoles/cm^2
C    REAL arr  RAAMNT  ref ammonia (NH3) amount    kiloMoles/cm^2
C    REAL arr  RWAMNT  reference water amount      kiloMoles/cm^2
C    REAL arr  SECANG  secant of path angle        none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  CO2MLT  CO2 multiplier              none
C    REAL arr  CPRED4  carbon monoxide pred set4   various
C    REAL arr  FIXMUL  fixed amount multiplier     none
C    REAL arr  FPRED1  fixed predictors set1       various
C    REAL arr  FPRED2  fixed predictors set2       various
C    REAL arr  FPRED3  fixed predictors set3       various
C    REAL arr  FPRED4  fixed predictors set4       various
C    REAL arr  FPRED5  fixed predictors set5       various
C    REAL arr  FPRED6  fixed predictors set6       various
C    REAL arr  FPRED7  fixed predictors set7       various
C    REAL arr  HNOMLT  HNO3 multiplier             none
C    REAL arr  MPRED3  methane predictors set3     various
C    REAL arr  N2OMLT  N2O multiplier              none
C    REAL arr  NH3MLT  NH3 multiplier              none
C    REAL arr  HDOMLT  HDO multiplier              none
C    REAL arr  OPRED1  ozone predictors set1       various
C    REAL arr  OPRED2  ozone predictors set2       various
C    REAL arr  OPRED4  ozone predictors set4       various
C    REAL arr  OPRED5  ozone predictors set5       various
C    REAL arr  OPRED6  ozone predictors set6       variou
C    REAL arr  OPRED7  ozone predictors set7       various
C    REAL arr  SO2MLT  SO2 multiplier              none
C    REAL arr  TRCPRD  trace gas pert predictors   various
C    REAL arr  WPRED1  water predictors set1       various
C    REAL arr  WPRED2  water predictors set2       various
C    REAL arr  WPRED3  water predictors set3       various
C    REAL arr  WPRED4  water predictors set4       various
C    REAL arr  WPRED5  water predictors set5       various
C    REAL arr  WPRED6  water predictors set6       various
C    REAL arr  WPRED7  water predictors set7       various


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
C    May 2008 version of the 100 layer AIRS Fast Transmittance
C    Code by L.L.Strow/S.Hannon.
C
C    Rapid transmittace algorithm predictors consisting of various gas
C    amount and temperature ratios and offsets relative to a reference
C    profile are calculated.
C
C    ===================================================================
C    The FTC profile variables computed for each layer are:
C
C    ---------------------------------
C    CONPRD: water continuum predictors (7 terms)
C       1) a*W/Tr^2    2) a*(W/Tr^2)^2   3) a*W/Tr  4) a*W^2/Tr
C       5) a*(W/Tr)^2  6) a*W/Tr^4       7) a*Wr
C
C    -------------------------------
C    Fixed predictors
C
C    FPRED1: FWO (8 terms):
C       1) a        2) a^2      3) a*Tr    4) a*Tr^2
C       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
C
C    FPRED2: FOW (8 terms):
C       1) a        2) a^2      3) a*Tr    4) a*Tr^2
C       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
C
C    FPRED3: FMW (8 terms):
C       1) a        2) a^2      3) a*Tr    4) a*Tr^2
C       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
C
C    FPRED4: FCOW (11 terms):
C       1) a        2) a^2      3) a*Tr    4) a*Tr^2
C       5) Tr       6) Tr^2     7) a*Trz   8) a^2*Trz
C       9) a^2*Tr  10) a^3     11) sqrt(a)
C
C    FPRED5: FWO (11 terms):
C       1) a        2) a^2      3) a*Tr    4) a*Tr^2
C       5) Tr       6) Tr^2     7) a*Trz   8) a*Trz/Tr
C       9) a^2*Tr  10) sqrt(a) 11) Trz
C
C    FPRED6: FWO (8 terms):
C       1) a        2) a^2      3) a*Tr    4) a*Tr^2
C       5) Tr       6) Tr^2     7) a*Trz   8) sqrt(a)
C
C    FPRED7: FWO (8 terms):
C       1) a        2) a^2      3) a*Tr    4) a*Tr^2
C       5) Tr       6) Tr^2     7) a*Trz   8) sqrt(a)
C
C    ---------------------------------
C    Water predictors
C
C    WPRED1: FWO (11 terms):
C       1) W*a           2) sqrt(W*a)       3) W*a*W/Wz
C       4) W*a*dT        5) (W*a)^2         6) sqrt(W*a)*dT
C       7) root^4(W*a)   8) sqrt(W*a)*W/Wz  9) (W*a)^3
C      10) W            11) W*a*dT*|dT|
C
C    WPRED2: FOW (11 terms):
C       1) W*a             2) sqrt(W*a)     3) W*a*dT
C       4) W*a*Ox*a        5) (W*a)^2       6) root^4(W*a)
C       7) sqrt(W*a)*dT    8) W*a*W/Wz      9) (W*a)^3
C      10) W*a*(Ox*a)^2   11) sqrt(W*a)*W/Wz
C
C    WPRED3: FMW (11 terms):
C       1) W*a             2) sqrt(W*a)     3) W*a*W/Wz
C       4) W*a*dT          5) (W*a)^2       6) sqrt(W*a)*dT
C       7) root^4(W*a)     8) (W*a)^3       9) W
C      10) sqrt(W*a)*W/Wz 11) sqrt(W*a)*Mz*a
C
C    WPRED4: FCOW (13 terms):
C       1) W*a             2) W             3) sqrt(W*a)
C       4) W*a*dT          5) (W*a)^2       6) sqrt(W*a)*dT
C       7) root^4(W*a)     8) W*a*W/Wz      9) W*a^2
C      10) (W*a)^3        11) W*a*Cz*a     12) sqrt(W*a)*W/Wz
C      13) W*a^2*dT
C
C    WPRED5: FWO bfsw (3 terms):
C       1) W*a           2) (W*a)^3/2       3) W*a*dT
C
C    WPRED6: FWO mfmw (7 terms):
C       1) W*a           2) (W*a)^3/2       3) W*a*dT
C       4) (W*a)^2       5) (W*a)^3/2*dT    6) (W*a)^3
C       7) W*a^2
C
C    WPRED7: FWO mfbw (13 terms):
C       1) W*a           2) (W*a)^3/2       3) W*a*dT
C       4) (W*a)^2       5) (W*a)^3/2*dT    6) (W*a)^3
C       7) W*a^2         8) W*a*W/Wz        9) (W*a)^3/2*W/Wz
C      10) (W*a)^5/4    11) (W*a)^2*W/Wz   12) W^2*a
C      13) (W*a)^7/4
C
C    ---------------------------
C    Ozone predictors
C
C    OPRED1: FWO (5 terms):
C       1) O*a             2) sqrt(O*a)     3) O*a*dT
C       4) (O*a)^2         5) sqrt(O*a)*dT
C
C    OPRED2: FOW (10 terms):
C       1) O*a             2) sqrt(O*a)     3) O*a*dT
C       4) (O*a)^2         5) sqrt(O*a)*dT  6) O*a*O/Ox
C       7) sqrt(O*a)*O/Ox  8) O*a*Oz/Ox     9) O*a*sqrt(Ox*a)
C      10) O*a*TOz*a
C
C    OPRED4: FCOW (3 terms):
C       1) O*a         2) sqrt(O*a)     3) O*a*dT
C
C    OPRED5: FWO bfsw (1 term):
C       1) O*a
C
C    OPRED6: FWO mfmw (1 term):
C       1) O*a
C
C    OPRED7: FWO mfbw (1 term):
C       1) O*a
C
C    ---------------------------
C    CPRED4: carbon monoxide predictors (11 terms):
C       1) C*a           2) sqrt(C*a)       3) C*a*dT
C       4) (C*a)^2       5) C*a*C/Cz        6) sqrt(C*a)*dT
C       7) root^4(C*a)   8) sqrt(C*a)*C/Cz  9) C
C
C    ---------------------------
C    MPRED3: methane predictors (9 terms):
C       1) M*a           2) sqrt(M*a)     3) M*a*dT
C       4) (M*a)^2       5) M*a^2         6) Mz*a
C       7) M*dT          8) TMz*a         9) sqrt(Mz*a)
C
C    ---------------------------
C    CO2PRD: CO2 perturbation coefs (4 terms):
C       1) a        2) Tr      3) a*Tr    4) a*Tr^2
C
C    -----
C    where:
C    "a" is the secant of the viewing angle SECANG
C    "Tr" is the temperature ratio PTEMP/RTEMP
C    "Trz" is the pressure weighted temperature ratio above, i.e.
C      the sum i=2 to i=L of { P(i) * ( P(i) -  P(i-1) )* Tr(i-1) }
C      where "P" is the pressure PRES and "L" is the layer number, and
C      Trz(L=1)=0
C    "W" is the water amount ratio PWAMNT/RWAMNT
C    "dT" is the temperature offset PTEMP-RTEMP
C    "Wz" is the pressure weighted water amount above ratio, the
C      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PWAMNT(i) },
C      divided by the same sum except using RWAMNT instead of PWAMNT.
C      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
C    "O" is the ozone amount ratio POAMNT/ROAMNT
C    "Oz" is the pressure weighted ozone amount above ratio, the
C      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * POAMNT(i) },
C      divided by the same sum except using ROAMNT instead of POAMNT.
C      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2)
C    "Ox" is the unweighted ozone amount above ratio, the
C      sum i=1 to i=L of { POAMNT(i) },
C      divided by the same sum except using ROAMNT instead of POAMNT.
C      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
C    "TOz" is the pressure and ozone weighted temperature ratio above,
C      sum i=2 to i=L of { P(i) * ( P(i)-P(i-1) )* dT(i-1) * O(i-1) }
C      and TOz(L=1)=0
C    "C" is the carbon monoxide amount ratio POAMNT/ROAMNT
C    "Cz" is the pressure weighted CO amount above ratio, the
C      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PCAMNT(i) },
C      divided by the same sum except using RCAMNT instead of PCAMNT.
C      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
C    "M" is the methane amount ratio PMAMNT/RMAMNT
C    "Mz" is the pressure weighted methane amount above ratio, the
C      sum i=1 to i=L of { P(i) * ( (P(i)-P(i-1) ) * PMAMNT(i) },
C      divided by the same sum except using RMAMNT instead of PMAMNT.
C      For these sums, term P(0) is defined as P(0)=2*P(1) - P(2).
C    "TMz" is the pressure and methane weighted temperature ratio above,
C      sum i=2 to i=L of { P(i) * ( P(i)-P(i-1) )* Tr(i-1) * M(i-1) }
C      and TMz(L=1)=0
C
C    -----------------------------------------------------
C    FIXMUL: the not-quite-fixed "fixed" amount multiplier.  The
C       value should be close to (within a few percent of) unity.
C       This term adjusts for the effects of water vapor displacement
C       and latitude dependent gravity.  The equations used below are
C       a combination of analytic adjustments for water and gravity,
C       as well as trial-and-error fudge factors to make it all work
C       accurately for any realistic surface pressure and altitude.
C    ===================================================================


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    Assumes the user has supplied vaguely realistic profile amounts
C    and temperatures.


!ROUTINE HISTORY:
C Date        Programmer     Comments
C ----------- -------------- --------------------------------------
C  1 Dec 1994 Scott Hannon   Created
C 10 Apr 1995 Scott Hannon   New header comments; redefined WZ;
C                               changed SECANG to array
C  6 Sep 1995 Scott Hannon   Correct WZ for top layer
C  3 Feb 1997 Scott Hannon   Re-wrote it for FWO+FOW+FMW+FCOW
C  7 Jul 1997 Scott Hannon   Re-wrote it for sets 1 thru 7
C 30 Sep 1997 Scott Hannon   Added CO2PRD
C  5 Mar 1998 Scott Hannon   Deleted water pred 12 & 13 of set 1 & 3
C 26 Aug 1998 Scott Hannon   Add LBOT to call; loop on LBOT instead
C                               of MAXLAY
C 31 Mar 2000 Scott Hannon   Change FIXMUL equation; add ALAT input
C                               var; add FX and DZREF data; add
C                               PWATER, PMULT, and GSCAL local vars.
C 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
C 17 Aug 2000 Scott Hannon   Add FX & DZREF input vars
C 12 Sep 2002 Scott Hannon   Add predictors 6 & 7 to H2O con
C 28 Jun 2005 Scott Hannon   "trace" version for CO2,SO2,HNO3,N2O
C 23 Jan 2008 Scott Hannon   Add LCO2,LN2O,LSO2,LHNO3 switches for
C                               perturbation multiplier calcs; add
C                               LCO2PM to allow CO2 ppmv profile
C 14 May 2008 Scott Hannon   Add no prof CO2MLT calc; add CO2TOP and
C                               CO2PPM to call; add CO2TOP calc
C 10 May 2018 C Hepplewhite  Add NH3
C  1 Feb 2019 C Hepplewhite  Add HDO

!END====================================================================

C      =================================================================
       SUBROUTINE YCALPAR_JAC ( LBOT,
     $    RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RSAMNT,RHAMNT,RNAMNT,
     $    RAAMNT,PTEMP,PFAMNT,PWAMNT,POAMNT,PCAMNT,PMAMNT,PSAMNT,PHAMNT,
     $    PNAMNT,PAAMNT,PRES,SECANG,  ALAT,    FX, DZREF,
     $    LCO2,  LN2O,  LSO2, LNH3,  LHDO, LHNO3,LCO2PM,
     $    CO2PPM,CO2TOP,FIXMUL,
     $    CONPRD,DPRED, 
     $    FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7,
     $    WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7,
     $    OPRED1,OPRED2,       OPRED4,OPRED5,OPRED6,OPRED7,
     $    MPRED3,CPRED4,TRCPRD,
     $    CO2MLT,SO2MLT,HNOMLT,N2OMLT,NH3MLT,HDOMLT,
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT)
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
       INTEGER   LBOT
       REAL  RTEMP(MAXLAY)
       REAL RFAMNT(MAXLAY)
       REAL RWAMNT(MAXLAY)
       REAL ROAMNT(MAXLAY)
       REAL RCAMNT(MAXLAY)
       REAL RMAMNT(MAXLAY)
       REAL RSAMNT(MAXLAY)
       REAL RHAMNT(MAXLAY)
       REAL RNAMNT(MAXLAY)
       REAL RAAMNT(MAXLAY)
       REAL  PTEMP(MAXLAY)
       REAL PFAMNT(MAXLAY)
       REAL PWAMNT(MAXLAY)
       REAL POAMNT(MAXLAY)
       REAL PCAMNT(MAXLAY)
       REAL PMAMNT(MAXLAY)
       REAL PSAMNT(MAXLAY)
       REAL PHAMNT(MAXLAY)
       REAL PNAMNT(MAXLAY)
       REAL PAAMNT(MAXLAY)
       REAL   PRES(MAXLAY)
       REAL SECANG(MAXLAY)
       REAL   ALAT
       REAL     FX(MAXLAY)
       REAL  DZREF(MAXLAY)
       LOGICAL LCO2
       LOGICAL LN2O
       LOGICAL LSO2
       LOGICAL LNH3
       LOGICAL LHDO
       LOGICAL LHNO3
       LOGICAL LCO2PM
       REAL CO2PPM
C
C      Output
       REAL CO2TOP
       REAL FIXMUL(MAXLAY), FIXMUL_T(MAXLAY), FIXMUL_1(MAXLAY), FIXMUL_3(MAXLAY)

       REAL CONPRD( N1CON,MAXLAY)
       REAL FPRED1( N1FIX,MAXLAY)
       REAL FPRED2( N2FIX,MAXLAY)
       REAL FPRED3( N3FIX,MAXLAY)
       REAL FPRED4( N4FIX,MAXLAY)
       REAL FPRED5( N5FIX,MAXLAY)
       REAL FPRED6( N6FIX,MAXLAY)
       REAL FPRED7( N7FIX,MAXLAY)
       REAL WPRED1( N1H2O,MAXLAY)
       REAL WPRED2( N2H2O,MAXLAY)
       REAL WPRED3( N3H2O,MAXLAY)
       REAL WPRED4( N4H2O,MAXLAY)
       REAL WPRED5( N5H2O,MAXLAY)
       REAL WPRED6( N6H2O,MAXLAY)
       REAL WPRED7( N7H2O,MAXLAY)
       REAL OPRED1(  N1O3,MAXLAY)
       REAL OPRED2(  N2O3,MAXLAY)
       REAL OPRED4(  N4O3,MAXLAY)
       REAL OPRED5(  N5O3,MAXLAY)
       REAL OPRED6(  N6O3,MAXLAY)
       REAL OPRED7(  N7O3,MAXLAY)
       REAL  DPRED(  NHDO,MAXLAY)
       REAL MPRED3( N3CH4,MAXLAY)
       REAL CPRED4(  N4CO,MAXLAY)
       REAL TRCPRD(NTRACE,MAXLAY)

       LOGICAL DOJAC
       INTEGER NWANTJ          ! number of wanted jacs (default 0=none)
       INTEGER  LISTJ(MAXPRO)  ! list of wanted channels
       !!! first index is the d/dT   next six are the d/dQi (i=GASID 1,2,3,4,5,6 9 12) ordered that way
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

       REAL CO2MLT(MAXLAY)
       REAL SO2MLT(MAXLAY)
       REAL HNOMLT(MAXLAY)
       REAL N2OMLT(MAXLAY)
       REAL NH3MLT(MAXLAY)
       REAL HDOMLT(MAXLAY)

       REAL CO2JACMLT(MAXLAY)
       REAL SO2JACMLT(MAXLAY)
       REAL HNOJACMLT(MAXLAY)
       REAL N2OJACMLT(MAXLAY)
       REAL NH3JACMLT(MAXLAY)
       REAL HDOJACMLT(MAXLAY)

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      L, IWHICHJAC
       REAL    PDP
       REAL  PNORM
       REAL     DT, DTZ_O, DT_T, DT_1, DT_3, DT_5, DT_6
       REAL     TR, TR_T, TR_1, TR_3,       TR_6
       REAL     TZ, TZ_T, TZ_1, TZ_3, TZ_N
       REAL    TRZ, TRZ_T, TRZ_1, TRZ_3, TRZ_N
       REAL    A_F, A_F_T, A_F_1, A_F_3
       REAL    A_W, A_W_T, A_W_1, A_W_3
       REAL  WZREF
       REAL     WZ, WZ_T,   WZ_1,   WZ_3
       REAL   AZ_W, AZ_W_T, AZ_W_1, AZ_W_3
       REAL    A_O, A_O_T, A_O_1, A_O_3
       REAL  XZREF
       REAL     XZ, XZ_T, XZ_1, XZ_3
       REAL   XZ_O, XZ_O_T, XZ_O_1, XZ_O_3
       REAL  OZREF
       REAL     OZ
       REAL   AZ_O, AZ_O_T, AZ_O_1, AZ_O_3
       REAL    TOZ, TOZ_T, TOZ_1, TOZ_3
       REAL  TAZ_O, TAZ_O_T, TAZ_O_1, TAZ_O_3
       REAL    A_C, A_C_T, A_C_1, A_C_3, A_C_5
       REAL     CZ
       REAL  CZREF
       REAL   AZ_C, AZ_C_T, AZ_C_1, AZ_C_3, AZ_C_5
       REAL    A_M, A_M_T, A_M_1, A_M_3, A_M_6
       REAL  MZREF
       REAL     MZ, MZ_T, MZ_1, MZ_3, MZ_6
       REAL   AZ_M, AZ_M_T, AZ_M_1, AZ_M_3, AZ_M_6
       REAL    TMZ, TMZ_T,   TMZ_1,   TMZ_3,   TMZ_6
       REAL  TAZ_M, TAZ_M_T, TAZ_M_1, TAZ_M_3, TAZ_M_6
       REAL TJUNKS, TJUNKS_T, TJUNKS_1, TJUNKS_3, TJUNKS_N
       REAL WJUNKA, WJUNKA_T, WJUNKA_1, WJUNKA_3
       REAL WJUNKR, WJUNKR_T, WJUNKR_1, WJUNKR_3
       REAL WJUNKS, WJUNKS_T, WJUNKS_1, WJUNKS_3
       REAL WJUNKZ, WJUNKZ_T, WJUNKZ_1, WJUNKZ_3
       REAL WJUNK4, WJUNK4_T, WJUNK4_1, WJUNK4_3
       REAL DJUNKA, DJUNKA_T, DJUNKA_1, DJUNKA_3
       REAL DJUNKR, DJUNKR_T, DJUNKR_1, DJUNKR_3
       REAL DJUNKS, DJUNKS_T, DJUNKS_1, DJUNKS_3
       REAL DJUNKZ, DJUNKz_T, DJUNKZ_1, DJUNKZ_3
       REAL DJUNK4, DJUNK4_T, DJUNK4_1, DJUNK4_3
       REAL OJUNKA, OJUNKA_T, OJUNKA_1, OJUNKA_3
       REAL OJUNKR, OJUNKR_T, OJUNKR_1, OJUNKR_3
       REAL OJUNKZ, OJUNKZ_T, OJUNKZ_1, OJUNKZ_3
       REAL OJUNKX, OJUNKX_T, OJUNKX_1, OJUNKX_3
       REAL CJUNKA, CJUNKA_T, CJUNKA_1, CJUNKA_3, CJUNKA_5  !! CO
       REAL CJUNKR, CJUNKR_T, CJUNKR_1, CJUNKR_3, CJUNKR_5  !! CO
       REAL CJUNKS, CJUNKS_T, CJUNKS_1, CJUNKS_3, CJUNKS_5  !! CO
       REAL CJUNKZ, CJUNKZ_T, CJUNKZ_1, CJUNKZ_3, CJUNKZ_5  !! CO
       REAL MJUNKA, MJUNKA_T, MJUNKA_1, MJUNKA_3, MJUNKA_6  !! CH4
       REAL MJUNKR, MJUNKR_T, MJUNKR_1, MJUNKR_3, MJUNKR_6  !! CH4
       REAL MJUNKZ, MJUNKZ_T, MJUNKZ_1, MJUNKZ_3, MJUNKZ_6  !! CH4
       REAL SUM_PDP_OVER_TREF
       INTEGER INTERSECT

C      Variables for fixed gases adjustment
       REAL PWATER, PWATER_T, PWATER_1, PWATER_3
       REAl  GSCAL
C      variables with DATA assignments
       REAL  PMULT
       REAL STDDEN
       REAL STDTMP
       REAL KMOLE
C      Data statments
       DATA PMULT /0.58/          ! fudge factor * (0.622=M_H2O/M_AIR)
       DATA STDDEN /2.6867E+19/   ! Loschmidt aka standard density
       DATA STDTMP /273.15/       ! Standard Temperature
       DATA KMOLE /6.022045E+26/  ! 1000 * Avagadro's Number


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
C      Calc the fixed gases gravity correction factor
       GSCAL =( 9.78050518 + 0.0518017*( COS( (ALAT - 90.0)*PI/180))**2 )/9.80683613
C
C      Initialize the sum terms to zero
       PNORM  = 0.0E+0
       TZ     = 0.0E+0
       WZREF  = 0.0E+0
       WZ     = 0.0E+0
       XZREF  = 0.0E+0
       XZ     = 0.0E+0
       OZREF  = 0.0E+0
       OZ     = 0.0E+0
       TOZ    = 0.0E+0
       CZREF  = 0.0E+0
       CZ     = 0.0E+0
       MZREF  = 0.0E+0
       MZ     = 0.0E+0
       TMZ    = 0.0E+0
       CO2TOP = 0.0E+0
       DTZ_O  = 0.0E+0
       SUM_PDP_OVER_TREF= 0.0E+0
C
       if (DEBUG) write(6,'(A,L3,ES11.3)') 'calpar: LCO2PM,CO2PPM ',LCO2PM,CO2PPM
C      --------------------
C      Loop over the layers
C      --------------------
       DO L = 1,LBOT
          include "ycalpar_INIT.f"
          include "ycalpar_predsINC.f"
          include "ycalpar_SWITCHES.f"    

          IF (DOJAC) THEN
            include "ycalpar_jacINIT.f"
            include "ycalpar_TjacpredsINC.f"
            include "ycalpar_WVjacpredsINC.f"
            include "ycalpar_OZjacpredsINC.f"
            include "ycalpar_GXjacpredsINC.f"
          END IF

       ENDDO
C      End loop over layers
C
C      Convert CO2TOP from sum to mean
       CO2TOP = CO2TOP/AMIN0(NTEBOT, LBOT)
C
C
       RETURN
       END
