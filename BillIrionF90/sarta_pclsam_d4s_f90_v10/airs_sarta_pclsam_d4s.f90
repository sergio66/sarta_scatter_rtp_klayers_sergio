!=======================================================================
!
!    University of Maryland Baltimore County (UMBC)
!
!    AIRS
!
!    SARTA_pclsam
!
!F77====================================================================


!ROUTINE NAME:
!    SARTA_pclsam


!ABSTRACT:
!    Program to quickly compute simulated AIRS radiances.
!
!    This variant of the code allows for the modelling of a
!    up to two scatter clouds using the PCLSAM method.


!CALL PROTOCOL
!    none (main program)


!INPUT PARAMETERS:
!    none


!OUTPUT PARAMETERS:
!    none


!INPUT/OUTPUT PARAMETERS:
!    none


!RETURN VALUES:
!    none


!PARENT(S)
!    none


!ROUTINES CALLED:
!    CALOWP : calc OPTRAN water predictors 
!    CALPAR : calculate a profile's predictors
!    CALRAD : calc radiance
!    CALT1  : calc effective layer trans for set1 (FWO)
!    CALT2  : calc effective layer trans for set2 (FOW)
!    CALT3  : calc effective layer trans for set3 (FMW)
!    CALT4  : calc effective layer trans for set4 (FCOW)
!    CALT5  : calc effective layer trans for set5 (FWO bfsw)
!    CALT6  : calc effective layer trans for set6 (FWO mfmw)
!    CALT7  : calc effective layer trans for set7 (FWO mfbw)
!    FAKETZ : calc a "fake" (rough approx) surface-to-space trans
!    RDCOEF : read the fast transmittance coefficients
!    RDPROF : read profile data
!    RDSUN  : read the solar radiance datafile
!    SUNPAR : calc a profile's predictors for sets4-7 (sun channels)
!    CALNTE : calc radiance contribution for non-LTE


!FILES ACCESSED:
!    incFTC.f : include file of parameter statements accessed during
!       compilation only.
!    unit IOUN: used by routines RDCOEF and RDPROF.
!    unit 6: USEFAST text messages to the screen
!    unit 5: USEFAST user input instructions, etc
!    unit 10: USEFAST output radiance, text file(s)


!COMMON BLOCKS
!    COMLEV : layer boundary pressure levels


!DESCRIPTION:
!    Dec 2005 version of the SARTA (Stand-Alone Rapid
!    Transmittance Algorith with RTP I/O) by
!    L.L.Strow, S.Hannon, and H.Mottler
!
!    Computes radiances for the layers profiles contained in the
!    input RTP file.  This is the main program, and consists
!    primarily of calls to the external routines to do most of the
!    computing.


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    This program is only intended as a demo of the fast model.


!ROUTINE HISTORY:
! Date         Programmer    Comments
! ----------- -------------- -------------------------------------------
! 01 Dec 1994 Scott Hannon   Created
! 10 Apr 1995 Scott Hannon   New header comments; added ALT; new
!                            external function VACONV; SECANG may
!                            vary with layer
! 03 Jul 1995 Scott Hannon   Add parameter DZ/RDZ to RDPROF call
! 03 Feb 1997 Scott Hannon   Re-written for FWO+FOW+FMW+FCOW
! 12 Sep 1997 Scott Hannon   Re-written for 7 sets and reflected sun
!                            and downwelling thermal
! 30 Sep 1997 Scott Hannon   Added variable CO2
! 27 Feb 1998 Scott Hannon   Added OPTRAN water
! 26 Aug 1998 Scott Hannon   Added LBOT to calls to CALPAR, CALOWP,
!                            and SUNPAR; rename TBOT to TSURF; calc
!                            fractional bottom layer temperature and
!                            put it in TEMPERATURE_PROFILE(LBOT)
! 15 Oct 1999 Scott Hannon   Add ANGMAX and re-arrange angle conv
! 31 Mar 2000 Scott Hannon   Redid calpar for FIXMUL and added getbot
! 15 Mar 2001 Scott Hannon   Major re-write for RTP
! 03 May 2001 Scott Hannon   Add COMLEV; add PLEV to getbot call
! 13 Sep 2001 Scott Hannon   Changes to check of FCHAN vs FREQ
! 01 Nov 2002 Scott Hannon   Added SATZEN & SALT to RDRTP call, and
!                            if valid use SATZEN rather than SATANG
! 03 Jan 2003 Scott Hannon   Delete SUNSEC, add XZ & SUNFDG & code
!                            to fudge large sun angles (previously
!                            sunang>80 were treated as no sun).
! 24 Jul 2003 Scott Hannon   Fix error in TEMPERATURE_PROFILE(LBOT) calc for
!                            bottom fractional layer; add PLAY
! 06 Feb 2004 Scott Hannon   Add call to TUNMLT; add call to MEAN_T
!                            and associated prep code; add PTYPE
!                            to OPNRTP call; add LRHOT to RDINFO,
!                            OPNRTP, & SETEMS calls.
! 20 Dec 2004 Scott Hanonn   Add NLAY to getbot.f call; add PTYPE
!                            to rdrtp_so2.f call
! 18 May 2005 Scott Hannon   Add HNO3 based on SO2 code
! 28 Jun 2005 Scott Hannon   "trace" version for CO2,SO2,HNO3,N2O
! 13 Oct 2005 Scott Hannon   Add non-LTE
! 08 Dec 2005 Scott Hannon   Update tunmlt call for non-LTE tuning
! 29 Mar 2006 Scott Hannon   Add clouds; change TAU from trans to od
! 26 Apr 2006 Scott Hannon   Add black clouds. Redo cloud fractions.
!                            Use RTP v1.06 cloud2 fields instead
!                            of udef(11-17). {Unfinished}
! 22 Dec 2006 Scott Hannon   New & revised code associated with the
!                            new & revised calrad* routines for more
!                            fexible cloud types including black
!                            clouds; changes to calt* calls and
!                            change TAUZ from (1 x n) to (m x n).
! 22 Jan 2007 Scott Hannon   Minor fix of cfrac checks
! 02 May 2007 Scott Hannon   Replace hardcoded default SALT value
!                            with XSALT from incFTC.
! 15 Nov 2007 Scott Hannon   Move most cloud prep to GETCLD & BKPREP.
! 31 Jan 2008 Scott Hannon   Add LCO2PM to allow CO2 profile in ppmv;
!                               add LCO2,LN2O,LSO2,LHNO3 switches
! 24 Mar 2008 Scott Hannon   Add COSDAZ
! 26 Nov 2008 Scott Hannon   Update for rtpV2101
! 01 Dec 2008 Scott Hannon   Add CSTMP1/2
! 05 Sep 2017 Bill Irion     Converted from F77 standalone program
!                              to F90 callable subroutine. Moved initial
!                              reads to airs_sarta_variables
! 01 Dec 2017 Bill Irion	 Added local TEMPERATURE_PROFILE copied from TEMPERATURE_PROFILE_IN so
!                            that TEMPERATURE_PROFILE can be modified without the changes
!                            being passed back to the calling program.
! 20 Dec 2017 Bill Irion 	 Modified some variable names to something approaching English.
!                            Added in IF statements to not call routines if the channel is
! 							 not in the set.

!END====================================================================

!=================================================================
SUBROUTINE AIRS_SARTA_PCLSAM_D4S( &
	PTYPE, & ! profile type (see airs_sarta_variables.f90)
	NCHAN, & ! number of channels to calculate
	SELCHAN, & ! Selected channels (unit offset)
	FREQCHAN, &
	LAT, & ! latitude of observation
	LON, & ! longitude of observation
	SATANG, & 
	SATZEN, &
	SUNANG, &
	SUN_AZIMUTH, &
	SAT_AZIMUTH, &
	PSURF, & ! surface pressure
	TSURF, &
	CLOUD_PCLSAM, &
	CLOUD_D4S_LAYER, &
	CLOUD_TRANSMISSIVITY_D4S, &
	TEMPERATURE_PROFILE_IN, &
	H2O_PROFILE, &
	CO2_PROFILE, &
	O3_PROFILE, &
	CO_PROFILE, &
	CH4_PROFILE, &
	SO2_PROFILE, &
	HNO3_PROFILE, &
	N2O_PROFILE, &
	ALT, &
	EMIS, &
	RHOTHR_IN, &
	RHOSUN_IN, &
	RAD, &
	BEGIN_LAYER, &		! optional
	SELECTED_SPECIES, & ! optional
	LEN_SELECTED_SPECIES)  	! optional

	!
	! INCLUDE FILES
	!
	USE AIRS_SARTA_VARIABLES  ! This implicitly imports INCFTC

	USE SELECT_SPECIES

	USE CALPAR
	USE SUNPAR

	USE CALT1
	USE CALT2
	USE CALT3
	USE CALT4
	USE CALT5
	USE CALT6
	USE CALT7

	USE CALT4_SOLAR
	USE CALT5_SOLAR
	USE CALT6_SOLAR
	USE CALT7_SOLAR

	USE GETBOT

	USE SECANG_MODULE

	IMPLICIT NONE

	!
	! EXTERNAL FUNCTIONS
	!
	REAL VACONV
	REAL SACONV

	!
	! ARGUMENTS
	!

	TYPE(CLOUD_TYPE_PCLSAM), INTENT(IN) :: CLOUD_PCLSAM 

	INTEGER, INTENT(IN) :: CLOUD_D4S_LAYER  ! layer number for D4S cloud
											! IF > 0, then use D4S and not PCLSAM
	REAL, INTENT(IN), DIMENSION(MXCHAN) :: CLOUD_TRANSMISSIVITY_D4S  ! transmissivity for D4S cloud

	INTEGER, INTENT(IN) :: PTYPE
	INTEGER, INTENT(IN) :: NCHAN
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: SELCHAN  ! selected channels
	REAL, INTENT(IN), DIMENSION(MXCHAN) :: FREQCHAN  ! selected channel frequencies
	REAL, INTENT(IN) :: LAT
	REAL, INTENT(IN) :: LON
	REAL, INTENT(IN) :: SATANG      ! input satellite scan angle (degrees)
	REAL, INTENT(IN) :: SATZEN      ! input satellite zenith angle (degrees)
	REAL, INTENT(IN) :: SUNANG		! solar zenith angle (at 0 altitude)
	REAL, INTENT(IN) :: SUN_AZIMUTH	! solar azimuth angle (degrees)
	REAL, INTENT(IN) :: SAT_AZIMUTH ! satellite azimuth angle (degrees)
	REAL, INTENT(IN) :: PSURF		! surface pressure (mb)
	REAL, INTENT(IN) :: TSURF		! surface temperature (K)
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: TEMPERATURE_PROFILE_IN ! layer temperatures (to remain unchanged from calling routine)
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: H2O_PROFILE ! water slab columns
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: O3_PROFILE ! ozone slab columns
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: CO_PROFILE ! carbon monoxide slab columns
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: CH4_PROFILE ! methane slab columns
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: CO2_PROFILE ! carbon dioxide slab columns
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: SO2_PROFILE ! sulphur dioxide slab columns
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: HNO3_PROFILE ! nitric acid slab columns
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: N2O_PROFILE ! nitrous oxide slab columns
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: ALT ! layer altitudes
	REAL, INTENT(IN), DIMENSION(MXCHAN) :: EMIS ! emissivity by channel
	REAL, INTENT(IN), DIMENSION(MXCHAN) :: RHOTHR_IN ! thermal reflectivity by channel (to remain unchanged from calling routine)
	REAL, INTENT(IN), DIMENSION(MXCHAN) :: RHOSUN_IN ! solar reflectivity by channel (to remain unchanged from calling routine)

	INTEGER, INTENT(IN), OPTIONAL :: BEGIN_LAYER
	INTEGER, INTENT(IN), OPTIONAL :: LEN_SELECTED_SPECIES
	CHARACTER (LEN = *), INTENT(IN), OPTIONAL :: SELECTED_SPECIES

	!-----------------------------------------------------------------------
	! OUTPUT VARIABLES
	!-----------------------------------------------------------------------
	
	REAL, INTENT(OUT), DIMENSION(MXCHAN) :: RAD ! output radiance	

	!-----------------------------------------------------------------------
	! LOCAL VARIABLES
	!-----------------------------------------------------------------------

	! local copy of TEMPERATURE_PROFILE_IN, so that the temperature can be modified without the
	! changes being passed back to the calling routine. ! FWI 12/1/17

	REAL, PARAMETER :: PI_OVER_180 = 0.01745329251994

	REAL, DIMENSION(MAXLAY) :: TEMPERATURE_PROFILE ! layer temperatures (can be modified)

	! local copies of RHOTHR and RHOSUN so that these can be modified without the
	! changes being passed back to the calling routine. ! FWI 12/15/17

	REAL, DIMENSION(MXCHAN) :: RHOTHR ! thermal reflectivity by channel (to remain unchanged from calling routine)
	REAL, DIMENSION(MXCHAN) :: RHOSUN ! solar reflectivity by channel (to remain unchanged from calling routine)

	! for surface
	INTEGER :: LBOT				! bottom layer index number
	INTEGER :: NEMIS			! # of emis pts
	REAL :: BLMULT				! bottom layer fractional multiplier

	! Other variables for the sun
	REAL COSDAZ         ! cosine(solazi - satazi) {COS Delta AZimuth}
	REAL SZALAY         ! solar zenith angle in some layer
	REAL SUNCOS         ! cosine of sun zenith angle
	REAL SCOS1          ! cosine of sun zenith angle at layer1
	REAL SUNFDG         ! fudge factor for large solar angles
	!REAL SECSUN(MAXLAY) ! secant of effective sun local path angle  ! moved to SECANG_MODULE ! FWI 1/17/18
	!REAL DISTES         ! distance of Earth from the sun
	REAL TAUZSN(MAXLAY,MXCHAN) ! sun space-to-surface-to-space OD
	!LOGICAL DOSUN       ! do sun calc? ! moved to SECANG_MODULE ! FWI 1/17/18

	! for satellite viewing angle
	REAL    SALT        ! input satellite altitude (kilometers)
	REAL    SVA         ! satellite viewing angle (degrees)

	! Predictors from CALPAR
	REAL, DIMENSION(:, :), POINTER :: H2O_CONTINUUM_PRED_PTR ! water continuum predictors
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

!
!	Basic cloud info  - moved to structured variable FWI 9/12/17
!	REAL XCEMI1(MXCHAN)    ! cloud1 emissivity
!	REAL XCEMI2(MXCHAN)    ! cloud2 emissivity
!	REAL XCRHO1(MXCHAN)    ! cloud1 reflectivity
!	REAL XCRHO2(MXCHAN)    ! cloud2 reflectivity
!	REAL CFRAC1            ! cloud1(total) fraction of FOV
!	REAL CFRAC2            ! cloud2(total) fraction of FOV
!	REAL CFRA1X            ! cloud1(exclusively) fraction of FOV
!	REAL CFRA2X            ! cloud2(exclusively) fraction of FOV
!	REAL CFRA12            ! cloud1+2(both) fraction of FOV
!	REAL CNGWA1            ! cloud1 non-gases water
!	REAL CNGWA2            ! cloud1 non-gases water
!	REAL CPRBO1            ! cloud1 bottom pressure
!	REAL CPRBO2            ! cloud2 bottom pressure
!	REAL CPRTO1            ! cloud1 top pressure
!	REAL CPRTO2            ! cloud2 top pressure
!	REAL CPSIZ1            ! cloud1 particle size
!	REAL CPSIZ2            ! cloud2 particle size
!	REAL CSTMP1            ! cloud1 top/surf temperature
!	REAL CSTMP2            ! cloud2 top/surf temperature
!	REAL FCLEAR            ! clear (no cloud) fraction of FOV
	REAL TEMPC1            ! cloud1 frac layer (above cloud) mean temp
	REAL TEMPC2            ! cloud2 frac layer (above cloud) mean temp
!	INTEGER CTYPE1         ! cloud1 type code number
!	INTEGER CTYPE2         ! cloud2 type code number

!	Exclusive cloud fractions
	
	REAL :: CFRA1X, CFRA2X
	REAL :: FCLEAR  ! Fraction of footprint that is clear

!
!	for GETMIE
!	LOGICAL LBLAC1  ! black cloud1? {Mie cloud if false}
!	LOGICAL LBLAC2  ! black cloud2? {Mie cloud if false}
	INTEGER INDMI1  ! index in MIETYP for CTYPE1
	INTEGER INDMI2  ! index in MIETYP for CTYPE2
	INTEGER  IERR1  ! error level of CTYPE1/MIETYP match
	INTEGER  IERR2  ! error level of CTYPE2/MIETYP match
!
!	for CCPREP cloud1
	INTEGER LCBOT1         ! layer containing cloud bottom
	INTEGER LCTOP1         ! layer containing cloud top
	REAL  CLRB1            ! frac of layer at bottom of cloud clear
	REAL  CLRT1            ! frac of layer at top of cloud clear
	REAL TCBOT1            ! temperature at cloud bottom
	REAL TCTOP1            ! temperature at cloud top
	REAL MASEC1            ! mean cloud view angle secant
	REAL MASUN1            ! mean cloud sun-only angle secant
	REAL CFRCL1(MAXLAY)    ! fraction of cloud in layer
	REAL G_ASY1(MXCHAN)    ! "g" asymmetry
	REAL NEXTO1(MXCHAN)    ! nadir extinction optical depth
	REAL NSCAO1(MXCHAN)    ! nadir scattering optical depth
!
!	for CCPREP cloud2
	INTEGER LCBOT2         ! layer containing cloud bottom
	INTEGER LCTOP2         ! layer containing cloud top
	REAL  CLRB2            ! frac of layer at bottom of cloud clear
	REAL  CLRT2            ! frac of layer at top of cloud clear
	REAL TCBOT2            ! temperature at cloud bottom
	REAL TCTOP2            ! temperature at cloud top
	REAL MASEC2            ! mean cloud view angle secant
	REAL MASUN2            ! mean cloud sun-only angle secant
	REAL CFRCL2(MAXLAY)    ! fraction of cloud in layer
	REAL G_ASY2(MXCHAN)    ! "g" asymmetry
	REAL NEXTO2(MXCHAN)    ! nadir extinction optical depth
	REAL NSCAO2(MXCHAN)    ! nadir scattering optical depth
!
!	used locally only
	INTEGER I      ! loop counter
	INTEGER L      ! loop counter
	REAL JUNK
!	INTEGER rtpclose    ! for call to RTP close interface routine
	REAL    EVA         ! (Earth) local viewing angle
	REAL ANGMAX         ! maximum allowed viewing angle
	REAL RJUNK1         ! junk/work
	REAL RJUNK2         ! another junk/work
!	REAL CO2PPM         ! Profile mean dry air CO2 mixing ratio
!	REAL PLAY(MAXLAY)   ! layer mean pressure
!	REAL C1V3           ! rad constant c1 times freq^3
!	REAL C2V            ! rad constant c2 times freq
	REAL VSTORE(6)      ! temporary storage for various variables
!
!	for function QIKEXP
	REAL QIKEXP


!
!      for RDCLDT
!       INTEGER MIENPS(NMIETY)            ! number of particle sizes
!       REAL  MIEPS(MXMIEA,NMIETY)        ! Mie particle size for table
!       REAL MIEABS(MXCHAN,MXMIEA,NMIETY) ! Mie absorption table
!       REAL MIEEXT(MXCHAN,MXMIEA,NMIETY) ! Mie extinction table
!       REAL MIEASY(MXCHAN,MXMIEA,NMIETY) ! Mie asymmetry table
       

! 	For now, do NOT turn off calculations for CO2, N2O, SO2, HNO3
! 	FWI 8/31/17
	LOGICAL, PARAMETER :: LCO2 = .TRUE.
	LOGICAL, PARAMETER :: LN2O = .TRUE.
	LOGICAL, PARAMETER :: LSO2 = .TRUE.
	LOGICAL, PARAMETER :: LHNO3 = .TRUE.
! 	Do not assume that CO2 input is in PPMV	
	LOGICAL, PARAMETER :: LCO2PM = .FALSE.

	! Copies of selected species (calculate new predictors for this species
	! and used saved predictors for others)

	CHARACTER (LEN = 20) :: SELECTED_SPECIES_COPY
	INTEGER :: START_LAYER 

	TYPE (DO_SPECIES_LOGICAL_ARRAY_TYPE) DO_SPECIES


!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none
!
!***********************************************************************
!***********************************************************************
!                    EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!

	! Check inputs
	IF (.FALSE.) THEN 
!		PRINT *, "PTYPE ", PTYPE
!		PRINT *, "NCHAN ", NCHAN
!		PRINT *, "SELCHAN ", SELCHAN(1:NCHAN)
!		PRINT *, "FREQCHAN ", FREQCHAN
!		PRINT *, "LAT ", LAT
!		PRINT *, "LON ", LON
!		PRINT *, "SATANG ", SATANG
!		PRINT *, "SATZEN ", SATZEN
!		PRINT *, "SUN_AZIMUTH ", SUN_AZIMUTH
!		PRINT *, "SAT_AZIMUTH ", SAT_AZIMUTH
!		PRINT *, "PSURF ", PSURF
!		PRINT *, "TSURF ", TSURF
		PRINT *, "CLOUD_PCLSAM%CTYPE1 ", CLOUD_PCLSAM%CTYPE1
		PRINT *, "CLOUD_PCLSAM%CFRAC1 ", CLOUD_PCLSAM%CFRAC1
		!PRINT *, "CLOUD_PCLSAM%CEMIS1 ", CLOUD_PCLSAM%CEMIS1
		!PRINT *, "CLOUD_PCLSAM%CRHO1 ", CLOUD_PCLSAM%CRHO1
		PRINT *, "CLOUD_PCLSAM%CPRTOP1 ", CLOUD_PCLSAM%CPRTOP1
		PRINT *, "CLOUD_PCLSAM%CPRBOT1 ", CLOUD_PCLSAM%CPRBOT1
		PRINT *, "CLOUD_PCLSAM%CNGWAT1 ", CLOUD_PCLSAM%CNGWAT1
		PRINT *, "CLOUD_PCLSAM%CPSIZE1 ", CLOUD_PCLSAM%CPSIZE1
		PRINT *, "CLOUD_PCLSAM%CSTEMP1 ", CLOUD_PCLSAM%CSTEMP1
		PRINT *, "CLOUD_PCLSAM%LBLAC1 ", CLOUD_PCLSAM%LBLAC1
		PRINT *, "CLOUD_PCLSAM%CTYPE2 ", CLOUD_PCLSAM%CTYPE2
		PRINT *, "CLOUD_PCLSAM%CFRAC2 ", CLOUD_PCLSAM%CFRAC2
		!PRINT *, "CLOUD_PCLSAM%CEMIS2 ", CLOUD_PCLSAM%CEMIS2
		!PRINT *, "CLOUD_PCLSAM%CRHO2 ", CLOUD_PCLSAM%CRHO2
		PRINT *, "CLOUD_PCLSAM%CPRTOP2 ", CLOUD_PCLSAM%CPRTOP2
		PRINT *, "CLOUD_PCLSAM%CPRBOT2 ", CLOUD_PCLSAM%CPRBOT2
		PRINT *, "CLOUD_PCLSAM%CNGWAT2 ", CLOUD_PCLSAM%CNGWAT2
		PRINT *, "CLOUD_PCLSAM%CPSIZE2 ", CLOUD_PCLSAM%CPSIZE2
		PRINT *, "CLOUD_PCLSAM%CSTEMP2 ", CLOUD_PCLSAM%CSTEMP2
		PRINT *, "CLOUD_PCLSAM%LBLAC2 ", CLOUD_PCLSAM%LBLAC2
		PRINT *, "CLOUD_D4S_LAYER ", CLOUD_D4S_LAYER
!		PRINT *, "CLOUD_TRANSMISSIVITY_D4S ", CLOUD_TRANSMISSIVITY_D4S
		PRINT *, "TEMPERATURE_PROFILE_IN ", TEMPERATURE_PROFILE_IN
!		PRINT *, "H2O_PROFILE ", H2O_PROFILE
!		PRINT *, "CO2_PROFILE ", CO2_PROFILE
!		PRINT *, "O3_PROFILE ", O3_PROFILE
!		PRINT *, "CO_PROFILE ", CO_PROFILE
!		PRINT *, "CH4_PROFILE ", CH4_PROFILE
!		PRINT *, "SO2_PROFILE ", SO2_PROFILE
!		PRINT *, "HNO3_PROFILE ", HNO3_PROFILE
!		PRINT *, "N2O_PROFILE ", N2O_PROFILE
!		PRINT *, "ALT ", ALT
!		PRINT *, "EMIS ", EMIS
!		PRINT *, "RHOTHR_IN ", RHOTHR_IN
!		PRINT *, "RHOSUN_IN ", RHOSUN_IN
!		PRINT *, "BEGIN_LAYER ", BEGIN_LAYER
!		PRINT *, "SELECTED_SPECIES ", SELECTED_SPECIES
!		PRINT *, "LEN_SELECTED_SPECIES ", LEN_SELECTED_SPECIES
	ENDIF

	SELECTED_SPECIES_COPY = ""
	SELECTED_SPECIES_COPY(1:LEN_SELECTED_SPECIES) = SELECTED_SPECIES(1:LEN_SELECTED_SPECIES)

	IF (PRESENT(BEGIN_LAYER)) THEN	
		START_LAYER = BEGIN_LAYER 
	ELSE
		START_LAYER = 1
	ENDIF

	IF (PRESENT(SELECTED_SPECIES)) THEN
		DO_SPECIES = SELECT_SPECIES__DETERMINE_DO_SPECIES_ARRAY(SELECTED_SPECIES_COPY)
	ELSE
		DO_SPECIES = SELECT_SPECIES__DO_ALL_SPECIES()
	ENDIF

	IF (SARTA_INITIALIZED .EQV. .TRUE.) THEN
		! Only check for channel changes if you're starting fresh.
		IF ( (DO_SPECIES % DO_ALL) .AND. (START_LAYER .EQ. 1) ) THEN
			! Check for changes in channel selection
			IF (AIRS_SARTA_VARIABLES__HAS_SELCHAN_CHANGED(NCHAN, SELCHAN)) THEN
				PRINT *, "Channels have changed -- reallocating arrays"
				CALL AIRS_SARTA_VARIABLES__INITIALIZE(NCHAN, SELCHAN)
				IF (SARTA_INITIALIZED .EQV. .TRUE.) THEN
					CALPAR_INITIALIZED = .FALSE.
					SUNPAR_INITIALIZED = .FALSE.
					CALL CALT1__ALLOCATE(NCHN1)
					CALL CALT2__ALLOCATE(NCHN2)
					CALL CALT3__ALLOCATE(NCHN3)
					CALL CALT4__ALLOCATE(NCHN4)
					CALL CALT4_SOLAR__ALLOCATE(NCHN4)
					CALL CALT5__ALLOCATE(NCHN5)
					CALL CALT5_SOLAR__ALLOCATE(NCHN5)
					CALL CALT6__ALLOCATE(NCHN6)
					CALL CALT6_SOLAR__ALLOCATE(NCHN6)
					CALL CALT7__ALLOCATE(NCHN7)
					CALL CALT7_SOLAR__ALLOCATE(NCHN7)
					PRINT *, "SARTA re-initialized for change in channel selection"
				ELSE
					PRINT *, "ERROR in re-initializing SARTA for change in channel selection"
					STOP
				ENDIF
			ENDIF
		ENDIF
	ELSE ! SARTA has not been initialized
		CALL AIRS_SARTA_VARIABLES__INITIALIZE(NCHAN, SELCHAN)
		IF (SARTA_INITIALIZED .EQV. .TRUE.) THEN
			CALL CALT1__ALLOCATE(NCHN1)
			CALL CALT2__ALLOCATE(NCHN2)
			CALL CALT3__ALLOCATE(NCHN3)
			CALL CALT4__ALLOCATE(NCHN4)
			CALL CALT4_SOLAR__ALLOCATE(NCHN4)
			CALL CALT5__ALLOCATE(NCHN5)
			CALL CALT5_SOLAR__ALLOCATE(NCHN5)
			CALL CALT6__ALLOCATE(NCHN6)
			CALL CALT6_SOLAR__ALLOCATE(NCHN6)
			CALL CALT7__ALLOCATE(NCHN7)
			CALL CALT7_SOLAR__ALLOCATE(NCHN7)
			PRINT *, "SARTA initialized"
		ELSE
			PRINT *, "ERROR in initializing SARTA"
			STOP
		ENDIF
	ENDIF

	IF (CALPAR_INITIALIZED .EQV. .FALSE.) THEN
		CALL CALPAR__INITIALIZE
		IF (CALPAR_INITIALIZED .EQV. .TRUE.) THEN
			PRINT *, "CALPAR initialized"
		ELSE
			PRINT *, "ERROR in initializing CALPAR"
			STOP
		ENDIF
	ENDIF

	IF (SUNPAR_INITIALIZED .EQV. .FALSE.) THEN
		CALL SUNPAR__INITIALIZE
		IF (SUNPAR_INITIALIZED .EQV. .TRUE.) THEN
			PRINT *, "SUNPAR initialized"
		ELSE
			PRINT *, "ERROR in initializing SUNPAR"
			STOP
		ENDIF
	ENDIF
	

	!--------------------------
	! Assign the I/O unit number
	!--------------------------
!	IOUN=11  ! not needed in this version FWI 4/26/18


	! Set satellite altitude to default.
	SALT = XSALT  	


	! This is done in AIRS_SARTA_VARIABLES and need not be repeated.  FWI 1/11/17
	!-----------------------------------------------
	! All channels from sets 1, 2, and 3 are to use a
	! fake effective sun angle layer-to-space trans
	!-----------------------------------------------
!	NFAKE=0

!	DO I=1,NCHN1
!		NFAKE=NFAKE + 1
!		INDFAK(NFAKE)=INDCHN( CLIST1(I) )
!	ENDDO

!	DO I=1,NCHN2
!		NFAKE=NFAKE + 1
!		INDFAK(NFAKE)=INDCHN( CLIST2(I) )
!	ENDDO
 
!	DO I=1,NCHN3
!		NFAKE=NFAKE + 1
!		INDFAK(NFAKE)=INDCHN( CLIST3(I) )
!	ENDDO


	!-------------------------------------
	! Determine bottom layer, CO2, & angles
	!-------------------------------------
	IF ( (DO_SPECIES % DO_ALL) .AND. (START_LAYER .EQ. 1) ) THEN
		CALL GETBOT__GETBOT(NLAY, PLEV, PSURF, LBOT, BLMULT)
	ELSE
		CALL GETBOT__GET_SAVED_GETBOT(LBOT, BLMULT)
	ENDIF

	IF ( (DO_SPECIES % DO_ALL) .AND. (START_LAYER .EQ. 1) ) THEN
		CALL SECANG_MODULE__CALCULATE_SECANG(SATZEN, SATANG, SUNANG, &
			SALT, ALT, LBOT)
	ENDIF

	!
	!! Calc the fractional bottom layer air temperature
	!!
	!! TEMPERATURE_PROFILE(LBOT)=TEMPERATURE_PROFILE(LBOT-1) + BLMULT*( TEMPERATURE_PROFILE(LBOT) - TEMPERATURE_PROFILE(LBOT-1) )
	!! Above line commented out & replaced by Scott Hannon, 24 July 2003.
	!! Mistakenly treats T at the center of the layer above as T at the
	!! bottom of the layer above.
	!!
	!
	!! CO2 profile switch
	! IF (ICO2 .LT. 1) THEN
	! 	LCO2=.FALSE.
	! ELSE
	!	LCO2=.TRUE.
	! ENDIF
	!! N2O profile switch
	! IF (IN2O .LT. 1) THEN
	!	LN2O=.FALSE.
	! ELSE
	!	LN2O=.TRUE.
	! ENDIF
	!! SO2 profile switch
	! IF (ISO2 .LT. 1) THEN
	!	LSO2=.FALSE.
	! ELSE
	!	LSO2=.TRUE.
	! ENDIF
	!! HNO3 profile switch
	! IF (IHNO3 .LT. 1) THEN
	!	LHNO3=.FALSE.
	! ELSE
	!	LHNO3=.TRUE.
	! ENDIF
	!

	! Make local copy of TEMPERATURE_PROFILE_IN that can be modified ! FWI 12/1/17
	TEMPERATURE_PROFILE(:) = TEMPERATURE_PROFILE_IN(:) 

	! Make local copy of RHOTHR and RHOSUN that can be modified ! FWI 12/1/17
	RHOTHR(:) = RHOTHR_IN(:)
	RHOSUN(:) = RHOSUN_IN(:)

	IF (PTYPE .EQ. AIRSLAY) THEN
		! Copy pseudo level temperatures to another array
		DO I = 1, LBOT
			TPSEUD(I)=TEMPERATURE_PROFILE(I)
		ENDDO
		! Convert temperatures
		CALL MEAN_T(LBOT, PLEV, PSURF, TPSEUD, TEMPERATURE_PROFILE)
	ELSE
		! Calc mean pressure for bottom fractional layer
		RJUNK1 = ( PSURF - PLEV(LBOT) )/LOG( PSURF/PLEV(LBOT) )
		! Do interpolation for fractional bottom layer mean temperature
		! assuming T is in linear in log(P)
		RJUNK2=( TEMPERATURE_PROFILE(LBOT) - TEMPERATURE_PROFILE(LBOT-1) ) / LOG( PLAY(LBOT)/PLAY(LBOT-1) ) ! slope
		TEMPERATURE_PROFILE(LBOT) = RJUNK2*LOG( RJUNK1/PLAY(LBOT-1) ) + TEMPERATURE_PROFILE(LBOT - 1)
	ENDIF


!	! Check satellite elevation
!	IF (SALT .GT. 0.0) THEN
!		! Warn and use default if invalid
!		IF (SALT .LT. XSALT-50 .OR. SALT .GT. XSALT+50) THEN
!			WRITE(IOINFO,1020) SALT, XSALT
!			1020 FORMAT('Warning! replacing invalid input satellite altitude ', 1PE11.4,' with default ',1PE11.4,' km')
!			SALT=XSALT
!		ENDIF
!	ELSE
!		SALT=XSALT
!	ENDIF
!
!	! Convert SATZEN or SATANG to viewing angle
!	IF (SATZEN .GE. 0.0 .AND. SATZEN .LT. 63.0) THEN
!		! Convert zenith angle at surface to view angle at satellite
!		SVA=SACONV( SATZEN, SALT*1000.0 )/CONV
!	ELSE
!		! Check if scan angle is valid
!		IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
!			! View angle should be within a few degrees of scan angle
!			SVA=ABS( SATANG )
!		ELSE
!			WRITE(IOERR,1030) SATZEN, SATANG
!			1030 FORMAT('Error! invalid angles for SATZEN ',1PE11.4,' and SATANG ',E11.4) 
!			STOP
!		ENDIF
!	ENDIF
!
!	ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
!	IF (SVA .GT. ANGMAX) THEN
!		! Truncate angle if too big
!		WRITE(IOINFO, 1040) SVA
!		1040 FORMAT('Warning! Profile: truncating view angle ',1PE11.4,' to 53 degrees')
!		SVA=ANGMAX
!	ENDIF
!
!	! Convert from satellite to earth viewing angle (in radians)
!	DO L=1,LBOT
!		EVA=VACONV(SVA, SALT, ALT(L))
!		SECANG(L)=1.0E+0/COS(EVA)
!		! for testing
!		! SECANG(L)=SVA
!	ENDDO
!
!	! Calc total sun angle secant
!	DOSUN=.FALSE.
!	IF (SUNANG .GE. 0.0 .AND. SUNANG .LT. 89.9) DOSUN=.TRUE.
!		IF (DOSUN) THEN
!			SUNCOS=COS(CONV*SUNANG)
!			SZALAY=SACONV(SUNANG,ALT(1))
!			SCOS1=COS(SZALAY)
!			RJUNK2=SECANG(LBOT) + 1.0/SUNCOS ! Total secant
!
!			! Calc non-unity fudge factor if total secant > 9
!			IF (RJUNK2 .GT. 9.0) THEN
!				! fudge factor = true_total_secant/calc_total_secant
!				SUNFDG=RJUNK2/9.0
!				! truncated solar angle to use to calc SECSUN
!				RJUNK1=ACOS( 1.0/(9.0 - SECANG(LBOT)) )/CONV
!			ELSE
!				SUNFDG=1.0
!				RJUNK1=SUNANG
!			ENDIF
!			! Should I change SUNFDG to SUNFDG(MAXLAY)?
!			DO L=1,LBOT
!				SZALAY=SACONV(RJUNK1,ALT(L))
!				SECSUN(L)=SECANG(L) + 1.0E+0/COS(SZALAY)
!			ENDDO
!	ENDIF


	!
	! Calculate the fast trans predictors
	!
	CALL CALPAR__CALCULATE_PREDICTORS(LBOT, &
		TEMPERATURE_PROFILE, CO2_PROFILE, H2O_PROFILE, O3_PROFILE, &
		CO_PROFILE, CH4_PROFILE, SO2_PROFILE, HNO3_PROFILE, N2O_PROFILE, &
		SECANG,   LAT,    FX,   &
		LCO2,  LN2O,  LSO2, LHNO3, LCO2PM, &  
		CO2PPM, CO2TOP, FIXMUL, &
		H2O_CONTINUUM_PRED_PTR, &
        FIXED_PRED1_PTR, FIXED_PRED2_PTR, FIXED_PRED3_PTR, FIXED_PRED4_PTR, &
		FIXED_PRED5_PTR, FIXED_PRED6_PTR, FIXED_PRED7_PTR, &
        H2O_PRED1_PTR,   H2O_PRED2_PTR,   H2O_PRED3_PTR,   H2O_PRED4_PTR, &
		H2O_PRED5_PTR,   H2O_PRED6_PTR,   H2O_PRED7_PTR, &
        O3_PRED1_PTR,    O3_PRED2_PTR,                     O3_PRED4_PTR, &
		O3_PRED5_PTR,    O3_PRED6_PTR,    O3_PRED7_PTR, &
        CH4_PRED3_PTR,   CO_PRED4_PTR,    TRACEGAS_PRED_PTR, &
		CO2MLT, SO2MLT, HNOMLT, N2OMLT, &
		START_LAYER, DO_SPECIES ) 

	!
	! Calculate the OPTRAN H2O predictors
	!
	CALL CALOWP ( LBOT, H2O_PROFILE, TEMPERATURE_PROFILE, SECANG, WAZOP, WAVGOP, &
		WAANG, LOPMIN, LOPMAX, LOPUSE, H2OPRD, LOPLOW, DAOP )

	!
	! Calculate the layer transmittances
	!

	! Calculate TAU for set 1 thru 7

	IF (NCHN1 .GT. 0) THEN
		CALL CALT1__CALCULATE_OD( INDCHN,  LBOT,   NCHN1, CLIST1,  COEF1, &
			FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED1_PTR, H2O_PRED1_PTR, O3_PRED1_PTR, TRACEGAS_PRED_PTR, &
			INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
			INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, &
			INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX, LOPLOW, &
			LOPUSE,   WAOP,   DAOP, WAANG,     TAU,   TAUZ, START_LAYER, DO_SPECIES)
	ENDIF


	IF (NCHN2 .GT. 0) THEN
		CALL CALT2__CALCULATE_OD( INDCHN,   LBOT,  NCHN2, CLIST2,  COEF2, &
			FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED2_PTR, O3_PRED2_PTR, H2O_PRED2_PTR, TRACEGAS_PRED_PTR, &
			INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
			INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, TAU, TAUZ, START_LAYER, DO_SPECIES)
	ENDIF

	IF (NCHN3 .GT. 0) THEN
		CALL CALT3__CALCULATE_OD( INDCHN,   LBOT,  NCHN3, CLIST3,  COEF3, &
			FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED3_PTR, CH4_PRED3_PTR, H2O_PRED3_PTR, TRACEGAS_PRED_PTR, &
			INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT, &
			INDN2O, COFN2O, N2OMLT, &
			INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX, LOPLOW, LOPUSE, &
			WAOP,   DAOP,    WAANG, &
			TAU, TAUZ, START_LAYER, DO_SPECIES)
	ENDIF

	IF (NCHN4 .GT. 0) THEN
		CALL CALT4__CALCULATE_OD( INDCHN,   LBOT,  NCHN4, CLIST4, COEF4, &
			FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED4_PTR, CO_PRED4_PTR, O3_PRED4_PTR, H2O_PRED4_PTR, &
			TRACEGAS_PRED_PTR, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
			TAU, TAUZ, START_LAYER, DO_SPECIES)
	ENDIF

	IF (NCHN5 .GT. 0) THEN
		CALL CALT5__CALCULATE_OD(  INDCHN,   LBOT,  NCHN5, CLIST5, COEF5, &
			FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED5_PTR, H2O_PRED5_PTR, O3_PRED5_PTR, &
			TRACEGAS_PRED_PTR, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
			TAU, TAUZ, START_LAYER, DO_SPECIES )
	ENDIF

	IF (NCHN6 .GT. 0) THEN
		CALL CALT6__CALCULATE_OD(  INDCHN,   LBOT,  NCHN6, CLIST6, COEF6, &
			FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED6_PTR, H2O_PRED6_PTR, O3_PRED6_PTR, TRACEGAS_PRED_PTR, &
			INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
			INDN2O, COFN2O, N2OMLT, &
			TAU, TAUZ, START_LAYER, DO_SPECIES )
	ENDIF

	IF (NCHN7 .GT. 0) THEN
		CALL CALT7__CALCULATE_OD(INDCHN,    LBOT,  NCHN7, CLIST7, COEF7, &
			FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED7_PTR, H2O_PRED7_PTR, O3_PRED7_PTR, &
			TRACEGAS_PRED_PTR, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
			TAU, TAUZ, START_LAYER, DO_SPECIES)
	ENDIF

	IF (DOSUN) THEN
		!
		! Calculate the fast trans predictors *for sun*
		!

		CALL SUNPAR__CALCULATE_PREDICTORS ( LBOT, &
			TEMPERATURE_PROFILE,  H2O_PROFILE,  O3_PROFILE,  CO_PROFILE, &
			SECSUN, H2O_CONTINUUM_PRED_PTR, & 
			FIXED_PRED4_PTR, FIXED_PRED5_PTR, FIXED_PRED6_PTR, FIXED_PRED7_PTR, &
			H2O_PRED4_PTR, H2O_PRED5_PTR, H2O_PRED6_PTR, H2O_PRED7_PTR, &
			O3_PRED4_PTR, O3_PRED5_PTR, O3_PRED6_PTR, O3_PRED7_PTR, &
			CO_PRED4_PTR, TRACEGAS_PRED_PTR, &
			START_LAYER, DO_SPECIES )

		!
		! Calculate the layer transmittances *for sun*
		!

		!  Calc fake TAUZSN for sets 1, 2, and 3

		IF ( (NCHN1 .GT. 0) .OR. (NCHN2 .GT. 0) .OR. (NCHN3 .GT. 0) ) THEN
			CALL FAKETZ( NFAKE, INDFAK, LBOT, TAUZ, SECANG, &
				SECSUN, TAUZSN)
		ENDIF


		! Calculate TAUZSN for sets 4 thru 7

		IF (NCHN4 .GT. 0) THEN
			CALL CALT4_SOLAR__CALCULATE_OD(INDCHN,   LBOT,  NCHN4, CLIST4, &
				 COEF4, FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED4_PTR, CO_PRED4_PTR, O3_PRED4_PTR, H2O_PRED4_PTR, &
				TRACEGAS_PRED_PTR, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
				TAUZSN, TAUZSN, START_LAYER, DO_SPECIES ) 
				!^^^^^^  ^^^^^^ 
				! dummy   actual 
		ENDIF

		IF (NCHN5 .GT. 0) THEN
			CALL CALT5_SOLAR__CALCULATE_OD(INDCHN,   LBOT,  NCHN5, CLIST5, &
				 COEF5, FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED5_PTR, H2O_PRED5_PTR, O3_PRED5_PTR, &
				TRACEGAS_PRED_PTR, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
				TAUZSN, TAUZSN, START_LAYER, DO_SPECIES ) 
		ENDIF
 
		IF (NCHN6 .GT. 0) THEN
			CALL CALT6_SOLAR__CALCULATE_OD(INDCHN,   LBOT,  NCHN6, CLIST6, &
				 COEF6, FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED6_PTR, H2O_PRED6_PTR, O3_PRED6_PTR, &
				TRACEGAS_PRED_PTR, INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT, &
				INDN2O, COFN2O, N2OMLT, TAUZSN, TAUZSN, START_LAYER, DO_SPECIES )  
 		ENDIF

		IF (NCHN7 .GT. 0) THEN
			CALL CALT7_SOLAR__CALCULATE_OD(INDCHN,   LBOT,  NCHN7, CLIST7, &
				 COEF7, FIXMUL, H2O_CONTINUUM_PRED_PTR, FIXED_PRED7_PTR, H2O_PRED7_PTR, O3_PRED7_PTR, &
				TRACEGAS_PRED_PTR, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT, &
				TAUZSN, TAUZSN, START_LAYER, DO_SPECIES  ) 
		ENDIF

		IF (SUNFDG .GT. 1.0001) THEN
			DO I=1,NCHAN
				DO L=1,LBOT
					TAUZSN(L,I)=TAUZSN(L,I)*SUNFDG
				ENDDO
			ENDDO
		ENDIF

	ENDIF

	SUNFAC=SUNCOS*PI*(RADSUN/DISTES)**2

	! Note: PI*(RADSUN/DISTES)^2 = solid angle [steradians] of
	! the sun as seen from Earth for the case DISTES >> RADSUN.

	!
	! Clear scene
	!
	IF ((CLOUD_D4S_LAYER .LT. 0) .AND. (CLOUD_PCLSAM%CFRAC1 .LE. 0.) .AND. (CLOUD_PCLSAM%CFRAC2 .LE. 0.) ) THEN

		DO I = 1, MXCHAN
			IF (INDCHN(I) .GT. 0) THEN
				! Radiation constants for current channel
				C1V3=C1*(FREQCHAN(I)**3)
				C2V=C2*FREQCHAN(I)

				! Calculate Planck & clear airs trans for full layers
				DO L=1, LBOT-1
					! Fix C1V3 and C2V so it doesn't have to be recomputed every time ?  FWI 9/1/17           
					RPLNCK(L)=C1V3/ ( EXP( C2V/TEMPERATURE_PROFILE(L) ) - 1.0 )
					TRANL(L)=QIKEXP( -TAU(L,I) )
				ENDDO
				! Note: TEMPERATURE_PROFILE(LBOT) already adjusted for bottom fractional layer
				RPLNCK(LBOT)=C1V3/( EXP( C2V/TEMPERATURE_PROFILE(LBOT) ) - 1.0 )

				! Calculate clear airs trans for bottom fractional layer

				!  Why are these lines repeated? What am I missing?
				RJUNK1=-TAU(LBOT,I)*BLMULT
				TRANL(LBOT)=QIKEXP( RJUNK1 )
				!TRANL(LBOT)=QIKEXP( RJUNK1 )
				TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LBOT-1,I) )
				TRANS(I)=QIKEXP( BLMULT*(TAUZSN(LBOT-1,I)-TAUZSN(LBOT,I)) - TAUZSN(LBOT-1,I) )

				! Planck for surface
				RSURFE=EMIS(I)*C1V3/( EXP( C2V/TSURF ) - 1.0 )
															
				CALL CALRAD0( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG, &
					 TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN, &
					RHOTHR, LABOVE, COEFF, RAD(I) ) 
			ENDIF
		ENDDO

		IF (DOSUN) THEN
			CALL CALNTE ( INDCHN, TEMPERATURE_PROFILE, SUNCOS, SCOS1, SECANG(1), &
				NCHNTE, CLISTN, COEFN, CO2TOP, RAD )
		ENDIF

		RETURN
	ENDIF
	


	! Note, GETCLD and SETEMS calls removed as cloud data and emissivity/reflectivity
	! are now input in the call to this subroutine.
	! FWI - 12/4/2107

	! Test each cloud modeling regimens.  If all show no cloud, then calculate for clear sky.
	! FWI 12/05/2017

	! If CLOUD_D4S_LAYER > 0, we're using pre-computed delta-4-stream optical depths
	! and transmissivity.  
	IF (CLOUD_D4S_LAYER .GE. 0) THEN

		DO I = 1, MXCHAN
			IF (INDCHN(I) .GT. 0) THEN
				!PRINT *, "MAIN I INDCHN(I) FREQCHAN(I) EMIS(I)", I, INDCHN(I), FREQCHAN(I), EMIS(I)

				! Radiation constants for current channel
				!C1V3=C1*(FREQ(I)**3)
				!C2V=C2*FREQ(I)
				C1V3=C1*(FREQCHAN(I)**3)
				C2V=C2*FREQCHAN(I)

				! Calculate Planck & clear airs trans for full layers
				DO L=1, LBOT-1
					! Fix C1V3 and C2V so it doesn't have to be recomputed every time ?  FWI 9/1/17           
					RPLNCK(L)=C1V3/ ( EXP( C2V/TEMPERATURE_PROFILE(L) ) - 1.0 )
					TRANL(L)=QIKEXP( -TAU(L, I) )
				ENDDO
				! Note: TEMPERATURE_PROFILE(LBOT) already adjusted for bottom fractional layer
				RPLNCK(LBOT)=C1V3/( EXP( C2V/TEMPERATURE_PROFILE(LBOT) ) - 1.0 )

				! Calculate clear airs trans for bottom fractional layer
				RJUNK1=-TAU(LBOT,I)*BLMULT
				TRANL(LBOT)=QIKEXP( RJUNK1 )
				TRANL(LBOT)=QIKEXP( RJUNK1 )
				TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LBOT-1,I) )
				TRANS(I)=QIKEXP( BLMULT*(TAUZSN(LBOT-1,I)-TAUZSN(LBOT,I)) - TAUZSN(LBOT-1,I) )

				! Planck for surface
				RSURFE=EMIS(I)*C1V3/( EXP( C2V/TSURF ) - 1.0 )
															
				CALL CALRAD_D4S( &
					 DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG, &
					 TRANL, TRANZ, CLOUD_TRANSMISSIVITY_D4S, CLOUD_D4S_LAYER, SUNFAC, HSUN, TRANS, RHOSUN, & 
					RHOTHR, LABOVE, COEFF, RAD(I) )
			ENDIF
		END DO

		!-----------------
		! Calculate non-LTE
		!-----------------
		! comment: the nonLTE calculation does not consider cloud effects,
		! but clouds are generally below the altitude where nonLTE occurs.
		IF (DOSUN) THEN
			CALL CALNTE ( INDCHN, TEMPERATURE_PROFILE, SUNCOS, SCOS1, SECANG(1), &
				NCHNTE, CLISTN, COEFN, CO2TOP, RAD )
		ENDIF

		RETURN		

	ENDIF


	! Set exclusive cloud fractions (moved from getcld.f)

	IF ((CLOUD_PCLSAM%CFRAC12 .LE. CLOUD_PCLSAM%CFRAC1) .AND. (CLOUD_PCLSAM%CFRAC12 .LE. CLOUD_PCLSAM%CFRAC2)) THEN
		CFRA1X=CLOUD_PCLSAM%CFRAC1 - CLOUD_PCLSAM%CFRAC12
		CFRA2X=CLOUD_PCLSAM%CFRAC2 - CLOUD_PCLSAM%CFRAC12
	ELSE
		CFRA1X=CLOUD_PCLSAM%CFRAC1 ! FWI 2/21/18
		CFRA2X=CLOUD_PCLSAM%CFRAC2 ! FWI 2/21/18
	ENDIF 

		
	FCLEAR = 1. - (CLOUD_PCLSAM%CFRAC1 + CLOUD_PCLSAM%CFRAC2 - CLOUD_PCLSAM%CFRAC12)  ! FWI 4/5/18
	
	IF (FCLEAR .LT. 0.001) FCLEAR = 0.
	IF (CFRA1X .LT. 0.001) CFRA1X = 0.
	IF (CFRA2X .LT. 0.001) CFRA2X = 0.

	!PRINT *, "CLOUD_PCLSAM%CFRAC1 ", CLOUD_PCLSAM%CFRAC1, "CFRAC2", CLOUD_PCLSAM%CFRAC2, "CFRAC12 ", CLOUD_PCLSAM%CFRAC12
	!PRINT *, "FCLEAR ", FCLEAR, "CFRA1X ", CFRA1X, "CFRA2X", CFRA2X
	!PRINT *, ""

	! Check and prepare (top) cloud1
	IF (CLOUD_PCLSAM%CFRAC1 .GT. 0.0) THEN
		IF (CLOUD_PCLSAM%LBLAC1 .EQ. 1) THEN
			CALL BKPREP(1, CLOUD_PCLSAM%CTYPE1, CLOUD_PCLSAM%CFRAC1, CLOUD_PCLSAM%CPRTOP1, &
				LBOT, PSURF, PLEV, PLAY, TEMPERATURE_PROFILE, LCTOP1, TCTOP1, &
				TEMPC1, CLRT1)
			IF (CLOUD_PCLSAM%CSTEMP1 .GT. 0.0) TCTOP1 = CLOUD_PCLSAM%CSTEMP1
		ELSE
			! Determine which lookup table to use
			CALL GETMIE(CLOUD_PCLSAM%CTYPE1, MIETYP, INDMI1, IERR1)
			! Prepare selected lookup table for given cpsize
			CALL CCPREP( NCHAN, INDCHN, LBOT, INDMI1, MIENPS, &
				CLOUD_PCLSAM%CNGWAT1, CLOUD_PCLSAM%CPSIZE1, &
				CLOUD_PCLSAM%CPRTOP1, CLOUD_PCLSAM%CPRBOT1, &
				PLEV, TEMPERATURE_PROFILE, SECANG, &
				SECSUN, MIEPS, MIEABS, MIEEXT, MIEASY, LCBOT1, LCTOP1, &
				CLRB1, CLRT1, TCBOT1, TCTOP1, MASEC1, MASUN1, &
				CFRCL1, G_ASY1, NEXTO1, NSCAO1 )
		ENDIF
	ENDIF
	
	! Check and prepare (bottom) cloud2
	IF (CLOUD_PCLSAM%CFRAC2 .GT. 0.0) THEN
		IF (CLOUD_PCLSAM%LBLAC2 .EQ. 1) THEN
			CALL BKPREP( 2, CLOUD_PCLSAM%CTYPE2, CLOUD_PCLSAM%CFRAC2, CLOUD_PCLSAM%CPRTOP2, &
				LBOT, PSURF, PLEV, PLAY, TEMPERATURE_PROFILE, LCTOP2, TCTOP2, &
				TEMPC2, CLRT2)
			IF (CLOUD_PCLSAM%CSTEMP2 .GT. 0.0) TCTOP2=CLOUD_PCLSAM%CSTEMP2
		ELSE
			! Determine which lookup table to use
			CALL GETMIE(CLOUD_PCLSAM%CTYPE2,MIETYP,INDMI2,IERR2)
			! Prepare lookup data for cloud2
			CALL CCPREP( NCHAN, INDCHN, LBOT, INDMI2, MIENPS, &
				CLOUD_PCLSAM%CNGWAT2, CLOUD_PCLSAM%CPSIZE2, &
				CLOUD_PCLSAM%CPRTOP2, CLOUD_PCLSAM%CPRBOT2, &
				PLEV, TEMPERATURE_PROFILE, SECANG, &
				SECSUN, MIEPS, MIEABS, MIEEXT, MIEASY, LCBOT2, LCTOP2, &
				CLRB2, CLRT2, TCBOT2, TCTOP2, MASEC2, MASUN2, &
				CFRCL2, G_ASY2, NEXTO2, NSCAO2 )
		ENDIF
	ELSE
		! Safe default for non-existant cloud2
		LCTOP2=1
	ENDIF

	!! this block for testing only
	!! PROF.udef(19)=TCTOP1
	!! PROF.udef(20)=TCTOP2
	!!ccccccccccccccccccccccccccccccccccc
	!----------------------
	! Loop over the channels
	! ----------------------

	COSDAZ = COS((SUN_AZIMUTH - SAT_AZIMUTH) * PI_OVER_180)  

	!DO I=1, NCHAN
	DO I=1, MXCHAN

		IF (INDCHN(I) .GT. 0) THEN
			! Radiation constants for current channel
			!C1V3=C1*(FREQ(I)**3)
			!C2V=C2*FREQ(I)
			C1V3=C1*(FREQCHAN(I)**3)
			C2V=C2*FREQCHAN(I)

			! Calculate Planck & clear airs trans for full layers
			DO L=1, LBOT-1
				! Fix C1V3 and C2V so it doesn't have to be recomputed every time ?  FWI 9/1/17           
				RPLNCK(L)=C1V3/ ( EXP( C2V/TEMPERATURE_PROFILE(L) ) - 1.0 )
				TRANL(L)=QIKEXP( -TAU(L,I) )
			ENDDO
			! Note: TEMPERATURE_PROFILE(LBOT) already adjusted for bottom fractional layer
			RPLNCK(LBOT)=C1V3/( EXP( C2V/TEMPERATURE_PROFILE(LBOT) ) - 1.0 )

			! Calculate clear airs trans for bottom fractional layer
			RJUNK1=-TAU(LBOT,I)*BLMULT
			TRANL(LBOT)=QIKEXP( RJUNK1 )
			TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LBOT-1,I) )
			TRANS(I)=QIKEXP( BLMULT*(TAUZSN(LBOT-1,I)-TAUZSN(LBOT,I)) - TAUZSN(LBOT-1,I) )

			! Planck for surface
			RSURFE=EMIS(I)*C1V3/( EXP( C2V/TSURF ) - 1.0 )  ! FWI 2/21/18

			! Calculate clear radiance
			IF (FCLEAR .GT. 0.0) THEN
				CALL CALRAD0( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG, &
					 TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN, &
					RHOTHR, LABOVE, COEFF, RAD0 ) 
			ELSE
				RAD0=0.0
			ENDIF

			! Store original values
			VSTORE(1)=TRANL(LCTOP2)
			VSTORE(2)=TRANZ(I)
			VSTORE(3)=TRANS(I)
			VSTORE(4)=RHOTHR(I)
			VSTORE(5)=RHOSUN(I)
			VSTORE(6)=RPLNCK(LCTOP2)
			! Updates for new surface if bottom cloud2 is black
			IF (CLOUD_PCLSAM%CFRAC2 .GT. 0.0 .AND. CLOUD_PCLSAM%LBLAC2 .EQ. 1) THEN
				RJUNK1=-TAU(LCTOP2,I)*CLRT2
				TRANL(LCTOP2)=QIKEXP( RJUNK1 )
				TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LCTOP2-1,I) )
				TRANS(I)=QIKEXP( CLRT2*(TAUZSN(LCTOP2-1,I)-  &
					TAUZSN(LCTOP2,I)) - TAUZSN(LCTOP2-1,I) )
!				RSURFC=CEMIS2(I)*C1V3/( EXP( C2V/TCTOP2 ) - 1.0 )
				RSURFC=CLOUD_PCLSAM%CEMIS2(I)*C1V3/( EXP( C2V/TCTOP2 ) - 1.0 )
!				RHOTHR(I)=CRHOT2(I)
!				RHOSUN(I)=CRHOS2(I)
				RHOTHR(I)=CLOUD_PCLSAM%CRHO2(I)
				RHOSUN(I)=CLOUD_PCLSAM%CRHO2(I)
				RPLNCK(LCTOP2)=C1V3/( EXP( C2V/TEMPC2 ) - 1.0 )
			ENDIF

			! Calculate bottom cloud2 radiance
			IF (CFRA2X .GT. 0.0) THEN
				IF (CLOUD_PCLSAM%LBLAC2 .EQ. 1) THEN
					CALL CALRAD0( DOSUN, I, LCTOP2, RPLNCK, RSURFC, SECANG, &
						 TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN, & 
						RHOTHR, LABOVE, COEFF, RADC2 )
				ELSE
					CALL CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG, &
							TAU,  TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN, & 
						 RHOTHR, LABOVE, COEFF, CFRCL2, MASEC2, MASUN2, COSDAZ, &
						 NEXTO2, NSCAO2, G_ASY2(I), LCTOP2, LCBOT2, RADC2 )
				 ENDIF
			ELSE
				RADC2=0.0
			ENDIF

			! Calculate combined cloud1+cloud2 radiance
			IF (CLOUD_PCLSAM%CFRAC12 .GT. 0.0) THEN
				IF (CLOUD_PCLSAM%LBLAC2 .EQ. 1) THEN
					CALL CALRAD1( DOSUN, I, LCTOP2, RPLNCK, RSURFC, SECANG, &
						TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN, &
						RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, COSDAZ, &
						NEXTO1, NSCAO1, G_ASY1(I), LCTOP1, LCBOT1, RADC12 )
				!ELSE  ! FWI 2/21/18
				ELSE IF (CLOUD_PCLSAM%CFRAC2 .GT. 0.0) THEN  ! FWI 2/21/18
					CALL CALRAD2( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG, &
						TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN, & 
						RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, NEXTO1, &
						NSCAO1, G_ASY1(I), LCTOP1, LCBOT1, CFRCL2, MASEC2, MASUN2, &
						COSDAZ, NEXTO2, NSCAO2, G_ASY2(I), LCTOP2, LCBOT2, RADC12 )
				ENDIF
			ELSE
				RADC12=0.0
			ENDIF

			! Restore original values
			TRANL(LCTOP2)=VSTORE(1)
			TRANZ(I)=VSTORE(2)
			TRANS(I)=VSTORE(3)
			RHOTHR(I)=VSTORE(4)
			RHOSUN(I)=VSTORE(5)
			RPLNCK(LCTOP2)=VSTORE(6)
			! Updates for new surface if top cloud1 is black
			IF (CLOUD_PCLSAM%CFRAC1 .GT. 0.0 .AND. CLOUD_PCLSAM%LBLAC1 .EQ. 1) THEN
				RJUNK1=-TAU(LCTOP1,I)*CLRT1
				TRANL(LCTOP1)=QIKEXP( RJUNK1 )
				TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LCTOP1-1,I) )
				TRANS(I)=QIKEXP( CLRT1*(TAUZSN(LCTOP1-1,I)- &
					TAUZSN(LCTOP1,I)) - TAUZSN(LCTOP1-1,I) )
				!RSURFC=CEMIS1(I)*C1V3/( EXP( C2V/TCTOP1 ) - 1.0 )
				RSURFC=CLOUD_PCLSAM%CEMIS1(I)*C1V3/( EXP( C2V/TCTOP1 ) - 1.0 )
!				RHOTHR(I)=CRHOT1(I)
!				RHOSUN(I)=CRHOS1(I)
				RHOTHR(I)=CLOUD_PCLSAM%CRHO1(I)
				RHOSUN(I)=CLOUD_PCLSAM%CRHO1(I)
				RPLNCK(LCTOP1)=C1V3/( EXP( C2V/TEMPC1 ) - 1.0 )
			ENDIF

			! Calculate top cloud1 radiance
			IF (CFRA1X .GT. 0.0) THEN
				IF (CLOUD_PCLSAM%LBLAC1 .EQ. 1) THEN
					CALL CALRAD0( DOSUN, I, LCTOP1, RPLNCK, RSURFC, SECANG, &
						TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN, &
						RHOTHR, LABOVE, COEFF, RADC1 )
				ELSE
					CALL CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG, &
						TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN, &
						RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, COSDAZ, &
						NEXTO1, NSCAO1, G_ASY1(I), LCTOP1, LCBOT1, RADC1 )
				ENDIF
			ELSE
				RADC1=0.0
			ENDIF
	
			! Total the clear & various cloudy radiances
			RAD(I)=RAD0*FCLEAR + RADC1*CFRA1X + RADC2*CFRA2X + RADC12*CLOUD_PCLSAM%CFRAC12
		ENDIF 
	ENDDO ! channels

	!-----------------
	! Calculate non-LTE
	!-----------------
	! comment: the nonLTE calculation does not consider cloud effects,
	! but clouds are generally below the altitude where nonLTE occurs.
	IF (DOSUN) THEN
		CALL CALNTE ( INDCHN, TEMPERATURE_PROFILE, SUNCOS, SCOS1, SECANG(1), &
			NCHNTE, CLISTN, COEFN, CO2TOP, RAD )
	ENDIF

	RETURN

END SUBROUTINE AIRS_SARTA_PCLSAM_D4S

END
