!=======================================================================
!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS (Atmospheric Infra-Red Sounder)
!
!    INCFTC
!
!F90====================================================================


! ROUTINE NAME:
!    incFTC (include file rewritten as module)


! ABSTRACT:
!    Include file consisting of parameter statements to size various
!    arrays in the USEFAST related routines source code.


! CALL PROTOCOL:
!    none (include file)


! INPUT PARAMETERS:
!    none


! OUTPUT PARAMETERS:
!    none


! INPUT/OUTPUT PARAMETERS:
!    none


! RETURN VALUES:
!    none


! PARENT(S):
!    CALOKW
!    CALOWP
!    CALPAR
!    CALRAD
!    CALT1
!    CALT2
!    CALT3
!    CALT4
!    CALT5
!    CALT6
!    CALT7
!    FAKETZ
!    RDCOEF
!    RDLIST
!    RDPROF
!    RDSUN
!    SUNPAR
!    SARTA


! ROUTINES CALLED:
!    none


! FILES ACCESSED:
!    none


! COMMON BLOCKS
!    none


! DESCRIPTION:
!    Include file for the December 2005 100 layer AIRS fast
!    Stand Alone RTA (SARTA) code by L.L.Strow/S.Hannon.
!
!    Parameter statements for the FTC routines.


! ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
! Date	 Programmer     Comments
! ----------- -------------- -------------------------------------------
!  1 Dec 1994 Scott Hannon   Created
! 31 Jan 1997 Scott Hannon   Re-wrote for FWO+FOW+FMW+FCOW=Feb97 FTC
!  3 Sep 1997 Scott Hannon   Re-wrote for sets 1 - 7
! 30 Sep 1997 Scott Hannon   Added NCO2 and MXCHNC
! 26 Feb 1998 Scott Hannon   Added OPTRAN variables for water, and
!				changed both N1H2O & N3H2O from 13 to 11
! 23 Sep 1999 Scott Hannon   Change number of channel dimensions for
!				new Sep99 version of fast model.
!  5 Apr 2000 Scott Hannon   Added MXEMIS
!  4 Aug 2000 Scott Hannon   Changes values for use with testfast
! 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
! 23 Jan 2001 Scott Hannon   Update values of C1 & C2
! 15 Feb 2001 Scott Hannon   Add MAXPRO, CO2STD, IOERR, IOINFO,
!				MXGAS, CSARTA, and all filenames
! 24 Apr 2001 Scott Hannon   Add MXMIEA and FNMIEA
! 14 Aug 2001 Scott Hannon   Add FNMIEE and FNMIEG
! 21 Nov 2001 Scott Hannon   Add VSARTA, VSCOEF, & VCLOUD; remove CSARTA
! 12 Sep 2002 Scott Hannon   Updated for m135f (-13.5 um with fringes)
! 17 Dec 2002 Scott Hannon   Updated for revised(Dec02) m135f 
!  3 Jan 2003 Scott Hannon   Updated VSARTA for version 1.04
!  3 Feb 2004 Scott Hannon   Updated for "cc2"
! 19 Feb 2004 Scott Hannon   Add FNTMLT & update VSARTA for v1.05
! 07 Apr 2005 Scott Hannon   NFCOEF increased from 5 to 6 for v1.06
! 18 May 2005 Scott Hannon   update for HNO3 version
! 29 Jun 2005 Scott Hannon   "trace" version v1.07 with CO2,SO2,HNO3,N2O
! 13 Oct 2005 Scott Hannon   Add variables for non-LTE
! 22 Nov 2005 Scott Hannon   Replace set1,set2,CO2 coefs for new M12
! 30 Mar 2006 Scott Hannon   Change from separate MI#* for each cloud
!				type to a single MIE* 3-D array for all
!				types, and add MIETYP & NMIETY.
! 02 May 2007 Scott Hannon   Added XSALT
! 12 May 2009 Scott Hannon   Add VTUNNG; delete VCLOUD

! 29 Aug 2017 FWI	     Conversion to F90 module
! 30 Aug 2017 FWI	     Added in CLOUD type
! 12 Sep 2017 FWI        Changed dimensioning of emissivity/reflectivity 
!                        of land and clouds from MXEMIS to MXCHAN. 
!                        Emissivities etc are input to SARTA on a channel-by-channel
!                        basis, and not interpolated.

!END====================================================================
!
!-----------------------------------------------------------------------
!      IMPLICIT NONE
!-----------------------------------------------------------------------
! Note: having an "implicit none" in both the include file & the main
! source code will cause some compilers to complain.
!	IMPLICIT NONE


!-----------------------------------------------------------------------
!      INCLUDE FILES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      EXECUTABLE CODE
!-----------------------------------------------------------------------
!      none

MODULE INCFTC

	!	-----------------------------------------------------------------
	!	Assign SARTA version strings
	!	-----------------------------------------------------------------
	!	The version strings consists of 3 parts: version number, date,
	!	and comment.  The version date should be updated to the
	!	current date whenever any portion of the code is updated.  The
	!	version number consists of two parts; a major version to the
	!	left of the decimal point, and a minor version to the right.
	!	The major number should be incremented only when major changes
	!	have been made to the overall SARTA code.  The minor number
	!	should be incremented only when minor but non-trivial changes
	!	are made to the code.  Bug fixes should generally be handled
	!	with the version date, but a fix for a serious bug may warrant
	!	a change to the minor version number.
	!	See the "Doc/last_update.txt" file for a description of the
	!	changes associated with every change of VSARTA.
	!
	CHARACTER*40 VSARTA  ! SARTA source code version
	CHARACTER*40 VSCOEF  ! SARTA coefficient version
	CHARACTER*40 VTUNNG  ! optical depth tuning version
	! version template    '#.## YYYY-MM-DD <--------comment------->'
	PARAMETER( VSARTA = '1.08 2010-09-14 rtpV201 PCLSAM slab HG3')
	PARAMETER( VSCOEF = 'AIRS 2008-04-30 m140x CO2=370' )
	PARAMETER( VTUNNG = 'v6 standard; refprof N2O x1/1.04')

	!
	! VARIABLES
	!
	! Note: these should not be changed by the user
	!
	!
	! Constants and other data
	!
	REAL, PARAMETER :: PI  = 3.1415926 ! pi, circle circumference/diameter (3.1415926)
	REAL, PARAMETER :: RADSUN = 6.956E+8 ! radius of the sun (6.956E+8 m)
	REAL, PARAMETER :: C1 = 1.191042722E-8 ! radiation constant c1 (1.1911E-8  W/(m2.st.(cm-1)4)
	REAL, PARAMETER :: C2 = 1.4387752  ! radiation constant c2 (1.4387863 K/cm-1)
 
	! Previously used values; agrees w/JPL pre-Dec2000
	! PARAMETER(    C1 = 1.1910439E-8)  ! JPL value is 1E+3 bigger
	! PARAMETER(    C2 = 1.4387687)
	!
	! Current values (CODATA98 from NIST); agrees w/JPL Dec2000
	! PARAMETER(  C1 = 1.191042722E-8)  ! JPL value is 1E+3 bigger
	! PARAMETER(  C2 = 1.4387752)
	
	! Profile type -- copied from RTP parameters, FWI, 8/31/17

	INTEGER, PARAMETER :: LEVPRO = 0
	INTEGER, PARAMETER :: LAYPRO = 1
	INTEGER, PARAMETER :: AIRSLAY = 2
	

	REAL, PARAMETER :: CO2STD = 385.0 ! standard CO2 PPMV mixing ratio ! m130, m140, m150
	! PARAMETER( CO2STD = 370.0 )  ! m130x, m140x

	REAL, PARAMETER :: XSALT = 705.0 ! expected nominal satellite altitude (km)

	! CONV = pi/180 = degrees to radians conversion factor
	REAL, PARAMETER :: CONV = 1.7453292E-02


	!
	! Channels and layers other variables
	!
	INTEGER, PARAMETER :: MAXLAY = 100 ! # of layers (100)
	INTEGER, PARAMETER :: MAXLEV = MAXLAY + 1 ! # of levels !FWI 8/29/17
	INTEGER, PARAMETER :: NSET = 7 ! # of coefficient data sets (7)
	INTEGER, PARAMETER :: MXCHAN = 2834 ! max total # of channels (2378)
	INTEGER, PARAMETER :: NFCOEF = 6 ! # of downwelling thermal "F" factor coefs 
	! INTEGER, PARAMETER :: MXEMIS = 100 ! max # of input emis/rho data points
	INTEGER, PARAMETER :: MAXPRO = 25 ! max # of user specified profiles
	INTEGER, PARAMETER ::  MXGAS = 44 ! max # of gases in user profile

	!
	!
	! Variables for the coefficient sets
	!
	!

	!	
	! For set1 = FWO
	!	
	! Used in part by modules: 12, 11, 10, 9, 8, 7, 6, 5, 3, 4b, 4a
	INTEGER, PARAMETER :: MXCHN1 = 1461 ! max # of channels for set1 = FWO (1461)
	INTEGER, PARAMETER :: N1CON = 7 ! # of water con predictors/coefs for set1 (5)
	INTEGER, PARAMETER :: N1FIX = 8 ! # of "fixed" predictors/coefs for set1 (8)
	INTEGER, PARAMETER :: N1H2O = 11 ! # of water predictors/coefs for set1 (13)
	INTEGER, PARAMETER :: N1O3 = 5 ! # of ozone predictors/coefs for set1 (5)
	INTEGER, PARAMETER :: N1COEF = N1CON + N1FIX + N1H2O + N1O3 ! total # of coefs for set1

	!	
	! For set2 = FOW
	!	
	! Used in part by modules: 6, 5
	INTEGER, PARAMETER :: MXCHN2 = 325 ! max # of channels for set2 = FOW  (325)
	INTEGER, PARAMETER :: N2CON = 7    ! # of water con predictors/coefs for set2 (5)
	INTEGER, PARAMETER :: N2FIX = 8    ! # of "fixed" predictors/coefs for set2 (8)
	INTEGER, PARAMETER :: N2O3 = 10    ! # of ozone predictors/coefs for set2 (10)
	INTEGER, PARAMETER :: N2H2O = 11   ! # of water predictors/coefs for set2 (11)
	INTEGER, PARAMETER :: N2COEF = N2CON + N2FIX + N2O3 + N2H2O ! total # of coefs for set2


	!
	! For set3 = FMW
	!
	! Used in part by modules: 4d, 4c, 3
	INTEGER, PARAMETER :: MXCHN3 = 396 	! max # of channels for set3 = FMW  (396)
	INTEGER, PARAMETER :: N3CON 	= 7	! # of water con predictors/coefs for set3 (5)
	INTEGER, PARAMETER :: N3FIX 	= 8	! # of "fixed" predictors/coefs for set3 (8)
	INTEGER, PARAMETER :: N3CH4	= 9 	! # of methane predictors/coefs for set3 (9)
	INTEGER, PARAMETER :: N3H2O	= 11 	! # of water predictors/coefs for set3 (13)
	INTEGER, PARAMETER :: N3COEF = N3CON + N3FIX + N3CH4 + N3H2O  ! total # of coefs for set3


	!
	! For set4 = sun FCOW
	!
	! Used in part by modules: 2b
	INTEGER, PARAMETER :: MXCHN4 = 85 ! max # of channels for set4 = FCOW (85)
	INTEGER, PARAMETER :: N4CON = 7 ! # of water con predictors/coefs for set4 (5)
	INTEGER, PARAMETER :: N4FIX = 11 ! # of "fixed" predictors/coefs for set4 (11)
	INTEGER, PARAMETER :: N4CO = 11 ! # of CO predictors/coefs for set4 (11)
	INTEGER, PARAMETER :: N4O3 = 3 ! # of ozone predictors/coefs for set4 (3)
	INTEGER, PARAMETER :: N4H2O = 13 ! # of water predictors/coefs for set4 (13)
	INTEGER, PARAMETER :: N4COEF = N4CON + N4FIX + N4CO + N4O3 + N4H2O ! total # of coefs for set4


	!
	! For set5 = sun BFSW
	!
	! Used in part by modules: 2b, 1b
	INTEGER, PARAMETER :: MXCHN5 = 210 ! max # of channels for set5 = BFSW (210)
	INTEGER, PARAMETER :: N5CON = 7 ! # of water con predictors/coefs for set5 (7)
	INTEGER, PARAMETER :: N5FIX = 11 ! # of "fixed" predictors/coefs for set5 (11)
	INTEGER, PARAMETER :: N5H2O = 3 ! # of water predictors/coefs for set5 (3)
	INTEGER, PARAMETER :: N5O3 = 1 ! # of ozone predictors/coefs for set5 (1)
	INTEGER, PARAMETER :: N5COEF = N5CON + N5FIX + N5H2O + N5O3 ! total # of coefs for set5
!
!
!	-----------------------
!	For set6 = sun MFMW
!	-----------------------
!	Used in part by modules: 1b, 2a
	INTEGER, PARAMETER :: MXCHN6 = 217 ! max # of channels for set6 = MFMW (217)
	INTEGER, PARAMETER :: N6CON = 7 ! # of water con predictors/coefs for set6 (7)
	INTEGER, PARAMETER :: N6FIX = 8 ! # of "fixed" predictors/coefs for set6 (8)
	INTEGER, PARAMETER :: N6H2O = 7 ! # of water predictors/coefs for set6 (7)
	INTEGER, PARAMETER :: N6O3 = 1 ! # of ozone predictors/coefs for set6 (1)
	INTEGER, PARAMETER :: N6COEF = N6CON + N6FIX + N6H2O + N6O3 ! total # of coefs for set6
!
!
!	-----------------------
!	For set7 = sun MFBW
!	-----------------------
!	Used in part by modules: 2a, 1a
	INTEGER, PARAMETER :: MXCHN7 = 140 ! max # of channels for set7 = MFBW (140)
	INTEGER, PARAMETER :: N7CON = 7 ! # of water con predictors/coefs for set7 (7)
	INTEGER, PARAMETER :: N7FIX = 8 ! # of "fixed" predictors/coefs for set7 (8)
	INTEGER, PARAMETER :: N7H2O = 13 ! # of water predictors/coefs for set7 (13)
	INTEGER, PARAMETER :: N7O3 = 1 ! # of ozone predictors/coefs for set7 (1)
	INTEGER, PARAMETER :: N7COEF = N7CON + N7FIX + N7H2O + N7O3 ! total # of coefs for set7
!
!
!	---------------
!	For trace gases predictors
!	---------------
	INTEGER, PARAMETER :: NTRACE = 7 ! number of trace gas perturbation predictors (7)	
!
!
!	----------------
!	For variable CO2
!	----------------
!	Used in part by modules: 12, 11, 10, 9, 7, 6, 5, 2b, 1b, 2a
	INTEGER, PARAMETER :: MXCHNC = 1082 ! max # of channels with CO2 pert coefs (1082)
	INTEGER, PARAMETER :: NCO2 = 5 ! number of CO2 coefficients
!
!
!	----------------
!	For variable SO2
!	----------------
	INTEGER, PARAMETER :: MXCHNS = 602 ! max # of channels with SO2 pert coefs (602)
	INTEGER, PARAMETER ::  NSO2 = 4 ! number of SO2 coefficients
!
!
!	-----------------
!	For variable HNO3
!	-----------------
	INTEGER, PARAMETER :: MXCHNH = 383 ! max # of channels with HNO3 pert coefs (383)
	INTEGER, PARAMETER :: NHNO3 = 4 ! number of HNO3 coefficients
!
!
!	-----------------
!	For variable N2O
!	-----------------
	INTEGER, PARAMETER :: MXCHNN = 586 ! max # of channels with N2O pert coefs (586)
	INTEGER, PARAMETER :: NN2O = 7 ! number of N2O coefficients
!
!
!	----------------------
!	For OPTRAN water coefs
!	----------------------
!	Used in part by modules:
	INTEGER, PARAMETER :: MXCHNW = 754 ! max # of channelss with OPTRAN H2O coefs (754)
	INTEGER, PARAMETER :: MXOWLY = 300 ! number of OPTRAN water layers
	INTEGER, PARAMETER :: NOWAVG = 4 ! # of OPTRAN water average profile values (4)
	INTEGER, PARAMETER :: NH2O = 9 ! number of OPTRAN H2O predictors/coefs (9)
!
!
!	-----------
!	For non-LTE
!	-----------
	INTEGER, PARAMETER :: MXCNTE = 203   ! max # of channels for non-LTE (203)
	INTEGER, PARAMETER :: NNCOEF = 7     ! # of coefs for non-LTE
	INTEGER, PARAMETER :: NTEBOT = 10    ! bottom layer for CO2TOP calc
	REAL, PARAMETER :: CO2NTE = 370.0 ! ref CO2 mixing ratio for non-LTE coefs (ppmv)
!
!
!      ---------
!      Filenames
!      ---------
!
	CHARACTER (LEN = *), PARAMETER :: FNCOF1='/asl/data/sarta_database/Data_AIRS_apr08/Coef/set1_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNCOF2='/asl/data/sarta_database/Data_AIRS_apr08/Coef/set2_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNCOF3='/asl/data/sarta_database/Data_AIRS_apr08/Coef/set3_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNCOF4='/asl/data/sarta_database/Data_AIRS_apr08/Coef/set4_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNCOF5='/asl/data/sarta_database/Data_AIRS_apr08/Coef/set5_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNCOF6='/asl/data/sarta_database/Data_AIRS_apr08/Coef/set6_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNCOF7='/asl/data/sarta_database/Data_AIRS_apr08/Coef/set7_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNCO2 ='/asl/data/sarta_database/Data_AIRS_apr08/Coef/CO2_5term_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNSO2 ='/asl/data/sarta_database/Data_AIRS_apr08/Coef/SO2_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNHNO3 ='/asl/data/sarta_database/Data_AIRS_apr08/Coef/HNO3_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNN2O ='/asl/data/sarta_database/Data_AIRS_apr08/Coef/N2O_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNOPTR='/asl/data/sarta_database/Data_AIRS_apr08/Coef/optran_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNTHER='/asl/data/sarta_database/Data_AIRS_apr08/Coef/therm_m140x370.dat'
	CHARACTER (LEN = *), PARAMETER :: FNCOFN = '/asl/data/sarta_database/Data_AIRS_apr08/Coef/nonLTE7_m140x.dat'
	CHARACTER (LEN = *), PARAMETER :: FNFX  ='/asl/data/sarta_database/Data_AIRS_apr08/Coef/fx.txt'
	CHARACTER (LEN = *), PARAMETER :: FNPREF='/asl/data/sarta_database/Data_AIRS_apr08/Coef/profref_trace370tuned'
	CHARACTER (LEN = *), PARAMETER :: FNSUN ='/asl/data/sarta_database/Data_AIRS_apr08/Solar/solar_m140x.txt'
!
!
! Mie lookup tables; also see "fnmie.f"
!
	INTEGER, PARAMETER :: MXMIEA = 10 ! max # of mie particle sizes ! ice aggregates=8, all others 10
	INTEGER, PARAMETER :: NMIETY = 3 ! number of mie particle types
!
!
! Tuning filename
	CHARACTER (LEN = *), PARAMETER :: FNTMLT='/asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_m140x.txt'
!	$ // 'tunmlt_ones.txt')
!	$ // 'tunmlt_wcononly.txt')

!
!
!	----------------
!	I/O unit numbers
!	----------------
!	Note: these units are not explicitly openned by the sarta code,
!	they should be set to standard I/O units for your compiler
	INTEGER, PARAMETER :: IOINFO = 6 ! unit number for non-error info messages (6)
	INTEGER, PARAMETER :: IOERR = 0 ! unit number for error messages (2 or 6)

	INTEGER, PARAMETER :: IOUNIT = 42 	! Unit number for self-contained readers 
										! Should not be passed to subroutines by calling routine
										! FWI 12/20/17

!
!	-----------------
!	Allowed input GUC (Gas Units Code number)
!	-----------------
	INTEGER GUCIN  ! The one & only allowed input GUC number
!	Note: GUCIN must be 1 or 2.  All gases in the input RTP
!	must be of this type.
	PARAMETER( GUCIN = 1 ) ! GUC number for:  molecules/cm^2
!	PARAMETER( GUCIN = 2 ) ! GUC number for:  kilomoles/cm^2

	!
	!	CLOUD structure (formerly in RTPDEFs)
	!	Added by Bill 
	!
	TYPE CLOUD_TYPE_PCLSAM
		! clear flag/code
		!INTEGER :: CLRFLAG		! clear flag/code  (was originally set as LOGICAL)
	
		! cloud1 data
		INTEGER :: CTYPE1 ! cloud type code
		
		REAL :: CFRAC1 ! cloud fraction 
		REAL :: CEMIS1(MXCHAN) ! cloud top emissivity
		REAL :: CRHO1(MXCHAN) ! cloud top reflectivity
		REAL :: CPRTOP1 ! cloud top pressure
		REAL :: CPRBOT1 ! cloud bottom pressure
		REAL :: CNGWAT1 ! cloud non-gas water
		REAL :: CPSIZE1 ! cloud particle size
		REAL :: CSTEMP1 ! cloud surface temperature		
		INTEGER :: LBLAC1  ! black cloud 1? Mie cloud if false (was originally set as LOGICAL)

		! cloud2 data
		INTEGER :: CTYPE2 ! cloud2 type code
		REAL :: CFRAC2 ! cloud2 fraction 
		REAL :: CEMIS2(MXCHAN) ! cloud2 top emissivity
		REAL :: CRHO2(MXCHAN) ! cloud2 top reflectivity
		REAL :: CPRTOP2 ! cloud2 top pressure
		REAL :: CPRBOT2 ! cloud2 bottom pressure
		REAL :: CNGWAT2 ! cloud2 non-gas water
		REAL :: CPSIZE2 ! cloud2 particle size
		REAL :: CSTEMP2 ! cloud2 surface temperature		
		INTEGER :: LBLAC2  ! black cloud 2? Mie cloud if false (was originally set as LOGICAL)

		REAL :: CFRAC12 ! cloud1+2 fraction

	END TYPE CLOUD_TYPE_PCLSAM


END MODULE INCFTC
