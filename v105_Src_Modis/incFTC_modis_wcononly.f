C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS (Atmospheric Infra-Red Sounder)
C
C    incFTC
C
!F77====================================================================


!ROUTINE NAME:
C    incFTC (include file)


!ABSTRACT:
C    Include file consisting of parameter statements to size various
C    arrays in the USEFAST related routines source code.


!CALL PROTOCOL:
C    none (include file)


!INPUT PARAMETERS:
C    none


!OUTPUT PARAMETERS:
C    none


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    CALOKW
C    CALOWP
C    CALPAR
C    CALRAD
C    CALT1
C    CALT2
C    CALT3
C    CALT4
C    CALT5
C    CALT6
C    CALT7
C    FAKETZ
C    RDCOEF
C    RDLIST
C    RDPROF
C    RDSUN
C    SUNPAR
C    SARTA


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    none


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    Include file for the January 2004 100 layer AIRS fast
C    Stand Alone RTA (SARTA) code by L.L.Strow/S.Hannon.
C
C    Parameter statements for the FTC routines.


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C Date        Programmer     Comments
C ----------- -------------- -------------------------------------------
C  1 Dec 1994 Scott Hannon   Created
C 31 Jan 1997 Scott Hannon   Re-wrote for FWO+FOW+FMW+FCOW=Feb97 FTC
C  3 Sep 1997 Scott Hannon   Re-wrote for sets 1 - 7
C 30 Sep 1997 Scott Hannon   Added NCO2 and MXCHNC
C 26 Feb 1998 Scott Hannon   Added OPTRAN variables for water, and
C                            changed both N1H2O & N3H2O from 13 to 11
C 23 Sep 1999 Scott Hannon   Change number of channel dimensions for
C                            new Sep99 version of fast model.
C  5 Apr 2000 Scott Hannon   Added MXEMIS
C  4 Aug 2000 Scott Hannon   Changes values for use with testfast
C 11 Aug 2000 Scott Hannon   Change from 4 to 5 term H2O continuum
C 23 Jan 2001 Scott Hannon   Update values of C1 & C2
C 15 Feb 2001 Scott Hannon   Add MAXPRO, CO2STD, IOERR, IOINFO,
C                            MXGAS, CSARTA, and all filenames
C 24 Apr 2001 Scott Hannon   Add MXMIEA and FNMIEA
C 14 Aug 2001 Scott Hannon   Add FNMIEE and FNMIEG
C 21 Nov 2001 Scott Hannon   Add VSARTA, VSCOEF, & VCLOUD; remove CSARTA
C 12 Sep 2002 Scott Hannon   Updated for m135f (-13.5 um with fringes)
C 17 Dec 2002 Scott Hannon   Updated for revised(Dec02) m135f 
C  3 Jan 2003 Scott Hannon   Updated VSARTA for version 1.04
C 06 Feb 2004 Scott Hannon   Add FNTMLT & update VSARTA for v1.05
C 15 Aug 2005 Scott Hannon   Add FNNU,NCOFNU, and COEFNU for MODIS channel
C                               freq coefs


!END====================================================================
C
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
C Note: having an "implicit none" in both the include file & the main
C source code will cause some compilers to complain.
c       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      EXECUTABLE CODE
C-----------------------------------------------------------------------
C      none

C      -----------------------------------------------------------------
C      Assign SARTA version strings
C      -----------------------------------------------------------------
C      The version strings consists of 3 parts: version number, date,
C      and comment.  The version date should be updated to the
C      current date whenever any portion of the code is updated.  The
C      version number consists of two parts; a major version to the
C      left of the decimal point, and a minor version to the right.
C      The major number should be incremented only when major changes
C      have been made to the overall SARTA code.  The minor number
C      should be incremented only when minor but non-trivial changes
C      are made to the code.  Bug fixes should generally be handled
C      with the version date, but a fix for a serious bug may warrant
C      a change to the minor version number.
C      See the "Doc/last_update.txt" file for a description of the
C      changes associated with every change of VSARTA.
C
       CHARACTER*40 VSARTA  ! SARTA source code version
       CHARACTER*40 VSCOEF  ! SARTA coefficient version
       CHARACTER*40 VCLOUD  ! cloud version
C      version template    '#.## YY-MM-DD <--------comment--------->'
       PARAMETER( VSARTA = '1.05 MODIS 05-08-26 w/chan freq shift' )
       PARAMETER( VSCOEF = '05-08-26 SRF=mean Aqua tun=AIRSjan04')
       PARAMETER( VCLOUD = 'no clouds' )

C      *********
C      VARIABLES
C      *********
C      Note: these should not be changed by the user
C
C      ------------------------
C      Constants and other data
C      ------------------------
       REAL     PI ! pi, circle circumference/diameter (3.1415926)
       REAL RADSUN ! radius of the sun (6.956E+8 m)
       REAL     C1 ! radiation constant c1 (1.1911E-8  W/(m2.st.(cm-1)4)
       REAL     C2 ! radiation constant c2 (1.4387863 K/cm-1)
       PARAMETER(    PI = 3.1415926)
       PARAMETER(RADSUN = 6.956E+8)
C
Ccc    Previously used values; agrees w/JPL pre-Dec2000
Ccc    PARAMETER(    C1 = 1.1910439E-8)  ! JPL value is 1E+3 bigger
Ccc    PARAMETER(    C2 = 1.4387687)
C
C      Current values (CODATA98 from NIST); agrees w/JPL Dec2000
       PARAMETER(  C1 = 1.191042722E-8)  ! JPL value is 1E+3 bigger
       PARAMETER(  C2 = 1.4387752)
C
       REAL CO2STD ! standard CO2 PPMV mixing ratio (370)
       PARAMETER( CO2STD = 370.0 )
C
C
C      -----------------------------------
C      Channels and layers other variables
C      -----------------------------------
       INTEGER MAXLAY ! # of layers (100)
       INTEGER   NSET ! # of coefficient data sets (7)
       INTEGER MXCHAN ! max total # of channels (2378)
       INTEGER NFCOEF ! # of downwelling thermal "F" factor coefs 
       INTEGER MXEMIS ! max # of input emis/rho data points
       INTEGER MAXPRO ! max # of user specified profiles
       INTEGER  MXGAS ! max # of gases in user profile
       INTEGER MXMIEA ! max # of mie particle sizes
       PARAMETER(MAXLAY = 100)
       PARAMETER(  NSET = 7)
       PARAMETER(MXCHAN = 36)
       PARAMETER(NFCOEF = 5)
       PARAMETER(MXEMIS = 100)
       PARAMETER(MAXPRO = 25)
       PARAMETER( MXGAS = 44)
       PARAMETER(MXMIEA = 10)
C
C***********************************************************************
C      Variables for the coefficient sets
C***********************************************************************
C
C      --------------
C      For set1 = FWO
C      -------------
C      Used in part by modules: 12, 11, 10, 9, 8, 7, 6, 5, 3, 4b, 4a
       INTEGER MXCHN1 ! max # of channels for set1 = FWO (1090)
       INTEGER  N1CON ! # of water con predictors/coefs for set1 (5)
       INTEGER  N1FIX ! # of "fixed" predictors/coefs for set1 (8)
       INTEGER  N1H2O ! # of water predictors/coefs for set1 (13)
       INTEGER   N1O3 ! # of ozone predictors/coefs for set1 (5)
       INTEGER N1COEF ! total # of coefs for set1
       PARAMETER(MXCHN1 = 8)
       PARAMETER( N1CON = 7)
       PARAMETER( N1FIX = 8)
       PARAMETER( N1H2O = 11)
       PARAMETER(  N1O3 = 5)
       PARAMETER(N1COEF = N1CON + N1FIX + N1H2O + N1O3 )
C
C
C      --------------
C      For set2 = FOW
C      --------------
C      Used in part by modules: 6, 5
       INTEGER MXCHN2 ! max # of channels for set2 = FOW  (382)
       INTEGER  N2CON ! # of water con predictors/coefs for set2 (5)
       INTEGER  N2FIX ! # of "fixed" predictors/coefs for set2 (8)
       INTEGER   N2O3 ! # of ozone predictors/coefs for set2 (10)
       INTEGER  N2H2O ! # of water predictors/coefs for set2 (11)
       INTEGER N2COEF ! total # of coefs for set2
       PARAMETER(MXCHN2 = 1)
       PARAMETER( N2CON = 7)
       PARAMETER( N2FIX = 8)
       PARAMETER(  N2O3 = 10)
       PARAMETER( N2H2O = 11)
       PARAMETER(N2COEF = N2CON + N2FIX + N2O3 + N2H2O )
C
C
C      --------------
C      For set3 = FMW
C      --------------
C      Used in part by modules: 4d, 4c, 3
       INTEGER MXCHN3 ! max # of channels for set3 = FMW  (392)
       INTEGER  N3CON ! # of water con predictors/coefs for set3 (5)
       INTEGER  N3FIX ! # of "fixed" predictors/coefs for set3 (8)
       INTEGER  N3CH4 ! # of methane predictors/coefs for set3 (9)
       INTEGER  N3H2O ! # of water predictors/coefs for set3 (13)
       INTEGER N3COEF ! total # of coefs for set3
       PARAMETER(MXCHN3 = 1)
       PARAMETER( N3CON = 7)
       PARAMETER( N3FIX = 8)
       PARAMETER( N3CH4 = 9)
       PARAMETER( N3H2O = 11)
       PARAMETER(N3COEF = N3CON + N3FIX + N3CH4 + N3H2O )
C
C
C      ---------------
C      For set4 = sun FCOW
C      ---------------
C      Used in part by modules: 2b
       INTEGER MXCHN4 ! max # of channels for set4 = FCOW (55)
       INTEGER  N4CON ! # of water con predictors/coefs for set4 (5)
       INTEGER  N4FIX ! # of "fixed" predictors/coefs for set4 (11)
       INTEGER   N4CO ! # of CO predictors/coefs for set4 (11)
       INTEGER   N4O3 ! # of ozone predictors/coefs for set4 (3)
       INTEGER  N4H2O ! # of water predictors/coefs for set4 (13)
       INTEGER N4COEF ! total # of coefs for set4
       PARAMETER(MXCHN4 = 1)
       PARAMETER( N4CON = 7)
       PARAMETER( N4FIX = 11)
       PARAMETER(  N4CO = 11)
       PARAMETER(  N4O3 = 3)
       PARAMETER( N4H2O = 13)
       PARAMETER(N4COEF = N4CON + N4FIX + N4CO + N4O3 + N4H2O )
C
C
C      -----------------------
C      For set5 = sun BFSW
C      -----------------------
C      Used in part by modules: 2b, 1b
       INTEGER MXCHN5 ! max # of channels for set5 = BFSW (202)
       INTEGER  N5CON ! # of water con predictors/coefs for set5 (5)
       INTEGER  N5FIX ! # of "fixed" predictors/coefs for set5 (11)
       INTEGER  N5H2O ! # of water predictors/coefs for set5 (3)
       INTEGER   N5O3 ! # of ozone predictors/coefs for set5 (1)
       INTEGER N5COEF ! total # of coefs for set5
       PARAMETER(MXCHN5 = 1)
       PARAMETER( N5CON = 7)
       PARAMETER( N5FIX = 11)
       PARAMETER( N5H2O = 3)
       PARAMETER(  N5O3 = 1)
       PARAMETER(N5COEF = N5CON + N5FIX + N5H2O + N5O3 )
C
C
C      -----------------------
C      For set6 = sun MFMW
C      -----------------------
C      Used in part by modules: 1b, 2a
       INTEGER MXCHN6 ! max # of channels for set6 = MFMW (174)
       INTEGER  N6CON ! # of water con predictors/coefs for set6 (5)
       INTEGER  N6FIX ! # of "fixed" predictors/coefs for set6 (8)
       INTEGER  N6H2O ! # of water predictors/coefs for set6 (7)
       INTEGER   N6O3 ! # of ozone predictors/coefs for set6 (1)
       INTEGER N6COEF ! total # of coefs for set6
       PARAMETER(MXCHN6 = 3)
       PARAMETER( N6CON = 7 )
       PARAMETER( N6FIX = 8 )
       PARAMETER( N6H2O = 7 )
       PARAMETER(  N6O3 = 1 )
       PARAMETER(N6COEF = N6CON + N6FIX + N6H2O + N6O3 )
C
C
C      -----------------------
C      For set7 = sun MFBW
C      -----------------------
C      Used in part by modules: 2a, 1a
       INTEGER MXCHN7 ! max # of channels for set7 = MFBW (83)
       INTEGER  N7CON ! # of water con predictors/coefs for set7 (5)
       INTEGER  N7FIX ! # of "fixed" predictors/coefs for set7 (8)
       INTEGER  N7H2O ! # of water predictors/coefs for set7 (13)
       INTEGER   N7O3 ! # of ozone predictors/coefs for set7 (1)
       INTEGER N7COEF ! total # of coefs for set7
       PARAMETER(MXCHN7 = 1)
       PARAMETER( N7CON = 7)
       PARAMETER( N7FIX = 8)
       PARAMETER( N7H2O = 13)
       PARAMETER(  N7O3 = 1)
       PARAMETER(N7COEF = N7CON + N7FIX + N7H2O + N7O3 )
C
C
C      ----------------
C      For variable CO2
C      ----------------
C      Used in part by modules: 12, 11, 10, 9, 7, 6, 5, 2b, 1b, 2a
       INTEGER MXCHNC ! max # of channels with CO2 pert coefs (868)
       INTEGER NCO2   ! number of CO2 pert predictors/coefs (4)
       PARAMETER(MXCHNC = 6)
       PARAMETER(  NCO2 = 4)
C
C
C      ----------------------
C      For OPTRAN water coefs
C      ----------------------
C      Used in part by modules:
       INTEGER MXCHNW ! max # of channelss with OPTRAN H2O coefs (637)
       INTEGER MXOWLY ! number of OPTRAN water layers
       INTEGER NOWAVG ! # of OPTRAN water average profile values (4)
       INTEGER NH2O   ! number of OPTRAN H2O predictors/coefs (9)
       PARAMETER(MXCHNW = 6)
       PARAMETER(MXOWLY = 300)
       PARAMETER(NOWAVG = 4)
       PARAMETER(  NH2O = 9)
C
C
C      ----------------------
C      For channel freq coefs
C      ----------------------
       INTEGER NCOFNU ! number of channel freq coefs (10)
       PARAMETER(NCOFNU = 10)

C      ---------
C      Filenames
C      ---------
       CHARACTER*80 FNCOF1 ! coef set1 
       CHARACTER*80 FNCOF2 ! coef set2 
       CHARACTER*80 FNCOF3 ! coef set3 
       CHARACTER*80 FNCOF4 ! coef set4 
       CHARACTER*80 FNCOF5 ! coef set5 
       CHARACTER*80 FNCOF6 ! coef set6 
       CHARACTER*80 FNCOF7 ! coef set7 
       CHARACTER*80 FNCO2  ! coef co2
       CHARACTER*80 FNOPTR ! coef optran
       CHARACTER*80 FNTHER ! coef therm
       CHARACTER*80 FNFX   ! coef fx
       CHARACTER*80 FNPREF ! ref prof
       CHARACTER*80 FNSUN  ! solar data
       CHARACTER*80 FNNU   ! channel freq coefs
C
       PARAMETER(FNCOF1=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/set1.dat')
       PARAMETER(FNCOF2=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/set2.dat')
       PARAMETER(FNCOF3=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/set3.dat')
       PARAMETER(FNCOF4=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/set4.dat')
       PARAMETER(FNCOF5=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/set5.dat')
       PARAMETER(FNCOF6=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/set6.dat')
       PARAMETER(FNCOF7=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/set7.dat')
       PARAMETER(FNCO2 =
     $ '/asl/data/sarta_database/Data_MODIS/Coef/co2.dat')
       PARAMETER(FNOPTR=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/optran.dat')
       PARAMETER(FNTHER=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/therm_pad.dat')
C
       PARAMETER(FNFX  =
     $ '/asl/data/sarta_database/Data_MODIS/Coef/fx.txt')
       PARAMETER(FNPREF=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/profref_truJul00')
       PARAMETER(FNSUN =
     $ '/asl/data/sarta_database/Data_MODIS/Solar/solar_data.txt')
       PARAMETER(FNNU =
     $ '/asl/data/sarta_database/Data_MODIS/Coef/nu10_pad.dat')

C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Tuning filename
       CHARACTER*80 FNTMLT ! tuning multiplier filename
C
       PARAMETER(FNTMLT=
     $ '/asl/data/sarta_database/Data_MODIS/Coef/'
     $ // 'tunmlt_wcononly_pad.txt')
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C      ----------------
C      I/O unit numbers
C      ----------------
C      Note: these units are not explicitly openned by the sarta code,
C      they should be set to standard I/O units for your compiler
       INTEGER IOINFO  ! unit number for non-error info messages (6)
       INTEGER IOERR   ! unit number for error messages (2 or 6)
       PARAMETER( IOINFO = 6 )
       PARAMETER( IOERR = 0 )
C
C
C      -----------------
C      Allowed input GUC (Gas Units Code number)
C      -----------------
       INTEGER GUCIN  ! The one & only allowed input GUC number
       PARAMETER( GUCIN = 1 ) ! GUC number for:  molecules/cm^2
c       PARAMETER( GUCIN = 2 ) ! GUC number for:  kilomoles/cm^2
C      Note: GUCIN must be 1 or 2.  All gases in the input RTP
C      must be of this type.

C      End of include file
