C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C    CrIS FSR using airs-oco2_pbl layering
C    incFTC

!F77====================================================================


!ROUTINE NAME:
C    incFTC (include file)


!ABSTRACT:
C    Include file consisting of parameter statements to size various
C    arrays in the SARTA related routines source code.


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
C    April 2009 version of the 100 layer CrIS fast model
C    code by L.L.Strow/S.Hannon.  This CrIS model
C    uses the same algorithm and source code (except for this
C    include file) as our AIRS fast model.
C
C    Parameter statements for the FTC routines.


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C Date        Programmer     Comments
C ----------- -------------- -------------------------------------------
C 27 Apr 2009 Scott Hannon   Created for CrIS April 2009 fast model
C 12 May 2009 Scott Hannon   Add VTUNNG string; delete VCLOUD
C 14 Sep 2018 C Hepplewhite  Updated IASI fast model
C 1  Mar 2019 C Hepplewhite  Added HDO
C 1  Jul 2019 C Hepplewhite  Added minor gas logicals
C    Mar 2021 C Hepplewhite  Added nh3, and other minor gases.
C    Jul 2022 C Hepplewhite  prod_2022, build: jul2022 HITRAN: 2020.

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
       LOGICAL DEBUG
       PARAMETER(DEBUG = .FALSE.)
C
       LOGICAL CFCO2
       LOGICAL CFHNO3
       LOGICAL CFN2O
       LOGICAL CFNH3
       LOGICAL CFSO2
       LOGICAL CFHDO
       LOGICAL CFTHER
       LOGICAL CFOPTR
       LOGICAL COFNTE
       PARAMETER(CFCO2  = .TRUE.)
       PARAMETER(CFHNO3 = .FALSE.)
       PARAMETER(CFN2O  = .FALSE.)
       PARAMETER(CFNH3  = .FALSE.)
       PARAMETER(CFSO2  = .FALSE.)
       PARAMETER(CFHDO  = .FALSE.)
c       PARAMETER(CFOPTR = .FALSE.)
       PARAMETER(CFOPTR = .TRUE.)
       PARAMETER(CFTHER = .TRUE.)
       PARAMETER(COFNTE = .FALSE.)

       ! PARAMETER(CFHNO3 = .TRUE.)
       ! PARAMETER(CFN2O  = .TRUE.)
       ! PARAMETER(CFNH3  = .TRUE.)
       ! PARAMETER(CFSO2  = .TRUE.)
       ! PARAMETER(CFHDO  = .FALSE.)
       ! PARAMETER(CFOPTR = .TRUE.)
       ! PARAMETER(CFTHER = .TRUE.)
       ! PARAMETER(COFNTE = .TRUE.)
C
       CHARACTER*40 VSARTA  ! SARTA source code version
       CHARACTER*40 VSCOEF  ! SARTA coefficient version
       CHARACTER*40 VTUNNG  ! optical depth tuning version
C      version template    '#.## YYYY-MM-DD <--------comment------->'
       PARAMETER( VSARTA = '2.02 2021-03-23. PBL prod_2025' )
       PARAMETER( VSCOEF = 'CrIS g4 0.8/0.8/0.8cm Hamming Jul 2022')
       PARAMETER( VTUNNG = 'none' )

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

c this is for Chou scaling
       REAL rXTang ! do we include Tang correction???? the paper says 0.5, but I found use 0.3 for ice, 0.1 for water clouds
                   ! use -1 for similarity (default before July 2023), -2 for chou, -3 for Maestri/Martinazzo
                   ! ccprep_slab.f used as ISCALING = ABS(FLOOR(rXTang)), then calrad1,calrad2 uses this as well
       PARAMETER(rXTang = 0.0000000)

Ccc    Previously used values; agrees w/JPL pre-Dec2000
Ccc    PARAMETER(    C1 = 1.1910439E-8)  ! JPL value is 1E+3 bigger
Ccc    PARAMETER(    C2 = 1.4387687)
C
C      Current values (CODATA98 from NIST); agrees w/JPL Dec2000
       PARAMETER(  C1 = 1.191042722E-8)  ! JPL value is 1E+3 bigger
       PARAMETER(  C2 = 1.4387752)
C
       REAL CO2STD ! standard CO2 PPMV mixing ratio (385)
       PARAMETER( CO2STD = 400.0 )
C
       REAL HDOSTD ! standard HDO depletion abundance (3.1069E-5)
       PARAMETER( HDOSTD = 0.00031069 )
C
       REAL HDOFCT ! vary proportion of HDO in H2O from std depletion
C                  ! (-1: 100% enhancement, 0:std HDO or zero depletion,
C                  1: 100% depleted))
       PARAMETER( HDOFCT = 0.00 )
C
       REAL  XSALT ! expected nominal satellite altitude (km)
c       PARAMETER( XSALT = 841.0 ) XXXXXX this was in earlier files
       PARAMETER( XSALT = 825.0 )
C
C      -----------------------------------
C      Channels and layers other variables
C      -----------------------------------
       INTEGER MAXLAY ! # of layers (100)
       INTEGER   NSET ! # of coefficient data sets (7)
       INTEGER MXCHAN ! max total # of channels (1305)
       INTEGER NFCOEF ! # of downwelling thermal "F" factor coefs 
       INTEGER MXEMIS ! max # of input emis/rho data points
       INTEGER MAXPRO ! max # of user specified profiles
       INTEGER  MXGAS ! max # of gases in user profile
       INTEGER MXMIEA ! max # of mie particle sizes (cloud code only)

       INTEGER MAXJAC    ! max # of jacs = 11 (since we have WV, CO2, O3, N2O, CO, CH4, SO2, NH3, HNO3, HDO, Tz .. WGT is separate)
       INTEGER OPTRANJAC ! max # of jacs = 2  (since we have WV, Tz)
       INTEGER CLDJAC    ! max # of jacs = 11  (since we have cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,cpsize1,cpsize2
       PARAMETER(MAXLAY = 100)
       PARAMETER(  NSET = 7)
       PARAMETER(MXCHAN = 2235)
       PARAMETER(NFCOEF = 6)
       PARAMETER(MXEMIS = 100)
       PARAMETER(MAXPRO = 25)
       PARAMETER( MXGAS = 44)
       PARAMETER(MXMIEA = 10)

       PARAMETER(MAXJAC = 11)
       PARAMETER(OPTRANJAC = 2)
       PARAMETER(CLDJAC = 12)  
C
C***********************************************************************
C      Variables for the coefficient sets
C***********************************************************************
C
C      --------------
C      For set1 = FWO
C      -------------
C      Used in part by modules: 12, 11, 10, 9, 8, 7, 6, 5, 3, 4b, 4a
       INTEGER MXCHN1 ! max # of channels for set1 = FWO (458)
       INTEGER  N1CON ! # of water con predictors/coefs for set1 (5)
       INTEGER  N1FIX ! # of "fixed" predictors/coefs for set1 (8)
       INTEGER  N1H2O ! # of water predictors/coefs for set1 (13)
       INTEGER   N1O3 ! # of ozone predictors/coefs for set1 (5)
       INTEGER N1COEF ! total # of coefs for set1
       PARAMETER(MXCHN1 = 493)
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
       INTEGER MXCHN2 ! max # of channels for set2 = FOW  (267)
       INTEGER  N2CON ! # of water con predictors/coefs for set2 (5)
       INTEGER  N2FIX ! # of "fixed" predictors/coefs for set2 (8)
       INTEGER   N2O3 ! # of ozone predictors/coefs for set2 (10)
       INTEGER  N2H2O ! # of water predictors/coefs for set2 (11)
       INTEGER N2COEF ! total # of coefs for set2
       PARAMETER(MXCHN2 = 228)
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
       INTEGER MXCHN3 ! max # of channels for set3 = FMW  (421)
       INTEGER  N3CON ! # of water con predictors/coefs for set3 (5)
       INTEGER  N3FIX ! # of "fixed" predictors/coefs for set3 (8)
       INTEGER  N3CH4 ! # of methane predictors/coefs for set3 (9)
       INTEGER  N3H2O ! # of water predictors/coefs for set3 (13)
       INTEGER N3COEF ! total # of coefs for set3
       PARAMETER(MXCHN3 = 873)
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
       INTEGER MXCHN4 ! max # of channels for set4 = FCOW (39)
       INTEGER  N4CON ! # of water con predictors/coefs for set4 (5)
       INTEGER  N4FIX ! # of "fixed" predictors/coefs for set4 (11)
       INTEGER   N4CO ! # of CO predictors/coefs for set4 (11)
       INTEGER   N4O3 ! # of ozone predictors/coefs for set4 (3)
       INTEGER  N4H2O ! # of water predictors/coefs for set4 (13)
       INTEGER N4COEF ! total # of coefs for set4
       PARAMETER(MXCHN4 = 126)
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
       INTEGER MXCHN5 ! max # of channels for set5 = BFSW (58)
       INTEGER  N5CON ! # of water con predictors/coefs for set5 (5)
       INTEGER  N5FIX ! # of "fixed" predictors/coefs for set5 (11)
       INTEGER  N5H2O ! # of water predictors/coefs for set5 (3)
       INTEGER   N5O3 ! # of ozone predictors/coefs for set5 (1)
       INTEGER N5COEF ! total # of coefs for set5
       PARAMETER(MXCHN5 = 75)
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
       INTEGER MXCHN6 ! max # of channels for set6 = MFMW (60)
       INTEGER  N6CON ! # of water con predictors/coefs for set6 (5)
       INTEGER  N6FIX ! # of "fixed" predictors/coefs for set6 (8)
       INTEGER  N6H2O ! # of water predictors/coefs for set6 (7)
       INTEGER   N6O3 ! # of ozone predictors/coefs for set6 (1)
       INTEGER N6COEF ! total # of coefs for set6
       PARAMETER(MXCHN6 = 415)
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
       INTEGER MXCHN7 ! max # of channels for set7 = MFBW (2)
       INTEGER  N7CON ! # of water con predictors/coefs for set7 (5)
       INTEGER  N7FIX ! # of "fixed" predictors/coefs for set7 (8)
       INTEGER  N7H2O ! # of water predictors/coefs for set7 (13)
       INTEGER   N7O3 ! # of ozone predictors/coefs for set7 (1)
       INTEGER N7COEF ! total # of coefs for set7
       PARAMETER(MXCHN7 = 25)
       PARAMETER( N7CON = 7)
       PARAMETER( N7FIX = 8)
       PARAMETER( N7H2O = 13)
       PARAMETER(  N7O3 = 1)
       PARAMETER(N7COEF = N7CON + N7FIX + N7H2O + N7O3 )
C
C
C      ---------------
C      For trace gases predictors
C      ---------------
       INTEGER NTRACE ! number of trace gas perturbation predictors (7)
       PARAMETER(NTRACE = 7)
C
C
C      ----------------
C      For variable CO2
C      ----------------
C      Used in part by modules: 12, 11, 10, 9, 7, 6, 5, 2b, 1b, 2a
       INTEGER MXCHNC ! max # of channels with CO2 pert coefs (was 567)
       INTEGER NCO2   ! number of CO2 pert predictors/coefs (5)
       PARAMETER(MXCHNC = 893)    ! placeholder
       PARAMETER(  NCO2 = 5)       ! was 4-term. otherwise 5-term XXXX before
c       PARAMETER(  NCO2 = 4)

C
C
C      ----------------
C      For variable SO2
C      ----------------
       INTEGER MXCHNS ! max # of channels with SO2 pert coefs (was 212)
       INTEGER   NSO2 ! number of SO2 coefficients
C       PARAMETER(MXCHNS = 295)   ! was 396?   XXXX 
       PARAMETER(MXCHNS = 396)    ! placeholder
       PARAMETER(  NSO2 = 4)
C
C
C      -----------------
C      For variable HNO3
C      -----------------
       INTEGER MXCHNH ! max # of channels with HNO3 pert coefs (was 253)
       INTEGER  NHNO3 ! number of HNO3 coefficients
C       PARAMETER(MXCHNH = 570)   ! was 570   XXXX
       PARAMETER(MXCHNH = 374)    ! placeholder
       PARAMETER( NHNO3 = 4)
C
C
C      -----------------
C      For variable N2O
C      -----------------
       INTEGER MXCHNN ! max # of channels with N2O pert coefs (was 181)
       INTEGER   NN2O ! number of N2O coefficients
C       PARAMETER(MXCHNN = 433)   ! was 289  XXXX
       PARAMETER(MXCHNN = 544)
       PARAMETER(  NN2O = 7)
C
C      -----------------
C      For variable NH3
C      -----------------
       INTEGER MXCHNA ! max # of channels with NH3 pert coefs (2075)
       INTEGER   NNH3           ! number of NH3 coefficients (4)
c       PARAMETER(MXCHNA = 1)        ! placeholder when not using this set XXXX was commented before
C       PARAMETER( NNH3 = 1)         ! placeholder when not using this set
C       PARAMETER(MXCHNA = 1422)   XXXXXX
       PARAMETER(MXCHNA = 605)
       PARAMETER(  NNH3 = 4)
C
C      -----------------
C      For variable HDO
C      -----------------
       INTEGER MXCHND ! max # of channels with HDO pert coefs (394)
       INTEGER   NHDO ! number of HDO coefficients (4)
       PARAMETER(MXCHND = 1)        ! placeholder when not using this set XXXX before 
C       PARAMETER(MXCHND = 1)        ! placeholder when not using this set
C       PARAMETER( NHDO = 1)         ! placeholder when not using this set
C       PARAMETER(MXCHND = 165)
       PARAMETER(  NHDO = 11)

C
C      ----------------------
C      For OPTRAN water coefs
C      ----------------------
C      Used in part by modules:
       INTEGER MXCHNW ! max # of channelss with OPTRAN H2O coefs (873)
       INTEGER MXOWLY ! number of OPTRAN water layers
       INTEGER NOWAVG ! # of OPTRAN water average profile values (4)
       INTEGER NH2O   ! number of OPTRAN H2O predictors/coefs (9)
C       PARAMETER(MXCHNW = 873)  XXXXX
       PARAMETER(MXCHNW = 708)
       PARAMETER(MXOWLY = 300)
       PARAMETER(NOWAVG = 4)
       PARAMETER(  NH2O = 9)
C
C      -----------
C      For non-LTE
C      -----------
       INTEGER MXCNTE ! max # of channels for non-LTE (264)
       INTEGER NNCOEF ! # of coefs for non-LTE (7)
       INTEGER NTEBOT ! bottom layer for CO2TOP calc
       REAL CO2NTE ! ref CO2 mixing ratio for non-LTE coefs (ppmv)
       PARAMETER(MXCNTE = 264)        ! was 264 placeholder
       PARAMETER(NNCOEF = 7)          ! was 7
       PARAMETER(NTEBOT = 10)
       PARAMETER(CO2NTE = 400.0)
C
C      ---------
C      Filenames
C      ---------
       CHARACTER*90 FNCOF1 ! coef set1 
       CHARACTER*90 FNCOF2 ! coef set2 
       CHARACTER*90 FNCOF3 ! coef set3 
       CHARACTER*90 FNCOF4 ! coef set4 
       CHARACTER*90 FNCOF5 ! coef set5 
       CHARACTER*90 FNCOF6 ! coef set6 
       CHARACTER*90 FNCOF7 ! coef set7 
       CHARACTER*90 FNCO2  ! coef CO2
       CHARACTER*90 FNSO2  ! coef SO2
       CHARACTER*90 FNHNO3 ! coef HNO3
       CHARACTER*90 FNN2O  ! coef N2O
       CHARACTER*90 FNNH3  ! coef NH3
       CHARACTER*90 FNHDO  ! coef HDO
       CHARACTER*90 FNOPTR ! coef optran
       CHARACTER*90 FNTHER ! coef therm
       CHARACTER*90 FNFX   ! coef fx
       CHARACTER*90 FNPREF ! ref prof
       CHARACTER*90 FNSUN  ! solar data
       CHARACTER*90 FNCOFN ! non-LTE
C
C
      PARAMETER(FNCOF1=
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/set1.dat')
      PARAMETER(FNCOF2=
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/set2.dat')
      PARAMETER(FNCOF3=
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/set3.dat')
      PARAMETER(FNCOF4=
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/set4.dat')
      PARAMETER(FNCOF5=
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/set5.dat')
      PARAMETER(FNCOF6=
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/set6.dat')
      PARAMETER(FNCOF7=
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/set7.dat')
      PARAMETER(FNOPTR=
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/optran.dat')
               
      PARAMETER(FNCO2 =
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/co2.dat')
      PARAMETER(FNSO2 =
     $     '/home/chepplew/data/sarta/prod_2021/cris_hr/dec2018/dbase/Coef/so2.dat')
      PARAMETER(FNHNO3 =
     $     '/home/chepplew/data/sarta/prod_2021/cris_hr/dec2018/dbase/Coef/hno3.dat')
      PARAMETER(FNN2O =
     $     '/home/chepplew/data/sarta/prod_2021/cris_hr/dec2018/dbase/Coef/n2o.dat')
      PARAMETER(FNNH3 =
     $     '/home/chepplew/data/sarta/prod_2021/cris_hr/dec2018/dbase/Coef/nh3.dat')
      PARAMETER(FNHDO =  
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/hdo.dat')
                
      PARAMETER(FNFX  =
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/fx_pbl.txt')
      PARAMETER(FNPREF =
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/refprof_400ppm_pbl')
      PARAMETER(FNSUN =
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Solar/sol.txt')
             
      PARAMETER(FNTHER =
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/refl_therm.dat')
      PARAMETER(FNCOFN =
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/xnte_2x7term.dat')
            
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Tuning filename
      CHARACTER (LEN=90) :: FNTMLT ! tuning multiplier filename
                
      PARAMETER(FNTMLT=
     $     '/home/chepplew/data/sarta/prod_2025/cris_hr_pbl/jan2025a/dbase/Coef/' 
     $   // 'tunmlt_ones.txt')
                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Mie lookup tables; also see "fnmie.f"
C
c       INTEGER MXMIEA  ! max # of mie particle sizes
c       PARAMETER(MXMIEA = 10) ! ice aggregates=8, all others 10
       INTEGER NMIETY  ! number of mie particle types
       PARAMETER(NMIETY = 3)
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
C      Note: GUCIN must be 1 or 2.  All gases in the input RTP
C      must be of this type.
       PARAMETER( GUCIN = 1 ) ! GUC number for:  molecules/cm^2
c       PARAMETER( GUCIN = 2 ) ! GUC number for:  kilomoles/cm^2

c
c rtpV201 compatibility
c       CHARACTER*40 VCLOUD
c       PARAMETER( VCLOUD = 'no clouds' )

C      End of include file
