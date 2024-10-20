C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              incLAY
C
!F77====================================================================

!ROUTINE NAME: incLAY

!ABSTRACT:
C    Include file for KLAYERS program and routines INTEG and INTLEV.

!CALL PROTOCOL: none

!INPUT PARAMETERS: none

!OUTPUT PARAMETERS: none

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S):
C    KLAYERS
C    INTEG
C    INTLEV

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C    Declare some array dimensions

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS: none

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C Mar  3 1995 Scott Hannon/UMBC created
C Jun 23 1995 Scott Hannon      Correct some comments
C Jul 30 1996 Scott Hannon      Increased MXIN from 300 to 500
C Jul 19 2000 Scott Hannon      Add NAMREF (used by layout_reg only)
C 15 Nov 2000 Scott Hannon      PLEV moved to here from klayers.f
C 21 Nov 2000 Scott Hannon      Add MRTH2O for water mixing ratio type
C  1 Feb 2001 Scott Hannon      Add NGIDS
C 21 Feb 2001 Scott Hannon      Add MAXGUC and IOERR & IOINFO
C 24 Jul 2001 Scott Hannon      Add CO2STD
C 30 Aug 2001 Scott Hannon      Add SATERR
C 12 Sep 2001 Scott Hannon      Add DFAFGL
C 20 Nov 2001 Scott Hannon      Remove NAMREF & add VKLAYE; v2.02
C  8 Mar 2002 Scott Hannon      Version to VKLAYE=2.03
C 29 Apr 2002 Scott Hannon      Version to VKLAYE=2.04
C 19 Feb 2007 Scott Hannon      Version to 2.05; add IDCLD,MXCLD, and
C                                  GUCCLD; change MXGAS to NGIDS+MXCLD

!END====================================================================

C      =================================================================
C      INCLUDE FILE incLAY
C      =================================================================
C
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
c       IMPLICIT NONE
c some stupid compilers choke on implicit none if it appears in both
c the include and the main code.

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
C      none (see description in body of code)

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
C      Declarations and parameters only

C      -----------------------------------------------------------------
C      Assign KLAYERS version string
C      -----------------------------------------------------------------
C      The version string consists of 3 parts: version number, date,
C      and a comment.  The version date should be updated to the
C      current date whenever any portion of the code is updated.  The
C      version number consists of two parts; a major version to the
C      left of the decimal point, and a minor version to the right.
C      The major number should be incremented only when major changes
C      have been made to the overall KLAYERS code.  The minor number
C      should be incremented only when minor but non-trivial changes
C      are made to the code.  Bug fixes should generally be handled
C      with the version date, but a fix for a serious bug may warrant
C      a change to the minor version number.
C      See the "Doc/last_update.txt" file for a description of the
C      changes associated with every change of VKLAYE.
C
       CHARACTER*40 VKLAYE
C      version template    '#.## YY-MM-DD <--------comment--------->'
       PARAMETER( VKLAYE = '2.05 11-08-30' )

C      -----------------------------------------------------------------
C      Assign generic variables
C      -----------------------------------------------------------------
       INTEGER MXIN         ! max number of input profile levels
       INTEGER NGRID        ! number of fine sub-grids per layer (10)
       INTEGER NGIDS        ! number of valid GAS IDs (see "cbgids.f")
       INTEGER MXCLD        ! max number of clouds
       INTEGER MXGAS        ! max number of gases+clouds for processing
C      Note: must have 1<=MXGAS<=NGIDS
C
       PARAMETER(MXIN = 600)
       PARAMETER(NGRID = 10)
       PARAMETER(NGIDS = 73)
       PARAMETER(MXCLD = 3)
       PARAMETER(MXGAS = NGIDS+MXCLD)


C      -----------------------------------------------------------------
C      Assign generic cloud "gas" ID number
C      -----------------------------------------------------------------
       INTEGER IDCLD1  ! cloud1 ID number
C      note: IDCLDn = IDCLD1 + n-1 for n=1 to MXCLD
       PARAMETER(IDCLD1 = 201)

C      -----------------------------------------------------------------
C      Assign how to interpret the water vapor mixing ratio
C      -----------------------------------------------------------------
C      Note: at runtime the code will inquire if the input profile
C      mixing ratios (for all gases) are for dry air (yes or no).  If
C      yes (ie dry), then the code will adjust/convert them to wet.
C      If no (ie wet), the code will not adjust the mixing ratios.
C      Variable MRTH2O may be used to force the water mixing ratio to
C      be dry or wet independent of what is specified for the other
C      gases, or to be of the same type as the other gases.
C
       INTEGER MRTH2O  ! code number for water mixing ratio type
C
C      Use the following code to set MRTH2O:
C         0 : same type as the other gases (recommended default)
C         1 : wet, ie MR will not be adjust (mimic old KLAYERS)
C         2 : dry, ie MR will be adjusted to convert to wet
C
       PARAMETER( MRTH2O = 1)

C      -----------------------------------------------------------------
C      Assign default file for AFGL models
C      -----------------------------------------------------------------
       CHARACTER*70 DFAFGL  ! Default file for AFGL models
C
       PARAMETER(DFAFGL=
c     $ '../Data/glatm.dat')
c     $ '/asl/packages/klayersV205/Data/glatm_9July2010.dat')
     $ '/asl/packages/klayersV205/Data/glatm.dat')

C      -----------------------------------------------------------------
C      Assign NBPLEV: number of boundary pressure levels defining layers
C      -----------------------------------------------------------------
       INTEGER NBPLEV  ! num of boundary pressure levels defining layers
C
C      Edit as needed
C      Must match PLEV asspecified in common block file cbplev.f
       PARAMETER(NBPLEV = 101)

C      -----------------------------------------------------------------
C      Assign parameters that depend on NBPLEV
C      -----------------------------------------------------------------
       INTEGER MYNLAY  ! number of layers = NLEV - 1
       INTEGER NSUBLV  ! max number of sublevels = MYNLAY*2*NGRID + 1
C
       PARAMETER(MYNLAY = NBPLEV - 1)
       PARAMETER(NSUBLV = MYNLAY*2*NGRID + 1)

C      -----------------------------------------------------------------
C      Assign max # of allowed Gas amount Units Code (GUC) numbers
C      -----------------------------------------------------------------
       INTEGER MAXGUC
       PARAMETER(MAXGUC = 9)
C      Must match GUCS assigned in cbgucs.f
C
       INTEGER GUCSTD  ! standard GUC for input
       INTEGER GUCOUT  ! gas units code for KLAYERS output
       PARAMETER(GUCSTD = 10) ! 10 = ppmv mixing ratio
c       PARAMETER(GUCOUT = 2)  !  2 = (Layer) kilomoles/cm^2
       PARAMETER(GUCOUT = 1)  !  1 = (Layer) molecules/cm^2
C      KLAYERS code must actually use these units!

       INTEGER GUCCLD ! gas units code number to use for cloud output
       PARAMETER(GUCCLD = 5)

C      -----------------------------------------------------------------
C      Assign standard CO2 mixing ratio
C      -----------------------------------------------------------------
       REAL CO2STD ! standard CO2 PPMV mixing ratio (385)
       PARAMETER( CO2STD = 385.0 )
C
C      -----------------------------------------------------------------
C      Assign minimum saturation error required for warning message
C      -----------------------------------------------------------------
       REAL SATERR ! percent min H2O over-saturation for warning
       PARAMETER( SATERR = 10.0 )

C      -----------------------------------------------------------------
C      Assign unit numbers for error and info/warning messages
C      -----------------------------------------------------------------
       INTEGER IOERR
       INTEGER IOINFO
       PARAMETER(IOERR = 6)
       PARAMETER(IOINFO = 6)
C
C      End of "incLAY.f" include file
