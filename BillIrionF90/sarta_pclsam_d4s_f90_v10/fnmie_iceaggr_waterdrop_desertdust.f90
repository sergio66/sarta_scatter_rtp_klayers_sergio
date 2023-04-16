!=======================================================================
!
!    University of Maryland Baltimore Country (UMBC)
!
!    AIRS
!
!    FNMIE
!
!F77====================================================================


!ROUTINE NAME:
!    FNMIE


!ABSTRACT:0
!    Set Mie table filenames; this is a work-around for FORTRAN 77's
!    inability to assign an array of filenames in an include file
!    using PARAMETER.  This routine is essentially the executable
!    equivant of an include file containing a list of filenames.
!
!    Create/edit a new version of this file for every set of Mie
!    files you wish to use.  The number of entries (NMIETY) must
!    match the value in the "incFTC.f" include file.


!CALL PROTOCOL
!    FNMIE( VCLOUD, MIETYP, FNMIEA, FNMIEE, FNMIEG )


!INPUT PARAMETERS:
!    none


!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    CHAR*240  VCLOUD  cloud version string        none
!    INT arr   MIETYP  particle type code number   none
!    CHAR arr  FNMIEA  absorption filenames        none
!    CHAR arr  FNMIEE  extinction filenames        none
!    CHAR arr  FNMIEG  "g" asymmetry filenames     none


!INPUT/OUTPUT PARAMETERS:
!    none


!RETURN VALUES:
!    none


!PARENT(S):
!    RDCOEF


!ROUTINES CALLED:
!    none


!FILES ACCESSED:
!    incFTC.f : include file of parameter statements accessed during
!       compilation only.


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    May 2009 version of the SARTA v1.08 code with PCLSAM Fast
!    Transmittance Code by S.Hannon.
!
!    Copies individual filenames into an array of filenames.
!
!    Recommended Mie particle type code numbers (08 Jan 2007)
!    Min - Max  Description
!    ---------  --------------------------------------------------------
!    000 - 099  Black clouds
!    100 - 199  Spherical liquid H2O droplets
!    200 - 299  Ice aggregates
!    300 - 399  Dust/mineral
!    400 - 499  Sea salt (liquid H2O + salt)
!    500 - 599  Smoke/soot
!    600 - 699  Sulfate/pollutants


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
! Date        Programmer     Comments
! ----------- -------------- -------------------------------------------
! 24 Apr 2006 Scott Hannon   Created
! 08 Jan 2007 Scott Hannon   Change MIETYP code number values
! 12 May 2009 Scott Hannon   Add VCLOUD & CLDSTR
! 05 Sep 2017 Bill Irion     Conversion to F90
! 05 Sep 2017 Bill Irion     Added variable for base directory of particle coefficients


!END====================================================================

!      =================================================================
       SUBROUTINE FNMIE ( VCLOUD, MIETYP, FNMIEA, FNMIEE, FNMIEG )

!      =================================================================

       USE INCFTC
       IMPLICIT NONE


!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
!      Output
       CHARACTER*240 VCLOUD
       INTEGER MIETYP(NMIETY)
       CHARACTER*79 FNMIEA(NMIETY)
       CHARACTER*79 FNMIEE(NMIETY)
       CHARACTER*79 FNMIEG(NMIETY)


!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER I
       INTEGER IS
       INTEGER IE
       INTEGER LENNB
       INTEGER N
       CHARACTER*40 CLDSTR(NMIETY)


!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none


!***********************************************************************
!***********************************************************************
!                    EXECUTABLE CODE
!***********************************************************************
!***********************************************************************
!
!      Make sure NMIETY is as this routine expects
       IF (NMIETY .NE. 3) THEN
          WRITE(IOERR,1010) NMIETY,3
 1010     FORMAT('incFTC.f NMIETY=',I3,' but fnmie.f expects ',I3)
       ENDIF

! Scott had these in /asl/data/sarta_database/Data_mie_apr08/

!      mie particles #1: Baran ice aggregates
!                 123456789012345678901234567890123456789
!                 <------------cloud string------------->
       CLDSTR(1)='201=ice_aggr 2008-06-13 Baran'
       MIETYP(1)=201
       FNMIEA(1)='/asl/data/sarta_database/Data_mie_apr08/ice_aggr_abs.dat'
       FNMIEE(1)='/asl/data/sarta_database/Data_mie_apr08/ice_aggr_ext.dat'
       FNMIEG(1)='/asl/data/sarta_database/Data_mie_apr08/ice_aggr_asy.dat'

!

!      mie particles #2: water drop
!                 <------------cloud string------------->
       CLDSTR(2)='101=waterdrop 2008-06-13'
       MIETYP(2)=101
       FNMIEA(2)='/asl/data/sarta_database/Data_mie_apr08/waterdrop_abs.dat'
       FNMIEE(2)='/asl/data/sarta_database/Data_mie_apr08/waterdrop_ext.dat'
       FNMIEG(2)='/asl/data/sarta_database/Data_mie_apr08/waterdrop_asy.dat'

!

!      mie particles #3: desert dust
!                 <------------cloud string------------->
       CLDSTR(3)='301=desertdust 2008-06-13'
       MIETYP(3)=301
       FNMIEA(3)='/asl/data/sarta_database/Data_mie_apr08/desertdust_abs.dat'
       FNMIEE(3)='/asl/data/sarta_database/Data_mie_apr08/desertdust_ext.dat'
       FNMIEG(3)='/asl/data/sarta_database/Data_mie_apr08/desertdust_asy.dat'
!

!      Assign VCLOUD string
       IE=0
       DO I=1,NMIETY
          N=LENNB(CLDSTR(I))
          IS=IE + 1
          IE=IS + N - 1
          VCLOUD(IS:IE)=CLDSTR(I)
          IE=IE + 1
          IF (I .EQ. NMIETY) THEN
             VCLOUD(IE:IE)=CHAR(0)
          ELSE
             VCLOUD(IE:IE)=';'
          ENDIF
       ENDDO
!
       RETURN
       END
