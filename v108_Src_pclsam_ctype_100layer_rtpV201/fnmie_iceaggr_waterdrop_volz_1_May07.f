C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore Country (UMBC)
C
C    AIRS
C
C    FNMIE
C
!F77====================================================================


!ROUTINE NAME:
C    FNMIE


!ABSTRACT:0
C    Set Mie table filenames; this is a work-around for FORTRAN 77's
C    inability to assign an array of filenames in an include file
C    using PARAMETER.  This routine is essentially the executable
C    equivant of an include file containing a list of filenames.
C
C    Create/edit a new version of this file for every set of Mie
C    files you wish to use.  The number of entries (NMIETY) must
C    match the value in the "incFTC.f" include file.


!CALL PROTOCOL
C    FNMIE( MIETYP, FNMIEA, FNMIEE, FNMIEG, VCLOUD )


!INPUT PARAMETERS:
C    none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INT arr   MIETYP  particle type code number   none
C    CHAR arr  FNMIEA  absorption filenames        none
C    CHAR arr  FNMIEE  extinction filenames        none
C    CHAR arr  FNMIEG  "g" asymmetry filenames     none
C    CHAR arr  VCLOUD  cloud version/comment       none


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    RDCOEF


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    April 2006 version of the SARTA v1.07 code with PCLSAM Fast
C    Transmittance Code by S.Hannon.
C
C    Copies individual filenames into an array of filenames.
C
C    Recommended Mie particle type code numbers (08 Jan 2007)
C    Min - Max  Description
C    ---------  --------------------------------------------------------
C    000 - 099  Black clouds
C    100 - 199  Spherical liquid H2O droplets
C    200 - 299  Ice aggregates
C    300 - 399  Dust/mineral
C    400 - 499  Sea salt (liquid H2O + salt)
C    500 - 599  Smoke/soot
C    600 - 699  Sulfate/pollutants


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date        Programmer     Comments
C    ----------- -------------- ----------------------------------------
C    24 Apr 2006 Scott Hannon   Created
C    08 Jan 2007 Scott Hannon   Change MIETYP code number values
C    23 Feb 2007 Scott Hannon   Add VCLOUD


!END====================================================================

C      =================================================================
       SUBROUTINE FNMIE ( MIETYP, FNMIEA, FNMIEE, FNMIEG, VCLOUD )

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
C      Output
       INTEGER MIETYP(NMIETY)
       CHARACTER*79 FNMIEA(NMIETY)
       CHARACTER*79 FNMIEE(NMIETY)
       CHARACTER*79 FNMIEG(NMIETY)
       CHARACTER*40 VCLOUD(NMIETY)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
C      none


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
C      Make sure NMIETY is as this routine expects
       IF (NMIETY .NE. 3) THEN
          WRITE(IOERR,1010) NMIETY,3
 1010     FORMAT('incFTC.f NMIETY=',I3,' but fnmie.f expects ',I3)
       ENDIF

C      mie particles #1: Baran ice aggregates
       MIETYP(1)=201
       FNMIEA(1)=
     $ '/asl/data/sarta_database/Data_mie_m135f/New/ice_aggr_abs.dat'
       FNMIEE(1)=
     $ '/asl/data/sarta_database/Data_mie_m135f/New/ice_aggr_ext.dat'
       FNMIEG(1)=
     $ '/asl/data/sarta_database/Data_mie_m135f/New/ice_aggr_asy.dat'
C        40 char <1234567890123456789012345678901234567890>
       VCLOUD(1)='201 ice_aggr 2004-02-18'
C      Recommended vcloud: MIETYP, name, date created
C
C      mie particles #2: water drop
       MIETYP(2)=101
       FNMIEA(2)=
     $ '/asl/data/sarta_database/Data_mie_m135f/New/waterdrop_abs.dat'
       FNMIEE(2)=
     $ '/asl/data/sarta_database/Data_mie_m135f/New/waterdrop_ext.dat'
       FNMIEG(2)=
     $ '/asl/data/sarta_database/Data_mie_m135f/New/waterdrop_asy.dat'
C        40 char <1234567890123456789012345678901234567890>
       VCLOUD(2)='101 waterdrop 2004-02-18'
C
C      mie particles #3: volz; correct lognormal, correct log(alpha) = log(2)
       MIETYP(3)=301
       FNMIEA(3)= 
     $ '/home/sergio/MATLABCODE/JIGOU/DATA/volz_1_0507_log2_abs.dat' 
       FNMIEE(3)= 
     $ '/home/sergio/MATLABCODE/JIGOU/DATA/volz_1_0507_log2_ext.dat' 
       FNMIEG(3)= 
     $ '/home/sergio/MATLABCODE/JIGOU/DATA/volz_1_0507_log2_asy.dat' 
C        40 char <1234567890123456789012345678901234567890>
       VCLOUD(3)='301 volz 2007-05-01'
C
       RETURN
       END
