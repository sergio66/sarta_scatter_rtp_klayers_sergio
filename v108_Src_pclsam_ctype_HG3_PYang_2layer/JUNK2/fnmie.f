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
C    FNMIE( MIETYP, FNMIEA, FNMIEE, FNMIEG )


!INPUT PARAMETERS:
C    none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INT arr   MIETYP  particle type code number   none
C    CHAR arr  FNMIEA  absorption filenames        none
C    CHAR arr  FNMIEE  extinction filenames        none
C    CHAR arr  FNMIEG  "g" asymmetry filenames     none


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
C    June 2008 version of the SARTA v1.08 code with PCLSAM Fast
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


!END====================================================================

C      =================================================================
       SUBROUTINE FNMIE ( MIETYP, FNMIEA, FNMIEE, FNMIEG )

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

c so this looks like ~/SARTA_CLOUDY/v107_Src_pclsam_ctype_HG3/
cC      mie particles #1: Baran ice aggregates
       MIETYP(1)=201
       FNMIEA(1)=
     $ '/asl/data/sarta_database/Data_mie_apr08/ice_aggr_abs.dat'
       FNMIEE(1)=
     $ '/asl/data/sarta_database/Data_mie_apr08/ice_aggr_ext.dat'
       FNMIEG(1)=
     $ '/asl/data/sarta_database/Data_mie_apr08/ice_aggr_asy.dat'
C
C      mie particles #2: water drop
       MIETYP(2)=101
       FNMIEA(2)=
     $ '/asl/data/sarta_database/Data_mie_apr08/waterdrop_abs.dat'
       FNMIEE(2)=
     $ '/asl/data/sarta_database/Data_mie_apr08/waterdrop_ext.dat'
       FNMIEG(2)=
     $ '/asl/data/sarta_database/Data_mie_apr08/waterdrop_asy.dat'
C

C      mie particles #3: volz CORRECT lognormal, with alpha = 2 
       MIETYP(3)=301 
       FNMIEA(3)= 
     $ '/home/sergio/MATLABCODE/JIGOU/AIRS2834/volz_1_0507_log2_abs.dat' 
       FNMIEE(3)= 
     $ '/home/sergio/MATLABCODE/JIGOU/AIRS2834/volz_1_0507_log2_ext.dat' 
       FNMIEG(3)= 
     $ '/home/sergio/MATLABCODE/JIGOU/AIRS2834/volz_1_0507_log2_asy.dat' 

C      mie particles #3: desertdust
c       MIETYP(3)=301
c       FNMIEA(3)=
c     $ '/asl/data/sarta_database/Data_mie_apr08/desertdust_abs.dat'
c       FNMIEE(3)=
c     $ '/asl/data/sarta_database/Data_mie_apr08/desertdust_ext.dat'
c       FNMIEG(3)=
c     $ '/asl/data/sarta_database/Data_mie_apr08/desertdust_asy.dat'
cC
 
       RETURN
       END
