C=======================================================================
C=======================================================================
c for the Feb97 fast model
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              FTCOUT
C
!F77====================================================================

!ROUTINE NAME: FTCOUT

!ABSTRACT:
C    Write profile averaged & integrated alt/pres/temp/amount data to an
C    output file in the format used by our 100 layer AIRS_FTC package.

!CALL PROTOCOL:
C    FTCOUT(IOUNIT, IFIX, IWAT, IOZO, ICO, ICH4, LAT, LON, ZO,
C       NAMIN, TITLE, LZ, LDRY, LSVP, PLAY, TLAY, ALAY, ZLAY, DZLAY)

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  ALAY    avg gas amount              k.mol/cm^2
C    REAL arr  DZLAY   layer thickness             m
c    INTEGER   ICH4    array position of CH4       none
c    INTEGER   iCO     array position of CO        none
C    INTEGER   IFIX    array position of "fixed"   none
C    INTEGER   IOUNIT  output file unit number     none
C    INTEGER   IOZO    array position of ozone     none
C    INTEGER   IWAT    array position of water     none
C    REAL      LAT     latitude                    degrees
C    LOGICAL   LDRY    input mix ratios dry air?   none
C    REAL      LON     longitude                   degrees
C    LOGICAL   LSVP    checked/fixed water sat?    none
C    LOGICAL   LZ      input altitudes supplied?   none
C    CHAR*70   NAMIN   name of input prof file     none
C    REAL arr  PLAY    avg air pressure            atm
C    REAL arr  TLAY    avg air temperature         K
C    CHAR*80   TITLE   profile title or comment    none
C    REAL      Z0      altitude of 1100 mb level   m
C    REAL arr  ZLAY    layer altitudes             m

!OUTPUT PARAMETERS:
C    none

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

ccc!PARENT(S): TESTLAY
!PARENT(S): KLAYERS

!ROUTINES CALLED: none

!FILES ACCESSED:
C    Output text file, unit IOUNIT

!COMMON BLOCKS: none

!DESCRIPTION:
C    Writes a profile's average layer altitude, pressure, temperature,
C    & layer integrated absorber amounts to an output text file in a
C    format useable with our 100 layer AIRS_FTC package.
C
C    First it writes a "header" that describes some of the features of
C    the raw-profile-to-AIRS-layer-profile conversion that were (or
C    were not) used, along with some info about the original input file
C    and a comment. This header section need not follow any exact
C    format for use with the AIRS_FTC package; the *only* requirement
C    is that each line starts with a "!" character.
C
C    Then it loops over the AIRS layers, writing out each layer's
C    alt/pres/temp/amnts in turn, starting closest to the ground and
C    looping *upward*. This data must conform to the format:
C       Counter, Altitude, Pressure, Temperature, Fixed, Water, Ozone
C    where:
C       Counter = loop counter over layers (1=lowest, MYNLAY=highest)
C       Altitude = average layer altitude (meters)
C       Pressure = average layer pressure (atmospheres)
C       Temperature = average layer temperature (Kelvin)
C       Fixed = integrated layer "fixed" absorber amount (k.mol/cm^2)
C       Water = integrated layer water amount (k.mol/cm^2)
C       Ozone = integrated layer ozone amount (k.mol/cm^2)
C    The actual FORTRAN FORMAT (ie the way the numbers will be written)
C    is not important.

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS:
C    none

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C Apr  6 1995 Scott Hannon/UMBC created
C Jun 23 1995 Scott Hannon      Correct some comments
C Jul  3 1995 Scott Hannon      Added parameter DZLAY = layer thickness
C May 15 1997 Scott Hannon      Added CO & CH4 for Feb97 fast model
C 16 Jun 1998 Scott Hannon      CHanged TITLE from 40 to 80 characters


!END====================================================================


C      =================================================================
       SUBROUTINE FTCOUT(IOUNIT, IFIX, IWAT, IOZO, iCO, ICH4, LAT,
     $    LON, Z0, NAMIN, TITLE, LZ, LDRY, LSVP, PLAY, TLAY, ALAY,
     $    ZLAY, DZLAY)
C      =================================================================


C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
      include 'incLAY.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input parameters:
       INTEGER IOUNIT, IFIX, IWAT, IOZO, ICO, ICH4
       REAL LAT, LON, Z0, PLAY(MYNLAY), TLAY(MYNLAY),
     $      ALAY(MYNLAY,MXGAS), ZLAY(MYNLAY), DZLAY(MYNLAY)
       CHARACTER*70 NAMIN
       CHARACTER*80 TITLE
       LOGICAL LZ, LDRY, LSVP


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER L
       CHARACTER*70 CZ, CDRY, CSVP


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
C      ----------------------------------
C      Assign output strings for logicals
C      ----------------------------------
       IF (LZ) THEN
          CZ='Input altitudes: supplied & interpolated to layers'
       ELSE
          CZ='Input altitudes: not supplied; calculated instead'
       ENDIF
C
       IF (LDRY) THEN
          CDRY='Input mixing ratios treated as relative to: dry air'
       ELSE
          CDRY='Input mixing ratios treated as relative to: wet air'
       ENDIF
C
       IF (LSVP) THEN
          CSVP='H2O saturation: checked & corrected if needed'
       ELSE
          CSVP='H2O saturation: not checked'
       ENDIF
C
C      ----------------------------
C      Write the output file header
C      ----------------------------
       WRITE(IOUNIT,1040) LAT, LON, Z0, NAMIN, CZ, CDRY, CSVP, TITLE
ccc 1040  FORMAT('! Profile created using TESTLAY',
 1040  FORMAT('! Profile created using KLAYERS',
     $      /,'! lat=',F7.3,', lon=',F7.3,', 1100mb alt(m)=',1PE12.5,
     $      /,'! input=',A70,
     $      /,'! ',A70,
     $      /,'! ',A70,
     $      /,'! ',A70,
     $      /,'! no altitude(m)  thickness(m) pressure(atm) temp(K) ',
     $        '"fixed" amnt  water amnt   ozone amnt    CO amnt   ',
     $        '   CH4 amnt',
     $      /,'! -- ------------ ------------ ------------ -------- ',
     $        '------------ ------------ ------------ ------------',
     $        ' ------------',/,A80)
C
C      -----------------------------
C      Loop over the 100 AIRS layers
C      -----------------------------
       DO L=1,MYNLAY
          WRITE(IOUNIT,2050) L, ZLAY(L), DZLAY(L), PLAY(L), TLAY(L),
     $       ALAY(L,IFIX), ALAY(L,IWAT), ALAY(L,IOZO), ALAY(L,ICO),
     $       ALAY(L,ICH4)
 2050     FORMAT(I4,3(' ',1PE12.5),' ',0PF8.3,5(' ',1PE12.5))
       ENDDO
C
       RETURN
       END
