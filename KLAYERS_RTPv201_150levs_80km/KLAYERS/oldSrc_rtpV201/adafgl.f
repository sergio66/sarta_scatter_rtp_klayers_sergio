C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              ADAFGL
C
!F77====================================================================


!ROUTINE NAME: ADAFGL


!ABSTRACT:
C    ADjust the AFGL mixing ratio if necessary.


!CALL PROTOCOL:
C    ADAFGL(NAFGL,IDAFGL,NLEVAF,MIXAF)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   NAFGL   current index in MIXAF      none
C    INT arr   IDAFGL  ID of gases in MIXAF        none
C    INTEGER   NLEVAF  number of levels in MIXAF   none


!OUTPUT PARAMETERS: none


!INPUT/OUTPUT PARAMETERS: none
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  MIXAF   AFGL mixing ratios          PPMV


!RETURN VALUES: none


!PARENT(S): RDINFO


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    none


!COMMON BLOCKS:
C    none


!DESCRIPTION:
C    Checks the current gas ID and do some simple adjustment of the
C    mixing ratio if needed.  This routine is bascially a kludge to
C    keep the AFGL profile database more-or-less up to date as some
C    of the gas concentrations change over the years.


!ALGORITHM REFERENCES: see DESCRIPTION
C    http://cdiac.esd.ornl.gov/pns/current_ghg.html

!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 24 Jul 2001 Scott Hannon      Created; previously some of these
C                                  adjustments were done by "rdafgl.f".
C 05 Nov 2003 Scott Hannon      Fix cfc-13 & 113 (were switched); add
C                                  adjustments for cfc-22 and SF6


!END====================================================================


C      =================================================================
       SUBROUTINE ADAFGL( NAFGL, IDAFGL, NLEVAF, MIXAF )
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
C
C      Input:
       INTEGER  NAFGL          ! index for current gas in MAXAF
       INTEGER IDAFGL(MXGAS)   ! ID of current gas
       INTEGER NLEVAF          ! number of AFGL levels
C
C      Input/Output:
       REAL  MIXAF(MXIN,MXGAS) ! AFGL mixing ratios


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER IGAS
       INTEGER I
       REAL MRMULT


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************
C
C      ID of current gas
       IGAS=IDAFGL(NAFGL)
C
       MRMULT=1.0
C
C      Gas multipliers appropriate for year ~2000/2001

       IF (IGAS .EQ. 2) THEN
C         AFGL CO2 is 330 ppm; want CO2STD
          MRMULT=CO2STD/330.0
C
       ELSEIF (IGAS .EQ. 6) THEN
C         AFGL CH4 is 1.7ppm; want 1.8 ppm
          MRMULT=1.8/1.7
C
       ELSEIF(IGAS .EQ. 30) THEN
C         AFGL SF6 is 1.42 ppt; want ~4.8 ppt = 3.4 * 1.4 ppt
          MRMULT=3.4
C
       ELSEIF (IGAS .EQ. 51) THEN
C         AFGL CFC-11 is 140 ppt; want ~260 ppt = 1.9 * 140 ppt
          MRMULT=1.9
C
       ELSEIF(IGAS .EQ. 52) THEN
C         AFGL CFC-12 is 240 ppt; want ~540 ppt = 2.3 * 240 ppt
          MRMULT=2.3
C
ccc old (wrong) adjustment for cfc-13
c       ELSEIF(IGAS .EQ. 53) THEN
cC         AFGL CFC-13 is 5 ppt; want ~82 ppt = 17 * 5 ppt
c          MRMULT=17.0
ccc
       ELSEIF(IGAS .EQ. 56) THEN
C         AFGL HCFC-22 is 60 ppt; want ~146 ppt = 2.4 * 60 ppt
          MRMULT=2.4
C
       ELSEIF(IGAS .EQ. 57) THEN
C         AFGL CFC-113 is 19 ppt; want ~82 ppt = 4.3 * 19 ppt
          MRMULT=4.3
C
       ELSEIF(IGAS .EQ. 60) THEN
C         AFGL CCl4 is 130 ppt; want ~95 ppt = 0.74 * 130 ppt
          MRMULT=0.74
       ENDIF
C

C      Adjust MR if MR multipler is not unity
       IF (MRMULT .NE. 1.0) THEN
          DO I=1,NLEVAF
             MIXAF(I,NAFGL)=MIXAF(I,NAFGL)*MRMULT
          ENDDO
       ENDIF
C

       RETURN
       END
