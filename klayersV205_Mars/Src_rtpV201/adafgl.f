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
C 14 Jul 2010 Scott Hannon      Add adjustments for CH3Cl, cfc-13,
C                                  hcfc-21, cfc-114, cfc-115; update
C                                  adjustments for SF6, hcfc-22,
C                                  cfc-113, cfc-12, cfc-11, CCl4
C 16 Jul 2010 Scott Hannon      Add xsec 81 for SF5 (same as gas 30)
C 16 Aug 2010 Scott Hannon      Adjust SF6 (30 & 81) for revised profile

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
C      Gas multipliers appropriate for year 2009 (updated from 2001)

CCCC UNCOMMENT 
CCCC UNCOMMENT        IF (IGAS .EQ. 2) THEN
CCCC UNCOMMENT C         AFGL CO2 is 330 ppm; want CO2STD
CCCC UNCOMMENT           MRMULT=CO2STD/330.0
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF (IGAS .EQ. 6) THEN
CCCC UNCOMMENT C         AFGL CH4 is 1.7ppm; want 1.8 ppm
CCCC UNCOMMENT           MRMULT=1.8/1.7
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 24) THEN
CCCC UNCOMMENT C         AFGL CH3Cl is 700 ppt; want ~540 ppt = 0.77 * 700 ppt
CCCC UNCOMMENT           MRMULT=0.77
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 30) THEN
CCCC UNCOMMENT ccc adjustment for old profile
CCCC UNCOMMENT cC         AFGL SF6 is 1.42 ppt; want ~6.6 ppt = 4.6 * 1.4 ppt
CCCC UNCOMMENT c          MRMULT=4.6  ! updated from 3.4
CCCC UNCOMMENT ccc
CCCC UNCOMMENT C         trop SF6 is 5.0 ppt; want ~6.6 ppt = 1.32 * 5.0 ppt
CCCC UNCOMMENT           MRMULT=1.32
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 51) THEN
CCCC UNCOMMENT C         AFGL CFC-11 is 140 ppt; want ~242 ppt = 1.7 * 140 ppt
CCCC UNCOMMENT           MRMULT=1.7  ! updated from 1.9
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 52) THEN
CCCC UNCOMMENT C         AFGL CFC-12 is 240 ppt; want ~536 ppt = 2.2 * 240 ppt
CCCC UNCOMMENT           MRMULT=2.2  ! updated from 2.3
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 53) THEN
CCCC UNCOMMENT C         AFGL CFC-13 is 5 ppt; want ~3.5 ppt = 0.7 * 5 ppt
CCCC UNCOMMENT           MRMULT=0.7
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 55) THEN
CCCC UNCOMMENT C         AFGL HCFC-21 is 1.6 ppt; want ~0.3 ppt = 0.2 * 1.6 ppt
CCCC UNCOMMENT           MRMULT=0.2
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 56) THEN
CCCC UNCOMMENT C         AFGL HCFC-22 is 60 ppt; want ~198 ppt = 3.3 * 60 ppt
CCCC UNCOMMENT           MRMULT=3.3  ! updated from 2.4
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 57) THEN
CCCC UNCOMMENT C         AFGL CFC-113 is 19 ppt; want ~76 ppt = 4.0 * 19 ppt
CCCC UNCOMMENT           MRMULT=4.0  ! updated from 4.3
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 58) THEN
CCCC UNCOMMENT C         AFGL CFC-114 is 12 ppt; want ~17 ppt = 1.4 * 12 ppt
CCCC UNCOMMENT           MRMULT=1.4
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 59) THEN
CCCC UNCOMMENT C         AFGL CFC-115 is 4 ppt; want ~8 ppt = 2.0 * 4 ppt
CCCC UNCOMMENT           MRMULT=2.0
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 60) THEN
CCCC UNCOMMENT C         AFGL CCl4 is 130 ppt; want ~88 ppt = 0.68 * 130 ppt
CCCC UNCOMMENT           MRMULT=0.68  ! updated from 0.74
CCCC UNCOMMENT C
CCCC UNCOMMENT        ELSEIF(IGAS .EQ. 81) THEN
CCCC UNCOMMENT C         trop SF6 is 5.0 ppt; want ~6.6 ppt = 1.32 * 5.0 ppt
CCCC UNCOMMENT           MRMULT=1.32
CCCC UNCOMMENT        ENDIF
CCCC UNCOMMENT C
CCCC UNCOMMENT 

C      Adjust MR if MR multipler is not unity
       IF (MRMULT .NE. 1.0) THEN
          DO I=1,NLEVAF
             MIXAF(I,NAFGL)=MIXAF(I,NAFGL)*MRMULT
          ENDDO
       ENDIF
C

       RETURN
       END
