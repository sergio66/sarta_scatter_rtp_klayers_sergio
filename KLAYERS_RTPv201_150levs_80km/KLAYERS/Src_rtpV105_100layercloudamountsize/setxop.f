C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              SETXOP
C
!F77====================================================================


!ROUTINE NAME: SETXOP


!ABSTRACT:
C    Set the user-to-model cross-over pressure.


!CALL PROTOCOL:
C    SETXOP(NIN, NGASI, NGASO, LGASI, LGASO, PSURF, PIN, PROF)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   NIN     # of user prof pres levels  none
C    INTEGER   NGASI   # of gases in user prof     none
C    INTEGER   NGASO   # of gases in output PROF   none
C    INT arr   LGASI   list of user gases          HITRAN ID numbers
C    INT arr   LGASO   list of output gases        HITRAN ID numbers
C    REAL      PSURF   surface pressure            millibar
C    REAL arr  PIN     user prof pressure levels   millibars


!OUTPUT PARAMETERS:
C    none


!INPUT/OUTPUT PARAMETERS: none
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    STRUCT    PROF    RTP profile structure       (see attributes)

!RETURN VALUES: none


!PARENT(S): KLAYERS


!ROUTINES CALLED: none


!FILES ACCESSED:
C    none


!COMMON BLOCKS: none


!DESCRIPTION:
C    Checks the xover fields in "prof".  If less then zero, assumes
C    nodata and assigns a value.  If NIN=0, xover=PSURF.  If NIN>0,
C    then if the gas was in the user profile then xover=PIN(NIN),
C    else xover=PSURF.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C  8 Feb 2002 Scott Hannon      created


!END====================================================================

C      =================================================================
       SUBROUTINE SETXOP(NIN, NGASI, NGASO, LGASI, LGASO, PSURF, PIN,
     $    PROF)
C      =================================================================


C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
      include 'incLAY.f'
      include 'rtpdefs.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input parameters:
       INTEGER    NIN        ! number of levels in user profile
       INTEGER  NGASI        ! number of gases in user profile
       INTEGER  NGASO        ! number of gases in output PROF
       INTEGER LGASI( MXGAS) ! list of gases in user profile
       INTEGER LGASO( MXGAS) ! list of gases in output PROF
       REAL  PSURF           ! surface pressure
       REAL    PIN(  MXIN)   ! user profile pressure levels

C      Input/Output paramters:
C      Profile data structure
       RECORD /RTPPROF/ PROF
C      note: the only fields that *might* be modified are
C         PROF.txover and PROF.gxover(1 to NGASO)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER     I  ! generic looping variable
       INTEGER     J  ! generic looping variable
       INTEGER     K  ! generic
       REAL    XOP    ! user-to-AFGL cross-over pressure


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************

C      Determine user-to-AFGL xover pressure
       IF (NIN .EQ. 0) THEN
          XOP=PSURF
       ELSE
          XOP=PIN(NIN)
       ENDIF

C      Set PROF.txover if needed
       IF (PROF.txover .LE. 0.0) PROF.txover = XOP

C      Set PROF.gxover if needed
       DO I=1,NGASO
          IF (PROF.gxover(I) .LE. 0.0) THEN
C            Look for current LGASO gas in LGASI
             K=0
             DO J=1,NGASI
                IF (LGASI(J) .EQ. LGASO(I)) K=J
             ENDDO
             IF (K .EQ. 0) THEN
                PROF.gxover(I) = PSURF
             ELSE
                PROF.gxover(I) = XOP
             ENDIF
          ENDIF
       ENDDO
C
       RETURN
       END
