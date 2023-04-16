C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              GETCLD
C
!F77====================================================================


!ROUTINE NAME: GETCLD


!ABSTRACT:
C    Get basic cloud info


!CALL PROTOCOL:
C    GETCLD( IPROF, HEAD, PROF,
C    LBLAC1, CTYPE1, CFRAC1, CPSIZ1, CPRTO1, CPRBO1, CNGWA1, CEMIS1, 
C    LBLAC2, CTYPE2, CFRAC2, CPSIZ2, CPRTO2, CPRBO2, CNGWA2, CEMIS2,
C    CFRA12, FCLEAR, CFRA1X, CFRA2X )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   IPROF   profile loop counter        none
C    STRUCT    HEAD    RTP header structure        various
C    STRUCT    PROF    RTP profile structure       various

!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    LOGICAL   LBLAC1  cloud1 is black?            none
C    INTEGER   CTYPE1  cloud type code number      none
C    REAL      CFRAC1  total cloud1 fraction       none (0.0 to 1.0)
C    REAL      CPSIZ1  particle size               um
C    REAL      CPRTO1  cloud top pressure          mb
C    REAL      CPRBO1  cloud bottom pressure       mb
C    REAL arr  CNGWA1  layer integrated profile    g/m^2
C    REAL      CEMIS1  emissivity                  none?
C    LOGICAL   LBLAC2  cloud2 is black?            none
C    INTEGER   CTYPE2  cloud type code number      none
C    REAL      CFRAC2  total cloud2 fraction       none (0.0 to 1.0)
C    REAL      CPSIZ2  particle size               um
C    REAL      CPRTO2  cloud top pressure          mb
C    REAL      CPRBO2  cloud bottom pressure       mb
C    REAL arr  CNGWA2  layer integrated profile    g/m^2
C    REAL      CEMIS2  emissivity                  none?
C    REAL      CFRA12  both clouds fraction        none (0.0 to 1.0)
C    REAL      FCLEAR  clear fraction              none (0.0 to 1.0)
C    REAL      CFRA1X  exclusive cloud1 fraction   none (0.0 to 1.0)
C    REAL      CFRA2X  exclusive cloud2 fraction   none (0.0 to 1.0)


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): SARTA


!ROUTINES CALLED: none


!FILES ACCESSED:
C    none


!COMMON BLOCKS: none


!DESCRIPTION:
C    Pulls out and checks cloud profile paramters.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 22 Feb 2007 Scott Hannon      created
C 10 Sep 2007 Scott Hannon      Mods for 100 layer particle size
C 14 Sep 2007 Scott Hannon      Add CPMIN check
C 15 Nov 2007 Scott Hannon      re-written for slab cloud

!END====================================================================

C      =================================================================
       SUBROUTINE GETCLD( IPROF, HEAD, PROF,
     $    LBLAC1, CTYPE1, CFRAC1, CPSIZ1, CPRTO1, CPRBO1, CNGWA1,CEMIS1,
     $    LBLAC2, CTYPE2, CFRAC2, CPSIZ2, CPRTO2, CPRBO2, CNGWA2,CEMIS2,
     $    CFRA12, FCLEAR, CFRA1X, CFRA2X )
C      =================================================================


C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
      include 'incFTC.f'
      include 'rtpdefs.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input parameters:
       INTEGER  IPROF          ! profile loop counter
       RECORD /RTPHEAD/ HEAD   ! header data
       RECORD /RTPPROF/ PROF   ! profile data

C      Output parameters:
C      cloud1
       LOGICAL LBLAC1      ! black cloud?
       INTEGER CTYPE1      ! cloud type code number
       REAL CFRAC1         ! cloud fraction
       REAL CPSIZ1         ! particle size (um)
       REAL CPRTO1         ! cloud top pressure (mb)
       REAL CPRBO1         ! cloud bottom pressure (mb)
       REAL CNGWA1         ! layer integrated cloud profile (g/m^2)
       REAL CEMIS1         ! emissivity
C      cloud2
       LOGICAL LBLAC2      ! black cloud?
       INTEGER CTYPE2      ! cloud type code number
       REAL CFRAC2         ! cloud fraction
       REAL CPSIZ2         ! particle size (um)
       REAL CPRTO2         ! cloud top pressure (mb)
       REAL CPRBO2         ! cloud bottom pressure (mb)
       REAL CNGWA2         ! layer integrated cloud profile (g/m^2)
       REAL CEMIS2         ! emissivity
C      other cloud fractions
       REAL CFRA12         ! both clouds fraction
       REAL FCLEAR         ! clear fraction
       REAL CFRA1X         ! exclusive cloud1 fraction
       REAL CFRA2X         ! exclusive cloud2 fraction

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       LOGICAL  LJUNK      ! generic

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************

C      cloud1 parameters
       CTYPE1=PROF.ctype
       CFRAC1=PROF.cfrac
       CPRTO1=PROF.cprtop
       CPRBO1=PROF.cprbot
       CNGWA1=PROF.cngwat
       CPSIZ1=PROF.cpsize
       CEMIS1=PROF.cemis
       IF (CFRAC1 .LT. 0.0 .OR. CTYPE1 .LT. 0) CFRAC1=0.0

C      cloud2 parameters
C      Note: currently uses generic RTP "udef" array variable
C      since RTP version 1.05 lacks standard cloud2 variables.
       CNGWA2=PROF.udef(11)
       CPSIZ2=PROF.udef(12)
       CPRTO2=PROF.udef(13)
       CPRBO2=PROF.udef(14)
       CFRAC2=PROF.udef(15)
       CFRA12=PROF.udef(16)
       CTYPE2=PROF.udef(17)
       CEMIS2=PROF.udef(18)
       IF (CFRAC2 .LT. 0.0 .OR. CTYPE2 .LT. 0) CFRAC2=0.0
       IF ( (CFRA12 .LT. 0.0) .OR.
     $    (CFRAC1 .EQ. 0.0 .OR. CFRAC2 .EQ. 0.0) ) CFRA12=0.0


C      Check cloud1 parameters
       IF (CFRAC1 .GT. 0.0) THEN
C
          IF (CFRAC1 .GT. 1.0) THEN
             WRITE(IOERR,1010) IPROF, 'CFRAC1', '0.0', '1.0'
 1010        FORMAT('Error! iprof=',I5,', ',A6,
     $       ' outside allowed ',A6,' to ',A6, ' range')
             STOP
          ENDIF
C
C         Check cprtop
          IF (CPRTO1 .LT. PROF.plevs(1) .OR.
     $        CPRTO1 .GT. PROF.spres) THEN
             WRITE(IOERR,1010) IPROF, 'CPRTO1', 'PLEVS1', 'SPRES'
             STOP
          ENDIF
C
          IF (CTYPE1 .LT. 100) THEN
C            Black cloud
             LBLAC1 = .TRUE.
             CPRBO1=CPRTO1*1.001 ! safe dummy value
             IF (CEMIS1 .LT. 0.0) CEMIS1=1.0  ! default is 1 not 0
C            Check cemis
             IF (CEMIS1 .GT. 1.0) THEN
                WRITE(IOERR,1010) IPROF, 'CEMIS1', '0.0', '1.0'
                STOP
             ENDIF
C
          ELSE
C            Slab cloud
             LBLAC1 = .FALSE.
C            Check cprbot, cpsize, & cngwat
             IF (CPRBO1 .LT. PROF.plevs(1) .OR.
     $           CPRBO1 .GT. PROF.spres) THEN
                WRITE(IOERR,1010) IPROF, 'CPRBO1', 'PLEVS1', 'SPRES'
                STOP
             ENDIF
             IF (CPRTO1 .GT. CPRBO1) THEN
                WRITE(IOERR,1020) IPROF, 'CPRTO1', 'CPRBO1'
 1020           FORMAT('Error! iprof=',I5,', ',A6,' > ',A6)
                STOP
             ENDIF
             IF (CPSIZ1 .LT. 0.0 .OR. CPSIZ1 .GT. 1E+3) THEN
                WRITE(IOERR,1010) IPROF, 'CPSIZ1', '0.0', '1E+3'
                STOP
             ENDIF
             IF (CNGWA1 .LT. 0.0 .OR. CNGWA1 .GT. 1E+6) THEN
                WRITE(IOERR,1010) IPROF, 'CNGWA1', '0.0', '1E+6'
                STOP
             ENDIF
          ENDIF

       ENDIF


C      Check cloud2 parameters
       IF (CFRAC2 .GT. 0.0) THEN
C
          IF (CFRAC2 .GT. 1.0) THEN
             WRITE(IOERR,1010) IPROF, 'CFRAC2', '0.0', '1.0'
             STOP
          ENDIF
C
C         Check cprtop
          IF (CPRTO2 .LT. PROF.plevs(1) .OR.
     $        CPRTO2 .GT. PROF.spres) THEN
             WRITE(IOERR,1010) IPROF, 'CPRTO2', 'PLEVS1', 'SPRES'
             STOP
          ENDIF
C
          IF (CTYPE2 .LT. 100) THEN
C            Black cloud
             LBLAC2 = .TRUE.
             CPRBO2=CPRTO2*1.001 ! safe dummy value
             IF (CEMIS2 .LT. 0.0) CEMIS2=1.0  ! default is 1 not 0
C            Check cemis
             IF (CEMIS2 .GT. 1.0) THEN
                WRITE(IOERR,1010) IPROF, 'CEMIS2', '0.0', '1.0'
                STOP
             ENDIF
C
          ELSE
C            Slab cloud
             LBLAC2 = .FALSE.
C            Check cprbot, cpsize, & cngwat
             IF (CPRBO2 .LT. PROF.plevs(1) .OR.
     $           CPRBO2 .GT. PROF.spres) THEN
                WRITE(IOERR,1010) IPROF, 'CPRBO2', 'PLEVS1', 'SPRES'
                STOP
             ENDIF
             IF (CPRTO2 .GT. CPRBO2) THEN
                WRITE(IOERR,1020) IPROF, 'CPRTO2', 'CPRBO2'
                STOP
             ENDIF
             IF (CPSIZ2 .LT. 0.0 .OR. CPSIZ2 .GT. 1E+3) THEN
                WRITE(IOERR,1010) IPROF, 'CPSIZ2', '0.0', '1E+3'
                STOP
             ENDIF
             IF (CNGWA2 .LT. 0.0 .OR. CNGWA2 .GT. 1E+6) THEN
                WRITE(IOERR,1010) IPROF, 'CNGWA2', '0.0', '1E+6'
                STOP
             ENDIF
          ENDIF
       ENDIF


C      Make sure top/only cloud uses cloud1 variables
       IF (CFRAC2 .GT. 0.0) THEN
          IF (CFRAC1 .EQ. 0.0) THEN
C            No cloud1 so shift cloud2 info to cloud1 variables
             CTYPE1=CTYPE2
             CFRAC1=CFRAC2
             CPRTO1=CPRTO2
             CPRBO1=CPRBO2
             CNGWA1=CNGWA2
             CPSIZ1=CPSIZ2
             CEMIS1=CEMIS2
             LBLAC1=LBLAC2
C            Clear CFRAC2
             CFRAC2=0.0
          ELSE
C            Two clouds
             IF (CPRTO1 .GT. CPRTO2) THEN
C               Switch variables so cloud1 is above cloud2
                CTYPE1=CTYPE2
                CFRAC1=CFRAC2
                CPRTO1=CPRTO2
                CPRBO1=CPRBO2
                CNGWA1=CNGWA2
                CPSIZ1=CPSIZ2
                CEMIS1=CEMIS2
                LJUNK=LBLAC1
                LBLAC1=LBLAC2
                CTYPE2=PROF.ctype
                CFRAC2=PROF.cfrac
                CPRTO2=PROF.cprtop
                CPRBO2=PROF.cprbot
                CNGWA2=PROF.cngwat
                CPSIZ2=PROF.cpsize
                CEMIS2=PROF.cemis
                LBLAC2=LJUNK
             ENDIF
          ENDIF
       ENDIF


C      Check CFRA12 and calc exclusive cloud fractions
       IF (CFRA12 .GT. 0.0) THEN
C
          IF (CFRA12 .GT. 1.0) THEN
             WRITE(IOERR,1010) IPROF, 'CFRA12', '0.0', '1.0'
             STOP
          ENDIF
C
          IF (LBLAC1) THEN
C            Note: If top cloud is black (trans=0) then must have cfra12=0
             WRITE(IOERR,1055) IPROF
 1055        FORMAT('Error! iprof=',I5,
     $       ' may not have CFRA12 > 0 with a top black cloud')
             STOP
          ELSE
c             IF (.NOT. LBLAC2 .AND. CPRBO1 .GT. CPRTO2) THEN
             IF (CPRBO1 .GT. CPRTO2) THEN
C               Bottom of top cloud is blow top of bottom cloud
                WRITE(IOERR,1020) IPROF, 'CPRBO1', 'CPRTO2'
             ENDIF
          ENDIF
C
          IF ((CFRA12 .LE. CFRAC1) .AND. (CFRA12 .LE. CFRAC2)) THEN
             CFRA1X=CFRAC1 - CFRA12
             CFRA2X=CFRAC2 - CFRA12
          ELSE
             WRITE(IOERR,1065) IPROF
 1065        FORMAT('Error! iprof=',I5,
     $       ' may not have cfra12 > cfrac1 or cfrac2')
          ENDIF
C
       ELSE
          CFRA1X=CFRAC1
          CFRA2X=CFRAC2
       ENDIF


C      Compute clear fraction
       FCLEAR=1.0 - (CFRA1X + CFRA2X + CFRA12)
       IF (FCLEAR .LT. 0.0) THEN
          IF (FCLEAR .GT. -1E-4) THEN
C            close enough to interpret as zero
             FCLEAR=0.0
          ELSE
             WRITE(IOERR,1070) IPROF, CFRA1X, CFRA2X, CFRA12, 1.0-FCLEAR
 1070        FORMAT('Error! iprof=',I5,', cfra1x=',F7.5,
     $       ' + cfra2x=',F7.5,' + cfra12=',F7.5,' =',F7.5,' > 1')
             STOP
          ENDIF
       ENDIF


       RETURN
       END
