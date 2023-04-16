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
C    GETCLD( IPROF, ICLD1, ICLD2, HEAD, PROF,
C    LBLAC1, CTYPE1, CFRAC1, CPSIZ1, CPRTO1, CEMIS1, CNGWA1,
C    LBLAC2, CTYPE2, CFRAC2, CPSIZ2, CPRTO2, CEMIS2, CNGWA2,
C    CFRA12, FCLEAR, CFRA1X, CFRA2X )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   IPROF   profile loop counter        none
C    INTEGER   ICLD1   cld1 index in HEAD.glist    none
C    INTEGER   ICLD2   cld2 index in HEAD.glist    none
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
C    REAL      CEMIS1  emissivity                  none?
C    REAL arr  CNGWA1  layer integrated profile    g/m^2
C    LOGICAL   LBLAC2  cloud2 is black?            none
C    INTEGER   CTYPE2  cloud type code number      none
C    REAL      CFRAC2  total cloud2 fraction       none (0.0 to 1.0)
C    REAL      CPSIZ2  particle size               um
C    REAL      CPRTO2  cloud top pressure          mb
C    REAL      CEMIS2  emissivity                  none?
C    REAL arr  CNGWA2  layer integrated profile    g/m^2
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


!END====================================================================

C      =================================================================
       SUBROUTINE GETCLD( IPROF, ICLD1, ICLD2, HEAD, PROF,
     $    LBLAC1, CTYPE1, CFRAC1, CPSIZ1, CPRTO1, CEMIS1, CNGWA1,
     $    LBLAC2, CTYPE2, CFRAC2, CPSIZ2, CPRTO2, CEMIS2, CNGWA2,
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
       INTEGER  IPROF                   ! profile loop counter
       INTEGER  ICLD1                   ! index of cloud1 in HEAD.glist
       INTEGER  ICLD2                   ! index of cloud2 in HEAD.glist
       RECORD /RTPHEAD/ HEAD            ! header data
       RECORD /RTPPROF/ PROF            ! profile data

C      Output parameters:
C      cloud1
       LOGICAL LBLAC1      ! black cloud?
       INTEGER CTYPE1      ! cloud type code number
       REAL CFRAC1         ! cloud fraction
       REAL CPSIZ1         ! particle size (um)
       REAL CPRTO1         ! cloud top pressure (mb)
       REAL CEMIS1         ! emissivity
       REAL CNGWA1(MAXLAY) ! layer integrated cloud profile (g/m^2)
C      cloud2
       LOGICAL LBLAC2      ! black cloud?
       INTEGER CTYPE2      ! cloud type code number
       REAL CFRAC2         ! cloud fraction
       REAL CPSIZ2         ! particle size (um)
       REAL CPRTO2         ! cloud top pressure (mb)
       REAL CEMIS2         ! emissivity
       REAL CNGWA2(MAXLAY) ! layer integrated cloud profile (g/m^2)
C      other cloud fractions
       REAL CFRA12         ! both clouds fraction
       REAL FCLEAR         ! clear fraction
       REAL CFRA1X         ! exclusive cloud1 fraction
       REAL CFRA2X         ! exclusive cloud2 fraction

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      I      ! dummy looping
       INTEGER      L      ! layer looping
       INTEGER     LR      ! reverse layer looping
       INTEGER   NLAY      ! number of layers to read
       INTEGER   NLEV      ! number of levels

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************

C      Number of layers
       NLEV=PROF.nlevs
       IF (HEAD.ptype .EQ. AIRSLAY) THEN
C         Special case for AIRS pseudo-levels
          NLAY=NLEV
       ELSE
          NLAY=NLEV - 1
       ENDIF

C      cloud1 parameters
       CFRAC1=PROF.cfrac
       CTYPE1=PROF.ctype
       IF (CFRAC1 .LT. CFRMIN) CFRAC1=0.0
       IF (CFRAC1 .GT. 1.0) THEN
          IF (CFRAC1 .LT. 1.00001) THEN
             CFRAC1=1.0
          ELSE
             WRITE(IOERR,1000) IPROF, 1
 1000        FORMAT('Error, prof=',I6,': cfrac',I1,' > 1.0')
             STOP
          ENDIF
       ENDIF
C

       IF (CTYPE1 .GE. 100) THEN
C         Transmissive cloud; load profile
C
          LBLAC1=.FALSE.
          CPSIZ1=PROF.cpsize
C
          IF (ICLD1 .LT. 1) THEN
C            Required cloud profile missing
             WRITE(IOERR,1010) IPROF, 1, CTYPE1, 1, IDCLD1
 1010        FORMAT('Error, prof=',I6,': ctype',I1,'=',I4,
     $       ' but HEAD.glist lacks cloud',I1,'=',I3)
             STOP
          ELSE
             IF (PROF.plevs(1) .LT. PROF.plevs(NLEV)) THEN
C               Prof is in top-down order
                DO L=1,NLAY
                   CNGWA1(L)=AMAX1(0.0,PROF.gamnt(L,ICLD1))
                ENDDO
             ELSE
C               Prof is in bottom-up order
                DO L=1,NLAY
                   LR=1 + NLAY - L  ! reversed layer index
                   CNGWA1(L)=AMAX1(0.0,PROF.gamnt(LR,ICLD2))
                ENDDO
             ENDIF
          ENDIF
       ELSE
C         Black cloud; load cloudtop surface
          LBLAC1=.TRUE.
          CPRTO1=PROF.cprtop
          DO I=1,PROF.nemis
            CEMIS1(I)=PROF.cemis(I)
          END DO
       ENDIF
C

C      cloud2 parameters
       CFRAC2=PROF.udef(15)
       CTYPE2=PROF.udef(17)
       IF (CFRAC2 .LT. CFRMIN) CFRAC2=0.0
       IF (CFRAC2 .GT. 1.0) THEN
          IF (CFRAC2 .LT. 1.0001) THEN
             CFRAC2=1.0
          ELSE
             WRITE(IOERR,1000) IPROF, 2
             STOP
          ENDIF
       ENDIF
C

       IF (CTYPE2 .GE. 100) THEN
C         Transmissive cloud; load profile
C
          LBLAC2=.FALSE.
          CPSIZ2=PROF.udef(12)
C
          IF (ICLD2 .LT. 1) THEN
C            Required cloud profile missing
             WRITE(IOERR,1010) IPROF, 2, CTYPE2, 2, IDCLD2
             STOP
          ELSE
             IF (PROF.plevs(1) .LT. PROF.plevs(NLEV)) THEN
C               Prof is in top-down order
                DO L=1,NLAY
                   CNGWA2(L)=AMAX1(0.0,PROF.gamnt(L,ICLD2))
                ENDDO
             ELSE
C               Prof is in bottom-up order
                DO L=1,NLAY
                   LR=1 + NLAY - L  ! reversed layer index
                   CNGWA2(L)=AMAX1(0.0,PROF.gamnt(LR,ICLD2))
                ENDDO
             ENDIF
          ENDIF

       ELSE
C         Black cloud; load cloudtop surface
          LBLAC2=.TRUE.
          CPRTO2=PROF.udef(13)
          DO I=1,PROF.nemis
            CEMIS2(I)=PROF.udef(18,I)
          END DO
       ENDIF
C

C      Both clouds fraction
       CFRA12=PROF.udef(16)
       IF (CFRA12 .LT. CFRMIN) CFRA12=0.0
       IF (CFRA12 .GT. 1.0) THEN
          IF (CFRA12 .LT. 1.0001) THEN
             CFRA12=1.0
          ELSE
             WRITE(IOERR,1030) IPROF
 1030        FORMAT('Error, prof=',I6,': cfra12 > 1.0')
             STOP
          ENDIF
       ENDIF
C

       IF (CFRA12 .GT. CFRAC1) THEN
          IF (CFRA12/CFRAC1 .LT. 1.0001) THEN
             CFRA12=CFRAC1
          ELSE
             WRITE(IOERR,1032) IPROF, 1
 1032        FORMAT('Error, prof=',I6,': cfra12 > cfrac',I1)
             STOP
          ENDIF
       ENDIF
C

       IF (CFRA12 .GT. CFRAC2) THEN
          IF (CFRA12/CFRAC2 .LT. 1.0001) THEN
             CFRA12=CFRAC2
          ELSE
             WRITE(IOERR,1030) IPROF, 2
          ENDIF
       ENDIF
C

       IF ((CFRA12 .GT. 0.0) .AND. (LBLAC1 .AND. LBLAC2)) THEN
          WRITE(IOERR,1035) IPROF
 1035     FORMAT('Error, iprof=',I5,
     $    ': may not have cfra12 > 0 with two black clouds')
          STOP
       ENDIF
C

C      Compute exclusive cloud fractions
       IF (CFRA12 .GT. 0.0) THEN
          CFRA1X=CFRAC1 - CFRA12
          CFRA2X=CFRAC2 - CFRA12
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
             WRITE(IOERR,1040) IPROF, CFRA1X, CFRA2X, CFRA12, 1.0-FCLEAR
 1040        FORMAT('Error! iprof=',I5,', cfra1x=',F7.5,
     $       ' + cfra2x=',F7.5,' + cfra12=',F7.5,' =',F7.5,' > 1')
             STOP
          ENDIF
       ENDIF
C
       RETURN
       END
