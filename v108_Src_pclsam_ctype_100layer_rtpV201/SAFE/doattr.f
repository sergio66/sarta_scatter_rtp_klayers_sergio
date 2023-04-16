c version for sarta
C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              DOATTR
C
!F77====================================================================


!ROUTINE NAME: DOATTR


!ABSTRACT:
C    Add/modify RTP HATT and PATT attribute strings


!CALL PROTOCOL:
C    DOATTR(HATT, PATT, LRHOT, VCLOUD)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    LOGICAL   LRHOT   force refl thermal rho?     none
C    CHAR arr  VCLOUD  cloud version string        none


!OUTPUT PARAMETERS:
C    none


!INPUT/OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    STRUCT    HATT    RTP header attributes       none
C    STRUCT    PATT    RTP profiles attributes     none


!RETURN VALUES: none


!PARENT(S): sarta


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    none


!COMMON BLOCKS: none


!DESCRIPTION:
C    Adds/changes RTP attribute strings


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 23 Feb 2007 Scott Hannon      Created (previously done in opnrtp.f)


!END====================================================================


C      =================================================================
       SUBROUTINE DOATTR(HATT, PATT, LRHOT, VCLOUD)
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
C      From "util.f"
C      function LENNB = length of string excluding trailing blanks


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input
       LOGICAL  LRHOT               ! force refl thermal rho=(1-e)/pi?
       CHARACTER*40 VCLOUD(NMIETY)  ! cloud version strings
C
C      Structures (see "rtpdefs.f")
       RECORD /RTPATTR/ HATT(MAXNATTR)  ! header attributes
       RECORD /RTPATTR/ PATT(MAXNATTR)  ! profile attributes


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      I
       INTEGER     IC
       INTEGER      J
       INTEGER      K
       INTEGER      N
       INTEGER  LENNB        ! for function LENNB
       INTEGER  NHATT        ! counter for # of header attributes
c       INTEGER  NPATT       ! counter for # of profile attributes
       CHARACTER*1      C1X  ! junk/work 1 char string
       CHARACTER*2      C2X  ! junk/work 2 char string
       CHARACTER*40   C40X1  ! junk/work 40 char string1
       CHARACTER*40   C40X2  ! junk/work 40 char string2
       CHARACTER*40   C40X3  ! junk/work 40 char string3
       CHARACTER*256  C256X  ! junk/work 256 char string
C

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
C      ---------------
C      HATT attributes
C      ---------------

C      Count the number of header attributes
       NHATT=1
       DO WHILE (ICHAR(HATT(I).fname) .NE. 0 .AND. I .LE. MAXNATTR)
          NHATT=NHATT + 1
       ENDDO

C      -----------------------------------------------------------------
C      "sarta" HATT entry
C
C      Assign the new "sarta" entry string
       C40X1=VSARTA   ! from incFTC.f
       C40X2=VSCOEF   ! from incFTC.f
       I=LENNB(C40X1)
       J=LENNB(C40X2)
       IF (LRHOT) THEN
          C1X='T'
       ELSE
          C1X='F'
       ENDIF
       C256X='src=' // C40X1(1:I) // '; coef=' // C40X2(1:J)
     $    // '; LRHOT=' // C1X
       J=LENNB(C256X)
C
C      Look for a previous "sarta" HATT entry
       IC=-1
       DO I=1,NHATT
          IF (HATT(I).aname(1:5) .EQ. 'sarta') IC=I
       ENDDO
       IF (IC .LT. 1) THEN
          NHATT=NHATT + 1
          IF (NHATT .GT. MAXNATTR) THEN
             WRITE(IOERR,1000)
 1000        FORMAT('MAXNATTR exceeded by doattr.f')
             STOP
          ENDIF
          IC=NHATT
       ENDIF
C
C      Plug in "sarta" HATT entry
       HATT(IC).fname='header'  // CHAR(0)
       HATT(IC).aname='sarta' // CHAR(0)
       HATT(IC).atext=C256X(1:J) // CHAR(0)
C

C      -----------------------------------------------------------------
C      Delete any/all old "vcloud<i>" entries
 10    IC=-1
       DO I=1,NHATT
          IF (HATT(I).aname(1:6) .EQ. 'vcloud') IC=I
       ENDDO
C
       IF (IC .GT. 0) THEN
          NHATT=NHATT - 1
          DO J=IC,NHATT
             HATT(J).fname=HATT(J+1).fname
             HATT(J).aname=HATT(J+1).aname
             HATT(J).atext=HATT(J+1).atext
          ENDDO
          GOTO 10
       ENDIF
C

C      -----------------------------------------------------------------
C      Add new "vcloud<i>" entries
C
       IF (NMIETY .GT. 99) THEN
          WRITE(IOERR,1010) NMIETY
 1010     FORMAT('Error: NMIETY=',I3,' > 99 max allowed by doattr.f')
       ENDIF
C      Loop over the cloud types
       DO K=1,NMIETY

C         "vcloud<i>" string
          IF (K .LT. 10) THEN
             WRITE(C1X,1021) K
 1021        FORMAT(I1)
             C40X3(1:7)='vcloud' // C1X
             N=7 ! length of C40X3
          ELSE
             WRITE(C2X,1022) K
 1022        FORMAT(I2)
             C40X3(1:8)='vcloud' // C2X
             N=8 ! length of C40X3
          ENDIF
C
C         Assign the new "vcloud<i>" entry string
          C40X1=VCLOUD(K)
          I=LENNB(C40X1)
          C256X='cloud table ' // C40X1(1:I)
C
C         Plug in "vcloud<i>" HATT entry
          NHATT=NHATT + 1
          IF (NHATT .GT. MAXNATTR) THEN
             WRITE(IOERR,1000)
             STOP
          ENDIF
          HATT(NHATT).fname='header'  // CHAR(0)
          HATT(NHATT).aname=C40X3(1:N) // CHAR(0)
          J=LENNB(C256X)
          HATT(NHATT).atext=C256X(1:J) // CHAR(0)
C
       ENDDO
C
C      -----------------------------------------------------------------
C
C      Add a char(0) to end of attributes if less than maxnattr
       IF (NHATT .LT. MAXNATTR) THEN
          HATT(NHATT + 1).fname=CHAR(0)
       ENDIF
C

C      ---------------
C      PATT attributes
C      ---------------
ccc do not bother to check profile attributes cccccccccccccccccccccccccc
C
       RETURN
       END
