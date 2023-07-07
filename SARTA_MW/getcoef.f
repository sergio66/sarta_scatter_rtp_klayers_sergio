      SUBROUTINE GETCOEF(LUNTR,IOS,FREQ)
C
C   NAME- GETCOEF   LANGUAGE- FORTRAN 77    PROGRAMMER- P. ROSENKRANZ
C
C   DATE- June 4, 2001
c         July 17, 2002 fixed initialization
C         Apr. 9, 2005 allow coefficients without profile
C         Mar. 29, 2006 maximum 8 coeff at each level
C         July 10, 2006 zero out any missing coeff's.
C
C   PURPOSE- READ TRANSMITTANCE COEFFICIENTS FOR MWTRAN INTO COMMON
c   The coefficient file should be opened and closed in the calling program.
C 
      IMPLICIT NONE
C   ARGUMENT SPECIFICATIONS-
      INTEGER LUNTR,IOS
      REAL FREQ(*)
C
C   NAME   I/O    DESCRIPTION
C
C   LUNTR   I     LOGICAL UNIT NUMBER TO READ TRANSMITTANCE COEFF FROM
C   IOS     O     FORTRAN ERROR CODE
C   FREQ    O     CHANNEL FREQUENCIES READ FROM FILE
C***********************************************************************
      INCLUDE 'mwtran.inc'
C
C   LOCAL VARIABLES
      CHARACTER*6 TYPEINST,DESCR,VERSION
      INTEGER NFREQ,NLEV,NCOEF,NDESCR,L,NCH,J,I,NCO1
      LOGICAL START
      DATA START/.TRUE./
      SAVE START
C
      IF(START) THEN
        NCHREAD = 0
        NTRLEV = 0
        START = .FALSE.
      ENDIF
C
C  READ HEADER AND SKIP OVER DESCRIPTION LINES FOR TRANS COEFF
      READ(LUNTR,4,IOSTAT=IOS,ERR=90) VERSION,TYPEINST,NFREQ,NLEV,
     & NCOEF,NDESCR
4     FORMAT(A5,A6,4I10)
      DO I=1,NDESCR
       READ(LUNTR,5,IOSTAT=IOS,ERR=90) DESCR
5      FORMAT(A6)
      END DO
C
      IF( NTRLEV.EQ.0 .AND. NLEV.GT.0 ) THEN 
!      READ STANDARD PROFILE FOR TRANSMITTANCE
       NTRLEV = NLEV
       DO L=1,NLEV
        READ(LUNTR,1,IOSTAT=IOS,ERR=90) PSTD(L),TSTD(L),PLB(L)
1       FORMAT(2F10.3,30X,F10.3)
       END DO
      ELSEIF( NLEV.GT.0 ) THEN 
!      SKIP STANDARD PROFILE, BECAUSE IT'S ALREADY IN COMMON
       DO L=1,NLEV
        READ(LUNTR,5,IOSTAT=IOS,ERR=90) DESCR
       END DO
      ENDIF
C
C  READ TRANSMITTANCE COEF
      DO 20 J=1,NFREQ
      READ(LUNTR,2,IOSTAT=IOS,ERR=90) FREQ(J),NCH,NLEV,NCOEF
2     FORMAT(F10.3,3I5)
      IF(NCHREAD.EQ.MAXCHAN) GOTO 90
      NCHREAD = NCHREAD + 1
      DO L=1,NLEV
       READ(LUNTR,3,IOSTAT=IOS,ERR=90) (TRCOEF(I,L,NCHREAD), I=1,NCOEF)
3      FORMAT(8E11.4)
       IF(NCOEF.LT.MAXCOEF) THEN
        NCO1 = NCOEF+1
        DO I=NCO1,MAXCOEF
          TRCOEF(I,L,NCHREAD) = 0.
        END DO
       ENDIF
      END DO
20    CONTINUE
C
90    RETURN
      END
