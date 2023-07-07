C  MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C  AIRS
C
C  MICROWAVE FIRST-GUESS RETRIEVAL
C
!ROUTINE NAME: TB10
!CALL INTERFACE:
      SUBROUTINE TB10 (TD,TR,E,ER,M,T,TRAN,TRANR)
C
!F77   LANGUAGE- FORTRAN 77
C                                                                       
!ABSTRACT: COMPUTES THE ATMOSPHERIC COMPONENTS OF BRIGHTNESS
C   TEMPERATURE
C
!ROUTINE HISTORY:
C  VERSION- 1 DATE- 4/5/94   PROGRAMMER-P.ROSENKRANZ
c    June 8, 2000 pwr - separated direct and reflected paths
C
!ARGUMENTS:
C  SPECIFICATIONS-
      IMPLICIT NONE
      INTEGER M
      REAL TD,TR,E,ER,T(M),TRAN(M),TRANR(M)
C  NAME     TYPE    I/O    DESCRIPTION
C
C  TD       R*4      O     UPWARD PROPAGATING ATMOSPHERIC EMISSION SEEN FROM
C                           LEVEL(1) ALONG DIRECT PATH, IN DEG K
C  TR       R*4      O     DOWNWARD PROPAGATING ATMOSPHERIC EMISSION SEEN FROM
C                           LEVEL(M) ALONG REFLECTED PATH, IN DEG K; 
C                           (NOT INCLUDING COSMIC BACKGROUND).
C  E        R*4      O     ONE WAY TRANSMITTANCE OF THE ATMOSPHERE, ALONG
C                           THE DIRECT PATH.
C  ER       R*4      O     ONE WAY TRANSMITTANCE OF THE ATMOSPHERE, ALONG
C                           THE REFLECTED PATH.
C  M        I*4      I     NUMBER OF ELEMENTS IN T,TRAN
C  T        R*4      I     VECTOR OF TEMPERATURES AT LEVELS (DEG K).
C  TRAN     R*4      I     LAYER TRANSMITTANCE BETWEEN LEVEL I-1 AND I ALONG
C                           THE DIRECT PATH
C  TRANR    R*4      I     LAYER TRANSMITTANCE BETWEEN LEVEL I-1 AND I ALONG
C                           THE REFLECTED PATH
C
!ROUTINES CALLED: none
!PARENT: 
!RETURN VALUES:
!FILES ACCESSED: none
!DESCRIPTION:
!KNOWN BUGS AND LIMITATIONS: 
!END HEADER*************************************************************
C  LOCAL VARIABLES:
      REAL EM,ERM,TAV
      INTEGER I
C
      E = TRAN(1)
      TD = T(1)*(1.-E) ! contrib from top layer
      DO I=2,M
        IF(E.LT.1.E-10) GOTO 50
        EM = E
        E = E*TRAN(I)
        TAV = (T(I) + T(I-1))/2.
        TD = TD + TAV*(EM-E)
      END DO
50    CONTINUE
      ER = 1.
      TR = 0.
      DO I=M,2,-1
        IF(ER.LT.1.E-10) GOTO 60
        ERM = ER
        ER = ER*TRANR(I)
        TAV = (T(I) + T(I-1))/2.
        TR = TR + TAV*(ERM-ER)
      END DO
      ERM = ER
      ER = ER*TRANR(1)
      TR = TR + T(1)*(ERM-ER) ! contrib from top layer
60    CONTINUE
      RETURN
      END
