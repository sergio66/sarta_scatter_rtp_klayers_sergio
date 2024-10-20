C
       SUBROUTINE WRTSPL(IOSPL,FNAME,NIN,PIN,TIN,NGASES,GASID,MRIN,
     $    NSIN,PSIN,TSIN,MRSIN,LAFGL,NMODEL)
C
      include 'incLAY.f'
C
C      Parameters
       INTEGER IOSPL,NIN,NGASES,GASID(MXGAS),NSIN,NMODEL
       REAL PIN(MXIN),TIN(MXIN),MRIN(MXIN,MXGAS),
     $    PSIN(MXIN),TSIN(MXIN),MRSIN(MXIN,MXGAS)
       CHARACTER*70 FNAME
       LOGICAL LAFGL
C
C      Local variables
       INTEGER I,J,K,II
C
C      Write the interpolated profile to the spline out file
       WRITE(IOSPL,1010)
 1010  FORMAT('! ----------------------------------------')
C
       IF (LAFGL) THEN
          WRITE(IOSPL,1020) NMODEL,FNAME,(GASID(I),I=1,NGASES)
 1020     FORMAT('! AFGL profile model ',I1,/,'! ',A70,/,
     $       '! gas IDs:',15(' ',I2))
       ELSE
          WRITE(IOSPL,1025) FNAME,(GASID(I),I=1,NGASES)
 1025     FORMAT('! User profile=',/,'! ',A70,/,
     $       '! gas IDs:',15(' ',I2))
       ENDIF
C
       WRITE(IOSPL,1030)
 1030  FORMAT('! If ?=1, level is interpolated',/,
     $       '! # ?  pressure     temperature  mixing ratios...')
C
       J=1
       DO I=1,NSIN
C
C         Determine whether or not the current level was interpolated
          IF (PSIN(I) .EQ. PIN(J)) THEN
C            Non-interpolated level
             K=0
             IF (J .LT. NIN) J=J+1
          ELSE
C            Interpolated level
             K=1
          ENDIF
C
C         Write data for this level to output file
          WRITE(IOSPL,1040) I, K, PSIN(I), TSIN(I),
     $       (MRSIN(I,II),II=1,NGASES)
 1040     FORMAT(I3,' ',I1,' ',1PE12.5,' ',0PF7.2,15(' ',1PE10.3))
C
       ENDDO
C
       RETURN
       END
