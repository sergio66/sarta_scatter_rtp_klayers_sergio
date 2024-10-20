C----------------------------------------------------------------------

C----------------------------------------------------------------------
C
       SUBROUTINE SPLINT(xa,ya,y2a,n,x,y)
C
C      Spline interpolation using the spline coefs determined by spline.
C      -----------------------------------------------------------------
C
      include 'incLAY.f'
C
       INTEGER n
       REAL x,y,xa(n),y2a(n),ya(n)
       INTEGER k,khi,klo
       REAL a,b,h
C
       klo=1
       khi=n
 1     IF (khi-klo .GT. 1) THEN
          k=(khi+klo)/2
          IF(xa(k) .GT. X) THEN
             khi=k
          ELSE
             klo=k
          ENDIF
          GOTO 1
       ENDIF
C
       h=xa(khi)-xa(klo)
       IF (h .EQ. 0.0) THEN
          WRITE(IOERR,1010)
 1010     FORMAT('Error! bad xa input in splint. Quitting')
          STOP
       ENDIF
C

       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*
     $    (h**2)/6.0
C
       RETURN
       END
