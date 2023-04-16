c this is copied from rad_misc.f in kCARTA

      SUBROUTINE GaussLegendre(nn,raX,raW) 
        
      IMPLICIT NONE 
 
      include 'incFTC.f'
      include 'pyang2cld.param'
 
      INTEGER nn 
      REAL raX(NANG_MAX),raW(NANG_MAX) 

      DOUBLE PRECISION daX(NANG_MAX),daW(NANG_MAX),dpi
 
      DOUBLE PRECISION daX1(2*NANG_MAX),daW1(2*NANG_MAX) 
      DOUBLE PRECISION x1,x2 
      INTEGER m,j,i,n 
      DOUBLE PRECISION z1,z,xm,xl,pp,p1,p2,p3,epss 
 
      epss = 3.0d-11 
      dpi  = pi * 1.0d0 
   
      x1 = -1.0D0 
      x2 = +1.0D0 
       
      IF ((nn .gt. NANG_MAX) .or. (nn .lt. 0)) THEN 
        write (*,*) 'need 0 < nn <= NANG_MAX'  
        Stop 
        END IF 
 
      n=nn*2 
 
      IF (MOD(n,2) .EQ. 1) THEN 
        write (*,*) 'need n to be even'  
        STOP 
        END IF 
 
      IF (x2 .LT. x1) THEN 
       xm = x1 
       x1 = x2 
       x2 = xm 
       END IF 
 
      m  = (n+1)/2 
 
      m  = n/2 

      m = n

      xm = 0.5*(x2+x1) 
      xl = 0.5*(x2-x1) 
 
      DO i = 1,m                    !loop over desired roots 

        z = cos(dpi*(i-0.25)/(n+0.5)) 
 20          CONTINUE 
        p1 = 1.0 
        p2 = 0.0 
        DO j = 1,n 
          p3 = p2 
          p2 = p1 
          p1 = ((2*j-1)*z*p2-(j-1)*p3)/j 
          END DO 
        pp = n*(z*p1-p2)/(z*z-1) 
        z1 = z 
        z  = z1-p1/pp 
        IF (abs(z-z1) .gt. epss) THEN 
          GOTO 20 
          END IF 
 
        daX(i)     = xm+xl*z 
        daW(i)     = 2*xl/((1-z*z)*pp*pp) 
 
        daX1(i)     = xm-xl*z 
        daX1(n+1-i) = xm+xl*z 
        daW1(i)     = 2*xl/((1-z*z)*pp*pp) 
        daW1(n+1-i) = daW(i) 
        END DO 
  
      DO i = 1,m 
        raX(i) = real(daX(i))
        raW(i) = real(daW(i))
        END DO

      RETURN 
      END 
 
