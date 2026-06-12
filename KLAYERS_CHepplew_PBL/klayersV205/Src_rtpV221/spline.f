C
       SUBROUTINE SPLINE(x,y,n,yp1,ypn,y2)
C
C      Calculate spline coefficients
C      from _Numerical Recipes in FORTRAN, second edition_
C
       INTEGER n,NMAX
       REAL yp1,ypn,x(n),y(n),y2(n)
       PARAMETER (NMAX=500)
       INTEGER i,k
       REAL p,qn,sig,un,u(NMAX)
C
C     ------------------------------------------------------------------
C
       IF (yp1 .GT. 0.99e30) THEN
          y2(1)=0.0
          u(1)=0.0
       ELSE
          y2(1)=-0.5
          u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       ENDIF
C
       DO 11 i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*y2(i-1)+2.0
          y2(i)=(sig-1.0)/p
          u(i)=(6.0*((y(i+1)-y(i))/(x(i+
     $       1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     $       u(i-1))/p
 11    CONTINUE
C
       IF (ypn .GT. 0.99E30) THEN
          qn=0.0
          un=0.0
       ELSE
         qn=0.5
         un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       ENDIF
C
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
C
       DO 12 k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
 12    CONTINUE
C
       RETURN
       END
