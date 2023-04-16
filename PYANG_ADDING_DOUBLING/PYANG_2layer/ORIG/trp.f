
      FUNCTION trp( A, !! c Input
                    N, !! c Input
                    X, !! c Input
                    Y ) ! c Input

        IMPLICIT NONE
       
c--------------
c  Arguments
c--------------

        REAL A
        INTEGER N
        REAL X(N)
        REAL Y(N)

c------------------
c Local variables
c------------------
        INTEGER  NM1
        INTEGER  IL
        INTEGER  IR
        INTEGER  I
        REAL     R

c------------
c Functions
c------------
        REAL trp

        NM1 = N-1
        IF ( A < X( 2 ) ) GO TO 50
        IF ( A >= X( NM1 ) ) GO TO 40
        IL = 2
        IR = NM1
c
c     BISECTION SEARCH
c
   10   I = ( IL + IR ) / 2
        IF ( I == IL ) GO TO 60
        IF ( A - X( I ) ) 20, 60, 30
   20   IR = I
        GO TO 10
   30   IL = I
        GO TO 10
c
c     A.LT.X(2) .OR. A.GE.X(N-1)
c
   40   I = NM1
        GO TO 60
   50   I = 1
c
c     EVALUATION
c
   60   R = ( A - X( I ) ) / ( X( I + 1 ) - X( I ) )
        TRP = Y( I ) + R * ( Y( I + 1 ) - Y( I ) )

      END FUNCTION trp
 
