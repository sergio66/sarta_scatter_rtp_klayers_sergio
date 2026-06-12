C
       SUBROUTINE SPLLEV(NIN, LNPIN, PIN, TIN, NGASES, MRIN,
     $    NSIN, LNPSIN, PSIN, TSIN, MRSIN, PMIN, PMAX)
C
C      Spline interpolation of input profile levels.
C
      include 'incLAY.f'
C
C      Parameters
       INTEGER NIN, NSIN,NGASES
       REAL LNPIN(MXIN), LNPSIN(MXIN), PSIN(MXIN), PIN(MXIN),
     $    TIN(MXIN), MRIN(MXIN,MXGAS), TSIN(MXIN), MRSIN(MXIN, MXGAS),
     $    PMIN, PMAX
C
C      Local variables
       INTEGER I,J,K,L,II, NLAY,NEWLL(MXIN),NEWLEV, NLEV,N,NOLD
       REAL COEF(5),X,DX,DLNP,NEWX(MXIN),OLDX(MXIN),YP1,YPN,
     $    OLDY(MXIN),NEWT(MXIN),TEMPMR(MXIN,MXGAS),XI,YI,Y2(MXIN)
C
C      -----------------------------------------------------------------
C
C      Layer thickness coefficients
       DATA (COEF(I),I=1,5)/
     $ -1.0999e-04, 1.9922e-03, -1.6689e-02, 9.6413e-02, -3.0577e-01/
C
C      -----------------------------------------------------------------
C
C      ----------------------------------
C      Loop over the input profile layers
C      ----------------------------------
       NLAY=0
       NSIN=0
       DO L=1,NIN-1
C
          NEWLL(L)=0
C
C         Check to see if the current level is in the desired range
          IF ( (PIN(L+1) .LT. PMAX) .OR. (PIN(L) .GT. PMIN)) THEN
C
C            Input profile layer thickness
             DLNP=LNPIN(L)-LNPIN(L+1)
C            Thickness for comparison
             X=-(COEF(5) + LNPIN(L)*( COEF(4) + LNPIN(L)*( COEF(3)
     $          + LNPIN(L)*( COEF(2) + LNPIN(L)*COEF(1) ) ) ))
C
             IF (DLNP .GT. X) THEN
C            --------------------------
C            Current layer is too thick
C            --------------------------
C               Add more levels
                NLAY=NLAY+1
                N=INT(DLNP/X + 0.5)
                NEWLL(L)=N
                DX=DLNP/(N+1)
                DO I=1,N
                   NEWLEV=NEWLEV+1
                   NEWX(NEWLEV)=LNPIN(L) - I*DX
                ENDDO
C
                IF (L. EQ. NIN-1 .AND. NLAY .GT. 3) THEN
C               -------------------------------------------------
C               If on last layer and NLAY>3, do spline interp now
C               -------------------------------------------------
C
C                  Do the mixing ratio spline interps
                   NOLD=NLAY+1
                   DO J=1,NGASES
                      DO I=1,NLAY
                         K=L - NLAY + I
                         II=NOLD + 1 - I
                         OLDX(II)=LNPIN(K)
                         OLDY(II)=MRIN(K,J)
                      ENDDO
                      OLDX(1)=LNPIN(L)
                      OLDY(1)=MRIN(L,J)
C
C                     Calc approx derivatives at 1st and last points
                      YP1=(OLDY(2)-OLDY(1))/(OLDX(2)-OLDX(1))
                      YPN=(OLDY(NOLD)-OLDY(NOLD-1))/
     $                   (OLDX(NOLD)-OLDX(NOLD-1))
C
                      CALL SPLINE(OLDX,OLDY,NOLD,YP1,YPN,Y2)
                      DO I=1,NEWLEV
                         XI=NEWX(I)
                         CALL SPLINT(OLDX,OLDY,Y2,NOLD,XI,YI)
                         TEMPMR(I,J)=YI
                      ENDDO
C
                   ENDDO
C
C                  Do the temperature spline interp
                   DO I=1,NLAY
                      K=L - NLAY + I
                      II=NOLD + 1 - I
                      OLDY(II)=TIN(K)
                   ENDDO
                   OLDY(1)=TIN(L)
C
C                  Calc approx derivatives at 1st and last points
                   YP1=(OLDY(2)-OLDY(1))/(OLDX(2)-OLDX(1))
                   YPN=(OLDY(NOLD)-OLDY(NOLD-1))/
     $                (OLDX(NOLD)-OLDX(NOLD-1))
C
                   CALL SPLINE(OLDX,OLDY,NOLD,YP1,YPN,Y2)
                   DO I=1,NEWLEV
                      XI=NEWX(I)
                      CALL SPLINT(OLDX,OLDY,Y2,NOLD,XI,YI)
                      NEWT(I)=YI
                   ENDDO
C
C                  ----------------------------
C                  Load up the output variables
C                  ----------------------------
                   NLEV=0
                   DO I=1,NLAY
                      K=L - NLAY + I
                      NSIN=NSIN+1
                      LNPSIN(NSIN)=LNPIN(K)
                      PSIN(NSIN)=PIN(K)
                      TSIN(NSIN)=TIN(K)
                      DO J=1,NGASES
                         MRSIN(NSIN,J)=MRIN(K,J)
                      ENDDO
C
                      DO II=1,NEWLL(K)
                         NSIN=NSIN+1
                         NLEV=NLEV+1
                         LNPSIN(NSIN)=NEWX(NLEV)
                         PSIN(NSIN)=EXP( NEWX(NLEV) )
                         TSIN(NSIN)=NEWT(NLEV)
                         DO J=1,NGASES
                            MRSIN(NSIN,J)=TEMPMR(NLEV,J)
                         ENDDO
                      ENDDO
                   ENDDO
C
                ENDIF
C
             ELSE
C            ------------------------------
C            Current layer is not too thick
C            ------------------------------
C
                IF (NLAY .GT. 3) THEN
C               ---------------------------------
C               Previous 3+ layers were too thick
C               ---------------------------------
C
C                  Do the mixing ratio spline interps
                   NOLD=NLAY+1
                   DO J=1,NGASES
                      DO I=1,NLAY
                         K=(L-1) - NLAY + I
                         II=NOLD + 1 - I
                         OLDX(II)=LNPIN(K)
                         OLDY(II)=MRIN(K,J)
                      ENDDO
                      OLDX(1)=LNPIN(L)
                      OLDY(1)=MRIN(L,J)
C
C                     Calc approx derivatives at 1st and last points
                      YP1=(OLDY(2)-OLDY(1))/(OLDX(2)-OLDX(1))
                      YPN=(OLDY(NOLD)-OLDY(NOLD-1))/
     $                   (OLDX(NOLD)-OLDX(NOLD-1))
C
                      CALL SPLINE(OLDX,OLDY,NOLD,YP1,YPN,Y2)
                      DO I=1,NEWLEV
                         XI=NEWX(I)
                         CALL SPLINT(OLDX,OLDY,Y2,NOLD,XI,YI)
                         TEMPMR(I,J)=YI
                      ENDDO
C
                   ENDDO
C
C                  Do the temperature spline interp
                   DO I=1,NLAY
                      K=(L-1) - NLAY + I
                      II=NOLD + 1 - I
                      OLDY(II)=TIN(K)
                   ENDDO
                   OLDY(1)=TIN(L)
C
C                  Calc approx derivatives at 1st and last points
                   YP1=(OLDY(2)-OLDY(1))/(OLDX(2)-OLDX(1))
                   YPN=(OLDY(NOLD)-OLDY(NOLD-1))/
     $                (OLDX(NOLD)-OLDX(NOLD-1))
C
                   CALL SPLINE(OLDX,OLDY,NOLD,YP1,YPN,Y2)
                   DO I=1,NEWLEV
                      XI=NEWX(I)
                      CALL SPLINT(OLDX,OLDY,Y2,NOLD,XI,YI)
                      NEWT(I)=YI
                   ENDDO
C
C                  ----------------------------
C                  Load up the output variables
C                  ----------------------------
                   NLEV=0
                   DO I=1,NLAY
                      K=(L-1) - NLAY + I
                      NSIN=NSIN+1
                      LNPSIN(NSIN)=LNPIN(K)
                      PSIN(NSIN)=PIN(K)
                      TSIN(NSIN)=TIN(K)
                      DO J=1,NGASES
                         MRSIN(NSIN,J)=MRIN(K,J)
                      ENDDO
C
                      DO II=1,NEWLL(K)
                         NSIN=NSIN+1
                         NLEV=NLEV+1
                         LNPSIN(NSIN)=NEWX(NLEV)
                         PSIN(NSIN)=EXP( NEWX(NLEV) )
                         TSIN(NSIN)=NEWT(NLEV)
                         DO J=1,NGASES
                            MRSIN(NSIN,J)=TEMPMR(NLEV,J)
                         ENDDO
                      ENDDO
                   ENDDO
C
                ELSEIF (NLAY .GT. 0) THEN
C               ----------------------------------------------
C               Update output variables of any previous layers
C               ----------------------------------------------
                   DO I=1,NLAY
                      K=(L-1) - NLAY + I
                      NSIN=NSIN+1
                      LNPSIN(NSIN)=LNPIN(K)
                      PSIN(NSIN)=PIN(K)
                      TSIN(NSIN)=TIN(K)
                      DO J=1,NGASES
                         MRSIN(NSIN,J)=MRIN(K,J)
                      ENDDO
                   ENDDO
                ENDIF
C
C               -----------------------------------------
C               Update output variables for current layer
C               -----------------------------------------
                NSIN=NSIN+1
                LNPSIN(NSIN)=LNPIN(L)
                PSIN(NSIN)=PIN(L)
                TSIN(NSIN)=TIN(L)
                DO J=1,NGASES
                   MRSIN(NSIN,J)=MRIN(L,J)
                ENDDO
C
                NLAY=0
                NEWLEV=0
             ENDIF
C
          ELSE
C         Currrent level is not in desired range
C
C            -----------------------------------------
C            Update output variables for current layer
C            -----------------------------------------
             NSIN=NSIN+1
             LNPSIN(NSIN)=LNPIN(L)
             PSIN(NSIN)=PIN(L)
             TSIN(NSIN)=TIN(L)
             DO J=1,NGASES
                MRSIN(NSIN,J)=MRIN(L,J)
             ENDDO
C
             NLAY=0
             NEWLEV=0
          ENDIF
C
       ENDDO
C
C      -------------------------------------
C      Update last level of output variables
C      -------------------------------------
       L=NIN
       NSIN=NSIN+1
       LNPSIN(NSIN)=LNPIN(L)
       PSIN(NSIN)=PIN(L)
       TSIN(NSIN)=TIN(L)
       DO J=1,NGASES
          MRSIN(NSIN,J)=MRIN(L,J)
       ENDDO
C
       RETURN
       END
