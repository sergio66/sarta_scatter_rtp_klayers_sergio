C
C         HDO terms
C         use water terms - but see below
C
C         ----------------------
C         Load up the predictors
C         ----------------------
C
C          -----
C          Fixed (for FWO, FOW, FMW, & FCOW)
C          recall d(u/v) = (v du - u dv)/v^2
C          -----

        IF (INTERSECT(2,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
          !!! CO2
          IWHICHJAC = 4 
          FJACPRED1(IWHICHJAC,1:8,L) = 0
          FJACPRED2(IWHICHJAC,1:8,L) = 0
          FJACPRED3(IWHICHJAC,1:8,L) = 0
          FJACPRED4(IWHICHJAC,1:11,L) = 0
          FJACPRED5(IWHICHJAC,1:11,L) = 0
          FJACPRED6(IWHICHJAC,1:8,L) = 0
          FJACPRED7(IWHICHJAC,1:8,L) = 0
          DJACPRED(IWHICHJAC,1:11,L) = 0
          TRCJACPRD(IWHICHJAC,1:4,L) = 0
        END IF
        IF (INTERSECT(4,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
          !!! N2O
          IWHICHJAC = 5
          FJACPRED1(IWHICHJAC,1:8,L) = 0
          FJACPRED2(IWHICHJAC,1:8,L) = 0
          FJACPRED3(IWHICHJAC,1:8,L) = 0
          FJACPRED4(IWHICHJAC,1:11,L) = 0
          FJACPRED5(IWHICHJAC,1:11,L) = 0
          FJACPRED6(IWHICHJAC,1:8,L) = 0
          FJACPRED7(IWHICHJAC,1:8,L) = 0
          DJACPRED(IWHICHJAC,1:11,L) = 0
          TRCJACPRD(IWHICHJAC,1:4,L) = 0
        END IF
        IF (INTERSECT(5,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
          !!! CO
          IWHICHJAC = 6
          FJACPRED1(IWHICHJAC,1:8,L) = 0
          FJACPRED2(IWHICHJAC,1:8,L) = 0
          FJACPRED3(IWHICHJAC,1:8,L) = 0
          FJACPRED4(IWHICHJAC,1:11,L) = 0
          FJACPRED5(IWHICHJAC,1:11,L) = 0
          FJACPRED6(IWHICHJAC,1:8,L) = 0
          FJACPRED7(IWHICHJAC,1:8,L) = 0
          DJACPRED(IWHICHJAC,1:11,L) = 0
          TRCJACPRD(IWHICHJAC,1:4,L) = 0
        END IF
        IF (INTERSECT(6,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
          !!! CH4
          IWHICHJAC = 7
          FJACPRED1(IWHICHJAC,1:8,L) = 0
          FJACPRED2(IWHICHJAC,1:8,L) = 0
          FJACPRED3(IWHICHJAC,1:8,L) = 0
          FJACPRED4(IWHICHJAC,1:11,L) = 0
          FJACPRED5(IWHICHJAC,1:11,L) = 0
          FJACPRED6(IWHICHJAC,1:8,L) = 0
          FJACPRED7(IWHICHJAC,1:8,L) = 0
          DJACPRED(IWHICHJAC,1:11,L) = 0
          TRCJACPRD(IWHICHJAC,1:4,L) = 0
        END IF
        IF (INTERSECT(9,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
          !!! SO2
          IWHICHJAC = 8
          FJACPRED1(IWHICHJAC,1:8,L) = 0
          FJACPRED2(IWHICHJAC,1:8,L) = 0
          FJACPRED3(IWHICHJAC,1:8,L) = 0
          FJACPRED4(IWHICHJAC,1:11,L) = 0
          FJACPRED5(IWHICHJAC,1:11,L) = 0
          FJACPRED6(IWHICHJAC,1:8,L) = 0
          FJACPRED7(IWHICHJAC,1:8,L) = 0
          DJACPRED(IWHICHJAC,1:11,L) = 0
          TRCJACPRD(IWHICHJAC,1:4,L) = 0
        END IF
        IF (INTERSECT(12,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
          !!! HNO3
          IWHICHJAC = 9
          FJACPRED1(IWHICHJAC,1:8,L) = 0
          FJACPRED2(IWHICHJAC,1:8,L) = 0
          FJACPRED3(IWHICHJAC,1:8,L) = 0
          FJACPRED4(IWHICHJAC,1:11,L) = 0
          FJACPRED5(IWHICHJAC,1:11,L) = 0
          FJACPRED6(IWHICHJAC,1:8,L) = 0
          FJACPRED7(IWHICHJAC,1:8,L) = 0
          DJACPRED(IWHICHJAC,1:11,L) = 0
          TRCJACPRD(IWHICHJAC,1:4,L) = 0
        END IF

        IF (INTERSECT(6,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
          !!! CH4
          IWHICHJAC = 7
C          -------
C          Methane for FMW = set3
C          -------
           MJUNKA_6=SECANG(L)*A_M_6
           MJUNKR_6=0.5/SQRT(MJUNKA)*MJUNKA_6
           MJUNKZ_6=SECANG(L)*AZ_M_6
           MJACPRED3(IWHICHJAC,1,L)=MJUNKA_6
           MJACPRED3(IWHICHJAC,2,L)=MJUNKR_6
           MJACPRED3(IWHICHJAC,3,L)=MJUNKA_6*DT + MJUNKA*DT_6
           MJACPRED3(IWHICHJAC,4,L)=2*MJUNKA*MJUNKA_6
           MJACPRED3(IWHICHJAC,5,L)=MJUNKA_6*SECANG(L)
           MJACPRED3(IWHICHJAC,6,L)=MJUNKZ_6
           MJACPRED3(IWHICHJAC,7,L)=A_M_6*DT + A_M*DT_6
           MJACPRED3(IWHICHJAC,8,L)=TAZ_M_6*SECANG(L)
           MJACPRED3(IWHICHJAC,9,L)=0.5/SQRT( MJUNKZ )*MJUNKZ_6
         END IF

        IF (INTERSECT(5,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
          !!! CO
          IWHICHJAC = 6
C          ---------------
C          Carbon monoxide for FCOW = set4
C          ---------------
           CJUNKA=SECANG(L)*A_C_5
           CJUNKR=0.5/SQRT( CJUNKA ) * CJUNKA_5
           CJUNKS=2*CJUNKA*CJUNKA_5
           CJUNKZ=CJUNKA_5*A_C/AZ_C + CJUNKA*(AZ_C*A_C_5 + A_C*AZ_C_5)/AZ_C/AZ_C
           CJACPRED4(IWHICHJAC,1,L)=CJUNKA_5
           CJACPRED4(IWHICHJAC,2,L)=CJUNKR_5
           CJACPRED4(IWHICHJAC,3,L)=CJUNKA_5*DT + CJUNKA*DT_5
           CJACPRED4(IWHICHJAC,4,L)=CJUNKS_5
           CJACPRED4(IWHICHJAC,5,L)=CJUNKZ_5
           CJACPRED4(IWHICHJAC,6,L)=CJUNKR_5*DT + CJUNKR*DT_5
           CJACPRED4(IWHICHJAC,7,L)=0.5/SQRT( CJUNKR )*CJUNKR_5
           CJACPRED4(IWHICHJAC,8,L)=(CJUNKR*CJUNKZ_5 -  CJUNKZ*CJUNKR_5)/CJUNKR/CJUNKR
           CJACPRED4(IWHICHJAC,9,L)=A_C_5
           CJACPRED4(IWHICHJAC,10,L)=CJUNKA_5*SECANG(L)
           CJACPRED4(IWHICHJAC,11,L)=CJUNKR_5*SECANG(L)
         END IF

        IF (INTERSECT(4,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
          !!! N2O
          IWHICHJAC = 5
C         ---------------
C         trace gas perturbation predictors
C         ---------------
C         The last 3 trace predictors are only used by N2O
          TRCJACPRD(IWHICHJAC,5,L)=0
          TRCJACPRD(IWHICHJAC,6,L)=0
          TRCJACPRD(IWHICHJAC,7,L)=0
        END IF
