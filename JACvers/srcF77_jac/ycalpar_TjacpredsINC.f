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
 
          TJUNKS_T = 2*TR*TR_T
          TJUNKS_1 = 2*TR*TR_1
          TJUNKS_3 = 2*TR*TR_3

          FIXMUL_T(L)=A_F_T
          FIXMUL_1(L)=A_F_1
          FIXMUL_3(L)=A_F_3
 
           FJACPRED1(1,1,L)=0
           FJACPRED1(1,2,L)=0
           FJACPRED1(1,3,L)=SECANG(L)*TR_T
           FJACPRED1(1,4,L)=SECANG(L)*TJUNKS_T
           FJACPRED1(1,5,L)=TR_T
           FJACPRED1(1,6,L)=TJUNKS_T
           FJACPRED1(1,7,L)=SECANG(L)*TRZ_T
           FJACPRED1(1,8,L)=SECANG(L)*(TR*TRZ_T - TRZ*TR_T)/TR/TR
 
           FJACPRED2(1,1,L)=FJACPRED1(1,1,L)
           FJACPRED2(1,2,L)=FJACPRED1(1,2,L)
           FJACPRED2(1,3,L)=FJACPRED1(1,3,L)
           FJACPRED2(1,4,L)=FJACPRED1(1,4,L)
           FJACPRED2(1,5,L)=FJACPRED1(1,5,L)
           FJACPRED2(1,6,L)=FJACPRED1(1,6,L)
           FJACPRED2(1,7,L)=FJACPRED1(1,7,L)
           FJACPRED2(1,8,L)=FJACPRED1(1,8,L)
 
           FJACPRED3(1,1,L)=FJACPRED1(1,1,L)
           FJACPRED3(1,2,L)=FJACPRED1(1,2,L)
           FJACPRED3(1,3,L)=FJACPRED1(1,3,L)
           FJACPRED3(1,4,L)=FJACPRED1(1,4,L)
           FJACPRED3(1,5,L)=FJACPRED1(1,5,L)
           FJACPRED3(1,6,L)=FJACPRED1(1,6,L)
           FJACPRED3(1,7,L)=FJACPRED1(1,7,L)
           FJACPRED3(1,8,L)=FJACPRED1(1,8,L)
 
           FJACPRED4(1,1,L)=0.0
           FJACPRED4(1,2,L)=0.0
           FJACPRED4(1,3,L)=SECANG(L)*TR_T
           FJACPRED4(1,4,L)=SECANG(L)*TJUNKS_T
           FJACPRED4(1,5,L)=TR_T
           FJACPRED4(1,6,L)=TJUNKS_T
           FJACPRED4(1,7,L)=SECANG(L)*TRZ_T
           FJACPRED4(1,8,L)=SECANG(L)*SECANG(L)*TRZ_T
           FJACPRED4(1,9,L)=SECANG(L)*SECANG(L)*TR_T
           FJACPRED4(1,10,L)=0
           FJACPRED4(1,11,L)=0
 
C          Fixed predictors for FWO sun bfsw = set5
           FJACPRED5(1,1,L)=0.0
           FJACPRED5(1,2,L)=0.0
           FJACPRED5(1,3,L)=SECANG(L)*TR_T
           FJACPRED5(1,4,L)=SECANG(L)*TJUNKS_T
           FJACPRED5(1,5,L)=TR_T
           FJACPRED5(1,6,L)=TJUNKS_T
           FJACPRED5(1,7,L)=SECANG(L)*TRZ_T
           FJACPRED5(1,8,L)=SECANG(L)*(TR*TRZ_T - TRZ*TR_T)/TR/TR
           FJACPRED5(1,9,L)=SECANG(L)*SECANG(L)*TR_T
           FJACPRED5(1,10,L)=0.0
           FJACPRED5(1,11,L)=TRZ_T
 
c          Fixed predictors for FWO sun mfmw = set6
           FJACPRED6(1,1,L)=0
           FJACPRED6(1,2,L)=0
           FJACPRED6(1,3,L)=SECANG(L)*TR_T
           FJACPRED6(1,4,L)=SECANG(L)*TJUNKS_T
           FJACPRED6(1,5,L)=TR_T
           FJACPRED6(1,6,L)=TJUNKS_T
           FJACPRED6(1,7,L)=SECANG(L)*TRZ_T
           FJACPRED6(1,8,L)=0.0
 
C          Fixed predictors for FWO sun mfbw = set7
           FJACPRED7(1,1,L)=0
           FJACPRED7(1,2,L)=0
           FJACPRED7(1,3,L)=SECANG(L)*TR_T
           FJACPRED7(1,4,L)=SECANG(L)*TJUNKS_T
           FJACPRED7(1,5,L)=TR_T
           FJACPRED7(1,6,L)=TJUNKS_T
           FJACPRED7(1,7,L)=SECANG(L)*TRZ_T
           FJACPRED7(1,8,L)=0
 
c          -----
C          Ozone
C          -----
           OJUNKA_T=SECANG(L)*A_O_T                     !!! A_O_T = 0!!!
           OJUNKR_T=SQRT(SECANG(L))*0.5/SQRT(A_O)*A_O_T !!! A_O_T = 0!!!
           OJUNKZ_T=(XZ_O*OJUNKA_T - OJUNKA*XZ_O_T)/XZ_O/XZ_O
           OJUNKZ_T=(XZ_O*OJUNKA_T - OJUNKA*XZ_O_T)/XZ_O/XZ_O
           OJUNKX_T=SECANG(L)*XZ_O_T

C          Ozone predictors for FWO = set1
           OJACPRED1(1,1,L)=OJUNKA_T
           OJACPRED1(1,2,L)=OJUNKR_T
           OJACPRED1(1,3,L)=OJUNKA_T*DT + OJUNKA*DT_T
           OJACPRED1(1,4,L)=2*OJUNKA*OJUNKA_T
           OJACPRED1(1,5,L)=OJUNKR_T*DT + OJUNKR*DT_T

C          ozone predictors for FOW = set2
           OJACPRED2(1, 1,L)=OJUNKA_T
           OJACPRED2(1, 2,L)=OJUNKR_T
           OJACPRED2(1, 3,L)=OJUNKA_T*DT + OJUNKA*DT_T
           OJACPRED2(1, 4,L)=2*OJUNKA*OJUNKA_T
           OJACPRED2(1, 5,L)=OJUNKR*DT_T + OJUNKR_T*DT
           OJACPRED2(1, 6,L)=OJUNKZ*A_O_T + OJUNKZ_T*A_O
           OJACPRED2(1, 7,L)=OJUNKR_T*A_O/XZ_O + OJUNKR*(XZ_O*A_O_T - A_O*XZ_O_T)/XZ_O/XZ_O 
           OJACPRED2(1, 8,L)=(OJUNKZ*AZ_O_T + OJUNKZ_T*AZ_O)
           OJACPRED2(1, 9,L)=(OJUNKA*0.5/SQRT( OJUNKX )*OJUNKX_T + OJUNKA_T*SQRT(OJUNKX))
           OJACPRED2(1,10,L)=(OJUNKA*TAZ_O_T + OJUNKA_T*TAZ_O)*SECANG(L)

C          There are no ozone predictors for set3 = FMW (the ozone
C          absorption in the region covered by FMW is negligible).
 
C          ozone predictors for FCOW = set4
           OJACPRED4(1,1,L)=OJUNKA_T
           OJACPRED4(1,2,L)=OJUNKR_T
           OJACPRED4(1,3,L)=OJUNKA_T*DT + OJUNKA*DT_T
C
C          ozone predictors for FWO sun bfsw = set5
           OJACPRED5(1,1,L)=OJUNKA_T
C
C          ozone predictors for FWO sun mfmw = set6
           OJACPRED6(1,1,L)=OJUNKA_T
C
C          ozone predictors for FWO sun mfbw = set7
           OJACPRED7(1,1,L)=OJUNKA_T
C
C
c STOPPED HERE STOPPED HERE STOPPED HERE 
c STOPPED HERE STOPPED HERE STOPPED HERE 
c STOPPED HERE STOPPED HERE STOPPED HERE 
C          -------
C          Methane for FMW = set3
C          -------
           MJUNKA=SECANG(L)*A_M
           MJUNKR=SQRT(MJUNKA)
           MJUNKZ=SECANG(L)*AZ_M
           MJACPRED3(1,1,L)=MJUNKA
           MJACPRED3(1,2,L)=MJUNKR
           MJACPRED3(1,3,L)=MJUNKA*DT
           MJACPRED3(1,4,L)=MJUNKA*MJUNKA
           MJACPRED3(1,5,L)=MJUNKA*SECANG(L)
           MJACPRED3(1,6,L)=MJUNKZ
           MJACPRED3(1,7,L)=A_M*DT
           MJACPRED3(1,8,L)=TAZ_M*SECANG(L)
           MJACPRED3(1,9,L)=SQRT( MJUNKZ )
c STOPPED HERE STOPPED HERE STOPPED HERE 
c STOPPED HERE STOPPED HERE STOPPED HERE 
c STOPPED HERE STOPPED HERE STOPPED HERE 

C
C
C          -----
C          Water
C          -----
           WJUNKA_T=SECANG(L)*A_W_T
           WJUNKR_T=SQRT(SECANG(L))*0.5/SQRT(A_W)*A_W_T
           WJUNKS_T=2*WJUNKA*WJUNKA_T
           WJUNKZ_T=WJUNKA_T*A_W/AZ_W + WJUNKA*(AZ_W*A_W_T - A_W*AZ_W_T)/AZ_W/AZ_W
           WJUNK4_T=(SECANG(L)**0.25)*0.25/(A_W**(0.75))*A_W_T

C          Water predictors for FWO = set1
           WJACPRED1(1, 1,L)=WJUNKA_T
           WJACPRED1(1, 2,L)=WJUNKR_T
           WJACPRED1(1, 3,L)=WJUNKZ_T
           WJACPRED1(1, 4,L)=WJUNKA_T*DT + WJUNKA*DT_T
           WJACPRED1(1, 5,L)=WJUNKS_T
           WJACPRED1(1, 6,L)=WJUNKR_T*DT + WJUNKR*DT_T
           WJACPRED1(1, 7,L)=WJUNK4_T
           WJACPRED1(1, 8,L)=(WJUNKR*WJUNKZ_T - WJUNKZ*WJUNKR_T)/WJUNKR/WJUNKR
           WJACPRED1(1, 9,L)=WJUNKS*WJUNKA_T + WJUNKS_T*WJUNKA 
           WJACPRED1(1,10,L)=A_W_T
           WJACPRED1(1,11,L)=WJUNKA*2*DT*DT_T + WJUNKA_T*(DT**2)

C          water predictors for FOW = set2
           WJACPRED2(1, 1,L)=WJUNKA_T
           WJACPRED2(1, 2,L)=WJUNKR_T
           WJACPRED2(1, 3,L)=WJUNKA_T*DT + WJUNKA*DT_T
           WJACPRED2(1, 4,L)=WJUNKA_T*OJUNKX + WJUNKA*OJUNKX_T
           WJACPRED2(1, 5,L)=WJUNKS_T
           WJACPRED2(1, 6,L)=WJUNK4_T
           WJACPRED2(1, 7,L)=WJUNKR_T*DT + WJUNKR*DT_T
           WJACPRED2(1, 8,L)=WJUNKZ_T
           WJACPRED2(1, 9,L)=WJUNKA_T*WJUNKS + WJUNKA*WJUNKS_T
           WJACPRED2(1,10,L)=WJUNKA_T*OJUNKX*OJUNKX + WJUNKA*OJUNKX_T*OJUNKX + WJUNKA*OJUNKX*OJUNKX_T
           WJACPRED2(1,11,L)=(WJUNKR*WJUNKZ_T - WJUNKZ*WJUNKR_T)/ WJUNKR/WJUNKR

C          water predictors for FMW = set3
           WJACPRED3(1, 1,L)=WJUNKA_T
           WJACPRED3(1, 2,L)=WJUNKR_T
           WJACPRED3(1, 3,L)=WJUNKZ_T
           WJACPRED3(1, 4,L)=WJUNKA_T*DT + WJUNKA*DT_T
           WJACPRED3(1, 5,L)=WJUNKS_T
           WJACPRED3(1, 6,L)=WJUNKR_T*DT + WJUNKR*DT_T
           WJACPRED3(1, 7,L)=WJUNK4_T
           WJACPRED3(1, 8,L)=WJUNKS_T*WJUNKA + WJUNKS*WJUNKA_T
           WJACPRED3(1, 9,L)=A_W_T
           WJACPRED3(1,10,L)=(WJUNKR*WJUNKZ_T - WJUNKZ*WJUNKR_T)/WJUNKR/WJUNKR
           WJACPRED3(1,11,L)=WJUNKR_T*MJUNKZ + WJUNKR*MJUNKZ_T

c STOPPPED HERE STOPPED HERE XXXXXXX
C          water predictors for FCOW = set4
           WJACPRED4(1, 1,L)=WJUNKA
           WJACPRED4(1, 2,L)=A_W
           WJACPRED4(1, 3,L)=WJUNKR
           WJACPRED4(1, 4,L)=WJUNKA*DT
           WJACPRED4(1, 5,L)=WJUNKS
           WJACPRED4(1, 6,L)=WJUNKR*DT
           WJACPRED4(1, 7,L)=WJUNK4
           WJACPRED4(1, 8,L)=WJUNKZ
           WJACPRED4(1, 9,L)=WJUNKA*SECANG(L)
           WJACPRED4(1,10,L)=WJUNKS*WJUNKA
           WJACPRED4(1,11,L)=WJUNKA*AZ_C*SECANG(L)
           WJACPRED4(1,12,L)=WJUNKZ/WJUNKR
           WJACPRED4(1,13,L)=WJUNKA*DT*SECANG(L)

C          Water predictors for FWO sun bfsw = set5
           WJACPRED5(1, 1,L)=WJUNKA
           WJACPRED5(1, 2,L)=WJUNKA*WJUNKR
           WJACPRED5(1, 3,L)=WJUNKA*DT

C          Water predictors for FWO sun mfmw = set6
           WJACPRED6(1, 1,L)=WJUNKA
           WJACPRED6(1, 2,L)=WJUNKA*WJUNKR
           WJACPRED6(1, 3,L)=WJUNKA*DT
           WJACPRED6(1, 4,L)=WJUNKS
           WJACPRED6(1, 5,L)=WJUNKA*WJUNKR*DT
           WJACPRED6(1, 6,L)=WJUNKA*WJUNKS
           WJACPRED6(1, 7,L)=WJUNKA*SECANG(L)

C          Water predictors for FWO sun mfbw = set7
           WJACPRED7(1, 1,L)=WJUNKA
           WJACPRED7(1, 2,L)=WJUNKA*WJUNKR
           WJACPRED7(1, 3,L)=WJUNKA*DT
           WJACPRED7(1, 4,L)=WJUNKS
           WJACPRED7(1, 5,L)=WJUNKA*WJUNKR*DT
           WJACPRED7(1, 6,L)=WJUNKA*WJUNKS
           WJACPRED7(1, 7,L)=WJUNKA*SECANG(L)
           WJACPRED7(1, 8,L)=WJUNKZ
           WJACPRED7(1, 9,L)=WJUNKZ*WJUNKR
           WJACPRED7(1,10,L)=WJUNKA*WJUNK4
           WJACPRED7(1,11,L)=WJUNKA*WJUNKZ
           WJACPRED7(1,12,L)=WJUNKA*A_W
           WJACPRED7(1,13,L)=WJUNKS/WJUNK4
c STOPPPED HERE STOPPED HERE XXXXXXX

C          ---------------
C          Water continuum (for FWO, FOW, FMW, FCOW)
C          ---------------
           CONJACPRD(1,1,L)=(TJUNKS*WJUNKA_T - WJUNKA*TJUNKS_T)/TJUNKS/TJUNKS
           CONJACPRD(1,2,L)=CONJACPRD(1,1,L)*A_W/TJUNKS + CONPRD(1,L)*(TJUNKS*A_W_T - A_W*TJUNKS_T)/TJUNKS/TJUNKS
           CONJACPRD(1,3,L)=(TR*WJUNKA_T - WJUNKA*TR_T)/TR/TR
           CONJACPRD(1,4,L)=CONJACPRD(1,3,L)*A_W + CONPRD(3,L)*A_W_T
           CONJACPRD(1,5,L)=CONJACPRD(1,1,L)*A_W + CONPRD(1,L)*A_W_T
           CONJACPRD(1,6,L)=(TJUNKS*CONJACPRD(1,1,L) - CONPRD(1,L)*TJUNKS_T)/TJUNKS/TJUNKS
           CONJACPRD(1,7,L)=WJUNKA_T
C         print *,'CALPAR CONJACPRD(1,1,L) = ',L,WJUNKA,TJUNKS,CONJACPRD(1,1,L)

C          ---------------
C          HDO
C          ---------------
           if (DEBUG) then
             IF(L .EQ. 96) write(6,'(A,X,I4,X,F6.2)') 'calpar: L,HDOFCT ',L,HDOFCT
           endif
           DJUNKA=SECANG(L)*A_W*(1 - HDOFCT)      ! *(1 - HDOFCT)
           DJUNKR=SQRT( DJUNKA )
           DJUNKS=DJUNKA*DJUNKA
           DJUNKZ=DJUNKA*A_W/AZ_W                 ! *(1 - HDOFCT)
           DJUNK4=SQRT( DJUNKR )

           DJACPRED(1, 1,L)=DJUNKA
           DJACPRED(1, 2,L)=DJUNKR
           DJACPRED(1, 3,L)=DJUNKZ
           DJACPRED(1, 4,L)=DJUNKA*DT
           DJACPRED(1, 5,L)=DJUNKS
           DJACPRED(1, 6,L)=DJUNKR*DT
           DJACPRED(1, 7,L)=DJUNK4
           DJACPRED(1, 8,L)=DJUNKZ/DJUNKR
           DJACPRED(1, 9,L)=DJUNKS*DJUNKA
           DJACPRED(1,10,L)=A_W
           DJACPRED(1,11,L)=DJUNKA*DT*ABS( DT )

C          ---------------
C          Carbon monoxide for FCOW = set4
C          ---------------
           CJUNKA=SECANG(L)*A_C
           CJUNKR=SQRT( CJUNKA )
           CJUNKS=CJUNKA*CJUNKA
           CJUNKZ=CJUNKA*A_C/AZ_C
           CJACPRED4(1,1,L)=CJUNKA
           CJACPRED4(1,2,L)=CJUNKR
           CJACPRED4(1,3,L)=CJUNKA*DT
           CJACPRED4(1,4,L)=CJUNKS
           CJACPRED4(1,5,L)=CJUNKZ
           CJACPRED4(1,6,L)=CJUNKR*DT
           CJACPRED4(1,7,L)=SQRT( CJUNKR )
           CJACPRED4(1,8,L)=CJUNKZ/CJUNKR
           CJACPRED4(1,9,L)=A_C
           CJACPRED4(1,10,L)=CJUNKA*SECANG(L)
           CJACPRED4(1,11,L)=CJUNKR*SECANG(L)
C
C         ---------------
C         trace gas perturbation predictors
C         ---------------
C         The first 4 trace predictors are used by all trace gases
          TRCJACPRD(1,1,L)=SECANG(L)
          TRCJACPRD(1,2,L)=TR
          TRCJACPRD(1,3,L)=SECANG(L)*TR
          TRCJACPRD(1,4,L)=SECANG(L)*TJUNKS
C         The last 3 trace predictors are only used by N2O
          TRCJACPRD(1,5,L)=SECANG(L)*SECANG(L)
          TRCJACPRD(1,6,L)=1.0
          TRCJACPRD(1,7,L)=SQRT( SECANG(L) )
