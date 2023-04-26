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
           FJACPRED1(1,7,L)=SECANG(L)*TRZ
           FJACPRED1(1,8,L)=SECANG(L)*TRZ/TR
 
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
 
           FJACPRED4(1,1,L)=SECANG(L)
           FJACPRED4(1,2,L)=SECANG(L)*SECANG(L)
           FJACPRED4(1,3,L)=SECANG(L)*TR
           FJACPRED4(1,4,L)=SECANG(L)*TJUNKS
           FJACPRED4(1,5,L)=TR
           FJACPRED4(1,6,L)=TJUNKS
           FJACPRED4(1,7,L)=SECANG(L)*TRZ
           FJACPRED4(1,8,L)=SECANG(L)*SECANG(L)*TRZ
           FJACPRED4(1,9,L)=SECANG(L)*SECANG(L)*TR
           FJACPRED4(1,10,L)=SECANG(L)*SECANG(L)*SECANG(L)
           FJACPRED4(1,11,L)=SQRT(SECANG(L))
 
C          Fixed predictors for FWO sun bfsw = set5
           FJACPRED5(1,1,L)=SECANG(L)
           FJACPRED5(1,2,L)=SECANG(L)*SECANG(L)
           FJACPRED5(1,3,L)=SECANG(L)*TR
           FJACPRED5(1,4,L)=SECANG(L)*TJUNKS
           FJACPRED5(1,5,L)=TR
           FJACPRED5(1,6,L)=TJUNKS
           FJACPRED5(1,7,L)=SECANG(L)*TRZ
           FJACPRED5(1,8,L)=SECANG(L)*TRZ/TR
           FJACPRED5(1,9,L)=SECANG(L)*SECANG(L)*TR
           FJACPRED5(1,10,L)=SQRT(SECANG(L))
           FJACPRED5(1,11,L)=TRZ
 
c          Fixed predictors for FWO sun mfmw = set6
           FJACPRED6(1,1,L)=SECANG(L)
           FJACPRED6(1,2,L)=SECANG(L)*SECANG(L)
           FJACPRED6(1,3,L)=SECANG(L)*TR
           FJACPRED6(1,4,L)=SECANG(L)*TJUNKS
           FJACPRED6(1,5,L)=TR
           FJACPRED6(1,6,L)=TJUNKS
           FJACPRED6(1,7,L)=SECANG(L)*TRZ
           FJACPRED6(1,8,L)=SQRT(SECANG(L))
 
C          Fixed predictors for FWO sun mfbw = set7
           FJACPRED7(1,1,L)=SECANG(L)
           FJACPRED7(1,2,L)=SECANG(L)*SECANG(L)
           FJACPRED7(1,3,L)=SECANG(L)*TR
           FJACPRED7(1,4,L)=SECANG(L)*TJUNKS
           FJACPRED7(1,5,L)=TR
           FJACPRED7(1,6,L)=TJUNKS
           FJACPRED7(1,7,L)=SECANG(L)*TRZ
           FJACPRED7(1,8,L)=SQRT(SECANG(L))
 
 
c          -----
C          Ozone
C          -----
           OJUNKA=SECANG(L)*A_O
           OJUNKR=SQRT( OJUNKA )
           OJUNKZ=OJUNKA/XZ_O
           OJUNKX=SECANG(L)*XZ_O

C          Ozone predictors for FWO = set1
           OJACPRED1(1,1,L)=OJUNKA
           OJACPRED1(1,2,L)=OJUNKR
           OJACPRED1(1,3,L)=OJUNKA*DT
           OJACPRED1(1,4,L)=OJUNKA*OJUNKA
           OJACPRED1(1,5,L)=OJUNKR*DT

C          ozone predictors for FOW = set2
           OJACPRED2(1, 1,L)=OJUNKA
           OJACPRED2(1, 2,L)=OJUNKR
           OJACPRED2(1, 3,L)=OJUNKA*DT
           OJACPRED2(1, 4,L)=OJUNKA*OJUNKA
           OJACPRED2(1, 5,L)=OJUNKR*DT
           OJACPRED2(1, 6,L)=OJUNKZ*A_O
           OJACPRED2(1, 7,L)=OJUNKR*A_O/XZ_O
           OJACPRED2(1, 8,L)=OJUNKZ*AZ_O
           OJACPRED2(1, 9,L)=OJUNKA*SQRT( OJUNKX )
           OJACPRED2(1,10,L)=OJUNKA*TAZ_O*SECANG(L)
C
C          There are no ozone predictors for set3 = FMW (the ozone
C          absorption in the region covered by FMW is negligible).
 
C          ozone predictors for FCOW = set4
           OJACPRED4(1,1,L)=OJUNKA
           OJACPRED4(1,2,L)=OJUNKR
           OJACPRED4(1,3,L)=OJUNKA*DT
C
C          ozone predictors for FWO sun bfsw = set5
           OJACPRED5(1,1,L)=OJUNKA
C
C          ozone predictors for FWO sun mfmw = set6
           OJACPRED6(1,1,L)=OJUNKA
C
C          ozone predictors for FWO sun mfbw = set7
           OJACPRED7(1,1,L)=OJUNKA
C
C
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
C
C
C          -----
C          Water
C          -----
           WJUNKA=SECANG(L)*A_W
           WJUNKR=SQRT( WJUNKA )
           WJUNKS=WJUNKA*WJUNKA
           WJUNKZ=WJUNKA*A_W/AZ_W
           WJUNK4=SQRT( WJUNKR )

C          Water predictors for FWO = set1
           WJACPRED1(1, 1,L)=WJUNKA
           WJACPRED1(1, 2,L)=WJUNKR
           WJACPRED1(1, 3,L)=WJUNKZ
           WJACPRED1(1, 4,L)=WJUNKA*DT
           WJACPRED1(1, 5,L)=WJUNKS
           WJACPRED1(1, 6,L)=WJUNKR*DT
           WJACPRED1(1, 7,L)=WJUNK4
           WJACPRED1(1, 8,L)=WJUNKZ/WJUNKR
           WJACPRED1(1, 9,L)=WJUNKS*WJUNKA
           WJACPRED1(1,10,L)=A_W
           WJACPRED1(1,11,L)=WJUNKA*DT*ABS( DT )

C          water predictors for FOW = set2
           WJACPRED2(1, 1,L)=WJUNKA
           WJACPRED2(1, 2,L)=WJUNKR
           WJACPRED2(1, 3,L)=WJUNKA*DT
           WJACPRED2(1, 4,L)=WJUNKA*OJUNKX
           WJACPRED2(1, 5,L)=WJUNKS
           WJACPRED2(1, 6,L)=WJUNK4
           WJACPRED2(1, 7,L)=WJUNKR*DT
           WJACPRED2(1, 8,L)=WJUNKZ
           WJACPRED2(1, 9,L)=WJUNKA*WJUNKS
           WJACPRED2(1,10,L)=WJUNKA*OJUNKX*OJUNKX
           WJACPRED2(1,11,L)=WJUNKZ/WJUNKR

C          water predictors for FMW = set3
           WJACPRED3(1, 1,L)=WJUNKA
           WJACPRED3(1, 2,L)=WJUNKR
           WJACPRED3(1, 3,L)=WJUNKZ
           WJACPRED3(1, 4,L)=WJUNKA*DT
           WJACPRED3(1, 5,L)=WJUNKS
           WJACPRED3(1, 6,L)=WJUNKR*DT
           WJACPRED3(1, 7,L)=WJUNK4
           WJACPRED3(1, 8,L)=WJUNKS*WJUNKA
           WJACPRED3(1, 9,L)=A_W
           WJACPRED3(1,10,L)=WJUNKZ/WJUNKR
           WJACPRED3(1,11,L)=WJUNKR*MJUNKZ

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

C          ---------------
C          Water continuum (for FWO, FOW, FMW, FCOW)
C          ---------------
           CONJACPRD(1,1,L)=WJUNKA/TJUNKS
           CONJACPRD(1,2,L)=CONJACPRD(1,1,L)*A_W/TJUNKS
           CONJACPRD(1,3,L)=WJUNKA/TR
           CONJACPRD(1,4,L)=CONJACPRD(1,3,L)*A_W
           CONJACPRD(1,5,L)=CONJACPRD(1,1,L)*A_W
           CONJACPRD(1,6,L)=CONJACPRD(1,1,L)/TJUNKS
           CONJACPRD(1,7,L)=WJUNKA
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
