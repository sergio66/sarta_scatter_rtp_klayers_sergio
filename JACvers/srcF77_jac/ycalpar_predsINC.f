C
C         HDO terms
C         use water terms - but see below
C
C         ----------------------
C         Load up the predictors
C         ----------------------
C
C         -----
C         Fixed (for FWO, FOW, FMW, & FCOW)
C         -----

          TJUNKS=TR*TR
          FIXMUL(L)=A_F

          FPRED1(1,L)=SECANG(L)
          FPRED1(2,L)=SECANG(L)*SECANG(L)
          FPRED1(3,L)=SECANG(L)*TR
          FPRED1(4,L)=SECANG(L)*TJUNKS
          FPRED1(5,L)=TR
          FPRED1(6,L)=TJUNKS
          FPRED1(7,L)=SECANG(L)*TRZ
          FPRED1(8,L)=SECANG(L)*TRZ/TR
C
          FPRED2(1,L)=FPRED1(1,L)
          FPRED2(2,L)=FPRED1(2,L)
          FPRED2(3,L)=FPRED1(3,L)
          FPRED2(4,L)=FPRED1(4,L)
          FPRED2(5,L)=FPRED1(5,L)
          FPRED2(6,L)=FPRED1(6,L)
          FPRED2(7,L)=FPRED1(7,L)
          FPRED2(8,L)=FPRED1(8,L)
C
          FPRED3(1,L)=FPRED1(1,L)
          FPRED3(2,L)=FPRED1(2,L)
          FPRED3(3,L)=FPRED1(3,L)
          FPRED3(4,L)=FPRED1(4,L)
          FPRED3(5,L)=FPRED1(5,L)
          FPRED3(6,L)=FPRED1(6,L)
          FPRED3(7,L)=FPRED1(7,L)
          FPRED3(8,L)=FPRED1(8,L)
C
          FPRED4(1,L)=SECANG(L)
          FPRED4(2,L)=SECANG(L)*SECANG(L)
          FPRED4(3,L)=SECANG(L)*TR
          FPRED4(4,L)=SECANG(L)*TJUNKS
          FPRED4(5,L)=TR
          FPRED4(6,L)=TJUNKS
          FPRED4(7,L)=SECANG(L)*TRZ
          FPRED4(8,L)=SECANG(L)*SECANG(L)*TRZ
          FPRED4(9,L)=SECANG(L)*SECANG(L)*TR
          FPRED4(10,L)=SECANG(L)*SECANG(L)*SECANG(L)
          FPRED4(11,L)=SQRT(SECANG(L))
C
C         Fixed predictors for FWO sun bfsw = set5
          FPRED5(1,L)=SECANG(L)
          FPRED5(2,L)=SECANG(L)*SECANG(L)
          FPRED5(3,L)=SECANG(L)*TR
          FPRED5(4,L)=SECANG(L)*TJUNKS
          FPRED5(5,L)=TR
          FPRED5(6,L)=TJUNKS
          FPRED5(7,L)=SECANG(L)*TRZ
          FPRED5(8,L)=SECANG(L)*TRZ/TR
          FPRED5(9,L)=SECANG(L)*SECANG(L)*TR
          FPRED5(10,L)=SQRT(SECANG(L))
          FPRED5(11,L)=TRZ
C
C         Fixed predictors for FWO sun mfmw = set6
          FPRED6(1,L)=SECANG(L)
          FPRED6(2,L)=SECANG(L)*SECANG(L)
          FPRED6(3,L)=SECANG(L)*TR
          FPRED6(4,L)=SECANG(L)*TJUNKS
          FPRED6(5,L)=TR
          FPRED6(6,L)=TJUNKS
          FPRED6(7,L)=SECANG(L)*TRZ
          FPRED6(8,L)=SQRT(SECANG(L))
C
C         Fixed predictors for FWO sun mfbw = set7
          FPRED7(1,L)=SECANG(L)
          FPRED7(2,L)=SECANG(L)*SECANG(L)
          FPRED7(3,L)=SECANG(L)*TR
          FPRED7(4,L)=SECANG(L)*TJUNKS
          FPRED7(5,L)=TR
          FPRED7(6,L)=TJUNKS
          FPRED7(7,L)=SECANG(L)*TRZ
          FPRED7(8,L)=SQRT(SECANG(L))
C
C
C         -----
C         Ozone
C         -----
          OJUNKA=SECANG(L)*A_O
          OJUNKR=SQRT( OJUNKA )
          OJUNKZ=OJUNKA/XZ_O
          OJUNKX=SECANG(L)*XZ_O
C
C         Ozone predictors for FWO = set1
          OPRED1(1,L)=OJUNKA
          OPRED1(2,L)=OJUNKR
          OPRED1(3,L)=OJUNKA*DT
          OPRED1(4,L)=OJUNKA*OJUNKA
          OPRED1(5,L)=OJUNKR*DT
C
C         ozone predictors for FOW = set2
          OPRED2( 1,L)=OJUNKA
          OPRED2( 2,L)=OJUNKR
          OPRED2( 3,L)=OJUNKA*DT
          OPRED2( 4,L)=OJUNKA*OJUNKA
          OPRED2( 5,L)=OJUNKR*DT
          OPRED2( 6,L)=OJUNKZ*A_O
          OPRED2( 7,L)=OJUNKR*A_O/XZ_O
          OPRED2( 8,L)=OJUNKZ*AZ_O
          OPRED2( 9,L)=OJUNKA*SQRT( OJUNKX )
          OPRED2(10,L)=OJUNKA*TAZ_O*SECANG(L)
C
C         There are no ozone predictors for set3 = FMW (the ozone
C         absorption in the region covered by FMW is negligible).

C         ozone predictors for FCOW = set4
          OPRED4(1,L)=OJUNKA
          OPRED4(2,L)=OJUNKR
          OPRED4(3,L)=OJUNKA*DT
C
C         ozone predictors for FWO sun bfsw = set5
          OPRED5(1,L)=OJUNKA
C
C         ozone predictors for FWO sun mfmw = set6
          OPRED6(1,L)=OJUNKA
C
C         ozone predictors for FWO sun mfbw = set7
          OPRED7(1,L)=OJUNKA
C
C
C         -------
C         Methane for FMW = set3
C         -------
          MJUNKA=SECANG(L)*A_M
          MJUNKR=SQRT(MJUNKA)
          MJUNKZ=SECANG(L)*AZ_M
          MPRED3(1,L)=MJUNKA
          MPRED3(2,L)=MJUNKR
          MPRED3(3,L)=MJUNKA*DT
          MPRED3(4,L)=MJUNKA*MJUNKA
          MPRED3(5,L)=MJUNKA*SECANG(L)
          MPRED3(6,L)=MJUNKZ
          MPRED3(7,L)=A_M*DT
          MPRED3(8,L)=TAZ_M*SECANG(L)
          MPRED3(9,L)=SQRT( MJUNKZ )
C
C
C         -----
C         Water
C         -----
          WJUNKA=SECANG(L)*A_W
          WJUNKR=SQRT( WJUNKA )
          WJUNKS=WJUNKA*WJUNKA
          WJUNKZ=WJUNKA*A_W/AZ_W
          WJUNK4=SQRT( WJUNKR )
C
C         Water predictors for FWO = set1
          WPRED1( 1,L)=WJUNKA
          WPRED1( 2,L)=WJUNKR
          WPRED1( 3,L)=WJUNKZ
          WPRED1( 4,L)=WJUNKA*DT
          WPRED1( 5,L)=WJUNKS
          WPRED1( 6,L)=WJUNKR*DT
          WPRED1( 7,L)=WJUNK4
          WPRED1( 8,L)=WJUNKZ/WJUNKR
          WPRED1( 9,L)=WJUNKS*WJUNKA
          WPRED1(10,L)=A_W
          WPRED1(11,L)=WJUNKA*DT*ABS( DT )
C
C         water predictors for FOW = set2
          WPRED2( 1,L)=WJUNKA
          WPRED2( 2,L)=WJUNKR
          WPRED2( 3,L)=WJUNKA*DT
          WPRED2( 4,L)=WJUNKA*OJUNKX
          WPRED2( 5,L)=WJUNKS
          WPRED2( 6,L)=WJUNK4
          WPRED2( 7,L)=WJUNKR*DT
          WPRED2( 8,L)=WJUNKZ
          WPRED2( 9,L)=WJUNKA*WJUNKS
          WPRED2(10,L)=WJUNKA*OJUNKX*OJUNKX
          WPRED2(11,L)=WJUNKZ/WJUNKR
C
C         water predictors for FMW = set3
          WPRED3( 1,L)=WJUNKA
          WPRED3( 2,L)=WJUNKR
          WPRED3( 3,L)=WJUNKZ
          WPRED3( 4,L)=WJUNKA*DT
          WPRED3( 5,L)=WJUNKS
          WPRED3( 6,L)=WJUNKR*DT
          WPRED3( 7,L)=WJUNK4
          WPRED3( 8,L)=WJUNKS*WJUNKA
          WPRED3( 9,L)=A_W
          WPRED3(10,L)=WJUNKZ/WJUNKR
          WPRED3(11,L)=WJUNKR*MJUNKZ
C
C         water predictors for FCOW = set4
          WPRED4( 1,L)=WJUNKA
          WPRED4( 2,L)=A_W
          WPRED4( 3,L)=WJUNKR
          WPRED4( 4,L)=WJUNKA*DT
          WPRED4( 5,L)=WJUNKS
          WPRED4( 6,L)=WJUNKR*DT
          WPRED4( 7,L)=WJUNK4
          WPRED4( 8,L)=WJUNKZ
          WPRED4( 9,L)=WJUNKA*SECANG(L)
          WPRED4(10,L)=WJUNKS*WJUNKA
          WPRED4(11,L)=WJUNKA*AZ_C*SECANG(L)
          WPRED4(12,L)=WJUNKZ/WJUNKR
          WPRED4(13,L)=WJUNKA*DT*SECANG(L)
C
C         Water predictors for FWO sun bfsw = set5
          WPRED5( 1,L)=WJUNKA
          WPRED5( 2,L)=WJUNKA*WJUNKR
          WPRED5( 3,L)=WJUNKA*DT
C
C         Water predictors for FWO sun mfmw = set6
          WPRED6( 1,L)=WJUNKA
          WPRED6( 2,L)=WJUNKA*WJUNKR
          WPRED6( 3,L)=WJUNKA*DT
          WPRED6( 4,L)=WJUNKS
          WPRED6( 5,L)=WJUNKA*WJUNKR*DT
          WPRED6( 6,L)=WJUNKA*WJUNKS
          WPRED6( 7,L)=WJUNKA*SECANG(L)
C
C         Water predictors for FWO sun mfbw = set7
          WPRED7( 1,L)=WJUNKA
          WPRED7( 2,L)=WJUNKA*WJUNKR
          WPRED7( 3,L)=WJUNKA*DT
          WPRED7( 4,L)=WJUNKS
          WPRED7( 5,L)=WJUNKA*WJUNKR*DT
          WPRED7( 6,L)=WJUNKA*WJUNKS
          WPRED7( 7,L)=WJUNKA*SECANG(L)
          WPRED7( 8,L)=WJUNKZ
          WPRED7( 9,L)=WJUNKZ*WJUNKR
          WPRED7(10,L)=WJUNKA*WJUNK4
          WPRED7(11,L)=WJUNKA*WJUNKZ
          WPRED7(12,L)=WJUNKA*A_W
          WPRED7(13,L)=WJUNKS/WJUNK4
C
C         ---------------
C         Water continuum (for FWO, FOW, FMW, FCOW)
C         ---------------
          CONPRD(1,L)=WJUNKA/TJUNKS
          CONPRD(2,L)=CONPRD(1,L)*A_W/TJUNKS
          CONPRD(3,L)=WJUNKA/TR
          CONPRD(4,L)=CONPRD(3,L)*A_W
          CONPRD(5,L)=CONPRD(1,L)*A_W
          CONPRD(6,L)=CONPRD(1,L)/TJUNKS
          CONPRD(7,L)=WJUNKA
C          print *,'CALPAR CONPRD(1,L) = ',L,WJUNKA,TJUNKS,CONPRD(1,L)
C
C         ---------------
C         HDO
C         ---------------
        if (DEBUG) then
          IF(L .EQ. 96) write(6,'(A,X,I4,X,F6.2)') 'calpar: L,HDOFCT ',L,HDOFCT
        endif
          DJUNKA=SECANG(L)*A_W*(1 - HDOFCT)      ! *(1 - HDOFCT)
          DJUNKR=SQRT( DJUNKA )
          DJUNKS=DJUNKA*DJUNKA
          DJUNKZ=DJUNKA*A_W/AZ_W                 ! *(1 - HDOFCT)
          DJUNK4=SQRT( DJUNKR )
C
          DPRED( 1,L)=DJUNKA
          DPRED( 2,L)=DJUNKR
          DPRED( 3,L)=DJUNKZ
          DPRED( 4,L)=DJUNKA*DT
          DPRED( 5,L)=DJUNKS
          DPRED( 6,L)=DJUNKR*DT
          DPRED( 7,L)=DJUNK4
          DPRED( 8,L)=DJUNKZ/DJUNKR
          DPRED( 9,L)=DJUNKS*DJUNKA
          DPRED(10,L)=A_W
          DPRED(11,L)=DJUNKA*DT*ABS( DT )
C
C         ---------------
C         Carbon monoxide for FCOW = set4
C         ---------------
          CJUNKA=SECANG(L)*A_C
          CJUNKR=SQRT( CJUNKA )
          CJUNKS=CJUNKA*CJUNKA
          CJUNKZ=CJUNKA*A_C/AZ_C
          CPRED4(1,L)=CJUNKA
          CPRED4(2,L)=CJUNKR
          CPRED4(3,L)=CJUNKA*DT
          CPRED4(4,L)=CJUNKS
          CPRED4(5,L)=CJUNKZ
          CPRED4(6,L)=CJUNKR*DT
          CPRED4(7,L)=SQRT( CJUNKR )
          CPRED4(8,L)=CJUNKZ/CJUNKR
          CPRED4(9,L)=A_C
          CPRED4(10,L)=CJUNKA*SECANG(L)
          CPRED4(11,L)=CJUNKR*SECANG(L)
C
C         ---------------
C         trace gas perturbation predictors
C         ---------------
C         The first 4 trace predictors are used by all trace gases
          TRCPRD(1,L)=SECANG(L)
          TRCPRD(2,L)=TR
          TRCPRD(3,L)=SECANG(L)*TR
          TRCPRD(4,L)=SECANG(L)*TJUNKS
C         The last 3 trace predictors are only used by N2O
          TRCPRD(5,L)=SECANG(L)*SECANG(L)
          TRCPRD(6,L)=1.0
          TRCPRD(7,L)=SQRT( SECANG(L) )
C
