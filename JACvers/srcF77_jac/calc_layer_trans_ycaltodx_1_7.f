      SUBROUTINE CALC_LAYER_TRANS_CALTODX_1_7(
     $ I,DOSUN, NCHAN, LSTCHN, FREQ, INDCHN, LBOT, 
     $ QUICKCLIST1, QUICKCLIST2, QUICKCLIST3, QUICKCLIST4, 
     $ QUICKCLIST5, QUICKCLIST6, QUICKCLIST7, 
     $ NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7, 
     $ CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, 
     $ FIXMUL, CONPRD, SUNCONPRD, TRCPRD, SUNTRCPRD, 
     $ DPRED, SUNDPRED, 
     $ INDCO2, COFCO2, INDSO2, COFSO2, INDHDO, COFHDO, 
     $ INDHNO, COFHNO, INDN2O, COFN2O, INDNH3, COFNH3,
     $ CO2MLT, SO2MLT, HNOMLT, N2OMLT, NH3MLT, HDOMLT, 
     $ INDH2O, H2OPRD, COFH2O,
     $ WAANG, LOPMIN, LOPMAX, LOPUSE, LOPLOW, DAOP, WAOP, 
     $ NFAKE, INDFAK, QUICKINDFAK, CO2TOP,
     $ COEF1, COEF2, COEF3, COEF4, COEF5, COEF6, COEF7, COEFF, 
     $ WPRED1, WPRED2, WPRED3, WPRED4, WPRED5, WPRED6, WPRED7, 
     $                     SUNWPRED4, SUNWPRED5, SUNWPRED6, SUNWPRED7, 
     $ OPRED1, OPRED2,         OPRED4, OPRED5, OPRED6, OPRED7, 
     $                     SUNOPRED4, SUNOPRED5, SUNOPRED6, SUNOPRED7, 
     $ FPRED1, FPRED2, FPRED3, FPRED4, FPRED5, FPRED6, FPRED7, 
     $                     SUNFPRED4, SUNFPRED5, SUNFPRED6, SUNFPRED7, 
     $ MPRED3, CPRED4, SUNCPRED4, SECANG, SECSUN, SUNFDG, SUNCOS,
     $ TAU, TAUZ, TAUZSN)

      IMPLICIT NONE
      INCLUDE "incFTC.f"

c output
       REAL    TAU(MAXLAY,MXCHAN) ! chan layer effective optical depth
       REAL   TAUZ(MAXLAY,MXCHAN) ! chan surface-to-space trans
       REAL TAUZSN(MAXLAY,MXCHAN) ! sun space-to-surface-to-space OD

c input
      INTEGER I               !! channel index loop
      LOGICAL DOSUN
      INTEGER NCHAN
      INTEGER LSTCHN(MXCHAN)  ! list of selected channels
      REAL   FREQ(MXCHAN)    ! chan center frequency

      REAL SECANG(MAXLAY)        ! local path angle secant
      REAL SECSUN(MAXLAY) ! secant of effective sun local path angle
      REAL SUNFDG         ! fudge factor for large solar angles
      REAL SUNCOS         ! cosine of sun zenith angle

      REAL MPRED3( N3CH4,MAXLAY) ! set3 methane predictors

      REAL CPRED4(  N4CO,MAXLAY) ! set4 carbon monoxide predictors
      REAL SUNCPRED4(  N4CO,MAXLAY) ! set4 carbon monoxide predictors
      REAL TRCPRD(NTRACE,MAXLAY) ! trace gas pert perdictors
      REAL SUNTRCPRD(NTRACE,MAXLAY) ! trace gas pert perdictors
      REAL FIXMUL(MAXLAY)        ! "fixed" amount multiplier (~1)
      REAL CONPRD( N1CON,MAXLAY) ! water continuum predictors
      REAL SUNCONPRD( N1CON,MAXLAY) ! water continuum predictors
      REAL  DPRED(  NHDO,MAXLAY) ! HDO perturbation predictors
      REAL  SUNDPRED(  NHDO,MAXLAY) ! HDO perturbation predictors

      REAL FPRED1( N1FIX,MAXLAY) ! set1 "fixed" predictors
      REAL FPRED2( N2FIX,MAXLAY) ! set2 "fixed" predictors
      REAL FPRED3( N3FIX,MAXLAY) ! set3 "fixed" predictors
      REAL FPRED4( N4FIX,MAXLAY) ! set4 "fixed" predictors
      REAL FPRED5( N5FIX,MAXLAY) ! set5 "fixed" predictors
      REAL FPRED6( N6FIX,MAXLAY) ! set6 "fixed" predictors
      REAL FPRED7( N7FIX,MAXLAY) ! set7 "fixed" predictors
      REAL SUNFPRED4( N4FIX,MAXLAY) ! set4 "fixed" predictors
      REAL SUNFPRED5( N5FIX,MAXLAY) ! set5 "fixed" predictors
      REAL SUNFPRED6( N6FIX,MAXLAY) ! set6 "fixed" predictors
      REAL SUNFPRED7( N7FIX,MAXLAY) ! set7 "fixed" predictors

      REAL OPRED1(  N1O3,MAXLAY) ! set1 ozone predictors
      REAL OPRED2(  N2O3,MAXLAY) ! set2 ozone predictors
      REAL OPRED4(  N4O3,MAXLAY) ! set4 ozone predictors
      REAL OPRED5(  N5O3,MAXLAY) ! set5 ozone predictors
      REAL OPRED6(  N6O3,MAXLAY) ! set6 ozone predictors
      REAL OPRED7(  N7O3,MAXLAY) ! set7 ozone predictors
      REAL SUNOPRED4( N4O3,MAXLAY) ! set4 "fixed" predictors
      REAL SUNOPRED5( N5O3,MAXLAY) ! set5 "fixed" predictors
      REAL SUNOPRED6( N6O3,MAXLAY) ! set6 "fixed" predictors
      REAL SUNOPRED7( N7O3,MAXLAY) ! set7 "fixed" predictors

      REAL WPRED1( N1H2O,MAXLAY) ! set1 water predictors
      REAL WPRED2( N2H2O,MAXLAY) ! set2 water predictors
      REAL WPRED3( N3H2O,MAXLAY) ! set3 water predictors
      REAL WPRED4( N4H2O,MAXLAY) ! set4 water predictors
      REAL WPRED5( N5H2O,MAXLAY) ! set5 water predictors
      REAL WPRED6( N6H2O,MAXLAY) ! set6 water predictors
      REAL WPRED7( N7H2O,MAXLAY) ! set7 water predictors
      REAL SUNWPRED4( N4H2O,MAXLAY) ! set4 "fixed" predictors
      REAL SUNWPRED5( N5H2O,MAXLAY) ! set5 "fixed" predictors
      REAL SUNWPRED6( N6H2O,MAXLAY) ! set6 "fixed" predictors
      REAL SUNWPRED7( N7H2O,MAXLAY) ! set7 "fixed" predictors

      REAL  COEF1(N1COEF,MAXLAY,MXCHN1) ! coefs for set1 chans
      REAL  COEF2(N2COEF,MAXLAY,MXCHN2) ! coefs for set2 chans
      REAL  COEF3(N3COEF,MAXLAY,MXCHN3) ! coefs for set3 chans
      REAL  COEF4(N4COEF,MAXLAY,MXCHN4) ! coefs for set4 chans
      REAL  COEF5(N5COEF,MAXLAY,MXCHN5) ! coefs for set5 chans
      REAL  COEF6(N6COEF,MAXLAY,MXCHN6) ! coefs for set6 chans
      REAL  COEF7(N7COEF,MAXLAY,MXCHN7) ! coefs for set7 chans
      REAL  COEFF(NFCOEF,MXCHAN)        ! coefs for chan "F" factor

      INTEGER  NFAKE              ! # of channels to "fake"
      INTEGER INDFAK(MXCHAN)      ! indices of channels to fake
      INTEGER QUICKINDFAK(MXCHAN) ! list of set1 channels

C     for CALOWP
      REAL  WAANG(MAXLAY)
      INTEGER LOPMIN
      INTEGER LOPMAX
      LOGICAL LOPUSE(MXOWLY)
      INTEGER LOPLOW(MAXLAY)

      REAL H2OPRD(  NH2O,MXOWLY)
      REAL  DAOP(MAXLAY)
      REAL   WAOP(MXOWLY)        ! OPTRAN abs coef scaling factor

      REAL CO2MLT(MAXLAY)        ! CO2 perturbation multiplier
      REAL SO2MLT(MAXLAY)        ! SO2 perturbation multiplier
      REAL HNOMLT(MAXLAY)        ! HNO3 perturbation multiplier
      REAL N2OMLT(MAXLAY)        ! N2O perturbation multiplier
      REAL NH3MLT(MAXLAY)        ! NH3 perturbation multiplier
      REAL HDOMLT(MAXLAY)        ! HDO perturbation multiplier
      REAL CO2TOP                ! top layers CO2 mixing ratio

      INTEGER INDCO2(MXCHAN)            ! chan indices for CO2 pert
      REAL COFCO2(  NCO2,MAXLAY,MXCHNC) ! coefs for CO2 pert
      INTEGER INDSO2(MXCHAN)            ! chan indices for SO2 pert
      REAL COFSO2(  NSO2,MAXLAY,MXCHNS) ! coefs for SO2 pert
      INTEGER INDHDO(MXCHAN)            ! chan indices for HDO pert
      REAL COFHDO(  NHDO,MAXLAY,MXCHND) ! coefs for HDO pert
      INTEGER INDHNO(MXCHAN)            ! chan indices for HNO3 pert
      REAL COFHNO( NHNO3,MAXLAY,MXCHNH) ! coefs for HNO3 pert
      INTEGER INDN2O(MXCHAN)            ! chan indices for N2O pert
      REAL COFN2O(  NN2O,MAXLAY,MXCHNN) ! coefs for N2O pert
      INTEGER INDNH3(MXCHAN)            ! chan indices for NH3 pert
      REAL COFNH3(  NNH3,MAXLAY,MXCHNA) ! coefs for NH3 pert
      INTEGER INDH2O(MXCHAN)            ! chan indices for OPTRAN H2O
      REAL COFH2O(  NH2O,MXOWLY,MXCHNW) ! coefs for OPTRAN H2O

      INTEGER QUICKCLIST1(MXCHAN) ! list of set1 channels
      INTEGER QUICKCLIST2(MXCHAN) ! list of set2 channels
      INTEGER QUICKCLIST3(MXCHAN) ! list of set3 channels
      INTEGER QUICKCLIST4(MXCHAN) ! list of set4 channels
      INTEGER QUICKCLIST5(MXCHAN) ! list of set5 channels
      INTEGER QUICKCLIST6(MXCHAN) ! list of set6 channels
      INTEGER QUICKCLIST7(MXCHAN) ! list of set7 channels
      INTEGER INDCHN(MXCHAN)  ! array indices for all channels
      INTEGER   LBOT             ! bottom layer index number

      INTEGER  NCHN1         ! # of set1 channels
      INTEGER  NCHN2         ! # of set2 channels
      INTEGER  NCHN3         ! # of set3 channels
      INTEGER  NCHN4         ! # of set4 channels
      INTEGER  NCHN5         ! # of set5 channels
      INTEGER  NCHN6         ! # of set6 channels
      INTEGER  NCHN7         ! # of set7 channels

      INTEGER CLIST1(MXCHN1) ! list of set1 channels
      INTEGER CLIST2(MXCHN2) ! list of set2 channels
      INTEGER CLIST3(MXCHN3) ! list of set3 channels
      INTEGER CLIST4(MXCHN4) ! list of set4 channels
      INTEGER CLIST5(MXCHN5) ! list of set5 channels
      INTEGER CLIST6(MXCHN6) ! list of set6 channels
      INTEGER CLIST7(MXCHN7) ! list of set7 channels

c local
      INTEGER III

C       ----------------------------------
C       Calculate the layer transmittances
C       ----------------------------------
C       Calculate TAU for set 1 thru 7

      IF (DEBUG) THEN
        DO III = 1,NCHAN
          !! print iI,h.ichan(iI),h.vchan(Ii)
          print *,III,LSTCHN(III),FREQ(III)
        END DO
      END IF

C     ---------------------------
C     Loop on channel (frequency)
C     eventually fills in matrix as    X=1,2,3,4,5,6,7
C      DO I=1,NCHNX
C         J=INDCHN( CLISTX(I) )
C         DO L = 1,100
C            Calc layer-to-space optical depth
C            KZ=KZ + KLAYER
C            TAUZ(ILAY,J)=KZ
C          END DO
C        END DO
C     ---------------------------

C     print *,INDCHN(CLIST1(1:NCHN1))
C     stop
!     DO II = 1,NCHN1
!       !! print iI,h.ichan(iI),h.vchan(Ii)
!       print *,II,LSTCHN(II),FREQ(II),CLIST1(II)
!     END DO
        
      III = QUICKCLIST1(I) ! III = intersect(I,INDCHN(CLIST1(1:NCHN1)), NCHN1)
      IF (III .GT. 0) THEN
        CALL YCALT1( INDCHN,  LBOT,   NCHN1, CLIST1,  COEF1,
     $       FIXMUL, CONPRD, FPRED1, WPRED1, DPRED, OPRED1, TRCPRD,
     $       INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $       INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT,
     $       INDNH3, COFNH3, NH3MLT, INDHDO, COFHDO, HDOMLT,
     $       INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX, LOPLOW,
     $       LOPUSE,   WAOP,   DAOP, WAANG,     TAU,   TAUZ,  III)
      END IF
c         print *,'oneNew  ',tau(10,1),tauz(10,1)

      III = QUICKCLIST2(I) ! III = intersect(I,INDCHN(CLIST2(1:NCHN2)), NCHN2)
      IF (III .GT. 0) THEN  
        CALL YCALT2( INDCHN, LBOT,   NCHN2, CLIST2,  COEF2,
     $      FIXMUL, CONPRD, FPRED2, OPRED2, WPRED2, DPRED, TRCPRD,
     $      INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $      INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, 
     $      INDNH3, COFNH3, NH3MLT, INDHDO, COFHDO, HDOMLT,TAU, TAUZ, 
     $      III)
      END IF
c         print *,'twoNew  ',tau(10,1),tauz(10,1)

      III = QUICKCLIST3(I) ! III = intersect(I,INDCHN(CLIST3(1:NCHN3)), NCHN3)
      IF (III .GT. 0) THEN  
        CALL YCALT3( INDCHN,   LBOT,  NCHN3, CLIST3,  COEF3,
     $       FIXMUL, CONPRD, FPRED3, MPRED3, WPRED3, DPRED, TRCPRD,
     $       INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT,
     $       INDN2O, COFN2O, N2OMLT, INDNH3, COFNH3, NH3MLT,
     $       INDHDO, COFHDO, HDOMLT, INDH2O, H2OPRD, COFH2O, 
     $       LOPMIN, LOPMAX, LOPLOW, LOPUSE,
     $         WAOP,   DAOP,  WAANG,    TAU,   TAUZ, III)
      END IF
c         print *,'threeNew  ',tau(10,1),tauz(10,1)

      III = QUICKCLIST4(I) !  III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4)
      IF (III .GT. 0) THEN  
        CALL YCALT4(INDCHN,   LBOT,  NCHN4, CLIST4,
     $       COEF4, FIXMUL, CONPRD, FPRED4, CPRED4, OPRED4, WPRED4,
     $       TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $       TAU,   TAUZ, III)
      END IF
c         print *,'fourNew  ',tau(10,1),tauz(10,1)

      III = QUICKCLIST5(I) ! III = intersect(I,INDCHN(CLIST5(1:NCHN5)), NCHN5)
      IF (III .GT. 0) THEN  
        CALL YCALT5(INDCHN,   LBOT,  NCHN5, CLIST5,
     $       COEF5, FIXMUL, CONPRD, FPRED5, WPRED5, OPRED5, 
     $       TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $       TAU,   TAUZ, III )
      END IF
c         print *,'fiveNew  ',tau(10,1),tauz(10,1)

      III = QUICKCLIST6(I) ! III = intersect(I,INDCHN(CLIST6(1:NCHN6)), NCHN6)
      IF (III .GT. 0) THEN  
        CALL YCALT6(INDCHN,   LBOT,  NCHN6, CLIST6,
     $       COEF6, FIXMUL, CONPRD, FPRED6, WPRED6, OPRED6, DPRED, TRCPRD,
     $      INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $      INDN2O, COFN2O, N2OMLT, INDHDO, COFHDO, HDOMLT,  TAU,  TAUZ, 
     $      III )
      END IF
c         print *,'sixNew  ',tau(10,1),tauz(10,1)

      III = QUICKCLIST7(I) ! III = intersect(I,INDCHN(CLIST7(1:NCHN7)), NCHN7)
      IF (III .GT. 0) THEN  
        CALL YCALT7(INDCHN,   LBOT,  NCHN7, CLIST7,
     $       COEF7, FIXMUL, CONPRD, FPRED7, WPRED7, OPRED7,
     $       TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $       TAU,   TAUZ, III )
      END IF
c         print *,'sevenNew  ',tau(10,1),tauz(10,1)

c***********************************************************************
      IF (DOSUN) THEN

C       Calc fake TAUZSN for sets 1, 2, and 3
c       III = intersect(I,INDFAK(1:NFAKE), NFAKE)  !! so I = INDFAK(III)
        III = QUICKINDFAK(I)
        IF (III .GT. 0) THEN
c          print *,'indfak',I,III,INDFAK(I),INDFAK(III),QUICKINDFAK(I)
          CALL FAKETZ( NFAKE, INDFAK, LBOT, TAUZ, SECANG,
     $        SECSUN, TAUZSN, III)
        END IF

        III = QUICKCLIST4(I) ! III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4) 
                               ! so I = INDCHN(CLIST4(III))
        IF (III .GT. 0) THEN  
          CALL YCALT4(INDCHN,   LBOT,  NCHN4, CLIST4,
     $         COEF4, FIXMUL, SUNCONPRD, SUNFPRED4, SUNCPRED4, SUNOPRED4, SUNWPRED4,
     $         SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $         TAUZSN, TAUZSN, III )
C       ^^^^^^  ^^^^^^
C           dummy   actual
        END IF

        III = QUICKCLIST5(I)  ! III = intersect(I,INDCHN(CLIST5(1:NCHN5)), NCHN5)
        IF (III .GT. 0) THEN  
          CALL YCALT5(INDCHN,   LBOT,  NCHN5, CLIST5,
     $         COEF5, FIXMUL, SUNCONPRD, SUNFPRED5, SUNWPRED5, SUNOPRED5,
     $         SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $         TAUZSN, TAUZSN, III )
        END IF

        III = QUICKCLIST6(I) !  III = intersect(I,INDCHN(CLIST6(1:NCHN6)), NCHN6)
        IF (III .GT. 0) THEN  
          CALL YCALT6(INDCHN,   LBOT,  NCHN6, CLIST6,
     $          COEF6, FIXMUL, SUNCONPRD, SUNFPRED6, SUNWPRED6, SUNOPRED6, DPRED,
     $          SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $          INDN2O, COFN2O, N2OMLT, INDHDO, COFHDO, HDOMLT, TAUZSN, 
     $          TAUZSN, III )
        END IF

        III = QUICKCLIST7(I) !  III = intersect(I,INDCHN(CLIST7(1:NCHN7)), NCHN7)
        IF (III .GT. 0) THEN  
          CALL YCALT7(INDCHN,   LBOT,  NCHN7, CLIST7,
     $          COEF7, FIXMUL, SUNCONPRD, SUNFPRED7, SUNWPRED7, SUNOPRED7,
     $          SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $          TAUZSN, TAUZSN, III )
        END IF
c       print *,'sun 7',TAU(LBOT,I),TAUZ(LBOT-1,I),TAUZSN(LBOT,I),TAUZSN(LBOT-1,I)

        IF (SUNFDG .GT. 1.0001) THEN
          TAUZSN(1:LBOT,1:NCHAN)=TAUZSN(1:LBOT,1:NCHAN)*SUNFDG
         ENDIF

      ELSE
C       DOSUN = 'FALSE'; No sun; set the sun surface-to-space trans to zero
        SUNCOS=0.0
        TAUZSN(1:LBOT,I) = 0.0
      ENDIF

      RETURN
      END
      
