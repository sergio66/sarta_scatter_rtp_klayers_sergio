      SUBROUTINE CALC_LAYER_TRANS_YCALTODX_1_7(
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
     $ DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED,H2OJACPRD,DAOPJAC,
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12,
     $ TAU, TAUZ, TAUZSN)

      IMPLICIT NONE
      INCLUDE "incFTC.f"

c output
       REAL    TAU(MAXLAY,MXCHAN) ! chan layer effective optical depth
       REAL   TAUZ(MAXLAY,MXCHAN) ! chan surface-to-space trans
       REAL TAUZSN(MAXLAY,MXCHAN) ! sun space-to-surface-to-space OD
       REAL DTAU_DTZ(MAXLAY,MXCHAN),DTAU_DG1(MAXLAY,MXCHAN),DTAU_DG2(MAXLAY,MXCHAN),DTAU_DG3(MAXLAY,MXCHAN)
       REAL DTAU_DG4(MAXLAY,MXCHAN),DTAU_DG5(MAXLAY,MXCHAN),DTAU_DG6(MAXLAY,MXCHAN),DTAU_DG9(MAXLAY,MXCHAN)
       REAL DTAU_DG12(MAXLAY,MXCHAN)

c input
       LOGICAL DOJAC
       INTEGER NWANTJ          ! number of wanted jacs (default 0=none)
       INTEGER  LISTJ(MAXPRO)  ! list of wanted channels
       REAL CO2JACMLT(MAXLAY)
       REAL SO2JACMLT(MAXLAY)
       REAL HNOJACMLT(MAXLAY)
       REAL N2OJACMLT(MAXLAY)
       REAL NH3JACMLT(MAXLAY)
       REAL HDOJACMLT(MAXLAY)
       !!! first index is the d/dT   second deriv is the d/dQ
       REAL CONJACPRD(MAXJAC, N1CON,MAXLAY)
       REAL FJACPRED1(MAXJAC, N1FIX,MAXLAY)
       REAL FJACPRED2(MAXJAC, N2FIX,MAXLAY)
       REAL FJACPRED3(MAXJAC, N3FIX,MAXLAY)
       REAL FJACPRED4(MAXJAC, N4FIX,MAXLAY)
       REAL FJACPRED5(MAXJAC, N5FIX,MAXLAY)
       REAL FJACPRED6(MAXJAC, N6FIX,MAXLAY)
       REAL FJACPRED7(MAXJAC, N7FIX,MAXLAY)
       REAL WJACPRED1(MAXJAC, N1H2O,MAXLAY)
       REAL WJACPRED2(MAXJAC, N2H2O,MAXLAY)
       REAL WJACPRED3(MAXJAC, N3H2O,MAXLAY)
       REAL WJACPRED4(MAXJAC, N4H2O,MAXLAY)
       REAL WJACPRED5(MAXJAC, N5H2O,MAXLAY)
       REAL WJACPRED6(MAXJAC, N6H2O,MAXLAY)
       REAL WJACPRED7(MAXJAC, N7H2O,MAXLAY)
       REAL OJACPRED1(MAXJAC,  N1O3,MAXLAY)
       REAL OJACPRED2(MAXJAC,  N2O3,MAXLAY)
       REAL OJACPRED4(MAXJAC,  N4O3,MAXLAY)
       REAL OJACPRED5(MAXJAC,  N5O3,MAXLAY)
       REAL OJACPRED6(MAXJAC,  N6O3,MAXLAY)
       REAL OJACPRED7(MAXJAC,  N7O3,MAXLAY)
       REAL  DJACPRED(MAXJAC,  NHDO,MAXLAY)
       REAL MJACPRED3(MAXJAC, N3CH4,MAXLAY)
       REAL CJACPRED4(MAXJAC,  N4CO,MAXLAY)
       REAL TRCJACPRD(MAXJAC,NTRACE,MAXLAY)
       REAL H2OJACPRD(OPTRANJAC,  NH2O,MXOWLY)    !! OPTRAN, used by ycalt1_od, ycalt3_od
       REAL  DAOPJAC(OPTRANJAC, MAXLAY)           !! OPTRAN, used by ycalt1_od, ycalt3_od

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
      INTEGER III,intersect

C       ----------------------------------
C       Calculate the layer transmittances
C       ----------------------------------
C       Calculate TAU for set 1 thru 7

c      so eg if you are doing IASI blah.rtp_2, then you should see    h.ichan,h.vchan
c           1        4232   1702.750
c           2        4233   1703.000
c           .......
c        4229        8460   2759.750
c        4230        8461   2760.000
c
c      IF (DEBUG) THEN
c        DO III = 1,NCHAN
c          !! print iI,h.ichan(iI),h.vchan(Ii)
c          print *,III,LSTCHN(III),FREQ(III)
c       END DO
c       stop
c      END IF

c      so eg if you are doing IASI blah.rtp_2, then you should see     iasi_rtp2_channels_set1_7
c      IF (DEBUG) THEN
c        DO III = 1,NCHAN
c          write(*,'(I6,I6,F12.5,7(I6))') III,LSTCHN(III),FREQ(III),QUICKCLIST1(III),QUICKCLIST2(III),
c     c            QUICKCLIST3(III),QUICKCLIST4(III),QUICKCLIST5(III),QUICKCLIST6(III),QUICKCLIST7(III)
c       END DO
c       stop
c      END IF

c************************************************************************
C this was the orig code for eg calt4.f
C     DO I=1,NCHN4
C         Index for TAU
C          J=INDCHN( CLIST4(I) )
c
C     this is the new code for eg ycalt4_od.f : IY = passed in using
c      III = QUICKCLIST4(I)      
c      III = QUICKCLIST4(LSTCHN(I))
c      III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4)
c      III = intersect(LSTCHN(I),INDCHN(CLIST4(1:NCHN4)), NCHN4)                       
C         I = IY     
C          J=INDCHN( CLIST4(I) )
c
c      
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
c
C     print *,INDCHN(CLIST1(1:NCHN1))
C     stop
c     DO II = 1,NCHN1
c       !! print iI,h.ichan(iI),h.vchan(Ii)
c       print *,II,LSTCHN(II),FREQ(II),CLIST1(II)
c     END DO
c************************************************************************

c        DO III = 1,NCHAN
c          write(*,'(I6,I6,F12.5,7(I6))') III,LSTCHN(III),FREQ(III),QUICKCLIST1(III),QUICKCLIST2(III),
c     c            QUICKCLIST3(III),QUICKCLIST4(III),QUICKCLIST5(III),QUICKCLIST6(III),QUICKCLIST7(III)
c       END DO

c III = QUICKCLISTN(I) is correct      
      III = QUICKCLIST1(I)            !III = intersect(I,INDCHN(CLIST1(1:NCHN1)), NCHN1)      
c III = QUICKCLISTN(I) is correct      
cc      III = QUICKCLIST1(LSTCHN(I))
cc      III = intersect(I,INDCHN(CLIST1(1:NCHN1)), NCHN1)
cc      III = intersect(LSTCHN(I),INDCHN(CLIST1(1:NCHN1)), NCHN1)                  
cc      print *, 'LIST 1 : I   III = BLAH(I) = ',I,III
cc      print *,'INDCHN = ',INDCHN
cc      print *,'LSTCHN = ',LSTCHN

c      IF (I .EQ. 254) print *,'moo moo moo moo QUICKLIST 1',I,III
      IF (III .GT. 0) THEN
        CALL YCALT1( INDCHN,  LBOT,   NCHN1, CLIST1,  COEF1,
     $       FIXMUL, CONPRD, FPRED1, WPRED1, DPRED, OPRED1, TRCPRD,
     $       INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $       INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT,
     $       INDNH3, COFNH3, NH3MLT, INDHDO, COFHDO, HDOMLT,
     $       INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX, LOPLOW,
     $       LOPUSE,   WAOP,   DAOP, WAANG,     TAU,   TAUZ,  III,
     $    DOJAC,LISTJ,NWANTJ,SECANG,CONJACPRD,DJACPRED,H2OJACPRD,DAOPJAC,
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
      END IF
c      print *,'sseerrggiioo oneNew  ',tau(10,1),tauz(10,4230)

c III = QUICKCLISTN(I) is correct      
      III = QUICKCLIST2(I)      ! III = intersect(I,INDCHN(CLIST2(1:NCHN2)), NCHN2)
c III = QUICKCLISTN(I) is correct      
cc      III = QUICKCLIST2(LSTCHN(I))
cc      III = intersect(I,INDCHN(CLIST2(1:NCHN2)), NCHN2)
cc      III = intersect(LSTCHN(I),INDCHN(CLIST2(1:NCHN2)), NCHN2)                  
cc      print *, 'LIST 2 : I   III = BLAH(I) = ',I,III      
cc      IF (I .EQ. 254) print *,'moo moo moo moo QUICKLIST 2',I,III
      IF (III .GT. 0) THEN  
        CALL YCALT2( INDCHN, LBOT,   NCHN2, CLIST2,  COEF2,
     $      FIXMUL, CONPRD, FPRED2, OPRED2, WPRED2, DPRED, TRCPRD,
     $      INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $      INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, 
     $      INDNH3, COFNH3, NH3MLT, INDHDO, COFHDO, HDOMLT,TAU, TAUZ, III,
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
      END IF
c      print *,'sseerrggiioo twoNew  ',tau(10,1),tauz(10,4230)
 
c III = QUICKCLISTN(I) is correct           
      III = QUICKCLIST3(I)      ! III = intersect(I,INDCHN(CLIST3(1:NCHN3)), NCHN3)
c III = QUICKCLISTN(I) is correct      
cc      III = QUICKCLIST3(LSTCHN(I))
cc      III = intersect(I,INDCHN(CLIST3(1:NCHN3)), NCHN3)
cc      III = intersect(LSTCHN(I),INDCHN(CLIST3(1:NCHN3)), NCHN3)                  
cc      print *, 'LIST 3 : I   III = BLAH(I) = ',I,III      
      IF (III .GT. 0) THEN  
        CALL YCALT3( INDCHN,   LBOT,  NCHN3, CLIST3,  COEF3,
     $       FIXMUL, CONPRD, FPRED3, MPRED3, WPRED3, DPRED, TRCPRD,
     $       INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT,
     $       INDN2O, COFN2O, N2OMLT, INDNH3, COFNH3, NH3MLT,
     $       INDHDO, COFHDO, HDOMLT, INDH2O, H2OPRD, COFH2O, 
     $       LOPMIN, LOPMAX, LOPLOW, LOPUSE,
     $         WAOP,   DAOP,  WAANG,    TAU,   TAUZ, III,
     $    DOJAC,LISTJ,NWANTJ,SECANG,CONJACPRD,DJACPRED,H2OJACPRD,DAOPJAC,
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
      END IF
c      print *,'sseerrggiioo threeNew  ',tau(10,1),tauz(10,4230)
      
c III = QUICKCLISTN(I) is correct      
      III = QUICKCLIST4(I)      !  III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4)
c III = QUICKCLISTN(I) is correct      
cc      III = QUICKCLIST4(LSTCHN(I))
cc      III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4)
cc      III = intersect(LSTCHN(I),INDCHN(CLIST4(1:NCHN4)), NCHN4)                  
cc      print *, 'LIST 4 : I   III = BLAH(I) = ',I,III      
      IF (III .GT. 0) THEN  
        CALL YCALT4(INDCHN,   LBOT,  NCHN4, CLIST4,
     $       COEF4, FIXMUL, CONPRD, FPRED4, CPRED4, OPRED4, WPRED4,
     $       TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $       TAU,   TAUZ, III,
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
      END IF
c      print *,'sseerrggiioo fourNew  ',tau(10,1),tauz(10,4230)
      
c III = QUICKCLISTN(I) is correct      
      III = QUICKCLIST5(I)      ! III = intersect(I,INDCHN(CLIST5(1:NCHN5)), NCHN5)
c III = QUICKCLISTN(I) is correct      
cc      III = QUICKCLIST5(LSTCHN(I))
cc      III = intersect(I,INDCHN(CLIST5(1:NCHN5)), NCHN5)
cc      III = intersect(LSTCHN(I),INDCHN(CLIST5(1:NCHN5)), NCHN5)                  
cc      print *, 'LIST 5 : I   III = BLAH(I) = ',I,III      
      IF (III .GT. 0) THEN  
        CALL YCALT5(INDCHN,   LBOT,  NCHN5, CLIST5,
     $       COEF5, FIXMUL, CONPRD, FPRED5, WPRED5, OPRED5, 
     $       TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $       TAU,   TAUZ, III,
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
      END IF
c      print *,'sseerrggiioo fiveNew  ',tau(10,1),tauz(10,4230)
      
c III = QUICKCLISTN(I) is correct      
      III = QUICKCLIST6(I)      ! III = intersect(I,INDCHN(CLIST6(1:NCHN6)), NCHN6)
c III = QUICKCLISTN(I) is correct      
cc      III = QUICKCLIST6(LSTCHN(I))
cc      III = intersect(I,INDCHN(CLIST6(1:NCHN6)), NCHN6)
cc      III = intersect(LSTCHN(I),INDCHN(CLIST6(1:NCHN6)), NCHN6)                  
cc      print *, 'LIST 6 : I   III = BLAH(I) = ',I,III      
      IF (III .GT. 0) THEN  
        CALL YCALT6(INDCHN,   LBOT,  NCHN6, CLIST6,
     $       COEF6, FIXMUL, CONPRD, FPRED6, WPRED6, OPRED6, DPRED, TRCPRD,
     $      INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $      INDN2O, COFN2O, N2OMLT, INDHDO, COFHDO, HDOMLT,  TAU,  TAUZ, III,
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
      END IF
c      print *,'sseerrggiioo sixNew  ',tau(10,1),tauz(10,4230)
      
c III = QUICKCLISTN(I) is correct      
      III = QUICKCLIST7(I)      ! III = intersect(I,INDCHN(CLIST7(1:NCHN7)), NCHN7)
c III = QUICKCLISTN(I) is correct      
cc      III = QUICKCLIST7(LSTCHN(I))
cc      III = intersect(I,INDCHN(CLIST7(1:NCHN7)), NCHN7)
cc      III = intersect(LSTCHN(I),INDCHN(CLIST7(1:NCHN7)), NCHN7)            
cc      print *, 'LIST 7 : I   III = BLAH(I) = ',I,III      
      IF (III .GT. 0) THEN  
        CALL YCALT7(INDCHN,   LBOT,  NCHN7, CLIST7,
     $       COEF7, FIXMUL, CONPRD, FPRED7, WPRED7, OPRED7,
     $       TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $       TAU,   TAUZ, III, 
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
      END IF
c      print *,'sseerrggiioo sevenNew  ',tau(10,1),tauz(10,4230)
      
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

        III = QUICKCLIST4(I)    ! III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4)
                                !                       so I = INDCHN(CLIST4(III))
c        III = QUICKCLIST4(LSTCHN(I))
c        III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4)
c        III = intersect(LSTCHN(I),INDCHN(CLIST4(1:NCHN4)), NCHN4)
c        print *, 'SUN LIST 4 : I   III = BLAH(I) = ',I,III        
        IF (III .GT. 0) THEN  
          CALL YCALT4(INDCHN,   LBOT,  NCHN4, CLIST4,
     $         COEF4, FIXMUL, SUNCONPRD, SUNFPRED4, SUNCPRED4, SUNOPRED4, SUNWPRED4,
     $         SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $         TAUZSN, TAUZSN, III,
C       ^^^^^^  ^^^^^^
C           dummy   actual
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
        END IF

        III = QUICKCLIST5(I)    ! III = intersect(I,INDCHN(CLIST5(1:NCHN5)), NCHN5)
c        III = QUICKCLIST5(LSTCHN(I))
c        III = intersect(I,INDCHN(CLIST5(1:NCHN5)), NCHN5)
c        III = intersect(LSTCHN(I),INDCHN(CLIST5(1:NCHN5)), NCHN5)
c        print *, 'SUN LIST 5 : I   III = BLAH(I) = ',I,III                
        IF (III .GT. 0) THEN  
          CALL YCALT5(INDCHN,   LBOT,  NCHN5, CLIST5,
     $         COEF5, FIXMUL, SUNCONPRD, SUNFPRED5, SUNWPRED5, SUNOPRED5,
     $         SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $         TAUZSN, TAUZSN, III, 
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
        END IF

        III = QUICKCLIST6(I)    !  III = intersect(I,INDCHN(CLIST6(1:NCHN6)), NCHN6)
c        III = QUICKCLIST7(LSTCHN(I))
c        III = intersect(I,INDCHN(CLIST6(1:NCHN6)), NCHN6)
c        III = intersect(LSTCHN(I),INDCHN(CLIST6(1:NCHN7)), NCHN6)
c        print *, 'SUN LIST 6 : I   III = BLAH(I) = ',I,III                
        IF (III .GT. 0) THEN  
          CALL YCALT6(INDCHN,   LBOT,  NCHN6, CLIST6,
     $          COEF6, FIXMUL, SUNCONPRD, SUNFPRED6, SUNWPRED6, SUNOPRED6, DPRED,
     $          SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $          INDN2O, COFN2O, N2OMLT, INDHDO, COFHDO, HDOMLT, 
     $          TAUZSN, TAUZSN, III,
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
        END IF

        III = QUICKCLIST7(I)    !  III = intersect(I,INDCHN(CLIST7(1:NCHN7)), NCHN7)
c        III = QUICKCLIST7(LSTCHN(I))
c        III = intersect(I,INDCHN(CLIST7(1:NCHN7)), NCHN7)
c        III = intersect(LSTCHN(I),INDCHN(CLIST7(1:NCHN7)), NCHN7)
c        print *, 'SUN LIST 7 : I   III = BLAH(I) = ',I,III                
        IF (III .GT. 0) THEN  
          CALL YCALT7(INDCHN,   LBOT,  NCHN7, CLIST7,
     $          COEF7, FIXMUL, SUNCONPRD, SUNFPRED7, SUNWPRED7, SUNOPRED7,
     $          SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $          TAUZSN, TAUZSN, III, 
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $    DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12)
        END IF
c       print *,'sun 7',TAU(LBOT,I),TAUZ(LBOT-1,I),TAUZSN(LBOT,I),TAUZSN(LBOT-1,I)

        IF (SUNFDG .GT. 1.0001) THEN
          TAUZSN(1:LBOT,I)=TAUZSN(1:LBOT,I)*SUNFDG
         ENDIF

      ELSE
C       DOSUN = 'FALSE'; No sun; set the sun surface-to-space trans to zero
        SUNCOS=0.0
        TAUZSN(1:LBOT,I) = 0.0
      ENDIF

c      print *,'sseerrggiioo calc_layer_trans_ycaltodx_1_7.f'
c      print *,'tau(:,1)    = ',tau(:,1)
c      print *,'tau(:,4230) = ',tau(:,4230)
c      DO III = 1,4230
c        print *,III,tau(20,III)
c      end do
c      print *,'calc_layer_trans_ycaltodx_1_7.f : stop'
c      stop
     
      RETURN
      END
      
