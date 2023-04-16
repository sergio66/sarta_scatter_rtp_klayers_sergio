!
!     This module contains local variables in sarta.f
!     Sung-Yung Lee
!     November 2006
!

MODULE AIRS_SARTA_VARIABLES

	USE INCFTC
	USE REF_PROFILES

	IMPLICIT NONE

	!
	! Variable used in AIRS SARTA package
	!

	LOGICAL :: SARTA_INITIALIZED  ! has sarta been initialized?
	INTEGER :: IOUN         ! I/O unit number

	! for FNMIE
	CHARACTER*240 VCLOUD        ! cloud version string
	INTEGER MIETYP(NMIETY)      ! mie type
	CHARACTER*79 FNMIEA(NMIETY) ! mie absorption filenames
	CHARACTER*79 FNMIEE(NMIETY) ! mie extinction filenames
	CHARACTER*79 FNMIEG(NMIETY) ! mie asymmetry filenames


	! for OPNRTP
	! INTEGER  PTYPE         ! profile type
	! INTEGER  NCHAN         ! # of selected channels ----- Taken out of this module by SYLee on Dec 14, 2007
	REAL     FCHAN(MXCHAN) ! chan center frequency
	INTEGER INDCHN(MXCHAN) ! array indices for all channels
	INTEGER LSTCHN(MXCHAN) ! list of selected channels
	INTEGER IH2O           ! index of H2O in gamnt
	INTEGER IO3            ! index of O3 in gamnt
	INTEGER ICO            ! index of CO in gamnt
	INTEGER ICH4           ! index of CH4 in gamnt
	INTEGER ICO2           ! index of CO2 in gamnt
	INTEGER ISO2           ! index of SO2 in gamnt
	INTEGER IHNO3          ! index of HNO3 in gamnt
	INTEGER IN2O           ! index of N2O in gamnt
	!INTEGER IOPCI          ! input RTP unit
	!INTEGER IOPCO          ! output RTP unit

	! for RDCLDT
	INTEGER MIENPS(NMIETY)            ! number of particle sizes
	REAL  MIEPS(MXMIEA,NMIETY)        ! Mie particle size for table
	REAL MIEABS(MXCHAN,MXMIEA,NMIETY) ! Mie absorption table
	REAL MIEEXT(MXCHAN,MXMIEA,NMIETY) ! Mie extinction table
	REAL MIEASY(MXCHAN,MXMIEA,NMIETY) ! Mie asymmetry table
   

	! for RDCOEF           ! Info for selected channels only
	INTEGER SETCHN(MXCHAN) ! set # for each channel
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
	INTEGER LABOVE(MXCHAN) ! chan downwelling thermal layer above
	REAL   FREQ(MXCHAN)    ! chan center frequency
	! REAL   C1V3(MXCHAN)    ! C1*(FREQ(I)**3)		Precomputed Plank function constant
	! REAL   C2V(MXCHAN)     ! C2*FREQ(I)
	REAL   C1V3    ! C1*(FREQ**3)		 Plank function constant
	REAL   C2V     ! C2*FREQ
	REAL  COEF1(N1COEF,MAXLAY,MXCHN1) ! coefs for set1 chans
	REAL  COEF2(N2COEF,MAXLAY,MXCHN2) ! coefs for set2 chans
	REAL  COEF3(N3COEF,MAXLAY,MXCHN3) ! coefs for set3 chans
	REAL  COEF4(N4COEF,MAXLAY,MXCHN4) ! coefs for set4 chans
	REAL  COEF5(N5COEF,MAXLAY,MXCHN5) ! coefs for set5 chans
	REAL  COEF6(N6COEF,MAXLAY,MXCHN6) ! coefs for set6 chans
	REAL  COEF7(N7COEF,MAXLAY,MXCHN7) ! coefs for set7 chans
	REAL  COEFF(NFCOEF,MXCHAN)        ! coefs for chan "F" factor
	INTEGER INDCO2(MXCHAN)            ! chan indices for CO2 pert
	REAL COFCO2(  NCO2,MAXLAY,MXCHNC) ! coefs for CO2 pert
	INTEGER INDSO2(MXCHAN)            ! chan indices for SO2 pert
	REAL COFSO2(  NSO2,MAXLAY,MXCHNS) ! coefs for SO2 pert
	INTEGER INDHNO(MXCHAN)            ! chan indices for HNO3 pert
	REAL COFHNO( NHNO3,MAXLAY,MXCHNH) ! coefs for HNO3 pert
	INTEGER INDN2O(MXCHAN)            ! chan indices for N2O pert
	REAL COFN2O(  NN2O,MAXLAY,MXCHNN) ! coefs for N2O pert
	INTEGER INDH2O(MXCHAN)            ! chan indices for OPTRAN H2O
	REAL   WAZOP(MXOWLY)              ! OPTRAN water l-to-s amounts
	REAL  WAVGOP(NOWAVG,MXOWLY)       ! OPTRAN raw predictor averages
	REAL COFH2O(  NH2O,MXOWLY,MXCHNW) ! coefs for OPTRAN H2O
	REAL     FX(MAXLAY)               ! fixed gases adjustment
	INTEGER NCHNTE                    ! number of non-LTE channels
	INTEGER CLISTN(MXCNTE)            ! non-LTE channel list
	REAL  COEFN(NNCOEF,MXCNTE)        ! non-LTE coefficients

	! for FAKETZ
	INTEGER  NFAKE         ! # of channels to "fake"
	INTEGER INDFAK(MXCHAN) ! indices of channels to fake

	! for RDPROF; reference profile  ! Variables now in REF_AMNTS module  ! FWI 12/20/17
	!	CHARACTER*40  RPNAM ! ref prof name/ID
	!	REAL   RALT(MAXLAY) ! ref prof layer altitude
	!	REAL    RDZ(MAXLAY) ! ref prof layer thickness
	!	REAL  RPRES(MAXLAY) ! ref prof layer average pressure
	!	REAL  RTEMP(MAXLAY) ! ref prof layer average temperature
	!	REAL RFAMNT(MAXLAY) ! ref prof layer "fixed" (CO2) amount
	!	REAL RWAMNT(MAXLAY) ! ref prof layer water (H2O) amount
	!	REAL ROAMNT(MAXLAY) ! ref prof layer ozone (O3) amount
	!	REAL RCAMNT(MAXLAY) ! ref prof layer carbon monoxide (CO) amount
	!	REAL RMAMNT(MAXLAY) ! ref prof layer methane (CH4) amount
	!	REAL RSAMNT(MAXLAY) ! ref prof layer sulfer dioxide (SO2) amount
	!	REAL RHAMNT(MAXLAY) ! ref prof layer nitric acid (HNO3) amount
	!	REAL RNAMNT(MAXLAY) ! ref prof layer nitrous oxide (N2O) amount

	! for RDRTP; profile to calculate
	INTEGER NLAY / 0 /  ! number of layers in profile
	!      REAL LAT            ! prof latitude
	!      REAL LON            ! prof longitude
	!      REAL    ALT(MAXLAY) ! prof layer altitudes
	!      REAL   TEMP(MAXLAY) ! prof layer average temperature
	!      REAL  WAMNT(MAXLAY) ! prof layer water (H2O) amount
	!      REAL  OAMNT(MAXLAY) ! prof layer ozone (O3) amount
	!      REAL  CAMNT(MAXLAY) ! prof layer carbon monoxide (CO) amount
	!      REAL  MAMNT(MAXLAY) ! prof layer methane (CH4) amount
	!      REAL  FAMNT(MAXLAY) ! prof layer CO2 amount
	!      REAL  SAMNT(MAXLAY) ! prof layer SO2 amount
	!      REAL  HAMNT(MAXLAY) ! prof layer HNO3 amount
	!      REAL  NAMNT(MAXLAY) ! prof layer N2O amount
	!
	!      for surface
	!      INTEGER   LBOT             ! bottom layer index number
	!      INTEGER  NEMIS             ! # of emis pts
	!      INTEGER   NRHO             ! # of rho pts
	!      REAL  PSURF                ! surface pressure
	!      REAL BLMULT                ! bottom layer fractional multiplier
	!      REAL  XEMIS(MXEMIS)        ! emis pts
	!      REAL  FEMIS(MXEMIS)        ! emis freq pts
	!      REAL   XRHO(MXEMIS)        ! reflec pts
	!      REAL   FRHO(MXEMIS)        ! reflec freq pts
	!

	! for MEAN_T
	REAL TPSEUD(MAXLAY)

	! for CALPAR
!	REAL SECANG(MAXLAY)        ! local path angle secant  ! Moved to SECANG_CALC module ! FWI 1/17/18
	REAL FIXMUL(MAXLAY)        ! "fixed" amount multiplier (~1)
!	REAL H2O_CONTINUUM_PRED( N1CON,MAXLAY) ! water continuum predictors
!!	REAL, DIMENSION(:, :), POINTER :: FIXED_PRED1_PTR ! set1 "fixed" predictors
!	REAL FIXED_PRED2( N2FIX,MAXLAY) ! set2 "fixed" predictors
!	REAL FIXED_PRED3( N3FIX,MAXLAY) ! set3 "fixed" predictors
!	REAL FIXED_PRED4( N4FIX,MAXLAY) ! set4 "fixed" predictors
!	REAL FIXED_PRED5( N5FIX,MAXLAY) ! set5 "fixed" predictors
!	REAL FIXED_PRED6( N6FIX,MAXLAY) ! set6 "fixed" predictors
!	REAL FIXED_PRED7( N7FIX,MAXLAY) ! set7 "fixed" predictors
!	REAL H2O_PRED1( N1H2O,MAXLAY) ! set1 water predictors
!	REAL H2O_PRED2( N2H2O,MAXLAY) ! set2 water predictors
!	REAL H2O_PRED3( N3H2O,MAXLAY) ! set3 water predictors
!	REAL H2O_PRED4( N4H2O,MAXLAY) ! set4 water predictors
!	REAL H2O_PRED5( N5H2O,MAXLAY) ! set5 water predictors
!	REAL H2O_PRED6( N6H2O,MAXLAY) ! set6 water predictors
!	REAL H2O_PRED7( N7H2O,MAXLAY) ! set7 water predictors
!	REAL O3_PRED1(  N1O3,MAXLAY) ! set1 ozone predictors
!	REAL O3_PRED2(  N2O3,MAXLAY) ! set2 ozone predictors
!	REAL O3_PRED4(  N4O3,MAXLAY) ! set4 ozone predictors
!	REAL O3_PRED5(  N5O3,MAXLAY) ! set5 ozone predictors
!	REAL O3_PRED6(  N6O3,MAXLAY) ! set6 ozone predictors
!	REAL O3_PRED7(  N7O3,MAXLAY) ! set7 ozone predictors
!	REAL CH4_PRED3( N3CH4,MAXLAY) ! set3 methane predictors
!	REAL CO_PRED4(  N4CO,MAXLAY) ! set4 carbon monoxide predictors
!	REAL TRACEGAS_PRED(NTRACE,MAXLAY) ! trace gas pert perdictors
	REAL CO2MLT(MAXLAY)        ! CO2 perturbation multiplier
	REAL SO2MLT(MAXLAY)        ! SO2 perturbation multiplier
	REAL HNOMLT(MAXLAY)        ! HNO3 perturbation multiplier
	REAL N2OMLT(MAXLAY)        ! N2O perturbation multiplier
	REAL CO2TOP                ! top layers CO2 mixing ratio


	! for CALOWP
	REAL  WAANG(MAXLAY)
	INTEGER LOPMIN
	INTEGER LOPMAX
	REAL H2OPRD(  NH2O,MXOWLY)
	LOGICAL LOPUSE(MXOWLY)
	INTEGER LOPLOW(MAXLAY)
	REAL  DAOP(MAXLAY)

	! for CALT
	REAL    TAU(MAXLAY,MXCHAN) ! chan layer effective trans
	REAL   TAUZ(MAXLAY,MXCHAN)        ! chan surface-to-space trans
	REAL   WAOP(MXOWLY)        ! OPTRAN abs coef scaling factor
	REAL     XZ                ! optical depth multiplier for TAUZ
	LOGICAL   LTAU             ! Calc all layer transmittances?

	! for SETEMS
!	REAL   EMIS(MXCHAN) ! chan surface emissivity
!	REAL RHOSUN(MXCHAN) ! chan reflectivity for sun
!	REAL RHOTHR(MXCHAN) ! chan reflectivity for downwelling thermal

!	REAL CEMIS1(MXCHAN) ! chan surface emissivity cloud1
!	REAL CRHOS1(MXCHAN) ! chan solar reflectivity cloud1
!	REAL CRHOT1(MXCHAN) ! chan thermal reflectivity cloud1
!	REAL CEMIS2(MXCHAN) ! chan surface emissivity cloud2
!	REAL CRHOS2(MXCHAN) ! chan solar reflectivity cloud2
!	REAL CRHOT2(MXCHAN) ! chan thermal reflectivity cloud2
	LOGICAL  LRHOT         ! force refl therm rho=(1-emis)/pi?

!
!      for CALRAD
   REAL SUNFAC         ! sun solid angles times cosine at surface
   REAL RPLNCK(MAXLAY) ! layer Planck
   REAL RSURFE         ! surface emission
   REAL RSURFC         ! black cloud surface emission
   REAL  TRANL(MAXLAY) ! clear air layer transmittance
   REAL  TRANZ(MXCHAN) ! clear air layer-to-space transmittance
   REAL  TRANS(MXCHAN) ! clear air total reflected solar trans
!	REAL    RAD(MXCHAN) ! chan radiance
!	For clear/cloudy radiances
   REAL   RAD0         ! radiance no clouds
   REAL  RADC1         ! radiance cloud1
   REAL  RADC2         ! radiance cloud2
   REAL RADC12         ! radiance cloud1+cloud2
!

!	for CALRAD
!	REAL  TSURF         ! surface temperature
!	REAL   EMIS(MXCHAN) ! chan surface emissivity
!	REAL RHOSUN(MXCHAN) ! chan reflectivity for sun
!	REAL RHOTHR(MXCHAN) ! chan reflectivity for downwelling thermal
!	REAL    RAD(MXCHAN) ! chan radiance
!	REAL     BT(MXCHAN) ! chan brightness temperature
!	REAL rad_atm_up(MXCHAN) ! Atmospheric component of upwelling radiance
!	REAL rad_atm_dn(MXCHAN) ! Atmospheric component of downwelling radiance at the surface
!

!	for RDSUN
	REAL   HSUN(MXCHAN) ! sun radiance (direct from sun)
!
!	Other variables for the sun
!	REAL SUNANG         ! solar zenith angle (at 0 altitude)34u
!	REAL SZALAY         ! solar zenith angle in some layer
!	REAL SUNCOS         ! cosine of sun zenith angle
!	REAL SCOS1          ! cosine of sun zenith angle at layer1
!	REAL SUNFDG         ! fudge factor for large solar angles
!	REAL SECSUN(MAXLAY) ! secant of effective sun local path angle
	REAL DISTES         ! distance of Earth from the sun
!	REAL TAUZSN(MXCHAN) ! chan eff sun angle surface-to-space trans
!	LOGICAL DOSUN       ! do sun calc?
!	
!	for satellite viewing angle
!	REAL    SATANG      ! input satellite scan angle (degrees)
!	REAL    SATZEN      ! input satellite zenith angle (degrees)
!	REAL    SALT        ! input satellite altitude (kilometers)
!	REAL    SVA         ! satellite viewing angle (degrees)
!
!	for RDRTP
!	INTEGER  IPROF      ! profile loop counter
!	LOGICAL  LWANT      ! do you want this profile?

!	used locally only
!	INTEGER      I      ! loop counter
!	INTEGER      L      ! loop counter
!	INTEGER rtpclose    ! for call to RTP close interface routine
!	REAL    EVA         ! (Earth) local viewing angle
!	REAL   CONV    / 1.7453292E-02 /    ! degrees to radians conversion factor
!	REAL ANGMAX         ! maximum allowed viewing angle
!	REAL RJUNK1         ! junk/work
!	REAL RJUNK2         ! another junk/work
	REAL CO2PPM         ! Profile mean dry air CO2 mixing ratio
	REAL PLAY(MAXLAY)   ! layer mean pressure
!
!	Profile data structure
!	INTEGER  ISTAT
!	RECORD /RTPPROF/ PROF            ! profile
!	RECORD /RTPHEAD/ HEAD            ! header data
!	RECORD /RTPATTR/ HATT(MAXNATTR)  ! header attributes
!	RECORD /RTPATTR/ PATT(MAXNATTR)  ! profile attributes
!
!	Boundary pressure levels
!	COMMON /COMLEV/ PLEV
	REAL :: PLEV(MAXLAY+1) = &
			(/  0.0050,    0.0161,    0.0384,    0.0769,    0.1370, &
				0.2244,    0.3454,    0.5064,    0.7140,    0.9753, &
				1.2972,    1.6872,    2.1526,    2.7009,    3.3398, &
				4.0770,    4.9204,    5.8776,    6.9567,    8.1655, &
				9.5119,   11.0038,   12.6492,   14.4559,   16.4318, &
			   18.5847,   20.9224,   23.4526,   26.1829,   29.1210, &
			   32.2744,   35.6505,   39.2566,   43.1001,   47.1882, &
			   51.5278,   56.1260,   60.9895,   66.1253,   71.5398, &
			   77.2396,   83.2310,   89.5204,   96.1138,  103.0172, &
			  110.2366,  117.7775,  125.6456,  133.8462,  142.3848, &
			  151.2664,  160.4959,  170.0784,  180.0183,  190.3203, &
			  200.9887,  212.0277,  223.4415,  235.2338,  247.4085, &
			  259.9691,  272.9191,  286.2617,  300.0000,  314.1369, &
			  328.6753,  343.6176,  358.9665,  374.7241,  390.8926, &
			  407.4738,  424.4698,  441.8819,  459.7118,  477.9607, &
			  496.6298,  515.7200,  535.2322,  555.1669,  575.5248, &
			  596.3062,  617.5112,  639.1398,  661.1920,  683.6673, &
			  706.5654,  729.8857,  753.6275,  777.7897,  802.3714, &
			  827.3713,  852.7880,  878.6201,  904.8659,  931.5236, &
			  958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, &
			 1100.0000 /)

	! channel selection of last initialization of AIRS_SARTA_VARIABLES
	INTEGER, DIMENSION(:), ALLOCATABLE :: LAST_SELCHAN
	! number of channels in last initialization of AIRS_SARTA_VARIABLES
	INTEGER :: LAST_NCHAN

CONTAINS

	FUNCTION AIRS_SARTA_VARIABLES__HAS_SELCHAN_CHANGED(NCHAN, SELCHAN) RESULT (SELCHAN_CHANGED)
		INTEGER, INTENT(IN) :: NCHAN ! number of selected channels ! FWI 12/15/17
		INTEGER, DIMENSION(MXCHAN), INTENT(IN) :: SELCHAN ! selected channel numbers (unit offset) ! FWI 12/15/17
		LOGICAL :: SELCHAN_CHANGED
 
		INTEGER :: ICHAN

		SELCHAN_CHANGED = .FALSE.

		! First check if there's a LAST_SELCHAN array at all
		IF (.NOT. ALLOCATED(LAST_SELCHAN)) THEN
			SELCHAN_CHANGED = .TRUE.
			RETURN
		ENDIF

		! Check NCHAN against LAST_NCHAN. If they're different, then the channel selection has changed
		IF (NCHAN .NE. LAST_NCHAN) THEN
			SELCHAN_CHANGED = .TRUE.
			RETURN
		ENDIF

		! Now do an element-by-element comparison of SELCHAN and LAST_SELCHAN. IF any are
		! different, then SELCHAN has changed.

		DO ICHAN = 1, NCHAN
			IF ((LAST_SELCHAN(ICHAN) - SELCHAN(ICHAN)) .NE. 0) THEN
				SELCHAN_CHANGED = .TRUE.
				EXIT
			ENDIF
		ENDDO

		RETURN

	END FUNCTION AIRS_SARTA_VARIABLES__HAS_SELCHAN_CHANGED

	SUBROUTINE AIRS_SARTA_VARIABLES__INITIALIZE(NCHAN, SELCHAN)

		INTEGER, INTENT(IN) :: NCHAN ! number of selected channels ! FWI 12/15/17
		INTEGER, DIMENSION(MXCHAN), INTENT(IN) :: SELCHAN ! selected channel numbers (unit offset) ! FWI 12/15/17

		INTEGER :: I
		INTEGER :: L	
		!INTEGER :: NCHAN

		SARTA_INITIALIZED = .FALSE.

		! commented out by FWI 12/15/17
		! Initialize INDCHN (channel indices of all channels)
		!DO I = 1, MXCHAN
		!	INDCHN(I) = I
		!ENDDO
		
		INDCHN(:) = 0
		DO I = 1, NCHAN
			INDCHN(SELCHAN(I)) = SELCHAN(I)
			LSTCHN(I) = SELCHAN(I)
		ENDDO

		! Mean layer pressure (KLAYERS definition)
		NLAY = MAXLAY
		DO L=1,MAXLAY
			PLAY(L) = ( PLEV(L+1) - PLEV(L) )/LOG( PLEV(L+1)/PLEV(L) )
		ENDDO
		ioun = 11

		ICO2 = 0			! No CO2 variability for now    SYLee Nov 2, 2007
		ICO2 = 1			! Have CO2 variability    ! FWI 12/19/17
		!CO2PPM = CO2STD ! FWI
		! SALT = XSALT             ! Nominal satelite height for now   SYLee Nov 2, 2007
		!-----------------------------
		! Read in the reference profile
		!-----------------------------
		!PRINT *, 'Reading reference profile ', FNPREF
		!CALL RDPROF(IOUN, FNPREF, RPNAM, RALT, RDZ, RPRES, RTEMP, &
		!	RFAMNT, RWAMNT, ROAMNT, RCAMNT, RMAMNT, RSAMNT,RHAMNT,RNAMNT)
		! PTYPE = 1


		! Read in and save the reference profiles
		CALL REF_PROFILES__READ_REF_PROFILES


		!
		! Get cloud table filenames
		!
		CALL FNMIE ( VCLOUD, MIETYP, FNMIEA, FNMIEE, FNMIEG )

		PRINT *, "VCLOUD ", VCLOUD
		PRINT *, "MIETYP ", MIETYP
		PRINT *, "FNMIEA ", FNMIEA

		!
		! Read cloud lookup tables
		!
		CALL RDCLDT( IOUN, INDCHN, MIETYP, FNMIEA, FNMIEE, FNMIEG, & 
			MIENPS, MIEPS, MIEABS, MIEEXT, MIEASY )       

		!
		! Set the variables set by opnrtp, which we are not calling
		!

		! Commented out by FWI 12/15/17
		!NCHAN = MXCHAN
		!DO L = 1, NCHAN
		!	INDCHN(L) = L		! Use all the channels
		!	LSTCHN(L) = L		! Use all the channels
		!END DO
		LRHOT = .true.
		IH2O = 0		! Temporarily use reference profile
		IO3 = 0			! Temporarily use reference profile
		ICO = 0			! Temporarily use reference profile
		ICH4 = 0		! Temporarily use reference profile
		ICO2 = 0		! Temporarily use reference profile
		ISO2 = 0		! Temporarily use reference profile
		IHNO3 = 0		! Temporarily use reference profile
		IN2O = 0		! Temporarily use reference profile


		!
		! Read the coef data files
		!
		PRINT *, 'Reading RTA coefficients'
		CALL RDCOEF( IOUN,   NCHAN,  INDCHN, SETCHN, &
			NCHN1,   NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7, &
			CLIST1,  CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
			COEF1,   COEF2,  COEF3,  COEF4,  COEF5,  COEF6,  COEF7, &
			FREQ,    LABOVE, COEFF,  INDCO2, COFCO2, INDSO2, COFSO2, &
			INDHNO,  COFHNO, INDN2O, COFN2O, &
			INDH2O,  WAZOP,  WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN )

		! Get and apply multipler tuning to coefficients {note: ignores HNO3}

		! Commented out by FWI for test  6/23/2016

		IF (.TRUE.) THEN
			CALL TUNMLT( IOUN, NCHAN, INDCHN, SETCHN, &
				NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7, &
			   CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
				COEF1,  COEF2,  COEF3,  COEF4,  COEF5,  COEF6,  COEF7, &
				 FREQ, LABOVE,  COEFF, INDCO2, COFCO2, INDSO2, COFSO2, &
			   INDHNO, COFHNO, INDN2O, COFN2O, &
			   INDH2O,  WAZOP, WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN )			   
		ELSE
			WRITE(*, *) "#### WARNING #### -- no tuning in SARTA"
		ENDIF


		! Calc OPTRAN absorption coefficient scaling factor WAOP
		WAOP(1)=WAZOP(1)
		DO L=2,MXOWLY
			WAOP(L)=WAZOP(L) - WAZOP(L-1)
		ENDDO

		!
		! Read in the solar radiance
		!
		CALL RDSUN(IOUN, INDCHN, HSUN)  ! added INDCHN ! FWI 12/19/17

		DISTES=1.496E+11  ! distance Earth to Sun
		!		--------------------
		!		Check FREQ and FCHAN
		!		--------------------
		!      Note: FREQ comes the coef data, while FCHAN comes
		!      from the input RTP file read by OPNRTP.  It is possible
		!      that FCHAN is "nodata", so we check the first element.
		!      IF (FCHAN(1) .GT. 640 .AND. FCHAN(1) .LT. 2670) THEN
		!         DO I=1,NCHAN
		!            RJUNK1=ABS(FREQ(I) - FCHAN(I))
		!            RJUNK2=0.01*FREQ(I)/1200.0   ! ~1% of a channel fullwidth
		!            IF (RJUNK1 .GT. RJUNK2) THEN
		!               WRITE(IOINFO,1010) I, LSTCHN(I), FREQ(I), FCHAN(I)
		!1010           FORMAT('Warning! index=',I4,', chan ID=',I4, &
		!               ', fastmodel freq=',F8.3,', RTP freq=',F8.3)
		!            ENDIF
		!            HEAD.vchan(I)=FREQ(I)
		!         ENDDO
		!      ELSE
		!         DO I=1,NCHAN
		!            HEAD.vchan(I)=FREQ(I)
		!         ENDDO
		!      ENDIF


		! We are not reading RTP file.   So copy frequency.   Nov 2, 2007    SYLee
		FCHAN(:) = FREQ(:)

		!      Precalculate Plank function related constants.  This will save CPU time.
		!      Sung-Yung Lee,    November 2007
		!       DO  I=1,MXCHAN
		!          C1V3(I)=C1*(FREQ(I)**3)
		!          C2V(I)=C2*FREQ(I)
		!       ENDDO

		!      -----------------------------------------------
		!      All channels from sets 1, 2, and 3 are to use a
		!      fake effective sun angle layer-to-space trans
		!      -----------------------------------------------

		NFAKE=0

		DO I=1,NCHN1
			NFAKE=NFAKE + 1
			INDFAK(NFAKE)=CLIST1(I)
		ENDDO

		DO I=1,NCHN2
			NFAKE=NFAKE + 1
			INDFAK(NFAKE)=CLIST2(I)
		ENDDO

		DO I=1,NCHN3
			NFAKE=NFAKE + 1
			INDFAK(NFAKE)=CLIST3(I)
		ENDDO

		LAST_NCHAN = NCHAN  ! number of channels at last initialization of AIRS_SARTA_VARIABLES
		IF (ALLOCATED(LAST_SELCHAN)) DEALLOCATE(LAST_SELCHAN)

		! Save channel selection for comparison to later calls
		ALLOCATE(LAST_SELCHAN(NCHAN))
		LAST_SELCHAN(1:NCHAN) = SELCHAN(1:NCHAN)


		SARTA_INITIALIZED = .TRUE.

	END SUBROUTINE AIRS_SARTA_VARIABLES__INITIALIZE

END MODULE AIRS_SARTA_VARIABLES

