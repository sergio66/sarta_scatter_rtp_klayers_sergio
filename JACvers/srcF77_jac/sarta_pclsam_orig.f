c This version with rtp input profiles and command-line arguments
C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County (UMBC)
C
C    AIRS
C
C    SARTA_pclsam
C
!F77====================================================================


!ROUTINE NAME:
C    SARTA_pclsam


!ABSTRACT:
C    Program to quickly compute simulated AIRS radiances.
C
C    This variant of the code allows for the modelling of a
C    up to two scatter clouds using the PCLSAM method.


!CALL PROTOCOL
C    none (main program)


!INPUT PARAMETERS:
C    none


!OUTPUT PARAMETERS:
C    none


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S)
C    none


!ROUTINES CALLED:
C    CALOWP  : calc OPTRAN water predictors 
C    CALPAR  : calculate a profile's predictors
C    CCPREP  : prepare lookup table for given cp-size.
C    BKPREP  : check cloud layer
C    GETMIE  : determine which scattering lookup table to use
C    GETCLD  : get basic cloud parameters
C    CALRAD0 : calc radiance
C    CALRAD1 : calc radiance
C    CALT1   : calc effective layer trans for set1 (FWO)
C    CALT2   : calc effective layer trans for set2 (FOW)
C    CALT3   : calc effective layer trans for set3 (FMW)
C    CALT4   : calc effective layer trans for set4 (FCOW)
C    CALT5   : calc effective layer trans for set5 (FWO bfsw)
C    CALT6   : calc effective layer trans for set6 (FWO mfmw)
C    CALT7   : calc effective layer trans for set7 (FWO mfbw)
C    FAKETZ  : calc a "fake" (rough approx) surface-to-space trans
C    RDCOEF  : read the fast transmittance coefficients
C    RDPROF  : read profile data
C    RDSUN   : read the solar radiance datafile
C    SUNPAR  : calc a profile's predictors for sets4-7 (sun channels)
C    CALNTE  : calc radiance contribution for non-LTE
C    OPNRTP  : open the RTP file
C    RDRTP   : read the RTP file
C    SETEMS  : sets surface parameters

!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.
C    unit IOUN: used by routines RDCOEF and RDPROF.
C    unit 6: USEFAST text messages to the screen
C    unit 5: USEFAST user input instructions, etc
C    unit 10: USEFAST output radiance, text file(s)


!COMMON BLOCKS
C    COMLEV : layer boundary pressure levels


!DESCRIPTION:
C    Dec 2005 version of the SARTA (Stand-Alone Rapid
C    Transmittance Algorith with RTP I/O) by
C    L.L.Strow, S.Hannon, and H.Mottler
C
C    Computes radiances for the layers profiles contained in the
C    input RTP file.  This is the main program, and consists
C    primarily of calls to the external routines to do most of the
C    computing.


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    This program is only intended as a demo of the fast model.


!ROUTINE HISTORY:
C Date         Programmer    Comments
C ----------- -------------- -------------------------------------------
C 01 Dec 1994 Scott Hannon   Created
C 10 Apr 1995 Scott Hannon   New header comments; added ALT; new
C                            external function VACONV; SECANG may
C                            vary with layer
C 03 Jul 1995 Scott Hannon   Add parameter DZ/RDZ to RDPROF call
C 03 Feb 1997 Scott Hannon   Re-written for FWO+FOW+FMW+FCOW
C 12 Sep 1997 Scott Hannon   Re-written for 7 sets and reflected sun
C                            and downwelling thermal
C 30 Sep 1997 Scott Hannon   Added variable CO2
C 27 Feb 1998 Scott Hannon   Added OPTRAN water
C 26 Aug 1998 Scott Hannon   Added LBOT to calls to CALPAR, CALOWP,
C                            and SUNPAR; rename TBOT to TSURF; calc
C                            fractional bottom layer temperature and
C                            put it in TEMP(LBOT)
C 15 Oct 1999 Scott Hannon   Add ANGMAX and re-arrange angle conv
C 31 Mar 2000 Scott Hannon   Redid calpar for FIXMUL and added getbot
C 15 Mar 2001 Scott Hannon   Major re-write for RTP
C 03 May 2001 Scott Hannon   Add COMLEV; add PLEV to getbot call
C 13 Sep 2001 Scott Hannon   Changes to check of FCHAN vs FREQ
C 01 Nov 2002 Scott Hannon   Added SATZEN & SALT to RDRTP call, and
C                            if valid use SATZEN rather than SATANG
C 03 Jan 2003 Scott Hannon   Delete SUNSEC, add XZ & SUNFDG & code
C                            to fudge large sun angles (previously
C                            sunang>80 were treated as no sun).
C 24 Jul 2003 Scott Hannon   Fix error in TEMP(LBOT) calc for
C                            bottom fractional layer; add PLAY
C 06 Feb 2004 Scott Hannon   Add call to TUNMLT; add call to MEAN_T
C                            and associated prep code; add PTYPE
C                            to OPNRTP call; add LRHOT to RDINFO,
C                            OPNRTP, & SETEMS calls.
C 20 Dec 2004 Scott Hanonn   Add NLAY to getbot.f call; add PTYPE
C                            to rdrtp_so2.f call
C 18 May 2005 Scott Hannon   Add HNO3 based on SO2 code
C 28 Jun 2005 Scott Hannon   "trace" version for CO2,SO2,HNO3,N2O
C 13 Oct 2005 Scott Hannon   Add non-LTE
C 08 Dec 2005 Scott Hannon   Update tunmlt call for non-LTE tuning
C 29 Mar 2006 Scott Hannon   Add clouds; change TAU from trans to od
C 26 Apr 2006 Scott Hannon   Add black clouds. Redo cloud fractions.
C                            Use RTP v1.06 cloud2 fields instead
C                            of udef(11-17). {Unfinished}
C 22 Dec 2006 Scott Hannon   New & revised code associated with the
C                            new & revised calrad* routines for more
C                            fexible cloud types including black
C                            clouds; changes to calt* calls and
C                            change TAUZ from (1 x n) to (m x n).
C 22 Jan 2007 Scott Hannon   Minor fix of cfrac checks
C 02 May 2007 Scott Hannon   Replace hardcoded default SALT value
C                            with XSALT from incFTC.
C 15 Nov 2007 Scott Hannon   Move most cloud prep to GETCLD & BKPREP.
C 31 Jan 2008 Scott Hannon   Add LCO2PM to allow CO2 profile in ppmv;
C                               add LCO2,LN2O,LSO2,LHNO3 switches
C 24 Mar 2008 Scott Hannon   Add COSDAZ
C 26 Nov 2008 Scott Hannon   Update for rtpV2101
C 01 Dec 2008 Scott Hannon   Add CSTMP1/2
C    Jul 2019 C Hepplewhite  Add NH3 and align with sarta clear updates

!END====================================================================

C      =================================================================
       PROGRAM SARTA
C      =================================================================


C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
       include 'incFTC.f'
       include 'rtpdefs.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
       REAL VACONV
       REAL SACONV


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      none (main program)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
C
       INTEGER   IOUN         ! I/O unit number
C
C      for RDINFO
       CHARACTER*120 FIN       ! input RTP filename
       CHARACTER*120 FOUT      ! output RTP filename
       LOGICAL  LRHOT         ! force refl therm rho=(1-emis)/pi?
       INTEGER NWANTP         ! number of wanted profiles (-1=all)
       INTEGER  LISTP(MAXPRO) ! list of wanted profiles
       INTEGER NWANTC         ! number of wanted channels (-1=all)
       INTEGER  LISTC(MAXPRO) ! list of wanted channels

C
C      for FNMIE
       CHARACTER*240 VCLOUD        ! cloud version string
       INTEGER MIETYP(NMIETY)      ! mie type
       CHARACTER*79 FNMIEA(NMIETY) ! mie absorption filenames
       CHARACTER*79 FNMIEE(NMIETY) ! mie extinction filenames
       CHARACTER*79 FNMIEG(NMIETY) ! mie asymmetry filenames
C
C      for OPNRTP
       INTEGER  PTYPE          ! profile type
       INTEGER  NCHAN          ! # of selected channels
       REAL     FCHAN(MXCHAN)  ! chan center frequency
       INTEGER LSTCHN(MXCHAN)  ! list of selected channels
       INTEGER INDCHN(MXCHAN)  ! array indices for all channels
       INTEGER RINDCHN(MXCHAN) ! list of locations of chans eg chID1291 = 1231 cm-1 but this is location 1520 in L1C
       INTEGER IH2O           ! index of H2O in gamnt
       INTEGER IO3            ! index of O3 in gamnt
       INTEGER ICO            ! index of CO in gamnt
       INTEGER ICH4           ! index of CH4 in gamnt
       INTEGER ICO2           ! index of CO2 in gamnt
       INTEGER ISO2           ! index of SO2 in gamnt
       INTEGER IHNO3          ! index of HNO3 in gamnt
       INTEGER IN2O           ! index of N2O in gamnt
       INTEGER INH3           ! index of NH3 in gamnt
       INTEGER IOPCI          ! input RTP unit
       INTEGER IOPCO          ! output RTP unit
       LOGICAL LCO2PM         ! CO2 profile in ppmv?
C
C      for RDCLDT
       INTEGER MIENPS(NMIETY)            ! number of particle sizes
       REAL  MIEPS(MXMIEA,NMIETY)        ! Mie particle size for table
       REAL MIEABS(MXCHAN,MXMIEA,NMIETY) ! Mie absorption table
       REAL MIEEXT(MXCHAN,MXMIEA,NMIETY) ! Mie extinction table
       REAL MIEASY(MXCHAN,MXMIEA,NMIETY) ! Mie asymmetry table
C
C      for RDCOEF             ! Info for selected channels only
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

       !!! this is to intersect the lists ie replace 
       !!! do I = 1 : NCHAN
       !!!   III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4)       
       !!!   IF III > 0
       !!!     call YCALT4(INDCHN,   LBOT,  NCHN4, CLIST4,...,III)
       !!!   END IF
       !!! end do
       INTEGER QUICKCLIST1(MXCHAN) ! list of set1 channels
       INTEGER QUICKCLIST2(MXCHAN) ! list of set2 channels
       INTEGER QUICKCLIST3(MXCHAN) ! list of set3 channels
       INTEGER QUICKCLIST4(MXCHAN) ! list of set4 channels
       INTEGER QUICKCLIST5(MXCHAN) ! list of set5 channels
       INTEGER QUICKCLIST6(MXCHAN) ! list of set6 channels
       INTEGER QUICKCLIST7(MXCHAN) ! list of set7 channels

       INTEGER LABOVE(MXCHAN) ! chan downwelling thermal layer above
       REAL   FREQ(MXCHAN)    ! chan center frequency
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
       INTEGER INDHDO(MXCHAN)            ! chan indices for HDO pert
       REAL COFHDO(  NHDO,MAXLAY,MXCHND) ! coefs for HDO pert
       INTEGER INDHNO(MXCHAN)            ! chan indices for HNO3 pert
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH) ! coefs for HNO3 pert
       INTEGER INDN2O(MXCHAN)            ! chan indices for N2O pert
       REAL COFN2O(  NN2O,MAXLAY,MXCHNN) ! coefs for N2O pert
       INTEGER INDNH3(MXCHAN)            ! chan indices for NH3 pert
       REAL COFNH3(  NNH3,MAXLAY,MXCHNA) ! coefs for NH3 pert
       INTEGER INDH2O(MXCHAN)            ! chan indices for OPTRAN H2O
       REAL   WAZOP(MXOWLY)              ! OPTRAN water l-to-s amounts
       REAL  WAVGOP(NOWAVG,MXOWLY)       ! OPTRAN raw predictor averages
       REAL COFH2O(  NH2O,MXOWLY,MXCHNW) ! coefs for OPTRAN H2O
       REAL     FX(MAXLAY)               ! fixed gases adjustment
       INTEGER NCHNTE                    ! number of non-LTE channels
       INTEGER CLISTN(MXCNTE)            ! non-LTE channel list
       REAL  COEFN(NNCOEF,MXCNTE)        ! non-LTE coefficients
C
C      for rtpopen
       INTEGER rtpopen
       CHARACTER*1 MODE
C
C      for FAKETZ
       INTEGER  NFAKE              ! # of channels to "fake"
       INTEGER INDFAK(MXCHAN)      ! indices of channels to fake
       INTEGER QUICKINDFAK(MXCHAN) ! list of set1 channels
C
C      for RDPROF; reference profile
       CHARACTER*40  RPNAM ! ref prof name/ID
       REAL   RALT(MAXLAY) ! ref prof layer altitude
       REAL    RDZ(MAXLAY) ! ref prof layer thickness
       REAL  RPRES(MAXLAY) ! ref prof layer average pressure
       REAL  RTEMP(MAXLAY) ! ref prof layer average temperature
       REAL RFAMNT(MAXLAY) ! ref prof layer "fixed" (CO2) amount
       REAL RWAMNT(MAXLAY) ! ref prof layer water (H2O) amount
       REAL ROAMNT(MAXLAY) ! ref prof layer ozone (O3) amount
       REAL RCAMNT(MAXLAY) ! ref prof layer carbon monoxide (CO) amount
       REAL RMAMNT(MAXLAY) ! ref prof layer methane (CH4) amount
       REAL RSAMNT(MAXLAY) ! ref prof layer sulfer dioxide (SO2) amount
       REAL RHAMNT(MAXLAY) ! ref prof layer nitric acid (HNO3) amount
       REAL RNAMNT(MAXLAY) ! ref prof layer nitrous oxide (N2O) amount
       REAL RAAMNT(MAXLAY) ! ref prof layer ammonia (NH3) amount
C
C      for RDRTP; profile to calculate
       INTEGER NLAY        ! number of layers in profile
       REAL LAT            ! prof latitude
       REAL LON            ! prof longitude
       REAL    ALT(MAXLAY) ! prof layer altitudes
       REAL   TEMP(MAXLAY) ! prof layer average temperature
       REAL  WAMNT(MAXLAY) ! prof layer water (H2O) amount
       REAL  OAMNT(MAXLAY) ! prof layer ozone (O3) amount
       REAL  CAMNT(MAXLAY) ! prof layer carbon monoxide (CO) amount
       REAL  MAMNT(MAXLAY) ! prof layer methane (CH4) amount
       REAL  FAMNT(MAXLAY) ! prof layer CO2 amount
       REAL  SAMNT(MAXLAY) ! prof layer SO2 amount
       REAL  HAMNT(MAXLAY) ! prof layer HNO3 amount
       REAL  NAMNT(MAXLAY) ! prof layer N2O amount
       REAL  AAMNT(MAXLAY) ! prof layer NH3 amount
C
C      for surface
       INTEGER   LBOT             ! bottom layer index number
       INTEGER  NEMIS             ! # of emis pts
       REAL  PSURF                ! surface pressure
       REAL BLMULT                ! bottom layer fractional multiplier
       REAL  FEMIS(MXEMIS)        ! emis freq pts
       REAL  XEMIS(MXEMIS)        ! emis pts
       REAL   XRHO(MXEMIS)        ! reflec pts
C
C      for MEAN_T
       REAL TPSEUD(MAXLAY)
C
C      for CALPAR
       LOGICAL   LCO2             ! CO2 profile switch
       LOGICAL   LN2O             ! N2O profile switch
       LOGICAL   LSO2             ! SO2 profile switch
       LOGICAL  LHNO3             ! HNO3 profile switch
       LOGICAL   LNH3             ! NH3 profile switch
       LOGICAL   LHDO             ! HDO profile switch
       REAL SECANG(MAXLAY)        ! local path angle secant
       REAL FIXMUL(MAXLAY)        ! "fixed" amount multiplier (~1)

       REAL CONPRD( N1CON,MAXLAY) ! water continuum predictors
       REAL SUNCONPRD( N1CON,MAXLAY) ! water continuum predictors

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

       REAL  DPRED(  NHDO,MAXLAY) ! HDO perturbation predictors

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

       REAL MPRED3( N3CH4,MAXLAY) ! set3 methane predictors

       REAL CPRED4(  N4CO,MAXLAY) ! set4 carbon monoxide predictors
       REAL SUNCPRED4(  N4CO,MAXLAY) ! set4 carbon monoxide predictors
       REAL TRCPRD(NTRACE,MAXLAY) ! trace gas pert perdictors
       REAL SUNTRCPRD(NTRACE,MAXLAY) ! trace gas pert perdictors

       REAL CO2MLT(MAXLAY)        ! CO2 perturbation multiplier
       REAL SO2MLT(MAXLAY)        ! SO2 perturbation multiplier
       REAL HNOMLT(MAXLAY)        ! HNO3 perturbation multiplier
       REAL N2OMLT(MAXLAY)        ! N2O perturbation multiplier
       REAL NH3MLT(MAXLAY)        ! NH3 perturbation multiplier
       REAL HDOMLT(MAXLAY)        ! HDO perturbation multiplier
       REAL CO2TOP                ! top layers CO2 mixing ratio
C
C      for CALOWP
       REAL  WAANG(MAXLAY)
       INTEGER LOPMIN
       INTEGER LOPMAX
       REAL H2OPRD(  NH2O,MXOWLY)
       LOGICAL LOPUSE(MXOWLY)
       INTEGER LOPLOW(MAXLAY)
       REAL  DAOP(MAXLAY)
C
C      for CALT
       REAL    TAU(MAXLAY,MXCHAN) ! chan layer effective optical depth
       REAL   TAUZ(MAXLAY,MXCHAN) ! chan surface-to-space trans
       REAL   WAOP(MXOWLY)        ! OPTRAN abs coef scaling factor
       REAL    XZ                 ! Optical depth multiplier for TAUZ
       LOGICAL  LTAU              ! calc all layer transmittances?
C
C      for SETEMS
       REAL   EMIS(MXCHAN) ! chan surface emissivity
       REAL RHOSUN(MXCHAN) ! chan reflectivity for sun
       REAL RHOTHR(MXCHAN) ! chan reflectivity for downwelling thermal
       REAL CEMIS1(MXCHAN) ! chan surface emissivity cloud1
       REAL CRHOS1(MXCHAN) ! chan solar reflectivity cloud1
       REAL CRHOT1(MXCHAN) ! chan thermal reflectivity cloud1
       REAL CEMIS2(MXCHAN) ! chan surface emissivity cloud2
       REAL CRHOS2(MXCHAN) ! chan solar reflectivity cloud2
       REAL CRHOT2(MXCHAN) ! chan thermal reflectivity cloud2
C
C      for CALRAD(0,1)
       REAL SUNFAC         ! sun solid angles times cosine at surface
       REAL RPLNCK(MAXLAY) ! layer Planck
       REAL RSURFE         ! surface emission
       REAL RSURFC         ! black cloud surface emission
       REAL  TRANL(MAXLAY) ! clear air layer transmittance
       REAL  TRANZ(MXCHAN) ! clear air layer-to-space transmittance
       REAL  TRANS(MXCHAN) ! clear air total reflected solar trans
       REAL  TSURF         ! surface temperature
       REAL    RAD(MXCHAN) ! chan radiance
C      For clear/cloudy radiances
       REAL   RAD0         ! radiance no clouds
       REAL  RADC1         ! radiance cloud1
       REAL  RADC2         ! radiance cloud2
       REAL RADC12         ! radiance cloud1+cloud2
C
C      for RDSUN
       REAL   HSUN(MXCHAN) ! sun radiance (direct from sun)
C      Other variables for the sun
       REAL SUNANG         ! solar zenith angle (at 0 altitude)
       REAL COSDAZ         ! cosine(solazi - satazi) {COS Delta AZimuth}
       REAL SZALAY         ! solar zenith angle in some layer
       REAL SUNCOS         ! cosine of sun zenith angle
       REAL SCOS1          ! cosine of sun zenith angle at layer1
       REAL SUNFDG         ! fudge factor for large solar angles
       REAL SECSUN(MAXLAY) ! secant of effective sun local path angle
       REAL DISTES         ! distance of Earth from the sun
       REAL TAUZSN(MAXLAY,MXCHAN) ! sun space-to-surface-to-space OD
       LOGICAL DOSUN       ! do sun calc?
C
C      for satellite viewing angle
       REAL    SATANG      ! input satellite scan angle (degrees)
       REAL    SATZEN      ! input satellite zenith angle (degrees)
       REAL    SALT        ! input satellite altitude (kilometers)
       REAL    SVA         ! satellite viewing angle (degrees)
C
C      for RDRTP
       INTEGER  IPROF      ! profile loop counter
       LOGICAL  LWANT      ! do you want this profile?
C
C      Basic cloud info
       REAL XCEMI1(MXEMIS)    ! cloud1 emissivity
       REAL XCEMI2(MXEMIS)    ! cloud2 emissivity
       REAL XCRHO1(MXEMIS)    ! cloud1 reflectivity
       REAL XCRHO2(MXEMIS)    ! cloud2 reflectivity
       REAL CFRAC1            ! cloud1(total) fraction of FOV
       REAL CFRAC2            ! cloud2(total) fraction of FOV
       REAL CFRA1X            ! cloud1(exclusively) fraction of FOV
       REAL CFRA2X            ! cloud2(exclusively) fraction of FOV
       REAL CFRA12            ! cloud1+2(both) fraction of FOV
       REAL CNGWA1            ! cloud1 non-gases water
       REAL CNGWA2            ! cloud1 non-gases water
       REAL CPRBO1            ! cloud1 bottom pressure
       REAL CPRBO2            ! cloud2 bottom pressure
       REAL CPRTO1            ! cloud1 top pressure
       REAL CPRTO2            ! cloud2 top pressure
       REAL CPSIZ1            ! cloud1 particle size
       REAL CPSIZ2            ! cloud2 particle size
       REAL CSTMP1            ! cloud1 top/surf temperature
       REAL CSTMP2            ! cloud2 top/surf temperature
       REAL FCLEAR            ! clear (no cloud) fraction of FOV
       REAL TEMPC1            ! cloud1 frac layer (above cloud) mean temp
       REAL TEMPC2            ! cloud2 frac layer (above cloud) mean temp
       INTEGER CTYPE1         ! cloud1 type code number
       INTEGER CTYPE2         ! cloud2 type code number
C
C      for GETMIE
       LOGICAL LBLAC1  ! black cloud1? {Mie cloud if false}
       LOGICAL LBLAC2  ! black cloud2? {Mie cloud if false}
       INTEGER INDMI1  ! index in MIETYP for CTYPE1
       INTEGER INDMI2  ! index in MIETYP for CTYPE2
       INTEGER  IERR1  ! error level of CTYPE1/MIETYP match
       INTEGER  IERR2  ! error level of CTYPE2/MIETYP match
C
C      for CCPREP cloud1
       INTEGER LCBOT1         ! layer containing cloud bottom
       INTEGER LCTOP1         ! layer containing cloud top
       REAL  CLRB1            ! frac of layer at bottom of cloud clear
       REAL  CLRT1            ! frac of layer at top of cloud clear
       REAL TCBOT1            ! temperature at cloud bottom
       REAL TCTOP1            ! temperature at cloud top
       REAL MASEC1            ! mean cloud view angle secant
       REAL MASUN1            ! mean cloud sun-only angle secant
       REAL CFRCL1(MAXLAY)    ! fraction of cloud in layer
       REAL G_ASY1(MXCHAN)    ! "g" asymmetry
       REAL NEXTO1(MXCHAN)    ! nadir extinction optical depth
       REAL NSCAO1(MXCHAN)    ! nadir scattering optical depth
C
C      for CCPREP cloud2
       INTEGER LCBOT2         ! layer containing cloud bottom
       INTEGER LCTOP2         ! layer containing cloud top
       REAL  CLRB2            ! frac of layer at bottom of cloud clear
       REAL  CLRT2            ! frac of layer at top of cloud clear
       REAL TCBOT2            ! temperature at cloud bottom
       REAL TCTOP2            ! temperature at cloud top
       REAL MASEC2            ! mean cloud view angle secant
       REAL MASUN2            ! mean cloud sun-only angle secant
       REAL CFRCL2(MAXLAY)    ! fraction of cloud in layer
       REAL G_ASY2(MXCHAN)    ! "g" asymmetry
       REAL NEXTO2(MXCHAN)    ! nadir extinction optical depth
       REAL NSCAO2(MXCHAN)    ! nadir scattering optical depth
C
C      used locally only
       INTEGER  I,II,III   ! loop counter
       INTEGER      L      ! loop counter
       INTEGER rtpclose    ! for call to RTP close interface routine
       REAL    EVA         ! (Earth) local viewing angle
       REAL   CONV         ! degrees to radians conversion factor
       REAL ANGMAX         ! maximum allowed viewing angle
       REAL RJUNK1         ! junk/work
       REAL RJUNK2         ! another junk/work
       REAL CO2PPM         ! Profile mean dry air CO2 mixing ratio
       REAL PLAY(MAXLAY)   ! layer mean pressure
       REAL C1V3           ! rad constant c1 times freq^3
       REAL C2V            ! rad constant c2 times freq
       REAL VSTORE(6)      ! temporary storage for various variables
C
C      Profile data structure
       INTEGER  ISTAT
       RECORD /RTPPROF/ PROF            ! profile
       RECORD /RTPHEAD/ HEAD            ! header data
       RECORD /RTPATTR/ HATT(MAXNATTR)  ! header attributes
       RECORD /RTPATTR/ PATT(MAXNATTR)  ! profile attributes
C
C      Boundary pressure levels
       COMMON /COMLEV/ PLEV
       REAL PLEV(MAXLAY+1)
C
C      for function QIKEXP
       REAL QIKEXP

C      for function intersect
       INTEGER intersect

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none
C
C***********************************************************************
C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
C      CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02

C      --------------------------
C      Assign the I/O unit number
C      --------------------------
       IOUN=11

C      ---------
C      Calc PLAY
C      ---------
C      Mean layer pressure (KLAYERS definition)
       DO L=1,MAXLAY
          PLAY(L) = ( PLEV(L+1) - PLEV(L) )/LOG( PLEV(L+1)/PLEV(L) )
       ENDDO

C      -----------------------------
C      Read in the reference profile
C      -----------------------------
       CALL RDPROF(IOUN, FNPREF, RPNAM, RALT, RDZ, RPRES, RTEMP,
     $    RFAMNT, RWAMNT, ROAMNT, RCAMNT, RMAMNT, RSAMNT,
     $    RHAMNT, RNAMNT, RAAMNT)

       if (DEBUG) then
          print*, 'sarta_cloudy: completed rdprof'
       endif
C      ---------------------
C      Get command-line info
C      ---------------------
       CALL RDINFO(FIN, FOUT, LRHOT, NWANTP, LISTP, NWANTC, LISTC)
ccc
       if (DEBUG) then
         print *, 'nwantp=', NWANTP
         print *, 'listp=', (LISTP(I),I=1,NWANTP)
         print *, 'FIN = ',FIN
         print *, 'FOUT = ',FOUT
       endif
ccc
       if (nwantp .gt. 0) then
         print *,nwantp
         print *,listp
       else
         print *,'nwantp = -1 (all profs)'
       end if

       if (nwantc .gt. 0) then
         print *,nwantc
         print *,listc
       else
         print *,'nwantc = -1 (all chans)'
       end if

C      -------------------------
C      Get cloud table filenames
C      -------------------------
       CALL FNMIE ( VCLOUD, MIETYP, FNMIEA, FNMIEE, FNMIEG )

C      ---------------------------
C      Open & check input RTP file
C      ---------------------------
       CALL OPNRTP(FIN, VCLOUD, LRHOT, PTYPE, NCHAN, FCHAN, LSTCHN,
     $    INDCHN, IH2O, IO3, ICO, ICH4, ICO2, ISO2, IHNO3, IN2O,
     $    INH3, IOPCI, HEAD, HATT, PATT, LCO2PM, NWANTC, LISTC, RINDCHN)
  
C       print*, 'sarta_cloudy: completed opnrtp'
ccc
C      ------------------------
C      Read cloud lookup tables
C      ------------------------
       CALL RDCLDT( IOUN, INDCHN, MIETYP, FNMIEA, FNMIEE, FNMIEG,
     $    MIENPS, MIEPS, MIEABS, MIEEXT, MIEASY )       

C      ------------------------
C      Read the coef data files
C      ------------------------
       CALL RDCOEF( IOUN, NCHAN, INDCHN, SETCHN,
     $  NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7,
     $ CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
     $  COEF1,  COEF2,  COEF3,  COEF4,  COEF5,  COEF6,  COEF7,
     $   FREQ, LABOVE,  COEFF, INDCO2, COFCO2, INDSO2, COFSO2,
     $ INDHNO, COFHNO, INDN2O, COFN2O, INDNH3, COFNH3, 
     $ INDHDO, COFHDO,
     $ INDH2O,  WAZOP, WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN)
C
C      Get and apply multipler tuning to coefficients {note: ignores HNO3}
       CALL TUNMLT( IOUN, NCHAN, INDCHN, SETCHN,
     $  NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7,
     $ CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
     $  COEF1,  COEF2,  COEF3,  COEF4,  COEF5,  COEF6,  COEF7,
     $   FREQ, LABOVE,  COEFF, INDCO2, COFCO2, INDSO2, COFSO2,
     $ INDHNO, COFHNO, INDN2O, COFN2O,
     $ INDH2O,  WAZOP, WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN )

       IF (NWANTC .GT. 0) THEN
         write(*,'(A)') '     ---------------------------------------------------------------------------------------------'
         write(*,'(A,I5,A)') 'after opnrtp, NCHAN = ',NCHAN,' .... here are the channels after opnrtp .... '
         write(*,'(A)') '              I          LSTCHN(I)      RINDCHN(I)    INDCHN(LSTCHN(I))  BREAKOUT     FCHAN(I) '
         write(*,'(A)') '                                                                     SETCHN(LSTCHN(I)          '
         write(*,'(A)') '     ---------------------------------------------------------------------------------------------'
         DO I = 1,NCHAN
           write(*,'(5(I15),F20.7)') I,LSTCHN(I),RINDCHN(I),INDCHN(LSTCHN(I)),SETCHN(LSTCHN(I)),FCHAN(I)
         END DO
         write(*,'(A)') '     ---------------------------------------------------------------------------------------------'
       END IF

      !!!! suppose this is AIRS and there are 2834 chans
      !!! see incFTC.f : so if you want all 2834 chans computed, then NCHN1 == MXCHN1 = 1461 for FWO, NCHN2 = MXCHN2 = 325 for FOW etc etc etc
      !!!      ie 1461   325     396     85    210    217    140     which sums to 2834
      !!!
      !!! but eg in my retrieval code CRODGERS_FAST_CLOUD, I have chose LW chnas only about 420 between 15 um, window, O3, WV
      !!! ie h.nchan = 426, h.ichan = >> h.ichan(1:15)' = 25 52 62 69 70 71 73 75 77 78 79 80 82 83 84 ...
      !!! then NCHNX will say of the ODs in setX, how many of these IDs should be computed!!!!!
      !!! but in the example from retrievals, this ends up being 277 56 69 67 7 0 0   which sums to 476     since no SW channels (set 6,7 = none chose) 
!      print *,NCHN1,NCHN2,NCHN3,NCHN4,NCHN5,NCHN6,NCHN7
      !!! and NCHAN witll b the sum of above == 476 = h.ichan
!      print *,NCHAN

      !!!! so there will be ZEROS where h.ichan does not exist, and the index where it does so for example this will look like
      !!!!                       0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           1           0           0           0           0     index 25
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           2           0     index 52
      !!!!           0           0           0           0           0           0
      !!!!           0           0           3           0           0           0     index 62
      !!!!           0           0           0           4           5           6     index 69 70 71
      !!!!           0           7           0           8           0           9     index 73 75 77
      !!!!          10          11          12           0          13          14     index 78 79 80 82 83 ... 

      QUICKCLIST1(INDCHN(CLIST1(1:NCHN1))) = (/(I,I=1,NCHN1)/)
      QUICKCLIST2(INDCHN(CLIST2(1:NCHN2))) = (/(I,I=1,NCHN2)/)
      QUICKCLIST3(INDCHN(CLIST3(1:NCHN3))) = (/(I,I=1,NCHN3)/)
      QUICKCLIST4(INDCHN(CLIST4(1:NCHN4))) = (/(I,I=1,NCHN4)/)
      QUICKCLIST5(INDCHN(CLIST5(1:NCHN5))) = (/(I,I=1,NCHN5)/)
      QUICKCLIST6(INDCHN(CLIST6(1:NCHN6))) = (/(I,I=1,NCHN6)/)
      QUICKCLIST7(INDCHN(CLIST7(1:NCHN7))) = (/(I,I=1,NCHN7)/)

       if (DEBUG) then
         print *,'NCHN4 = ',NCHN4
         print *,'SETCHN(CLIST4(1:NCHN4)) = ',SETCHN(CLIST4(1:NCHN4))
         print *,'CLIST4(1:NCHN4) = ',CLIST4(1:NCHN4)
         print *,'INDCHN(CLIST4(1:NCHN4)) = ',INDCHN(CLIST4(1:NCHN4))
         print *,'QLIST4(INDCHN(CLIST4(1:NCHN4))) = ',QUICKCLIST4(INDCHN(CLIST4(1:NCHN4)))
         do I = 1,NCHN4
           print *,I,CLIST4(I),INDCHN(CLIST4(I)),QUICKCLIST4(INDCHN(CLIST4(I)))
         end do
       end if

      if (DEBUG) then
        print *,NCHN1,NCHN2,NCHN3,NCHN4,NCHN5,NCHN6,NCHN7
        print *,NCHAN
        print *,'INDCHN(1:NCHAN) = ',INDCHN(1:NCHAN)
        DO I = 1,NCHAN
           !! print i,h.ichan(i),h.vchan(i)
           print *,I,LSTCHN(I),SETCHN(I),FREQ(I)
         END DO
         print *,'CLIST1 = ',CLIST1
         print *,'CLIST2 = ',CLIST2
!        print *,'CLIST3 = ',CLIST3
!        print *,'CLIST4 = ',CLIST4
!        print *,'CLIST5 = ',CLIST5
!        print *,'CLIST6 = ',CLIST6
!        print *,'CLIST7 = ',CLIST7
       END IF

       if (DEBUG) then
         print *,'CLIST2 = ',CLIST2(1:NCHN2)
         print *,'INDCHN(CLIST2) = ',INDCHN(CLIST2(1:NCHN2))
         DO I = 1,NCHAN
           III = intersect(I,INDCHN(CLIST2(1:NCHN2)), NCHN2)
           IF (III .GT. 0) print *,I,III,INDCHN(CLIST2(III))
         END DO
         STOP
       end if

      !!! of the above INDCHN, which match up with set OD1
      !!! CLIST1 = 25 52 62 69 70 71 73 75 77 78 79 80 82 83 84 .... 1785
      !!! CLIST2 = 149 150 196 203 206 207 217 220 221 223 224 ... 1221 1236

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    SO print *,INDCHN(CLIST1) should be INDCHN([25 52 62 69 70 71 73 75 77 ...]) = 1 2 3 4 5 6 ....
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       print *,INDCHN(CLIST2(1:NCHN2))
!         63  64  110 117 120 121
!         131 134 135 137 138 139
!         142 143 144 145 146 147
!         148 149 150 164 166 167
!         168 169 170 172 179 180
!         242 243 244 245 246 247
!         248 249 250 251 252 253
!         254 255 256 257 259 260
!         261 262 263 264 265 266
!         267 268            

C      Calc OPTRAN absorption coefficient scaling factor WAOP
       WAOP(1)=WAZOP(1)
       DO L=2,MXOWLY
          WAOP(L)=WAZOP(L) - WAZOP(L-1)
       ENDDO

C      --------------------------
C      Read in the solar radiance
C      --------------------------
       CALL RDSUN(IOUN, INDCHN, HSUN)
C
       DISTES=1.496E+11  ! distance Earth to Sun

C      --------------------
C      Check FREQ and FCHAN
C      --------------------
C      Note: FREQ comes the coef data, while FCHAN comes
C      from the input RTP file read by OPNRTP.  It is possible
C      that FCHAN is "nodata", so we check the first element.
       IF (FCHAN(1) .GT. 640.0 .AND. FCHAN(1) .LT. 2670.0) THEN
          DO I=1,NCHAN
             RJUNK1=ABS(FREQ(I) - FCHAN(I))
             RJUNK2=0.01*FREQ(I)/1200.0   ! ~1% of a channel fullwidth
             IF (RJUNK1 .GT. RJUNK2) THEN
                WRITE(IOINFO,1010) I, LSTCHN(I), FREQ(I), FCHAN(I)
 1010           FORMAT('Warning! index=',I4,', chan ID=',I4,
     $          ', fastmodel freq=',F8.3,', RTP freq=',F8.3)
             ENDIF
             HEAD%vchan(I)=FREQ(I)
          ENDDO
       ELSE
          DO I=1,NCHAN
             HEAD%vchan(I)=FREQ(I)
          ENDDO
       ENDIF

C      ------------------------
C      Open the output RTP file
C      ------------------------
       MODE='c'
       ISTAT=rtpopen(FOUT, MODE, HEAD, HATT, PATT, IOPCO)

ccc
       if (DEBUG)  print *, 'read open status = ', ISTAT
ccc

C      -----------------------------------------------
C      All channels from sets 1, 2, and 3 are to use a
C      fake effective sun angle layer-to-space trans
c      at the end, NFAKE = NCHN1 + NCHN2 + NCHN3
C      -----------------------------------------------
       NFAKE=0
C 
       DO I=1,NCHN1
          NFAKE=NFAKE + 1
          INDFAK(NFAKE)=INDCHN( CLIST1(I) )
C          QUICKINDFAK(NFAKE)=INDCHN( CLIST1(I) )
       ENDDO
C
       DO I=1,NCHN2
          NFAKE=NFAKE + 1
          INDFAK(NFAKE)=INDCHN( CLIST2(I) )
C          QUICKINDFAK(NFAKE)=INDCHN( CLIST2(I) )
       ENDDO
C
       DO I=1,NCHN3
          NFAKE=NFAKE + 1
          INDFAK(NFAKE)=INDCHN( CLIST3(I) )
C          QUICKINDFAK(NFAKE)=INDCHN( CLIST3(I) )
       ENDDO

      DO I = 1,NFAKE
        III = intersect(I,INDFAK(1:NFAKE), NFAKE)  !! so I = INDFAK(III)
        QUICKINDFAK(I) = III
      END DO

C      ---------------------------
C      Start of loop over profiles
C      ---------------------------
       IPROF=1  ! initialize profile counter
C      Do you want this profile?
 10    LWANT=.TRUE.
       IF (NWANTP .GT. 0) THEN
C         Look for current profile on list of wanted profiles
          LWANT=.FALSE.
          DO I=1,NWANTP
             IF (IPROF .EQ. LISTP(I)) LWANT=.TRUE.
          ENDDO
       ENDIF

C      --------------
C      Read input RTP
C      --------------
       CALL RDRTP( LWANT, IPROF, IOPCI,
     $    IH2O, IO3, ICO, ICH4, ICO2, ISO2, IHNO3, IN2O, PTYPE,
     $    RALT, LCO2PM, NLAY, NEMIS, LAT, LON, SATANG, SATZEN,
     $    SALT, SUNANG, COSDAZ, PSURF, TSURF, CO2PPM,
     $    FEMIS, XEMIS, XRHO,
     $    TEMP, WAMNT, OAMNT, CAMNT, MAMNT, FAMNT, SAMNT, HAMNT, NAMNT,
     $     ALT, PROF, ISTAT )
C
       IF (ISTAT .EQ. -1) GOTO 9999  ! reached End Of File
C
       IF (.NOT. LWANT) THEN
C         Skip this profile
          IPROF=IPROF+ 1 
          GOTO 10
       ENDIF

C      -------------------------------------
C      Determine bottom layer, CO2, & angles
C      -------------------------------------
       CALL GETBOT(NLAY, PLEV, PSURF, LBOT, BLMULT)

C      Calc the fractional bottom layer air temperature
ccc
c       TEMP(LBOT)=TEMP(LBOT-1) + BLMULT*( TEMP(LBOT) - TEMP(LBOT-1) )
c Above line commented out & replaced by Scott Hannon, 24 July 2003.
c Mistakenly treats T at the center of the layer above as T at the
c bottom of the layer above.
ccc

C      CO2 profile switch
       IF (ICO2 .LT. 1) THEN
          LCO2=.FALSE.
       ELSE
          LCO2=.TRUE.
       ENDIF
C      N2O profile switch
       IF (IN2O .LT. 1) THEN
          LN2O=.FALSE.
       ELSE
          LN2O=.TRUE.
       ENDIF
C      SO2 profile switch
       IF (ISO2 .LT. 1) THEN
          LSO2=.FALSE.
       ELSE
          LSO2=.TRUE.
       ENDIF
C      NH3 profile switch
       IF (INH3 .LT. 1) THEN
          LNH3=.FALSE.
       ELSE
          LNH3=.TRUE.
       ENDIF
C      HNO3 profile switch
       IF (IHNO3 .LT. 1) THEN
          LHNO3=.FALSE.
       ELSE
          LHNO3=.TRUE.
       ENDIF
C      HDO switch (default .TRUE. from water)
       LHDO=.FALSE.
C
C
       IF (PTYPE .EQ. AIRSLAY) THEN
C         Copy pseudo level temperatures to another array
          DO I=1,LBOT
             TPSEUD(I)=TEMP(I)
          ENDDO
C         Convert temperatures
          CALL MEAN_T(LBOT, PLEV, PSURF, TPSEUD, TEMP)
C
       ELSE
C         Calc mean pressure for bottom fractional layer
          RJUNK1 = ( PSURF - PLEV(LBOT) )/LOG( PSURF/PLEV(LBOT) )
C         Do interpolation for fractional bottom layer mean temperature
C         assuming T is in linear in log(P)
          RJUNK2=( TEMP(LBOT) - TEMP(LBOT-1) )/
     $       LOG( PLAY(LBOT)/PLAY(LBOT-1) )             ! slope
          TEMP(LBOT)=RJUNK2*LOG( RJUNK1/PLAY(LBOT-1) ) + TEMP(LBOT - 1)
       ENDIF
C
C      Check satellite elevation
       IF (SALT .GT. 0.0) THEN
C         Warn and use default if invalid
          IF (SALT .LT. XSALT-50 .OR. SALT .GT. XSALT+50) THEN
             WRITE(IOINFO,1020) IPROF, SALT, XSALT
 1020        FORMAT('Warning! Profile',I5,
     $          ': replacing invalid input satellite altitude ',
     $          1PE11.4,' with default ',1PE11.4,' km')
             SALT=XSALT
          ENDIF
       ELSE
          SALT=XSALT
       ENDIF
C
C      Convert SATZEN or SATANG to viewing angle
       IF (SATZEN .GE. 0.0 .AND. SATZEN .LT. 63.0) THEN
C         Convert zenith angle at surface to view angle at satellite
          SVA=SACONV( SATZEN, SALT*1000.0 )/CONV
       ELSE
C         Check if scan angle is valid
          IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
C            View angle should be within a few degrees of scan angle
             SVA=ABS( SATANG )
          ELSE
             WRITE(IOERR,1030) IPROF, SATZEN, SATANG
 1030        FORMAT('Error! Profile',I5,
     $          ': invalid angles for SATZEN ',1PE11.4,
     $          ' and SATANG ',E11.4) 
             STOP
          ENDIF
       ENDIF

       ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
       IF (SVA .GT. ANGMAX) THEN
C         Truncate angle if too big
          WRITE(IOINFO,1040) IPROF, SVA
 1040     FORMAT('Warning! Profile',I5,': truncating view angle ',
     $       1PE11.4,' to 53 degrees')
          SVA=ANGMAX
       ENDIF

C      Convert from satellite to earth viewing angle (in radians)
       DO L=1,LBOT
             EVA=VACONV(SVA, SALT, ALT(L))
             SECANG(L)=1.0E+0/COS(EVA)
ccccccccccccc
c            for testing
c             SECANG(L)=SVA
ccccccccccccc
       ENDDO

C      Calc total sun angle secant
       DOSUN=.FALSE.
       SECSUN = 1.0
       SUNFDG = 1.0
       IF (SUNANG .GE. 0.0 .AND. SUNANG .LT. 89.9) DOSUN=.TRUE.
       IF (DOSUN) THEN
          SUNCOS=COS(CONV*SUNANG)
          SZALAY=SACONV(SUNANG,ALT(1))
          SCOS1=COS(SZALAY)
          RJUNK2=SECANG(LBOT) + 1.0/SUNCOS ! Total secant

C         Calc non-unity fudge factor if total secant > 9
          IF (RJUNK2 .GT. 9.0) THEN
C            fudge factor = true_total_secant/calc_total_secant
             SUNFDG=RJUNK2/9.0
C            truncated solar angle to use to calc SECSUN
             RJUNK1=ACOS( 1.0/(9.0 - SECANG(LBOT)) )/CONV
          ELSE
             SUNFDG=1.0
             RJUNK1=SUNANG
          ENDIF
c Should I change SUNFDG to SUNFDG(MAXLAY)?
C
          DO L=1,LBOT
             SZALAY=SACONV(RJUNK1,ALT(L))
             SECSUN(L)=SECANG(L) + 1.0E+0/COS(SZALAY)
          ENDDO

       ENDIF

C      -----------------------------------
C      Calculate the fast trans predictors
C      -----------------------------------
C
       CALL CALPAR (LBOT,
     $    RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RSAMNT,RHAMNT,RNAMNT,
     $    RAAMNT, TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, 
     $    NAMNT, AAMNT, RPRES,SECANG,   LAT,    FX,   RDZ,
     $     LCO2,  LN2O,  LSO2, LNH3, LHDO, LHNO3,LCO2PM,CO2PPM,CO2TOP,
     $   FIXMUL,CONPRD, DPRED,
     $   FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7,
     $   WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7,
     $   OPRED1,OPRED2,       OPRED4,OPRED5,OPRED6,OPRED7,
     $   MPRED3,CPRED4,TRCPRD,CO2MLT,SO2MLT,HNOMLT,N2OMLT,NH3MLT,HDOMLT)

         IF (DOSUN) THEN
           CALL SUNPAR ( LBOT,
     $       RTEMP, RWAMNT, ROAMNT, RCAMNT,
     $        TEMP,  WAMNT,  OAMNT,  CAMNT,
     $       RPRES,  SECSUN, SUNCONPRD,
     $       SUNFPRED4, SUNFPRED5, SUNFPRED6, SUNFPRED7,
     $       SUNWPRED4, SUNWPRED5, SUNWPRED6, SUNWPRED7,
     $       SUNOPRED4, SUNOPRED5, SUNOPRED6, SUNOPRED7,
     $       SUNCPRED4, SUNTRCPRD )
         END IF

C      -----------------------------------
C      Calculate the OPTRAN H2O predictors
C      -----------------------------------
       CALL CALOWP ( LBOT, WAMNT, RPRES, TEMP, SECANG, WAZOP, WAVGOP,
     $    WAANG, LOPMIN, LOPMAX, LOPUSE, H2OPRD, LOPLOW, DAOP )

C        Get basic cloud parameters from input RTP
         CALL GETCLD( IPROF, HEAD, PROF,
     $    LBLAC1, CTYPE1, CFRAC1, CPSIZ1, CPRTO1, CPRBO1, CNGWA1,
     $    XCEMI1, XCRHO1, CSTMP1,
     $    LBLAC2, CTYPE2, CFRAC2, CPSIZ2, CPRTO2, CPRBO2, CNGWA2,
     $    XCEMI2, XCRHO2, CSTMP2, CFRA12, FCLEAR, CFRA1X, CFRA2X )
         if (DEBUG) then
           print *,'getcld ',IPROF,CTYPE1, CFRAC1, CPSIZ1, CPRTO1,
     $                          CPRBO1, CNGWA1,CFRA1X     
         endif

C        ---------------------------------------------------
C        Set the emissivity & reflectivity for every channel
C        ---------------------------------------------------
         CALL SETEMS( NCHAN, NEMIS, FREQ, FEMIS, XEMIS, XRHO,
     $    XCEMI1, XCRHO1, XCEMI2, XCRHO2, LRHOT,
     $    EMIS, RHOSUN, RHOTHR, CEMIS1, CRHOS1, CRHOT1,
     $    CEMIS2, CRHOS2, CRHOT2,I) 
c        print *,CFRAC1,CFRAC2,CFRA12,LBLAC1,LBLAC2

C        Check and prepare (top) cloud1
         IF (CFRAC1 .GT. 0.0) THEN
           IF (LBLAC1) THEN
             CALL BKPREP(IPROF, 1, CTYPE1, CFRAC1, CPRTO1,
     $          LBOT, PSURF, PLEV, PLAY, TEMP, LCTOP1, TCTOP1,
     $          TEMPC1, CLRT1)
             IF (CSTMP1 .GT. 0.0) TCTOP1=CSTMP1
           ELSE
C             Determine which lookup table to use
              CALL GETMIE(CTYPE1,MIETYP,INDMI1,IERR1)
C             Prepare selected lookup table for given cpsize
              CALL CCPREP( NCHAN, LBOT, INDMI1, MIENPS,
     $          CNGWA1, CPSIZ1, CPRTO1, CPRBO1, PLEV, TEMP, SECANG,
     $          SECSUN, MIEPS, MIEABS, MIEEXT, MIEASY, LCBOT1, LCTOP1,
     $          CLRB1, CLRT1, TCBOT1, TCTOP1, MASEC1, MASUN1,
     $          CFRCL1, G_ASY1, NEXTO1, NSCAO1 )
            ENDIF
         ENDIF

C        Check and prepare (bottom) cloud2
         IF (CFRAC2 .GT. 0.0) THEN
            IF (LBLAC2) THEN
               CALL BKPREP(IPROF, 2, CTYPE2, CFRAC2, CPRTO2,
     $          LBOT, PSURF, PLEV, PLAY, TEMP, LCTOP2, TCTOP2,
     $          TEMPC2, CLRT2)
             IF (CSTMP2 .GT. 0.0) TCTOP2=CSTMP2
            ELSE
C            Determine which lookup table to use
             CALL GETMIE(CTYPE2,MIETYP,INDMI2,IERR2)
C            Prepare lookup data for cloud2
             CALL CCPREP( NCHAN, LBOT, INDMI2, MIENPS,
     $          CNGWA2, CPSIZ2, CPRTO2, CPRBO2, PLEV, TEMP, SECANG,
     $          SECSUN, MIEPS, MIEABS, MIEEXT, MIEASY, LCBOT2, LCTOP2,
     $          CLRB2, CLRT2, TCBOT2, TCTOP2, MASEC2, MASUN2,
     $          CFRCL2, G_ASY2, NEXTO2, NSCAO2 )
           ENDIF
         ELSE
C           Safe default for non-existant cloud2
            LCTOP2=1
         ENDIF

         SUNFAC=SUNCOS*PI*(RADSUN/DISTES)**2
C        Note: PI*(RADSUN/DISTES)^2 = solid angle [steradians] of
C        the sun as seen from Earth for the case DISTES >> RADSUN.

C------------------------------------------------------------------------
C      ----------------------
C      Loop over the channels
C      ----------------------
       DO I=1,NCHAN

C        ----------------------------------
C        Calculate the layer transmittances
C        ----------------------------------
C        Calculate TAU for set 1 thru 7

         IF (DEBUG) THEN
           DO II = 1,NCHAN
             !! print iI,h.ichan(iI),h.vchan(Ii)
             print *,II,LSTCHN(II),FREQ(II)
           END DO
         END IF

C      ---------------------------
C      Loop on channel (frequency)
C      eventually fills in matrix as    X=1,2,3,4,5,6,7
C       DO I=1,NCHNX
C          J=INDCHN( CLISTX(I) )
C          DO L = 1,100
C             Calc layer-to-space optical depth
C             KZ=KZ + KLAYER
C             TAUZ(ILAY,J)=KZ
C           END DO
C         END DO
C      ---------------------------

C        print *,INDCHN(CLIST1(1:NCHN1))
C        stop
!       DO II = 1,NCHN1
!         !! print iI,h.ichan(iI),h.vchan(Ii)
!         print *,II,LSTCHN(II),FREQ(II),CLIST1(II)
!       END DO

c         III = intersect(I,INDCHN(CLIST1(1:NCHN1)), NCHN1)
         III = QUICKCLIST1(I)
         IF (III .GT. 0) THEN
           CALL YCALT1( INDCHN,  LBOT,   NCHN1, CLIST1,  COEF1,
     $       FIXMUL, CONPRD, FPRED1, WPRED1, DPRED, OPRED1, TRCPRD,
     $       INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $       INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT,
     $       INDNH3, COFNH3, NH3MLT, INDHDO, COFHDO, HDOMLT,
     $       INDH2O, H2OPRD, COFH2O, LOPMIN, LOPMAX, LOPLOW,
     $       LOPUSE,   WAOP,   DAOP, WAANG,     TAU,   TAUZ,  III)
         END IF

c         III = intersect(I,INDCHN(CLIST2(1:NCHN2)), NCHN2)
         III = QUICKCLIST2(I)
         IF (III .GT. 0) THEN  
           CALL YCALT2( INDCHN, LBOT,   NCHN2, CLIST2,  COEF2,
     $      FIXMUL, CONPRD, FPRED2, OPRED2, WPRED2, DPRED, TRCPRD,
     $      INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $      INDHNO, COFHNO, HNOMLT, INDN2O, COFN2O, N2OMLT, 
     $      INDNH3, COFNH3, NH3MLT, INDHDO, COFHDO, HDOMLT,TAU, TAUZ, 
     $      III)
         END IF

c         III = intersect(I,INDCHN(CLIST3(1:NCHN3)), NCHN3)
         III = QUICKCLIST3(I)
         IF (III .GT. 0) THEN  
           CALL YCALT3( INDCHN,   LBOT,  NCHN3, CLIST3,  COEF3,
     $       FIXMUL, CONPRD, FPRED3, MPRED3, WPRED3, DPRED, TRCPRD,
     $       INDSO2, COFSO2, SO2MLT, INDHNO, COFHNO, HNOMLT,
     $       INDN2O, COFN2O, N2OMLT, INDNH3, COFNH3, NH3MLT,
     $       INDHDO, COFHDO, HDOMLT, INDH2O, H2OPRD, COFH2O, 
     $       LOPMIN, LOPMAX, LOPLOW, LOPUSE,
     $         WAOP,   DAOP,  WAANG,    TAU,   TAUZ, III)
          END IF

c         III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4)
         III = QUICKCLIST4(I)
         IF (III .GT. 0) THEN  
c           print *,'I,III,QUICKCLIST4(I) = ',I,III
c           print *,'I,III,QUICKCLIST4(I) = ',I,III,QUICKCLIST4(I)
           CALL YCALT4(INDCHN,   LBOT,  NCHN4, CLIST4,
     $       COEF4, FIXMUL, CONPRD, FPRED4, CPRED4, OPRED4, WPRED4,
     $       TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $       TAU,   TAUZ, III)
         END IF

c         III = intersect(I,INDCHN(CLIST5(1:NCHN5)), NCHN5)
         III = QUICKCLIST5(I)
         IF (III .GT. 0) THEN  
           CALL YCALT5(INDCHN,   LBOT,  NCHN5, CLIST5,
     $       COEF5, FIXMUL, CONPRD, FPRED5, WPRED5, OPRED5, 
     $       TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $       TAU,   TAUZ, III )
         END IF

c         III = intersect(I,INDCHN(CLIST6(1:NCHN6)), NCHN6)
         III = QUICKCLIST6(I)
         IF (III .GT. 0) THEN  
           CALL YCALT6(INDCHN,   LBOT,  NCHN6, CLIST6,
     $       COEF6, FIXMUL, CONPRD, FPRED6, WPRED6, OPRED6, DPRED, TRCPRD,
     $      INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $      INDN2O, COFN2O, N2OMLT, INDHDO, COFHDO, HDOMLT,  TAU,  TAUZ, 
     $      III )
         END IF

c         III = intersect(I,INDCHN(CLIST7(1:NCHN7)), NCHN7)
         III = QUICKCLIST7(I)
         IF (III .GT. 0) THEN  
           CALL YCALT7(INDCHN,   LBOT,  NCHN7, CLIST7,
     $       COEF7, FIXMUL, CONPRD, FPRED7, WPRED7, OPRED7,
     $       TRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $       TAU,   TAUZ, III )
         END IF

c************************************************************************
         IF (DOSUN) THEN
c           CALL SUNPAR ( LBOT,
c     $       RTEMP, RWAMNT, ROAMNT, RCAMNT,
c     $        TEMP,  WAMNT,  OAMNT,  CAMNT,
c     $       RPRES,  SECSUN, SUNCONPRD,
c     $       SUNFPRED4, SUNFPRED5, SUNFPRED6, SUNFPRED7,
c     $       SUNWPRED4, SUNWPRED5, SUNWPRED6, SUNWPRED7,
c     $       SUNOPRED4, SUNOPRED5, SUNOPRED6, SUNOPRED7,
c     $       SUNCPRED4, SUNTRCPRD )

C           Calc fake TAUZSN for sets 1, 2, and 3
c           III = intersect(I,INDFAK(1:NFAKE), NFAKE)  !! so I = INDFAK(III)
           III = QUICKINDFAK(I)
           IF (III .GT. 0) THEN
c              print *,'indfak',I,III,INDFAK(I),INDFAK(III),QUICKINDFAK(I)
              CALL FAKETZ( NFAKE, INDFAK, LBOT, TAUZ, SECANG,
     $         SECSUN, TAUZSN, III)
           END IF

c           III = intersect(I,INDCHN(CLIST4(1:NCHN4)), NCHN4) !! so I = INDCHN(CLIST4(III))
           III = QUICKCLIST4(I)
           IF (III .GT. 0) THEN  
             CALL YCALT4(INDCHN,   LBOT,  NCHN4, CLIST4,
     $         COEF4, FIXMUL, SUNCONPRD, SUNFPRED4, SUNCPRED4, SUNOPRED4, SUNWPRED4,
     $         SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $         TAUZSN, TAUZSN, III )
           END IF
C            ^^^^^^  ^^^^^^
C            dummy   actual

c           III = intersect(I,INDCHN(CLIST5(1:NCHN5)), NCHN5)
           III = QUICKCLIST5(I)
           IF (III .GT. 0) THEN  
             CALL YCALT5(INDCHN,   LBOT,  NCHN5, CLIST5,
     $         COEF5, FIXMUL, SUNCONPRD, SUNFPRED5, SUNWPRED5, SUNOPRED5,
     $         SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $         TAUZSN, TAUZSN, III )
           END IF

c           III = intersect(I,INDCHN(CLIST6(1:NCHN6)), NCHN6)
           III = QUICKCLIST6(I)
           IF (III .GT. 0) THEN  
             CALL YCALT6(INDCHN,   LBOT,  NCHN6, CLIST6,
     $          COEF6, FIXMUL, SUNCONPRD, SUNFPRED6, SUNWPRED6, SUNOPRED6, DPRED,
     $          SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDSO2, COFSO2, SO2MLT,
     $          INDN2O, COFN2O, N2OMLT, INDHDO, COFHDO, HDOMLT, TAUZSN, 
     $          TAUZSN, III )
           END IF

c           III = intersect(I,INDCHN(CLIST7(1:NCHN7)), NCHN7)
           III = QUICKCLIST7(I)
           IF (III .GT. 0) THEN  
             CALL YCALT7(INDCHN,   LBOT,  NCHN7, CLIST7,
     $          COEF7, FIXMUL, SUNCONPRD, SUNFPRED7, SUNWPRED7, SUNOPRED7,
     $          SUNTRCPRD, INDCO2, COFCO2, CO2MLT, INDN2O, COFN2O, N2OMLT,
     $          TAUZSN, TAUZSN, III )
           END IF
c           print *,'sun 7',TAU(LBOT,I),TAUZ(LBOT-1,I),TAUZSN(LBOT,I),TAUZSN(LBOT-1,I)

            IF (SUNFDG .GT. 1.0001) THEN
               DO II=1,NCHAN
                  DO L=1,LBOT
                     TAUZSN(L,II)=TAUZSN(L,II)*SUNFDG
                  ENDDO
               ENDDO
            ENDIF
         ELSE
C           DOSUN = 'FALSE'; No sun; set the sun surface-to-space trans to zero
            SUNCOS=0.0
C            DO II=1,NCHAN
C              DO L=1,LBOT
C                TAUZSN(L,II)=0.0
C              ENDDO
C            ENDDO

C            DO L=1,LBOT      !!! outer loop is slowest
C              DO II=1,NCHAN  !!! inner loop is fastest
C                TAUZSN(L,II)=0.0
C              ENDDO
C            ENDDO
C           TAUZSN(1:LBOT,1:NCHAN) = 0.0
           TAUZSN(1:LBOT,I) = 0.0

         ENDIF
c************************************************************************

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        Calculate cloudy radiance

C         Radiation constants for current channel
          C1V3=C1*(FREQ(I)**3)
          C2V=C2*FREQ(I)

c          IF (NWANTC .GT. 0) THEN
c            write(*,'(A,I5,7(F12.4))') 'start rads/tau',I,FREQ(I),TSURF,TEMP(LBOT),
c     $             TAU(LBOT,I),TAUZ(LBOT-1,I),TAUZSN(LBOT,I),TAUZSN(LBOT-1,I)
c          END IF

C         Calculate Planck & clear airs trans for full layers
          DO L=1,LBOT-1
             RPLNCK(L)=C1V3/( EXP( C2V/TEMP(L) ) - 1.0 )
             TRANL(L)=QIKEXP( -TAU(L,I) )
          ENDDO
C         Note: TEMP(LBOT) already adjusted for bottom fractional layer
          RPLNCK(LBOT)=C1V3/( EXP( C2V/TEMP(LBOT) ) - 1.0 )

C         Calculate clear airs trans for bottom fractional layer
          RJUNK1=-TAU(LBOT,I)*BLMULT
          TRANL(LBOT)=QIKEXP( RJUNK1 )
          TRANL(LBOT)=QIKEXP( RJUNK1 )
          TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LBOT-1,I) )
          TRANS(I)=QIKEXP( BLMULT*(TAUZSN(LBOT-1,I)-TAUZSN(LBOT,I)) -
     $       TAUZSN(LBOT-1,I) )

C         Planck for surface
          RSURFE=EMIS(I)*C1V3/( EXP( C2V/TSURF ) - 1.0 )

C         Calculate clear radiance
          IF (FCLEAR .GT. 0.0) THEN
             CALL CALRAD0( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $       TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $       RHOTHR, LABOVE, COEFF, RAD0 )
          ELSE
             RAD0=0.0
          ENDIF

C         Store original values
          VSTORE(1)=TRANL(LCTOP2)
          VSTORE(2)=TRANZ(I)
          VSTORE(3)=TRANS(I)
          VSTORE(4)=RHOTHR(I)
          VSTORE(5)=RHOSUN(I)
          VSTORE(6)=RPLNCK(LCTOP2)
C         Updates for new surface if bottom cloud2 is black
          IF (CFRAC2 .GT. 0.0 .AND. LBLAC2) THEN
             RJUNK1=-TAU(LCTOP2,I)*CLRT2
             TRANL(LCTOP2)=QIKEXP( RJUNK1 )
             TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LCTOP2-1,I) )
             TRANS(I)=QIKEXP( CLRT2*(TAUZSN(LCTOP2-1,I)-
     $          TAUZSN(LCTOP2,I)) - TAUZSN(LCTOP2-1,I) )
             RSURFC=CEMIS2(I)*C1V3/( EXP( C2V/TCTOP2 ) - 1.0 )
             RHOTHR(I)=CRHOT2(I)
             RHOSUN(I)=CRHOS2(I)
             RPLNCK(LCTOP2)=C1V3/( EXP( C2V/TEMPC2 ) - 1.0 )
c             RSURFC=C1V3/( EXP( C2V/TEMPC2 ) - 1.0 )
          ENDIF

C         Calculate bottom cloud2 radiance
          IF (CFRA2X .GT. 0.0) THEN
             IF (LBLAC2) THEN
                CALL CALRAD0( DOSUN, I, LCTOP2, RPLNCK, RSURFC, SECANG,
     $          TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, RADC2 )
             ELSE
                CALL CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $          TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, CFRCL2, MASEC2, MASUN2, COSDAZ,
     $          NEXTO2, NSCAO2, G_ASY2, LCTOP2, LCBOT2, RADC2 )
             ENDIF
          ELSE
             RADC2=0.0
          ENDIF

C         Calculate combined cloud1+cloud2 radiance
          IF (CFRA12 .GT. 0.0) THEN
             IF (LBLAC2) THEN
                CALL CALRAD1( DOSUN, I, LCTOP2, RPLNCK, RSURFC, SECANG,
     $          TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, COSDAZ,
     $          NEXTO1, NSCAO1, G_ASY1, LCTOP1, LCBOT1, RADC12 )
             ELSE
                CALL CALRAD2( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $          TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, NEXTO1,
     $          NSCAO1, G_ASY1, LCTOP1, LCBOT1, CFRCL2, MASEC2, MASUN2,
     $          COSDAZ, NEXTO2, NSCAO2, G_ASY2, LCTOP2, LCBOT2, RADC12 )
             ENDIF
          ELSE
             RADC12=0.0
          ENDIF

C         Restore original values
          TRANL(LCTOP2)=VSTORE(1)
          TRANZ(I)=VSTORE(2)
          TRANS(I)=VSTORE(3)
          RHOTHR(I)=VSTORE(4)
          RHOSUN(I)=VSTORE(5)
          RPLNCK(LCTOP2)=VSTORE(6)
C         Updates for new surface if top cloud1 is black
          IF (CFRAC1 .GT. 0.0 .AND. LBLAC1) THEN
             RJUNK1=-TAU(LCTOP1,I)*CLRT1
             TRANL(LCTOP1)=QIKEXP( RJUNK1 )
             TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LCTOP1-1,I) )
             TRANS(I)=QIKEXP( CLRT1*(TAUZSN(LCTOP1-1,I)-
     $          TAUZSN(LCTOP1,I)) - TAUZSN(LCTOP1-1,I) )
             RSURFC=CEMIS1(I)*C1V3/( EXP( C2V/TCTOP1 ) - 1.0 )
             RHOTHR(I)=CRHOT1(I)
             RHOSUN(I)=CRHOS1(I)
             RPLNCK(LCTOP1)=C1V3/( EXP( C2V/TEMPC1 ) - 1.0 )
c             RSURFC=C1V3/( EXP( C2V/TEMPC1 ) - 1.0 )
          ENDIF

C         Calculate top cloud1 radiance
          IF (CFRA1X .GT. 0.0) THEN
             IF (LBLAC1) THEN
                CALL CALRAD0( DOSUN, I, LCTOP1, RPLNCK, RSURFC, SECANG,
     $          TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, RADC1 )
             ELSE
                CALL CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $          TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, COSDAZ,
     $          NEXTO1, NSCAO1, G_ASY1, LCTOP1, LCBOT1, RADC1 )
             ENDIF
          ELSE
             RADC1=0.0
          ENDIF

C         Total the clear & various cloudy radiances
          RAD(I)=RAD0*FCLEAR + RADC1*CFRA1X + RADC2*CFRA2X +
     $       RADC12*CFRA12

          IF (NWANTC .GT. 0) THEN
            write(*,'(A,I5,13(F12.4))') 'rads',I,FREQ(I),EMIS(I),TSURF,FCLEAR,CFRA1X,CFRA2X,CFRA12,RSURFE,RAD0,RADC1,RADC2,RADC12
          END IF

ccc this block for testing
       IF (I .EQ. 1291) THEN
c         print *,'chan1291 : iPROF,rad0,radc1,radc2,radc12,FINAL=',
c     $      IPROF,RAD0,RADC1,RADC2,RADC12,RAD(I)
C         print *,'chan1291 : IPROF,rad0,FCLEAR,CFRA1X,CFRA2X,CFRA12=',
C     $      IPROF,RAD0,FCLEAR,CFRA1X,CFRA2X,CFRA12
c         PRINT *,'CLOUD1 emis,temp = ',CEMIS1(I),TCTOP1
c         PRINT *,'CLOUD2 emis,temp = ',CEMIS2(I),TCTOP2
       endif
ccc

       ENDDO ! channels


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      -----------------
C      Calculate non-LTE
C      -----------------
C      comment: the nonLTE calculation does not consider cloud effects,
C      but clouds are generally below the altitude where nonLTE occurs.
       IF (DOSUN) THEN
          CALL CALNTE ( INDCHN, TEMP, SUNCOS, SCOS1, SECANG(1),
     $       NCHNTE, CLISTN, COEFN, CO2TOP, RAD )
       ENDIF
C

C      -------------------
C      Output the radiance
C      -------------------
       CALL WRTRTP(IPROF, IOPCO, NCHAN, RAD, PROF, NWANTC, RINDCHN)
C

C      ----------------------
C      End loop over profiles
C      ----------------------
       IPROF=IPROF + 1  ! increment profile counter
       GOTO 10
C

C      -------------------
C      Close the RTP files
C      -------------------
 9999  ISTAT=rtpclose(IOPCI)
       ISTAT=rtpclose(IOPCO)
C
       STOP
       END
