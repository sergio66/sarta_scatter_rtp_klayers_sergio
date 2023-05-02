c     This version with rtp input profiles and command-line arguments
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
C
c*********************************************************************************************************************
C SSM : H2OJACPRD,DAOPJAC typically size(OPTRANJAC,MXLAY,MXCHAN) where the 2 derivatives are T,G1 (OPTRAN)
C SSM : XJACPRD           typically   size(MAXJAC, MXLAY,MXCHAN) where the   derivatives are T,G1,G2,G3,G4,G5,G6,G9,G12
C
C --------------------------------------------
C  IWHICHJAC  NAME           ID
C --------------------------------------------
C    1        T              100
C    2        WV             1
C    3        O3             3
C    4        CO2            2
C    5        N2O            4
C    6        CO             5
C    7        CH4            6
C    8        SO2            9
C    9        HNO3           12
C    10       NH3            11   eventually
C    11       HDO            103  eventually
C --------------------------------------------
C
c*********************************************************************************************************************

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
C       REAL VACONV
C       REAL SACONV


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
       LOGICAL  LRHOT          ! force refl therm rho=(1-emis)/pi?
       INTEGER NUMPROF,NUMCHAN ! number of channels, profiles in rtp file
       INTEGER NWANTP          ! number of wanted profiles (default -1=all)
       INTEGER  LISTP(MAXPRO)  ! list of wanted profiles
       INTEGER NWANTC          ! number of wanted channels (default -1=all)
       INTEGER  LISTC(MAXPRO)  ! list of wanted channels
       INTEGER NWANTJ          ! number of wanted jacs (default 0=none)
       INTEGER  LISTJ(MAXPRO)  ! list of wanted channels

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
       REAL COFH2O(  NH2O,MXOWLY,MXCHNW) ! coefs for OPTRAN H2O
       REAL     FX(MAXLAY)               ! fixed gases adjustment
       REAL   WAZOP(MXOWLY)              ! OPTRAN water l-to-s amounts
       REAL  WAVGOP(NOWAVG,MXOWLY)       ! OPTRAN raw predictor averages
       
       REAL RADNTE                       ! temporary NLTE rad                        
       INTEGER QUICKINDNTE(MXCHAN)       ! list of non-LTE channels
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
c       REAL TPSEUD(MAXLAY)
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
       REAL  SUNDPRED(  NHDO,MAXLAY) ! HDO perturbation predictors

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
C FOR JACS
       REAL CO2JACMLT(MAXLAY)
       REAL SO2JACMLT(MAXLAY)
       REAL HNOJACMLT(MAXLAY)
       REAL N2OJACMLT(MAXLAY)
       REAL NH3JACMLT(MAXLAY)
       REAL HDOJACMLT(MAXLAY)
       LOGICAL DOJAC
       REAL DTAU_DTZ(MAXLAY,MXCHAN),DTAU_DG1(MAXLAY,MXCHAN),DTAU_DG2(MAXLAY,MXCHAN),DTAU_DG3(MAXLAY,MXCHAN)
       REAL DTAU_DG4(MAXLAY,MXCHAN),DTAU_DG5(MAXLAY,MXCHAN),DTAU_DG6(MAXLAY,MXCHAN),DTAU_DG9(MAXLAY,MXCHAN)
       REAL DTAU_DG12(MAXLAY,MXCHAN)
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
c       REAL RPLNCK(MAXLAY) ! layer Planck
c       REAL RSURFE         ! surface emission
c       REAL RSURFC         ! black cloud surface emission
c       REAL  TRANL(MAXLAY) ! clear air layer transmittance
c       REAL  TRANZ(MXCHAN) ! clear air layer-to-space transmittance
c       REAL  TRANS(MXCHAN) ! clear air total reflected solar trans
       REAL  TSURF         ! surface temperature
       REAL    RAD(MXCHAN) ! chan radiance
C      For clear/cloudy radiances
c       REAL   RAD0         ! radiance no clouds
c       REAL  RADC1         ! radiance cloud1
c       REAL  RADC2         ! radiance cloud2
c       REAL RADC12         ! radiance cloud1+cloud2
C
C      for RDSUN
       REAL   HSUN(MXCHAN) ! sun radiance (direct from sun)
C      Other variables for the sun
       REAL SUNANG         ! solar zenith angle (at 0 altitude)
       REAL COSDAZ         ! cosine(solazi - satazi) {COS Delta AZimuth}
c       REAL SZALAY         ! solar zenith angle in some layer
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
c       INTEGER INDMI1  ! index in MIETYP for CTYPE1
c       INTEGER INDMI2  ! index in MIETYP for CTYPE2
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
c       REAL    EVA         ! (Earth) local viewing angle
c       REAL   CONV         ! degrees to radians conversion factor
c       REAL ANGMAX         ! maximum allowed viewing angle
       REAL RJUNK1         ! junk/work
       REAL RJUNK2         ! another junk/work
       REAL CO2PPM         ! Profile mean dry air CO2 mixing ratio
       REAL PLAY(MAXLAY)   ! layer mean pressure
c       REAL C1V3           ! rad constant c1 times freq^3
c       REAL C2V            ! rad constant c2 times freq
c       REAL VSTORE(6)      ! temporary storage for various variables

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

C      for jacobians
       REAL  DAOPJAC(OPTRANJAC, MAXLAY)        !!!! OPTRAN
       REAL  H2OJACPRD(OPTRANJAC,NH2O,MXOWLY)  !!!! OPTRAN
       REAL JAC_ST_C(MXCHAN),         JAC_ST_1(MXCHAN),         JAC_ST_2(MXCHAN),         JAC_ST_12(MXCHAN)
       REAL JAC_TZ_C(MAXLAY,MXCHAN),  JAC_TZ_1(MAXLAY,MXCHAN),  JAC_TZ_2(MAXLAY,MXCHAN),  JAC_TZ_12(MAXLAY,MXCHAN)
       REAL JAC_G1_C(MAXLAY,MXCHAN),  JAC_G1_1(MAXLAY,MXCHAN),  JAC_G1_2(MAXLAY,MXCHAN),  JAC_G1_12(MAXLAY,MXCHAN)
       REAL JAC_G2_C(MAXLAY,MXCHAN),  JAC_G2_1(MAXLAY,MXCHAN),  JAC_G2_2(MAXLAY,MXCHAN),  JAC_G2_12(MAXLAY,MXCHAN)
       REAL JAC_G3_C(MAXLAY,MXCHAN),  JAC_G3_1(MAXLAY,MXCHAN),  JAC_G3_2(MAXLAY,MXCHAN),  JAC_G3_12(MAXLAY,MXCHAN)
       REAL JAC_G4_C(MAXLAY,MXCHAN),  JAC_G4_1(MAXLAY,MXCHAN),  JAC_G4_2(MAXLAY,MXCHAN),  JAC_G4_12(MAXLAY,MXCHAN)
       REAL JAC_G5_C(MAXLAY,MXCHAN),  JAC_G5_1(MAXLAY,MXCHAN),  JAC_G5_2(MAXLAY,MXCHAN),  JAC_G5_12(MAXLAY,MXCHAN)
       REAL JAC_G6_C(MAXLAY,MXCHAN),  JAC_G6_1(MAXLAY,MXCHAN),  JAC_G6_2(MAXLAY,MXCHAN),  JAC_G6_12(MAXLAY,MXCHAN)
       REAL JAC_G9_C(MAXLAY,MXCHAN),  JAC_G9_1(MAXLAY,MXCHAN),  JAC_G9_2(MAXLAY,MXCHAN),  JAC_G9_12(MAXLAY,MXCHAN)
       REAL JAC_G12_C(MAXLAY,MXCHAN), JAC_G12_1(MAXLAY,MXCHAN), JAC_G12_2(MAXLAY,MXCHAN), JAC_G12_12(MAXLAY,MXCHAN)
       REAL JAC_WGT_C(MAXLAY,MXCHAN), JAC_WGT_1(MAXLAY,MXCHAN), JAC_WGT_2(MAXLAY,MXCHAN), JAC_WGT_12(MAXLAY,MXCHAN)
       REAL TAU4(4,MAXLAY,MXCHAN) !  chan layer effective optical depth for CLR,CLD1,CLD2,CLD12
       REAL RAD4(4,MAXLAY,MXCHAN) ! -chan radiance + planck(TL)         for CLR,CLD1,CLD2,CLD12
       REAL L2S4(4,MAXLAY,MXCHAN),WGT4(4,MAXLAY,MXCHAN)
       REAL DBTDT(MAXLAY,MXCHAN)  ! dBT(T,L)/dT
       INTEGER IOUNTZ,IOUNG1,IOUNG2,IOUNG3,IOUNG4,IOUNG5,IOUNG6,IOUNG9,IOUNG11,IOUNG12,IOUNG103,IOUNWGT,iFileErr
       INTEGER JAC_OUTPUT_UNITS           ! 0 for drad/dT and drad/dq, 1 for dBT/dT and dBT/d(log q) = q dBT/dq
       CHARACTER*180 caJacTZ,caJACWGT,caJACG1,caJACG2,caJACG3,caJACG4,caJACG5,
     $               caJACG6,caJACG9,caJACG11,caJACG12,caJACG103

       INTEGER ijunk
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
       CALL RDINFO(FIN, FOUT, LRHOT, NWANTP, LISTP, NWANTC, LISTC, 
     $             NWANTJ, LISTJ, NUMCHAN, NUMPROF, 
     $     caJacTZ,caJACWGT,caJACG1,caJACG2,caJACG3,caJACG4,caJACG5,
     $     caJACG6,caJACG9,caJACG11,caJACG12,caJACG103)
ccc
       if (DEBUG) then
         print *, 'nwantp=', NWANTP
         print *, 'listp=', (LISTP(I),I=1,NWANTP)
         print *, 'FIN = ',FIN
         print *, 'FOUT = ',FOUT
       endif
ccc
       if (nwantp .gt. 0) then
         print *,'nwantp = ',nwantp,'followed by list ....'
         print *,listp(1:nwantp)
       else
         print *,'nwantp = -1 (all profs)'
       end if

       if (nwantc .gt. 0) then
         print *,'nwantc = ',nwantc,'followed by list ....'
         print *,listc(1:nwantc)
       else
         print *,'nwantc = -1 (all chans)'
       end if

       DOJAC = .FALSE.
       IF (NWANTJ .GT. 0) THEN
         DOJAC = .TRUE.
         print *,'for NUMPROF = ',numprof,' profiles with NUMCHAN = ',numchan,' channels '
         print *, 'want this # jacs nwantj = ',nwantj,' followed by list (100=ST/T, 200=WGT, 1,3 = WV/OZ)....'
         print *,listj(1:nwantj)
         IOUNWGT = 400
         IOUNTZ  = 200
         IOUNG1  = 21
         IOUNG2  = 22
         IOUNG3  = 23
         IOUNG4  = 24
         IOUNG5  = 25
         IOUNG6  = 26
         IOUNG9  = 29
         IOUNG11 = 211
         IOUNG12 = 212
         IOUNG103 = 2103

         IF (INTERSECT(1,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG1 = ',caJacG1
           OPEN(UNIT=IOUNG1,FILE=caJacG1,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG1) NUMPROF
           WRITE(IOUNG1) NUMCHAN
         END IF
         IF (INTERSECT(2,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG2 = ',caJacG2
           OPEN(UNIT=IOUNG2,FILE=caJacG2,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG2) NUMPROF
           WRITE(IOUNG2) NUMCHAN
         END IF
         IF (INTERSECT(3,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG3 = ',caJacG3
           OPEN(UNIT=IOUNG3,FILE=caJacG3,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG3) NUMPROF
           WRITE(IOUNG3) NUMCHAN
         END IF
         IF (INTERSECT(4,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG4 = ',caJacG4
           OPEN(UNIT=IOUNG4,FILE=caJacG4,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG4) NUMPROF
           WRITE(IOUNG4) NUMCHAN
         END IF
         IF (INTERSECT(5,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG5 = ',caJacG5
           OPEN(UNIT=IOUNG5,FILE=caJacG5,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG5) NUMPROF
           WRITE(IOUNG5) NUMCHAN
         END IF
         IF (INTERSECT(6,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG6 = ',caJacG6
           OPEN(UNIT=IOUNG6,FILE=caJacG6,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG6) NUMPROF
           WRITE(IOUNG6) NUMCHAN
         END IF
         IF (INTERSECT(9,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG9 = ',caJacG9
           OPEN(UNIT=IOUNG9,FILE=caJacG9,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG9) NUMPROF
           WRITE(IOUNG9) NUMCHAN
         END IF
         IF (INTERSECT(11,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG11 = ',caJacG11
           OPEN(UNIT=IOUNG11,FILE=caJacG11,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG11) NUMPROF
           WRITE(IOUNG11) NUMCHAN
         END IF
         IF (INTERSECT(12,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG12 = ',caJacG12
           OPEN(UNIT=IOUNG12,FILE=caJacG12,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG12) NUMPROF
           WRITE(IOUNG12) NUMCHAN
         END IF
         IF (INTERSECT(103,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacG103 = ',caJacG103
           OPEN(UNIT=IOUNG103,FILE=caJacG103,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNG103) NUMPROF
           WRITE(IOUNG103) NUMCHAN
         END IF

         IF (INTERSECT(100,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacTZ = ',caJacTZ
           OPEN(UNIT=IOUNTZ,FILE=caJacTZ,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNTZ) NUMPROF
           WRITE(IOUNTZ) NUMCHAN
         END IF
         IF (INTERSECT(200,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
           write(*,'(A,A)') 'opening jac file caJacWGT = ',caJacWGT
           OPEN(UNIT=IOUNWGT,FILE=caJacWGT,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
           WRITE(IOUNWGT) NUMPROF
           WRITE(IOUNWGT) NUMCHAN
         END IF
       END IF

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

C      ------------------------
C      Read cloud lookup tables
C      ------------------------
       CALL RDCLDT( IOUN, INDCHN, MIETYP, FNMIEA, FNMIEE, FNMIEG,
     $    MIENPS, MIEPS, MIEABS, MIEEXT, MIEASY )       

C      --------------------------------------------------------------
C      Read the coef data files and apply multiplier tuning to coeffs
C      --------------------------------------------------------------
       CALL RDCOEF_TUNMLT (IOUN, NCHAN, INDCHN, SETCHN, FCHAN, 
     $  NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7,
     $ CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
     $  COEF1,  COEF2,  COEF3,  COEF4,  COEF5,  COEF6,  COEF7,
     $   FREQ, LABOVE,  COEFF, INDCO2, COFCO2, INDSO2, COFSO2,
     $ INDHNO, COFHNO, INDN2O, COFN2O, INDNH3, COFNH3, 
     $ INDHDO, COFHDO,
     $ INDH2O,  WAZOP, WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN,
     $ QUICKCLIST1, QUICKCLIST2, QUICKCLIST3, QUICKCLIST4, QUICKCLIST5, 
     $ QUICKCLIST6, QUICKCLIST7,
     $ NWANTC, LISTC, LSTCHN, RINDCHN)
       IF (NCHAN .EQ. 2645) THEN
         write(*,'(A)') 'SET(X) NCNHN(X)  IMIN(X) IMAX(X)  FMIN(X)  FMAX(X)'
         write(*,'(A)') '--------------------------------------------------'
         write(*,'(4(I8),2(F8.2))') 1,NCHN1,CLIST1(1),CLIST1(NCHN1),
     $                                      MINVAL(FREQ(INDCHN(CLIST1(1:NCHN1)))),MAXVAL(FREQ(INDCHN(CLIST1(1:NCHN1))))
         write(*,'(4(I8),2(F8.2))') 2,NCHN2,CLIST2(1),CLIST2(NCHN2),
     $                                      MINVAL(FREQ(INDCHN(CLIST2(1:NCHN2)))),MAXVAL(FREQ(INDCHN(CLIST2(1:NCHN2))))
         write(*,'(4(I8),2(F8.2))') 3,NCHN3,CLIST3(1),CLIST3(NCHN3),
     $                                      MINVAL(FREQ(INDCHN(CLIST3(1:NCHN3)))),MAXVAL(FREQ(INDCHN(CLIST3(1:NCHN3))))
         write(*,'(4(I8),2(F8.2))') 4,NCHN4,CLIST4(1),CLIST4(NCHN4),
     $                                      MINVAL(FREQ(INDCHN(CLIST4(1:NCHN4)))),MAXVAL(FREQ(INDCHN(CLIST4(1:NCHN4))))
         write(*,'(4(I8),2(F8.2))') 5,NCHN5,CLIST5(1),CLIST5(NCHN5),
     $                                      MINVAL(FREQ(INDCHN(CLIST5(1:NCHN5)))),MAXVAL(FREQ(INDCHN(CLIST5(1:NCHN5))))
         write(*,'(4(I8),2(F8.2))') 6,NCHN6,CLIST6(1),CLIST6(NCHN6),
     $                                      MINVAL(FREQ(INDCHN(CLIST6(1:NCHN6)))),MAXVAL(FREQ(INDCHN(CLIST6(1:NCHN6))))
         write(*,'(4(I8),2(F8.2))') 7,NCHN7,CLIST7(1),CLIST7(NCHN7),
     $                                      MINVAL(FREQ(INDCHN(CLIST7(1:NCHN7)))),MAXVAL(FREQ(INDCHN(CLIST7(1:NCHN7))))

cc       write(*,'(4(I8),2(F8.2))') 1,NCHN1,CLIST1(1),CLIST1(NCHN1),FREQ(INDCHN(CLIST1(1))),FREQ(INDCHN(CLIST1(NCHN1)))
cc       write(*,'(4(I8),2(F8.2))') 2,NCHN2,CLIST2(1),CLIST2(NCHN2),FREQ(INDCHN(CLIST2(1))),FREQ(INDCHN(CLIST2(NCHN2)))
cc       write(*,'(4(I8),2(F8.2))') 3,NCHN3,CLIST3(1),CLIST3(NCHN3),FREQ(INDCHN(CLIST3(1))),FREQ(INDCHN(CLIST3(NCHN3)))
cc       write(*,'(4(I8),2(F8.2))') 4,NCHN4,CLIST4(1),CLIST4(NCHN4),FREQ(INDCHN(CLIST4(1))),FREQ(INDCHN(CLIST4(NCHN4)))
cc       write(*,'(4(I8),2(F8.2))') 5,NCHN5,CLIST5(1),CLIST5(NCHN5),FREQ(INDCHN(CLIST5(1))),FREQ(INDCHN(CLIST5(NCHN5)))
cc       write(*,'(4(I8),2(F8.2))') 6,NCHN6,CLIST6(1),CLIST6(NCHN6),FREQ(INDCHN(CLIST6(1))),FREQ(INDCHN(CLIST6(NCHN6)))
cc       write(*,'(4(I8),2(F8.2))') 7,NCHN7,CLIST7(1),CLIST7(NCHN7),FREQ(INDCHN(CLIST7(1))),FREQ(INDCHN(CLIST7(NCHN7)))
       END IF

C      Calc OPTRAN absorption coefficient scaling factor WAOP
       WAOP(1)=WAZOP(1)
       DO L=2,MXOWLY
          WAOP(L)=WAZOP(L) - WAZOP(L-1)
       ENDDO

C      --------------------------
C      Read in the solar radiance
C      --------------------------
       CALL RDSUN_SET_NFAKE(IOUN, INDCHN, HSUN, 
     $         NCHN1,NCHN2,NCHN3,CLIST1,CLIST2,CLIST3,
     $         NFAKE,INDFAK,QUICKINDFAK)

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

C************************************************************************
C   this is end of generic stuff eg read in HEAD = nchan,ichan,vchan, predictors
C************************************************************************

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

       IF (ISTAT .EQ. -1) GOTO 9999  ! reached End Of File
C
       IF (.NOT. LWANT) THEN
C         Skip this profile
          IPROF=IPROF+ 1 
          GOTO 10
       ENDIF

       CALL get_lbot_fix_salt_temp(
     $   NLAY, PLEV, PLAY, PSURF, LBOT, BLMULT,                 !!! for getbot
     $   TEMP, SALT, IPROF, AIRSLAY)
 
       CALL CALC_SVA_SECANG_SECSUN(SALT,SATZEN,SATANG,SUNANG,LBOT,ALT,
     $      IPROF,LSTCHN,NCHNTE,CLISTN,FREQ,INDCHN,
     $      SVA,SECANG,SECSUN,SUNFDG,SUNCOS,SCOS1,DOSUN,QUICKINDNTE)

c************************************************************************

C      ---------------------------------------------------------
C      ---------------------------------------------------------
C      Calculate the fast trans predictors and OPTRAN predictors
C      which directly use T(z),WV(z),O3(z) 
C      ---------------------------------------------------------
C      ---------------------------------------------------------
C
       IF (DOJAC .EQV. .FALSE.) THEN
         CALL CALPAR (LBOT,
     $      RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RSAMNT,RHAMNT,RNAMNT,
     $      RAAMNT, TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, 
     $      NAMNT, AAMNT, RPRES,SECANG,   LAT,    FX,   RDZ,
     $      LCO2,  LN2O,  LSO2, LNH3, LHDO, LHNO3,LCO2PM,
     $      CO2PPM,CO2TOP,FIXMUL,
     $      CONPRD, DPRED,
     $     FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7,
     $     WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7,
     $     OPRED1,OPRED2,       OPRED4,OPRED5,OPRED6,OPRED7,
     $     MPRED3,CPRED4,TRCPRD,CO2MLT,SO2MLT,HNOMLT,N2OMLT,NH3MLT,HDOMLT)
      ELSE
        CALL YCALPAR_JAC ( LBOT,
     $      RTEMP,RFAMNT,RWAMNT,ROAMNT,RCAMNT,RMAMNT,RSAMNT,RHAMNT,RNAMNT,
     $      RAAMNT, TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, 
     $      NAMNT, AAMNT, RPRES,SECANG,   LAT,    FX,   RDZ,
     $    LCO2,  LN2O,  LSO2, LNH3,  LHDO, LHNO3,LCO2PM,
     $    CO2PPM,CO2TOP,FIXMUL,
     $    CONPRD,DPRED, 
     $    FPRED1,FPRED2,FPRED3,FPRED4,FPRED5,FPRED6,FPRED7,
     $    WPRED1,WPRED2,WPRED3,WPRED4,WPRED5,WPRED6,WPRED7,
     $    OPRED1,OPRED2,       OPRED4,OPRED5,OPRED6,OPRED7,
     $    MPRED3,CPRED4,TRCPRD,
     $    CO2MLT,SO2MLT,HNOMLT,N2OMLT,NH3MLT,HDOMLT,
     $    DOJAC,LISTJ,NWANTJ,CONJACPRD,DJACPRED, 
     $    FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $    WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $    OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $    MJACPRED3,CJACPRED4,TRCJACPRD,
     $    CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT)
      END IF

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
       CALL YCALOWP ( LBOT, WAMNT, RPRES, TEMP, SECANG, WAZOP, WAVGOP,
     $    WAANG, LOPMIN, LOPMAX, LOPUSE, H2OPRD, LOPLOW, DAOP, 
     $    DOJAC, H2OJACPRD, DAOPJAC )

c************************************************************************

C      -----------------------------------
C      -----------------------------------
C      Calculate the CLOUD, EMIS, REFL param
C      -----------------------------------
C      -----------------------------------

C      ---------------------------------------------------
C      Set the emissivity & reflectivity for every channel
C      ---------------------------------------------------
       CALL SETEMS( NCHAN, NEMIS, FREQ, FEMIS, XEMIS, XRHO,
     $    XCEMI1, XCRHO1, XCEMI2, XCRHO2, LRHOT,
     $    EMIS, RHOSUN, RHOTHR, CEMIS1, CRHOS1, CRHOT1,
     $    CEMIS2, CRHOS2, CRHOT2) 
c        print *,CFRAC1,CFRAC2,CFRA12,LBLAC1,LBLAC2

C      --------------------------------------
C      Set the cloud fractions and mie tables
C      --------------------------------------
       CALL prepare_clds(
     $    LBLAC1, CTYPE1, CFRAC1, CPSIZ1, CPRTO1, CPRBO1, CNGWA1,  ! from GETCLD
     $    XCEMI1, XCRHO1, CSTMP1,                                  ! from GETCLD
     $    LBLAC2, CTYPE2, CFRAC2, CPSIZ2, CPRTO2, CPRBO2, CNGWA2,  ! from GETCLD
     $    XCEMI2, XCRHO2, CSTMP2, CFRA12, FCLEAR, CFRA1X, CFRA2X,  ! from GETCLD
     $    NCHAN, IPROF, LBOT, PSURF, PLEV, PLAY, TEMP, SECANG, SECSUN,
     $    MIETYP, MIENPS, MIEPS, MIEABS, MIEEXT, MIEASY,
     $    LCBOT1, LCTOP1, CLRB1, CLRT1, TCBOT1, TCTOP1, MASEC1, MASUN1,
     $    CFRCL1, G_ASY1, NEXTO1, NSCAO1, TEMPC1, 
     $    LCBOT2, LCTOP2, CLRB2, CLRT2, TCBOT2, TCTOP2, MASEC2, MASUN2,
     $    CFRCL2, G_ASY2, NEXTO2, NSCAO2, TEMPC2
     $    )

       SUNFAC=SUNCOS*PI*(RADSUN/DISTES)**2
C      Note: PI*(RADSUN/DISTES)^2 = solid angle [steradians] of
C      the sun as seen from Earth for the case DISTES >> RADSUN.

c************************************************************************
c************************************************************************
c************************************************************************
C      ----------------------
C      Loop over the channels
C      ----------------------
       DO I=1,NCHAN

C        compute OD : indirectly uses T(z),WV(z),O3(z) through the PREDS, to get TAU, TAUZ, TAUZSN
         CALL CALC_LAYER_TRANS_YCALTODX_1_7(
     $     I,DOSUN, NCHAN, LSTCHN, FREQ, INDCHN, LBOT, 
     $     QUICKCLIST1, QUICKCLIST2, QUICKCLIST3, QUICKCLIST4, 
     $     QUICKCLIST5, QUICKCLIST6, QUICKCLIST7, 
     $     NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7, 
     $     CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, 
     $     FIXMUL, CONPRD, SUNCONPRD, TRCPRD, SUNTRCPRD, 
     $     DPRED, SUNDPRED, 
     $     INDCO2, COFCO2, INDSO2, COFSO2, INDHDO, COFHDO, 
     $     INDHNO, COFHNO, INDN2O, COFN2O, INDNH3, COFNH3,
     $     CO2MLT, SO2MLT, HNOMLT, N2OMLT, NH3MLT, HDOMLT, 
     $     INDH2O, H2OPRD, COFH2O,
     $     WAANG, LOPMIN, LOPMAX, LOPUSE, LOPLOW, DAOP, WAOP, 
     $     NFAKE, INDFAK, QUICKINDFAK, CO2TOP,
     $     COEF1, COEF2, COEF3, COEF4, COEF5, COEF6, COEF7, COEFF, 
     $     WPRED1, WPRED2, WPRED3, WPRED4, WPRED5, WPRED6, WPRED7, 
     $                         SUNWPRED4, SUNWPRED5, SUNWPRED6, SUNWPRED7, 
     $     OPRED1, OPRED2,         OPRED4, OPRED5, OPRED6, OPRED7, 
     $                         SUNOPRED4, SUNOPRED5, SUNOPRED6, SUNOPRED7, 
     $     FPRED1, FPRED2, FPRED3, FPRED4, FPRED5, FPRED6, FPRED7, 
     $                         SUNFPRED4, SUNFPRED5, SUNFPRED6, SUNFPRED7, 
     $     MPRED3, CPRED4, SUNCPRED4, SECANG, SECSUN, SUNFDG, SUNCOS,
     $       DOJAC,LISTJ,NWANTJ, CONJACPRD, DJACPRED, H2OJACPRD, DAOPJAC, 
     $       FJACPRED1,FJACPRED2,FJACPRED3,FJACPRED4,FJACPRED5,FJACPRED6,FJACPRED7,
     $       WJACPRED1,WJACPRED2,WJACPRED3,WJACPRED4,WJACPRED5,WJACPRED6,WJACPRED7,
     $       OJACPRED1,OJACPRED2,       OJACPRED4,OJACPRED5,OJACPRED6,OJACPRED7,
     $       MJACPRED3,CJACPRED4,TRCJACPRD,
     $       CO2JACMLT,SO2JACMLT,HNOJACMLT,N2OJACMLT,NH3JACMLT,HDOJACMLT,
     $       DTAU_DTZ,DTAU_DG1,DTAU_DG2,DTAU_DG3,DTAU_DG4,DTAU_DG5,DTAU_DG6,DTAU_DG9,DTAU_DG12,
     $     TAU, TAUZ, TAUZSN)

C        Calculate cloudy radiance; also no NLTE if needed
         CALL docloudyTwoSlab_RT(I, FREQ, LBOT, NWANTC, INDCHN, 
     $      TEMP,TSURF,TAU,TAUZ, TAUZSN, 
     $      BLMULT, EMIS, FCLEAR, COSDAZ, SECANG, SECSUN, DOSUN, SUNFAC, HSUN, RHOSUN, 
     $      RHOTHR, LABOVE, COEFF, LCTOP1, LCBOT1, LBLAC1, LCTOP2, LCBOT2, LBLAC2,
     $      TCBOT1, TCTOP1, TCBOT2, TCTOP2, CFRAC1, CFRAC2, CLRB1, CLRT1, CLRB2, CLRT2,
     $      CFRA12, CFRA1X, CFRA2X, 
     $      CEMIS1, CRHOS1, CRHOT1, CEMIS2, CRHOS2, CRHOT2, TEMPC1, TEMPC2,
     $      MASEC1, MASUN1, CFRCL1, G_ASY1, NEXTO1, NSCAO1, 
     $      MASEC2, MASUN2, CFRCL2, G_ASY2, NEXTO2, NSCAO2,
     $      QUICKINDNTE, NCHNTE, CLISTN, COEFN, SUNCOS, SCOS1, CO2TOP,
     $      RAD, DOJAC, TAU4, RAD4, DBTDT)

       ENDDO ! channels


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       IF (DOJAC) THEN
         CALL L2Scalc(DOJAC,DOSUN,SECANG,SECSUN,NLAY,NCHAN,TAU4,L2S4,WGT4)
       END IF

C      ----------------------------
C      Output the radiance and jacs
C      ----------------------------
       CALL WRTRTP(IPROF, IOPCO, NCHAN, RAD, PROF, NWANTC, RINDCHN)
       IF (DOJAC) THEN 
         JAC_OUTPUT_UNITS = 1
         include "writeout_jacs.f"         
       END IF

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

       IF (DOJAC) THEN
         IF (INTERSECT(  1,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) CLOSE(IOUNTZ)
         IF (INTERSECT(  3,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) CLOSE(IOUNG3)
         IF (INTERSECT(100,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) CLOSE(IOUNTZ)
         IF (INTERSECT(200,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) CLOSE(IOUNWGT)
       END IF
C
       STOP
       END
