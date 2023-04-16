c version for sarta
C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              RDRTP
C
!F77====================================================================


!ROUTINE NAME: RDRTP


!ABSTRACT:
C    Read a profile from a previously openned RTP file


!CALL PROTOCOL:
C    RDRTP( LWANT, IPROF, IOPCI, IH2O, IO3, ICO, ICH4, PTYPE, RALT,
C       NLAY, NEMIS, NRHO, LAT, LON, SATANG, SATZEN, ZSAT, SUNANG,
C       PSURF, TSURF, CO2PPM, EMIS, FEMIS, RHO, FRHO,
C       TEMP, WAMNT, OAMNT, CAMNT, MAMNT, ALT, PROF, ISTAT )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    LOGICAL   LWANT   do you want this profile?   none
C    INTEGER   IPROF   profile number              none
C    INTEGER   IOPCI   input RTP file I/O number   none
C    INTEGER   IH2O    index of H2O in gamnt       none
C    INTEGER   IO3     index of O3 in gamnt        none
C    INTEGER   ICO     index of CO in gamnt        none
C    INTEGER   ICH4    index of CH4 in gamnt       none
C    INTEGER   PTYPE   profile type code number    none
C    REAL arr  RALT    ref prof layer altitudes    meters


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   NLAY    number of used layers       none
C    INTEGER   NEMIS   number of used emis pts     none
C    INTEGER   NRHO    number of used rho pts      none
C    REAL      LAT     latitude                    degrees
C    REAL      LON     longitude                   degrees
C    REAL      SATANG  satellite scan angle        degrees
C    REAL      SATZEN  satellite zenith angle      degrees
C    REAL      ZSAT    satellite altitude          kilometer
C    REAL      SUNANG  sun zenith angle            degrees
C    REAL      PSURF   surface pressure            millibars
C    REAL      TSURF   surface skin temperature    Kelvin
C    REAL      CO2PPM  mean trop-strat CO2 mix rat PPMV
C    REAL arr  EMIS    emis points                 none (0 to 1)
C    REAL arr  FEMIS   emis freq points            cm^-1
C    REAL arr  RHO     rho points                  none (0 to 1/pi)
C    REAL arr  FRHO    rho freq points             cm^-1
C    REAL arr  TEMP    layer temperature           Kelvin
C    REAL arr  WAMNT   layer Water vapor amount    kilomoles/cm^2
C    REAL arr  OAMNT   layer Ozone amount          kilomoles/cm^2
C    REAL arr  CAMNT   layer CO amount             kilomoles/cm^2
C    REAL arr  MAMNT   layer Methane amount        kilomoles/cm^2
ccc
cC    REAL arr  WAMNT   layer Water vapor amount    molecules/cm^2
cC    REAL arr  OAMNT   layer Ozone amount          molecules/cm^2
cC    REAL arr  CAMNT   layer CO amount             molecules/cm^2
cC    REAL arr  MAMNT   layer Methane amount        molecules/cm^2
ccc
C    REAL arr  ALT     layer average altitude      meters
C    STRUCT    PROF    RTP profile structure       various
C    INTEGER   ISTAT   I/O status                  none


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): KLAYERS


!ROUTINES CALLED: none


!FILES ACCESSED:
C    Input RTP file withI/O number IOPCI
C    unit IOERR: error messages
C    unit IOINFO: info/warning messages


!COMMON BLOCKS: none


!DESCRIPTION:
C    Reads a single profile from a previously openned RTP file.
C    The routine expects to find the data specified in the header
C    of the input RTP file.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 14 Feb 2001 Scott Hannon      created based on klayers version
C 13 Sep 2001 Scott Hannon      Added checks of PSURF & TSURF
C 31 Oct 2002 Scott Hannon      Add output vars SATZEN and ZSAT
C 06 Oct 2004 Scott Hannon      Add error trap for LAT
C 17 Dec 2004 Scott Hannon      Add PTYPE to call; fix error in NLAY
C                                  when PTYPE=AIRSLAY


!END====================================================================

C      =================================================================
       SUBROUTINE RDRTP(LWANT, IPROF, IOPCI, IH2O, IO3, ICO, ICH4,
     $    PTYPE, RALT, NLAY, NEMIS, NRHO, LAT, LON, SATANG, SATZEN,
     $    ZSAT, SUNANG, PSURF, TSURF, CO2PPM, EMIS, FEMIS, RHO, FRHO,
     $    TEMP, WAMNT, OAMNT, CAMNT, MAMNT, ALT, PROF, ISTAT )
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
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input parameters:
       LOGICAL LWANT      ! do we want this profile?
       INTEGER IPROF      ! number of current profile
       INTEGER IOPCI      ! input RTP unit
       INTEGER IH2O       ! index of H2O in gamnt
       INTEGER IO3        ! index of O3 in gamnt
       INTEGER ICO        ! index of CO in gamnt
       INTEGER ICH4       ! index of CH4 in gamnt
       INTEGER PTYPE      ! profile type code number
       REAL RALT(MAXLAY)  ! ref prof layer average altitudes


C      Output parameters:
       INTEGER   NLAY
       INTEGER  NEMIS
       INTEGER   NRHO
       REAL    LAT
       REAL    LON
       REAL SATANG
       REAL SATZEN
       REAL   ZSAT
       REAL SUNANG
       REAL  PSURF
       REAL  TSURF
       REAL CO2PPM
       REAL   EMIS(MXEMIS)
       REAL  FEMIS(MXEMIS)
       REAL    RHO(MXEMIS)
       REAL   FRHO(MXEMIS)
       REAL   TEMP(MAXLAY)
       REAL  WAMNT(MAXLAY)
       REAL  OAMNT(MAXLAY)
       REAL  CAMNT(MAXLAY)
       REAL  MAMNT(MAXLAY)
       REAL    ALT(MAXLAY)
C
C      Profile data structure
       RECORD /RTPPROF/ PROF
C
       INTEGER ISTAT


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER I        ! generic
       INTEGER L        ! layer looping
       INTEGER LR       ! reversed layer looping
       INTEGER NLEV     ! number of levels
       INTEGER rtpread  ! for calling read rtp interface routine
       REAL ZSURF       ! surface altitude (read but ignored for now)
       REAL RJUNK1      ! generic junk/work
       REAL RJUNK2      ! generic junk/work

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************

C      ------------------------
C      Read the current profile
C      ------------------------
       ISTAT=rtpread(IOPCI, PROF)
C
       IF (ISTAT .EQ. -1) GOTO 9999  ! reached end of file
C
       IF (.NOT. LWANT) GOTO 9999    ! skip prof if not wanted
C
C      --------------------
C      Pull out needed data
C      --------------------
C      Latitude & longitude
       LAT=PROF.plat
       LON=PROF.plon
C
C      Can not correctly compute gravity if LAT is bad
       IF (ABS(LAT) .GT. 90.01) THEN
          WRITE(IOERR,1005) IPROF, LAT
 1005     FORMAT('ERROR! input profile PROF(',I4,').plat=',1PE11.4,
     $    ' is out of range -90 to +90')
          STOP
       ENDIF
C      Note: LON is not currently used and thus need not be checked
C
C      Number of levels
       NLEV=PROF.nlevs
       IF (NLEV .LT. 2 .OR. NLEV .GT. MAXLAY+1) THEN
          WRITE(IOERR,1010) IPROF, NLEV, MAXLAY+1
 1010     FORMAT('ERROR! input profile PROF(',I4,').nlevs=',I4,
     $    ' is out of range 2 to ',I3)
          STOP
       ENDIF

C      Number of layers
       IF (PTYPE .EQ. AIRSLAY) THEN
C         Special case for AIRS pseudo-levels
          NLAY=NLEV
       ELSE
          NLAY=NLEV - 1
       ENDIF
C
C      CO2 mixing ratio
       CO2PPM=PROF.co2ppm
       IF (CO2PPM .LT. -998) CO2PPM=CO2STD
       RJUNK1=0.8*CO2PPM
       RJUNK2=1.2*CO2PPM
       IF (CO2PPM .LT. RJUNK1 .OR. CO2PPM .GT. RJUNK2) THEN
          WRITE(IOERR,1015) IPROF, CO2PPM, RJUNK1, RJUNK2
 1015     FORMAT('Error! PROF(',I4,').co2ppm=',1PE10.3,
     $    ' is outside allowed range ',0PF5.1,' to ',F5.1)
       ENDIF
C
C      Angles
       SUNANG=PROF.solzen
       SATANG=PROF.scanang
       SATZEN=PROF.satzen
C
C      Satellite altitude above ellipsoid surface (convert m to km)
       ZSAT=PROF.zobs/1000
C
C      Surface
       PSURF=PROF.spres
       TSURF=PROF.stemp
       ZSURF=PROF.salti  ! note: ZSURF is currently ignored
       IF (PSURF .LE. 0) THEN
          WRITE(IOERR,1017) IPROF
 1017     FORMAT('ERROR! Prof(',I4,') has no surface pressure')
          STOP
       ENDIF
       IF (TSURF .LE. 0) THEN
          WRITE(IOERR,1018) IPROF
 1018     FORMAT('ERROR! Prof(',I4,') has no surface temperature')
          STOP
       ENDIF
C
C      Emissivity (range 0 to 1)
       NEMIS=PROF.nemis
       IF (NEMIS .EQ. 0) THEN
          WRITE(IOERR,1020) IPROF
 1020     FORMAT('ERROR! PROF(',I4,') has no emissivity')
          STOP
       ENDIF
       DO I=1,NEMIS
          EMIS(I)=PROF.emis(I)
          FEMIS(I)=PROF.efreq(I)
       ENDDO
C
C      Reflectivity (range 0 to 1/pi)
       NRHO=PROF.nrho
       DO I=1,NRHO
          RHO(I)=PROF.rho(I)
          FRHO(I)=PROF.rfreq(I)
       ENDDO
       IF (NRHO .EQ. 0) THEN
          WRITE(IOINFO,1022) IPROF
 1022     FORMAT('Note! PROF(',I4,') has no reflectivity',
     $    '; will use (1-emis)/pi')
          DO I=1,NEMIS
             RHO(I)=(1 - EMIS(I))/PI
             FRHO(I)=FEMIS(I)
          ENDDO
       ENDIF
C
C      ----------------------------------
C      Get layer temperature & gas amount
C      ----------------------------------
       IF (GUCIN .EQ. 1) THEN
C         Input gas units are molecules/cm^2; convert to kilomoles/cm^2
          IF (PROF.plevs(1) .LT. PROF.plevs(NLEV)) THEN
C            Prof is in top-down order
             DO L=1,NLAY
                TEMP(L)=PROF.ptemp(L)
                WAMNT(L)=PROF.gamnt(L,IH2O)/6.02214199E+26
                OAMNT(L)=PROF.gamnt(L,IO3 )/6.02214199E+26
                CAMNT(L)=PROF.gamnt(L,ICO )/6.02214199E+26
                MAMNT(L)=PROF.gamnt(L,ICH4)/6.02214199E+26
                ALT(L)=0.5*( PROF.palts(L) + PROF.palts(L+1) )
             ENDDO
          ELSE
C            Prof is in bottom-up order
             DO L=1,NLAY
                LR=1 + NLAY - L  ! reversed layer index
                TEMP(L)=PROF.ptemp(LR)
                WAMNT(L)=PROF.gamnt(LR,IH2O)/6.02214199E+26
                OAMNT(L)=PROF.gamnt(LR,IO3 )/6.02214199E+26
                CAMNT(L)=PROF.gamnt(LR,ICO )/6.02214199E+26
                MAMNT(L)=PROF.gamnt(LR,ICH4)/6.02214199E+26
                ALT(L)=0.5*( PROF.palts(LR) + PROF.palts(LR+1) )            
             ENDDO
          ENDIF
C
       ELSEIF (GUCIN .EQ. 2) THEN
C         Input gas units are kilomoles/cm^2
          IF (PROF.plevs(1) .LT. PROF.plevs(NLEV)) THEN
C            Prof is in top-down order
             DO L=1,NLAY
                TEMP(L)=PROF.ptemp(L)
                WAMNT(L)=PROF.gamnt(L,IH2O)
                OAMNT(L)=PROF.gamnt(L,IO3)
                CAMNT(L)=PROF.gamnt(L,ICO)
                MAMNT(L)=PROF.gamnt(L,ICH4)
                ALT(L)=0.5*( PROF.palts(L) + PROF.palts(L+1) )
             ENDDO
          ELSE
C            Prof is in bottom-up order
             DO L=1,NLAY
                LR=1 + NLAY - L  ! reversed layer index
                TEMP(L)=PROF.ptemp(LR)
                WAMNT(L)=PROF.gamnt(LR,IH2O)
                OAMNT(L)=PROF.gamnt(LR,IO3)
                CAMNT(L)=PROF.gamnt(LR,ICO)
                MAMNT(L)=PROF.gamnt(LR,ICH4)
                ALT(L)=0.5*( PROF.palts(LR) + PROF.palts(LR+1) )            
             ENDDO
          ENDIF
       ELSE
       ENDIF
C

C      -----------------
C      Default altitudes
C      -----------------
       IF (PROF.palts(1) .LT. -998 .OR. PTYPE .EQ. AIRSLAY) THEN
C         No altitudes, use reference
          DO L=1,NLAY
             ALT(L)=RALT(L)
          ENDDO
       ENDIF
c
c for now ignore clouds
c

C
 9999  RETURN
       END
