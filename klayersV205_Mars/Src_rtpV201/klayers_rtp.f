c This version for use with RTP files
C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              KLAYERS
C
!F77====================================================================


!ROUTINE NAME: KLAYERS


!ABSTRACT: Profile integration/averaging program for use with kCARTA


!CALL PROTOCOL: main program


!INPUT PARAMETERS: none


!OUTPUT PARAMETERS: none


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: main program


!PARENT(S): main program


!ROUTINES CALLED:
C    INTEG  : integrates sublayers to make averaged layers.
C    INTLEV : interpolates a profile onto sublevels.
C    RDINFO : reads info about a user supplied profile
C    LAYOUT : writes output in format required for use with kCARTA.
C    MERGE  : combine a user and AFGL profile.
C    SPLLEV : spline interpolation of profile levels.
C    WRTSPL : writes spline interpolated profile to a output file.
C    RDUSER : read a user profile in the LAYERS user profile format.
C    RDAFGL : read a AFGL model profile in the LAYERS AFGL format.


!FILES ACCESSED:
C    unit 5: standard input/keyboard (filenames, yes/no, etc)
C    IOERR : error messages
C    IOOUT : unit 30 : "fort.30" main output file
C    IOIN  : unit  9 : profile input
C    IOSPL : unit 11 : spline output file


!COMMON BLOCKS:
C    COMLEV : file "cbplev.f"; assigns values to array PLEV


!DESCRIPTION:
C    Driver program for the routines.
C
C    RDINFO: reads in profile info (filename, various yes/no, etc)
C
C    RDUSER: reads in the user supplied profile
C
C    RDAFGL: reads in an AFGL model.
C
C    MERGE: merges the user and AFGL profiles.
C
C    Adjust the profile mixing ratios for "wet air" conditions if
C    instructed.
C
C    INTLEV: interpolate the profile onto a set of sub-levels.
C
C    INTEG: integrate the sub-levels to get the averaged and
C       integrated layer values.
C
C    LAYOUT: write the integrated/averaged profile to output (fort.30)


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    None


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------


!END====================================================================

C      =================================================================
       PROGRAM KLAYERS
C      =================================================================


C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
      include 'incLAY.f'
      include 'rtpdefs.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      none (main program)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER   IOIN
       INTEGER  IOOUT
       INTEGER      I
       INTEGER      J
       INTEGER      L
       INTEGER    NIN
       INTEGER   NSUB(MYNLAY)
       INTEGER  NFINE
       INTEGER     JO
       INTEGER  GASID( MXGAS)
       INTEGER NGASES
       INTEGER   IMAX
       INTEGER   IFIX
       INTEGER   IWAT
       INTEGER  NWANT
       INTEGER  LISTG( MXGAS)
       INTEGER  NAFGL
       INTEGER IDAFGL( MXGAS)
       INTEGER NLEVAF
       INTEGER MNAFGL
       INTEGER GORDER( MXGAS)
       INTEGER   NXIN
       INTEGER  IOSPL
       INTEGER   NINX
       INTEGER     IP              ! profile loop counter
       INTEGER LAYBOT              ! bottom layer number
       REAL  PSURF                 ! surface pressure
       REAL  ZSURF                 ! surface altitude
       REAL     PB(MYNLAY+1)
       REAL   LNPB(MYNLAY+1)
       REAL    LAT
       REAL    LON
       REAL  ZPMAX
       REAL    PIN(  MXIN)
       REAL  LNPIN(  MXIN)
       REAL    TIN(  MXIN)
       REAL   MRIN(  MXIN, MXGAS)
       REAL  PFINE(NSUBLV)
       REAL  TFINE(NSUBLV)
       REAL  TGLAY(MYNLAY, MXGAS)
       REAL    ZIN(  MXIN)
       REAL MRFINE(NSUBLV, MXGAS)
       REAL  TGSUB(NSUBLV-1, MXGAS)
       REAL  DZLAY(MYNLAY)
       REAL  ZFINE(NSUBLV)
       REAL   PSUB(NSUBLV-1)
       REAL   TSUB(NSUBLV-1)
       REAL   ZSUM
       REAL   PLAY(MYNLAY)
       REAL   TLAY(MYNLAY)
       REAL   ALAY(MYNLAY, MXGAS)
       REAL   TOFF
       REAL  MRSUB(NSUBLV-1, MXGAS)
       REAL  DZSUB(NSUBLV-1)
       REAL   ZLAY(MYNLAY)
       REAL  LATAF
       REAL    HAF(  MXIN)
       REAL    PAF(  MXIN)
       REAL  LNPAF(  MXIN)
       REAL    TAF(  MXIN)
       REAL    DAF(  MXIN)
       REAL  MIXAF(  MXIN, MXGAS)
       REAL HBOUND(MYNLAY+1)
       REAL  PARTP( MXGAS,MYNLAY)
       REAL   DELT(MYNLAY)
       REAL LNPXIN(  MXIN)
       REAL   PXIN(  MXIN)
       REAL   TXIN(  MXIN)
       REAL  MRXIN(  MXIN, MXGAS)
       REAL   PMIN
       REAL   PMAX
       REAL  MASSW( MXGAS)  ! mass of all wanted gases
       REAL CO2MLT          ! CO2 mixing ratio multiplier
C
       CHARACTER*70 FNAFGL
C
       LOGICAL LZ
       LOGICAL LSVP
       LOGICAL LTAIR
       LOGICAL LDRY
       LOGICAL LSPLIN
       LOGICAL LPMAX
       LOGICAL LAFGL
C
C      for rdinfo
       REAL SCALEH             ! scale height for default ZSURF calc
       CHARACTER*80 FIN        ! input RTP filename
       CHARACTER*80 FOUT       ! output RTP filename
       CHARACTER*256 COMMNT    ! klayers info comment string
C
C      for opnrtp
       INTEGER  NGASI          ! number of gases in input profile
       INTEGER GLISTI( MXGAS)  ! list of all gases in input file
       INTEGER  NGASF          ! # of wanted gases found in input file
       INTEGER GLISTF( MXGAS)  ! list of found wanted gases
       INTEGER GUNITF( MXGAS)  ! gas units code#'s of found gases
       INTEGER INDEXF( MXGAS)  ! indices in RTP of found wanted gases
       REAL  MASSF( MXGAS)     ! mass of found wanted gases
       INTEGER IOPCI           ! input RTP file I/O channel
       INTEGER IOPCO           ! output RTP file I/O channel
C
C      for rdrtp, wrtrtp, & rtpclose
       INTEGER ISTAT
       INTEGER rtpclose
C
C      Boundary pressure levels
       COMMON /COMLEV/ NAMGRD, PLEV
       CHARACTER*40 NAMGRD
       REAL PLEV(NBPLEV)
C
C      Allowed gas IDs
       COMMON /COMGID/ GIDS, GMASS
       INTEGER GIDS(NGIDS)
       REAL GMASS(NGIDS)
C
       RECORD /RTPPROF/ PROF

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
C      ------------------------------------
C      Assign input and output unit numbers
C      ------------------------------------
C      Input files
       IOIN=9
C      Output file
       IOOUT=30
C      Output file for spline interpolated profiles
       IOSPL=11
C

C      ----------------------------------
C      One air temperature for all gases?
C      ----------------------------------
C      Use air temperature rather than individual gas temperature T/F
       LTAIR=.TRUE.
C

C      ----------------------------
C      Assign PB and LNPB from PLEV
C      ----------------------------
C      Min and max pressure levels (see "cbplev.f")
       PMIN=PLEV(NBPLEV)
       PMAX=PLEV(1)
C
       DO I=1,NBPLEV
          PB(I)=PLEV(I)
          LNPB(I)=LOG(PB(I))
C
C         Assign maximum index in PLEV for sat vap pres check
C         SVP equations seem to break-down below ~100 mb
          IF (PLEV(I) .GT. 2.0) IMAX=I
       ENDDO
C

C      -----------------------------------------------
C      Get filenames & other program control variables
C      -----------------------------------------------
       CALL RDINFO(FIN, FOUT, COMMNT, FNAFGL, MNAFGL, NWANT, LISTG,
     $   MASSW, TOFF, SCALEH, LDRY, LSVP, LSPLIN, LPMAX)
C      note: new rtp version of rdinfo; parameter list is slightly
C      different than old rdinfo
C
       IF (LSPLIN) THEN
C         Open spline interpolation output file
          OPEN(UNIT=IOSPL,FILE='spline.out',STATUS='NEW',
     $       FORM='FORMATTED')
       ENDIF
C

C      ------------------------
C      Read in the AFGL profile 
C      ------------------------
       CALL RDAFGL(FNAFGL, IOIN, MNAFGL, NWANT, LISTG, NLEVAF, NAFGL,
     $    IDAFGL, LATAF, HAF, PAF, TAF, DAF, MIXAF)
C
       DO I=1,NLEVAF
          LNPAF(I)=LOG(PAF(I))
       ENDDO
C

C      -----------------------------
C      Open RTP input & output files
C      -----------------------------
       CALL OPNRTP(FIN, FOUT, COMMNT, NWANT, LISTG, MASSW,
     $    PMIN, PMAX, NGASI, GLISTI, NGASF, GLISTF, GUNITF, MASSF,
     $    INDEXF, IOPCI, IOPCO)
C

       IF (LSPLIN) THEN
C         ---------------------------------------
C         Do spline interpolation of AFGL profile
C         ---------------------------------------
C         Load up copy of input data
          NXIN=NLEVAF
          DO I=1,NLEVAF
             PXIN(I)=PAF(I)
             LNPXIN(I)=LNPAF(I)
             TXIN(I)=TAF(I)
             DO J=1,NAFGL
                MRXIN(I,J)=MIXAF(I,J)
             ENDDO
          ENDDO
C
C         Call spline interpolation subroutine
C         Note: the call is made so that the spline interpolated
C            data comes back in the NLEVAF/PAF/TAF/MIXAF variables
          CALL SPLLEV(NXIN,LNPXIN,PXIN,TXIN,NAFGL,MRXIN,
     $       NLEVAF,LNPAF, PAF, TAF, MIXAF, PMIN, PMAX)
C
C         Write the interpolated profile to the spline out file
          LAFGL=.TRUE.
          CALL WRTSPL(IOSPL, FNAFGL, NXIN, PXIN, TXIN, NAFGL, IDAFGL,
     $       MRXIN, NLEVAF, PAF, TAF, MIXAF, LAFGL, MNAFGL)
       ENDIF

C      -----------------------------------------------------------------
C      Start of while loop over profiles
       IP=0
C      ------------------------
C      Read the current profile
C      ------------------------
 10    CALL RDRTP(ISTAT, IOPCI, NGASF, GLISTF, INDEXF, LAT, LON,
     $    NIN, PIN, LNPIN, TIN, MRIN, ZIN, LZ, PSURF, ZSURF, PROF)

c       DO J=1,NIN
c         print *,'just after rdrtp',J,PIN(J),TIN(J)
c       END DO

C
ccc
c for testing only
c      LON=0.0
ccc
       IF (ISTAT .EQ. -1) GOTO 99  ! reached End Of File
C
       IP=IP + 1  ! increment profile counter
C
       IF (LPMAX) THEN
C         Force profile to extend to PMAX
          LAYBOT=1

C         Calc Z at Pmax using kludged adjustment to scale height
          ZPMAX=-1E+3*(SCALEH+LOG(PMAX+1)/5)*LOG(PMAX/10.00)  !! use 10.00 insted of 1013.25
C
C         Adjust ZPMAX if profile has any altitude info
          IF (LZ) THEN
             IF (PSURF .GT. PIN(1) .AND. ZSURF .GT. -999.89) THEN
C               Adjust ZPMAX based on PSURF & ZSURF
                ZSUM=ZSURF +
     $             1E+3*(SCALEH+LOG(PSURF+1)/5)*LOG(PSURF/10.00) !! use 10.00 insted of 1013.25
             ELSE
C               Adjust ZPMAX based on PIN(1) & ZIN(1)
                ZSUM=ZIN(1) +
     $             1E+3*(SCALEH+LOG(PIN(1)+1)/5)*LOG(PIN(1)/10.00) !! use 10.00 insted of 1013.25
             ENDIF
             ZPMAX=ZPMAX + ZSUM
          ELSE
             IF (PSURF .GT. 0.0 .AND. ZSURF .GT. -998.9) THEN
C               Calc offset between specified Zsurf and default value
                ZSUM=ZSURF +
     $             1E+3*(SCALEH+LOG(PSURF+1)/5)*LOG(PSURF/10.00) !! use 10.00 insted of 1013.25
                ZPMAX=ZPMAX + ZSUM
             ENDIF
          ENDIF
C             
          PSURF=PMAX
          ZSURF=ZPMAX
C
       ELSE
C         Check PSURF and calc LAYBOT
          IF (PSURF .LT. PLEV(MYNLAY)) THEN
             WRITE(IOERR,1020) IP, PSURF
 1020        FORMAT('ERROR! profile ',I5,' has bad surface pressure=',
     $          F10.4)
          ELSE
C            Calc bottom layer number counting upward (ie decreasing P)
             LAYBOT=1
 20          IF (PLEV(LAYBOT+1) .GE. PSURF) THEN
                LAYBOT=LAYBOT+1
                GOTO 20
             ENDIF
          ENDIF
C
C         Check Zsurf; calc default if nodata
          IF (ZSURF .LT. -998.9) THEN
C            Zsurf is probably -999 or -9999 nodate flag; calc default
             ZSURF=-1E+3*(SCALEH+LOG(PSURF+1)/5)*LOG(PSURF/10.00) !! use 10.00 insted of 1013.25
C
             IF (LZ) THEN
C               Find highest user level at or below ZSURF
                J=1
 30             IF (PIN(J) .GT. PSURF .AND. J .LT. NIN) THEN
                   J=J + 1
                   GOTO 30
                ENDIF
                ZSUM=ZIN(J) +
     $             1E+3*(SCALEH+LOG(PIN(J)+1)/5)*LOG(PIN(J)/10.00) !! use 10.00 insted of 1013.25
                ZSURF=ZSURF + ZSUM
             ENDIF
C
             WRITE(IOINFO,1030) IP, ZSURF, PSURF
 1030        FORMAT('Note: profile ',I5,' using default calc zsurf=',
     $       F7.1,'m at psurf=',F8.3,'mb')
          ENDIF
       ENDIF
C
C      Reset GASID; it gets re-written every loop by MERGE
       DO I=1,NGASF
          GASID(I)=GLISTF(I)
       ENDDO

C      --------------------------------------
C      If necessary, convert gas amount units
C      --------------------------------------
       CALL TOPPMV(NGASF, GLISTF, GUNITF, NIN, MASSF, PIN, TIN, MRIN)
C

C      ---------------------------------------
C      Do spline interpolation of user profile
C      ---------------------------------------
       IF (LSPLIN .AND. .NOT. LZ) THEN
C         Load up copy of input data
          NXIN=NIN
          DO I=1,NIN
             PXIN(I)=PIN(I)
             LNPXIN(I)=LNPIN(I)
             TXIN(I)=TIN(I)
             DO J=1,NGASF
                MRXIN(I,J)=MRIN(I,J)
             ENDDO
          ENDDO
C
C         Call spline interpolation subroutine
C         Note: the call is made so that the spline interpolated
C            data comes back in the NIN/PIN/TIN/MRIN variables
          CALL SPLLEV(NXIN,LNPXIN,PXIN,TXIN,NGASF,MRXIN,
     $       NIN,LNPIN, PIN, TIN, MRIN, PMIN, PMAX)
C
C         Write the interpolated profile to the spline out file
          LAFGL=.FALSE.
          CALL WRTSPL(IOSPL, FIN, NXIN, PXIN, TXIN, NGASF, GLISTF,
     $       MRXIN, NIN, PIN, TIN, MRIN, LAFGL, MNAFGL)
C
       ENDIF
C

C      --------------
C      Merge profiles
C      --------------
       CO2MLT=1.0
       IF (PROF.co2ppm .GT. 0 .AND. PROF.co2ppm .LE. 1E+6) THEN
          CO2MLT=PROF.co2ppm/CO2STD
       ENDIF
C
       CALL MERGE(LZ, NGASES, PMIN, NINX, LISTG, GORDER,
     $    NIN,    ZIN, PIN, LNPIN, NGASF, GASID,  TIN, MRIN,
     $    NLEVAF, HAF, PAF, LNPAF, NAFGL, IDAFGL, TAF, MIXAF, CO2MLT)
C
C      Check the number of gases
       IF (NGASES .NE. NWANT) THEN
          WRITE(IOERR,1010) NWANT, NGASES
 1010     FORMAT('ERROR! Unexpected number of gases after merging',
     $       ' the input and AFGL',/,'profiles.  NWANT=',I2,
     $       ' but after merging have NGASES=',I2)
          WRITE(IOERR,1011) (LISTG(I), I=1,NWANT)
 1011     FORMAT('Wanted gases LISTG=',44(1X,I2))
          WRITE(IOERR,1012) (GASID(I), I=1,NGASES)
 1012     FORMAT('Merged gases GASID=',44(1X,I2))
          STOP
       ENDIF
C
C      Note on MERGE: on return NGASES should be the same as
C      NWANT.  On sending GASID should be the same as GLISTF,
C      but on return GASID will be extended to length NGASES
C      by adding on AFGL model values for the gases not found
C      in the input file.  GORDER is the order of gases in
C      the returned GASID relative to LISTG.

C      ------------------
C      Get index of water
C      ------------------
       IWAT=0
       DO J=1,NGASES
          IF (GASID(J) .EQ. 1) IWAT=J
       ENDDO
C
C      -----------------------------------------------------
C      Adjust dry air mixing ratios to wet air mixing ratios
C      -----------------------------------------------------
C      First do water
       IF (( MRTH2O .EQ. 2) .OR. (LDRY .AND. MRTH2O .EQ. 0)) THEN
C         Water is dry or other gases are dry and water is same
          DO I=1,NINX
             MRIN(I,IWAT)=1.0E+6*MRIN(I,IWAT)/(MRIN(I,IWAT)+ 1.0E+6)
          ENDDO
       ENDIF
C
C      Now do all the other gases
       IF (LDRY) THEN
          DO J=1,NGASES
             IF (J .NE. IWAT) THEN
                DO I=1,NINX
                   MRIN(I,J)=MRIN(I,J)*(1.0E+6 - MRIN(I,IWAT))/
     $             1.0E+6
                ENDDO
             ENDIF
          ENDDO
       ENDIF
C
C      --------------------------------------------------
C      Interpolate the input profile onto fine sub-levels
C      -------------------------------------------------
       CALL INTLEV(IP, LZ, LSVP, IMAX, NINX, NGASES, PIN, LNPIN, TIN,
     $    MRIN, ZIN, LAT, LON, PB, LNPB, LAYBOT, PSURF, ZSURF,
     $    NSUB, NFINE, PFINE, TFINE,
     $    MRFINE, ZFINE, IWAT, PSUB, TSUB, TGSUB, MRSUB, DZSUB, ZPMAX)
C

ccccccccccc
c for testing
c       do J=1,NFINE
c       print *, J, PFINE(J), TFINE(J), ZFINE(J), MRFINE(J,IWAT)
c       enddo
ccccccccccc

C      -------------------------------------
C      Calculate layer heights & thicknesses
C      -------------------------------------
       HBOUND(LAYBOT)=ZPMAX
       J=0
       DO I=LAYBOT,MYNLAY
          JO=J
          J=J+NSUB(I)
          HBOUND(I+1)=ZFINE(J+1)
          DELT(I)=-(TFINE(J+1)-TFINE(JO+1))
       ENDDO
C
C      ----------------------------------------------------
C      Integrate the sub-levels to form the averaged layers
C      ----------------------------------------------------
       CALL INTEG(LTAIR, LAYBOT, NSUB, NFINE, NGASES,
     $    PFINE, PSUB, TSUB,
     $    TGSUB, MRSUB, DZSUB, PLAY, TLAY, TGLAY, ALAY, DZLAY)
C
C      -----------------------------------
C      Calculate the gas partial pressures
C      -----------------------------------
       ZSUM=ZPMAX
       DO L=LAYBOT,MYNLAY
          DO I=1,NGASES
C            denav=1.2027E-12 (see integ.f)
             PARTP(I,L)=ALAY(L,I)*TLAY(L)/
     $          (DZLAY(L)*1.2027E-12*1.0E+6*1013.25)
          ENDDO
          ZLAY(L)=ZSUM + (DZLAY(L)/2.0)
          ZSUM=ZSUM + DZLAY(L)
       ENDDO
C

C      --------------------------------------------
C      Set user-to-AFGL profile cross-over pressure
C      --------------------------------------------
       CALL SETXOP(NIN, NGASF, NWANT, GLISTF, LISTG, PSURF, PIN, PROF)


C      --------------------------------
C      Write profile to RTP output file
C      --------------------------------
       CALL WRTRTP(IP, ISTAT, IOPCO, LAYBOT, NGASES, LISTG, GORDER,
     $    HBOUND, PLEV, PLAY, TLAY, ALAY, TOFF, PROF)
C
       IF (ISTAT .EQ. -1) GOTO 99  ! reached End Of File on output!
C

C      --------------------
C      Loop to next profile
C      --------------------
       GOTO 10

 99    CONTINUE  ! end of while loop over profiles
C      -----------------------------------------------------------------

C      -----------
C      Close files
C      -----------
C      Close RTP files
       ISTAT=rtpclose(IOPCI)
       ISTAT=rtpclose(IOPCO)
C
C      Close spline file
       IF (LSPLIN) CLOSE(IOSPL)

C
       STOP
       END
