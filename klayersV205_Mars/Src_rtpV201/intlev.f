C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              INTLEV
C
!F77====================================================================

!ROUTINE NAME: INTLEV

!ABSTRACT:
C    Interpolate a profile onto a set of AIRS sub-levels.

!CALL PROTOCOL:
C    INTLEV( IP, LZ, LSVP, IMAX, NIN, NGASES, PIN, LNPIN, TIN,
C          MRIN, ZIN, LAT, LON, PB, LNPB, LAYBOT, PSURF, ZSURF, 
C          NSUB, NFINE, PFINE, TFINE, MRFINE, ZFINE, WATID,
C          PSUB, TSUB, TGSUB, MRSUB, DZSUB, ZPMAX )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   IP      profile number              none
C    INTEGER   IMAX    Max index of PB for SVP     none
C    REAL      LAT     Latitude                    degrees
C    INTEGER   LAYBOT  bottom layer number         none
C    REAL arr  LNPB    Log(AIRS pressure levels)   PB in mb
C    REAL arr  LNPIN   Log(profile pressure levs)  Pin in mb
C    REAL      LON     Longitude                   degrees
C    LOGICAL   LSVP    Check H2O sat vapor pres?   none
C    LOGICAL   LZ      Prof altitudes supplied?    none
C    REAL arr  MRIN    Profile Mixing Ratios       ppmv
C    INTEGER   NGASES  Number of gases             none
C    INTEGER   NIN     Number of profile levels    none
C    REAL arr  PB      AIRS pressure levels        mb
C    REAL arr  PIN     Profile pressure levels     mb
C    REAL      PSURF   surface pressure            mb
C    REAL arr  TIN     Profile temperatures        K
C    INTEGER   WATID   Water index in MRIN         none
C    REAL arr  ZIN     Profile altitudes           m
C    REAL      ZSURF   surface altitude            m


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  DZSUB   Sub-layer thickness         m
C    REAL arr  MRFINE  Sub-level mixing ratios     ppmv
C    REAL arr  MRSUB   Sub-layer avg mixing ratio  ppmv
C    INTEGER   NFINE   Number of sub-levels        none
C    INT arr   NSUB    # of sub-layers per layer   none
C    REAL arr  PFINE   Sub-level pressures         mb
C    REAL arr  PSUB    Sub-layer average presure   mb
C    REAL arr  TFINE   Sub-level temperatures      K
C    REAL arr  TGSUB   Sub-layer gas avg temp      K
C    REAL arr  TSUB    Sub-layer air avg temp      K
C    REAL arr  ZFINE   Sub-level altitudes         m
C    REAL      ZPMAX   bottom-most used Plev alt   m

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S): KLAYERS


!ROUTINES CALLED:
C    GRAV (function)  : computes Earth's gravity
C    WEXSVP (function): computes saturation vapor pressure
C    ZBOT (function) : computes altitude at bottom of bottom layer


!FILES ACCESSED: none

!COMMON BLOCKS: none


!DESCRIPTION:
C    Interpolates an input profile onto a set of sub-levels. The
C    sub-levels are the layer boundary levels divided into finer
C    sub-layers. The layer sub-division is adjusted as needed to ensure
C    that each input profile data point (in the range covered by the
C    layers) is accounted for in the sub-levels.
C
C    The interpolations of Temperature and Mixing Ratio are done by
C    assuming a relationship linear with the log of Pressure.
C
C    If profile altitudes are supplied along with the pressure,
C    temperature, and mixing ratios, this is also interpolated assuming
C    a linear in ln(P) dependence.
C
C    If profile altitudes are not supplied, they are calculated using
C    a finite slab approximation of the hydrostatic equation,
C       delta_z = - delta_P * Tavg /( gavg * Cavg * Pavg)
C    where Cavg is the average mass of one mole of "air" divided by
C    the gas constat R, and Tavg, gavg, and Pavg are the average
C    temperature, gravity, and pressure. Since this equation supplies
C    only layer thickness (delta_z), it is necessary to supply a
C    surface altitude as well as surface pressure, and this is then
C    used to calc the altitude at the bottom of the the bottom-most
C    AIRS pressure level.
C 
C    Since "air" mass-per-mole varies primarily due to variations in
C    the water mixing ratio, it is necessary to know which of the
C    NGASES mixing ratios supplied corresponds to water (WATID).
C
C    If desired, the routine will check the water profile to see if
C    it is over-saturated and truncate it to 100% humidity.


!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS:
C    Assumes the input profile is reasonable and can be interpolated
C    to span the entire pressure range of the layers.


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C Apr  4 1995 Scott Hannon/UMBC created
C Apr 26 1995 Scott Hannon      fixed array(I-1) for case I=1
C May  8 1995 Scott Hannon      modified code for case WATID=0
C Jun 23 1995 Scott Hannon      Correct some comments
C Sep 12 1995 Scott Hannon      Correct math error in dry air CFINE
C 21 Feb 2001 Scott Hannon      Warning messages now use unit IOINFO
C 27 Mar 2001 Scott Hannon      Add IP,LAYBOT,PSURF,ZSURF to call;
C                               replace input Z0 with output ZPMAX;
C                               add local GSURF; change loops from
C                               1 to MYNLAY to LAYBOT to MYNLAY; add
C                               call to function ZB0T; revise warning
C                               messages.
C 15 May 2001 Scott Hannon      Remove RJUNK2 from ZBOT call.
C 30 Aug 2001 Scott Hannon      Reverse "small"/"large" in format 1030;
C                               Add SATERR to over-saturation warning


!END====================================================================

C      =================================================================
       SUBROUTINE INTLEV(IP, LZ, LSVP, IMAX, NIN, NGASES, PIN, LNPIN,
     $    TIN, MRIN, ZIN, LAT, LON, PB, LNPB, LAYBOT, PSURF, ZSURF,
     $    NSUB, NFINE, PFINE, TFINE,
     $    MRFINE, ZFINE, WATID, PSUB, TSUB, TGSUB, MRSUB, DZSUB, ZPMAX)
C      =================================================================
C
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
      include 'incLAY.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
       REAL GRAV
       REAL WEXSVP
       REAL ZBOT


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input:
       INTEGER     IP
       LOGICAL     LZ
       LOGICAL   LSVP
       INTEGER   IMAX
       INTEGER    NIN
       INTEGER NGASES
       REAL    PIN(MXIN)
       REAL  LNPIN(MXIN)
       REAL    TIN(MXIN)
       REAL   MRIN(MXIN,MXGAS)
       REAL    ZIN(MXIN)
       REAL    LAT
       REAL    LON
       REAL     PB(MYNLAY+1)
       REAL   LNPB(MYNLAY+1)
       INTEGER LAYBOT
       REAL  PSURF
       REAL  ZSURF
C
C      Output:
       INTEGER   NSUB(MYNLAY)
       INTEGER  NFINE
       REAL  PFINE(NSUBLV)
       REAL  TFINE(NSUBLV)
       REAL  MRFINE(NSUBLV,MXGAS)
       REAL   ZFINE(NSUBLV)
       INTEGER  WATID
       REAL    PSUB(NSUBLV-1)
       REAL    TSUB(NSUBLV-1)
       REAL   TGSUB(NSUBLV-1,MXGAS)
       REAL   MRSUB(NSUBLV-1,MXGAS)
       REAL   DZSUB(NSUBLV-1)
       REAL   ZPMAX


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      I
       INTEGER      J
       INTEGER      K
       REAL     AM(MXGAS)
       REAL     AT
       REAL     AZ
       REAL     BM(MXGAS)
       REAL     BT
       REAL     BZ
       REAL   CAVG(NSUBLV-1)
       REAL  CFINE(NSUBLV)
       REAL     D1
       REAL     D2
       REAL  GSURF
       REAL LNPFIN(NSUBLV)
       REAL  RJUNK
       REAL RJUNK2


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
c       print *,'profile ',IP,' IMAX = ',IMAX,' LAYBOT = ',LAYBOT
c       print *,'ZSURF = ',ZSURF,' PSURF = ',PSURF
c       print *,'wrt 101 layers : ',exp(LNPB(LAYBOT)) ,exp(LNPB(IMAX)) 
c       DO I = 1,LAYBOT
c         print *,'inside intlev',I,PIN(I),TIN(I),MRIN(I,1),MRIN(I,2)
c       END DO

C      ------------------------------------------------
C      Determine the sub-layer grid NSUB for each layer
C      ------------------------------------------------
C      Loop over the layers
       J=1
       DO I=LAYBOT,MYNLAY
C
C         Assign a default value to NSUB
          NSUB(I)=NGRID
C
C         Change NSUB if needed to ensure that there will be no more
C         than one PIN data point inside any AIRS sub-layer.
C
C         Find the first profile pressure level in PIN less than the
C         current AIRS pressure level PB(I). Assumes PIN is sorted
C         biggest to smallest pressure level.
 5        IF (PIN(J) .GT. PB(I) .AND. J .LT. NIN-1) THEN
             J=J+1
             GOTO 5
          ENDIF
C
C         Check to see if the current PIN point and the next one
C         both lie within the current AIRS layer
          IF (PIN(J) .LE. PB(I) .AND. PIN(J+1) .GE. PB(I+1)) THEN
C
C            Calculate PIN point spacing
             RJUNK=LNPIN(J)-LNPIN(J+1)
C
C            Calculate the minimum number of grids needed to at least
C            match the PIN spacing.
             K=1 + INT( (LNPB(I)-LNPB(I+1))/RJUNK )
C
             NSUB(I)=MAX(NSUB(I),K)
C
C            Check for any other pairs of PIN points in this layer
             IF (J .LT. NIN-1) THEN
                J=J+1
                GOTO 5
             ENDIF
C
          ENDIF
       ENDDO
C
C      ---------------------------------------------------------------
C      Calculate the fine grid points to define the sub-layers used in
C      calc'ing the layer weighted temperatures and integrated amounts
C      ---------------------------------------------------------------
C      Loop over the layers
       NFINE=1
       DO I=LAYBOT,MYNLAY
C
C         Bottom sub-level is same as layer lower boundary level
          PFINE(NFINE)=PB(I)
          LNPFIN(NFINE)=LNPB(I)
C
C         Calc grid spacing in terms of ln(P)
          RJUNK2=NSUB(I)
          RJUNK=(LNPB(I)-LNPB(I+1))/RJUNK2
C
C         Load up the other sub-levels
          DO K=1,NSUB(I)-1
             RJUNK2=K
             LNPFIN(NFINE+K)=LNPB(I) - (RJUNK2*RJUNK)
             PFINE(NFINE+K)=EXP( LNPFIN(NFINE+K) )
          ENDDO
C
C         Increment the counter of the total number of sub-levels
          NFINE=NFINE+NSUB(I)
C
       ENDDO
C
C      Topmost sub-level is the same as the top layer upper boundary
       PFINE(NFINE)=PB(MYNLAY+1)
       LNPFIN(NFINE)=LNPB(MYNLAY+1)
C
C      -----------------------------------------------------------------
C      Now that the sub-level pressure values are known, interpolate
C      the input profile temperature and mixing ratio (and Z, if
C      given) onto them in terms of ln(P).
C
C      To help prevent unexpected errors, the program checks to be
C      sure any interpolations below the lowest altitude pressure
C      (which will depend upon the two lowest altitude data points)
C      does not go too crazy. Here, we truncate any temperature
C      change to +-10K and any mixing ratio change to a factor of 2
C      and write a warning to the screen.
C      -----------------------------------------------------------------
C      Assign initital interpolation coefficient values
       J=2
       IF (LZ) THEN
          CALL INTERP( LNPIN(J),   ZIN(J), LNPIN(J-1),   ZIN(J-1),
     $       AZ, BZ )
       ENDIF
       CALL INTERP( LNPIN(J),  TIN(J), LNPIN(J-1),  TIN(J-1), AT,BT )
       DO K=1,NGASES
          CALL INTERP( LNPIN(J), MRIN(J,K), LNPIN(J-1), MRIN(J-1,K),
     $    AM(K), BM(K) )
       ENDDO
C
C      Loop over the output sub-levels 
       DO I=1,NFINE
C
C         Update the interpolation coefficients if necessary
 10       IF (PIN(J) .GT. PFINE(I) .AND. J .LT. NIN) THEN
             J=J+1
             CALL INTERP( LNPIN(J),  TIN(J), LNPIN(J-1),  TIN(J-1),
     $          AT,BT )
             DO K=1,NGASES
                CALL INTERP( LNPIN(J), MRIN(J,K), LNPIN(J-1),
     $             MRIN(J-1,K), AM(K), BM(K) )
             ENDDO
             IF (LZ) THEN
                CALL INTERP( LNPIN(J),   ZIN(J), LNPIN(J-1),
     $               ZIN(J-1), AZ, BZ )
             ENDIF
C
             GOTO 10
C
          ENDIF
C
C         Use the linear interpolation coefficients to calc T & MR
C         (and Z, if given, else calc the water MR)
C
          TFINE(I) =(LNPFIN(I)*AT) + BT
C
C         Check bottom-most temp interpolations for crazy behavior
          IF (PFINE(I) .GT. PIN(1)) THEN
             IF (TFINE(I) .LT. TIN(1)-10.0) THEN
                WRITE(IOINFO,1010) IP, I, PFINE(I), 'small', TFINE(I)
 1010           FORMAT('Warning: prof=',I5, ', sublay=',I4,
     $             ', pres=',F8.3,': truncating ',A5,' Temp=',F7.2)
                TFINE(I)=TIN(1)-10.0
             ENDIF
             IF (TFINE(I) .GT. TIN(1)+10.0) THEN
                WRITE(IOINFO,1010) IP, I, PFINE(I), 'large', TFINE(I)
                TFINE(I)=TIN(1)+10.0
             ENDIF
          ENDIF
C
          DO K=1,NGASES
             MRFINE(I,K)=(LNPFIN(I)*AM(K)) + BM(K)
C
C            Check bottom-most mixing ratio interps for crazy behavior
             IF (PFINE(I) .GT. PIN(1)) THEN
                IF (MRFINE(I,K) .LT. MRIN(1,K)/2.0) THEN
                   WRITE(IOINFO,1030) IP, I, K, PFINE(I), 'small',
     $                MRFINE(I,K)
 1030              FORMAT('Warning: prof=',I5,', sublay=',I4,
     $                ', gas=',I2,', pres=',F8.3,
     $                ': truncating ',A5,' mixing ratio=',1PE10.3)
                   MRFINE(I,K)=MRIN(1,K)/2.0
                ENDIF
                IF (MRFINE(I,K) .GT. MRIN(1,K)*2.0) THEN
                   WRITE(IOINFO,1030) IP, I, K, PFINE(I), 'large',
     $                MRFINE(I,K)
                   MRFINE(I,K)=MRIN(1,K)*2.0
                ENDIF
             ENDIF
C
C            Check for and correct water vapor over-saturation
C            in the lower atmosphere
             IF (LSVP .AND. (K .EQ. WATID) .AND.
     $       (PFINE(I) .GE. PB(IMAX))) THEN
C               Saturation mixing ratio in ppmv = (Psat/Pair)*10^6
                D1=WEXSVP( TFINE(I) )*1.0E+6/PFINE(I)
C
                IF (D1 .LT. MRFINE(I,K)) THEN
C                  Compute percent over-saturated
                   D2=((MRFINE(I,K)/D1) - 1.0)*100.0
C                  Write warning to screen
                   IF (D2 .GT. SATERR) THEN
                      WRITE(IOINFO,1050) IP, PFINE(I), D2
 1050                 FORMAT('! Warning: profile',I5,', ',F8.3,
     $                'mb, corrected',F8.2,'% over-saturation')
                   ENDIF
C                  Truncate mixing ratio to saturation value
                   MRFINE(I,K)=D1
                ENDIF
             ENDIF
          ENDDO
C
          IF (WATID .NE. 0) THEN
C            Note: 1.2028E-2 = (1E-3 kg/g) * (1E+2 mb/Pa) / R
C            with R = (8.314 J/Mole/K)
             CFINE(I)=1.2028E-2*( 28.97*(1.0 - (MRFINE(I,WATID)/
     $          1.0E+6)) + (18.01*MRFINE(I,WATID)/1.0E+6) )
          ELSE
C            Use dry air value = 1.2028E-2 * 28.97
             CFINE(I)=3.4845E-1
          ENDIF
C
          IF (I .GT. 1) THEN
             D2=CFINE(I)*PFINE(I)/TFINE(I)
             D1=CFINE(I-1)*PFINE(I-1)/TFINE(I-1)
             PSUB(I-1)=(PFINE(I)-PFINE(I-1))/LOG(PFINE(I)/PFINE(I-1))
             TSUB(I-1)=( (D1*TFINE(I-1)) + (D2*TFINE(I)) )/(D1+D2)
             DO K=1,NGASES
                TGSUB(I-1,K)=( (D1*TFINE(I-1)*MRFINE(I-1,K)) +
     $             (D2*TFINE(I)*MRFINE(I,K)) )/
     $             ( (D1*MRFINE(I-1,K)) + (D2*MRFINE(I,K)) ) 
                MRSUB(I-1,K)=0.5*(MRFINE(I,K)+MRFINE(I-1,K))
             ENDDO
          ENDIF
C
          IF (LZ) THEN
             ZFINE(I)=(LNPFIN(I)*AZ) + BZ
          ELSE
C            Calculate the layer heights if LZ is false
             IF (I .EQ. 1) THEN
C               Special calcs for bottom-most fine sublayer
C
C               Calc gravity at the surface
                GSURF=GRAV(ZSURF,0,0,LAT,LON)
C
C               Calc altitude at the bottom of the bottom-most
C               used layer
                IF (WATID .NE. 0) THEN
                   RJUNK=MRFINE(1,WATID)
                ELSE
                   RJUNK=0
                ENDIF
                ZFINE(1) = ZBOT( RJUNK, TFINE(1), GSURF, PFINE(1),
     $             PFINE(2), PSURF, ZSURF )
C
             ELSE
                CAVG(I-1)=( (D1*CFINE(I-1)) + (D2*CFINE(I)) )/(D1+D2)
                RJUNK=TSUB(I-1)*( PFINE(I)-PFINE(I-1) )/( CAVG(I-1)*
     +             PSUB(I-1) )
C
C               Convert RJUNK in K.mb/kg to K.Pa/kg
                RJUNK=RJUNK*100.0
C
C               Calculate gravity in m/s^2
C               Currently assumes no wind
                RJUNK2=-GRAV(ZFINE(I-1),0,0,LAT,LON)
C
                ZFINE(I)=ZFINE(I-1) + (RJUNK/RJUNK2)
C
             ENDIF
          ENDIF
C
          IF (I .GT. 1) DZSUB(I-1)=ZFINE(I)-ZFINE(I-1)
C
       ENDDO
C
C      Assign altitude at bottom of bottom-most used layer
       ZPMAX=ZFINE(1)

ccc
c      print *, 'zpmax=',ZPMAX
ccc

C
       RETURN
       END

C-----------------------------------------------------------------------
C
       SUBROUTINE INTERP(X1,Y1,X2,Y2,A,B)
C
C      Calculate linear interpolation coefficients A & B for Y = A*X + B
C      when given points (X1,Y1) and (X2,Y2).
C
C      Input: two data point pairs (x,y)
       REAL X1
       REAL Y1
       REAL X2
       REAL Y2
C
C      Output: constants for line y=a*x + b
       REAL A  ! slope
       REAL B  ! y-intercept
C
C      -----------------------------------------------------------------
C
       A=(Y2-Y1)/(X2-X1)
       B=Y1 - (A*X1)
C
       RETURN
       END
