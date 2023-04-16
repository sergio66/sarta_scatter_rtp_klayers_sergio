C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALRAD
C
!F77====================================================================


!ROUTINE NAME:
C    CALRAD


!ABSTRACT:
C    Calculate a profile's radiance.


!CALL PROTOCOL:
C    CALRAD ( NCHAN, FREQ, TAU, TP, TBOT, EBOT, LBOT,
C             SUNCOS, RHOSUN, DISTES, HSUN, TAUZSN,
C             SECANG, RHOTHR, LABOVE, COEFF, COEFNU, TAUZ, RAD, BT)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  COEFF   thermal F factor coefs      various
C    REAL arr  COEFNU  channl freq coefs           various
C    REAL      DISTES  Earth-sun distance          meters
C    REAL arr  EBOT    bottom surface emissivity   none
C    REAL arr  FREQ    channel frequencies         cm^-1
C    INT arr   LABOVE  layer-above for thermal     none
C    INTEGER   LBOT    bottom layer                none
C    INTEGER   NCHAN   number of channels          none
C    REAL arr  RHOSUN  reflectivity for solar      1/steradian
C    REAL arr  RHOTHR  reflectivity for thermal    1/steradian
C    REAL arr  SECANG  path angle secant           none
C    REAL      SUNCOS  sun angle cosine            none
C    REAL arr  TAU     effective layer trans       none
C    REAL arr  TAUZ    layer-to-space trans        none
C    REAL arr  TAUZSN  eff sun angle l-to-s trans  none
C    REAL      TBOT    bottom surface temperature  Kelvin
C    REAL arr  TP      temperature profile         Kelvin


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  BT      brightness temperature      Kelvin
C    REAL arr  RAD     radiance                    W/(m^2.str.cm^-1)


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    USEFAST


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    March 1998 version of the 100 layer AIRS Fast Transmittance
C    Code by L.L.Strow/S.Hannon.
C
C    The radiance is calculated for each channels in turn.  The rad
C    is a sum of four components: surface, upwelling (non-reflected)
C    atmospheric, reflected downwelling atmospheric thermal, and
C    reflected solar.  No scattering.
C
C    Comment: this routine could easily be re-written to use layer-to-
C    space transmittances rather than layer transittances.  Currently
C    the CALT# routines compute layer transmittances (since it's a bit
C    faster and more accurate when using the QIKEXP function.
C
C    ===================================================================
C    Computes black body emissions for each layer using the Planck
C    equation:
C       planck = c1*v^3/( exp(c2*v/T) - 1 )
C    where c1 and c2 are the radiation constants, T is the temperature
C    TP, and v is the frequency FREQ.
C
C    We assume the layers emit radiances of
C       rad_layer = (1 - tau)*planck
C    where tau is the layer transmittance TAU.
C
C    The total radiance leaving the bottom surface and going upward
C    is the surface emission and reflected solar & thermal. The
C    reflected solar and thermal are handled as seperate terms added
C    to radiance arriving at the satellite.
C       rad_surface =  e*planck
C    where e is the bottom surface emissivity EBOT, and the surface is
C    at temperature TBOT.
C
C    We trace the upward radiance thru the atmosphere and determine
C    the total radiance leaving the top layer (and then reaching the
C    satellite) is:
C       the sum L=L_bot downto 1 of { rad(L-1)*tau(L) + rad(L) }
C    where rad(L_bot-1) = rad_surface, and rad(1) = RAD.
C
C    The reflected solar term is based on an approximation suggested
C    by J.Susskind et al.  The reflected solar radiance reaching the
C    satellite is given by
C       Rsun = rho * omega * TAUZSN * Hsun
C    where omega is the solid angle of the sun as seen from Earth,
C    Hsun is the (non-reflected) solar radiance at the top of the
C    atmosphere, and TAUZSN is (surface) layer-to-space transmittance
C    of a path along an effective total angle defined as
C       secant_eff = secant_view + secant_sun
C    Note that this requires a seperate transmittance calculation
C    at the effective sun angle.  Hsun is passed to this routine (it
C    is close to planck for 5800 K), while omega is computed using
C    the distance of the Earth from the sun DISTES
C       omega = pi * ( radius_sun / distance_Earth_sun )^2
C
C    The reflected downwelling thermal is another approximation based
C    on a method suggested by Susskind et al.  It uses the viewing
C    angle layer-to-space transmittance TAUZ, the radiance of a
C    single layer somewhere above the surface, and a parameterized "F"
C    factor (which is sort of a fudge factor determined by regression).
C       Rtherm = rho * pi * F * planck * TAUZ*(1-TAUZ)
C    where the planck radiance is computed for layer L = LBOT - LABOVE.
C
C    For convenience we also output brightness temperature, which is
C    related to the radiance by inverting the planck equation:
C       BT = c2*v/ln( 1 + c1*v^3/RAD )
C    ===================================================================


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    Currently this routine does not handle scattering or clouds.


!ROUTINE HISTORY:
C    Date        Programmer     Comments
C    ----------- -------------- ----------------------------------------
C     2 Sep 1997 Scott Hannon   Created from an extensive re-write of
C                               our Feb97 CALRAD routine for Mar98 FTC.
C    12 Aug 2005 Scott Hannon   Added COEFNU and change scaler SEC to
C                               vector SECANG; change loop over layer
C                               to downward.
C    26 Aug 2005 Scott Hannon   Revised fshift predictors with fmin &
C                               fmax checks


!END====================================================================

C      =================================================================
       SUBROUTINE CALRAD ( NCHAN, FREQ, TAU, TP, TBOT, EBOT, LBOT,
     $    SUNCOS, RHOSUN, DISTES, HSUN, TAUZSN,
     $    SECANG, RHOTHR, LABOVE, COEFF, COEFNU, TAUZ, RAD, BT)
C      =================================================================

C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
       include 'incFTC.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
       INTEGER NCHAN
       REAL   FREQ(MXCHAN)
       REAL    TAU(MAXLAY,MXCHAN)
       REAL     TP(MAXLAY)
       REAL   TBOT
       REAL   EBOT(MXCHAN)
       INTEGER LBOT
       REAL SUNCOS
       REAL RHOSUN(MXCHAN)
       REAL DISTES
       REAL   HSUN(MXCHAN)
       REAL TAUZSN(MXCHAN)
       REAL SECANG(MAXLAY)
       REAL RHOTHR(MXCHAN)
       INTEGER LABOVE(MXCHAN)
       REAL  COEFF(NFCOEF,MXCHAN)
       REAL COEFNU(NCOFNU,MAXLAY+1,MXCHAN)
       REAL   TAUZ(MXCHAN)
       REAL    RAD(MXCHAN)
       REAL     BT(MXCHAN)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      I
       INTEGER      J
       INTEGER      L
       INTEGER     LS
       INTEGER LTHERM
       REAL   C1V3
       REAL    C2V
       REAL      F
       REAL FSHIFT
       REAL FPLNCK
       REAL NUPRED(8)
       REAL  RJUNK
       REAL RPLNCK(MAXLAY)
       REAL   RSUN
       REAL RTHERM
       REAL SUNFAC
       REAL     T0
       REAL     TZ
       REAL   TZM1


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
C      Assign NU shift predictor elements
       T0 = 273.15
       NUPRED(1) = 1.0

C      Note: on average, DISTES = 1.496E+11 m.  The exact value varies
C      with time since the Earth's orbit is slightly elliptical.
       SUNFAC=SUNCOS*PI*(RADSUN/DISTES)**2
C      Note: PI*(RADSUN/DISTES)^2 = omega = solid angle [steradians] of
C      the sun as seen from Earth.  The above equation is actually an
C      approximation for the case DISTES >> RADSUN.  The exact equation
C      is omega = 2*pi*(1 - DISTES/sqrt(DISTES^2 + RADSUN^2)).
C
C      ---------------------------
C      Loop on channel (frequency)
C      ---------------------------
       DO 210 I=1,NCHAN
C
C         ---------------------------
C         Calc the up going radiance
C         (excluding reflected terms)
C         ---------------------------
C
C         Loop downward over the layers
          RAD(I) = 0.0
          TZ = 1.0
          DO L=1,LBOT
C
C            Calculate the layer-to-space trans
             TZM1 = TZ
             TZ = TAU(L,I)*TZM1
C
C            Assign NU shift predictors
             NUPRED(2) = TP(L)/T0
             NUPRED(3) = SECANG(L)
             NUPRED(4) = 1.0/SECANG(L)
             NUPRED(5) = T0/TP(L)
             NUPRED(6) = SECANG(L)*TP(L)/T0
             NUPRED(7) = TZ
             IF (COEFNU(8,L,I) .EQ. 0.0) THEN
                NUPRED(8) = 0.0
             ELSE
                IF (TZ .GE. 1.0) THEN
C                  Error trap for TZ=1
                   NUPRED(8) = 0.0
                ELSE
                   IF (TZ .LT. 4.5E-5) THEN
C                     Error trap for TZ=0
                      NUPRED(8) = 10.0
                   ELSE
                      NUPRED(8) = -LOG(TZ)
                   ENDIF
                ENDIF
             ENDIF

C            Calculate the effective planck freq
             FSHIFT =
     $       NUPRED(1)*COEFNU(1,L,I) + NUPRED(2)*COEFNU(2,L,I) +
     $       NUPRED(3)*COEFNU(3,L,I) + NUPRED(4)*COEFNU(4,L,I) +
     $       NUPRED(5)*COEFNU(5,L,I) + NUPRED(6)*COEFNU(6,L,I) +
     $       NUPRED(7)*COEFNU(7,L,I) + NUPRED(8)*COEFNU(8,L,I)
C
             FPLNCK = FREQ(I) + FSHIFT
ccc for testing
c      FPLNCK=FREQ(I)
ccc
C
C            Enforce range limits on FPLNCK
             IF (FPLNCK .LT. COEFNU( 9,L,I)) FPLNCK = COEFNU( 9,L,I)
             IF (FPLNCK .GT. COEFNU(10,L,I)) FPLNCK = COEFNU(10,L,I)
C
C            Calc c1*v^3 and c2*v
             C1V3=C1*(FPLNCK**3)
             C2V=C2*FPLNCK
C
C            Calculate the Planck function for this layer
             RPLNCK(L)=C1V3/( EXP( C2V/TP(L) ) - 1.0 )
C
C            Calc the radiance contribution from this layer
             RAD(I)=RAD(I) + (TZM1 - TZ)*RPLNCK(L)

          ENDDO
C
C         Add on the bottom surface emission
          NUPRED(2) = TBOT/T0
          NUPRED(3) = SECANG(LBOT)
          NUPRED(4) = 1.0/SECANG(LBOT)
          NUPRED(5) = T0/TBOT
          NUPRED(6) = SECANG(LBOT)*TBOT/T0
          NUPRED(7) = TAUZ(I)
          IF (TAUZ(I) .LT. 4.5E-5) THEN
C            Error trap for TZ=0
             NUPRED(8) = 10.0
          ELSE
             NUPRED(8) = -LOG(TAUZ(I))
          ENDIF
          LS = MAXLAY + 1  ! layer index for fshift surface coeffs
          FSHIFT =
     $       NUPRED(1)*COEFNU(1,LS,I) + NUPRED(2)*COEFNU(2,LS,I) +
     $       NUPRED(3)*COEFNU(3,LS,I) + NUPRED(4)*COEFNU(4,LS,I) +
     $       NUPRED(5)*COEFNU(5,LS,I) + NUPRED(6)*COEFNU(6,LS,I) +
     $       NUPRED(7)*COEFNU(7,LS,I) + NUPRED(8)*COEFNU(8,LS,I)
          FPLNCK = FREQ(I) + FSHIFT
ccc for testing
c      FPLNCK=FREQ(I)
ccc
          IF (FPLNCK .LT. COEFNU( 9,LS,I)) FPLNCK = COEFNU( 9,LS,I)
          IF (FPLNCK .GT. COEFNU(10,LS,I)) FPLNCK = COEFNU(10,LS,I)
          C1V3=C1*(FPLNCK**3)
          C2V=C2*FPLNCK
          RAD(I)=RAD(I) + TAUZ(I)*EBOT(I)*C1V3/( EXP( C2V/TBOT ) - 1.0 )

C         --------------------------
C         Calc the reflected solar
C         rad reaching the satellite
C         --------------------------
          RSUN=RHOSUN(I)*SUNFAC*HSUN(I)*TAUZSN(I)
C
C         ----------------------------------
C         Calc the reflected downwelling
C         thermal rad reaching the satellite
C         ----------------------------------
          IF (LABOVE(I) .GT. 0 .AND. TAUZ(I) .GT. 0.0001) THEN
             LTHERM=LBOT - LABOVE(I)
C            Note: fractions of a layer may be ignored for LTHERM
             F=   COEFF(1,I) +
     $          ( COEFF(2,I)/SECANG(LBOT) ) +
     $          ( COEFF(3,I)*RPLNCK(LTHERM) ) +
     $          ( COEFF(4,I)*RPLNCK(LTHERM)/SECANG(LBOT) ) +
     $          ( COEFF(5,I)*RPLNCK(LBOT)/RPLNCK(LTHERM) )
             RTHERM=RHOTHR(I)*PI*RPLNCK(LTHERM)*F*TAUZ(I)*
     $          (1.0E+0 - TAUZ(I))
          ELSE
             RTHERM=0.0E+0
          ENDIF
C
C         --------------------------------------------------
C         Add on the reflected solar and downwelling thermal
C         --------------------------------------------------
ccc
c for testing
c      RTHERM=0.0
c      RSUN=0.0
ccc
          RAD(I)=RAD(I) + RSUN + RTHERM
C
C         ------------------------------------------
C         Convert radiance to brightness temperature
C         ------------------------------------------
          BT(I)=C2V/LOG( 1.0 + C1V3/RAD(I) )
C
 210   CONTINUE
C      End loops on channel number (frequency)
C
       RETURN
       END
