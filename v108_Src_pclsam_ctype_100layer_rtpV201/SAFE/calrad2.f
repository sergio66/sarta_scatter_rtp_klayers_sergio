C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALRAD2
C
C    This version for 100 layer clouds using the PCLSAM method
C       PCLSAM = Parameterization for Cloud Longwave Scattering
C          for Atmospheric Models
C
!F77====================================================================


!ROUTINE NAME:
C    CALRAD2


!ABSTRACT:
C    Calculate the channel radiance for an atmosphere with two
C    complex clouds.


!CALL PROTOCOL:
C    CALRAD2( DOSUN, LCL, I, LBOT, RPLNCK, RSURFE, COSANG, SECANG,
C       SCOSL, SSECL, ODL, TRANL, TRANZ, SUNFAC, HSUN, TRANS,
C       RHOSUN, RHOTHR, COEFF,
C       CNGWA1, NEXTB1, NSCAB1, G_ASY1,
C       CNGWA2, NEXTB2, NSCAB2, G_ASY2, LTOPC ,LBOTC, RAD12 )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    LOGICAL   DOSUN   do sun radiance calcs?      true/false
C    LOG arr   LCL     complex cloud this layer?   none
C    INTEGER   I       channel index               none
C    INTEGER   LBOT    bottom layer                none
C    REAL arr  RPLNCK  Planck function             mW/(m^2 cm^-1 sterad)
C    REAL arr  RSURFE  surface emission            mW/(m^2 cm^-1 sterad)
C    REAL arr  COSANG  view angle cosine           none
C    REAL arr  SECANG  view angle secant           none
C    REAL arr  SCOSL   sun angle cosine            none
C    REAL arr  SSECL   sun angle secant            none
C    REAL arr  ODL     layer air optical depth     none
C    REAL arr  TRANL   layer air transmittances    none
C    REAL arr  TRANZ   surface-to-space air trans  none
C    REAL      SUNFAC  sun solid angle * cosine    sterad
C    REAL arr  HSUN    solar irradiance at TOA     mW/(m^2 cm^-1 sterad)?
C    REAL arr  TRANS   total solar air path trans  none
C    REAL      RHOSUN  solar surface reflectivity  none
C    REAL      RHOTHR  down thermal surf refl      none
C    REAL arr  COEFF   "F" factor coefficients     various
C    REAL arr  CNGWA1  cloud1 amount               g/m^2
C    REAL arr  NEXTB1  cloud1 extinction           per g/m^2
C    REAL arr  NSCAB1  cloud1 scattering           per g/m^2
C    REAL arr  G_ASY1  cloud1 asymmetry parameter  none
C    REAL arr  CNGWA2  cloud2 amount               g/m^2
C    REAL arr  NEXTB2  cloud2 extinction           per g/m^2
C    REAL arr  NSCAB2  cloud2 scattering           per g/m^2
C    REAL arr  G_ASY2  cloud2 asymmetry parameter  none
C    INTEGER   LTOPC   cloud top layer index       none
C    INTEGER   LBOTC   cloud bottom layer index    none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  RAD12   radiance                    mW/(m^2 cm^-1 sterad)


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    sarta


!ROUTINES CALLED:
C    function QIKEXP = EXP(X) calculation, faster than EXP(X) if X is small.
C    function HG2 = HG phase function


!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    Calculates the radiance for an atmosphere with two combined
C    clouds with the PCLSAM method.


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    The temperature is treated a constant within each layer (ie
C    no adjustments for temperature gradiants).


!ROUTINE HISTORY:
C    Date        Programmer     Comments
C    ----------- -------------- ----------------------------------------
C    22 Feb 2007 Scott Hannon   created (based on calrad1.f)
C    02 Mar 2007 Scott Hannon   Add LCL and revise "do" loops to use it.

!END====================================================================

C      =================================================================
       SUBROUTINE CALRAD2( DOSUN, LCL, I, LBOT, RPLNCK, RSURFE,
     $    COSANG, SECANG, SCOSL, SSECL, ODL, TRANL, TRANZ,
     $    SUNFAC, HSUN, TRANS, RHOSUN, RHOTHR, COEFF,
     $    CNGWA1, NEXTB1, NSCAB1, G_ASY1,
     $    CNGWA2, NEXTB2, NSCAB2, G_ASY2, LTOPC, LBOTC, RAD12 )
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
C      QIKEXP  : see file 'qikexp.f'
C      HG2     : see file 'hg2.f'

C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input
       LOGICAL  DOSUN      ! do sun radiance calcs?
       LOGICAL    LCL(MAXLAY) ! complex cloud this layer?
       INTEGER      I      ! channel index
       INTEGER   LBOT      ! bottom layer of atmosphere
       REAL RPLNCK(MAXLAY) ! layer Planck function
       REAL RSURFE         ! surface emission
       REAL COSANG(MAXLAY) ! viewing angle cosine
       REAL SECANG(MAXLAY) ! viewing angle secant
       REAL  SCOSL(MAXLAY) ! actual solar angle cosine
       REAL  SSECL(MAXLAY) ! actual solar angle secant
       REAL    ODL(MAXLAY,MXCHAN) ! clear air layer optical depth
       REAL  TRANL(MAXLAY) ! clear air layer transmittance
       REAL  TRANZ(MXCHAN) ! clear air surface-to-space transmittance
C      Sun info
       REAL SUNFAC         ! sun solid angle times cosine at surface
       REAL   HSUN(MXCHAN) ! irradiance from Sun at top of atmosphere
       REAL  TRANS(MXCHAN) ! up plus down clear air solar transmittance
       REAL RHOSUN         ! surface reflectivity for solar
C      Downwelling thermal info
       REAL RHOTHR         ! surface reflectivity for downwelling thermal
       REAL  COEFF(NFCOEF,MXCHAN) ! "F" factor coefficients
C      Cloud info
       REAL CNGWA1(MAXLAY) ! cloud1 amount (g/m^2)
       REAL NEXTB1(MXCHAN) ! cloud1 extinction (per g/m^2)
       REAL NSCAB1(MXCHAN) ! cloud1 scattering (per g/m^2)
       REAL G_ASY1(MXCHAN) ! cloud1 asymmetry
       REAL CNGWA2(MAXLAY) ! cloud2 amount (g/m^2)
       REAL NEXTB2(MXCHAN) ! cloud2 extinction (per g/m^2)
       REAL NSCAB2(MXCHAN) ! cloud2 scattering (per g/m^2)
       REAL G_ASY2(MXCHAN) ! cloud2 asymmetry
       INTEGER LTOPC       ! cloud top layer index
       INTEGER LBOTC       ! cloud bottom layer index
C
C      Output
       REAL  RAD12         ! upwelling radiance at satellite

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       LOGICAL DOSUNL(MAXLAY) ! layer solar scattering true/false
       INTEGER      L      ! layer index
       REAL      F         ! reflected therm "F" (fudge) factor
       REAL     GL(MAXLAY) ! layer scattering asymmetry
       REAL     K1         ! cloud optical depth
       REAL K1SUMS         ! total cloud optical depth along view angle
       REAL K1SUMV         ! total cloud optical depth along sun angle
       REAL   KAIR         ! layer nadir air (no cloud) optical depth
       REAL NEXTO1         ! nadir extinction optical depth cloud1
       REAL NSCAO1         ! nadir scattering optical depth cloud1
       REAL NEXTO2         ! nadir extinction optical depth cloud2
       REAL NSCAO2         ! nadir scattering optical depth cloud2
       REAL NEXTOD         ! nadir extinction optical depth both clouds
       REAL NSCAOD         ! nadir scattering optical depth both clouds
       REAL  ODSUM         ! sum of optical depth
       REAL ODTOTL(MAXLAY) ! total nadir layer optical depth
       REAL ODTOTZ(MAXLAY) ! tot nadir layer-above-to-space OD no scat adjust
       REAL  RADUP         ! upward radiance
       REAL   RSUN         ! reflected solar radiance
       REAL RSUNSC         ! scatter solar radiance
       REAL RTHERM         ! reflected downwelling thermal radiance
       REAL TRANLX(MAXLAY) ! layer transmittance
       REAL TRANCU         ! cloud transmittance at upward view angle
       REAL TRANCD         ! cloud trans at down sun angle no scat adjust
       REAL PI4INV         ! 1/4pi
C
       REAL XFUDGE(MAXLAY) ! Sergio's fudged optical depth for RSUNSC
       REAL WTILDE(MAXLAY) ! single scattering albedo including KAIR
C
C      Downwelling atmospheric thermal emission terms
       REAL TDOWNN ! "near-side" layer-to-surface trans
       REAL TDOWNF ! "far-side" layer-to-surface trans
       REAL  RDOWN ! downward radiance
C
C      for function QIKEXP
       REAL QIKEXP

C      for function HG2
       REAL HG2


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C***********************************************************************

       PI4INV = 1.0/(4.0*PI)
C
C      -----------------------------------------------------------------
C      Loop downward over the layers
C      -----------------------------------------------------------------
       ODSUM =0.0
       K1SUMV=0.0
       K1SUMS=0.0

C      Loop down over layers above cloud
       DO L=1,LBOT
C
          IF (LCL(L)) THEN
             DOSUNL(L)=DOSUN
             ODTOTZ(L)=ODSUM
             KAIR=ODL(L,I)/SECANG(L)
C
             NEXTO1=NEXTB1(I)*CNGWA1(L)
             NSCAO1=NSCAB1(I)*CNGWA1(L)
             NEXTO2=NEXTB2(I)*CNGWA2(L)
             NSCAO2=NSCAB2(I)*CNGWA2(L)
             NEXTOD=NEXTO1 + NEXTO2
             NSCAOD=NSCAO1 + NSCAO2
             GL(L)=( G_ASY1(I)*NSCAO1 + G_ASY2(I)*NSCAO2 )/NSCAOD
C
C            Cloud nadir optical depth including scattering adjustment
             K1=NEXTOD - NSCAOD*(1.0+GL(L))/2.0
             K1SUMV=K1SUMV + K1*SECANG(L)
             K1SUMS=K1SUMS + K1*SSECL(L)
C
             ODTOTL(L)=KAIR + K1
             ODSUM=ODSUM + KAIR + NEXTOD
             XFUDGE(L)=KAIR + NEXTOD
             WTILDE(L)=NSCAOD / XFUDGE(L)
C
             TRANLX(L)=QIKEXP( -ODTOTL(L)*SECANG(L) )
C
          ELSE
             DOSUNL(L)=.FALSE.
             ODTOTZ(L)=ODSUM         ! note ODTOTZ excludes current layer
             KAIR=ODL(L,I)/SECANG(L) ! clear air nadir optical depth
C
             ODTOTL(L)=KAIR
             ODSUM=ODSUM + KAIR
             TRANLX(L)=TRANL(L)
          ENDIF
       ENDDO
C
C      Calc the surface-to-space transmittance thru the cloud(only)
       TRANCU=QIKEXP( -K1SUMV ) ! upward path along view angle
C

C      -----------------------------------------------------------------
C      Loop upward over layers
C      -----------------------------------------------------------------
       RADUP=RSURFE
       RDOWN=0.0
       TDOWNN=1.0
       DO L=LBOT,1,-1

          IF (DOSUNL(L)) THEN
C            Scattered solar
             RSUNSC=(SCOSL(L)/(COSANG(L)+SCOSL(L)))*PI4INV*WTILDE(L)*
     $          HG2( -SCOSL(L), COSANG(L), GL(L) )*SUNFAC*HSUN(I)*
     $          QIKEXP( -ODTOTZ(L)*SSECL(L) )*
ccc fudged PCLSAM equation uses XFUDGE instead of ODTOTL
     $          (1.0 - QIKEXP( -XFUDGE(L)*(SECANG(L)+SSECL(L)) ))
ccc standard PCLSAM equation
c     $          (1.0 - QIKEXP( -ODTOTL(L)*(SECANG(L)+SSECL(L)) ))
C comment: According to Sergio Machado, the standard PCLSAM equation
C under-estimates the amount of solar radianced scattered into the
C view angle (RSUNSC).  If ODTOTL is replaced by XFUDGE the solar
C scattering term is increased.

          ELSE
             RSUNSC=0.0
          ENDIF
          RADUP=RADUP*TRANLX(L) + RPLNCK(L)*(1.0 - TRANLX(L)) + RSUNSC

C         Calc the downward radiance from this layer
          TDOWNF=TDOWNN*TRANLX(L)
          RDOWN = RDOWN + ( RPLNCK(L)*(TDOWNN - TDOWNF) )
          TDOWNN=TDOWNF

       ENDDO
C

C      ------------------------
C      Reflected solar radiance
C      ------------------------
C      Calc the reflected solar reaching the satellite
       IF (DOSUN) THEN
          TRANCD=QIKEXP( -K1SUMS ) ! downward path along sun angle
          RSUN=RHOSUN*SUNFAC*HSUN(I)*TRANS(I)*TRANCU*TRANCD
       ELSE
          RSUN=0.0
       ENDIF
C

C      --------------------------------------
C      Reflected downwelling thermal radiance
C      --------------------------------------
       F=1.0
       IF (TRANZ(I) .GT. 0.0005) THEN
          F=   COEFF(1,I) +
     $       ( COEFF(2,I)/SECANG(LBOT) ) +
     $       ( COEFF(3,I)*TRANZ(I) ) +
     $       ( COEFF(4,I)*TRANZ(I)*TRANZ(I) ) +
     $       ( COEFF(5,I)*TRANZ(I)/SECANG(LBOT) ) +
     $       ( COEFF(6,I)*TRANZ(I)/RDOWN )
C         Truncate F at limits as needed
          F = AMAX1( AMIN1(F,2.09), 0.696 )
       ENDIF
       RTHERM=RHOTHR*PI*RDOWN*F*TRANZ(I)*TRANCU
C

C      --------------
C      Total radiance
C      --------------
       RAD12=RADUP + RSUN + RTHERM
C
       RETURN
       END
