C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALRAD1
C    Calculate the channel radiance for an atmosphere with one cloud
C    This version of CALRAD uses "Parameterization for Cloud
C    Longwave Scattering for Atmospheric Models" (PCLSAM)
C
!F77====================================================================


!ROUTINE NAME:
C    CALRAD1


!ABSTRACT:
C    Calculate the channel radiance for an atmosphere with one cloud.


!CALL PROTOCOL:
C    CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
C       ODL, TAUL, TAUZ, SUNFAC, HSUN, TAUZSN, RHOSUN,
C       RHOTHR, LABOVE, COEFF,
C       CFRCL1, MASEC1, MASUN1, COSDAZ,
C       NEXTO1, NSCAO1, G_ASY1, LCTOP1, LCBOT1, RAD1, POLYNOM_BACKSCAT1, ISCALING1,
C       DOJAC, CLDTAU, RADLAY, RTHERM, CLDEFFOD )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    LOGICAL   DOSUN   do sun radiance calcs?      true/false
C    INTEGER   I       channel index               none
C    INTEGER   LBOT    bottom layer                none
C    REAL arr  RPLNCK  Planck function             mW/(m^2 cm^-1 sterad)
C    REAL arr  RSURFE  surface emission            mW/(m^2 cm^-1 sterad)
C    REAL arr  SECANG  path secant angles          none
C    REAL arr  TAU     layer transmittances        none
C    REAL arr  TAUZ    surface-to-space trans      none
C    REAL      SUNFAC  sun solid angle * cosine    sterad
C    REAL arr  HSUN    solar irradiance at TOA     mW/(m^2 cm^-1 sterad)?
C    REAL arr  TAUZ    surface-to-space trans      none
C    REAL arr  TAUZSN  surface-to-space trans      none
C    REAL arr  RHOSUN  solar surface reflectivity  none
C    REAL arr  RHOTHR  down thermal surf refl      none
C    INT arr   LABOVE  down therm layer above      none
C    REAL arr  COEFF   "F" factor coefficients     various
C    REAL arr  CFRCL1  fraction of cloud in layer  none
C    REAL      MASEC1  mean view secant in cloud   none
C    REAL      MASUN1  mean sun-only sec in cloud  none
C    REAL      COSDAZ  cosine(satazi-solazi)       none
C    REAL arr  NEXTO1  cloud extinction opt depth  none
C    REAL arr  NSCAO1  cloud scattering opt depth  none
C    REAL arr  G_ASY1  cloud asymmetry parameter   none
C    INTEGER   LCTOP1  cloud top layer index       none
C    INTEGER   LCBOT1  cloud bottom layer index    none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  RAD1    radiance                    mW/(m^2 cm^-1 sterad)
C    REAL arr  CLDTAU  ODs with clouds             none
C    REAL arr  RADLAY  rads at each layer          mW/(m^2 cm^-1 sterad)


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    sarta


!ROUTINES CALLED:
C    function QIKEXP = EXP(X) calculation, faster than EXP(X) if X is small.
C    function HG3 = HG phase function


!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    Calculates the radiance for an atmosphere with one cloud
C    with the PCLSAM method.


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    The temperature is treated a constant within each layer (ie
C    no adjustments for temperature gradiants).
C
C    The reflected downwelling emission term is adjusted for the
C    clouds after it is reflected from the surface, but the downward
C    radiance hitting the surface is for a clear sky.



!ROUTINE HISTORY:
C    Date        Programmer     Comments
C    ----------- -------------- ----------------------------------------
C    13 Jan 2006 Scott Hannon   Created
C    01 Feb 2006 S.Machado      Correct solar scattering by 2pi (replace
C       0.5 by PI4INV)
C    03 Feb 2006 S.Hannon       Remove scattering adjustments from some
C       downward optical depths to be consistenct with kcarta. Correct
C       ODSUM=0 initialization by moving it outside layer loop.
C    08 Feb 2006 S.Hannon       Correct error in RSUNSC so HSUN is
C       multiplied by SUNFAC not SCOSL
C    28 Mar 2006 S.Hannon       Implement Sergio's RSUNSC fudge
C       (non-standard PCLSAM) to improve agreement with other codes.
C       Change ODTOTZ used by RSUNSC to layer-above-to-space (was
C       (layer-to-space).
C    29 Mar 2006 S.Hannon       Updated RTHERM for sartaV107
C    03 Apr 2006 S.Hannon/S.Machado  Add missing w_tilde to RSUNSC
C    24 Mar 2008 S.Hannon       Add COSDAZ and use HG3 instead of HG2


!END====================================================================

C      =================================================================
       SUBROUTINE CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $    ODL, TAUL, TAUZ, SUNFAC, HSUN, TAUZSN, RHOSUN,
     $    RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, COSDAZ,
     $    NEXTO1, NSCAO1, G_ASY1, LCTOP1,LCBOT1, 
     $    POLYNOM_BACKSCAT1, ISCALING1,
     $    RAD1, DOJAC, CLDTAU, RADLAY, RTHERM, CLDEFFOD )
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
C      HG3     : see file 'hg3.f'

C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input
       LOGICAL  DOSUN      ! do sun radiance calcs?
       INTEGER      I      ! channel index
       INTEGER   LBOT      ! bottom layer of atmosphere
       REAL RPLNCK(MAXLAY) ! layer Planck function
       REAL RSURFE         ! surface emission
       REAL SECANG(MAXLAY) ! viewing angle secant
       REAL    ODL(MAXLAY,MXCHAN) ! clear air layer optical depth
       REAL   TAUL(MAXLAY) ! clear air layer transmittance
       REAL   TAUZ(MXCHAN) ! clear air surface-to-space transmittance
C      Sun info
       REAL SUNFAC         ! sun solid angle times cosine at surface
       REAL   HSUN(MXCHAN) ! irradiance from Sun at top of atmosphere
       REAL TAUZSN(MXCHAN) ! up plus down clear air solar transmittance
       REAL RHOSUN(MXCHAN) ! surface reflectivity for solar
C      Downwelling thermal info
       REAL RHOTHR(MXCHAN) ! surface reflectivity for downwelling thermal
       INTEGER LABOVE(MXCHAN) ! representative layer above surface
       REAL  COEFF(NFCOEF,MXCHAN) ! "F" factor coefficients
C      Cloud1 info
       REAL CFRCL1(MAXLAY) ! fraction of cloud in layer
       REAL MASEC1         ! mean view secant in cloud
       REAL MASUN1         ! mean sun-only secant in cloud
       REAL COSDAZ         ! cosine of delta azimuth angles
       REAL POLYNOM_BACKSCAT1(4) ! polynom backscat coeffs
       INTEGER ISCALING1         !type of scho.similarity,maestri scaling
       REAL NEXTO1(MXCHAN) ! cloud nadir extinction optical depth
       REAL NSCAO1(MXCHAN) ! cloud nadir scattering optical depth
       REAL G_ASY1(MXCHAN) ! cloud asymmetry
       INTEGER LCTOP1      ! cloud top layer index
       INTEGER LCBOT1      ! cloud bottom layer index
       LOGICAL DOJAC
C
C      Output
       REAL   RAD1         ! upwelling radiance at satellite
       REAL CLDTAU(MAXLAY)        ! chan layer effective optical depth for CLDONLY
       REAL RADLAY(MAXLAY)        ! chan layer radiance                for CLDONLY
       REAL RTHERM         ! reflected downwelling thermal radiance, used to be local param
       REAL CLDEFFOD       ! total cloud effect OD (from PCLSAM)
C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       LOGICAL DOSUNL(MAXLAY) ! layer solar scattering true/false
       INTEGER      L      ! layer index
       INTEGER LTHERM      ! layer for RTHERM calc
       REAL      F         ! reflected therm "F" (fudge) factor
       REAL     GL(MAXLAY) ! layer scattering asymmetry
       REAL     K1         ! cloud1 optical depth
       REAL   KAIR         ! layer nadir air (no cloud) optical depth
       REAL  ODSUM         ! sum of optical depth
       REAL ODTOTL(MAXLAY) ! total nadir layer optical depth
       REAL ODTOTZ(MAXLAY) ! tot nadir layer-above-to-space OD no scat adjust
       REAL RADUP          ! upward radiance
       REAL RSUN           ! reflected solar radiance
       REAL RSUNSC         ! scatter solar radiance
       REAL  SSECL(MAXLAY) ! solar angle secant
       REAL  SCOSL(MAXLAY) ! solar angle cosine
       REAL  TAULX(MAXLAY) ! layer transmittance
       REAL TAUZCU         ! cloud transmittance at upward view angle
       REAL TAUZCD         ! cloud trans at down sun angle no scat adjust
       REAL  VCOSL(MAXLAY) ! view angle cosine
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

C      for function HG3
       REAL HG3

c      for Tang correction to Chou calcs, rXTang is in incFTC.f
       REAL raTangCorrect(MAXLAY)
       REAL raTauL2S(MAXLAY),Bscatter

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C***********************************************************************

       IF (DOJAC) THEN 
         CLDTAU = 0.0
         RADLAY = 0.0
       END IF

       PI4INV = 1.0/(4.0*PI)
C
C      Optical depth of cloud1 including scattering adjustment
C      this is similarity scaling, see eg Chou 1999 or Martinazzo/Maestri JQSRT 2021
C      abs(floor([0 0.2 0.5 0.7 -1 -2 -3])) = 0     0     0     0     1     2     3
       IF (ISCALING1 .LE. 1) THEN
         !! similarity scaling
         K1 = NEXTO1(I) - NSCAO1(I)*(1.0+G_ASY1(I))/2.0   
       ELSE
         !! this is Chou scaling using eqn 7 in Martinazzo/Maestri JQSRT 2021
         K1 = 1 - (POLYNOM_BACKSCAT1(1) + POLYNOM_BACKSCAT1(2)*G_ASY1(I) + 
     $             POLYNOM_BACKSCAT1(3)*G_ASY1(I)**2 + POLYNOM_BACKSCAT1(4)*G_ASY1(I)**3)
         K1 = NEXTO1(I)-NSCAO1(I) + K1*NSCAO1(I)
       END IF

c      print *,'CALRAD1 -- ISCALING1',ISCALING1
c      print *,'CALRAD1 -- POLYNOM_BACKSCAT1 = ',POLYNOM_BACKSCAT1

       CLDEFFOD = K1

c      prep for Tang correction
       IF (rXTang .GT. 0.0) THEN
         raTangCorrect = 0.0 
         RDOWN = 0.0
         TDOWNN=1.0
       END IF

C      -----------------------------------------------------------------
C      Loop downward over the layers
C      -----------------------------------------------------------------
       ODSUM=0.0
       DO L=1,LBOT
          DOSUNL(L)=.FALSE.
          KAIR=ODL(L,I)/SECANG(L)
c added 28 Mar 2006; layer-above-to-space
          ODTOTZ(L)=ODSUM
          WTILDE(L)=0.0
C
          IF (CFRCL1(L) .GT. 0.0) THEN
             SSECL(L)=MASUN1       ! note: if no sun, this is garbage
             SCOSL(L)=1.0/SSECL(L) ! note: if no sun, this is garbage
             VCOSL(L)=1.0/SECANG(L)
             GL(L)=G_ASY1(I)
             DOSUNL(L)=DOSUN
             ODTOTL(L)=KAIR + K1*CFRCL1(L)
C replaced 03Feb2006          ODSUM=ODSUM + ODTOTL(L)
c these are used in SW when sun = on (if DOSUN(L)), else not really used in LW
c so derivatives wrt CFRCL1(L) should not be impacted
             ODSUM=ODSUM + KAIR + NEXTO1(I)*CFRCL1(L)
             XFUDGE(L)=KAIR + NEXTO1(I)*CFRCL1(L)
             WTILDE(L)=CFRCL1(L)*NSCAO1(I) / XFUDGE(L)
          ELSE
             XFUDGE(L)=KAIR
             ODTOTL(L)=KAIR
             ODSUM=ODSUM + KAIR
          ENDIF

c Tang
          IF ((rXTang .GT. 0.0) .AND. (CFRCL1(L) .GT. 0.0)) THEN
            raTauL2S(L) = ODSUM
            TDOWNF=TDOWNN*exp(-ODSUM)
            RDOWN = RDOWN + ( RPLNCK(L)*(TDOWNN - TDOWNF) )
            TDOWNN=TDOWNF
 
c           Tang uses factor of 0.5; see /home/sergio/KCARTA/SRCv1.22_f90/Revisions_cloud.txt where I say 0.10 for water clouds, 0.30 for ice clouds
            Bscatter    = 1-GL(L)
            raTangCorrect(L) = (RDOWN-RPLNCK(L))*(1-ODSUM*ODSUM)*TDOWNN
            raTangCorrect(L) = 0.20*WTILDE(L)*Bscatter/(1-WTILDE(L)*(1-Bscatter))*raTangCorrect(L)
         END IF

         IF (DOJAC) CLDTAU(L) = ODTOTL(L)*SECANG(L)
C
c removed 28 Mar 2006; layer-to-space
c          ODTOTZ(L)=ODSUM
C
          IF (L .LT. LCTOP1 .OR. L .GT. LCBOT1) THEN
             TAULX(L)=TAUL(L)
          ELSE
             TAULX(L)=QIKEXP( -ODTOTL(L)*SECANG(L) )
          ENDIF

       ENDDO ! downward loop over layers

C      Calc the surface-to-space transmittance thru the clouds(only)
       TAUZCU=QIKEXP(-K1*MASEC1) ! upward path
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
             RSUNSC=(SCOSL(L)/(VCOSL(L)+SCOSL(L)))*PI4INV*WTILDE(L)*
     $          HG3(-SCOSL(L),VCOSL(L),COSDAZ,GL(L))*SUNFAC*HSUN(I)*
     $          QIKEXP( -ODTOTZ(L)*SSECL(L) )*
ccc fudged PCLSAM equation uses XFUDGE instead of ODTOTL
     $          (1.0 - QIKEXP( -XFUDGE(L)*(SECANG(L)+SSECL(L)) ))
ccc standard PCLSAM equation
c     $          (1.0 - QIKEXP( -ODTOTL(L)*(SECANG(L)+SSECL(L)) ))
C comment: According to Sergio Machado, the standard PCLSAM equation
C under-estimates the amount of solar radianced scattered into the
C view angle (RSUNSC).  If ODTOTL is replaced by XFUDGE the solar
C scattering term is increased.

c check for infinity
            IF ( (ABS(RSUNSC) .GT. 1000) .OR. ( (ABS(RSUNSC)-1) .EQ. ABS(RSUNSC)))  RSUNSC = 0.0
          ELSE
             RSUNSC=0.0
          ENDIF
          RADUP=RADUP*TAULX(L) + RPLNCK(L)*(1.0 - TAULX(L)) + RSUNSC
          IF (DOJAC) RADLAY(L) = RADUP

C         Calc the downward radiance from this layer
          TDOWNF=TDOWNN*TAULX(L)
          RDOWN = RDOWN + ( RPLNCK(L)*(TDOWNN - TDOWNF) )
          TDOWNN=TDOWNF

       ENDDO
C

C      ------------------------
C      Reflected solar radiance
C      ------------------------
C      Calc the reflected solar reaching the satellite
       IF (DOSUN) THEN
C replaced 03Feb2006          TAUZCD=QIKEXP( -MASUN1*K1 ) ! downward path
          TAUZCD=QIKEXP( -MASUN1*NEXTO1(I) ) ! downward path
          RSUN=RHOSUN(I)*SUNFAC*HSUN(I)*TAUZSN(I)*TAUZCU*TAUZCD
c check for infinity
          IF ( (ABS(RSUN) .GT. 1000) .OR. ( (ABS(RSUN)-1) .EQ. ABS(RSUN)))  RSUN = 0.0
       ELSE
          RSUN=0.0
       ENDIF
C

C      --------------------------------------
C      Reflected downwelling thermal radiance
C      --------------------------------------
       F=1.0
       IF (TAUZ(I) .GT. 0.0005) THEN
          F=   COEFF(1,I) +
     $       ( COEFF(2,I)/SECANG(LBOT) ) +
     $       ( COEFF(3,I)*TAUZ(I) ) +
     $       ( COEFF(4,I)*TAUZ(I)*TAUZ(I) ) +
     $       ( COEFF(5,I)*TAUZ(I)/SECANG(LBOT) ) +
     $       ( COEFF(6,I)*TAUZ(I)/RDOWN )
C         Truncate F at limits as needed
          F = MAX( MIN(F,2.09), 0.696 )
       ENDIF
       RTHERM=RHOTHR(I)*PI*RDOWN*F*TAUZ(I)*TAUZCU
C

C      --------------
C      Total radiance
C      --------------

c      raTauL2S(L) = ODSUM
       RAD1=RADUP + RSUN + RTHERM
c       IF (I .EQ. 1100) print *,'calrad1 final ',RADUP,RSUN,RTHERM,RAD1,sum(raTangCorrect)*rXTang

       IF (rXTang .GT. 0.0) THEN
         RAD1 = RAD1 + sum(raTangCorrect)*rXTang
       ENDIF

       IF ((DOJAC) .AND. (DOSUN)) RTHERM = RSUN + RTHERM

       RETURN
       END
