C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALRAD0
C    Calculate channel radiance for a clear atmosphere above a surface
C
!F77====================================================================


!ROUTINE NAME:
C    CALRAD0


!ABSTRACT:
C    Calculate channel radiance for a clear atmosphere above a surface.


!CALL PROTOCOL:
C    CALRAD0(DOSUN, I, LBOT, RPLNCK, TSURF, RSURFE, SECANG,
C       TAUL, TAUZ, SUNFAC, HSUN, TAUZSN, RHOSUN,
C       RHOTHR, LABOVE, COEFF, RAD0, DOJAC, CLDTAU, RADLAY, RTHERM )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    LOGICAL   DOSUN   do sun radiance calcs?      true/false
C    INTEGER   I       channel index               none
C    INTEGER   LBOT    bottom layer                none
C    REAL arr  RPLNCK  Planck function             mW/(m^2 cm^-1 sterad)
C    REAL      TSURF   surface temperature         K, to help debug N. Nalli emiss
C    REAL      RSURFE  surface emission            mW/(m^2 cm^-1 sterad)
C    REAL arr  SECANG  path secant angles          none
C    REAL arr  TAUL    layer transmittance         none
C    REAL arr  TAUZ    surface-to-space trans      none
C    REAL      SUNFAC  sun solid angle * cosine    sterad
C    REAL arr  HSUN    solar irradiance at TOA     mW/(m^2 cm^-1 sterad)?
C    REAL arr  TAUZ    surface-to-space trans      none
C    REAL arr  TAUZSN  surface-to-space trans      none
C    REAL arr  RHOSUN  solar surface reflectivity  none
C    REAL arr  RHOTHR  down thermal surf refl      none
C    INT arr   LABOVE  down therm layer above      none
C    REAL arr  COEFF   "F" factor coefficients     various


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  RAD0    radiance                    mW/(m^2 cm^-1 sterad)
C    REAL arr  CLDTAU  ODs with clouds             none
C    REAL arr  RADLAY  rads at each layer          mW/(m^2 cm^-1 sterad)


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    sarta


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    Calculates the radiance for an atmosphere with no clouds.


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    The temperature is treated a constant within each layer (ie
C    no adjustments for temperature gradiants).



!ROUTINE HISTORY:
C    Date        Programmer     Comments
C    ----------- -------------- ----------------------------------------
C    13 Jan 2006 Scott Hannon   Created
C    29 Mar 2006 Scott Hannon   Updated RTHERM for sartaV107


!END====================================================================

C      =================================================================
       SUBROUTINE CALRAD0( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $    TAUL, TAUZ, SUNFAC, HSUN, TAUZSN, RHOSUN,
     $    RHOTHR, LABOVE, COEFF, RAD0, DOJAC, CLDTAU, RADLAY, RTHERM,
     $    EMIS, TSURF, IPROF, FREQ )
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
C      Input
       LOGICAL  DOSUN      ! do sun radiance calcs?
       INTEGER      I      ! channel index
       INTEGER   LBOT      ! bottom layer of atmosphere
       REAL RPLNCK(MAXLAY) ! layer Planck function
       REAL RSURFE         ! surface emission
       REAL SECANG(MAXLAY) ! viewing angle secant
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
       LOGICAL DOJAC

c for Nick Nalli
       REAL EMIS(MXCHAN)   ! surface emissivity
       REAL TSURF          ! surface temperature
       INTEGER IPROF
       REAL FREQ(MXCHAN)   ! channel center

C
C      Output
       REAL   RAD0                ! upwelling radiance at satellite
       REAL CLDTAU(MAXLAY)        ! chan layer effective optical depth for CLDONLY
       REAL RADLAY(MAXLAY)        ! chan layer radiance                for CLDONLY
       REAL RTHERM                ! reflected downwelling thermal radiance, now an output

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      L      ! layer index
       INTEGER LTHERM      ! layer for RTHERM calc
       REAL      F         ! reflected therm "F" (fudge) factor
       REAL RADUP          ! upward radiance
       REAL RSUN           ! reflected solar radiance
       INTEGER iNalli      ! usual or nalli refl thermal

C      Downwelling atmospheric thermal emission terms
       REAL TDOWNN ! "near-side" layer-to-surface trans
       REAL TDOWNF ! "far-side" layer-to-surface trans
       REAL  RDOWN ! downward radiance

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
         !CLDTAU = -log(TAUL)
       END IF


C      -----------------------------------------------------------------
C      Loop upward over layers
C      -----------------------------------------------------------------
       RADUP=RSURFE
c       write(*,'(A,I,3F12.5)'),'sseerrggiioo',I,FREQ(I),TSURF,RADUP
       RDOWN=0.0
       TDOWNN=1.0
       DO L=LBOT,1,-1
          RADUP=RADUP*TAUL(L) + RPLNCK(L)*(1.0 - TAUL(L))
          IF (DOJAC) RADLAY(L) = RADUP

C         Calc the downward radiance from this layer
          TDOWNF=TDOWNN*TAUL(L)
          RDOWN = RDOWN + ( RPLNCK(L)*(TDOWNN - TDOWNF) )
          TDOWNN=TDOWNF
c          write(*,'(A,I,4F12.5)') 'sseerrggiioo',L,TAUL(L),RPLNCK(L),RADUP,RDOWN
       ENDDO
C
c       IF (I .EQ. 1291) THEN
c         PRINT *,'FINAL',RSURFE,L,RPLNCK(L),RADUP
c         PRINT *,'sergio stop A'
c       END IF

C      ------------------------
C      Reflected solar radiance
C      ------------------------
       IF (DOSUN) THEN
          RSUN=RHOSUN(I)*SUNFAC*HSUN(I)*TAUZSN(I)
c check for infinity
          IF ( (ABS(RSUN) .GT. 1000) .OR. ( (ABS(RSUN)-1) .EQ. ABS(RSUN)))  RSUN = 0.0
       ELSE
          RSUN=0.0
       ENDIF
C

C      --------------------------------------
C      Reflected downwelling thermal radiance
C      --------------------------------------

       iNalli = +1
       iNalli = -1   !! default Scott

       IF (iNalli .LT. 0) THEN      
         !!! this is default, Scott's version
         F=1.0
         IF (TAUZ(I) .GT. 0.0005) THEN
            F =   COEFF(1,I) +
     $         ( COEFF(2,I)/SECANG(LBOT) ) +
     $         ( COEFF(3,I)*TAUZ(I) ) +
     $         ( COEFF(4,I)*TAUZ(I)*TAUZ(I) ) +
     $         ( COEFF(5,I)*TAUZ(I)/SECANG(LBOT) ) +
     $         ( COEFF(6,I)*TAUZ(I)/RDOWN )
            !  Truncate F at limits as needed
            F = MAX( MIN(F,2.09), 0.696 )
         ENDIF
         RTHERM = RHOTHR(I)*PI*RDOWN*F*TAUZ(I)
       ELSE
         !  ccc uncomment these next few lines for nick nalli 2023 IEEE paper
         !  ccc see setems.f               RHOTHR(I)=(1.0 - EMIS(I))/PI
         !  ccc    so (1.0 - EMIS(I)) = RHOTHR(I)*PI
         !        ccc   and F = 1.0     !! for nick nalli
         !  ccc uncomment these next few lines for nick nalli 2023 IEEE paper
         F = 1.0
         RTHERM = RHOTHR(I)*PI*RDOWN*TAUZ(I)  !!! uncomment this for nick nalli
         IF (I .EQ. 1520) THEN
           ! Basically I would want the variables used in the RTE (Ts, satzen, Rdown, Tauz, F, Emis, etc.). 
           write(*,'(2(I8),11(F12.5,1X))') IPROF,I,FREQ(I),TSURF,EMIS(I),SECANG(LBOT),RHOTHR(I),RDOWN,F,TAUZ(I),RTHERM,RADUP,RSUN
         END IF
       END IF

C      --------------
C      Total radiance
C      --------------
       RAD0=RADUP + RSUN + RTHERM
c       print *,'sseerrggiioo ',RADUP,RSUN,RTHERM,RAD0
c       STOP

c     IF (I .EQ. 1100) print *,'calrad0 ',RADUP,RSUN,RTHERM,RAD0

       IF ((DOJAC) .AND. (DOSUN)) RTHERM = RSUN + RTHERM
C
C
       RETURN
       END
