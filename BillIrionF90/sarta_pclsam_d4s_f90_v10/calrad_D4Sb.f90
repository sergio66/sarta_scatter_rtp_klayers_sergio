!=======================================================================
!
!    University of Maryland Baltimore County [UMBC]
!
!    AIRS
!
!    CALRAD0
!    Calculate channel radiance for a clear atmosphere above a surface
!
!F77====================================================================


!ROUTINE NAME:
!    CALRAD_D4S


!ABSTRACT:
!    Calculate channel radiance for a clear atmosphere above a surface.


!CALL PROTOCOL:
!    CALRAD_D4S(DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
!       TAUL, TAUZ, SUNFAC, HSUN, TAUZSN, RHOSUN,
!       RHOTHR, LABOVE, COEFF, RPLNCK, RAD0 )


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    LOGICAL   DOSUN   do sun radiance calcs?      true/false
!    INTEGER   I       channel index               none
!    INTEGER   LBOT    bottom layer                none
!    REAL arr  RPLNCK  Planck function             mW/(m^2 cm^-1 sterad)
!    REAL      RSURFE  surface emission            mW/(m^2 cm^-1 sterad)
!    REAL arr  SECANG  path secant angles          none
!    REAL arr  TAUL    layer transmittance         none
!    REAL arr  TAUZ    surface-to-space trans      none
!    REAL      SUNFAC  sun solid angle * cosine    sterad
!    REAL arr  HSUN    solar irradiance at TOA     mW/(m^2 cm^-1 sterad)?
!    REAL arr  TAUZ    surface-to-space trans      none
!    REAL arr  TAUZSN  surface-to-space trans      none
!    REAL arr  RHOSUN  solar surface reflectivity  none
!    REAL arr  RHOTHR  down thermal surf refl      none
!    INT arr   LABOVE  down therm layer above      none
!    REAL arr  COEFF   "F" factor coefficients     various

!    REAL arr  TAU_D4S Transmissivity of cloud     none
!    INT       D4S_LAYER Layer number for cloud    none

!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL arr  RAD0    radiance                    mW/(m^2 cm^-1 sterad)


!INPUT/OUTPUT PARAMETERS:
!    none


!RETURN VALUES:
!    none


!PARENT(S):
!    sarta


!ROUTINES CALLED:
!    none


!FILES ACCESSED:
!    incFTC.f : include file of parameter statements accessed during
!       compilation only.


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    Calculates the radiance for an atmosphere with no clouds.


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    The temperature is treated a constant within each layer (ie
!    no adjustments for temperature gradiants).



!ROUTINE HISTORY:
!    Date        Programmer     Comments
!    ----------- -------------- ----------------------------------------
!    13 Jan 2006 Scott Hannon   Created
!    29 Mar 2006 Scott Hannon   Updated RTHERM for sartaV107
!    04 Dec 2017 Bill Irion     Converted from CALRAD0 to use opt depths
!                               from Delta-4-stream cloud calculations


SUBROUTINE CALRAD_D4S( &
	DOSUN, &
	I, &
	LBOT, &
	RPLNCK, &
	RSURFE, &
	SECANG, &
	TAUL, &
	TAUZ, &
	TAU_D4S, &
	D4S_LAYER, &
	SUNFAC, &
	HSUN, &
	TAUZSN, &
	RHOSUN, &
	RHOTHR, &
	LABOVE, &
	COEFF, &
	RAD0 )

	USE INCFTC
	IMPLICIT NONE


	!
	! ARGUMENTS
	!
	! Input
	LOGICAL, INTENT(IN) ::  DOSUN      ! do sun radiance calcs?
	INTEGER, INTENT(IN) ::      I      ! channel index
	INTEGER, INTENT(IN) ::   LBOT      ! bottom layer of atmosphere
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: RPLNCK ! layer Planck function
	REAL, INTENT(IN) :: RSURFE         ! surface emission
	REAL, INTENT(IN), DIMENSION(MAXLAY) :: SECANG ! viewing angle secant
	REAL, INTENT(IN), DIMENSION(MAXLAY) ::   TAUL ! clear air layer transmittance
	REAL, INTENT(IN), DIMENSION(MXCHAN) ::   TAUZ ! clear air surface-to-space transmittance
	!      Sun info
	REAL, INTENT(IN) :: SUNFAC         ! sun solid angle times cosine at surface
	REAL, INTENT(IN), DIMENSION(MXCHAN) :: HSUN   ! irradiance from Sun at top of atmosphere
	REAL, INTENT(IN), DIMENSION(MXCHAN) :: TAUZSN ! up plus down clear air solar transmittance
	REAL, INTENT(IN), DIMENSION(MXCHAN) :: RHOSUN ! surface reflectivity for solar
	!      Downwelling thermal info
	REAL, INTENT(IN), DIMENSION(MXCHAN) :: RHOTHR! surface reflectivity for downwelling thermal
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: LABOVE ! representative layer above surface
	REAL, INTENT(IN), DIMENSION(NFCOEF, MXCHAN) :: COEFF ! "F" factor coefficients

	REAL, INTENT(IN), DIMENSION(MXCHAN) :: TAU_D4S ! transmittance of cloud
	INTEGER, INTENT(IN) :: D4S_LAYER  ! layer number of cloud

	!
	!      Output
	REAL   RAD0                ! upwelling radiance at satellite

	!
	! LOCAL VARIABLES
	!
	INTEGER :: L           ! layer index
	INTEGER :: LTHERM      ! layer for RTHERM calc
	REAL :: F              ! reflected therm "F" (fudge) factor
	REAL :: RADUP          ! upward radiance
	REAL :: RSUN           ! reflected solar radiance
	REAL :: RTHERM         ! reflected downwelling thermal radiance


	! Copies are made of TAUZ and TAUL so as not to change them from the 
	! calling program
	REAL, DIMENSION(MAXLAY) :: TAUL_CLOUDY !  air layer transmittance with cloudy layer 
	REAL, DIMENSION(MXCHAN) :: TAUZ_CLOUDY !  air surface-to-space transmittance with cloudy layer
	REAL, DIMENSION(MXCHAN) :: TAUZSN_CLOUDY ! up plus down air solar transmittance with cloudy layer


	! Downwelling atmospheric thermal emission terms
	REAL TDOWNN ! "near-side" layer-to-surface trans
	REAL TDOWNF ! "far-side" layer-to-surface trans
	REAL  RDOWN ! downward radiance

	!
	!
	! EXECUTABLE CODE
	!
	!

	! Adjust transmissivity for cloudy layer
	TAUL_CLOUDY(:) = TAUL(:)
	TAUZ_CLOUDY(:) = TAUZ(:)
	TAUZSN_CLOUDY(:) = TAUZSN(:)

	TAUL_CLOUDY(D4S_LAYER) = TAUL(D4S_LAYER) * TAU_D4S(I)
	TAUZ_CLOUDY(I) = TAUZ(I) * TAU_D4S(I)
	TAUZSN_CLOUDY(I) = TAUZSN(I) * TAU_D4S(I)

	! 
	! Loop upward over layers
	!
	RADUP = RSURFE
	RDOWN = 0.0
	TDOWNN = 1.0
	DO L = LBOT, 1, -1
		RADUP = RADUP * TAUL_CLOUDY(L) + RPLNCK(L) * (1.0 - TAUL_CLOUDY(L))

		! Calc the downward radiance from this layer
		TDOWNF = TDOWNN * TAUL_CLOUDY(L)
		RDOWN = RDOWN + (RPLNCK(L) * (TDOWNN - TDOWNF) )
		TDOWNN = TDOWNF
	ENDDO
	
	!
	! Reflected solar radiance
	!
	IF (DOSUN) THEN
		RSUN=RHOSUN(I)*SUNFAC*HSUN(I)*TAUZSN_CLOUDY(I)
	ELSE
		RSUN=0.0
	ENDIF

	!
	! Reflected downwelling thermal radiance
	!
	F=1.0
	IF (TAUZ_CLOUDY(I) .GT. 0.0005) THEN
		F = &
			COEFF(1,I) +  &
			( COEFF(2,I)/SECANG(LBOT) ) + &
			( COEFF(3,I)*TAUZ_CLOUDY(I) ) + &
			( COEFF(4,I)*TAUZ_CLOUDY(I)*TAUZ_CLOUDY(I) ) + &
			( COEFF(5,I)*TAUZ_CLOUDY(I)/SECANG(LBOT) ) + &
			( COEFF(6,I)*TAUZ_CLOUDY(I)/RDOWN )
		! Truncate F at limits as needed
		F = MAX( MIN(F,2.09), 0.696 )
	ENDIF
	RTHERM=RHOTHR(I)*PI*RDOWN*F*TAUZ_CLOUDY(I)


	! 
	! Total radiance
	! 
	RAD0=RADUP + RSUN + RTHERM


	RETURN

END SUBROUTINE CALRAD_D4S
