MODULE SECANG_MODULE
	
	USE INCFTC
	IMPLICIT NONE

	REAL, DIMENSION(MAXLAY) :: SECANG
	REAL, DIMENSION(MAXLAY) :: SECSUN
	LOGICAL :: DOSUN

CONTAINS

	SUBROUTINE SECANG_MODULE__CALCULATE_SECANG( &
		SATZEN, &
		SATANG, &		
		SUNANG, &
		SALT, &
		ALT, &
		LBOT)
		
		REAL, INTENT(IN) :: SATZEN
		REAL, INTENT(IN) :: SATANG
		REAL, INTENT(IN) :: SUNANG
		REAL, INTENT(IN) :: SALT
		REAL, INTENT(IN), DIMENSION(MAXLAY) :: ALT
		INTEGER, INTENT(IN) :: LBOT

		REAL, PARAMETER :: CONV = 1.7453292E-02  ! pi/180 = degrees to radians conversion factor
		REAL, PARAMETER :: ANGMAX = 53.  ! max satellite view angle (49.5 scan + 3.5 spacecraft)

		REAL :: SVA
		REAL :: EVA
		REAL :: SZALAY
		REAL :: SCOS1
		REAL :: SUNCOS
		REAL :: RJUNK1, RJUNK2
		REAL :: SUNFDG
		INTEGER :: L 
		
		REAL :: VACONV
		REAL :: SACONV


		! Convert SATZEN or SATANG to viewing angle
		IF (SATZEN .GE. 0.0 .AND. SATZEN .LT. 63.0) THEN
			! Convert zenith angle at surface to view angle at satellite
			SVA = SACONV( SATZEN, SALT*1000.0 )/ CONV
		ELSE
			! Check if scan angle is valid
			IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
				! View angle should be within a few degrees of scan angle
				SVA=ABS( SATANG )
			ELSE
				WRITE(IOERR,1030) SATZEN, SATANG
				1030 FORMAT('Error! invalid angles for SATZEN ',1PE11.4,' and SATANG ',E11.4) 
				STOP
			ENDIF
		ENDIF

		!ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
		IF (SVA .GT. ANGMAX) THEN
			! Truncate angle if too big
			WRITE(IOINFO, 1040) SVA
			1040 FORMAT('Warning! Profile: truncating view angle ',1PE11.4,' to 53 degrees')
			SVA=ANGMAX
		ENDIF

		! Convert from satellite to earth viewing angle (in radians)
		DO L=1,LBOT
			EVA=VACONV(SVA, SALT, ALT(L))
			SECANG(L)=1.0E+0/COS(EVA)
			! for testing
			! SECANG(L)=SVA
		ENDDO

		! Calc total sun angle secant
		DOSUN=.FALSE.
		IF (SUNANG .GE. 0.0 .AND. SUNANG .LT. 89.9) DOSUN=.TRUE.
			IF (DOSUN) THEN
				SUNCOS=COS(CONV*SUNANG)
				SZALAY=SACONV(SUNANG,ALT(1))
				SCOS1=COS(SZALAY)
				RJUNK2=SECANG(LBOT) + 1.0/SUNCOS ! Total secant

				! Calc non-unity fudge factor if total secant > 9
				IF (RJUNK2 .GT. 9.0) THEN
					! fudge factor = true_total_secant/calc_total_secant
					SUNFDG=RJUNK2/9.0
					! truncated solar angle to use to calc SECSUN
					RJUNK1=ACOS( 1.0/(9.0 - SECANG(LBOT)) )/CONV
				ELSE
					SUNFDG=1.0
					RJUNK1=SUNANG
				ENDIF
				! Should I change SUNFDG to SUNFDG(MAXLAY)?
				DO L=1,LBOT
					SZALAY=SACONV(RJUNK1,ALT(L))
					SECSUN(L)=SECANG(L) + 1.0E+0/COS(SZALAY)
				ENDDO
		ENDIF

	END SUBROUTINE SECANG_MODULE__CALCULATE_SECANG

END MODULE SECANG_MODULE
