!=======================================================================
!
!    University of Maryland Baltimore Country (UMBC)
!
!    AIRS
!
!    MIE_PARAM (was RDCLDT)
!
!F90====================================================================


!ROUTINE NAME:
!    RDCLDT


!ABSTRACT:
!    Read in the cloud lookup tables and save the parameters


!CALL PROTOCOL
!    RDCLDT ( IOUN, INDCHN, MIETYP, FNMIEA, FNMIEE, FNMIEG,
!       MIENPS, MIEPS, MIEABS, MIEEXT, MIEASY )


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INTEGER   IOUN    I/O unit number             none
!    INT arr   INDCHN  indices of channels         none
!    INT arr   MIETYP  Mie particle type code #    none
!    CHAR arr  FNMIEA  Mie absorption filenames    none
!    CHAR arr  FNMIEE  Mie extinction filenames    none
!    CHAR arr  FNMIEG  Mie asymmetry filenames     none

!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INT arr   MIENPS  Mie # of particle sizes     none
!    REAL arr  MIEPS   Mie cloud particle size     um
!    REAL arr  MIEABS  Mie cloud absorption        m^2/g
!    REAL arr  MIEEXT  Mie cloud extinction        m^2/g
!    REAL arr  MIEASY  Mie cloud asymmetry         none


!INPUT/OUTPUT PARAMETERS:
!    none


!RETURN VALUES:
!    none


!PARENT(S):
!    USEFAST


!ROUTINES CALLED:
!    none


!FILES ACCESSED:
!    incFTC.f : include file of parameter statements accessed during
!       compilation only.
!    unit IOUN : input file, binary FORTRAN data file. The file is
!       opened, read, and closed. This is done 3*NMIETY times.


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    May 2009 version of the cloudy SARTA code by L.Strow/S.Hannon.
!
!    Cloud lookup tables of absorption, extinction, and asymmetry
!    are read in from FORTRAN binary data files.


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
! Date        Programmer     Comments
! ----------- -------------- -------------------------------------------
! 12 May 2009 Scott Hannon   Created from rdcoef_pclsam.f
! 30 Aug 2017 Bill Irion     Converted from F77 subroutine to F90 module


!END====================================================================


SUBROUTINE RDCLDT( IOUN, INDCHN, MIETYP, FNMIEA, FNMIEE, FNMIEG, &
	MIENPS, MIEPS, MIEABS, MIEEXT, MIEASY )

	USE INCFTC
	IMPLICIT NONE

	! input
	INTEGER, INTENT(IN) ::   IOUN            	        ! I/O unit number
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDCHN	! channel use/index
	INTEGER, INTENT(IN), DIMENSION(NMIETY) :: MIETYP	! particle type code number
	CHARACTER (LEN = 79), INTENT(IN), DIMENSION(NMIETY) :: FNMIEA       ! absorption filenames
	CHARACTER (LEN = 79), INTENT(IN), DIMENSION(NMIETY) :: FNMIEE       ! extinction filenames
	CHARACTER (LEN = 79), INTENT(IN), DIMENSION(NMIETY) :: FNMIEG       ! asymmetry filenames

	! Output
	INTEGER, INTENT(OUT) :: MIENPS(NMIETY)            ! number of particle sizes
	REAL, INTENT(OUT), DIMENSION(MXMIEA,NMIETY) ::  MIEPS       ! particle size
	REAL, INTENT(OUT), DIMENSION(MXCHAN,MXMIEA,NMIETY) :: MIEABS ! absorption
	REAL, INTENT(OUT), DIMENSION(MXCHAN,MXMIEA,NMIETY) :: MIEEXT ! extinction
	REAL, INTENT(OUT), DIMENSION(MXCHAN,MXMIEA,NMIETY) :: MIEASY ! asymmetry

	! local
	REAL ::  RJUNK
	REAL ::  XJUNK(MXMIEA)
	INTEGER ::      I
	INTEGER ::     IC
	INTEGER ::   IERR
	INTEGER ::     IL
	INTEGER ::      J
	INTEGER ::      K

	! Read Mie lookup tables
	WRITE(6, *) "Reading Mie lookup tables"
	DO K=1, NMIETY ! NMIETY == number of mie particle types
		!
		! Read Mie absorption table
		!
		WRITE(6, *) "MIETYP ", MIETYP(K)
		WRITE(6, *) "READING ", FNMIEA(K) 
		OPEN(UNIT=IOUN, FILE=FNMIEA(K), FORM='UNFORMATTED', CONVERT = "BIG_ENDIAN", STATUS='OLD', IOSTAT=IERR)
		IF (IERR .NE. 0) THEN
			WRITE(6, 1020) IERR, FNMIEA(K)
1020		FORMAT('Error ',I5,' opening file:',/,A80, " in subroutine RDCLDT")
			STOP
		ENDIF
!
!		Read the number of channels and mie types
		READ(IOUN) I, IC
		IF (I .NE. MXCHAN) THEN
			WRITE(6,1071) FNMIEA(K), I, MXCHAN	
1071		FORMAT('ERROR! unexpected number of channels in mie file', A80,'  File has', I10,' channels, but MXCHAN=',I5)
			STOP
		ENDIF
		IF (IC .GT. MXMIEA) THEN
			WRITE(6,1072) FNMIEA(K), IC, MXMIEA
1072		FORMAT('ERROR! too many particle sizes in mie file', A80,'  File has', I5,' sizes, but MXMIEA=',I3)
			STOP
		ENDIF
		IF (IC .LT. 1) THEN
			WRITE(6,1073) FNMIEA(K), IC
1073		FORMAT('ERROR! invalid # of particle sizes in mie file', A80,' File has',I12,' sizes')
			STOP
		ENDIF
		MIENPS(K)=IC

		! Read mie particle sizes
		READ(IOUN) (MIEPS(IL,K), IL=1, MIENPS(K))
		
		! Read mie abs data for required channels
		
		J = 0  ! FWI 2/20/18
		DO I=1, MXCHAN
			! INDCHN(I) = I
			IF (INDCHN(I) .NE. 0) THEN
				! IF (.TRUE.) THEN
				! first index is channel. Second is for particle size. Third is lookup table.
				!PRINT *, "RDCLDT I, INDCHN(I) ", I, INDCHN(I) 
				READ(IOUN) (MIEABS(I, IL , K), IL=1, MIENPS(K))
				!READ(IOUN) (MIEABS(J, IL , K), IL=1, MIENPS(K))
				!PRINT *, "RDCLDT I, INDCHN(I) K MIEABS(INDCHN(I), IL , K)", I, INDCHN(I), K, (MIEABS(J, IL , K), IL=1, MIENPS(K))
				!J = J + 1  ! FWI 2/20/18
			ELSE
				READ(IOUN) (XJUNK(IL), IL=1, MIENPS(K))
			ENDIF
		ENDDO

		CLOSE(IOUN)

		!
		! Read Mie extinction table
		!
		WRITE(6, *) "READING ", FNMIEE(K) 
		OPEN(UNIT=IOUN, FILE=FNMIEE(K), FORM='UNFORMATTED', CONVERT = "BIG_ENDIAN", STATUS='OLD', IOSTAT=IERR)
		IF (IERR .NE. 0) THEN
			WRITE(6,1020) IERR, FNMIEE(K)
 			STOP
		ENDIF

		! Read the number of channels and mie types
		READ(IOUN) I, IC
		IF (I .NE. MXCHAN) THEN
			WRITE(6,1071) FNMIEE(K), I, MXCHAN
			STOP
		ENDIF
		IF (IC .NE. MIENPS(K)) THEN
			WRITE(6,1074) FNMIEE(K), IC, MIENPS(K)
1074		FORMAT('ERROR! unexpected # of particle sizes in mie file', A80, ' File has',I12,' sizes, but MIENPS=',I3)
			STOP
		ENDIF
		!
		! Read mie particle sizes
		!
		READ(IOUN) (XJUNK(I),I=1,MIENPS(K))
		! Check particle sizes are consistent
		DO I=1,MIENPS(K)
			RJUNK=ABS( MIEPS(I,K) - XJUNK(I) )
			IF (RJUNK .GT. 1E-4*MIEPS(I,K)) THEN
				WRITE(6,1075) FNMIEE(K)
1075			FORMAT('ERROR! unexpected particle sizes in mie file', /,A80)
				STOP
			ENDIF
		ENDDO
		!
		! Read mie ext data for required channels
		!
!		J=1 ! FWI 2/20/18
		DO I=1,MXCHAN
			IF (INDCHN(I) .NE. 0) THEN
				READ(IOUN) (MIEEXT(INDCHN(I),IL,K),IL=1,MIENPS(K)) ! FWI 2/20/18
				!READ(IOUN) (MIEEXT(J,IL,K),IL=1,MIENPS(K)) ! FWI 2/20/18
				!J = J + 1  ! FWI 2/20/18
			ELSE
				READ(IOUN) (XJUNK(IL),IL=1,MIENPS(K))
			ENDIF
		ENDDO

		CLOSE(IOUN)


		!
		! Read Mie asymmetry ("g") table
		!
		WRITE(6, *) "READING ", FNMIEG(K) 
		OPEN(UNIT=IOUN, FILE=FNMIEG(K), FORM='UNFORMATTED', CONVERT = "BIG_ENDIAN", STATUS='OLD', IOSTAT=IERR)
		IF (IERR .NE. 0) THEN
			WRITE(6,1020) IERR, FNMIEG(K)
			STOP
		ENDIF

		! Read the number of channels and mie types
		READ(IOUN) I, IC
		IF (I .NE. MXCHAN) THEN
			WRITE(6,1071) FNMIEG(K), I, MXCHAN
			STOP
		ENDIF
		IF (IC .NE. MIENPS(K)) THEN
			WRITE(6,1074) FNMIEG(K), IC, MIENPS(K)
			STOP
		ENDIF

		! Read mie particle sizes
		READ(IOUN) (XJUNK(I),I=1,MIENPS(K))
		! Check particle sizes are consistent
		DO I=1,MIENPS(K)
			RJUNK=ABS( MIEPS(I,K) - XJUNK(I) )
			IF (RJUNK .GT. 1E-4*MIEPS(I,K)) THEN
				WRITE(6,1075) FNMIEG(K)
				STOP
			ENDIF
		ENDDO

		! Read mie data for required channels
!		J=1 ! FWI 2/20/18
		DO I=1,MXCHAN
			IF (INDCHN(I) .NE. 0) THEN
				READ(IOUN) (MIEASY(I,IL,K),IL=1,MIENPS(K)) ! FWI 2/20/18
				!READ(IOUN) (MIEASY(J,IL,K),IL=1,MIENPS(K)) ! FWI 2/20/18
				!J = J + 1 ! FWI 2/20/18
			ELSE
				READ(IOUN) (XJUNK(IL),IL=1,MIENPS(K))
			ENDIF
		ENDDO

		CLOSE(IOUN)

	ENDDO

	RETURN
       
END SUBROUTINE RDCLDT
