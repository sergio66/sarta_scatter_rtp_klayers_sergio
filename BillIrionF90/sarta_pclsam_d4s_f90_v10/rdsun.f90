! version2 with hardcoded include filename instead of prompt
!=======================================================================
!=======================================================================
!
!    University of Maryland Baltimore Country (UMBC)
!
!    AIRS
!
!    RDSUN
!
!F77====================================================================


!ROUTINE NAME:
!    RDSUN


!ABSTRACT:
!    Read in the AIRS solar radiance data


!CALL PROTOCOL
!    RDSUN ( IOUN, INDCHN, HSUN )


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INT arr   INDCHN  indices of channels         none
!    INTEGER   IOUN    I/O unit number             none


!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL arr  HSUN    solar radiance              W/(m2.str.cm-1)


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
!    unit IOUN : input file, ASCII text file. The file is opened,
!       read, and closed.


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    August 2000 version of the 100 layer AIRS Fast Transmittance
!    Code by L.Strow/S.Hannon.
!
!    Reads in a text file with solar radiance data for each AIRS
!    channel.  This is the solar rad direct from the sun at the top
!    of Earth's atmosphere.


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
!    Date        Programmer     Comments
!    ----------- -------------- ----------------------------------------
!    12 Sep 1997 Scott Hannon   Created
!    12 Feb 2001 Scott Hannon   hardcoded filename instead of prompt
!    19 Dec 2017 Bill Irion     Conversion to F90
!    19 Dec 2017 Bill Irion     Removed GOTO's

!END====================================================================

SUBROUTINE RDSUN ( IOUN, INDCHN, HSUN )

	USE INCFTC
	IMPLICIT NONE

	!
	! ARGUMENTS
	!
	! Input
	INTEGER, INTENT(IN) ::   IOUN
	INTEGER, INTENT(IN), DIMENSION(MXCHAN) :: INDCHN
	!
	! Output
	REAL, INTENT(OUT), DIMENSION(MXCHAN) :: HSUN


	!
	! LOCAL VARIABLES
	!
	CHARACTER*80 :: CLINE
	INTEGER :: I
	INTEGER :: IERR
	INTEGER :: ICHAN
	REAL :: FRQCHN
	REAL :: SUNCHN

	!
	!
	! EXECUTABLE CODE
	!

	!
	!
	! Open the solar radiance file
	!
	WRITE(*, *) "READING: ", FNSUN 
	OPEN(UNIT=IOUN,FILE=FNSUN,FORM='FORMATTED',STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNSUN
1020 	FORMAT('Error ',I5,' openning file:',/,A80)
		STOP
	ENDIF

	! Initialize the channel counter
	I=0

	!
	! Read the solar rad file
	!
	DO
		! Read a line of text from the file
		READ(IOUN, "(A80)", END=910) CLINE
	
		! Determine if the text line is data or a comment
		IF (CLINE(1:1) .NE. '!') THEN

			! It's data, so increment the channel counter
			I=I+1

			! Read the data from the text line
			READ(CLINE,*)  ICHAN, FRQCHN, SUNCHN

			! Check to be sure the channel value is OK
			IF ((ICHAN .LT. 1) .OR. (ICHAN .GT. MXCHAN)) THEN
				WRITE(6,1040) MXCHAN, ICHAN
1040			FORMAT('Error! Channel number is out of range. Range is 1 to ',I4,', but channel list has ',I7,'.')
				STOP
			ENDIF

			! Keep the data if the current channel is on the list
			IF (INDCHN(ICHAN) .NE. 0) THEN
				HSUN( INDCHN(ICHAN) )=SUNCHN
			ENDIF

		ENDIF

		

		! Exit if the next line is > MXCHAN
		IF (I .GE. MXCHAN) EXIT
		! Note: this routine expects data for every channel
	END DO

910 CLOSE(IOUN)

	RETURN

END SUBROUTINE RDSUN
