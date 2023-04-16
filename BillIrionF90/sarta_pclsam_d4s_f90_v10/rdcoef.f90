!=======================================================================
!
!    University of Maryland Baltimore Country (UMBC)
!
!    AIRS
!
!    RDCOEF version with trace gases CO2, SO2, & HNO3
!
!F77====================================================================


!ROUTINE NAME:
!    RDCOEF


!ABSTRACT:
!    Read in the AIRS fast transmittance coefficients.


!CALL PROTOCOL
!    RDCOEF ( IOUN, NCHAN, INDCHN, SETCHN,
!       NCHN1, NCHN2, NCHN3, NCHN4, NCHN5, NCHN6, NCHN7,
!       CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
!       COEF1, COEF2, COEF3, COEF4, COEF5, COEF6, COEF7,
!       FREQ, LABOVE, COEFF, INDCO2, COFCO2, INDSO2, COFSO2,
!       INDHNO, COFHNO, INDN2O, COFN2O,
!       INDH2O, WAZOP, WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN )


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INT arr   INDCHN  indices of channels         none
!    INTEGER   IOUN    I/O unit number             none
!    INTEGER   NCHAN   number of channels          none


!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INT arr   CLIST1  set1 channel list           none
!    INT arr   CLIST2  set2 channel list           none
!    INT arr   CLIST3  set3 channel list           none
!    INT arr   CLIST4  set4 channel list           none
!    INT arr   CLIST5  set5 channel list           none
!    INT arr   CLIST6  set6 channel list           none
!    INT arr   CLIST7  set7 channel list           none
!    INT arr   CLISTN  non-LTE channel list        none
!    REAL arr  COEF1   set1 fast trans coefs       various
!    REAL arr  COEF2   set2 fast trans coefs       various
!    REAL arr  COEF3   set3 fast trans coefs       various
!    REAL arr  COEF4   set4 fast trans coefs       various
!    REAL arr  COEF5   set5 fast trans coefs       various
!    REAL arr  COEF6   set6 fast trans coefs       various
!    REAL arr  COEF7   set7 fast trans coefs       various
!    REAL arr  COEFF   thermal "F" factor coefs    various
!    REAL arr  COEFN   non-LTE coefficients        various
!    REAL arr  COFCO2  CO2 perturbation coefs      various
!    REAL arr  COFSO2  SO2 perturbation coefs      various
!    REAL arr  COFHNO  HNO3 perturbation coefs     various
!    REAL arr  COFN2O  N2O perturbation coefs      various
!    REAL arr  COFH2O  OPTRAN H2O coefs            various
!    REAL arr  FREQ    channel freqs               cm-1
!    REAL arr  FX      fixed gases adjustment      none
!    INT arr   INDCO2  CO2 pert channel indices    none
!    INT arr   INDSO2  SO2 pert channel indices    none
!    INT arr   INDHNO  HNO3 pert channel indices   none
!    INT arr   INDN2O  N2O pert channel indices    none
!    INT arr   INDH2O  OPTRAN H2O channel indices  none
!    INT arr   LABOVE  layer above for thermal     none
!    INTEGER   NCHN1   set1 number of channels     none
!    INTEGER   NCHN2   set2 number of channels     none
!    INTEGER   NCHN3   set3 number of channels     none
!    INTEGER   NCHN4   set4 number of channels     none
!    INTEGER   NCHN5   set5 number of channels     none
!    INTEGER   NCHN6   set6 number of channels     none
!    INTEGER   NCHN7   set7 number of channels     none
!    INTEGER   NCHNTE  non-LTE number of channels  none
!    REAL arr  WAZOP   OPTRAN water grid           kiloMoles/cm^2
!    REAL arr  WAVGOP  OPTRAN water pred averges   various
!    INT arr   SETCHN  set# (1-7) chan belongs to  none (integer, 1 - 7)


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
!       opened, read, and closed. This is done 10 times, once per
!       each of the 7 coef sets, and once each for the variable CO2,
!       OPTRAN water, and thermal F factor coefs.


!COMMON BLOCKS
!    none


!DESCRIPTION:
!    June 2005 version of the 100 layer AIRS Fast Transmittance
!    Code by L.Strow/S.Hannon.
!
!    Seven sets of binary data files containing the main fast
!    transmittance coefficients are opened and read one channel at
!    a time.  The seven sets of coefs are each stored in their own
!    arrays.  Next, preturbation coefficients for four trace gases,
!    (CO2, SO2, HNO3, & N2O) are read in from four binary files.
!    Next, OPTRAN water fast trans coefs for some channels are read
!    in from a binary file file.  The header of the OPTRAN file
!    specifies 300 OPTRAN water levels, and also the mean value of
!    4 predictor terms for each of the levels.  Next, comes the
!    read of the binary file with the reflected downwelling thermal
!    radiance "F factor" coefficients.  Next is a read of the "FX"
!    fixed gases adjustment term from an ASCII text file.  Lastly
!    comes the read of the non-LTE coefficients from a binary file.


!ALGORITHM REFERENCES:
!    none


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
!    Date        Programmer     Comments
!    ----------- -------------- ----------------------------------------
!    Dec  1 1994 Scott Hannon   Created
!    Dec 21 1994 Scott Hannon   Fixed error with IOPF (now assigned)
!     5 Feb 1997 Scott Hannon   Re-wrote for FWO+FOW+FMW+FCOW.
!    28 Aug 1997 Scott Hannon   Re-wrote for sets 1 - 7 and thermal
!    30 Sep 1997 Scott Hannon   Added COFCO2 and INDCO2
!    27 Feb 1998 Scott Hannon   Added COFH2O, INDH2O, WAZOP, & WAVGOP
!    17 Aug 2000 Scott Hannon   Add FX
!    12 Feb 2001 Scott Hannon   hardcoded filenames instead of prompts
!    18 May 2005 Scott Hannon   Add HNO3 based on SO2 code
!    28 Jun 2005 Scott Hannon   "trace" version for CO2,SO2,HNO3,N2O
!    13 Oct 2005 Scott Hannon   Add non-LTE variables
!    15 Sep 2107 Bill Irion     Conversion to F90

!END====================================================================

!=================================================================
SUBROUTINE RDCOEF ( &
	  IOUN,  NCHAN, INDCHN, SETCHN, &
	 NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7, &
    CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7, &
	 COEF1,  COEF2,  COEF3,  COEF4,  COEF5,  COEF6,  COEF7, &
	  FREQ, LABOVE,  COEFF, INDCO2, COFCO2, INDSO2, COFSO2, &
	INDHNO, COFHNO, INDN2O, COFN2O, &
	INDH2O,  WAZOP, WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN )
!=================================================================


	USE INCFTC
	IMPLICIT NONE

!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
!	Input
	INTEGER, INTENT(IN) :: IOUN
	INTEGER, INTENT(IN) ::  NCHAN
	INTEGER, INTENT(IN) :: INDCHN(MXCHAN)


	INTEGER ::  SETCHN(MXCHAN) ! set # for each channel
	INTEGER ::   NCHN1         ! # of set1 channels
	INTEGER ::   NCHN2         ! # of set2 channels
	INTEGER ::   NCHN3         ! # of set3 channels
	INTEGER ::   NCHN4         ! # of set4 channels
	INTEGER ::   NCHN5         ! # of set5 channels
	INTEGER ::   NCHN6         ! # of set6 channels
	INTEGER ::   NCHN7         ! # of set7 channels
	INTEGER ::  CLIST1(MXCHN1) ! list of set1 channels
	INTEGER ::  CLIST2(MXCHN2) ! list of set2 channels
	INTEGER ::  CLIST3(MXCHN3) ! list of set3 channels
	INTEGER ::  CLIST4(MXCHN4) ! list of set4 channels
	INTEGER ::  CLIST5(MXCHN5) ! list of set5 channels
	INTEGER ::  CLIST6(MXCHN6) ! list of set6 channels
	INTEGER ::  CLIST7(MXCHN7) ! list of set7 channels
	REAL ::   COEF1(N1COEF,MAXLAY,MXCHN1) ! coefs for set1 chans
	REAL ::   COEF2(N2COEF,MAXLAY,MXCHN2) ! coefs for set2 chans
	REAL ::   COEF3(N3COEF,MAXLAY,MXCHN3) ! coefs for set3 chans
	REAL ::   COEF4(N4COEF,MAXLAY,MXCHN4) ! coefs for set4 chans
	REAL ::   COEF5(N5COEF,MAXLAY,MXCHN5) ! coefs for set5 chans
	REAL ::   COEF6(N6COEF,MAXLAY,MXCHN6) ! coefs for set6 chans
	REAL ::   COEF7(N7COEF,MAXLAY,MXCHN7) ! coefs for set7 chans
	REAL ::    FREQ(MXCHAN)    ! chan center frequency
	INTEGER ::  LABOVE(MXCHAN) ! chan downwelling thermal layer above
	REAL ::   COEFF(NFCOEF,MXCHAN)        ! coefs for chan "F" factor
	INTEGER ::  INDCO2(MXCHAN)            ! chan indices for CO2 pert
	REAL ::  COFCO2(  NCO2,MAXLAY,MXCHNC) ! coefs for CO2 pert
	INTEGER ::  INDSO2(MXCHAN)            ! chan indices for SO2 pert
	REAL ::  COFSO2(  NSO2,MAXLAY,MXCHNS) ! coefs for SO2 pert
	INTEGER ::  INDHNO(MXCHAN)            ! chan indices for HNO3 pert
	REAL ::  COFHNO( NHNO3,MAXLAY,MXCHNH) ! coefs for HNO3 pert
	INTEGER ::  INDN2O(MXCHAN)            ! chan indices for N2O pert
	REAL ::  COFN2O(  NN2O,MAXLAY,MXCHNN) ! coefs for N2O pert
	INTEGER ::  INDH2O(MXCHAN)            ! chan indices for OPTRAN H2O
	REAL ::    WAZOP(MXOWLY)              ! OPTRAN water l-to-s amounts
	REAL ::   WAVGOP(NOWAVG,MXOWLY)       ! OPTRAN raw predictor averages
	REAL ::  COFH2O(  NH2O,MXOWLY,MXCHNW) ! coefs for OPTRAN H2O
	REAL ::      FX(MAXLAY)               ! fixed gases adjustment
	INTEGER ::  NCHNTE                    ! number of non-LTE channels
	INTEGER ::  CLISTN(MXCNTE)            ! non-LTE channel list
	REAL ::   COEFN(NNCOEF,MXCNTE)        ! non-LTE coefficients


!-----------------------------------------------------------------------
!	LOCAL VARIABLES
!-----------------------------------------------------------------------
	CHARACTER*80  CLINE
	REAL FRQCHN
	REAL  FCHAN(NFCOEF)
	REAL  RJUNK
	INTEGER      I
	INTEGER     IC
	INTEGER  ICHAN
	INTEGER   IERR
	INTEGER     IL
	INTEGER      J
	INTEGER ICOUNT
	INTEGER LACHAN


	!
	!
	! EXECUTABLE CODE
	!
	!
	!

	! Initialize "set"-independent index arrays
	! Removed do-loop  FWI 12/15/17

	! Trace gases
	INDCO2(:)=0
	INDSO2(:)=0
	INDHNO(:)=0
	INDN2O(:)=0
	! OPTRAN water
	INDH2O(:)=0


	!----------
	! Read set 1
	!----------
	OPEN(UNIT=IOUN, FILE=FNCOF1, FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNCOF1
1020		FORMAT('Error ',I5,' opening file:',/,A80)
			STOP
	ELSE
		WRITE(6, *) "READING: ", FNCOF1
	ENDIF

	J=1
	DO I=1, MXCHN1
		! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COEF1(IC, IL, J), IC=1, N1COEF), IL=1, MAXLAY)

		! Set the set number for this channel
		SETCHN(ICHAN) = 1

		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			CLIST1(J) = ICHAN
			FREQ(INDCHN(ICHAN)) = FRQCHN
			J = J + 1
		ENDIF
	ENDDO
	NCHN1 = J - 1

	CLOSE(IOUN)

	!----------
	! Read set 2
	!----------
	OPEN(UNIT=IOUN, FILE=FNCOF2, FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNCOF2
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNCOF2
	ENDIF

	J = 1
	DO I = 1, MXCHN2
		! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COEF2(IC, IL, J), IC = 1, N2COEF), IL=1, MAXLAY)

		SETCHN(ICHAN) = 2

		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			CLIST2(J) = ICHAN
			FREQ(INDCHN(ICHAN)) = FRQCHN
			J = J + 1
		ENDIF
	ENDDO
	NCHN2 = J - 1

	CLOSE(IOUN)

	!----------
	! Read set 3
	! ----------
	OPEN(UNIT=IOUN, FILE=FNCOF3, FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNCOF3
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNCOF3
	ENDIF

	J = 1
	DO I = 1, MXCHN3
	! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COEF3(IC, IL, J), IC=1, N3COEF), IL=1, MAXLAY)

		SETCHN(ICHAN)=3

		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			CLIST3(J) = ICHAN
			FREQ(INDCHN(ICHAN)) = FRQCHN
			J=J + 1
		ENDIF
	ENDDO
	NCHN3 = J - 1

	CLOSE(IOUN)

	!----------
	! Read set 4
	!----------
	OPEN(UNIT=IOUN, FILE=FNCOF4, FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNCOF4
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNCOF4
	ENDIF

	J=1
	DO I=1,MXCHN4
		! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COEF4(IC,IL,J),IC=1,N4COEF), IL=1,MAXLAY)
		!PRINT *, ICHAN, INDCHN(ICHAN)

		SETCHN(ICHAN)=4
		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			CLIST4(J)=ICHAN
			FREQ( INDCHN(ICHAN) )=FRQCHN
			J=J + 1
		ENDIF
	ENDDO
	NCHN4=J - 1

	CLOSE(IOUN)

	!----------
	! Read set 5
	!----------
	OPEN(UNIT=IOUN, FILE=FNCOF5, FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNCOF5
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNCOF5
	ENDIF

	J=1
	DO I=1,MXCHN5
	! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COEF5(IC,IL,J),IC=1,N5COEF), IL=1,MAXLAY)
		!PRINT *, ICHAN, INDCHN(ICHAN)

		SETCHN(ICHAN)=5
		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			CLIST5(J)=ICHAN
			FREQ( INDCHN(ICHAN) )=FRQCHN
			J=J + 1
		ENDIF
	ENDDO
	NCHN5=J - 1

	CLOSE(IOUN)
	
	!----------
	! Read set 6
	!----------
	OPEN(UNIT=IOUN, FILE=FNCOF6, FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNCOF6
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNCOF6
	ENDIF

	J=1
	DO I=1, MXCHN6
		! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COEF6(IC,IL,J),IC=1,N6COEF), IL=1,MAXLAY)
		!PRINT *, ICHAN, INDCHN(ICHAN)
		SETCHN(ICHAN)=6
		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			CLIST6(J)=ICHAN
			FREQ( INDCHN(ICHAN) )=FRQCHN
			J=J + 1
		ENDIF
	ENDDO
	NCHN6=J - 1

	CLOSE(IOUN)

	!----------
	! Read set 7
	!----------
	OPEN(UNIT=IOUN, FILE=FNCOF7, FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNCOF7
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNCOF7
	ENDIF

	J=1
	DO I=1,MXCHN7
		! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COEF7(IC,IL,J),IC=1,N7COEF), IL=1,MAXLAY)
		!PRINT *, ICHAN, INDCHN(ICHAN)
		SETCHN(ICHAN)=7

		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			CLIST7(J)=ICHAN
			FREQ( INDCHN(ICHAN) )=FRQCHN
			J=J + 1
		ENDIF
	ENDDO
	NCHN7 = J - 1

	CLOSE(IOUN)


	!---------------------------
	! Read CO2 perturbation coefs
	!---------------------------
	OPEN(UNIT=IOUN,FILE=FNCO2, FORM = 'UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNCO2
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNCO2
	ENDIF

	J=1
	DO I=1,MXCHNC
	! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COFCO2(IC,IL,J),IC=1,NCO2), IL=1,MAXLAY)
		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			INDCO2(ICHAN)=J
			J=J + 1
		ENDIF
	ENDDO

	CLOSE(IOUN)

	!---------------------------
	! Read SO2 perturbation coefs
	!---------------------------
	OPEN(UNIT=IOUN,FILE=FNSO2,FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNSO2
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNSO2
	ENDIF

	J=1
	DO I=1,MXCHNS
		! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COFSO2(IC,IL,J),IC=1,NSO2), IL=1,MAXLAY)
		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			INDSO2(ICHAN)=J
			J=J + 1
		ENDIF
	ENDDO

	CLOSE(IOUN)

	! ---------------------------
	! Read HNO3 perturbation coefs
	! ---------------------------
	OPEN(UNIT=IOUN,FILE=FNHNO3,FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNHNO3
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNHNO3
	ENDIF

	J=1
	DO I=1,MXCHNH
		! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COFHNO(IC,IL,J),IC=1,NHNO3), IL=1,MAXLAY)

		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			INDHNO(ICHAN)=J
			J=J + 1
		ENDIF
	ENDDO

	CLOSE(IOUN)


	! ---------------------------
	! Read N2O perturbation coefs
	! ---------------------------
	OPEN(UNIT=IOUN,FILE=FNN2O,FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNN2O
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNN2O
	ENDIF

	J=1
	DO I=1,MXCHNN
	! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COFN2O(IC,IL,J),IC=1,NN2O), IL=1,MAXLAY)
		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			INDN2O(ICHAN)=J
			J=J + 1
		ENDIF
	ENDDO

	CLOSE(IOUN)

	! ---------------------
	! Read OPTRAN H2O coefs
	! ---------------------
	OPEN(UNIT=IOUN,FILE=FNOPTR,FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNOPTR
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNOPTR
	ENDIF

	READ(IOUN) (WAZOP(IL),IL=1,MXOWLY)
	DO IC=1,NOWAVG
	! Read the header section
		READ(IOUN) (WAVGOP(IC,IL),IL=1,MXOWLY)
	ENDDO

	J=1
	DO I=1,MXCHNW
		! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, ((COFH2O(IC,IL,J),IC=1,NH2O), IL=1,MXOWLY)
		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			INDH2O(ICHAN)=J
			J=J + 1
		ENDIF
	ENDDO

	CLOSE(IOUN)

	!-----------------------------------------------
	! Read the downward thermal F factor coefficients
	!-----------------------------------------------
	OPEN(UNIT=IOUN,FILE=FNTHER,FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNTHER
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNTHER
	ENDIF

	DO I=1,MXCHAN
		! Read data for this frequency/channel
		!cc changed 18 May 2005
		! 	READ(IOUN) ICHAN, FRQCHN, LACHAN, (FCHAN(IC),IC=1,NFCOEF)
		READ(IOUN) ICHAN, FRQCHN, (FCHAN(IC),IC=1,NFCOEF)
		LACHAN=-1   ! assign dummy value
		!cc

		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			LABOVE( INDCHN(ICHAN) )=LACHAN
			DO IC=1,NFCOEF
				COEFF(IC,INDCHN(ICHAN))=FCHAN(IC)
			ENDDO
		ENDIF
	ENDDO

	CLOSE(IOUN)


	!-------
	! Read FX
	!-------
	OPEN(UNIT=IOUN,FILE=FNFX,FORM='FORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNFX
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNFX
	ENDIF

	! Read the file
	ICOUNT=0
10	READ(IOUN,9000,END=99) CLINE
9000  FORMAT(A80)
	IF (CLINE(1:1) .NE. '!') THEN
	! Note: fx file format is:  layer_number  fx_value 
		READ(CLINE,*) IC, RJUNK
		ICOUNT=ICOUNT + 1
		FX(IC)=RJUNK
	ENDIF
	GOTO 10

99	CLOSE(IOUN)

	IF (ICOUNT .NE. MAXLAY) THEN
		WRITE(6,1047) MAXLAY, ICOUNT
1047	FORMAT('Error! Unexpected number of layers in fx file.',/, 'Expected fx to have ',I4,' layers, but found ',I4)
	ENDIF


	!------------
	! Read non-LTE
	!------------
	OPEN(UNIT=IOUN,FILE=FNCOFN,FORM='UNFORMATTED', CONVERT = 'BIG_ENDIAN', STATUS='OLD', IOSTAT=IERR)
	IF (IERR .NE. 0) THEN
		WRITE(6,1020) IERR, FNCOFN
		STOP
	ELSE
		WRITE(6, *) "READING: ", FNCOFN
	ENDIF

	J=1
	DO I=1,MXCNTE
		! Read data for this frequency/channel
		READ(IOUN) ICHAN, FRQCHN, (COEFN(IC,J),IC=1,NNCOEF)
		! Keep the data if the current channel is on the list
		IF (INDCHN(ICHAN) .NE. 0) THEN
			CLISTN(J)=ICHAN
			J=J + 1
		ENDIF
	ENDDO
	NCHNTE=J - 1

	CLOSE(IOUN)

	!---------------------------------------------
	! Make sure all channels on the list were found
	!---------------------------------------------
	ICOUNT=NCHN1 + NCHN2 + NCHN3 + NCHN4 + NCHN5 + NCHN6 + NCHN7
	IF (ICOUNT .NE. NCHAN) THEN
		WRITE(6,1050) NCHAN, ICOUNT
1050 	FORMAT('Error! Unexpected number of channels found.',/, 'The channel list had ',I4,' channels, but found ',I4)
	ENDIF

	!----------------------------
	! Show summary of channel sets
	!----------------------------
	!cc
	!       WRITE(6,1060) 1, NCHN1
	! 1060  FORMAT('Number of channels for set',I1,' = ',I4)
	!       WRITE(6,1060) 2, NCHN2
	!       WRITE(6,1060) 3, NCHN3
	!       WRITE(6,1060) 4, NCHN4
	!       WRITE(6,1060) 5, NCHN5
	!       WRITE(6,1060) 6, NCHN6
	!       WRITE(6,1060) 7, NCHN7
	!cc
	!
	RETURN
	
END SUBROUTINE RDCOEF
