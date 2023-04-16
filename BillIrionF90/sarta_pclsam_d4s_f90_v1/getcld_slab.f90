!=======================================================================
!
!            University of Maryland Baltimore County [UMBC]
!
!            AIRS
!
!            GETCLD
!
!F90====================================================================


!ROUTINE NAME: GETCLD


!ABSTRACT:
!    Get basic cloud info

!CALL PROTOCOL:
!    GETCLD( IPROF, HEAD, PROF,
!    LBLAC1, CTYPE1, CFRAC1, CPSIZ1, CPRTO1, CPRBO1, CNGWA1, ...
!    XCEMI1, XCRHO1, CSTMP1,
!    LBLAC2, CTYPE2, CFRAC2, CPSIZ2, CPRTO2, CPRBO2, CNGWA2,
!    XCEMI2, XCRHO2, CSTMP2, CFRA12, FCLEAR, CFRA1X, CFRA2X )


!INPUT PARAMETERS:
!    type      name    purpose                 units
!    --------  ------  --------------------------  ---------------------
!    INTEGER   IPROF   profile loop counter        none
!    STRUCT    HEAD    RTP header structure        various
!    STRUCT    PROF    RTP profile structure       various

!OUTPUT PARAMETERS:
!    type      name    purpose                 units
!    --------  ------  --------------------------  ---------------------
!    LOGICAL   LBLAC1  cloud1 is black?          none
!    INTEGER   CTYPE1  cloud type code number      none
!    REAL      CFRAC1  total cloud1 fraction       none (0.0 to 1.0)
!    REAL      CPSIZ1  particle size             um
!    REAL      CPRTO1  cloud top pressure        mb
!    REAL      CPRBO1  cloud bottom pressure       mb
!    REAL      CNGWA1  layer integrated profile    g/m^2
!    REAL arr  XCEMI1  emissivity cloud1         none (0.0 to 1.0)
!    REAL arr  XCRHO1  reflectivity cloud1         none
!    REAL      CSTMP1  cloud top temperature       Kelvin
!    LOGICAL   LBLAC2  cloud2 is black?          none
!    INTEGER   CTYPE2  cloud type code number      none
!    REAL      CFRAC2  total cloud2 fraction       none (0.0 to 1.0)
!    REAL      CPSIZ2  particle size             um
!    REAL      CPRTO2  cloud top pressure        mb
!    REAL      CPRBO2  cloud bottom pressure       mb
!    REAL      CNGWA2  layer integrated profile    g/m^2
!    REAL arr  XCEMI2  emissivity cloud2         none (0.0 to 1.0)
!    REAL arr  XCRHO2  reflectivity cloud2         none
!    REAL      CSTMP2  cloud top temperature       Kelvin
!    REAL      CFRA12  both clouds fraction        none (0.0 to 1.0)
!    REAL      FCLEAR  clear fraction            none (0.0 to 1.0)
!    REAL      CFRA1X  exclusive cloud1 fraction   none (0.0 to 1.0)
!    REAL      CFRA2X  exclusive cloud2 fraction   none (0.0 to 1.0)


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): SARTA


!ROUTINES CALLED: none


!FILES ACCESSED:
!    none


!COMMON BLOCKS: none


!DESCRIPTION:
!    Pulls out and checks cloud profile paramters.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
!    Date     Programmer        Comments
!------------ ----------------- ----------------------------------------
! 22 Feb 2007 Scott Hannon      created
! 10 Sep 2007 Scott Hannon      Mods for 100 layer particle size
! 14 Sep 2007 Scott Hannon      Add CPMIN check
! 15 Nov 2007 Scott Hannon      re-written for slab cloud
! 26 Nov 2008 Scott Hannon      Partial re-write for rtpV201 including
!                         spectral XCEMI1/2 & XCRHO1/2
! 01 Dec 2008 Scott Hannon      Add CSTMP1/2
! 07 May 2009 Scott Hannon      Bug fix: initialize CFRA1X=CFRA2X=0
! 30 Aug 2017 Bill Irion        Changed references to HEAD and PROF records 
!                               Converted from F77 to F90
! 12 Sep 2017 Bill Irion        Removed non-structured variables where appropriate
!                               Routine does not copy to a non-structured scalar not
!                               needed.

!END====================================================================
!-----------------------------------------------------------------------
!      INCLUDE FILES
!-----------------------------------------------------------------------
      !include 'rtpdefs.f'

!=================================================================
SUBROUTINE GETCLD( &
	CLOUD,  &  ! cloud input parameters formerly in PROF structure
	NEMIS,  &  ! length of emissivity array -- formerly in PROF structure
	PLEVS,  &  ! pressure levels -- formerly in PROF structure
	PSURF,  &  ! surface pressure -- formerly in PROF structure
	LBLAC1, &
!	CTYPE1, &
!	CFRAC1, &
!	CPSIZ1, &
!	CPRTO1, &
!	CPRBO1, &
!	CNGWA1, &
!	XCEMI1, &
!	XCRHO1, &
!	CSTMP1, &
!	LBLAC2, &
!	CTYPE2, &
!	CFRAC2, &
!	CPSIZ2, &
!	CPRTO2, &
!	CPRBO2, &
!	CNGWA2, &
!	XCEMI2, &
!	XCRHO2, &
!	CSTMP2, &
!	CFRA12, &
	FCLEAR, &
	CFRA1X, &
	CFRA2X)
!=================================================================

	USE INCFTC

	IMPLICIT NONE

!-----------------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!	none


!-----------------------------------------------------------------------
!	ARGUMENTS
!-----------------------------------------------------------------------
!	Input parameters:
!	RECORD /RTPHEAD/ HEAD   ! header data
!	RECORD /RTPPROF/ PROF   ! profile data
!	Variables formerly in PROF structure
	TYPE(CLOUD_TYPE), INTENT(IN) :: CLOUD
	INTEGER, INTENT(IN) :: NEMIS  ! length of emissivity array
	REAL, INTENT(IN) :: PSURF   ! surface pressure
	REAL, DIMENSION(:), INTENT(IN) :: PLEVS  ! pressure levels
!	Output parameters:
!	cloud1
	LOGICAL, INTENT(OUT) ::   LBLAC1      ! black cloud?
!	INTEGER, INTENT(OUT) ::   CTYPE1      ! cloud type code number
!	REAL, INTENT(OUT) ::   CFRAC1         ! cloud fraction
!	REAL, INTENT(OUT) ::   CPSIZ1         ! particle size (um)
!	REAL, INTENT(OUT) ::   CPRTO1         ! cloud top pressure (mb)
!	REAL, INTENT(OUT) ::   CPRBO1         ! cloud bottom pressure (mb)
!	REAL, INTENT(OUT) ::   CNGWA1         ! cloud amount (g/m^2)
!	REAL, INTENT(OUT) :: XCEMI1(MXCHAN) ! emissivity
!	REAL, INTENT(OUT) :: XCRHO1(MXCHAN) ! reflectivity
!	REAL, INTENT(OUT) :: CSTMP1         ! cloud/surf temperature (K)
!	cloud2
	LOGICAL, INTENT(OUT) :: LBLAC2      ! black cloud?
!	INTEGER, INTENT(OUT) :: CTYPE2      ! cloud type code number
!	REAL, INTENT(OUT) :: CFRAC2         ! cloud fraction
!	REAL, INTENT(OUT) ::CPSIZ2         ! particle size (um)
!	REAL, INTENT(OUT) ::CPRTO2         ! cloud top pressure (mb)
!	REAL, INTENT(OUT) ::CPRBO2         ! cloud bottom pressure (mb)
!	REAL, INTENT(OUT) ::CNGWA2         ! cloud amount (g/m^2)
!	REAL, INTENT(OUT) ::XCEMI2(MXCHAN) ! emissivity
!	REAL, INTENT(OUT) ::XCRHO2(MXCHAN) ! reflectivity
!	REAL, INTENT(OUT) ::CSTMP2         ! cloud/surf temperature (K)
!	other cloud fractions
!	REAL, INTENT(OUT) ::CFRA12         ! both clouds fraction
	REAL, INTENT(OUT) ::FCLEAR         ! clear fraction
	REAL, INTENT(OUT) ::CFRA1X         ! exclusive cloud1 fraction
	REAL, INTENT(OUT) ::CFRA2X         ! exclusive cloud2 fraction

!-----------------------------------------------------------------------
!	LOCAL VARIABLES
!-----------------------------------------------------------------------
	INTEGER ::     I      ! generic integer
	INTEGER :: ICASE      ! case number for clouds

!-----------------------------------------------------------------------
!	SAVE STATEMENTS
!-----------------------------------------------------------------------
!	none


!***********************************************************************
!***********************************************************************
!      EXECUTABLE CODE begins below
!***********************************************************************
!***********************************************************************

!	Determine which of 5 possible cloud cases is in effect:
!		1: cfrac1=0 and cfrac2=0
!		2: cfrac1>0 and cfrac2=0
!		3: cfrac1=0 and cfrac2>0
!		4: cfrac1>0 and cfrac2>0 and cprtop1 <= cprtop2
!		5: cfrac1>0 and cfrac2>0 and cprtop1 > cprtop2
!
!	CFRAC1=CLOUD%cfrac
!	CFRAC2=CLOUD%cfrac2
!	CPRTO1=CLOUD%cprtop
!	CPRTO2=CLOUD%cprtop2
!
!	Initialize with default values
	FCLEAR=1.0
!	CFRA12=0.0
	CFRA1X=0.0
	CFRA2X=0.0
	LBLAC1=.TRUE.
	LBLAC2=.TRUE.
	ICASE=0
!	
	IF (CLOUD%CFRAC1 .LE. 0.0) THEN
		IF (CLOUD%CFRAC2 .LE. 0.0) THEN
			ICASE=1
		ELSE
			ICASE=3
		ENDIF
	ELSE
		IF (CLOUD%CFRAC2 .LE. 0.0) THEN
			ICASE=2
		ELSE
			IF (CLOUD%CPRTOP1 .LE. CLOUD%CPRTOP2) THEN
				ICASE=4
			ELSE
				ICASE=5
			ENDIF
		ENDIF
	ENDIF
!

!cc
!      print *, 'icase=',ICASE
!cc

!	Assign remaining cloud variables if icase > 1
	IF (ICASE .GT. 1) THEN
!
		IF (ICASE .EQ. 2 .OR. ICASE .EQ. 4) THEN
!			cloud1 RTP fields to cloud1 variables
!			CTYPE1=CLOUD%ctype
!			CPRTO1=CLOUD%cprtop
!			CPRBO1=CLOUD%cprbot
!			CNGWA1=CLOUD%cngwat
!			CPSIZ1=CLOUD%cpsize
			IF (CLOUD%CTYPE1 .LT. 100) THEN
!			WARNING! does not check if cemis is ok
!			CSTMP1=CLOUD%cstemp
				DO I=1, NEMIS
!					XCEMI1(I)=CLOUD%cemis(I)
!					XCRHO1(I)=CLOUD%crho(I)
					IF (CLOUD.CRHO1(I) .LT. 0.0) THEN
						CLOUD.CRHO1(I)=(1 - CLOUD%CEMIS1(I))/PI
					ENDIF
				ENDDO
			ENDIF
!
			IF (ICASE .EQ. 4) THEN
!				cloud2 RTP fields to cloud2 variables
				CTYPE2=CLOUD%ctype2
				CPRTO2=CLOUD%cprtop2
				CPRBO2=CLOUD%cprbot2
				CNGWA2=CLOUD%cngwat2
				CPSIZ2=CLOUD%cpsize2
				IF (CTYPE2 .LT. 100) THEN
!					WARNING! does not check if cemis2 is ok
					CSTMP2=CLOUD%cstemp2
					DO I=1,NEMIS
						XCEMI2(I)=CLOUD%cemis2(I)
						XCRHO2(I)=CLOUD%crho2(I)
						IF (XCRHO2(I) .LT. 0.0) THEN
							XCRHO2(I)=(1 - XCEMI2(I))/PI
						ENDIF
					ENDDO
				ENDIF
			ENDIF
!
		ELSE ! icase=3 or 5
!			cloud2 RTP fields into cloud1 variables
			CFRAC1=CFRAC2
			CFRAC2=0.0
			CTYPE1=CLOUD%ctype2
			CPRTO1=CLOUD%cprtop2
			CPRBO1=CLOUD%cprbot2
			CNGWA1=CLOUD%cngwat2
			CPSIZ1=CLOUD%cpsize2
           IF (CTYPE1 .LT. 100) THEN
!             WARNING! does not check if cemis2 is ok
              CSTMP1=CLOUD%cstemp2
              DO I=1, NEMIS
                 XCEMI1(I)=CLOUD%cemis2(I)
                 XCRHO1(I)=CLOUD%crho2(I)
                 IF (XCRHO1(I) .LT. 0.0) THEN
                  XCRHO1(I)=(1 - XCEMI1(I))/PI
                 ENDIF
              ENDDO
           ENDIF
!
           IF (ICASE .EQ. 5) THEN
!             cloud1 RTP fields to cloud2 variables
              CFRAC2=CLOUD%cfrac
              CTYPE2=CLOUD%ctype
              CPRTO2=CLOUD%cprtop
              CPRBO2=CLOUD%cprbot
              CNGWA2=CLOUD%cngwat
              CPSIZ2=CLOUD%cpsize
              IF (CTYPE2 .LT. 100) THEN
!                WARNING! does not check if cemis is ok
                 CSTMP2=CLOUD%cstemp
                 DO I=1, NEMIS
                  XCEMI2(I)=CLOUD%cemis(I)
                  XCRHO2(I)=CLOUD%crho(I)
                  IF (XCRHO2(I) .LT. 0.0) THEN
                     XCRHO2(I)=(1 - XCEMI2(I))/PI
                  ENDIF
                 ENDDO
              ENDIF
           ENDIF
        ENDIF
!

!         Check cloud1 varibles
        IF (CFRAC1 .GT. 1.0) THEN
           WRITE(IOERR,1010) , 'CFRAC1', '0.0', '1.0'
 1010        FORMAT('Error!', A6, ' outside allowed ',A6,' to ',A6, ' range')
           STOP
        ENDIF
        IF (CPRTO1 .LT. PLEVS(1) .OR. CPRTO1 .GT. PSURF) THEN
           WRITE(IOERR,1010) 'CPRTO1', 'PLEVS1', 'SPRES'
           STOP
        ENDIF
        IF (CTYPE1 .LT. 100) THEN
!          Black cloud
           LBLAC1 = .TRUE.
           CPRBO1=CPRTO1*1.001 ! safe dummy value
        ELSE
!          Slab cloud
           LBLAC1 = .FALSE.
!          Check cprbot, cpsize, & cngwat
           IF (CPRBO1 .LT. PLEVS(1) .OR. CPRBO1 .GT. PSURF) THEN
              WRITE(IOERR,1010)  'CPRBO1', 'PLEVS1', 'SPRES'
              STOP
           ENDIF
           IF (CPRTO1 .GT. CPRBO1) THEN
              WRITE(IOERR,1020)  'CPRTO1', 'CPRBO1'
 1020         FORMAT('Error! ',A6,' > ',A6)
              STOP
           ENDIF
           IF (CPSIZ1 .LT. 0.0 .OR. CPSIZ1 .GT. 1E+3) THEN
              WRITE(IOERR,1010) 'CPSIZ1', '0.0', '1E+3'
              STOP
           ENDIF
           IF (CNGWA1 .LT. 0.0 .OR. CNGWA1 .GT. 1E+6) THEN
              WRITE(IOERR,1010) 'CNGWA1', '0.0', '1E+6'
              STOP
           ENDIF
        ENDIF ! end check of cloud1 variables
!

!         Check cloud2 variables if needed
        IF (CFRAC2 .GT. 0.0) THEN
           IF (CFRAC2 .GT. 1.0) THEN
              WRITE(IOERR,1010)  'CFRAC2', '0.0', '1.0'
              STOP
           ENDIF
           IF (CPRTO1 .LT. PLEVS(1) .OR. CPRTO1 .GT. PSURF) THEN
              WRITE(IOERR,1010)  'CPRTO2', 'PLEVS1', 'SPRES'
              STOP
           ENDIF
           IF (CTYPE2 .LT. 100) THEN
!             Black cloud
              LBLAC2 = .TRUE.
              CPRBO2=CPRTO2*1.001 ! safe dummy value
           ELSE
!             Slab cloud
              LBLAC2 = .FALSE.
!             Check cprbot, cpsize, & cngwat
              IF (CPRBO2 .LT. PLEVS(1) .OR. CPRBO2 .GT. PSURF) THEN
                 WRITE(IOERR,1010)  'CPRBO2', 'PLEVS1', 'SPRES'
                 STOP
              ENDIF
              IF (CPRTO2 .GT. CPRBO2) THEN
                 WRITE(IOERR,1020)  'CPRTO2', 'CPRBO2'
                 STOP
              ENDIF
              IF (CPSIZ2 .LT. 0.0 .OR. CPSIZ2 .GT. 1E+3) THEN
                 WRITE(IOERR,1010)  'CPSIZ2', '0.0', '1E+3'
                 STOP
              ENDIF
              IF (CNGWA2 .LT. 0.0 .OR. CNGWA2 .GT. 1E+6) THEN
                 WRITE(IOERR,1010)  'CNGWA2', '0.0', '1E+6'
                 STOP
              ENDIF
           ENDIF
        ENDIF ! end check of cloud2 variables
!

!         Compute exclusive cloud fractions
        IF (CFRA12 .GT. 0.0) THEN
!          First check two-cloud sanity
           IF (CFRA12 .GT. 1.0) THEN
              WRITE(IOERR,1010)  'CFRA12', '0.0', '1.0'
              STOP
           ENDIF
           IF (LBLAC1) THEN
!             If top cloud is black then must have cfra12=0
              WRITE(IOERR,1055) 
 1055         FORMAT('Error! may not have CFRA12 > 0 with a black top cloud')
              STOP
           ELSE
              IF (CPRBO1 .GT. CPRTO2) THEN
!                Bottom of top cloud is blow top of bottom cloud
                 WRITE(IOERR,1020)  'CPRBO1', 'CPRTO2'
              ENDIF
           ENDIF
!          Now compute exclusive cloud fractions
           IF ((CFRA12 .LE. CFRAC1) .AND. (CFRA12 .LE. CFRAC2)) THEN
              CFRA1X=CFRAC1 - CFRA12
              CFRA2X=CFRAC2 - CFRA12
           ELSE
              WRITE(IOERR,1065) 
 1065         FORMAT('Error!  may not have cfra12 > cfrac1 or cfrac2')
!       write(6,*) 'cfrac1,cfrac2,cfra12=',CFRAC1,CFRAC2,CFRA12
           ENDIF
        ELSE
           CFRA1X=CFRAC1
           CFRA2X=CFRAC2
        ENDIF
!

!         Compute clear fraction
        FCLEAR=1.0 - (CFRA1X + CFRA2X + CFRA12)
        IF (FCLEAR .LT. 0.0) THEN
           IF (FCLEAR .GT. -1E-4) THEN
!             close enough to interpret as zero
              FCLEAR=0.0
           ELSE
              WRITE(IOERR,1070) CFRA1X,CFRA2X,CFRA12,1.0-FCLEAR
 1070         FORMAT('Error! cfra1x=',F7.5, ' + cfra2x=',F7.5,' + cfra12=',F7.5,' =',F7.5,' > 1')
              STOP
           ENDIF
        ENDIF
!
       ENDIF ! end icase > 1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      uncomment for cloud summary 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       WRITE(6,2010) IPROF
! 2010  FORMAT('-------- Cloud Summary for profile',I5' --------')
!       WRITE(6,2011) FCLEAR
! 2011  FORMAT('fclear=',F6.3)
!C
!       IF (CFRAC1 .GT. 0.0) THEN
!        WRITE(6,2021) CTYPE1, CFRAC1, CFRA1X, CPRTO1
! 2021   FORMAT('Cloud1: ctype1=',I5,', cfrac1=',F6.3,', cfra1x=',F6.3,
!     $  ', cprto1=',F8.3)
!        IF (LBLAC1) THEN
!         WRITE(6,2022) CSTMP1
! 2022      FORMAT('Black cloud1: cstmp1=',F9.3)
!        ELSE
!         WRITE(6,2023) CPRBO1, CNGWA1, CPSIZ1
! 2023      FORMAT('Complex cloud1: cprbo1=',F8.3,', cngwa1=',F7.3,
!     $     ', cpsiz1=',F7.2)
!        ENDIF
!C
!        IF (CFRAC2 .GT. 0.0) THEN
!         WRITE(6,2034) CTYPE2, CFRAC2, CFRA2X, CPRTO2
! 2034    FORMAT('Cloud2: ctype2=',I5,', cfrac2=',F6.3,', cfra2x=',F6.3,
!     $   ', cprto2=',F8.3)
!         IF (LBLAC2) THEN
!          WRITE(6,2035) CSTMP2
! 2035       FORMAT('Black cloud2: cstmp2=',F9.3)
!         ELSE
!          WRITE(6,2036) CPRBO2, CNGWA2, CPSIZ2
! 2036       FORMAT('Complex cloud2: cprbo2=',F8.3,', cngwa2=',F7.3,
!     $      ', cpsiz2=',F7.2)
!         ENDIF
!         WRITE(6,2037) CFRA12
! 2037    FORMAT('cfra12=',F6.3)
!        ENDIF ! cfrac2 > 0
!       ENDIF ! cfrac1 > 0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
       RETURN
       END SUBROUTINE GETCLD
