!=======================================================================
!
!              University of Maryland Baltimore County [UMBC]
!
!              AIRS
!
!              BKPREP
!
!F77====================================================================


!ROUTINE NAME: BKPREP


!ABSTRACT:
!    Prepare black cloud parameters


!CALL PROTOCOL:
!    BKPREP( IPROF, ICLOUD, CTYPE, CFRAC, CPRTOP, LBOT, PSURF,
!       PLEV, PLAY, TEMP, LTOPC, TTOPC, TLAC, FLAC )


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INTEGER   IPROF   profile number              none
!    INTEGER   ICLOUD  cloud number                none
!    INTEGER   CTYPE   cloud type code number      none
!    REAL      CFRAC   cloud fraction              none
!    REAL      CPRTOP  cloud top pressure          mb
!    INTEGER   LBOT    bottom layer                none
!    REAL      PSURF   surface pressure            mb
!    REAL arr  PLEV    pressure level grid         mb
!    REAL arr  PLAY    layer mean pressure         mb
!    REAL arr  TEMP    layer mean temperature      Kelvin


!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    INTEGER   LTOPC   layer containing cloud top  none
!    REAL      TTOPC   cloud top temperature       Kelvin
!    REAL      TLAC    T in frac layer above cloud Kelvin
!    REAL      FLAC    fractional layer above cld  none


!RETURN VALUES: none


!PARENT(S): SARTA


!ROUTINES CALLED: none


!FILES ACCESSED:
!    none


!COMMON BLOCKS: none


!DESCRIPTION:
!    Prepare black cloud parameters required for RTA calculations.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
!    none


!ROUTINE HISTORY:
!    Date     Programmer        Comments
!------------ ----------------- ----------------------------------------
! 21 Feb 2007 Scott Hannon      created
! 10 Sep 2007 Scott Hannon      Set CRHO (was unassigned) and maybe
!                                  reset CEMIS depending on CTYPE
! 14 Nov 2007 Scott Hannon      Added CTYPE=2 (emis=1,rho=(1-cemis)/pi)
! 25 Nov 2008 Scott Hannon      Remove CEMIS & CRHO; ctype not used

!END====================================================================

!      =================================================================
       SUBROUTINE BKPREP(ICLOUD, CTYPE, CFRAC, CPRTOP, &
         LBOT, PSURF, PLEV, PLAY, TEMP, LTOPC, TTOPC, TLAC, FLAC)
!      =================================================================


       USE INCFTC
       IMPLICIT NONE


!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none


!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
!      Input parameters:
       INTEGER  IPROF        ! profile number
       INTEGER ICLOUD        ! cloud number
       INTEGER  CTYPE        ! cloud type code number
       REAL  CFRAC           ! cloud fraction
       REAL CPRTOP           ! cloud top pressure (mb)
       INTEGER   LBOT        ! bottom layer
       REAL  PSURF           ! surface pressure (mb)
       REAL   PLEV(MAXLAY+1) ! pressure level grid (mb)
       REAL   PLAY(MAXLAY)   ! layer mean pressure (mb)
       REAL   TEMP(MAXLAY)   ! layer mean temperature (K)

!      Output parameters:
       INTEGER  LTOPC        ! layer containing cloud top
       REAL  TTOPC           ! cloud top temperature (K)
       REAL  TLAC            ! temp of fractional layer above cloud (K)
       REAL  FLAC            ! fractional layer above cloud (K)



!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
       INTEGER  L      ! looping
       REAL RJUNK1     ! generic
       REAL RJUNK2     ! generic


!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none


!***********************************************************************
!***********************************************************************
!      EXECUTABLE CODE begins below
!***********************************************************************
!***********************************************************************

!      Check cloud top pressure is OK
       IF (CPRTOP .LT. PLEV(1)) THEN
          WRITE(IOERR,1010) ICLOUD, CPRTOP, PLEV(1)
 1010     FORMAT('Error!  cprto',I1,'=',1PE12.5,' < Plev(1) =',1PE12.5)
          STOP
       ENDIF
       IF (CPRTOP .GT. PSURF) THEN
          WRITE(IOERR,1012) ICLOUD, CPRTOP, PSURF
 1012     FORMAT('Error!: cprto',I1,'=',1PE12.5, ' > Psurf =',1PE12.5)
          STOP
       ENDIF

!      Determine cloud top layer
       DO L=1,LBOT
          IF (PLEV(L) .LE. CPRTOP) LTOPC=L
       ENDDO

!      Calc fraction of layer at top of cloud (clear of this cloud)
       FLAC=(CPRTOP - PLEV(LTOPC))/(PLEV(LTOPC+1) - PLEV(LTOPC))

!      Calc fractional layer mean P (rjunk1)
       RJUNK2= CPRTOP/PLEV(LTOPC)
       IF (ABS(1.0 - RJUNK2) .LT. 1E-3) THEN
!         Use linear approx to avoid divide by zero
          RJUNK1=CPRTOP + 0.5*(CPRTOP - PLEV(LTOPC))
       ELSE
          RJUNK1=(CPRTOP - PLEV(LTOPC))/LOG( RJUNK2 )
       ENDIF

!      Linearly interpolate layer mean T to frac layer mean P
       RJUNK2=( TEMP(LTOPC) - TEMP(LTOPC-1) )/ &
         LOG( PLAY(LTOPC)/PLAY(LTOPC-1) )        ! slope
       TLAC=RJUNK2*LOG( RJUNK1/PLAY(LTOPC-1) ) + TEMP(LTOPC-1)

!      Linearly interpolate layer mean T to cloud top P
       L=LTOPC
       RJUNK1=(PLEV(L+1) - PLEV(L))/LOG( PLEV(L+1)/PLEV(L) )
       IF (RJUNK1 .GT. CPRTOP .OR. L .EQ. LBOT) THEN
          RJUNK2=(PLEV(L) - PLEV(L-1))/LOG( PLEV(L)/PLEV(L-1) )
          TTOPC=TEMP(L) + LOG(CPRTOP/RJUNK1)* &
            (TEMP(L-1) - TEMP(L))/LOG( RJUNK2/RJUNK1 )
       ELSE
          RJUNK2=(PLEV(L+2) - PLEV(L+1))/LOG( PLEV(L+2)/PLEV(L+1) )
          TTOPC=TEMP(L) + LOG(CPRTOP/RJUNK1)* &
            (TEMP(L+1) - TEMP(L))/LOG( RJUNK2/RJUNK1 )
       ENDIF
!
       RETURN
       END
