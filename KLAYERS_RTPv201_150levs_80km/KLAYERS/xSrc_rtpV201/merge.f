C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              MERGE
C
!F77====================================================================

!ROUTINE NAME: MERGE

!ABSTRACT:
C    Merge the AFGL profile with the user profile.


!CALL PROTOCOL:
C    MERGE(LZ, NGASES, PMIN, NINX, LISTG, GORDER,
C          NIN,    ZIN, PIN, LNPIN, NUSER, GASID,  TIN, MRIN,
C          NLEVAF, ZAF, PAF, LNPAF, NAFGL, IDAFGL, TAF, MIXAF, CO2MLT)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      CO2MLT  mult for AFGL CO2 mix ratio none
C    INT arr   IDAFGL  AFGL gas IDs (HITRAN)       none
C    INT arr   LISTG   list of desired gases       none
C    REAL arr  LNPAF   log of AFGL pres levels     log(atm)
C    LOGICAL   LZ      user profile altitudes      none
C    REAL arr  MIXAF   AFGL gas mixing ratios      ppmv
C    INTEGER   NAFGL   number of AFGL prof gases   none
C    INTEGER   NLEVAF  number of AFGL prof levels  none
C    INTEGER   NIN     number of user prof levels  none
C    INTEGER   NUSER   number of user prof gases   none
C    REAL arr  PAF     AFGL prof pressure levels   atm
C    REAL arr  PMIN    minimum AIRS pressure lev   atm
C    REAL arr  TAF     AFGL prof temperatures      K
C    REAL arr  ZAF     AFGL altitudes              meters


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INT arr   GORDER  GASID order rel to LISTG    none
C    INTEGER   NGASES  number of gases             none
C    INTEGER   NINX    # extended user prof levs   none


!INPUT/OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INT arr   GASID   gas IDs (HITRAN)            none
C    REAL arr  LNPIN   log of user pres levels     log(atm)
C    REAL arr  MRIN    gas mixing ratios           ppmv
C    REAL arr  PIN     user prof pressure levels   atm
C    REAL arr  TIN     temperatures                K
C    REAL arr  ZIN     user altitudes              meters


!RETURN VALUES: none

!PARENT(S): KLAYERS

!ROUTINES CALLED: none

!FILES ACCESSED:
C    Input unit IOUNIT

!COMMON BLOCKS: none

!DESCRIPTION:
C    If necessary, extend the user profile with the AFGL model.

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C Apr 21 1995 Scott Hannon/UMBC created
C Jun 23 1995 Scott Hannon      Correct some comments
C Oct  5 1998 Scott Hannon      Add MRAFGL to test AFGL interpolation
C 21 Feb 2001 Scott Hannon      Warning messages now use unit IOINFO
C 24 Jul 2001 Scott Hannon      Add CO2MLT and XBOTH & XONLY
C  8 Feb 2002 Scott Hannon      Add IF THEN ELSE to handle NIN=0
C 29 Apr 2002 Scott Hannon      Add LZ, ZIN, & ZAF to call to allow
C                                  extending profile with AFGL alt.
C 16 Feb 2007 Scott Hannon      Correct IDAFGL dim to MXGAS (was MXIN)


!END====================================================================

C      =================================================================
       SUBROUTINE MERGE(LZ, NGASES, PMIN, NINX, LISTG, GORDER,
     $    NIN,    ZIN, PIN, LNPIN, NUSER, GASID,  TIN, MRIN,
     $    NLEVAF, ZAF, PAF, LNPAF, NAFGL, IDAFGL, TAF, MIXAF, CO2MLT)
C      =================================================================

C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE

C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
      include 'incLAY.f'

C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input parameters:
       LOGICAL     LZ           ! user profile uses altitudes
       INTEGER    NIN           ! number of user profile levels
       INTEGER  NUSER           ! number of user profile gases
       INTEGER  NAFGL           ! number of AFGL gases
       INTEGER NLEVAF           ! number of AFGL levels
c       INTEGER IDAFGL(  MXIN)   ! AFGL gas IDs
       INTEGER IDAFGL( MXGAS)   ! AFGL gas IDs
       INTEGER  LISTG( MXGAS)   ! user gas IDs
       REAL   PMIN                 ! min pressure  
       REAL    PIN(  MXIN)         ! user profile pressures
       REAL  LNPIN(  MXIN)         ! log of user prof press
       REAL    ZAF(  MXIN)         ! AFGL altitudes *kilometers*
       REAL    PAF(  MXIN)         ! AFGL pressures
       REAL  LNPAF(  MXIN)         ! log of AFGL press
       REAL    TAF(  MXIN)         ! AFGL temperatures
       REAL  MIXAF(  MXIN, MXGAS)  ! AFGL mixing ratios
       REAL CO2MLT                 ! AFGL CO2 mixing ratio multiplier
C
C      Output parameters:
       INTEGER NGASES          ! Total number of gases
       INTEGER   NINX          ! number of extended profile levels
       INTEGER GORDER( MXGAS)  ! order of gases in gasid
C
C      Input/Output parameters
       INTEGER  GASID( MXGAS)      ! user/extended profile gas IDs
       REAL    ZIN(  MXIN)         ! user/extended profile alts *meters*
       REAL    TIN(  MXIN)         ! user/extended profile temperatures
       REAL   MRIN(  MXIN, MXGAS)  ! user/extended profile mixing ratios


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      I                 ! generic looping
       INTEGER  IBOTH( MXGAS,     2)  ! ID & index of "both" gases
       INTEGER   ILEV                 ! level looping
       INTEGER   IOLD
       INTEGER  IONLY( MXGAS)         ! IDs of AFGL-only gases
       INTEGER      J                 ! generic looping
       INTEGER  NBOTH                 ! # "both" (user and AFGL) gases
       INTEGER NFIRST                 ! first AFGL level for extending
       INTEGER  NLAST                 ! last AFGL level for extending
       INTEGER  NONLY                 ! # gases AFGL-only
C
       REAL      A          ! linear interpolation slope
       REAL      B          ! linear interpolation intercept
       REAL MRAFGL
       REAL  XAFGL
       REAL  XBOTH( MXGAS)  ! mixing ratio mult for "both" gases
       REAL  XONLY( MXGAS)  ! mixing ratio mult for "only" gases
       REAL  XUSER
       REAL ZOFFST          ! altitude offset for extended levels


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
C      -----------------------------------------------------
C      Look for gases in both user & AFGL profs or AFGL only
C      -----------------------------------------------------
       NBOTH=0
       NONLY=0
       DO J=1,NAFGL
          IOLD=NBOTH
          DO I=1,NUSER
C
C            Look for matches in both profiles
             IF (GASID(I) .EQ. IDAFGL(J)) THEN
                NBOTH=NBOTH+1
                IBOTH(NBOTH,1)=I
                IBOTH(NBOTH,2)=J
C               Set "both" AFGL mixing ratio multiplier
                IF (IDAFGL(J) .EQ. 2) THEN
                   XBOTH(NBOTH)=CO2MLT
                ELSE
                   XBOTH(NBOTH)=1.0
                ENDIF
             ENDIF
          ENDDO
C
C         If no match in both, then AFGL only
          IF (IOLD .EQ. NBOTH) THEN
             NONLY=NONLY+1
             IONLY(NONLY)=J
C               Set "only" AFGL mixing ratio multiplier
             IF (IDAFGL(J) .EQ. 2) THEN
                XONLY(NONLY)=CO2MLT
             ELSE
                XONLY(NONLY)=1.0
             ENDIF
          ENDIF
       ENDDO
C
C      ------------------
C      Extend the profile
C      ------------------
       NINX=NIN
       XAFGL=PAF(NLEVAF)
       IF (NIN .EQ. 0) THEN
          XUSER=1.1*PAF(1)  ! arbitrary value > max(PAF)
       ELSE
          XUSER=PIN(NIN)
       ENDIF
C
       IF (XUSER .GT. XAFGL .AND. XUSER .GT. PMIN) THEN
C
C         Find first point in AFGL prof to add
          NFIRST=1
 10       IF (PAF(NFIRST) .GT. 0.999*XUSER) THEN
             NFIRST=NFIRST+1
             GOTO 10
          ENDIF
C
C         Find last point in AFGL prof to add
          IF (PMIN .GT. XAFGL) THEN
             NLAST=NFIRST
 20          IF (PAF(NLAST) .GT. PMIN) THEN
                NLAST=NLAST+1
                GOTO 20
             ENDIF
          ELSE
             NLAST=NLEVAF
          ENDIF
C
C         Calculate offset for extended altitudes
          ZOFFST=0.0
          IF (LZ) THEN
             IF (NFIRST .EQ. 1) THEN
                I=NFIRST
             ELSE
                I=NFIRST - 1
             ENDIF
             A=1000.0*( ZAF(I+1) - ZAF(I) )/( LNPAF(I+1) - LNPAF(I) )
             B=1000.0*ZAF(I) - LNPAF(I)*A
             ZOFFST=ZIN(NIN) - ( LNPIN(NIN)*A + B )
          ENDIF
C
C         Tack on the temp & mix ratio extensions
          DO I=NFIRST,NLAST
C
             NINX=NIN+I+1-NFIRST
             ZIN(NINX)=ZAF(I)*1000.0 + ZOFFST
             PIN(NINX)=PAF(I)
             LNPIN(NINX)=LNPAF(I)
             TIN(NINX)=TAF(I)
C
C            Mixing ratios for extended levels
             DO J=1,NBOTH
                MRIN(NINX,IBOTH(J,1))=MIXAF(I,IBOTH(J,2))*XBOTH(J)
             ENDDO
             DO J=1,NONLY
                MRIN(NINX,NUSER+J)=MIXAF(I,IONLY(J))*XONLY(J)
             ENDDO
          ENDDO
C
       ENDIF
C
C      --------------------------------------------------
C      Interpolate AFGL-only mixing ratios to user levels
C      --------------------------------------------------
       IF (NONLY .GT. 0) THEN
C         Loop over the PIN levels
          ILEV=2
          DO I=1,NIN
 30          IF (PAF(ILEV) .GT. PIN(I) .AND. ILEV .LT. NLEVAF) THEN
                ILEV=ILEV+1
                GOTO 30
             ENDIF
             DO J=1,NONLY
                CALL INTERP(LNPAF(ILEV-1),MIXAF(ILEV-1,IONLY(J)),
     $                      LNPAF(ILEV),  MIXAF(ILEV,IONLY(J)),A,B)
C
                MRAFGL=( A*LNPIN(I) + B )*XONLY(J)
C               Test that the interpolated value is non-negative
                IF (MRAFGL .LT. 0) THEN
                   WRITE(IOINFO,1010) IDAFGL( IONLY(J) ), I
 1010              FORMAT(/,'Warning! Interpolated AFGL mixing ratio',
     $             ' in MERGE was negative for',/,'gas ID',I3,
     $             ', user prof level',I4,
     $             '.  Reset mixing ratio to zero.')
                   MRAFGL=0.0E+0
                ENDIF
                MRIN(I,NUSER+J)=MRAFGL
             ENDDO
          ENDDO
       ENDIF
C
       NGASES=NUSER + NONLY
       DO I=1,NONLY
          GASID(NUSER+I)=IDAFGL(IONLY(I))
       ENDDO
C
C      ---------------------------------------------------------
C      Determine array order of gases in GASID relative to LISTG
C      ---------------------------------------------------------
       DO I=1,NGASES
          DO J=1,NGASES
             IF (LISTG(I) .EQ. GASID(J)) GORDER(I)=J
          ENDDO
       ENDDO
C
       RETURN
       END
