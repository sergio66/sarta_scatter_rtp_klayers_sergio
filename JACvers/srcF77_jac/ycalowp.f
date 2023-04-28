c This version uses amounts at center of layer, not lower boundary
C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALOWP
C
!F77====================================================================


!ROUTINE NAME:
C    CALOWP


!ABSTRACT:
C    Calculate the OPTRAN water (H2O) predictors for a profile.


!CALL PROTOCOL:
C    CALOWP ( LBOT, WAMNT, P, T, SECANG, WAZOP, WAVGOP,
C       WAANG, LOPMIN, LOPMAX, LOPUSE, H2OPRD, LOPLOW, DAOP, DOJAC, H2OJACPRD, DAOPJAC )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   LBOT    bottom layer number         none
C    REAL arr  WAMNT   profile layer water         kiloMoles/cm^2
C    REAL arr  P       layer pressures             atmospheres
C    REAL arr  T       profile temperature         K
C    REAL arr  SECANG  secant of path angle        none
C    REAL arr  WAZOP   OPTRAN l-to-s water grid    kiloMoles/cm^2
C    REAL arr  WAVGOP  OPTRAN average preds        various
C    LOGICAL   DOJAC

!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  WAANG   water amount in layer       kiloMoles/cm^2
C    INTEGER   LOPMIN  min OPTRAN level to use     none
C    INTEGER   LOPMAX  max OPTRAN level to use     none
C    LOG arr   LOPUSE  OPTRAN level needed?        none
C    REAL arr  H2OPRD  OPTRAN predictors           various
C    INTEGER   LOPLOW  low bracketing OPTRAN lev   none
C    REAL arr  DAOP    OPTRAN-to-AIRS interp frac  none
C    REAL arr  H2OJACPRD  OPTRAN jac predictors           various
C    REAL arr  DAOPJAC    OPTRAN-to-AIRS interp frac  none


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    USEFAST


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    March 1998 version of the 100 layer AIRS Fast Transmittance
C    Code by L.L.Strow/S.Hannon.


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    Assumes the user has supplied vaguely realistic profile amounts
C    and temperatures.


!ROUTINE HISTORY:
C    Date         Programmer      Comments
C    -----------  --------------  --------------------------------------
C    27 Feb 1998  Scott Hannon    Created
C    26 Aug 1998  Scott Hannon    Add LBOT to call; loop on LBOT instead
C                                 of MAXLAY


!END====================================================================

C      =================================================================
       SUBROUTINE YCALOWP ( LBOT, WAMNT, P, T, SECANG, WAZOP, WAVGOP,
     $    WAANG, LOPMIN, LOPMAX, LOPUSE, H2OPRD, LOPLOW, DAOP, 
     $    DOJAC, H2OJACPRD, DAOPJAC )
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
       INTEGER   LBOT
       REAL  WAMNT(MAXLAY)
       REAL      P(MAXLAY)
       REAL      T(MAXLAY)
       REAL SECANG(MAXLAY)
       REAL  WAZOP(MXOWLY)
       REAL WAVGOP(NOWAVG,MXOWLY)
       LOGICAL DOJAC
C
C      Output
       REAL  WAANG(MAXLAY)
       INTEGER LOPMIN
       INTEGER LOPMAX
       REAL H2OPRD(  NH2O,MXOWLY)
       LOGICAL LOPUSE(MXOWLY)
       INTEGER LOPLOW(MAXLAY)
       REAL  DAOP(MAXLAY)
       REAL  DAOPJAC(2, MAXLAY)
       REAL H2OJACPRD(2, NH2O,MXOWLY)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      L
       INTEGER     LL
       INTEGER    LOP
       INTEGER   LOPL
       INTEGER   LOPU
       INTEGER     LU
       REAL    WAZ(MAXLAY), WAZ_T(MAXLAY), WAZ_1(MAXLAY)
       REAL WAZSUM, WAZSUM_T, WAZSUM_1
       REAL WPZSUM, WPZSUM_T, WPZSUM_1
       REAL WTZSUM, WTZSUM_T, WTZSUM_1
       REAL     PZ(MAXLAY), PZ_T(MAXLAY), PZ_1(MAXLAY)
       REAL     TZ(MAXLAY), TZ_T(MAXLAY), TZ_1(MAXLAY)
       REAL     DA, DA_T, DA_1
       REAL    POP, POP_T, POP_1
       REAL    TOP, TOP_T, TOP_1
       REAL   PZOP, PZOP_T, PZOP_1
       REAL   TZOP, TZOP_T, TZOP_1
       REAL  ANGOP, ANGOP_T, ANGOP_1
       REAL  WAZOP_T(MXOWLY), WAZOP_1(MXOWLY)
       REAL WAVGOP_T(NOWAVG,MXOWLY), WAVGOP_1(NOWAVG,MXOWLY)
       REAL  WAANG_T(MAXLAY), WAANG_1(MAXLAY)
       REAL JUNK, JUNK_T, JUNK_1
       LOGICAL   LAST


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
C      Initialize amount above sums
       WAZSUM=0.0E+0
       WPZSUM=0.0E+0
       WTZSUM=0.0E+0
C
C      ---------------------------------------
C      Calculate raw predictors for all layers
C      ---------------------------------------
       DO L=1,LBOT
C
C         Layer amount*angle
          WAANG(L)=WAMNT(L)*SECANG(L)
C
C         Center-of-layer-to-space amount*angle
C         Note: do this before updating AZSUM
          WAZ(L)=5.0E-1*WAANG(L) + WAZSUM
C
C         Bottom-of-layer-to-space amount sum
          WAZSUM=WAANG(L) + WAZSUM
C
C         Pressure above sum
          WPZSUM=WAANG(L)*P(L) + WPZSUM
          PZ(L)=WPZSUM/WAZSUM
C
C         Temperature above sum
          WTZSUM=WAANG(L)*T(L) + WTZSUM
          TZ(L)=WTZSUM/WAZSUM
C
c          write(6,'(A,I3,X,E11.4,X,E11.4)') 'calowp:WAZ(L),WAZSUM : ',L,WAZ(L),WAZSUM
           IF (DOJAC) THEN
             WAANG_T(L) = 0
             WAANG_1(L) = SECANG(L)
             WAZ_T(L)   = 0
             WAZ_1(L)   = 0.5*SECANG(L)
             WAZSUM_T   = 0
             WAZSUM_1   = SECANG(L)
             WPZSUM_T   = 0
             WPZSUM_1   = SECANG(L)*P(L)
             WTZSUM_T   = WAANG(L)
             WTZSUM_1   = SECANG(L)*T(L)
             PZ_T       = (WAZSUM*WPZSUM_T -  WPZSUM*WAZSUM_T)/WAZSUM/WAZSUM
             PZ_1       = (WAZSUM*WPZSUM_1 -  WPZSUM*WAZSUM_1)/WAZSUM/WAZSUM
             TZ_T       = (WAZSUM*WTZSUM_T -  WTZSUM*WAZSUM_T)/WAZSUM/WAZSUM
             TZ_1       = (WAZSUM*WTZSUM_1 -  WTZSUM*WAZSUM_1)/WAZSUM/WAZSUM
           END IF
       ENDDO
C
       if (DEBUG) print*,'calowp: completed raw predictors'
C      --------------------------------------------------
C      Find the max OPTRAN level that is less than WAZ(1)
C      --------------------------------------------------
       LOPMIN=1
       if (DEBUG) WRITE(6,'(A, E11.4)') 'calowp: WAZ(1) = ',WAZ(1)
 30    IF (WAZOP(LOPMIN+1) .LT. WAZ(1)) THEN
          LOPMIN=LOPMIN + 1
          GOTO 30
       ENDIF
C
C      Initialize the upper and lower (pressure) layer index
       LL=1
       LU=2
       LAST=.FALSE.
C
       if (DEBUG) print*,'calowp: completed find max optran level'
C      ----------------------------------------
C      Loop over the OPTRAN layers (while loop)
C      ----------------------------------------
       LOP=LOPMIN
 10    IF (LOP .LE. MXOWLY) THEN
C
C         --------------------------------------------------------
C         Find the two pressure layers closest to the OPTRAN layer
C         --------------------------------------------------------
 20       IF (WAZ(LU) .LT. WAZOP(LOP)) THEN
             IF (LU .LT. LBOT) THEN
                LL=LU
                LU=LU + 1
                GOTO 20
             ELSE
                LAST=.TRUE.
             ENDIF
          ENDIF
C
C         Compute the interpolation fractor
          DA=(WAZOP(LOP) - WAZ(LL))/(WAZ(LU) - WAZ(LL))
C
C         Do the interpolation
          POP=( DA*(  P(LU) -  P(LL) ) +  P(LL) )/WAVGOP(1,LOP)
          TOP=( DA*(  T(LU) -  T(LL) ) +  T(LL) )/WAVGOP(2,LOP)
          PZOP=( DA*( PZ(LU) - PZ(LL) ) + PZ(LL) )/WAVGOP(3,LOP)
          TZOP=( DA*( TZ(LU) - TZ(LL) ) + TZ(LL) )/WAVGOP(4,LOP)
          ANGOP=DA*( SECANG(LU) - SECANG(LL) ) + SECANG(LL)
C
C         Assign the predictors
          H2OPRD(1,LOP)=1.0E+0
          H2OPRD(2,LOP)=POP
          H2OPRD(3,LOP)=TOP
          H2OPRD(4,LOP)=SQRT( POP )
          H2OPRD(5,LOP)=TOP**2
          H2OPRD(6,LOP)=POP*TOP
          H2OPRD(7,LOP)=ANGOP
          H2OPRD(8,LOP)=PZOP
          H2OPRD(9,LOP)=TZOP

          IF (DOJAC) THEN
            DA_T = 1/(WAZ(LU) - WAZ(LL)) * WAZ_T(LL)
            DA_1 = 1/(WAZ(LU) - WAZ(LL)) * WAZ_1(LL)

            JUNK = P(LU) -  P(LL)
            POP_T = DA_T*junk/WAVGOP(1,LOP) - (DA*junk + P(LL))/WAVGOP(1,LOP)/WAVGOP(1,LOP)*WAVGOP_T(1,LOP)
            POP_1 = DA_1*junk/WAVGOP(1,LOP) - (DA*junk + P(LL))/WAVGOP(1,LOP)/WAVGOP(1,LOP)*WAVGOP_1(1,LOP)

            JUNK = T(LU) -  T(LL)
            TOP_T = DA_T*junk/WAVGOP(2,LOP) - (DA*junk + T(LL))/WAVGOP(2,LOP)/WAVGOP(2,LOP)*WAVGOP_T(2,LOP)
            TOP_1 = DA_1*junk/WAVGOP(2,LOP) - (DA*junk + T(LL))/WAVGOP(2,LOP)/WAVGOP(2,LOP)*WAVGOP_1(2,LOP)

            !!! assume dPZ/dT = dPZ/dQ = 0
            JUNK = PZ(LU) -  PZ(LL)
            PZOP_T = DA_T*junk/WAVGOP(3,LOP) - (DA*junk + PZ(LL))/WAVGOP(3,LOP)/WAVGOP(3,LOP)*WAVGOP_T(3,LOP)
            PZOP_1 = DA_1*junk/WAVGOP(3,LOP) - (DA*junk + PZ(LL))/WAVGOP(3,LOP)/WAVGOP(3,LOP)*WAVGOP_1(3,LOP)

            !!! assume dTZ/dT = dTZ/dQ = 0
            JUNK = TZ(LU) -  TZ(LL)
            TZOP_T = DA_T*junk/WAVGOP(4,LOP) - (DA*junk + TZ(LL))/WAVGOP(4,LOP)/WAVGOP(4,LOP)*WAVGOP_T(4,LOP)
            TZOP_1 = DA_1*junk/WAVGOP(4,LOP) - (DA*junk + TZ(LL))/WAVGOP(4,LOP)/WAVGOP(4,LOP)*WAVGOP_1(4,LOP)

            ANGOP_T = DA_T*( SECANG(LU) - SECANG(LL) )
            ANGOP_1 = DA_1*( SECANG(LU) - SECANG(LL) )

            !! TJAC 
            H2OJACPRD(1,1,LOP)=0
            H2OJACPRD(1,2,LOP)=POP_T
            H2OJACPRD(1,3,LOP)=TOP_T
            H2OJACPRD(1,4,LOP)=0.5/SQRT(POP)*POP_T
            H2OJACPRD(1,5,LOP)=2*TOP*TOP_T
            H2OJACPRD(1,6,LOP)=POP*TOP_T + TOP*POP_T
            H2OJACPRD(1,7,LOP)=ANGOP_T
            H2OJACPRD(1,8,LOP)=PZOP_T
            H2OJACPRD(1,9,LOP)=TZOP_T

            !! WJAC 
            H2OJACPRD(2,1,LOP)=0
            H2OJACPRD(2,2,LOP)=POP_1
            H2OJACPRD(2,3,LOP)=TOP_1
            H2OJACPRD(2,4,LOP)=0.5/SQRT(POP)*POP_1
            H2OJACPRD(2,5,LOP)=2*TOP*TOP_1
            H2OJACPRD(2,6,LOP)=POP*TOP_1 + TOP*POP_1
            H2OJACPRD(2,7,LOP)=ANGOP_1
            H2OJACPRD(2,8,LOP)=PZOP_1
            H2OJACPRD(2,9,LOP)=TZOP_1
          END IF

C
C         Update LOP and loop
          IF (LAST .EQV. .TRUE.) THEN
             LOPMAX=LOP
C            Set LOP > MXOWLY to exit loop over LOP
             LOP=MXOWLY + 1
          ELSE
             LOP=LOP + 1
          ENDIF
          GOTO 10
C
       ENDIF
C      End while loop over LOP
C
C      -----------------
C      Initialize LOPUSE
C      -----------------
       DO LOP=1,MXOWLY
          LOPUSE(LOP)=.FALSE.
       ENDDO
C
C      ---------------------------------------
C      Determine what OPTRAN layers are needed
C      ---------------------------------------
C      Initialize LOPL and LOPU
       LOPL=LOPMIN
       LOPU=LOPMIN + 1
C
C      Loop over the AIRS pressure layers
       DO L=1,LBOT
C         Find the two OPTRAN levels that bracket the AIRS layer
 40       IF (WAZOP(LOPU) .LT. WAZ(L) .AND. LOPU .LT. LOPMAX) THEN
             LOPL=LOPU
             LOPU=LOPU + 1
             GOTO 40
          ENDIF
C
          LOPUSE(LOPL)=.TRUE.
          LOPUSE(LOPU)=.TRUE.
C         Assign the lower OPTRAN level
          LOPLOW(L)=LOPL
C         Assign the interpolation fraction
          DAOP(L)=(WAZ(L) - WAZOP(LOPL))/(WAZOP(LOPU) - WAZOP(LOPL))

          IF (DOJAC) THEN
            DAOPJAC(1,L) = WAZ_T(L)/(WAZOP(LOPU) - WAZOP(LOPL))
            DAOPJAC(2,L) = WAZ_1(L)/(WAZOP(LOPU) - WAZOP(LOPL))
          END IF

       ENDDO
C
       RETURN
       END
