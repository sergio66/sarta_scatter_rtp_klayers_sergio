C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CALOKW
C
!F77====================================================================


!ROUTINE NAME:
C    CALOKW


!ABSTRACT:
C    Calculate the OPTRAN derived water pressure layer effective
C    optical depth for a single channel.


!CALL PROTOCOL:
C    YCALOKW ( LBOT, ICHAN, LOPMIN, LOPMAX, LOPLOW, LOPUSE,
C       H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW, 
C       DOJAC, H2OJACPRD,KW_T,KW_1,SECANG )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   LBOT    bottom pres layer number    none
C    INTEGER   ICHAN   OPTRAN water channel index  none
C    INTEGER   LOPMIN  min OPTRAN level to use     none
C    INTEGER   LOPMAX  max OPTRAN level to use     none
C    INTEGER   LOPLOW  low OPTRAN bracketing lev   none
C    LOG arr   LOPUSE  Need this OPTRAN level?     none
C    REAL arr  SECANG  secant of path angle        none
C    REAL arr  H2OPRD  OPTRAN water predictors     various
C    REAL arr  COFH2O  OPTRAN H2O fast trans coef  various
C    REAL arr  WAOP    OPTRAN layer water amounts  kiloMoles/cm^2
C    REAL arr  DAOP    OPTRAN-to-AIRS interp fact  none
C    REAL arr  WAANG   AIRS layer water amounts    kiloMoles/cm^2
c
C    LOGICAL   DOJAC
C    REAL matr H2OJACPRD  predictor derivatives
C    REAL arr  KW_T    dK/dT
C    REAL arr  KW_1    dK/dG1

!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  KW      AIRS H2O layer eff op dep   none


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    CALT1


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
C
C    The OPTRAN predictors and fast transmittance coefficients are
C    used to calculate the water effective absorption coefficient on
C    on the OPTRAN level grid.  Only the OPTRAN levels actually needed
C    are calculated.
C    Note: the COFH2O*H2OPRD result must be divided by WAOP, a scaling
C    factor which was originally applied during the fast transmittance
C    coefficient regression.
C    The OPTRAN absorption coefficients are then interpolated onto the
C    100 AIRS layers and multiplied by the AIRS layer water amount (to
C    convert absorption coefficient into optical depth).


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date         Programmer      Comments
C    -----------  --------------  --------------------------------------
C    27 Feb 1998  Scott Hannon    Created
C    26 Aug 1998  Scott Hannon    Add LBOT to call; loop on LBOT instead
C                                 of MAXLAY


!END====================================================================

C      =================================================================
       SUBROUTINE YCALOKW ( LBOT, ICHAN, LOPMIN, LOPMAX, LOPLOW, LOPUSE,
     $    H2OPRD, COFH2O, WAOP, DAOP, WAANG, KW, 
     $    DOJAC, SECANG, H2OJACPRD, DAOPJAC, KW_T, KW_1 )
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
       INTEGER  ICHAN
       INTEGER LOPMIN
       INTEGER LOPMAX
       REAL SECANG(MAXLAY)
       INTEGER LOPLOW(MAXLAY)
       LOGICAL LOPUSE(MXOWLY)
       REAL   DAOP(MAXLAY)
       REAL   WAANG(MAXLAY)
       REAL   WAOP(MXOWLY)
       REAL  H2OPRD(  NH2O,MXOWLY)
       REAL  COFH2O(  NH2O,MXOWLY,MXCHNW)
       LOGICAL DOJAC
       REAL   H2OJACPRD(OPTRANJAC,NH2O,MXOWLY)    !! OPTRAN, used by ycalt1_od, ycalt3_od
       REAL  DAOPJAC(OPTRANJAC, MAXLAY)           !! OPTRAN, used by ycalt1_od, ycalt3_od

C      Output
       REAL  KW(MAXLAY)
       REAL  KW_T(MAXLAY)
       REAL  KW_1(MAXLAY)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      L
       INTEGER    LOP
       REAL   KWOP(MXOWLY)
       REAL   KWOP_T(MXOWLY)
       REAL   KWOP_1(MXOWLY)


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
C      ---------------------------------
C      Loop over the OPTRAN water levels
C      ---------------------------------
C$$$       IF (DEBUG)  write(6,'(A,I4, X, I4)') 'calokw: LOPMIN, LOPMAX ', 
C$$$     $      LOPMIN, LOPMAX
C      Only do calc for OPTRAN levels that are needed
       DO LOP=LOPMIN,LOPMAX
          IF (LOPUSE(LOP)) THEN
             KWOP(LOP)=
     $          COFH2O(1,LOP,ICHAN)*H2OPRD(1,LOP) +
     $          COFH2O(2,LOP,ICHAN)*H2OPRD(2,LOP) +
     $          COFH2O(3,LOP,ICHAN)*H2OPRD(3,LOP) +
     $          COFH2O(4,LOP,ICHAN)*H2OPRD(4,LOP) +
     $          COFH2O(5,LOP,ICHAN)*H2OPRD(5,LOP) +
     $          COFH2O(6,LOP,ICHAN)*H2OPRD(6,LOP) +
     $          COFH2O(7,LOP,ICHAN)*H2OPRD(7,LOP) +
     $          COFH2O(8,LOP,ICHAN)*H2OPRD(8,LOP) +
     $          COFH2O(9,LOP,ICHAN)*H2OPRD(9,LOP)
C            Remove WAOP scaling factor
             KWOP(LOP)=KWOP(LOP)/WAOP(LOP)
C            Check for negative value
             IF (KWOP(LOP) .LT. 0.0E+0) KWOP(LOP)=0.0E+0
          ENDIF
       ENDDO
C
C      -------------------------
C      Loop over the AIRS layers
C      -------------------------
       DO L=1,LBOT
C
CCC       catch bug: KWOP(0)
          IF (LOPLOW(L) .GT. 0.0E+0) THEN
          
C         Interpolate abs coef and convert to optical depth
c    ycalowp.f shows WAANG(L)=WAMNT(L)*SECANG(L)
          KW(L)=( DAOP(L)*( KWOP(LOPLOW(L) + 1) -
     $       KWOP(LOPLOW(L)) ) + KWOP(LOPLOW(L)) )*WAANG(L)
          IF (KW(L) .LT. 0.0E+0) KW(L)=0.0E+0
C
         ENDIF
       ENDDO
C
c************************************************************************
       IF (DOJAC) THEN
cbaba TJAC
       DO LOP=LOPMIN,LOPMAX
          IF (LOPUSE(LOP)) THEN
             KWOP_T(LOP)=
     $          COFH2O(1,LOP,ICHAN)*H2OJACPRD(1,1,LOP) +
     $          COFH2O(2,LOP,ICHAN)*H2OJACPRD(1,2,LOP) +
     $          COFH2O(3,LOP,ICHAN)*H2OJACPRD(1,3,LOP) +
     $          COFH2O(4,LOP,ICHAN)*H2OJACPRD(1,4,LOP) +
     $          COFH2O(5,LOP,ICHAN)*H2OJACPRD(1,5,LOP) +
     $          COFH2O(6,LOP,ICHAN)*H2OJACPRD(1,6,LOP) +
     $          COFH2O(7,LOP,ICHAN)*H2OJACPRD(1,7,LOP) +
     $          COFH2O(8,LOP,ICHAN)*H2OJACPRD(1,8,LOP) +
     $          COFH2O(9,LOP,ICHAN)*H2OJACPRD(1,9,LOP)
C            Remove WAOP scaling factor
C            WAOP depends on WAZOP(L)-WAZOP(L-1), and WAZOP are from a file so constant!
             KWOP_T(LOP)=KWOP_T(LOP)/WAOP(LOP)
C            Check for negative value
c             IF (KWOP_T(LOP) .LT. 0.0E+0) KWOP_T(LOP)=0.0E+0
          ENDIF
       ENDDO
C
C      -------------------------
C      Loop over the AIRS layers
C      -------------------------
       DO L=1,LBOT
C
CCC       catch bug: KWOP_T(0)
          IF (LOPLOW(L) .GT. 0.0E+0) THEN
          
C         Interpolate abs coef and convert to optical depth
c    ycalowp.f shows WAANG(L)=WAMNT(L)*SECANG(L)
            KW_T(L)=( DAOP(L)*( KWOP_T(LOPLOW(L) + 1) -
     $       KWOP_T(LOPLOW(L)) ) + KWOP_T(LOPLOW(L)) )*WAANG(L)

c    ycalowp.f shows DAOP depends on T and WV
            KW_T(L)=KW_T(L) + ( DAOPJAC(1,L)*( KWOP(LOPLOW(L) + 1) -
     $       KWOP(LOPLOW(L)) ))*WAANG(L)

c            IF (KW_T(L) .LT. 0.0E+0) KW_T(L)=0.0E+0
C
         ENDIF
       ENDDO
cbaba TJAC

cbaba QJAC
       DO LOP=LOPMIN,LOPMAX
          IF (LOPUSE(LOP)) THEN
             KWOP_1(LOP)=
     $          COFH2O(1,LOP,ICHAN)*H2OJACPRD(2,1,LOP) +
     $          COFH2O(2,LOP,ICHAN)*H2OJACPRD(2,2,LOP) +
     $          COFH2O(3,LOP,ICHAN)*H2OJACPRD(2,3,LOP) +
     $          COFH2O(4,LOP,ICHAN)*H2OJACPRD(2,4,LOP) +
     $          COFH2O(5,LOP,ICHAN)*H2OJACPRD(2,5,LOP) +
     $          COFH2O(6,LOP,ICHAN)*H2OJACPRD(2,6,LOP) +
     $          COFH2O(7,LOP,ICHAN)*H2OJACPRD(2,7,LOP) +
     $          COFH2O(8,LOP,ICHAN)*H2OJACPRD(2,8,LOP) +
     $          COFH2O(9,LOP,ICHAN)*H2OJACPRD(2,9,LOP)
C            Remove WAOP scaling factor, which depends on wazop (sarta_pclsam.f) which is read in by rdcoef.f
C            WAOP depends on WAZOP(L)-WAZOP(L-1), and WAZOP are from a file so constant!
             KWOP_1(LOP) = KWOP_1(LOP)/WAOP(LOP)
C            Check for negative value
c             IF (KWOP_1(LOP) .LT. 0.0E+0) KWOP_1(LOP)=0.0E+0
          ENDIF
       ENDDO
C
C      -------------------------
C      Loop over the AIRS layers
C      -------------------------
       DO L=1,LBOT
C
CCC       catch bug: KWOP_1(0)
          IF (LOPLOW(L) .GT. 0.0E+0) THEN
          
C         Interpolate abs coef and convert to optical depth
c    ycalowp.f shows WAANG(L)=WAMNT(L)*SECANG(L)
            KW_1(L)=( DAOP(L)*( KWOP_1(LOPLOW(L) + 1) -
     $       KWOP_1(LOPLOW(L)) ) + KWOP_1(LOPLOW(L)) )*WAANG(L)

c    ycalowp.f shows DAOP depends on T and WV
            KW_1(L)=KW_1(L) + ( DAOPJAC(2,L)*( KWOP(LOPLOW(L) + 1) -
     $       KWOP(LOPLOW(L)) ) + KWOP(LOPLOW(L)) )*WAANG(L)

c    ycalowp.f shows WAANG(L)=WAMNT(L)*SECANG(L) so d WAANG(L)/dQ = SECANG(L)
            KW_1(L)=KW_1(L) + ( DAOP(L)*( KWOP(LOPLOW(L) + 1) -
     $       KWOP(LOPLOW(L)) ) + KWOP(LOPLOW(L)) )*SECANG(L)

c            IF (KW_T(L) .LT. 0.0E+0) KW(L)=0.0E+0
C
         ENDIF
       ENDDO
cbaba TJAC
       END IF

       RETURN
       END
