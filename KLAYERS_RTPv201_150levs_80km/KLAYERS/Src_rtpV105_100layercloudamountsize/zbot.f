C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              ZBOT
C
!F77====================================================================

!ROUTINE NAME: ZBOT

!ABSTRACT: Function for computing altitude at bottom of bottom layer

!CALL PROTOCOL:
C ZBOT(WMRAVG, TAVG, GAVG, PBOT, PTOP, PSURF, ZSURF)

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      WMRAVG  wet air H2O avg mix ratio   ppmv
C    REAL      TAVG    layer average temperature   Kelvin
C    REAL      GAVG    layer average gravity       m/s^2
C    REAL      PBOT    pres at layer bot boundary  mb
C    REAL      PTOP    pres at layer top boundary  mb
C    REAL      PSURF   surface pressure            mb
C    REAL      ZSURF   surface altitude            meters


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      ZBOT    bottom boundary altitutde   meters


!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S):
C    INTLEV

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C    Function to calculate altitude at the bottom of the bottom layer.
C    This calc is approximate, but should be fairly accurate.  First
C    it calcs the thickness of the bottom layer.  Then it assumes
C    altitude varies linearly with ln(pressure), and uses this to
C    figure out the altitude at the bottom of the layer.
C
C      Layer thickness is approximately
C         dz = (dP * T)/(g * M/R * P)
C      where
C         dz = layer thickness (meters)
C         dP = difference in layer boundary pressures (Pascals)
C         T  = layer average temperature (Kelvin)
C         g  = layer average gravity (meters per second^2)
C         M  = layer average molar mass of air molecules (kilograms)
C         R  = molar gas constant (8.314472 Joules per mole per K)
C         P  = layer average pressure (Pascals)
C
C      With the thickness determined, it is now a simple interpolation
C      to figure out the pressure at the bottom of the layer.  We
C      assume ln(pressure) varies linearly with altitude, that is
C         Z = a * ln(pressure) + b where a and b are unknown constants.
C      Therefore the ratio
C         ( ln(Psurf) - ln(Pbot) )/( ln(Ptop) - ln(Pbot) )
C      is the same as
C         ( Zsurf - Zbot )/(Ztop - Zbot)
C      where Ztop - Zbot is the layer thickness.  Thus we can solve
C      for the one unknown, Zbot. 
C

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS: none

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 23 Mar 2001 Scott Hannon      created


!END====================================================================

C      =================================================================
       REAL FUNCTION ZBOT(WMRAVG, TAVG, GAVG, PBOT, PTOP, PSURF, ZSURF)
C      =================================================================
C
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE

C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input:
       REAL WMRAVG ! wet air water layer average mixing ratio (ppmv)
       REAL   TAVG ! layer average temperature (Kelvin)
       REAL   GAVG ! layer average gravity (m/s^2)
       REAL   PBOT ! pressure at layer bottom boundary level (mb)
       REAL   PTOP ! pressure at layer top boundary level (mb)
       REAL  PSURF ! surface pressure (mb)
       REAL  ZSURF ! surface altitude (meters)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       REAL R     ! molar gas constant (joules per mole per Kelvin)
       REAL M_DRY ! mass of 1 mole of dry air (kilograms)
       REAL M_WAT ! mass of 1 mole of water vapor (kilograms)
       REAL M_AIR ! mass of 1 mole of wet air (kilograms)
       REAL DELP  ! delta-P (Pascals)
       REAL PAVG  ! average P (Pascals)
       REAL DELZ  ! layer thickness (meters)


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************
C
       R = 8.314472
       M_DRY = 28.966E-3
       M_WAT= 18.016E-3
C
C      Convert pressure from mb=hPa to Pa
       DELP = (PBOT - PTOP)*100
       PAVG = 100*( PBOT - PTOP ) / LOG( PBOT/PTOP )
C
C      Average wet air mass
       M_AIR = (1 - (WMRAVG/1E+6))*M_DRY + (WMRAVG/1E+6)*M_WAT
C
C      Calc the layer thickness
       DELZ = DELP * TAVG  / (GAVG * (M_AIR/R) * PAVG )
C
C      Solve for ZBOT
       ZBOT=ZSURF - DELZ*( LOG(PSURF/PBOT) / LOG(PTOP/PBOT) )
C
       RETURN
       END
