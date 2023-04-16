C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              TOPPMV
C
!F77====================================================================


!ROUTINE NAME: TOPPMV


!ABSTRACT:
C    Convert gas amount profile from whatever units to PPMV


!CALL PROTOCOL:
C    TOPPMV(NGASF, GLISTF, GUNITF, NLEV, MASSF, PIN, TIN, MRIN)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   NGASF   number of found gases       none
C    INT arr   GLISTF  list of found gases         HITRAN gas ID #'s
C    INT arr   GUNITF  gas amount unit code #'s    "GUC" code numbers
C    INTEGER   NLEV    number of profile levels    none
C    REAL arr  MASSF   molecular masses            AMU
C    REAL arr  PIN     pressure levels             mb or hPa
C    REAL arr  TIN     temperature                 Kelvin


!OUTPUT PARAMETERS:
C    none


!INPUT/OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  MRIN    mixing ratios               GUNITF in; PPMV out

!RETURN VALUES: none


!PARENT(S): klayers_rtp.f


!ROUTINES CALLED:
C    FUNCTION WEXSVP : calc H2O saturation vapor pressure


!FILES ACCESSED:
C    none


!DESCRIPTION:
C    Converts gas amount units from <whatever> to PPMV.  The
C    <whatever> units must be oneof the allowed code numbers
C    found in GUCS (see cbgucs.f)!


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 27 Feb 2001 Scott Hannon      Created


!END====================================================================


C      =================================================================
       SUBROUTINE TOPPMV(NGASF, GLISTF, GUNITF, NLEV, MASSF, PIN, TIN,
     $    MRIN)
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
C      Input:
       INTEGER  NGASF           ! # of wanted gases found in input file
       INTEGER GLISTF( MXGAS)   ! list of found wanted gases
       INTEGER GUNITF( MXGAS)   ! gas units code#'s of found gases
       INTEGER   NLEV           ! number of levels in profile
       REAL     MASSF( MXGAS)   ! molecular masses of found gases
       REAL       PIN(  MXIN)   ! molecular masses of found gases
       REAL       TIN(  MXIN)   ! molecular masses of found gases

C
C      Input/Output:
       REAL MRIN(  MXIN, MXGAS) ! gas amount profiles
C
C      Output: none


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER ICODE
       INTEGER IG
       INTEGER IL
       REAL MDAIR
       REAL RJUNK
C
       REAL WEXSVP  ! real function


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      RTP code numbers for various gas amount units
C      ---------------------------------------------
C
C      This is the recommended list of code numbers for various
C      gas amount units.  The list is not complete; it includes
C      only the most common units.  Additional units can be added
C      where ever there is a free number.
C
C      An attempt has been made to organize the gas amount units
C      into classes so that similar units have similar numbers.
C      This is merely an attempt to make it easier for users to
C      search the list for some particular unit.
C      The gas amount unit classes are as follows:
C         -9999 = no data
C             0 = unknown units
C        1 -  9 = integrated column density
C       10 - 19 = volume mixing ratio
C       20 - 29 = mass mixing ratio
C       30 - 39 = partial pressure
C       40 - 49 = water vapor humidity
C       50 - 59 = none-of-the-above oddball units
C       Negative numbers are not recommended as standard code numbers.
C       However, the user might possibly find it useful to use the
C       negative value of a code number to indicate some sort of
C       alternate interpretation of the units represented by the
C       positive valued number.  For example, the meaning of "air"
C       in mixing ratios might mean "dry air" if positive and "wet
C       air" if negative.
C
C      Note: currently the KLAYERS program only uses input=PPMV
C      and output=kilomoles/cm^2
C
C      Code#   Gas amount Units and comments
C      -----   ---------------------------------------------------------
C      -9999   no data
C
C          0   unknown units (?)
C              It's something but you don't know what
C
C              ------- 1-9 = integrated column density -----------------
C              Note: The gas amount has been integrated over some finite
C              thickness vertical slab.  This "thickness" may correspond
C              to a single slab layer, or it might be the total column
C              from the current pressure level to the top of the
C              atmosphere.  This implicit "thickness" is a separate
C              piece of info which is required to correctly interpret
C              the meaning of the gas amounts.
C
C          1   molecules per square centimeter (molecules/cm^2)
C              Layer
C
C          2   kilomoles per square centimeter (kilomoles/cm^2)
C              Layer
C
C          3   molecules per square centimeter (molecules/cm^2)
C              Total column (aka layer-to-space)
C
C          4   kilomoles per square centimeter (kilomoles/cm^2)
C              Total column (aka layer-to-space)
C
C          5   [not assigned]
C          6   [not assigned]
C          7   [not assigned
C          8   [not assigned]
C          9   [not assigned]
C
C              ------- 10-19 = volume mixing ratio ---------------------
C              Note: the meaning of "air" can vary.  There are two
C              common usages with slightly different meanings.  The most
C              obvious (but less often encountered) meaning is to refer
C              to the air as it actually exists (aka "wet air") at the
C              corresponding pressure level.  But the most common
C              meaning is to refer to some simple idealized model of
C              air.  Typically this "dry air" is assumed to be composed
C              entirely of imaginary "air" molecules with an atomic
C              mass of roughly 28.966 AMU.  This distinction between
C              "wet" and "dry" air can affect the mixing ratio in the
C              lower troposphere by a few percent.
C              
C         10   parts per million volume mixing ratio (ppmv)
C              Number of gas X molecules per 1E+6 "air" molecules
C
C         11   parts per billion volume mixing ratio (ppbv)
C              Number of gas X molecules per 1E+9 "air" molecules
C
C         12   volume mixing ratio (vmr, unitless fraction)
C              Number of gas X molecules per "air" molecule
C
C         13   [not assigned]
C         14   [not assigned]
C         15   [not assigned]
C         16   [not assigned
C         17   [not assigned]
C         18   [not assigned]
C         19   [not assigned]
C
C              ------- 20-29 = mass mixing ratio -----------------------
C
C         20   mass mixing ratio in (g/kg)
C              Grams of gas X per kilogram of "air"
C
C         21   mass mixing ratio in (g/g)
C              Grams of gas X per gram of "air"
C
C         22   [not assigned]
C         23   [not assigned]
C         24   [not assigned]
C         25   [not assigned]
C         26   [not assigned]
C         27   [not assigned]
C         28   [not assigned]
C         29   [not assigned]
C
C              ------- 30-39 = partial pressure ------------------------
C
C         30   partial pressure in millibars (mb) Note: mb=hPa
C              Pressure of gas X as it exists in the atmosphere.
C              To clarify, this means at the corresponding profile
C              temperature and pressure level total pressure.
C
C         31   partial pressure in atmospheres (atm)
C              Pressure of gas X as it exists in the atmosphere
C
C         32   [not assigned]
C         33   [not assigned]
C         34   [not assigned]
C         35   [not assigned]
C         36   [not assigned]
C         37   [not assigned]
C         38   [not assigned]
C         39   [not assigned]
C
C              ------- 40-49 = water vapor humidity units --------------
C
C         40   relative humidity in (percent)
C              100 times actual vapor pressure divided by saturation
C              vapor pressure
C
C         41   relative humidity (unitless fraction)
C              Actual vapor pressure divided by saturation vapor
C              pressure
C
C         42   dew point temperature (Kelvin)
C              Temperaure at which vapor will start to condense out
C
C         43   dew point temperature (Celcius)
C              Temperaure at which vapor will start to condense out
C
C              Possible additional units might be
C                 inches or centimenters of water vapor
C                 grams of water vapor
C
C         44   [not assigned]
C         45   [not assigned]
C         46   [not assigned]
C         47   [not assigned
C         48   [not assigned]
C         49   [not assigned]
C
C              ------- 50-59 = oddball  --------------------------------
C
C         50   Dobson Units (DU)
C              A Dobson Unit is a special type of integrated column
C              density often used for ozone profiles.  The meaning is
C              a bit obtuse:  1000 times the height (in centimeters)
C              of some column of ozone if it were squashed down to STP.
C
C         51   [not assigned]
C         52   [not assigned]
C         53   [not assigned]
C         54   [not assigned]
C         55   [not assigned]
C         56   [not assigned]
C         57   [not assigned
C         58   [not assigned]
C         59   [not assigned]
C
C              ------- 60 and above = whatever you want ----------------
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
       MDAIR=28.966  ! molecular mass of dry air
C
C      Loop over the gases
       DO IG=1,NGASF
C         Note: This routine is for GUCSTD = 10 (ppmv) only!
          ICODE=GUNITF(IG)
          IF (ICODE .NE. GUCSTD) THEN
C
C            *****************************
C            Convert profile to PPMV units
C            *****************************
C
C            parts per billion volume mixing ratio
             IF (ICODE .EQ. 11) THEN
C               PPMV = PPBV*1E-3
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*1E-3
                ENDDO
C
C            volume mixing ratio
             ELSEIF (ICODE .EQ. 12) THEN
C               PPMV = VMR*1E+6
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*1E+6
                ENDDO
C
C            mass mixing ratio in g/kg
             ELSEIF (ICODE .EQ. 20) THEN
C               PPMV = MRgkg*((MDAIR*1E-3)/MASSF)*1E+6
                RJUNK=1E+3*MDAIR/MASSF(IG)
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*RJUNK
                ENDDO
C
C            mass mixing ratio in g/g
             ELSEIF (ICODE .EQ. 21) THEN
C               PPMV = MRgg*(MDAIR/MASSF)*1E+6
                RJUNK=1E+6*MDAIR/MASSF(IG)
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*RJUNK
                ENDDO
C
C            partial pressure in mb
             ELSEIF (ICODE .EQ. 30) THEN
C               PPMV = (PPmb/PIN)*1E+6
                DO IL=1,NLEV
                   MRIN(IL,IG)=(MRIN(IL,IG)/PIN(IL))*1E+6
                ENDDO
C
C            partial pressure in atm
             ELSEIF (ICODE .EQ. 31) THEN
C               PPMV = (PPatm*1013.25/PIN)*1E+6
                DO IL=1,NLEV
                   MRIN(IL,IG)=(MRIN(IL,IG)*1013.25/PIN(IL))*1E+6
                ENDDO
C
C            relative humidy in percent
             ELSEIF (ICODE .EQ. 40 .AND. GLISTF(IG) .EQ. 1) THEN
C               PPMV = (RH%/100)*(SVP/PIN)*1E+6
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*
     $                (WEXSVP( TIN(IL) )/PIN(IL))*1E+4
                ENDDO
C
C            relative humidity (fraction)
             ELSEIF (ICODE .EQ. 41 .AND. GLISTF(IG) .EQ. 1) THEN
C               PPMV = RH*(SVP/PIN)*1E+6
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*
     $                (WEXSVP( TIN(IL) )/PIN(IL))*1E+6
                ENDDO
C
             ENDIF
C
          ENDIF
C
       ENDDO
C
       RETURN
       END

