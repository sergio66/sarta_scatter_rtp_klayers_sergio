C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              CBGUCS
C
!F77====================================================================

!ROUTINE NAME: CBGUCS

!ABSTRACT:
C    Common Block for allowed Gas amount Units Code numberS

!CALL PROTOCOL: none

!INPUT PARAMETERS: none

!OUTPUT PARAMETERS: none

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S):
C    KLAYERS

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C    Assign values to COMGUC

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS: none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 21 Feb 2001 Scott Hannon      Created
C  1 Mar 2001 Scott Hannon      Add 8 more code numbers
C 19 Feb 2007 Scott Hannon      Update some comments


!END====================================================================

C      =================================================================
       BLOCK DATA CBGUCS
C      =================================================================
C
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
C      none

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      EXECUTABLE CODE
C-----------------------------------------------------------------------
C      none; common block data
C
C      Edit as needed, but GUCS must match the units code used by
C      whoever assembled the input RTP files.  And of course KLAYERS
C      must contain the code needed to deal with the all the units.
C
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
C
C          2   kilomoles per square centimeter (kilomoles/cm^2)
C
C          3   
C
C          4   
C
C          5   grams per square meter (g/m^2)
C              Note: commonly used for clouds
C
C          6   [not assigned]
C          7   [not assigned]
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
C         12   volume mixing ratio (unitless fraction)
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
C              a bit obtuse:  1000 times the height in centimeters
C              of a gas column if it were squashed down to STP.
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

C***********************************************************************
       COMMON /COMGUC/ GUCS
C
       INTEGER GUCS(MAXGUC) ! size must match data statements below
C
C      List of all allowed input profile gas amount unit code numbers
C      for use with KLAYERS.
C      Currently only the following are allowed: (see toppmv.f)
       DATA GUCS( 1) / 10 /   ! PPMV (parts per million volume)
       DATA GUCS( 2) / 11 /   ! PPBV (parts per billion volume)
       DATA GUCS( 3) / 12 /   ! volume mixing ratio (fraction)
       DATA GUCS( 4) / 20 /   ! mass mixing ratio in g/kg 
       DATA GUCS( 5) / 21 /   ! mass mixing ratio in g/g 
       DATA GUCS( 6) / 30 /   ! partial pressure in millibars or hPa
       DATA GUCS( 7) / 31 /   ! partial pressure in atm
       DATA GUCS( 8) / 40 /   ! relative humidity in percent
       DATA GUCS( 9) / 41 /   ! relative humidity (fraction)
C
ccc
C      Easy to convert to ppmv
c         Parts per billion volume
c         Volume mixing ratio (ie PPMV without the million)
c         Partial pressure in millibar or atmospheres
c
C      Not too difficult to convert to ppmv
c         Mass mixing ratio in gram/gram (need masses for all gases)
c         Relative humidity (water only & need SVP code)
C
C      Difficult to convert to ppmv
c         Any sort of column density such as water in cm or gram or
c            ozone in Dobson units.
ccc

C
       END
