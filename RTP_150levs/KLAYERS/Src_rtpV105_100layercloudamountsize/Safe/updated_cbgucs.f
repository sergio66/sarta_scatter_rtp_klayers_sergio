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
C 20 Feb 2007 Scott Hannon      Add 5 more code numbers
C 22 Aug 2007 Scott Hannon      Add code number 55(microns) for cloud
C                                  particle size; add 5 & 55 to GUCS


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
C      Last update: 20 February 2007, Scott Hannon
C
C      Note: aka = Also Known As
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
C
C      The gas amount unit classes are as follows:
C         -9999 = no data
C             0 = unknown units
C        1 -  9 = integrated column density (aka cross-sectional density)
C       10 - 19 = volume mixing ratio
C       20 - 29 = mass mixing ratio
C       30 - 39 = partial pressure
C       40 - 49 = water vapor humidity
C       50 - 59 = none-of-the-above oddball units
C
C      Note: The Earth's "air" is a composite of various molecules
C      and atoms whos proportions vary in time and space.  There
C      is an implicit assumed composition of "air" in some of the
C      gas amount units, and not all gas units are based on the same
C      assumption.  There are two commonly used air types which
C      I will refer to as "dry air" and "wet air". Below about
C      100 km altitude the Earth's air has a fairly uniform bulk
C      composition, and the air can be approximated as composed
C      of imaginary "dry air" molecules (with mass 28.97 AMU) plus
C      minor trace constituents.  Over 99.9% of "dry air" is
C      composed of N2, O2, and argon, all of which are nearly
C      inactive in the infra-red portion of the spectrum.  The
C      only IR active gas which typically varies enough to affect
C      the overall bulk composition of air by more than a tiny
C      fraction of a percent is water vapor.  The distinction
C      between "wet" and "dry" air can affect the mixing ratio in
C      the lower troposphere by a few percent.
C
C      Code#   Gas amount Units and comments
C      -----   ---------------------------------------------------------
C      -9999   no data
C
C          0   unknown units (?)
C              It's something but you don't know what
C
C              ------- 1-9 = integrated column density -----------------
C              AKA cross-sectional density.
C              The gas amount in a particular volume of air has been
C              integrated over its vertical thickness.  This implicit
C              "thickness" is a separate piece of info which is required
C              to correctly interpret the meaning of the gas amount.
C
C          1   molecules per square centimeter (molecules/cm^2)
C
C          2   kilomoles per square centimeter (kilomoles/cm^2)
C
C          3   [not assigned]
C          4   [not assigned]
C
C          5   grams per square meter (g/m^2)
C              Note: typically used for clouds
C
C          6   [not assigned]
C          7   [not assigned
C          8   [not assigned]
C          9   [not assigned]
C
C              ------- 10-19 = volume mixing ratio ---------------------
C
C         10   Parts Per Million Volume mixing ratio (ppmv), dry air
C              Number of gas X molecules per 1E+6 "dry air" molecules
C
C         11   Parts Per Billion Volume mixing ratio (ppbv), dry air
C              Number of gas X molecules per 1E+9 "dry air" molecules
C
C         12   Volume Mixing Ratio (vmr), dry air
C              Fraction of a gas X molecule per "dry air" molecule
C
C         13   [not assigned]
C         14   [not assigned]
C
C         15   Parts Per Million Volume mixing ratio (ppmv), wet air
C         16   Parts Per Billion Volume mixing ratio (ppmv), wet air
C         17   Volume mixing Ratio (vmr), wet air
C
C         18   [not assigned]
C         19   [not assigned]
C
C              ------- 20-29 = mass mixing ratio -----------------------
C              AKA specific humidity (which is usually but not always
C              defined as a wet air mixing ratio).
C
C         20   mass mixing ratio in (g/kg), dry air
C              Grams of gas X per kilogram of "dry air"
C
C         21   mass mixing ratio in (g/g) or (kg/kg), dry air
C              Grams of gas X per gram of "dry air"
C
C         22   [not assigned]
C         23   [not assigned]
C         24   [not assigned]
C
C         25   mass mixing ratio in (g/kg), wet air
C         26   mass mixing ratio in (g/g) or (kg/kg), wet air
C
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
C
C         55   micron (um) {10^-6 meters}
C              Note: typically used for cloud particle size
C
C         56   [not assigned]
C         57   [not assigned]
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
       DATA GUCS( 1) /  5 /   ! g/m^2 {cloud amount}
       DATA GUCS( 2) / 10 /   ! dry air PPMV (parts per million volume)
       DATA GUCS( 3) / 11 /   ! dry air PPBV (parts per billion volume)
       DATA GUCS( 4) / 12 /   ! dry air volume mixing ratio (fraction)
       DATA GUCS( 5) / 15 /   ! wet air PPMV (parts per million volume)
       DATA GUCS( 6) / 16 /   ! wet air PPBV (parts per billion volume)
       DATA GUCS( 7) / 17 /   ! wet air volume mixing ratio (fraction)
       DATA GUCS( 8) / 20 /   ! dry air mass mixing ratio in g/kg 
       DATA GUCS( 9) / 21 /   ! dry air mass mixing ratio in g/g 
       DATA GUCS(10) / 25 /   ! wet air mass mixing ratio in g/kg 
       DATA GUCS(11) / 26 /   ! wet air mass mixing ratio in g/g 
       DATA GUCS(12) / 30 /   ! partial pressure in millibars or hPa
       DATA GUCS(13) / 31 /   ! partial pressure in atm
       DATA GUCS(14) / 40 /   ! relative humidity in percent
       DATA GUCS(15) / 41 /   ! relative humidity (fraction)
       DATA GUCS(16) / 55 /   ! microns {cloud particle size}
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
