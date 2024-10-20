C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              WEXSVP
C
!F77====================================================================

!ROUTINE NAME: WEXSVP

!ABSTRACT: Hyland-WEXler's Saturation Vapor Pressure (over water)

!CALL PROTOCOL:
C    WEXSVP( T )

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      T       temperature                 Kelvin

!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      WEXSVP  saturation vapor pressure   millibar = hPa

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S):
C    INTLEV

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C      Hyland & Wexler, 1983, equation for SVP over water
C      (equation taken from Holger Voemel's web page at:
C        http://cires.colorado.edu/~voemel/vp.html )
C
C      log Pw =  -0.58002206E+4/T
C               + 0.13914993E+1
C               - 0.48640239E-1 T
C               + 0.41764768E-4 T^2
C               - 0.14452093E-7 T^3
C               + 0.65459673E+1 log(T)
C      with T in [K] and Pw in [Pa]

!ALGORITHM REFERENCES: see DESCRIPTION
C         Hyland, R. W. and A. Wexler,
C         Formulations for the Thermodynamic Properties of the
C            saturated Phases of H2O from 173.15K to 473.15K,
C         ASHRAE Trans, 89(2A), 500-519, 1983.

!KNOWN BUGS AND LIMITATIONS:
C      Temperature limited to 173-473K

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 23 Jul 2003 Scott Hannon/UMBC created to replace P.Flatau's polynominal
C                               fit of this equation which breaks down
C                               below 188K.

!END====================================================================

C      =================================================================
       REAL FUNCTION WEXSVP(T)
C      =================================================================

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
       REAL T

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       DOUBLE PRECISION LOGP
       DOUBLE PRECISION TEMP

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
C      Do the math using all doubles until the last step
       TEMP = DBLE(T)
       LOGP = -0.58002206D+4 / TEMP
     $       + 0.13914993D+1
     $       - 0.48640239D-1 * TEMP
     $       + 0.41764768D-4 * TEMP**2
     $       - 0.14452093D-7 * TEMP**3
     $       + 0.65459673D+1 * DLOG(TEMP)
       WEXSVP = SNGL( 1.0D-2 * DEXP( LOGP ) )
C
C      Hyland & Wexler, 1983, equation for SVP over ice
C           log Pi =  -0.56745359E+4 / T
C                    + 0.63925247E+1
C                    - 0.96778430E-2  T
C                    + 0.62215701E-6  T^2
C                    + 0.20747825E-8  T^3
C                    - 0.94840240E-12 T^4
C                    + 0.41635019E+1 log(T)
C       with T in [K] and Pi in [Pa]
C
       RETURN
       END
