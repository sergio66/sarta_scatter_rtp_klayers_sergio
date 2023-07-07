C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              MARS
C
C              GRAV
C
!F77====================================================================

!ROUTINE NAME: GRAV

!ABSTRACT: Function for computing Mars's gravity.

!CALL PROTOCOL:
C    GRAV(Z, WINDE, WINDN, LAT, LON)

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      LAT     latitude                    degrees          
C    REAL      LON     longitude                   degrees
C    REAL      WINDE   wind velecity east          m/s
C    REAL      WINDN   wind velocity north         m/s
C    REAL      Z       altitude                    m

!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      GRAV    Mars gravity                 m/s^2

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S):
C    INTLEV

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C    Function to calculate Mars gravity (gravitation plus
C    centripetal acceleration) for points in the atmosphere.
C
C    It calculates surface gravity using an equation given in
C    "American Institute of Physics Handbook", 1963, page 2-102.
C    This equation is essentially a variant of the International
C    Gravity Formula with an extra term for longitude (which is
C    very minor).
C
C    Centripetal acceleration is tangental velocity squared over the
C    radius.
C
C    Gravitation is the gravitational constant times the Mars's mass
C    divided by the square of the radius.
C
C    Gravity at any point in the atmosphere is simply surface gravity
C    plus the change in gravitation and centripetal acceleration at
C    that point compared to the surface.

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS: none

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C Mar  8 1995 Scott Hannon/UMBC created
C Jun 23 1995 Scott Hannon      Correct some comments
C 28 Sep 2007 Scott Hannon      Add more centripetel term comments
c May 2021    Sergio Machado    Mars

!END====================================================================

C      =================================================================
       REAL FUNCTION GRAV(Z, WINDE, WINDN, LAT, LON)
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

       REAL Z, WINDE, WINDN, LAT, LON

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
C
       REAL G_SUR, R, COSLT, COSLT2, SINLT2, SIN2LT, COSLON,
     $    LTRAD, C_SUR, C_Z, RTOT, GRAVZ
C
c https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
C      Constants for (1 - b^2/a^2) with
C      a = 3,3962+6 m = equatorial radius, and
C      b = 3.3762+6 m = polar radius.
       REAL B2, ABTERM
C
C      Constants for normal gravity equation
C      (see "American Institute of Physics Handbook", 1963, pg 2-102) 
       REAL G0
       REAL C1, C2, C3
C
C      Constants for pi/180, 2*pi, and Mars's rotational speed
C      of w=1/86400 rev/s iday = 24.6597
       REAL PI180,PI2, W
C 
       DATA B2, ABTERM /1.1399e+13, 1.1743e-02/
       DATA G0 /3.721/
       DATA C1,C2,C3 /0.0, 0.0, 0.0/
       DATA PI180,PI2,W /1.7453293E-2, 6.28318531, 1.1264e-05/

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
C      Calculate longitude term
C      Add offset of 18 degrees, double it, convert to radians, and
C      take the cosine
       COSLON=COS( PI180*2.0*( LON + 18.0 ) )
C
C      Calculate the latitude terms
C      Convert Latitude into radians
       LTRAD=PI180*LAT
C      Calculate sine and cosine terms
       COSLT = COS(LTRAD)
       COSLT2 = COSLT**2
       SINLT2 = ( SIN(LTRAD ) )**2
       SIN2LT = ( SIN( 2.0*LTRAD ) )**2
C
C      Calculate the Mars's radius at this latitude
       R = SQRT( B2/( 1.0 - COSLT2*ABTERM ) )
C
C      Calculate total distance from Mars's center
       RTOT = R + Z
C
C      Calculate gravity at the Mars's surface
c      https://en.wikipedia.org/wiki/Theoretical_gravity
c       G_SUR = G0*( 1.0 + C1*SINLT2 + C2*SIN2LT + C3*COSLT2*COSLON )
       G_SUR = G0

C
C      Calculate the centripetal term at the Mars's surface
C      Note: the centripetal acceleration due to Mars's rotation
C      is in a direction perpendicular to the Mars's rational axis,
C      but for this gravity function we are only interested in the
C      component parallel to the radial direction.
       C_SUR = COSLT2*R*(PI2*W)**2
C
C      Calculate the centripetal term at altitude z (with wind)
       C_Z = ( ( PI2*RTOT*COSLT*W + WINDE )**2 + (WINDN)**2 )/RTOT
C
C      Calculate the change in gravitation with altitude
       GRAVZ=(G_SUR + C_SUR)*(1.0 - R**2/RTOT**2)
C
       GRAV = G_SUR + (C_SUR - C_Z) - GRAVZ
C
       RETURN
       END
