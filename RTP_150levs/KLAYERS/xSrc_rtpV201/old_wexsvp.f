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

!ABSTRACT: Function for WEXler's Saturation Vapor Pressure expression.

!CALL PROTOCOL:
C    WEXSVP( T )

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      T       temperature                 K

!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      WEXSVP  saturation vapor pressure   mb

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S):
C    INTLEV

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C      Function based upon the article:
C      "Polynomial Fits to Saturation Vapor Pressure"
C      by P. J. Flatau, R. L. Walko, and W. R. Cotton
C      Journal of Applied Meteorology, v.31 Dec 1992, pg 1507  
C
C      Uses an 8th-order polynomial fit to Wexler's equations of
C      the form:
C         e_sat=sum{i=1 to 9} a(i)*T^(i-1)
C      where T is in C.
C      Over water (-85 to +70 degrees C) with relative error norm.

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS:
C      I'm not sure this works right for pressures less than ~100 mb.

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C Mar  8 1995 Scott Hannon/UMBC created
C Jun 23 1995 Scott Hannon      Correct some comments

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
       INTEGER I, J
       DOUBLE PRECISION RJ1, RJ2, WCOEF(9)
C      Over water (-85 to +70 degrees C); relative error norm
       DATA (WCOEF(I),I=1,9)/
     +  6.11583699E+00, 4.44606896E-01, 1.43177157E-02, 2.64224321E-04,
     +  2.99291081E-06, 2.03154182E-08, 7.02620698E-11, 3.79534310E-14,
     + -3.21582393E-16/
C
ccc
cC      Over ice (-90 to 0 degrees); relative error norm
c       DATA (ICOEF(I),I=1,9)/
c     +  6.09868993E+00, 4.99320233E-01, 1.84672631E-02, 4.02737184E-04,
c     +  5.65392987E-06, 5.21693933E-08, 3.07839583E-10, 1.05785160E-12,
c     +  1.61444444E-15/
ccc

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE
C***********************************************************************
C***********************************************************************

C      Convert T in K to C
       RJ2=T-273.15
C
       RJ1=WCOEF(1)
       DO I=2,9
          J=I-1
          RJ1=RJ1 + ( WCOEF(I)*(RJ2**J) )
       ENDDO
C
       WEXSVP=SNGL(RJ1)
C
       RETURN
       END

