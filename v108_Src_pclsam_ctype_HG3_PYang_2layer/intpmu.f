
c------------------------------------------------------------
c  PROGRAM    INTPMU SUBROUTINE
c------------------------------------------------------------
c  PURPOSE    SUBROUTINE TO MAKE R ! T INTERPOLATION
c             FOR CORRESPONDING REQUIRED WAVENUMBER
c------------------------------------------------------------
c  VERSION    1.0      JIANGUO NIU      2004/08/04
c------------------------------------------------------------
c  DESCRIPTION 
c------------------------------------------------------------
c  VARIABLE   DESCRIPTION
c  nu         WAVENUMBER GRID IN CLD DATABASE
c  inu1       Start wavenumber index
c  inu2       End wavenumber index
c  NUM        TOTAL WAVENUMBERS IN R&T DATABASE
c  RQWV       REQUIRED WAVENUMBER GRID FOR INTERPOLATION
c  R1         ORIGINAL R DATA FOR WAVENUMBER INTERPOLATION
c  T1         ORIGINAL T DATA FOR WAVENUMBER INTERPOLATION
c  R2         ORIGINAL R DATA FOR WAVENUMBER INTERPOLATION
c  T2         ORIGINAL T DATA FOR WAVENUMBER INTERPOLATION
c  R          OUTPUT R ON THE REQUIRED WAVENUMBER GRID
c  T          OUTPUT T ON THE REQUIRED WAVENUMBER GRID
c------------------------------------------------------------
       SUBROUTINE INTPMU( nu, inu1, inu2, R1, T1, R2, T2, RQWV, R, T)

          IMPLICIT NONE
          INCLUDE 'pyang2cld.param'
c--------------
c Arguments
c--------------

          REAL nu(*)
          INTEGER inu1
          INTEGER inu2
          REAL R1
          REAL T1
          REAL R2
          REAL T2
          REAL RQWV
          REAL R
          REAL T

c------------------
c Local variables
c------------------
          
          INTEGER i, jT    

c
c Don't bother to interpolate if input wavennumber occurs at either endpoint of lookup table
c
          IF ( nu( inu1 ) .EQ. RQWV ) THEN
             R = R1
             T = T1
             RETURN
          END IF
          IF ( nu( inu2 ) .EQ. RQWV ) THEN
             R = R2
             T = T2
             RETURN
          END IF

          IF ( inu2 .NE. inu1 ) THEN
            R = R1 + ( R2 - R1 ) * 
     &  ( RQWV - nu( inu1 ) ) / ( nu( inu2 ) - nu( inu1 ) )
            T = T1 + ( T2 - T1 ) * 
     &  ( RQWV - nu( inu1 ) ) / ( nu( inu2 ) - nu( inu1 ) )
          ELSE
            R = R2
            T = T2
          END IF

       RETURN
       END
