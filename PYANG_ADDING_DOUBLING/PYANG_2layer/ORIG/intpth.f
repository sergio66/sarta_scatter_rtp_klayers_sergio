
c----------------------------------------------------------------
c  PROGRAM    INTPTH  SUBROUTINE
c----------------------------------------------------------------
c  PURPOSE    VIEWING ANGLE INTERPOLATION FOR R&T OF
c             ICE ! WATER CLOUD IN DATABASE
c----------------------------------------------------------------
c  VERSION    1.0 JIANGUO NIU  2004/AUGUST/02
c----------------------------------------------------------------
c  VARIABLE   DESCRIPTION
c----------------------------------------------------------------
c  THETA      INPUT, SATELITTE VIEWING ANGLE
c  RT         INPUT, R&T ARRAY
c  R          OUTPUT, REFLECTANCES AT REQUIRED THETA
c  T          OUTPUT, TREANMITTANCES AT REQUIRED THETA
c----------------------------------------------------------------
       SUBROUTINE INTPTH( THETA,   ! c Input
                          RT,      ! c Input
                          R,       ! c Output
                          T )     !  c Ouput

          IMPLICIT NONE
          
c-----------
c Arguments 
c-----------           
          REAL THETA
          REAL RT( 2 * THETMAX )
          REAL R
          REAL T
          
c------------------
c Local variables
c------------------
          INTEGER i, k
          REAL SLP

          DO i = 1, THETMAX
            IF ( THET( i ) == THETA ) THEN
              THET_idx  = i
              THET_idx_ = i
              GOTO 210
            ELSE IF ( THET( i ) > THETA ) THEN
              THET_idx  = i
              THET_idx_ = i - 1
              GOTO 210
            END IF
          END DO

210       IF ( THET_idx == THET_idx_ ) THEN
            R = RT( THET_idx )
            T = RT( THETMAX + THET_idx )
          ELSE

c This code does not interpolate properly - TG 5/23/2006
c            R = ( THETA - THET( THET_idx_ ) ) * ( RT( THET_idx ) &
c                      - RT( THET_idx_ ) ) / ( THET( THET_idx ) - THETA &
c                      - THET( THET_idx_ ) ) + RT( THET_idx_ )
c This code does not interpolate properly - TG 5/23/2006
c            T = ( THETA - THET( THET_idx_ ) ) * ( RT( THETMAX + THET_idx ) &
c                      - RT( THETMAX + THET_idx_ ) ) / ( THET( THET_idx ) - THETA &
c                      - THET( THET_idx_ ) ) + RT( THETMAX + THET_idx_ )
            
            SLP = ( THETA - THET( THET_idx_ ) ) / ( THET( THET_idx ) - THET( THET_idx_ ) )
            R = ( 1.0_fp_kind - SLP ) * RT( THET_idx_ ) + SLP * RT( THET_idx )
            T = ( 1.0_fp_kind - SLP ) * RT( THETMAX + THET_idx_) + SLP * RT( THETMAX + THET_idx )

c    print *, ' R = ', r
c    print *, ' rt = ', rt(thet_idx), rt(thet_idx_)
c    print *, ' T = ', t
c    print *, ' THETA = ', theta, thet(thet_idx), thet(thet_idx_)
c    print *, ' rt = ', rt(thetmax+thet_idx), rt(thetmax+thet_idx_)

          END IF

       END SUBROUTINE intpth

