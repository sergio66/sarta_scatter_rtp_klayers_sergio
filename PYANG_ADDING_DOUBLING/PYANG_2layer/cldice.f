
c--------------------------------------------------------
c PROGRAM      CLDICE SUBROUTINE
c--------------------------------------------------------
c PURPOSE      Reading ice cloud reflectance, transmittance
c              and interpolating them if necessary
c--------------------------------------------------------
c VERSION:     1.0  JIANGUO NIU  2004/07/12
c--------------------------------------------------------
c DISCRIPTION  CLD LOOKUP TABLE IS COMPUTED IN THE 
c              WAVENUMBER REGION FROM 587.30 cm-1 
c              TO 2349.50 cm-1 WITH THE INTERVAL OF
c              0.6 cm-1, 2938 data points.
c--------------------------------------------------------
c SUBROUTINES  INTPTH INTPMU
c--------------------------------------------------------
c INPUT PARAMETERS:
c RQWV         REQUIRED WAVENUMBER TO INTERPOLATE 
c TAU_vis      CLD VISIBLE OPTICAL THICKNESS
c De           CLD PARTICLE EFFECTIVE SIZE
c THETA        SATELITTE ZENITH ANGLE
       
c OUTPUTS:
c R            CLD REFLECTANCE
c T            CLD TRANSMITTANCE
       
c NUM          TOTAL WAVENUMBERS IN R&T LOOKUP TABLE
c TAUMAX       TOTAL DIFFERENT TAUvis
c DMAX         TOTAL DIFFERENT De
c THETMAX      TOTAL DIFFERENT THETA ANGLES
c--------------------------------------------------------

       SUBROUTINE CLDICE( RQWV,TAU_vis,De,THETA,R,T)

          IMPLICIT NONE
          INCLUDE 'pyang2cld.param'

c-------------
c  Arguments
c-------------
          REAL RQWV
          REAL TAU_vis
          REAL De
          REAL THETA
          REAL R
          REAL T

c-----------------
c Local variables
c-----------------
          INTEGER indexT
          INTEGER indexT_
          INTEGER indexD
          INTEGER indexD_
          INTEGER i,j
          INTEGER index
          REAL TAU_v(TAUMAX_ICE)
          REAL D(DMAX_ICE)
          REAL nu(NUM_ICE)
          REAL    RTT1( 2 * THETMAX )
          REAL    RTT2( 2 * THETMAX )
          REAL    RTD1( 2 * THETMAX )
          REAL    RTD2( 2 * THETMAX )
          REAL RD1
          REAL TD1
          REAL RD2
          REAL TD2
          REAL RT1
          REAL TT1
          REAL RT2
          REAL TT2
          REAL Rtemp1
          REAL Ttemp1
          REAL Rtemp2
          REAL Ttemp2
          REAL RT( 2 * THETMAX)
          REAL R1
          REAL T1
          REAL R2
          REAL T2
          INTEGER index_wn
          INTEGER index_wn1
          INTEGER index_wn2
          LOGICAL first

         INTEGER TAU_idx,TAU_idx_,D_idx,D_idx_,THET_idx,THET_idx_
         REAL THET(THETMAX)
         REAL RT_ICE( 2 * THETMAX, NUM_ICE, TAUMAX_ICE * DMAX_ICE )
         REAL RT_LIQ( 2 * THETMAX, NUM_LIQ, TAUMAX_LIQ * DMAX_LIQ )
         COMMON /common1/TAU_idx,TAU_idx_,D_idx,D_idx_,THET_idx,THET_idx_
         COMMON /common2/THET,RT_ICE,RT_LIQ

c----------
c  Data
c----------
          DATA D / 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150.,
     &             160., 170., 180. /
          DATA TAU_v / 0.01, 0.05, 0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
     &                 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0 / 
          DATA first / .TRUE. /
          SAVE


c Compute wavenumber grid used in tables first time through only
          IF ( first ) THEN
            DO i = 1, NUM_ICE
              nu( i ) = 100.0 + 1.0 * ( i - 1 )
            END DO
            DO i = 1, THETMAX
              THET( i ) =  10.0 * ( i - 1 ) 
            END DO
            first = .FALSE.
          END IF
  
          IF ( TAU_vis < TAU_v( 1 ) ) THEN
            PRINT*, 'INPUT TAU SMALLER THAN LIMIT (0.01)'
            GOTO 1000
          END IF
           
          IF ( TAU_vis > TAU_v( TAUMAX_ICE ) ) THEN
            TAU_idx = TAUMAX_ICE
            TAU_idx_= TAUMAX_ICE
            GOTO 100
          END IF
          DO i = 1, TAUMAX_ICE
            IF (TAU_v( i ) == TAU_vis ) THEN
              TAU_idx = i
              TAU_idx_= i
              GOTO 100    
            ELSE IF ( TAU_v( i ) > TAU_vis ) THEN
              TAU_idx = i
              TAU_idx_= i-1
              GOTO 100
            END IF
          END DO
         
 100       IF ( De > D( DMAX_ICE ) ) THEN  !  c Check added by TG 3/1/06
            D_idx = DMAX_ICE              
            D_idx_= DMAX_ICE              
            GOTO 110                      
          END IF                          
          DO i = 1, DMAX_ICE
            IF ( D( i ) == De ) THEN
              D_idx = i
              D_idx_= i
              GOTO 110
            ELSE IF ( D( i ) > De ) THEN
              D_idx = i
              D_idx_= i-1
              GOTO 110
            END IF
          END DO
       
 110      CONTINUE

c-------------------------
c Get wavenumber indices
c-------------------------
c          index_wn = int( ( RQWV - 587.3 ) / 0.6 + 1.5 )
          index_wn = int( ( RQWV - 100.0 ) + 1.5 )
          IF ( RQWV > nu( NUM_ICE ) ) index_wn = NUM_ICE
          index_wn1 = index_wn
          IF ( index_wn == NUM_ICE ) THEN
            index_wn2 = index_wn1
          ELSE
            index_wn2 = index_wn1 + 1
          END IF

c--------------------------------------------------------
c  TAUvis and De are on the grid of cld table
c--------------------------------------------------------
          IF ( ( TAU_idx == TAU_idx_ ) .AND. ( D_idx == D_idx_ ) ) THEN
            index = DMAX_ICE * ( TAU_idx - 1 ) + D_idx
            DO j = 1,  2 * THETMAX
              RT(j) = RT_ICE( j, index_wn1, index )
              END DO
            CALL INTPTH( THETA, RT, R1, T1 )
            DO j = 1,  2 * THETMAX
              RT(j) = RT_ICE( j, index_wn2, index )
              END DO
            CALL INTPTH( THETA, RT, R2, T2 )
            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF
c--------------------------------------------------------
c TAUvis on the grid of table, De not on
c--------------------------------------------------------
          IF ( (TAU_idx == TAU_idx_ ) .AND. ( D_idx /= D_idx_ ) ) THEN
            indexD_ = DMAX_ICE * ( TAU_idx - 1 ) + D_idx_
            indexD  = DMAX_ICE * ( TAU_idx - 1 ) + D_idx
            DO j = 1,  2 * THETMAX
              RTD1(j) = RT_ICE( j, index_wn1, indexD_ )
              END DO
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            DO j = 1,  2 * THETMAX
              RTD2(j) = RT_ICE( j, index_wn1, indexD )
              END DO
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            R1 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            T1 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1

            DO j = 1,  2 * THETMAX
              RTD1(j) = RT_ICE( j, index_wn2, indexD_ )
              END DO
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            DO j = 1,  2 * THETMAX
              RTD2(j) = RT_ICE( j, index_wn2, indexD )
              END DO
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            R2 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            T2 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1

            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF
c--------------------------------------------------------
c BOTH TAUvis and De not on the grid of cld table
c--------------------------------------------------------
          IF ( ( TAU_idx /= TAU_idx_ ) .AND. ( D_idx /= D_idx_ ) ) THEN
            indexT_= DMAX_ICE * ( TAU_idx_ - 1 ) + D_idx_
            indexT = DMAX_ICE * ( TAU_idx_ - 1 ) + D_idx
            indexD_= DMAX_ICE * ( TAU_idx - 1 ) + D_idx_
            indexD = DMAX_ICE * ( TAU_idx - 1 ) + D_idx
            DO j = 1,  2 * THETMAX
              RTT1(j) = RT_ICE( j, index_wn1, indexT_ )
              END DO
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            DO j = 1,  2 * THETMAX
              RTT2(j) = RT_ICE( j, index_wn1, indexT )
              END DO
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            Rtemp1 = ( De - D( D_idx_ ) ) * ( RT2 - RT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RT1
            Ttemp1 = ( De - D( D_idx_ ) ) * ( TT2 - TT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TT1
            DO j = 1,  2 * THETMAX
              RTD1(j) = RT_ICE( j, index_wn1, indexD_ )
              END DO
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            DO j = 1,  2 * THETMAX
              RTD2(j) = RT_ICE( j, index_wn1, indexD )
              END DO
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            Rtemp2 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            Ttemp2 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1
            R1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Rtemp2 - Rtemp1 ) /
     &                  ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Rtemp1
            T1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Ttemp2 - Ttemp1 ) /  
     &                  ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Ttemp1
            DO j = 1,  2 * THETMAX
              RTT1(j) = RT_ICE( j, index_wn2, indexT_ )
              END DO
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            DO j = 1,  2 * THETMAX
              RTT2(j) = RT_ICE( j, index_wn2, indexT )
              END DO
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            Rtemp1 = ( De - D( D_idx_ ) ) * ( RT2 - RT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RT1
            Ttemp1 = ( De - D( D_idx_ ) ) * ( TT2 - TT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TT1
            DO j = 1,  2 * THETMAX
              RTD1(j) = RT_ICE( j, index_wn2, indexD_ )
              END DO
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            DO j = 1,  2 * THETMAX
              RTD2(j) = RT_ICE( j, index_wn2, indexD )
              END DO
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            Rtemp2 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            Ttemp2 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1
            R2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Rtemp2 - Rtemp1 ) / 
     &                  ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Rtemp1
            T2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Ttemp2 - Ttemp1 ) /  
     &                  ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Ttemp1
            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R , T )

          END IF
c--------------------------------------------------------
c TAUvis not on the grid, De on the grid of cld table
c--------------------------------------------------------
          IF ( ( TAU_idx /= TAU_idx_ ) .AND. ( D_idx == D_idx_ ) ) THEN
            indexT_= DMAX_ICE * ( TAU_idx_ - 1 ) + D_idx
            indexT = DMAX_ICE * ( TAU_idx  - 1 ) + D_idx
            DO j = 1,  2 * THETMAX
              RTT1(j) = RT_ICE( j, index_wn1, indexT_ )
              END DO
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            DO j = 1,  2 * THETMAX
              RTT2(j) = RT_ICE( j, index_wn1, indexT )
              END DO
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            R1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * (RT2 - RT1 ) / 
     &                  ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + RT1
            T1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( TT2 - TT1 ) / 
     &                  ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + TT1

            DO j = 1,  2 * THETMAX
              RTT1(j) = RT_ICE( j, index_wn2, indexT_ )
              END DO
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            DO j = 1,  2 * THETMAX
              RTT2(j) = RT_ICE( j, index_wn2, indexT )
              END DO
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            R2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * (RT2 - RT1 ) / 
     &                   ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + RT1
            T2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( TT2 - TT1 ) / 
     &                  ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + TT1

            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF

 1000  RETURN
       END
