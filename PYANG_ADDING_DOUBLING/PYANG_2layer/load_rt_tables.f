
      SUBROUTINE load_rt_tables( input_filename_ice,input_filename_liq  )

      IMPLICIT NONE
      INCLUDE 'pyang2cld.param'

c------------
c Arguments
c------------
      CHARACTER(*) input_filename_ice
      CHARACTER(*) input_filename_liq

         INTEGER TAU_idx,TAU_idx_,D_idx,D_idx_,THET_idx,THET_idx_
         REAL THET(THETMAX)
         REAL RT_ICE( 2 * THETMAX, NUM_ICE, TAUMAX_ICE * DMAX_ICE )
         REAL RT_LIQ( 2 * THETMAX, NUM_LIQ, TAUMAX_LIQ * DMAX_LIQ )
         COMMON /common1/TAU_idx,TAU_idx_,D_idx,D_idx_,THET_idx,THET_idx_
         COMMON /common2/THET,RT_ICE,RT_LIQ

c------------------
c Local variables
c------------------
      INTEGER i, j, k, IOSTAT

      IOSTAT = 0

      OPEN( 16, FILE = input_filename_ice, RECL = NUM_ICE * THETMAX * 2 * 4, ACCESS = 'DIRECT', 
     &       FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = IOSTAT )
      IF ( IOSTAT > 0 ) THEN
        print *, ' ** OPEN ERROR in load_rt_tables for file ', input_filename_ice, ' *** '
        STOP
      END IF 

      DO i = 1, TAUMAX_ICE * DMAX_ICE
        DO j = 1, NUM_ICE
          DO k = 1, 2 *THETMAX
            READ( 16, REC = i, ERR = 99 ) RT_ICE( k, j, i )
            END DO
          END DO
      END DO      
      CLOSE( 16 )

      OPEN( 16, FILE = input_filename_liq, RECL = NUM_LIQ * THETMAX * 2 * 4, ACCESS = 'DIRECT',
     &       FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = IOSTAT )
      IF ( IOSTAT > 0 ) THEN
        print *, ' ** OPEN ERROR in load_rt_tables for file ', input_filename_liq, ' *** '
        STOP
      END IF 

      DO i = 1, TAUMAX_LIQ * DMAX_LIQ
        DO j = 1, NUM_ICE
          DO k = 1, 2 *THETMAX
            READ( 16, REC = i, ERR = 99 ) RT_LIQ( k, j, i )
            END DO
          END DO
      END DO      
      CLOSE( 16 )

 99   PRINT *, ' ** READ ERROR in load_rt_tables ** '
      STOP

      RETURN
      END
