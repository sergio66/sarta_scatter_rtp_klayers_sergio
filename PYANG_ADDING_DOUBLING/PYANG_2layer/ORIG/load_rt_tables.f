
   SUBROUTINE load_rt_tables( input_filename_ice, ! c Input
                              input_filename_liq  ) c Input

      IMPLICIT NONE
      
c------------
c Arguments
c------------
      CHARACTER(*) input_filename_ice
      CHARACTER(*) input_filename_liq

c------------------
c Local variables
c------------------
      INTEGER i, IOSTAT

      IOSTAT = 0

      OPEN( 16, FILE = input_filename_ice, RECL = NUM_ICE * THETMAX * 2 * 4, ACCESS = 'DIRECT', &
            FORM = TRIM( FILE_FORMAT ), STATUS = 'OLD', IOSTAT = IOSTAT )
      IF ( IOSTAT > 0 ) THEN
        print *, ' ** OPEN ERROR in load_rt_tables for file ', trim(input_filename_ice), ' *** '
        STOP
      END IF 

      DO i = 1, TAUMAX_ICE * DMAX_ICE
        READ( 16, REC = i, ERR = 99 ) RT_ICE( :, :, i )
      END DO      
      CLOSE( 16 )

      OPEN( 16, FILE = input_filename_liq, RECL = NUM_LIQ * THETMAX * 2 * 4, ACCESS = 'DIRECT', &
            FORM = TRIM( FILE_FORMAT ), STATUS = 'OLD', IOSTAT = IOSTAT )
      IF ( IOSTAT > 0 ) THEN
        print *, ' ** OPEN ERROR in load_rt_tables for file ', trim(input_filename_liq), ' *** '
        STOP
      END IF 

      DO i = 1, TAUMAX_LIQ * DMAX_LIQ
        READ( 16, REC = i, ERR = 99 ) RT_LIQ( :, :, i )
      END DO      
      CLOSE( 16 )

      RETURN
      
 99   PRINT *, ' ** READ ERROR in load_rt_tables ** '
      STOP

   END SUBROUTINE load_rt_tables

