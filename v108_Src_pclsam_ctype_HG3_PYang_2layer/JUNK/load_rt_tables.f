
      SUBROUTINE load_rt_tables( input_filename_ice,input_filename_liq )

      IMPLICIT NONE
      INCLUDE 'pyang2cld.param'

c------------
c Arguments
c------------
      CHARACTER*256 input_filename_ice
      CHARACTER*256 input_filename_liq

         INTEGER TAU_idx,TAU_idx_,D_idx,D_idx_,THET_idx,THET_idx_
         REAL THET(THETMAX)
         REAL RT_ICE( 2 * THETMAX, NUM_ICE, TAUMAX_ICE * DMAX_ICE )
         REAL RT_LIQ( 2 * THETMAX, NUM_LIQ, TAUMAX_LIQ * DMAX_LIQ )
         COMMON /common1/TAU_idx,TAU_idx_,D_idx,D_idx_,
     $               THET_idx,THET_idx_
         COMMON /common2/THET,RT_ICE,RT_LIQ

c------------------
c Local variables
c------------------
      INTEGER i, j, k, ISTAT

      ISTAT = 0

      OPEN( 16, FILE = input_filename_ice, 
     &       RECL = NUM_ICE * THETMAX * 2 * 4, ACCESS = 'DIRECT', 
     &       FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = ISTAT )
      IF ( ISTAT > 0 ) THEN
        print *, ' ** OPEN ERROR for file ', input_filename_ice, ' *** '
        STOP
      END IF 

      DO i = 1, TAUMAX_ICE * DMAX_ICE
        READ( 16, REC = i, ERR = 98 ) 
     $     ((RT_ICE(k,j,i),k=1,2 *THETMAX),j=1,NUM_ICE)
        END DO      
      CLOSE( 16 )

      OPEN( 16, FILE = input_filename_liq, 
     &       RECL = NUM_LIQ * THETMAX * 2 * 4, ACCESS = 'DIRECT',
     &       FORM = 'UNFORMATTED', STATUS = 'OLD', IOSTAT = ISTAT )
      IF ( ISTAT > 0 ) THEN
        print *, ' ** OPEN ERROR for file ', input_filename_liq, ' *** '
        STOP
      END IF 

      DO i = 1, TAUMAX_LIQ * DMAX_LIQ
        READ( 16, REC = i, ERR = 99 ) 
     $    ((RT_LIQ(k,j,i),k=1,2 *THETMAX),j=1,NUM_LIQ)
        END DO
      CLOSE( 16 )

c      k = 1     !!! angle
c      j = 1500  !!! wavenumber
c      i = 100   !!! OD * SZ
c      print *,'ice ',k,j,i,THETMAX,NUM_ICE,RT_ICE(k,j,i)
c      print *,'liq ',k,j,i,THETMAX,NUM_LIQ,RT_LIQ(k,j,i)
c      print *,(RT_ICE(k,j,i),k=1,2*THETMAX)
c      print *,(RT_LIQ(k,j,i),k=1,2*THETMAX)
c      print *,(RT_ICE(k,j,i),j=1,NUM_ICE)
c      print *,(RT_LIQ(k,j,i),j=1,NUM_LIQ)
c      print *,(RT_ICE(k,j,i),j=1,25)
c      print *,(RT_LIQ(k,j,i),j=1,25)

      RETURN

 98   PRINT *, ' ** READ ERROR in load_rt_tables (ICE) ** ',
     $           i,TAUMAX_ICE * DMAX_ICE
 99   PRINT *, ' ** READ ERROR in load_rt_tables (LIQ) ** ',
     $           i,TAUMAX_LIQ * DMAX_LIQ

      END

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE CAT12(f1,f2,fx)

      character*128 f1,f2
      character*256 fx

      integer i,i1,i2

      i1 = 128
 10   CONTINUE
      IF (f1(i1:i1) .EQ. ' ') THEN
        i1 = i1 - 1
        GOTO 10
        END IF

      i2 = 128
 20   CONTINUE
      IF (f2(i2:i2) .EQ. ' ') THEN
        i2 = i2 - 1
        GOTO 20
        END IF
      
      DO i = 1,256
        fx(i:i) = ' '
        END DO
      DO i = 1,i1
        fx(i:i) = f1(i:i)
        END DO
      DO i = 1,i2
        fx(i+i1:i+i1) = f2(i:i)
        END DO

      RETURN
      END 

c************************************************************************
