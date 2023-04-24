      SUBROUTINE WRTJAC_T(IOUNTZ,IPROF,NLAY,NUMCHAN,JAC_ST_C,JAC_TZ_C)

      IMPLICIT NONE
      
      include 'incFTC.f'

      INTEGER IOUNTZ   !!! I/O unit ID
      INTEGER IPROF    !!! current profile number
      INTEGER NLAY     !!! number of layers to be written out
      INTEGER NUMCHAN  !!! number of chan to be written out
      REAL    JAC_ST_C(MXCHAN)        !!! ST jacobian
      REAL    JAC_TZ_C(MXCHAN,MAXLAY) !!! TZ jacobian

      INTEGER iC,iL

      WRITE(IOUNTZ) IPROF,NLAY+1,NUMCHAN  !! prof number, number of layers to write, number of chans to write
      DO iL = 1,NLAY
        WRITE(IOUNTZ) (JAC_TZ_C(iC,IL),iC=1,NUMCHAN)
      END DO
      WRITE(IOUNTZ) (JAC_ST_C(iC),iC=1,NUMCHAN)

      RETURN 
      END

c************************************************************************
      SUBROUTINE WRTJAC_GAS(IOUNG1,IPROF,NLAY,NUMCHAN,iGASID,JAC_G1_C)

      IMPLICIT NONE

      include 'incFTC.f'

      INTEGER iGasID   !!! gasID
      INTEGER IOUNG1  !!! I/O unit ID
      INTEGER IPROF    !!! current profile number
      INTEGER NLAY     !!! number of layers to be written out
      INTEGER NUMCHAN  !!! number of chan to be written out
      REAL    JAC_G1_C(MXCHAN,MAXLAY) !!! GN jacobian

      INTEGER iC,iL

      WRITE(IOUNG1) IPROF,NLAY,NUMCHAN  !! prof number, number of layers to write, number of chans to write
      DO iL = 1,NLAY
        WRITE(IOUNG1) (JAC_G1_C(iC,IL),iC=1,NUMCHAN)
      END DO

      RETURN 
      END

c************************************************************************
