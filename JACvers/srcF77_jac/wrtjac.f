      SUBROUTINE WRTJAC_T(IOUNTZ,IPROF,NLAY,NCHAN,FREQ,RAD,JAC_OUTPUT_UNITS,JAC_ST_C,JAC_TZ_C)

      IMPLICIT NONE
      
      include 'incFTC.f'

      INTEGER IOUNTZ   !!! I/O unit ID
      INTEGER IPROF    !!! current profile number
      INTEGER NLAY     !!! number of layers to be written out
      INTEGER NCHAN    !!! number of chan to be written out
      REAL    JAC_ST_C(MXCHAN)        !!! ST jacobian
      REAL    JAC_TZ_C(MAXLAY,MXCHAN) !!! TZ jacobian

      REAL RAD(MXCHAN),FREQ(MXCHAN)
      INTEGER JAC_OUTPUT_UNITS

      INTEGER iC,iL
      REAL     raDeriv_Rad(MXCHAN),a(MXCHAN),b(MXCHAN),c(MXCHAN),d(MXCHAN)

      !! prof number, number of layers to write, number of chans to write,100 = TZ/ST JAC
      WRITE(IOUNTZ) IPROF,NLAY+1,NCHAN,100  

      IF (JAC_OUTPUT_UNITS .EQ. 0) THEN
        !! JAC_OUTPUT_UNITS = 0  ==> output drad/dT
        DO iL = 1,NLAY
          WRITE(IOUNTZ) (JAC_TZ_C(IL,iC),iC=1,NCHAN)
        END DO
        WRITE(IOUNTZ) (JAC_ST_C(iC),iC=1,NCHAN)
      ELSEIF (JAC_OUTPUT_UNITS .EQ. 1) THEN
        !! JAC_OUTPUT_UNITS = 1  ==> output dBT/dT
        a = c1*c2*(FREQ(1:NCHAN)**4)
        b = c1*(FREQ(1:NCHAN)**3)/RAD(1:NCHAN) + 1
        c = (log(b))**2
        d = (RAD(1:NCHAN))**2
        raDeriv_Rad(1:NCHAN) = a/(b*c*d)
        DO iL = 1,NLAY
          WRITE(IOUNTZ) (raDeriv_rad(IC) * JAC_TZ_C(IL,iC),iC=1,NCHAN)
        END DO
        WRITE(IOUNTZ) (raDeriv_rad(IC) * JAC_ST_C(iC),iC=1,NCHAN)      
      END IF

      RETURN 
      END

c************************************************************************
      SUBROUTINE WRTJAC_GAS(IOUNG1,IPROF,NLAY,NCHAN,iGASID,FREQ,RAD,JAC_OUTPUT_UNITS,GAMNT,JAC_G1_C)

      IMPLICIT NONE

      include 'incFTC.f'

      INTEGER iGasID   !!! gasID
      INTEGER IOUNG1   !!! I/O unit ID
      INTEGER IPROF    !!! current profile number
      INTEGER NLAY     !!! number of layers to be written out
      INTEGER NCHAN    !!! number of chan to be written out
      REAL    JAC_G1_C(MAXLAY,MXCHAN) !!! GN jacobian

      REAL GAMNT(MAXLAY)
      REAL RAD(MXCHAN),FREQ(MXCHAN)
      INTEGER JAC_OUTPUT_UNITS

      INTEGER iC,iL
      REAL     raDeriv_Rad(MXCHAN),a(MXCHAN),b(MXCHAN),c(MXCHAN),d(MXCHAN)

      !! prof number, number of layers to write, number of chans to write, what jac being done
      WRITE(IOUNG1) IPROF,NLAY,NCHAN,iGASID  
      IF ((IGASID .EQ. 200) .OR. (JAC_OUTPUT_UNITS .EQ. 0)) THEN
        !! IGASID == 200 ==> weighting function = no units
        !! JAC_OUTPUT_UNITS = 0  ==> output drad/dq
        DO iL = 1,NLAY
          WRITE(IOUNG1) (JAC_G1_C(IL,iC),iC=1,NCHAN)
        END DO
      ELSEIF ((IGASID .LE. 12) .AND. (JAC_OUTPUT_UNITS .EQ. 1)) THEN
        !! JAC_OUTPUT_UNITS = 1  ==> output q dBT/dq = dBT/dlog(q)) 
        a = c1*c2*(FREQ(1:NCHAN)**4)
        b = c1*(FREQ(1:NCHAN)**3)/RAD(1:NCHAN) + 1
        c = (log(b))**2
        d = (RAD(1:NCHAN))**2
        raDeriv_Rad(1:NCHAN) = a/(b*c*d)
        DO iL = 1,NLAY
          WRITE(IOUNG1) (raDeriv_rad(IC) * JAC_G1_C(IL,iC) * GAMNT(IL),iC=1,NCHAN)
        END DO
      END IF

      RETURN 
      END

c************************************************************************
