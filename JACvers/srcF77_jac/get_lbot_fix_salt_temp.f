      SUBROUTINE get_lbot_fix_salt_temp(
     $   NLAY, PLEV, PLAY, PSURF, LBOT, BLMULT,                 !!! for getbot
     $   TEMP, SALT, IPROF, AIRSLAY)

      IMPLICIT NONE
      
      include "incFTC.f"

c input
       REAL PLEV(MAXLAY+1)
       REAL PLAY(MAXLAY)   ! layer mean pressureo
       INTEGER NLAY        ! number of layers in profile
       REAL  PSURF         ! surface pressure
       INTEGER  PTYPE      ! profile type
       REAL   TEMP(MAXLAY) ! prof layer average temperature
       REAL    SALT        ! input satellite altitude (kilometers)
       INTEGER IPROF       ! profile loop
       INTEGER AIRSLAY
c output
       INTEGER   LBOT             ! bottom layer index number
       REAL BLMULT                ! bottom layer fractional multiplier

c local
       INTEGER I
       REAL    RJUNK1,RJUNK2
C      for MEAN_T
       REAL TPSEUD(MAXLAY)

C      -------------------------------------
C      Determine bottom layer, CO2, & angles
C      -------------------------------------
       CALL GETBOT(NLAY, PLEV, PSURF, LBOT, BLMULT)

C      Calc the fractional bottom layer air temperature
ccc
c       TEMP(LBOT)=TEMP(LBOT-1) + BLMULT*( TEMP(LBOT) - TEMP(LBOT-1) )
c Above line commented out & replaced by Scott Hannon, 24 July 2003.
c Mistakenly treats T at the center of the layer above as T at the
c bottom of the layer above.
cc
C
C
       IF (PTYPE .EQ. AIRSLAY) THEN
C         Copy pseudo level temperatures to another array
          DO I=1,LBOT
             TPSEUD(I)=TEMP(I)
          ENDDO
C         Convert temperatures
          CALL MEAN_T(LBOT, PLEV, PSURF, TPSEUD, TEMP)
C
       ELSE
C         Calc mean pressure for bottom fractional layer
          RJUNK1 = ( PSURF - PLEV(LBOT) )/LOG( PSURF/PLEV(LBOT) )
C         Do interpolation for fractional bottom layer mean temperature
C         assuming T is in linear in log(P)
          RJUNK2=( TEMP(LBOT) - TEMP(LBOT-1) )/
     $       LOG( PLAY(LBOT)/PLAY(LBOT-1) )             ! slope
          TEMP(LBOT)=RJUNK2*LOG( RJUNK1/PLAY(LBOT-1) ) + TEMP(LBOT - 1)
       ENDIF
C
C      Check satellite elevation
       IF (SALT .GT. 0.0) THEN
C         Warn and use default if invalid
          IF (SALT .LT. XSALT-50 .OR. SALT .GT. XSALT+50) THEN
             WRITE(IOINFO,1020) IPROF, SALT, XSALT
 1020        FORMAT('Warning! Profile',I5,
     $          ': replacing invalid input satellite altitude ',
     $          1PE11.4,' with default ',1PE11.4,' km')
             SALT=XSALT
          ENDIF
       ELSE
          SALT=XSALT
       ENDIF

      RETURN
      END
