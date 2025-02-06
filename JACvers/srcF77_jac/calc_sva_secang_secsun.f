      SUBROUTINE CALC_SVA_SECANG_SECSUN(SALT,SATZEN,SATANG,SUNANG,LBOT,ALT,
     $      IPROF,LSTCHN,NCHNTE,CLISTN,FREQ,INDCHN,
     $      SVA,SECANG,SECSUN,SUNFDG,SUNCOS,SCOS1,DOSUN,QUICKINDNTE)

       IMPLICIT NONE

       include 'incFTC.f'

c input
C      for satellite viewing angle
       REAL    SATANG      ! input satellite scan angle (degrees)
       REAL    SATZEN      ! input satellite zenith angle (degrees)
       REAL    SALT        ! input satellite altitude (kilometers)
       REAL    ALT(MAXLAY) ! prof layer altitudes
       REAL SUNANG         ! solar zenith angle (at 0 altitude)
       INTEGER   LBOT      ! bottom layer index number
       INTEGER  IPROF      ! profile loop counter

       INTEGER LSTCHN(MXCHAN)  ! list of selected channels
       INTEGER INDCHN(MXCHAN)  ! array indices for all channels
       REAL   FREQ(MXCHAN)    ! chan center frequency
       INTEGER NCHNTE                    ! number of non-LTE channels
       INTEGER CLISTN(MXCNTE)            ! non-LTE channel list

c output
       REAL    SVA         ! satellite viewing angle (degrees)
       REAL SECANG(MAXLAY)        ! local path angle secant
       REAL SUNFDG         ! fudge factor for large solar angles
       REAL SECSUN(MAXLAY) ! secant of effective sun local path angle
       REAL SUNCOS         ! cosine of sun zenith angle
       LOGICAL DOSUN       ! do sun calc?
       INTEGER QUICKINDNTE(MXCHAN)       ! list of non-LTE channels
       REAL SCOS1          ! cosine of sun zenith angle at layer1

c      local
       INTEGER I,III,L
       REAL   CONV         ! degrees to radians conversion factor
       REAL ANGMAX         ! maximum allowed viewing angle
       REAL    EVA         ! (Earth) local viewing angle
       REAL SZALAY         ! solar zenith angle in some layer
       REAL RJUNK2,RJUNK1

C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
       REAL VACONV
       REAL SACONV
       INTEGER INTERSECT

c************************************************************************

C      CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02

C
C      Convert SATZEN or SATANG to viewing angle
       IF (SATZEN .GE. 0.0 .AND. SATZEN .LT. 63.0) THEN
C         Convert zenith angle at surface to view angle at satellite
          SVA=SACONV( SATZEN, SALT*1000.0 )/CONV
       ELSE
C         Check if scan angle is valid
          IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
C            View angle should be within a few degrees of scan angle
             SVA=ABS( SATANG )
          ELSE
             WRITE(IOERR,1030) IPROF, SATZEN, SATANG
 1030        FORMAT('Error! Profile',I5,
     $          ': invalid angles for SATZEN ',1PE11.4,
     $          ' and SATANG ',E11.4) 
             STOP
          ENDIF
       ENDIF

       ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
       IF (SVA .GT. ANGMAX) THEN
C         Truncate angle if too big
          WRITE(IOINFO,1040) IPROF, SVA
 1040     FORMAT('Warning! Profile',I5,': truncating view angle ',
     $       1PE11.4,' to 53 degrees')
          SVA=ANGMAX
       ENDIF

C      Convert from satellite to earth viewing angle (in radians)
       DO L=1,LBOT
         EVA=VACONV(SVA, SALT, ALT(L))
         SECANG(L)=1.0E+0/COS(EVA)
ccccccccccccc
c        for testing
c        SECANG(L)=SVA
ccccccccccccc
       ENDDO

C      Calc total sun angle secant
       DOSUN=.FALSE.
       SECSUN = 1.0
       SUNFDG = 1.0
       IF (SUNANG .GE. 0.0 .AND. SUNANG .LT. 89.9) DOSUN=.TRUE.
       IF (DOSUN) THEN
          SUNCOS=COS(CONV*SUNANG)
          SZALAY=SACONV(SUNANG,ALT(1))
          SCOS1=COS(SZALAY)
          RJUNK2=SECANG(LBOT) + 1.0/SUNCOS ! Total secant

C         Calc non-unity fudge factor if total secant > 9
          IF (RJUNK2 .GT. 9.0) THEN
C            fudge factor = true_total_secant/calc_total_secant
             SUNFDG=RJUNK2/9.0
C            truncated solar angle to use to calc SECSUN
             RJUNK1=ACOS( 1.0/(9.0 - SECANG(LBOT)) )/CONV
          ELSE
             SUNFDG=1.0
             RJUNK1=SUNANG
          ENDIF
c Should I change SUNFDG to SUNFDG(MAXLAY)?
C
          DO L=1,LBOT
             SZALAY=SACONV(RJUNK1,ALT(L))
             SECSUN(L)=SECANG(L) + 1.0E+0/COS(SZALAY)
          ENDDO
 
          QUICKINDNTE = 0
          !! LSTCHN = h.ichan
          !! so eg if I = 1520, h.ichan(1520) = listchn(1520) = 1291; h.vchan(1520) = freq(1520) = 1231.3 cm-1
          !! so eg if I = 2600, h.ichan(2600) = listchn(2600) = 2333; h.vchan(2600) = freq(2600) = 2616.4 cm-1
          !! so eg if I = 2371, h.ichan(2371) = listchn(2371) = 2100; h.vchan(2371) = freq(2371) = 2379.4 cm-1

          IF (DEBUG) THEN
            do I = 1,MXCHAN
              if (LSTCHN(I) .GT. 0) THEN
                write(*,'(A,I5,I5,F12.5,F12.5)') 'NTE CHECK(1) : I,LSTCHN(I)=h.ichan,FREQ(I) = ',
     $                I,LSTCHN(I),FREQ(I),FREQ(INDCHN(LSTCHN(I)))
              end if
            end do
            print *,'  '
  
            do I = 1,NCHNTE
              write(*,'(A,I5,I5,I5,F12.4)') 
     $           'NTE CHECK(2) : I,CLISTN(I)=h.ichan --> LSTCHN(I),FREQ(CLISTN(I))',
     $           I,CLISTN(I),INDCHN(CLISTN(I)),FREQ(INDCHN(CLISTN(I)))
            end do
            print *,'  '
          END IF
  
          DO I = 1,NCHNTE
            QUICKINDNTE(INDCHN(CLISTN(I))) = I
c            print *,'MIAOW MIAOW',I,III,INDCHN(CLISTN(I)),QUICKINDNTE(INDCHN(CLISTN(I)))
          END DO
       ENDIF

       RETURN
       END
