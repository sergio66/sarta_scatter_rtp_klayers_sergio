      
      IF (INTERSECT(100,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_ST_C(1:NCHAN)  = EMIS(1:NCHAN) * DBTDT(NLAY+1,1:NCHAN)
        JAC_ST_1(1:NCHAN)  = EMIS(1:NCHAN) * DBTDT(NLAY+1,1:NCHAN)
        JAC_ST_2(1:NCHAN)  = EMIS(1:NCHAN) * DBTDT(NLAY+1,1:NCHAN)
        JAC_ST_12(1:NCHAN) = EMIS(1:NCHAN) * DBTDT(NLAY+1,1:NCHAN)
!        !!!! this gives right answers compared to sarta finite diff + kcarta  so DO NOT CHANGE THIS     USE L2S4, USE L2S4, USE L2S4, USE L2S4, USE L2S4, USE L2S4
!        JAC_ST_C = FCLEAR*JAC_ST_C(1:NCHAN)  * L2S4(1,NLAY,1:NCHAN) + 
!     $             CFRA1X*JAC_ST_1(1:NCHAN)  * L2S4(2,NLAY,1:NCHAN) +  
!     $             CFRA2X*JAC_ST_2(1:NCHAN)  * L2S4(3,NLAY,1:NCHAN) + 
!     $             CFRA12*JAC_ST_12(1:NCHAN) * L2S4(4,NLAY,1:NCHAN)
        !!!! this gives right answers compared to sarta finite diff + kcarta so DO NOT CHANGE THIS     USE L2S4ABOVE(NLAY+1)
        JAC_ST_C(1:NCHAN) = FCLEAR*JAC_ST_C(1:NCHAN)  * L2S4above(1,NLAY+1,1:NCHAN) + 
     $                      CFRA1X*JAC_ST_1(1:NCHAN)  * L2S4above(2,NLAY+1,1:NCHAN) +  
     $                      CFRA2X*JAC_ST_2(1:NCHAN)  * L2S4above(3,NLAY+1,1:NCHAN) + 
     $                      CFRA12*JAC_ST_12(1:NCHAN) * L2S4above(4,NLAY+1,1:NCHAN)

        JAC_TZ_C(1:NLAY,1:NCHAN)  = PLANCK_RAD4(1,1:NLAY,1:NCHAN) * DTAU_DTZ(1:NLAY,1:NCHAN)
        JAC_TZ_1(1:NLAY,1:NCHAN)  = PLANCK_RAD4(2,1:NLAY,1:NCHAN) * DTAU_DTZ(1:NLAY,1:NCHAN)
        JAC_TZ_2(1:NLAY,1:NCHAN)  = PLANCK_RAD4(3,1:NLAY,1:NCHAN) * DTAU_DTZ(1:NLAY,1:NCHAN)
        JAC_TZ_12(1:NLAY,1:NCHAN) = PLANCK_RAD4(4,1:NLAY,1:NCHAN) * DTAU_DTZ(1:NLAY,1:NCHAN)

        JAC_TZ_C(1:NLAY,1:NCHAN)  = JAC_TZ_C(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * dTAU_DTZ(1:NLAY,1:NCHAN)
        JAC_TZ_1(1:NLAY,1:NCHAN)  = JAC_TZ_1(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(2,1:NLAY,1:NCHAN) * dTAU_DTZ(1:NLAY,1:NCHAN)
        JAC_TZ_2(1:NLAY,1:NCHAN)  = JAC_TZ_2(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(3,1:NLAY,1:NCHAN) * dTAU_DTZ(1:NLAY,1:NCHAN)
        JAC_TZ_12(1:NLAY,1:NCHAN) = JAC_TZ_12(1:NLAY,1:NCHAN) - RTHERM4_SOLAR4(4,1:NLAY,1:NCHAN) * dTAU_DTZ(1:NLAY,1:NCHAN)

        JAC_TZ_C(1:NLAY,1:NCHAN)  = JAC_TZ_C(1:NLAY,1:NCHAN)  + DBTDT(1:NLAY,1:NCHAN) * (1 - EXP(-TAU4(1,1:NLAY,1:NCHAN)))
        JAC_TZ_1(1:NLAY,1:NCHAN)  = JAC_TZ_1(1:NLAY,1:NCHAN)  + DBTDT(1:NLAY,1:NCHAN) * (1 - EXP(-TAU4(2,1:NLAY,1:NCHAN)))
        JAC_TZ_2(1:NLAY,1:NCHAN)  = JAC_TZ_2(1:NLAY,1:NCHAN)  + DBTDT(1:NLAY,1:NCHAN) * (1 - EXP(-TAU4(3,1:NLAY,1:NCHAN)))
        JAC_TZ_12(1:NLAY,1:NCHAN) = JAC_TZ_12(1:NLAY,1:NCHAN) + DBTDT(1:NLAY,1:NCHAN) * (1 - EXP(-TAU4(4,1:NLAY,1:NCHAN)))

        JAC_TZ_C(1:NLAY,1:NCHAN) = FCLEAR*JAC_TZ_C(1:NLAY,1:NCHAN)  * L2S4above(1,1:NLAY,1:NCHAN) + 
     $                             CFRA1X*JAC_TZ_1(1:NLAY,1:NCHAN)  * L2S4above(2,1:NLAY,1:NCHAN) + 
     $                             CFRA2X*JAC_TZ_2(1:NLAY,1:NCHAN)  * L2S4above(3,1:NLAY,1:NCHAN) + 
     $                             CFRA12*JAC_TZ_12(1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN)

        JAC_TZ_C(1:5,1:NCHAN) = JAC_TZ_C(1:5,1:NCHAN) + NLTEJACPRED5T(1:5,1:NCHAN) * L2S4above(1,1:5,1:NCHAN)

c        DO IJUNK = 1,100
c          write(*,'(A,I5,10(ES15.6))') 'TADA ',IJUNK,FCLEAR,CFRA1X,CFRA2X,CFRA12,PLANCK_RAD4(1,IJUNK,254),DTAU_DTZ(IJUNK,254),
c                                 DBTDT(IJUNK,254),TAU4(1,IJUNK,254),L2S4(1,IJUNK,254),JAC_TZ_C(IJUNK,254)
c        END DO

        CALL WRTJAC_T(IOUNTZ,IPROF,NLAY,NCHAN,FREQ,RAD,JAC_OUTPUT_UNITS,JAC_ST_C,JAC_TZ_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(200,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_WGT_C(1:NLAY,1:NCHAN)  = WGT4(1,1:NLAY,1:NCHAN) * L2S4above(1,1:NLAY,1:NCHAN)
        JAC_WGT_1(1:NLAY,1:NCHAN)  = WGT4(2,1:NLAY,1:NCHAN) * L2S4above(2,1:NLAY,1:NCHAN)
        JAC_WGT_2(1:NLAY,1:NCHAN)  = WGT4(3,1:NLAY,1:NCHAN) * L2S4above(3,1:NLAY,1:NCHAN)
        JAC_WGT_12(1:NLAY,1:NCHAN) = WGT4(4,1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN)

        JAC_WGT_C(1:NLAY,1:NCHAN) = FCLEAR*JAC_WGT_C(1:NLAY,1:NCHAN) + 
     $                              CFRA1X*JAC_WGT_1(1:NLAY,1:NCHAN) + CFRA2X*JAC_WGT_2(1:NLAY,1:NCHAN) + 
     $                              CFRA12*JAC_WGT_12(1:NLAY,1:NCHAN)
        CALL WRTJAC_GAS(IOUNWGT,IPROF,NLAY,NCHAN,200,FREQ,RAD,JAC_OUTPUT_UNITS,WAMNT,JAC_WGT_C)
      END IF   
