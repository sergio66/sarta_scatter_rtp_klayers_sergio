c testing
c       FCLEAR = 0.0; CFRA1X = 1.0; CFRA2X = 0.0; CFRA12 = 0.0
c       FCLEAR = 0.0; CFRA1X = 0.0; CFRA2X = 1.0; CFRA12 = 0.0
c       FCLEAR = 0.0; CFRA1X = 0.0; CFRA2X = 0.0; CFRA12 = 1.0
c       FCLEAR = 1.0; CFRA1X = 0.0; CFRA2X = 0.0; CFRA12 = 0.0

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(1,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_G1_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * dTAU_DG1(1:NLAY,1:NCHAN)
        JAC_G1_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * dTAU_DG1(1:NLAY,1:NCHAN)
        JAC_G1_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * dTAU_DG1(1:NLAY,1:NCHAN)
        JAC_G1_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * dTAU_DG1(1:NLAY,1:NCHAN)

        JAC_G1_C = FCLEAR*JAC_G1_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
     $             CFRA1X*JAC_G1_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
     $             CFRA2X*JAC_G1_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
     $             CFRA12*JAC_G1_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)
        CALL WRTJAC_GAS(IOUNG1,IPROF,NLAY,NCHAN,1,FREQ,RAD,JAC_OUTPUT_UNITS,WAMNT,JAC_G1_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(2,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_G2_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * dTAU_DG2(1:NLAY,1:NCHAN)
        JAC_G2_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * dTAU_DG2(1:NLAY,1:NCHAN)
        JAC_G2_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * dTAU_DG2(1:NLAY,1:NCHAN)
        JAC_G2_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * dTAU_DG2(1:NLAY,1:NCHAN)

        JAC_G2_C = FCLEAR*JAC_G2_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
     $             CFRA1X*JAC_G2_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
     $             CFRA2X*JAC_G2_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
     $             CFRA12*JAC_G2_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)

        CALL WRTJAC_GAS(IOUNG2,IPROF,NLAY,NCHAN,2,FREQ,RAD,JAC_OUTPUT_UNITS,FAMNT,JAC_G2_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(3,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_G3_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * dTAU_DG3(1:NLAY,1:NCHAN)
        JAC_G3_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * dTAU_DG3(1:NLAY,1:NCHAN)
        JAC_G3_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * dTAU_DG3(1:NLAY,1:NCHAN)
        JAC_G3_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * dTAU_DG3(1:NLAY,1:NCHAN)

        JAC_G3_C = FCLEAR*JAC_G3_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
     $             CFRA1X*JAC_G3_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
     $             CFRA2X*JAC_G3_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
     $             CFRA12*JAC_G3_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)

        CALL WRTJAC_GAS(IOUNG3,IPROF,NLAY,NCHAN,3,FREQ,RAD,JAC_OUTPUT_UNITS,OAMNT,JAC_G3_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(4,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_G4_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * dTAU_DG4(1:NLAY,1:NCHAN)
        JAC_G4_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * dTAU_DG4(1:NLAY,1:NCHAN)
        JAC_G4_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * dTAU_DG4(1:NLAY,1:NCHAN)
        JAC_G4_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * dTAU_DG4(1:NLAY,1:NCHAN)

        JAC_G4_C = FCLEAR*JAC_G4_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
     $             CFRA1X*JAC_G4_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
     $             CFRA2X*JAC_G4_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
     $             CFRA12*JAC_G4_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)

        CALL WRTJAC_GAS(IOUNG4,IPROF,NLAY,NCHAN,4,FREQ,RAD,JAC_OUTPUT_UNITS,NAMNT,JAC_G4_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(5,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_G5_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * dTAU_DG5(1:NLAY,1:NCHAN)
        JAC_G5_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * dTAU_DG5(1:NLAY,1:NCHAN)
        JAC_G5_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * dTAU_DG5(1:NLAY,1:NCHAN)
        JAC_G5_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * dTAU_DG5(1:NLAY,1:NCHAN)

        JAC_G5_C = FCLEAR*JAC_G5_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
     $             CFRA1X*JAC_G5_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
     $             CFRA2X*JAC_G5_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
     $             CFRA12*JAC_G5_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)

        CALL WRTJAC_GAS(IOUNG5,IPROF,NLAY,NCHAN,5,FREQ,RAD,JAC_OUTPUT_UNITS,CAMNT,JAC_G5_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(6,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_G6_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * dTAU_DG6(1:NLAY,1:NCHAN)
        JAC_G6_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * dTAU_DG6(1:NLAY,1:NCHAN)
        JAC_G6_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * dTAU_DG6(1:NLAY,1:NCHAN)
        JAC_G6_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * dTAU_DG6(1:NLAY,1:NCHAN)

        JAC_G6_C = FCLEAR*JAC_G6_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
     $             CFRA1X*JAC_G6_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
     $             CFRA2X*JAC_G6_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
     $             CFRA12*JAC_G6_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)

        CALL WRTJAC_GAS(IOUNG6,IPROF,NLAY,NCHAN,6,FREQ,RAD,JAC_OUTPUT_UNITS,MAMNT,JAC_G6_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(9,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_G9_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * dTAU_DG9(1:NLAY,1:NCHAN)
        JAC_G9_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * dTAU_DG9(1:NLAY,1:NCHAN)
        JAC_G9_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * dTAU_DG9(1:NLAY,1:NCHAN)
        JAC_G9_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * dTAU_DG9(1:NLAY,1:NCHAN)

        JAC_G9_C = FCLEAR*JAC_G9_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
     $             CFRA1X*JAC_G9_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
     $             CFRA2X*JAC_G9_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
     $             CFRA12*JAC_G9_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)

        CALL WRTJAC_GAS(IOUNG9,IPROF,NLAY,NCHAN,9,FREQ,RAD,JAC_OUTPUT_UNITS,SAMNT,JAC_G9_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
c      IF (INTERSECT(11,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
c        JAC_G11_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * dTAU_DG11(1:NLAY,1:NCHAN)
c        JAC_G11_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * dTAU_DG11(1:NLAY,1:NCHAN)
c        JAC_G11_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * dTAU_DG11(1:NLAY,1:NCHAN)
c        JAC_G11_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * dTAU_DG11(1:NLAY,1:NCHAN)

c        JAC_G11_C = FCLEAR*JAC_G11_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
c     $              CFRA1X*JAC_G11_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
c     $              CFRA2X*JAC_G11_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
c     $              CFRA11*JAC_G11_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)

c        CALL WRTJAC_GAS(IOUNG11,IPROF,NLAY,NCHAN,11,FREQ,RAD,JAC_OUTPUT_UNITS,AAMNT,JAC_G11_C)
c      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(12,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_G12_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * dTAU_DG12(1:NLAY,1:NCHAN)
        JAC_G12_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * dTAU_DG12(1:NLAY,1:NCHAN)
        JAC_G12_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * dTAU_DG12(1:NLAY,1:NCHAN)
        JAC_G12_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * dTAU_DG12(1:NLAY,1:NCHAN)

        JAC_G12_C = FCLEAR*JAC_G12_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
     $              CFRA1X*JAC_G12_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
     $              CFRA2X*JAC_G12_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
     $              CFRA12*JAC_G12_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)

        CALL WRTJAC_GAS(IOUNG12,IPROF,NLAY,NCHAN,12,FREQ,RAD,JAC_OUTPUT_UNITS,HAMNT,JAC_G12_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (INTERSECT(100,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_ST_C(1:NCHAN)  = EMIS(1:NCHAN) * DBTDT(NLAY+1,1:NCHAN)
        JAC_ST_1(1:NCHAN)  = EMIS(1:NCHAN) * DBTDT(NLAY+1,1:NCHAN)
        JAC_ST_2(1:NCHAN)  = EMIS(1:NCHAN) * DBTDT(NLAY+1,1:NCHAN)
        JAC_ST_12(1:NCHAN) = EMIS(1:NCHAN) * DBTDT(NLAY+1,1:NCHAN)
        JAC_ST_C = FCLEAR*JAC_ST_C(1:NCHAN)  * L2S4(1,NLAY,1:NCHAN) + 
     $             CFRA1X*JAC_ST_1(1:NCHAN)  * L2S4(2,NLAY,1:NCHAN) +  
     $             CFRA2X*JAC_ST_2(1:NCHAN)  * L2S4(3,NLAY,1:NCHAN) + 
     $             CFRA12*JAC_ST_12(1:NCHAN) * L2S4(4,NLAY,1:NCHAN)

        JAC_TZ_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * DTAU_DTZ(1:NLAY,1:NCHAN)
        JAC_TZ_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * DTAU_DTZ(1:NLAY,1:NCHAN)
        JAC_TZ_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * DTAU_DTZ(1:NLAY,1:NCHAN)
        JAC_TZ_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * DTAU_DTZ(1:NLAY,1:NCHAN)

        JAC_TZ_C(1:NLAY,1:NCHAN)  = JAC_TZ_C(1:NLAY,1:NCHAN)  + DBTDT(1:NLAY,1:NCHAN) * (1 - EXP(-TAU4(1,1:NLAY,1:NCHAN)))
        JAC_TZ_1(1:NLAY,1:NCHAN)  = JAC_TZ_1(1:NLAY,1:NCHAN)  + DBTDT(1:NLAY,1:NCHAN) * (1 - EXP(-TAU4(2,1:NLAY,1:NCHAN)))
        JAC_TZ_2(1:NLAY,1:NCHAN)  = JAC_TZ_2(1:NLAY,1:NCHAN)  + DBTDT(1:NLAY,1:NCHAN) * (1 - EXP(-TAU4(3,1:NLAY,1:NCHAN)))
        JAC_TZ_12(1:NLAY,1:NCHAN) = JAC_TZ_12(1:NLAY,1:NCHAN) + DBTDT(1:NLAY,1:NCHAN) * (1 - EXP(-TAU4(4,1:NLAY,1:NCHAN)))

        JAC_TZ_C = FCLEAR*JAC_TZ_C(1:NLAY,1:NCHAN)  * L2S4(1,1:NLAY,1:NCHAN) + 
     $             CFRA1X*JAC_TZ_1(1:NLAY,1:NCHAN)  * L2S4(2,1:NLAY,1:NCHAN) + 
     $             CFRA2X*JAC_TZ_2(1:NLAY,1:NCHAN)  * L2S4(3,1:NLAY,1:NCHAN) + 
     $             CFRA12*JAC_TZ_12(1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)

c        DO IJUNK = 1,100
c          write(*,'(A,I5,10(ES15.6))') 'TADA ',IJUNK,FCLEAR,CFRA1X,CFRA2X,CFRA12,RAD4(1,IJUNK,254),DTAU_DTZ(IJUNK,254),
c                                 DBTDT(IJUNK,254),TAU4(1,IJUNK,254),L2S4(1,IJUNK,254),JAC_TZ_C(IJUNK,254)
c        END DO

        CALL WRTJAC_T(IOUNTZ,IPROF,NLAY,NCHAN,FREQ,RAD,JAC_OUTPUT_UNITS,JAC_ST_C,JAC_TZ_C)
      END IF   

c%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (INTERSECT(200,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        JAC_WGT_C  = WGT4(1,1:NLAY,1:NCHAN) * L2S4(1,1:NLAY,1:NCHAN)
        JAC_WGT_1  = WGT4(2,1:NLAY,1:NCHAN) * L2S4(2,1:NLAY,1:NCHAN)
        JAC_WGT_2  = WGT4(3,1:NLAY,1:NCHAN) * L2S4(3,1:NLAY,1:NCHAN)
        JAC_WGT_12 = WGT4(4,1:NLAY,1:NCHAN) * L2S4(4,1:NLAY,1:NCHAN)
        JAC_WGT_C = FCLEAR*JAC_WGT_C + CFRA1X*JAC_WGT_1 + CFRA2X*JAC_WGT_2 + CFRA12*JAC_WGT_12
        CALL WRTJAC_GAS(IOUNWGT,IPROF,NLAY,NCHAN,200,FREQ,RAD,JAC_OUTPUT_UNITS,WAMNT,JAC_WGT_C)
      END IF   
c%%%%%%%%%%%%%%%%%%%%%%%%%