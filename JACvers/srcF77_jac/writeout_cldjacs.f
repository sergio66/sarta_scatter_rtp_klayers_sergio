      IF (INTERSECT(300,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
        
        !!! would be cfrac1,cfrac2,cfrac12,amt1,amt2,sze1,sze2

        JAC_CLD_OUT = 0  
        
        !! ** <> ** <> ** <> ** <> ** !!
        !! figure out cloud top and cloud bottom layers
        IFOUND1 = -1
        IEFFCLD_TOP1 = 0
        DO WHILE ((IFOUND1 .LT. 0) .AND. (IEFFCLD_TOP1 .LT. LBOT))
          IEFFCLD_TOP1 = IEFFCLD_TOP1 + 1
          IF (CFRCL1(IEFFCLD_TOP1) .GE. 1e-5) IFOUND1 = +1
        END DO
        IF (IFOUND1 .GT. 0) THEN 
          IFOUND2 = -1
          IEFFCLD_BOT1 = LBOT+1
          DO WHILE ((IFOUND2 .LT. 0) .AND. (IEFFCLD_BOT1 .GT. 1))
            IEFFCLD_BOT1 = IEFFCLD_BOT1 - 1
            IF (CFRCL1(IEFFCLD_BOT1) .GE. 1e-5) IFOUND2 = +1
          END DO
          print *,'found CLD1 between layers ',IEFFCLD_TOP1,IEFFCLD_BOT1
        ELSE
          IEFFCLD_TOP1 = -1
          IEFFCLD_BOT1 = -1
        END IF

        IFOUND1 = -1
        IEFFCLD_TOP2 = 0
        DO WHILE ((IFOUND1 .LT. 0) .AND. (IEFFCLD_TOP2 .LT. LBOT))
          IEFFCLD_TOP2 = IEFFCLD_TOP2 + 1
          IF (CFRCL2(IEFFCLD_TOP2) .GE. 1e-5) IFOUND1 = +1
        END DO
        IF (IFOUND1 .GT. 0) THEN 
          IFOUND2 = -1
          IEFFCLD_BOT2 = LBOT+1
          DO WHILE ((IFOUND2 .LT. 0) .AND. (IEFFCLD_BOT2 .GT. 1))
            IEFFCLD_BOT2 = IEFFCLD_BOT2 - 1
            IF (CFRCL2(IEFFCLD_BOT2) .GE. 1e-5) IFOUND2 = +1
          END DO
          print *,'found CLD2 between layers ',IEFFCLD_TOP2,IEFFCLD_BOT2
        ELSE
          IEFFCLD_TOP2 = -1
          IEFFCLD_BOT2 = -1
        END IF 

        !! ** <> ** <> ** <> ** <> ** !!
        !! cld frac jacs cfrac1,cfrac2, cfrac12
        JAC_CLD_OUT(1,1:NCHAN) = -RAD4(1,1,1:NCHAN) + RAD4(2,1,1:NCHAN)  !!! dr/d cfrac1 = -rclr + r1
        JAC_CLD_OUT(2,1:NCHAN) = -RAD4(1,1,1:NCHAN) + RAD4(3,1,1:NCHAN)  !!! dr/d cfrac1 = -rclr + r2
        JAC_CLD_OUT(3,1:NCHAN) = +RAD4(1,1,1:NCHAN) + RAD4(4,1,1:NCHAN)  !!! dr/d cfrac12 = rclr + r12 - r1 - r2
     $                          - RAD4(2,1,1:NCHAN) - RAD4(3,1,1:NCHAN)  

        !! ** <> ** <> ** <> ** <> ** !!
        !! cld amt jac, cld 1
        RAMU = 0
        JAC_CLD_C  = 0  !! no clouds
        JAC_CLD_1  = 0
        JAC_CLD_2  = 0
        JAC_CLD_12 = 0
        IF ((IEFFCLD_BOT1 .GT. 0)    .AND. (IEFFCLD_TOP1 .GT. 0) .AND. 
     $      (IEFFCLD_BOT1 .LE. LBOT) .AND. (IEFFCLD_TOP1 .LE. LBOT)) THEN     
            
          DO IJUNK = IEFFCLD_TOP1,IEFFCLD_BOT1
             RAMU(IJUNK,1:NCHAN) = JACA_FINAL_1(1:NCHAN)
          END DO          
          JAC_CLD_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = JAC_CLD_C(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = JAC_CLD_1(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = JAC_CLD_2(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = JAC_CLD_12(1:NLAY,1:NCHAN) - RTHERM4_SOLAR4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = FCLEAR*JAC_CLD_C(1:NLAY,1:NCHAN)  * L2S4above(1,1:NLAY,1:NCHAN)
          JAC_CLD_1(1:NLAY,1:NCHAN)  = CFRA1X*JAC_CLD_1(1:NLAY,1:NCHAN)  * L2S4above(2,1:NLAY,1:NCHAN) 
          JAC_CLD_2(1:NLAY,1:NCHAN)  = CFRA2X*JAC_CLD_2(1:NLAY,1:NCHAN)  * L2S4above(3,1:NLAY,1:NCHAN) 
          JAC_CLD_12(1:NLAY,1:NCHAN) = CFRA12*JAC_CLD_12(1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN)                      

          JAC_CLD_C(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) = JAC_CLD_C(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_1(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_2(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_12(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN)

          JAC_CLD_OUT(4,1:NCHAN) = 0
          DO IJUNK = IEFFCLD_TOP1,IEFFCLD_BOT1
            JAC_CLD_OUT(4,1:NCHAN) = JAC_CLD_OUT(4,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN)
          END DO
        END IF 

        !! ** <> ** <> ** <> ** <> ** !!
        !! cld amt jac, cld 2
        RAMU = 0
        JAC_CLD_C  = 0  !! no clouds
        JAC_CLD_1  = 0
        JAC_CLD_2  = 0
        JAC_CLD_12 = 0
        IF ((IEFFCLD_BOT2 .GT. 0)    .AND. (IEFFCLD_TOP2 .GT. 0) .AND. 
     $      (IEFFCLD_BOT2 .LE. LBOT) .AND. (IEFFCLD_TOP2 .LE. LBOT)) THEN     
            
          DO IJUNK = IEFFCLD_TOP2,IEFFCLD_BOT2
             RAMU(IJUNK,1:NCHAN) = JACA_FINAL_2(1:NCHAN)
          END DO          
          JAC_CLD_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = JAC_CLD_C(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = JAC_CLD_1(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = JAC_CLD_2(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = JAC_CLD_12(1:NLAY,1:NCHAN) - RTHERM4_SOLAR4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = FCLEAR*JAC_CLD_C(1:NLAY,1:NCHAN)  * L2S4above(1,1:NLAY,1:NCHAN)
          JAC_CLD_1(1:NLAY,1:NCHAN)  = CFRA1X*JAC_CLD_1(1:NLAY,1:NCHAN)  * L2S4above(2,1:NLAY,1:NCHAN) 
          JAC_CLD_2(1:NLAY,1:NCHAN)  = CFRA2X*JAC_CLD_2(1:NLAY,1:NCHAN)  * L2S4above(3,1:NLAY,1:NCHAN) 
          JAC_CLD_12(1:NLAY,1:NCHAN) = CFRA12*JAC_CLD_12(1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN)                      

          JAC_CLD_C(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) = JAC_CLD_C(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_1(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_2(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_12(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN)

          JAC_CLD_OUT(5,1:NCHAN) = 0
          DO IJUNK = IEFFCLD_TOP2,IEFFCLD_BOT2
            JAC_CLD_OUT(5,1:NCHAN) = JAC_CLD_OUT(5,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN)
          END DO
        END IF 

        !! ** <> ** <> ** <> ** <> ** !!
        !! cld sze jac, cld 1
        RAMU = 0
        JAC_CLD_C  = 0  !! no clouds
        JAC_CLD_1  = 0
        JAC_CLD_2  = 0
        JAC_CLD_12 = 0
        IF ((IEFFCLD_BOT1 .GT. 0)    .AND. (IEFFCLD_TOP1 .GT. 0) .AND. 
     $      (IEFFCLD_BOT1 .LE. LBOT) .AND. (IEFFCLD_TOP1 .LE. LBOT)) THEN     
            
          DO IJUNK = IEFFCLD_TOP1,IEFFCLD_BOT1
             RAMU(IJUNK,1:NCHAN) = JACS_FINAL_1(1:NCHAN)
          END DO          
          JAC_CLD_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = JAC_CLD_C(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = JAC_CLD_1(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = JAC_CLD_2(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = JAC_CLD_12(1:NLAY,1:NCHAN) - RTHERM4_SOLAR4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = FCLEAR*JAC_CLD_C(1:NLAY,1:NCHAN)  * L2S4above(1,1:NLAY,1:NCHAN)
          JAC_CLD_1(1:NLAY,1:NCHAN)  = CFRA1X*JAC_CLD_1(1:NLAY,1:NCHAN)  * L2S4above(2,1:NLAY,1:NCHAN) 
          JAC_CLD_2(1:NLAY,1:NCHAN)  = CFRA2X*JAC_CLD_2(1:NLAY,1:NCHAN)  * L2S4above(3,1:NLAY,1:NCHAN) 
          JAC_CLD_12(1:NLAY,1:NCHAN) = CFRA12*JAC_CLD_12(1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN)                      

          JAC_CLD_C(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) = JAC_CLD_C(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_1(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_2(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_12(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN)

          JAC_CLD_OUT(6,1:NCHAN) = 0
          DO IJUNK = IEFFCLD_TOP1,IEFFCLD_BOT1
            JAC_CLD_OUT(6,1:NCHAN) = JAC_CLD_OUT(6,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN)
          END DO
        END IF 

        !! ** <> ** <> ** <> ** <> ** !!
        !! cld sze jac, cld 2
        RAMU = 0
        JAC_CLD_C  = 0  !! no clouds
        JAC_CLD_1  = 0
        JAC_CLD_2  = 0
        JAC_CLD_12 = 0
        IF ((IEFFCLD_BOT2 .GT. 0)    .AND. (IEFFCLD_TOP2 .GT. 0) .AND. 
     $      (IEFFCLD_BOT2 .LE. LBOT) .AND. (IEFFCLD_TOP2 .LE. LBOT)) THEN     
            
          DO IJUNK = IEFFCLD_TOP2,IEFFCLD_BOT2
             RAMU(IJUNK,1:NCHAN) = JACS_FINAL_2(1:NCHAN)
          END DO          
          JAC_CLD_C(1:NLAY,1:NCHAN)  = RAD4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = RAD4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = RAD4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = RAD4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = JAC_CLD_C(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = JAC_CLD_1(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = JAC_CLD_2(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = JAC_CLD_12(1:NLAY,1:NCHAN) - RTHERM4_SOLAR4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = FCLEAR*JAC_CLD_C(1:NLAY,1:NCHAN)  * L2S4above(1,1:NLAY,1:NCHAN)
          JAC_CLD_1(1:NLAY,1:NCHAN)  = CFRA1X*JAC_CLD_1(1:NLAY,1:NCHAN)  * L2S4above(2,1:NLAY,1:NCHAN) 
          JAC_CLD_2(1:NLAY,1:NCHAN)  = CFRA2X*JAC_CLD_2(1:NLAY,1:NCHAN)  * L2S4above(3,1:NLAY,1:NCHAN) 
          JAC_CLD_12(1:NLAY,1:NCHAN) = CFRA12*JAC_CLD_12(1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN)                      

          JAC_CLD_C(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) = JAC_CLD_C(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_1(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_2(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_12(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN)

          JAC_CLD_OUT(7,1:NCHAN) = 0
          DO IJUNK = IEFFCLD_TOP2,IEFFCLD_BOT2
            JAC_CLD_OUT(7,1:NCHAN) = JAC_CLD_OUT(7,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN)
          END DO
        END IF 

        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write it out!!!
        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write it out!!!
        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write it out!!!
        CALL WRTJAC_CLD(IOUNCLD,IPROF,7,NCHAN,300,FREQ,RAD,JAC_OUTPUT_UNITS,JAC_CLD_OUT)
        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write it out!!!
        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write it out!!!
        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write it out!!!

      END IF   
c%%%%%%%%%%%%%%%%%%%%%%%%%
