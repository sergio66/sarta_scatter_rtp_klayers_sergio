      IF (INTERSECT(300,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN

        IF (NCHAN .EQ. 1) THEN
          print *,'JACA_FINAL_1(1) = ',JACA_FINAL_1(1)        
          print *,'JACA_FINAL_2(1) = ',JACA_FINAL_2(1)        
          print *,'JACS_FINAL_1(1) = ',JACS_FINAL_1(1)        
          print *,'JACS_FINAL_2(1) = ',JACS_FINAL_2(1)        
        END IF

        !!! would be cfrac1,cfrac2,cfrac12,amt1,amt2,sze1,sze2

        JAC_CLD_OUT = 0  
        ICLDJAC = 0
        CLDJACFAKEAMT(1) = CFRAC1
        CLDJACFAKEAMT(2) = CNGWA1
        CLDJACFAKEAMT(3) = CPSIZ1
        CLDJACFAKEAMT(4) = CFRAC2
        CLDJACFAKEAMT(5) = CNGWA2
        CLDJACFAKEAMT(6) = CPSIZ2
        CLDJACFAKEAMT(7) = CFRA12
        print *,'CLDJACFAKEAMT = ',CLDJACFAKEAMT
        
        !! ** <> ** <> ** <> ** <> ** !!
        !! figure out cloud top and cloud bottom layers
        IEFFCLD_TOP1 = -1
        IEFFCLD_BOT1 = -1
        IF ((CFRAC1 .GT. 0) .AND. (CNGWA1 .GT. 0)) THEN
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
            print *,CFRCL1(IEFFCLD_TOP1-1:IEFFCLD_BOT1+1)
          ELSE
            IEFFCLD_TOP1 = -1
            IEFFCLD_BOT1 = -1
          END IF
        END IF

        IEFFCLD_TOP2 = -1
        IEFFCLD_BOT2 = -1
        IF ((CFRAC2 .GT. 0) .AND. (CNGWA2 .GT. 0)) THEN
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
            print *,CFRCL2(IEFFCLD_TOP2-1:IEFFCLD_BOT2+1)
          ELSE
            IEFFCLD_TOP2 = -1
            IEFFCLD_BOT2 = -1
          END IF 
        END IF

c%%%%%%%%%%%%%%%%%%%%%%%%%        

        !! ** <> ** <> ** <> ** <> ** !!
        !! cld frac jacs cfrac1,cfrac2, cfrac12
        IF ((CFRAC1 .GT. 0) .AND. (CNGWA1 .GT. 0) .AND.
     $      (IEFFCLD_BOT1 .GT. 0)    .AND. (IEFFCLD_TOP1 .GT. 0) .AND. 
     $      (IEFFCLD_BOT1 .LE. LBOT) .AND. (IEFFCLD_TOP1 .LE. LBOT)) THEN

          ICLDJAC = ICLDJAC + 1 !!!! 1
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = -PURE_RAD4(1,1,1:NCHAN) + PURE_RAD4(2,1,1:NCHAN)  !!! dr/d cfrac1 = -rclr + r1

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld amt jac, cld 1
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
              
          DO IJUNK = IEFFCLD_TOP1,IEFFCLD_BOT1
             RAMU(IJUNK,1:NCHAN) = JACA_FINAL_1(1:NCHAN)
          END DO          

          JAC_CLD_C(1:NLAY,1:NCHAN)  = PLANCK_RAD4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = PLANCK_RAD4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = PLANCK_RAD4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = PLANCK_RAD4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = JAC_CLD_C(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = JAC_CLD_1(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = JAC_CLD_2(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = JAC_CLD_12(1:NLAY,1:NCHAN) - RTHERM4_SOLAR4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = FCLEAR*JAC_CLD_C(1:NLAY,1:NCHAN)  * L2S4above(1,1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = CFRA1X*JAC_CLD_1(1:NLAY,1:NCHAN)  * L2S4above(2,1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = CFRA2X*JAC_CLD_2(1:NLAY,1:NCHAN)  * L2S4above(3,1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = CFRA12*JAC_CLD_12(1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN) * 1                     

          JAC_CLD_C(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) = JAC_CLD_C(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_1(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_2(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_12(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 2
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = IEFFCLD_TOP1,IEFFCLD_BOT1
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL1(IJUNK)
          END DO

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld sze jac, cld 1
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = IEFFCLD_TOP1,IEFFCLD_BOT1
             RAMU(IJUNK,1:NCHAN) = JACS_FINAL_1(1:NCHAN)
          END DO          
          JAC_CLD_C(1:NLAY,1:NCHAN)  = PLANCK_RAD4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = PLANCK_RAD4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = PLANCK_RAD4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = PLANCK_RAD4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = JAC_CLD_C(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = JAC_CLD_1(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = JAC_CLD_2(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = JAC_CLD_12(1:NLAY,1:NCHAN) - RTHERM4_SOLAR4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = FCLEAR*JAC_CLD_C(1:NLAY,1:NCHAN)  * L2S4above(1,1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = CFRA1X*JAC_CLD_1(1:NLAY,1:NCHAN)  * L2S4above(2,1:NLAY,1:NCHAN) * 1
          JAC_CLD_2(1:NLAY,1:NCHAN)  = CFRA2X*JAC_CLD_2(1:NLAY,1:NCHAN)  * L2S4above(3,1:NLAY,1:NCHAN) * 0
          JAC_CLD_12(1:NLAY,1:NCHAN) = CFRA12*JAC_CLD_12(1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN) * 1                     

          JAC_CLD_C(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) = JAC_CLD_C(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_1(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_2(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN) + 
     $                                                   JAC_CLD_12(IEFFCLD_TOP1:IEFFCLD_BOT1,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 3
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = IEFFCLD_TOP1,IEFFCLD_BOT1
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL1(IJUNK)
          END DO
        END IF 
  
c%%%%%%%%%%%%%%%%%%%%%%%%%        

        IF ((CFRAC2 .GT. 0) .AND. (CNGWA2 .GT. 0) .AND.
     $      (IEFFCLD_BOT2 .GT. 0)    .AND. (IEFFCLD_TOP2 .GT. 0) .AND. 
     $      (IEFFCLD_BOT2 .LE. LBOT) .AND. (IEFFCLD_TOP2 .LE. LBOT)) THEN

          ICLDJAC = ICLDJAC + 1 !!!! 4
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = -PURE_RAD4(1,1,1:NCHAN) + PURE_RAD4(3,1,1:NCHAN)  !!! dr/d cfrac1 = -rclr + r2

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld amt jac, cld 2
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = IEFFCLD_TOP2,IEFFCLD_BOT2
             RAMU(IJUNK,1:NCHAN) = JACA_FINAL_2(1:NCHAN)
          END DO          
          JAC_CLD_C(1:NLAY,1:NCHAN)  = PLANCK_RAD4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = PLANCK_RAD4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = PLANCK_RAD4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = PLANCK_RAD4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = JAC_CLD_C(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = JAC_CLD_1(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = JAC_CLD_2(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = JAC_CLD_12(1:NLAY,1:NCHAN) - RTHERM4_SOLAR4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = FCLEAR*JAC_CLD_C(1:NLAY,1:NCHAN)  * L2S4above(1,1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = CFRA1X*JAC_CLD_1(1:NLAY,1:NCHAN)  * L2S4above(2,1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = CFRA2X*JAC_CLD_2(1:NLAY,1:NCHAN)  * L2S4above(3,1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = CFRA12*JAC_CLD_12(1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN) * 1                     

          JAC_CLD_C(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) = JAC_CLD_C(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_1(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_2(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_12(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 5
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = IEFFCLD_TOP2,IEFFCLD_BOT2
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL2(IJUNK)
          END DO

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld sze jac, cld 2
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = IEFFCLD_TOP2,IEFFCLD_BOT2
             RAMU(IJUNK,1:NCHAN) = JACS_FINAL_2(1:NCHAN)
          END DO          
          JAC_CLD_C(1:NLAY,1:NCHAN)  = PLANCK_RAD4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = PLANCK_RAD4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = PLANCK_RAD4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = PLANCK_RAD4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = JAC_CLD_C(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = JAC_CLD_1(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(2,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = JAC_CLD_2(1:NLAY,1:NCHAN)  - RTHERM4_SOLAR4(3,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = JAC_CLD_12(1:NLAY,1:NCHAN) - RTHERM4_SOLAR4(4,1:NLAY,1:NCHAN) * RAMU(1:NLAY,1:NCHAN) * 1
  
          JAC_CLD_C(1:NLAY,1:NCHAN)  = FCLEAR*JAC_CLD_C(1:NLAY,1:NCHAN)  * L2S4above(1,1:NLAY,1:NCHAN) * 0
          JAC_CLD_1(1:NLAY,1:NCHAN)  = CFRA1X*JAC_CLD_1(1:NLAY,1:NCHAN)  * L2S4above(2,1:NLAY,1:NCHAN) * 0
          JAC_CLD_2(1:NLAY,1:NCHAN)  = CFRA2X*JAC_CLD_2(1:NLAY,1:NCHAN)  * L2S4above(3,1:NLAY,1:NCHAN) * 1
          JAC_CLD_12(1:NLAY,1:NCHAN) = CFRA12*JAC_CLD_12(1:NLAY,1:NCHAN) * L2S4above(4,1:NLAY,1:NCHAN) * 1                     

          JAC_CLD_C(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) = JAC_CLD_C(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_1(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_2(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN) + 
     $                                                   JAC_CLD_12(IEFFCLD_TOP2:IEFFCLD_BOT2,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 6
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = IEFFCLD_TOP2,IEFFCLD_BOT2
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL2(IJUNK)
          END DO
        END IF 

c%%%%%%%%%%%%%%%%%%%%%%%%%
        IF ((CFRAC1 .GT. 0) .AND. (CNGWA1 .GT. 0) .AND. (CFRAC2 .GT. 0) .AND. (CNGWA2 .GT. 0)) THEN
          ICLDJAC = ICLDJAC + 1 !!!! 7
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = +PURE_RAD4(1,1,1:NCHAN) + PURE_RAD4(4,1,1:NCHAN)  !!! dr/d cfrac12 = rclr + r12 - r1 - r2
     $                                  - PURE_RAD4(2,1,1:NCHAN) - PURE_RAD4(3,1,1:NCHAN)  
        END IF

c%%%%%%%%%%%%%%%%%%%%%%%%%

        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write it out!!!
!        CALL WRTJAC_CLD(IOUNCLD,IPROF,7,NCHAN,300,FREQ,RAD,JAC_OUTPUT_UNITS*0,CLDJACFAKEAMT,JAC_CLD_OUT)    !!! test raw rad diffs
        CALL WRTJAC_CLD(IOUNCLD,IPROF,CLDJAC,NCHAN,300,FREQ,RAD,JAC_OUTPUT_UNITS*1,CLDJACFAKEAMT,JAC_CLD_OUT)
        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write it out!!!

      END IF   
c%%%%%%%%%%%%%%%%%%%%%%%%%
