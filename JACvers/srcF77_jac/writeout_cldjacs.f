      IF (INTERSECT(300,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN

c        IF (NCHAN .EQ. 1) THEN
c          print *,'JACA_FINAL_1(1) = ',JACA_FINAL_1(1)        
c          print *,'JACA_FINAL_2(1) = ',JACA_FINAL_2(1)        
c          print *,'JACS_FINAL_1(1) = ',JACS_FINAL_1(1)        
c          print *,'JACS_FINAL_2(1) = ',JACS_FINAL_2(1)        
c        END IF

        !!! would be [cfrac1,amnt1,sze1,ctop1,cbot1],[cfrac2,amt2,sze2,ctop2,cbot2],cfrac12

        JAC_CLD_OUT = 0  
        ICLDJAC = 0

        CLDJACFAKEAMT( 1) = CFRAC1
        CLDJACFAKEAMT( 2) = CNGWA1
        CLDJACFAKEAMT( 3) = CPSIZ1
        CLDJACFAKEAMT( 4) = CPRTO1
        CLDJACFAKEAMT( 5) = CPRBO1

        CLDJACFAKEAMT( 6) = CFRAC2
        CLDJACFAKEAMT( 7) = CNGWA2
        CLDJACFAKEAMT( 8) = CPSIZ2
        CLDJACFAKEAMT( 9) = CPRTO2
        CLDJACFAKEAMT(10) = CPRBO2

        CLDJACFAKEAMT(11) = CFRA12

        CLDJACFAKEAMT(12) = 1      !! this is stemp so don't worry

c        print *,'CLDJACFAKEAMT = ',CLDJACFAKEAMT
c        IEFFCLD_TOP1 = LCTOP1
c        IEFFCLD_TOP2 = LCTOP2
c        IEFFCLD_BOT1 = LCBOT1
c        IEFFCLD_BOT2 = LCBOT2

c %%%%%%%%%%%%%%%%%%%%%%%%%        
c cld1 : cfrac1, cngwat1, cpsize1,cprtop1,cprbot1
        !! ** <> ** <> ** <> ** <> ** !!
        IF ((CFRAC1 .GT. 0) .AND. (CNGWA1 .GT. 0) .AND.
     $      (LCBOT1 .GT. 0)    .AND. (LCTOP1 .GT. 0) .AND. 
     $      (LCBOT1 .LE. LBOT) .AND. (LCTOP1 .LE. LBOT)) THEN

          ICLDJAC = ICLDJAC + 1 !!!! 1
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = -PURE_RAD4(1,1,1:NCHAN) + PURE_RAD4(2,1,1:NCHAN)  !!! dr/d cfrac1 = -rclr + r1

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld amt jac, cld 1
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
              
          DO IJUNK = LCTOP1,LCBOT1
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

          JAC_CLD_C(LCTOP1:LCBOT1,1:NCHAN) = JAC_CLD_C(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_1(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_2(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_12(LCTOP1:LCBOT1,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 2
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = LCTOP1,LCBOT1
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL1(IJUNK)
          END DO

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld sze jac, cld 1
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = LCTOP1,LCBOT1
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

          JAC_CLD_C(LCTOP1:LCBOT1,1:NCHAN) = JAC_CLD_C(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_1(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_2(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_12(LCTOP1:LCBOT1,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 3
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = LCTOP1,LCBOT1
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL1(IJUNK)
          END DO

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld top jac, cld 1
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = LCTOP1,LCBOT1
             RAMU(IJUNK,1:NCHAN) = JACTOP_CFRCL1_v(IJUNK,1:NCHAN)
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

          JAC_CLD_C(LCTOP1:LCBOT1,1:NCHAN) = JAC_CLD_C(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_1(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_2(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_12(LCTOP1:LCBOT1,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 4
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = LCTOP1,LCBOT1
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL1(IJUNK)
          END DO

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld bot jac, cld 1
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = LCTOP1,LCBOT1
             RAMU(IJUNK,1:NCHAN) = JACBOT_CFRCL1_v(IJUNK,1:NCHAN)
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

          JAC_CLD_C(LCTOP1:LCBOT1,1:NCHAN) = JAC_CLD_C(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_1(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_2(LCTOP1:LCBOT1,1:NCHAN) + 
     $                                                   JAC_CLD_12(LCTOP1:LCBOT1,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 5
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = LCTOP1,LCBOT1
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL1(IJUNK)
          END DO
        ELSE
          ICLDJAC = ICLDJAC + 5    !!! skipped 5 cld1 jacs
        END IF 
  
c %%%%%%%%%%%%%%%%%%%%%%%%%        
c cld1 : cfrac1, cngwat1, cpsize1,cprtop1,cprbot1
        IF ((CFRAC2 .GT. 0) .AND. (CNGWA2 .GT. 0) .AND.
     $      (LCBOT2 .GT. 0)    .AND. (LCTOP2 .GT. 0) .AND. 
     $      (LCBOT2 .LE. LBOT) .AND. (LCTOP2 .LE. LBOT)) THEN

          ICLDJAC = ICLDJAC + 1 !!!! 6
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = -PURE_RAD4(1,1,1:NCHAN) + PURE_RAD4(3,1,1:NCHAN)  !!! dr/d cfrac1 = -rclr + r2

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld amt jac, cld 2
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = LCTOP2,LCBOT2
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

          JAC_CLD_C(LCTOP2:LCBOT2,1:NCHAN) = JAC_CLD_C(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_1(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_2(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_12(LCTOP2:LCBOT2,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 7
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = LCTOP2,LCBOT2
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL2(IJUNK)
          END DO

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld sze jac, cld 2
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = LCTOP2,LCBOT2
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

          JAC_CLD_C(LCTOP2:LCBOT2,1:NCHAN) = JAC_CLD_C(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_1(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_2(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_12(LCTOP2:LCBOT2,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 8
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = LCTOP2,LCBOT2
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL2(IJUNK)
          END DO

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld top jac, cld 2
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = LCTOP2,LCBOT2
             RAMU(IJUNK,1:NCHAN) = JACTOP_CFRCL2_v(IJUNK,1:NCHAN)
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
          JAC_CLD_C(LCTOP2:LCBOT2,1:NCHAN) = JAC_CLD_C(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_1(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_2(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_12(LCTOP2:LCBOT2,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 9
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = LCTOP2,LCBOT2
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL2(IJUNK)
          END DO

          !! ** <> ** <> ** <> ** <> ** !!
          !! cld bot jac, cld 2
          RAMU = 0
          JAC_CLD_C  = 0  !! no clouds
          JAC_CLD_1  = 0
          JAC_CLD_2  = 0
          JAC_CLD_12 = 0
            
          DO IJUNK = LCTOP2,LCBOT2
             RAMU(IJUNK,1:NCHAN) = JACBOT_CFRCL2_v(IJUNK,1:NCHAN)
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

          JAC_CLD_C(LCTOP2:LCBOT2,1:NCHAN) = JAC_CLD_C(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_1(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_2(LCTOP2:LCBOT2,1:NCHAN) + 
     $                                                   JAC_CLD_12(LCTOP2:LCBOT2,1:NCHAN)

          ICLDJAC = ICLDJAC + 1 !!!! 10
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = 0
          DO IJUNK = LCTOP2,LCBOT2
            JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_CLD_OUT(ICLDJAC,1:NCHAN) + JAC_CLD_C(IJUNK,1:NCHAN) * CFRCL2(IJUNK)
          END DO
        ELSE
          ICLDJAC = ICLDJAC + 5    !!! skipped 5 cld2 jacs        
        END IF 

c %%%%%%%%%%%%%%%%%%%%%%%%%
c cld1+cld2 : cfrac12 
        IF ((CFRAC1 .GT. 0) .AND. (CNGWA1 .GT. 0) .AND. (CFRAC2 .GT. 0) .AND. (CNGWA2 .GT. 0)) THEN
          ICLDJAC = ICLDJAC + 1 !!!! 11
          JAC_CLD_OUT(ICLDJAC,1:NCHAN) = +PURE_RAD4(1,1,1:NCHAN) + PURE_RAD4(4,1,1:NCHAN)  !!! dr/d cfrac12 = rclr + r12 - r1 - r2
     $                                  - PURE_RAD4(2,1,1:NCHAN) - PURE_RAD4(3,1,1:NCHAN)  
        ELSE
          ICLDJAC = ICLDJAC + 1    !!! skipped 1 cld1,2 jacs        
        END IF

c %%%%%%%%%%%%%%%%%%%%%%%%%
c stemp jac for fun
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
        ICLDJAC = ICLDJAC + 1 !!!! 12
        JAC_CLD_OUT(ICLDJAC,1:NCHAN) = JAC_ST_C(1:NCHAN)

c %%%%%%%%%%%%%%%%%%%%%%%%%
c write it all out
        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write out the [5 cld1, 5 cld2, 1 cld12, 1 stemp jacs ] !!!
!        CALL WRTJAC_CLD(IOUNCLD,IPROF,11,NCHAN,300,FREQ,RAD,JAC_OUTPUT_UNITS*0,CLDJACFAKEAMT,JAC_CLD_OUT)    !!! test raw rad diffs
        CALL WRTJAC_CLD(IOUNCLD,IPROF,CLDJAC,NCHAN,300,FREQ,RAD,JAC_OUTPUT_UNITS*1,CLDJACFAKEAMT,JAC_CLD_OUT)
        !! ** <> ** <> ** <> ** <> ** !!
        !! finally, write it out!!!

      END IF   
c %%%%%%%%%%%%%%%%%%%%%%%%%
