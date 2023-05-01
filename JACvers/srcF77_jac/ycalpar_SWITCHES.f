        IF (LCO2) THEN
           IF (LCO2PM) THEN
              CO2MLT(L)=100.0*(PFAMNT(L) - CO2STD)/(3.0*CO2STD)
           ELSE
C             CO2 mult=1 when prof amount = 1.03 * ref amount
              CO2MLT(L)=33.3333*( PFAMNT(L) - FIXMUL(L)*RFAMNT(L) )/RFAMNT(L)
C             Ignore changes in CO2 of less than ~0.03%
              IF (ABS(CO2MLT(L)) .LT. 1E-2) CO2MLT(L)=0.0
           ENDIF
        ELSE
           CO2MLT(L)=100.0*(CO2PPM - CO2STD)/(3.0*CO2STD)
        ENDIF
        IF (L .LE. NTEBOT) THEN
           CO2TOP=CO2TOP + CO2STD*(1.0 + CO2MLT(L)*3.0E-2)
        ENDIF
C       if (DEBUG) write(6,'(a,X,I4,3(X,ES11.3E3))') 'calpar: L,CO2MLT(L) ',L,PFAMNT(L),RFAMNT(L),CO2MLT(L)
C
        IF (LN2O) THEN
C          N2O mult=-1 when prof amount = 0.75 * ref amount
           N2OMLT(L)=4.0*( PNAMNT(L) - FIXMUL(L)*RNAMNT(L) )/RNAMNT(L)
C          Ignore changes in N2O less than ~0.3%
           IF (ABS(N2OMLT(L)) .LT. 1E-2) N2OMLT(L)=0.0
        ELSE
           N2OMLT(L)=0.0
        ENDIF
C
        IF (LSO2) THEN
C          SO2 mult=1 when prof amount = 1000 * ref amount
           SO2MLT(L)=1.0010E-3*( PSAMNT(L) - FIXMUL(L)*RSAMNT(L) )/ RSAMNT(L)
C          Ignore changes in SO2 of less than ~10%
           IF (ABS(SO2MLT(L)) .LT. 1E-4) SO2MLT(L)=0.0
        ELSE
           SO2MLT(L)=0.0
        ENDIF
C
        IF (LNH3) THEN
C          NH3 mult=1 when prof amount = 100 * ref amount
           NH3MLT(L)=1.0101E-2*( PAAMNT(L) - FIXMUL(L)*RAAMNT(L) )/RAAMNT(L)
C          Ignore changes in NH3 of less than ~10%
           IF (ABS(NH3MLT(L)) .LT. 1E-4) NH3MLT(L)=0.0
        ELSE
           NH3MLT(L)=0.0
        ENDIF
C
        IF (LHNO3) THEN
C          HNO3 mult=1 when prof amount = 2 * ref amount
           HNOMLT(L)=( PHAMNT(L) - FIXMUL(L)*RHAMNT(L) )/ RHAMNT(L)
C          Ignore changes in HNO3 less than ~1%
           IF (ABS(HNOMLT(L)) .LT. 1E-2) HNOMLT(L)=0.0
        ELSE
           HNOMLT(L)=0.0
        ENDIF
C
        IF (LHDO) THEN
C          HDO mult=1 when no depletion/enhancement of HDO
           HDOMLT(L)=( 1 - HDOFCT )
C          Ignore changes in HDO of less than ~1%
C           IF (ABS(HDOMLT(L)) .LT. 1E-5) HDOMLT(L)=0.0
        ELSE
           HDOMLT(L)=0.0
        ENDIF

