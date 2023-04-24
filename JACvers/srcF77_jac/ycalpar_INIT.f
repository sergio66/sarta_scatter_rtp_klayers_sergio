C
C     ---------------------------
C     Calculate the basic profile
C     dependent predictors.
C     ---------------------------
C
      IF (L .EQ. 1) THEN
         PDP = PRES(1)*( PRES(2) - PRES(1))
         TRZ = 0.0E+0
         TAZ_O = 0.0E+0
         TAZ_M = 0.0E+0
      ELSE
         PDP = PRES(L)*( PRES(L) - PRES(L-1) )
         PNORM = PNORM + PDP
C
C        Note: TRZ, TOZ, and TMZ use layer-above terms
         TZ = TZ + PDP*TR
         TRZ = TZ/PNORM
C
         TOZ = TOZ + PDP*DT*A_O
         TAZ_O = TOZ/PNORM
C
         TMZ = TMZ + PDP*TR*A_M
         TAZ_M = TMZ/PNORM
      ENDIF
C
C     Temperature terms
      DT = PTEMP(L) - RTEMP(L)
      TR = PTEMP(L)/RTEMP(L)
C
C     Calc the fixed gases correction term for this layer
      PWATER = KMOLE*PWAMNT(L)*PTEMP(L)/(STDDEN*STDTMP*100*DZREF(L))
      A_F = ( 1 - PMULT*PWATER/PRES(L) )/( FX(L)*GSCAL )
ccc
c for testing
c      A_F = 1.0
ccc
C
C     Water terms
      WZREF = WZREF + PDP*RWAMNT(L)
      WZ    = WZ    + PDP*PWAMNT(L)
      A_W  = PWAMNT(L)/RWAMNT(L)
      AZ_W = WZ/WZREF
C
C     Ozone terms
      A_O = POAMNT(L)/ROAMNT(L)
      XZREF = XZREF + ROAMNT(L)
      XZ    = XZ   + POAMNT(L)
      XZ_O  = XZ/XZREF
      OZREF = OZREF + PDP*ROAMNT(L)
      OZ    = OZ    + PDP*POAMNT(L)
      AZ_O  = OZ/OZREF
C
C     Carbon monoxide terms
      A_C   = PCAMNT(L)/RCAMNT(L)
      CZREF = CZREF + PDP*RCAMNT(L)
      CZ    = CZ + PDP*PCAMNT(L)
      AZ_C  = CZ/CZREF
C
C     Methane terms
      A_M   = PMAMNT(L)/RMAMNT(L)
      MZREF = MZREF + PDP*RMAMNT(L)
      MZ    = MZ + PDP*PMAMNT(L)
      AZ_M  = MZ/MZREF

