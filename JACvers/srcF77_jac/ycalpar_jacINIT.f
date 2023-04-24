CjacXC
CjacXC     ---------------------------
CjacXC     Calculate the basic profile
CjacXC     dependent predictors DERIVATIVES
CjacXC     ---------------------------
CjacXC
CjacX      IF (L .EQ. 1) THEN
CjacX         PDP = PRES(1)*( PRES(2) - PRES(1))
CjacX         TRZ = 0.0E+0
CjacX         TAZ_O = 0.0E+0
CjacX         TAZ_M = 0.0E+0
CjacX      ELSE
CjacX         PDP = PRES(L)*( PRES(L) - PRES(L-1) )
CjacX         PNORM = PNORM + PDP
CjacXC
CjacXC        Note: TRZ, TOZ, and TMZ use layer-above terms
CjacX         TZ = TZ + PDP*TR
CjacX         TRZ = TZ/PNORM
CjacXC
CjacX         TOZ = TOZ + PDP*DT*A_O
CjacX         TAZ_O = TOZ/PNORM
CjacXC
CjacX         TMZ = TMZ + PDP*TR*A_M
CjacX         TAZ_M = TMZ/PNORM
CjacX      ENDIF
CjacXC
CjacXC     Temperature terms
      DT = PTEMP(L) - RTEMP(L)
      TR = PTEMP(L)/RTEMP(L)
CjacXC
CjacXC     Calc the fixed gases correction term for this layer
      PWATER = KMOLE*PWAMNT(L)*PTEMP(L)/(STDDEN*STDTMP*100*DZREF(L))
      A_F = ( 1 - PMULT*PWATER/PRES(L) )/( FX(L)*GSCAL )
CCjacXccc
CCjacXc for testing
CCjacXc      A_F = 1.0
CCjacXccc
CjacXC
CjacXC     Water terms
      A_W = PWAMNT(L)/RWAMNT(L)
CjacX      WZREF = WZREF + PDP*RWAMNT(L)
      WZ = WZ + PDP*PWAMNT(L)
      AZ_W = WZ/WZREF
CjacXC
CjacXC     Ozone terms
      A_O = POAMNT(L)/ROAMNT(L)
CjacX      XZREF = XZREF + ROAMNT(L)
      XZ = XZ + POAMNT(L)
      XZ_O = XZ/XZREF
CjacX      OZREF = OZREF + PDP*ROAMNT(L)
      OZ = OZ + PDP*POAMNT(L)
      AZ_O = OZ/OZREF
CjacXC
CjacXC     Carbon monoxide terms
      A_C = PCAMNT(L)/RCAMNT(L)
CjacX      CZREF = CZREF + PDP*RCAMNT(L)
      CZ = CZ + PDP*PCAMNT(L)
      AZ_C = CZ/CZREF
CjacXC
CjacXC     Methane terms
      A_M = PMAMNT(L)/RMAMNT(L)
CjacX      MZREF = MZREF + PDP*RMAMNT(L)
      MZ = MZ + PDP*PMAMNT(L)
      AZ_M = MZ/MZREF

