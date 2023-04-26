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
      DT_T = 1.0
      TR_T = 1/RTEMP(L)
CjacXC     water terms
      DT_1 = 0.0
      TR_1 = 0.0
CjacXC     ozone terms
      DT_3 = 0.0
      TR_3 = 0.0
CjacXC
CjacXC     Calc the fixed gases correction term for this layer
      PWATER = KMOLE*PWAMNT(L)*PTEMP(L)/(STDDEN*STDTMP*100*DZREF(L))
      PWATER_T = PWATER/PTEMP(L)
      PWATER_1 = PWATER/PWAMNT(L)
      PWATER_3 = 0.0
      A_F = ( 1 - PMULT*PWATER/PRES(L) )/( FX(L)*GSCAL )
      A_F_T = 1/(FX(L)*GSCAL)/PRES(L)*(-PMULT*PWATER_T)
      A_F_1 = 1/(FX(L)*GSCAL)/PRES(L)*(-PMULT*PWATER_1)
      A_F_3 = 0.0
CCjacXccc
CCjacXc for testing
CCjacXc      A_F = 1.0
CCjacXccc
CjacXC
CjacXC     Water terms
      A_W = PWAMNT(L)/RWAMNT(L)
      A_W_T = 0.0
      A_W_1 = 1/RWAMNT(L)
      A_W_3 = 0.0
CjacX      WZREF = WZREF + PDP*RWAMNT(L)
CjacX      WZ = WZ + PDP*PWAMNT(L)
      AZ_W = WZ/WZREF
      AZ_W_T = 0.0
      AZ_W_1 = 0.0
      AZ_W_1 = AZ_W_1 + PDP
      AZ_W_3 = 0.0
CjacXC
CjacXC     Ozone terms
      A_O = POAMNT(L)/ROAMNT(L)
      A_O_T = 0.0
      A_O_1 = 0.0
      A_O_3 = 1/ROAMNT(L)
CjacX      XZREF = XZREF + ROAMNT(L)
CjacX      XZ = XZ + POAMNT(L)
      XZ_O = XZ/XZREF
      XZ_O_T = 0.0
      XZ_O_1 = 0.0
      XZ_O_3 = 1/XZREF
CjacX      OZREF = OZREF + PDP*ROAMNT(L)
CjacX      OZ = OZ + PDP*POAMNT(L)
      AZ_O = OZ/OZREF
      AZ_O_T = 0.0
      AZ_O_1 = 0.0
      AZ_O_3 = 1/OZREF

CjacXC
CjacXC     Carbon monoxide terms
      A_C = PCAMNT(L)/RCAMNT(L)
CjacX      CZREF = CZREF + PDP*RCAMNT(L)
CjacX      CZ = CZ + PDP*PCAMNT(L)
      AZ_C = CZ/CZREF
CjacXC
CjacXC     Methane terms
      A_M = PMAMNT(L)/RMAMNT(L)
CjacX      MZREF = MZREF + PDP*RMAMNT(L)
CjacX      MZ = MZ + PDP*PMAMNT(L)
      AZ_M = MZ/MZREF

