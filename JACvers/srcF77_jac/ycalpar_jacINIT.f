c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%% this is from ycalpar_SWITCHES.f %%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IF (DOJAC) THEN
          IF (LCO2) THEN
            CO2JACMLT(L)=33.3333*(1)/RFAMNT(L)
          ELSE !! if LCO2ppm
            CO2JACMLT(L) = 100.0/(3.0*CO2STD)*1
          END IF

          IF (LN2O) THEN
            N2OJACMLT(L) = 4.0/RNAMNT(L)
          ELSE
            N2OJACMLT(L) = 0.0
          END IF

          IF (LSO2) THEN
            SO2JACMLT(L) = 1.0010E-3/RSAMNT(L)
          ELSE
            SO2JACMLT(L) = 0.0
          END IF

          IF (LNH3) THEN
            NH3JACMLT(L) = 1.0101E-2/RAAMNT(L)
          ELSE
            NH3JACMLT(L) = 0.0
          END IF

          IF (LHNO3) THEN
            HNOJACMLT(L) = 1/RHAMNT(L)
          ELSE
            HNOJACMLT(L) = 0.0
          END IF

          IF (LHDO) THEN
            HDOJACMLT(L) = -1
          ELSE
            HDOJACMLT(L) = 0.0
          END IF
        END IF

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%% this is from ycalpar_SWITCHES.f %%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CjacXC
CjacXC     ---------------------------
CjacXC     Calculate the basic profile
CjacXC     dependent predictors DERIVATIVES
CjacXC     ---------------------------
CjacXC
           IF (L .EQ. 1) THEN
              PDP = PRES(1)*( PRES(2) - PRES(1))
CjacX         TRZ = 0.0E+0
              TRZ_T = 0
              TRZ_1 = 0
              TRZ_3 = 0
CjacX         TAZ_O = 0.0E+0
              TAZ_O_T = 0
              TAZ_O_1 = 0
              TAZ_O_3 = 0
CjacX         TAZ_M = 0.0E+0
              TAZ_M_T = 0
              TAZ_M_1 = 0
              TAZ_M_3 = 0

              AZ_W_T = 0
              AZ_W_1 = 0
              AZ_W_3 = 0

              DT_T = 1
              DT_1 = 0
              DT_3 = 0
              DT_5 = 0

              TZ_T = PDP/RTEMP(L)
              TZ_1 = 0.0
              TZ_3 = 0.0

              TRZ_T = 0.0
              TRZ_1 = 0.0
              TRZ_3 = 0.0

              TOZ_T = PDP*(DT_T*A_O + DT*A_O_T)
              TOZ_1 = PDP*(DT_1*A_O + DT*A_O_1)
              TOZ_3 = PDP*(DT_3*A_O + DT*A_O_3)

CjacX         TAZ_O = TOZ/PNORM
!              TAZ_O_T = TOZ_T/PNORM
!              TAZ_O_1 = TOZ_1/PNORM
!              TAZ_O_3 = TOZ_3/PNORM
              TAZ_O_T = 0
              TAZ_O_1 = 0
              TAZ_O_3 = 0
CjacXC
CjacX         TMZ = TMZ + PDP*TR*A_M
CjacX         TAZ_M = TMZ/PNORM
!              TAZ_M_T = 1/PNORM*PDP*(TR_T*A_M + TR*A_M_T)
!              TAZ_M_1 = 1/PNORM*PDP*(TR_1*A_M + TR*A_M_1)
!              TAZ_M_3 = 1/PNORM*PDP*(TR_3*A_M + TR*A_M_3)
!              TAZ_M_6 = 1/PNORM*PDP*(TR_6*A_M + TR*A_M_6)
              TAZ_M_T = 0
              TAZ_M_1 = 0
              TAZ_M_3 = 0
              TAZ_M_6 = 0

            ELSE
              PDP = PRES(L)*( PRES(L) - PRES(L-1) )
CjacX         PNORM = PNORM + PDP
CjacXC
CjacXC        Note: TRZ, TOZ, and TMZ use layer-above terms
CjacX         TZ = TZ + PDP*TR
              TZ_T = PDP/RTEMP(L)
              !! TZ_T = PDP*PTEMP(L)/TZ  !! tobias wehr HIRS writeup
              TZ_1 = 0.0
              TZ_3 = 0.0
CjacX         TRZ = TZ/PNORM
              TRZ_T = TZ_T/PNORM
              TRZ_1 = 0
              TRZ_3 = 0
CjacXC
CjacX         TOZ = TOZ + PDP*DT*A_O
              TOZ_T = PDP*(DT_T*A_O + DT*A_O_T)
              TOZ_1 = PDP*(DT_1*A_O + DT*A_O_1)
              TOZ_3 = PDP*(DT_3*A_O + DT*A_O_3)

CjacX         TAZ_O = TOZ/PNORM
              TAZ_O_T = TOZ_T/PNORM
              TAZ_O_1 = TOZ_1/PNORM
              TAZ_O_3 = TOZ_3/PNORM
CjacXC
CjacX         TMZ = TMZ + PDP*TR*A_M
CjacX         TAZ_M = TMZ/PNORM
              TAZ_M_T = 1/PNORM*PDP*(TR_T*A_M + TR*A_M_T)
              TAZ_M_1 = 1/PNORM*PDP*(TR_1*A_M + TR*A_M_1)
              TAZ_M_3 = 1/PNORM*PDP*(TR_3*A_M + TR*A_M_3)
              TAZ_M_6 = 1/PNORM*PDP*(TR_6*A_M + TR*A_M_6)
            ENDIF
CjacXC
CjacXC     Temperature terms
      DT_T = 1.0
      DT_1 = 0.0
      DT_3 = 0.0
      DT_6 = 0.0

      TR_T = 1/RTEMP(L)
      TR_1 = 0.0
      TR_3 = 0.0
      TR_6 = 0.0
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
      AZ_W_1 = AZ_W_1 + PDP/WZREF
      AZ_W_1 = PDP/WZREF   !! tobias Wehr HIRS write up
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
      AZ_O_3 = PDP/OZREF

CjacXC
CjacXC     Carbon monoxide terms
      A_C = PCAMNT(L)/RCAMNT(L)
      A_C_T = 0
      A_C_1 = 0
      A_C_3 = 0
      A_C_5 = 1/RCAMNT(L)
CjacX      CZREF = CZREF + PDP*RCAMNT(L)
CjacX      CZ = CZ + PDP*PCAMNT(L)
      AZ_C = CZ/CZREF
      AZ_C_T = 0.0
      AZ_C_1 = 0.0
      AZ_C_3 = 0.0
      AZ_C_5 = PDP
CjacXC
CjacXC     Methane terms
      A_M = PMAMNT(L)/RMAMNT(L)
      A_M_T = 0
      A_M_1 = 0
      A_M_3 = 0
      A_M_6 = 1/RMAMNT(L)
CjacX      MZREF = MZREF + PDP*RMAMNT(L)
CjacX      MZ = MZ + PDP*PMAMNT(L)
      MZ_T   = 0
      MZ_1   = 0
      MZ_3   = 0
      MZ_6   = PDP
      AZ_M = MZ/MZREF
      AZ_M_T = MZ_T/MZREF
      AZ_M_1 = MZ_1/MZREF
      AZ_M_3 = MZ_3/MZREF
      AZ_M_6 = MZ_6/MZREF

