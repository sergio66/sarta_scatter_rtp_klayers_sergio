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

              WZ_T = 0
              WZ_1 = PDP
              WZ_3 = 0

c              A_O = POAMNT(L)/ROAMNT(L)
               A_O_T = 0.0
               A_O_1 = 0.0
               A_O_3 = 1/ROAMNT(L)

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
!! but PNORM == 0 at the beginning!!!!!
c              TAZ_M_T = 1/PNORM*PDP*(TR_T*A_M + TR*A_M_T)
c              TAZ_M_1 = 1/PNORM*PDP*(TR_1*A_M + TR*A_M_1)
c              TAZ_M_3 = 1/PNORM*PDP*(TR_3*A_M + TR*A_M_3)
c              TAZ_M_6 = 1/PNORM*PDP*(TR_6*A_M + TR*A_M_6)
              TAZ_M_T = 0
              TAZ_M_1 = 0
              TAZ_M_3 = 0
              TAZ_M_6 = 0
              TAZ_M_T = (TR_T*A_M + TR*A_M_T)
              TAZ_M_1 = (TR_1*A_M + TR*A_M_1)
              TAZ_M_3 = (TR_3*A_M + TR*A_M_3)
              TAZ_M_6 = (TR_6*A_M + TR*A_M_6)

            ELSE
              PDP = PRES(L)*( PRES(L) - PRES(L-1) )
CjacX         PNORM = PNORM + PDP
CjacXC
CjacXC        Note: TRZ, TOZ, and TMZ use layer-above terms
CjacX         TZ = TZ + PDP*TR
              TZ_T = PDP/RTEMP(L)
              TZ_1 = 0.0
              TZ_3 = 0.0

CjacX         TRZ = TZ/PNORM
              TRZ_T = TZ_T/PNORM   !!! this is just Layer L  being perturbed
              TRZ_T = SUM_PDP_OVER_TREF/PNORM !!! this is all layers above being perturbed, see notes book 47,
              TRZ_1 = 0
              TRZ_3 = 0
CjacXC
CjacX         A_O = POAMNT(L)/ROAMNT(L)
              A_O_T = 0.0
              A_O_1 = 0.0
              A_O_3 = 1/ROAMNT(L)
CjacXC
CjacX         TOZ = TOZ + PDP*DT*A_O
              TOZ_T = PDP*(DT_T*A_O + DT*A_O_T)
              TOZ_1 = PDP*(DT_1*A_O + DT*A_O_1)
              TOZ_3 = PDP*(DT_3*A_O + DT*A_O_3)

CjacX         TAZ_O = TOZ/PNORM
              !!! this is just Layer L  being perturbed
              TAZ_O_T = TOZ_T/PNORM
              TAZ_O_1 = TOZ_1/PNORM
              TAZ_O_3 = TOZ_3/PNORM
              !!! this is all layers above being perturbed, see notes book 47, BUT KILLS O3 JACS so forget it
!              TAZ_O_T = OZ/PNORM
!              TAZ_O_1 = 0
!              TAZ_O_3 = DTZ_O/PNORM/POAMNT(L) 

cSSM STOPPED HERE
CjacXC
CjacX         TMZ = TMZ + PDP*TR*A_M
CjacX         TAZ_M = TMZ/PNORM
              TMZ_T = PDP*(TR_T*A_M + TR*A_M_T)
              TMZ_1 = PDP*(TR_1*A_M + TR*A_M_1)
              TMZ_3 = PDP*(TR_3*A_M + TR*A_M_3)
              TMZ_6 = PDP*(TR_6*A_M + TR*A_M_6)

              TAZ_M_T = 1/PNORM*TMZ_T
              TAZ_M_1 = 1/PNORM*TMZ_1
              TAZ_M_3 = 1/PNORM*TMZ_3
              TAZ_M_6 = 1/PNORM*TMZ_6
cSSM STOPPED HERE

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
C      PWATER = KMOLE*PWAMNT(L)*PTEMP(L)/(STDDEN*STDTMP*100*DZREF(L))
      PWATER_T = PWATER/PTEMP(L)
      PWATER_1 = PWATER/PWAMNT(L)
      PWATER_3 = 0.0
c      A_F = ( 1 - PMULT*PWATER/PRES(L) )/( FX(L)*GSCAL )
      A_F_T = 1/(FX(L)*GSCAL)/PRES(L)*(-PMULT*PWATER_T)
      A_F_1 = 1/(FX(L)*GSCAL)/PRES(L)*(-PMULT*PWATER_1)
      A_F_3 = 0.0
CCjacXccc
CCjacXc for testing
CCjacXccc
CjacXC
CjacXC     Water terms
c     WZ  = WZ + PDP*PWANT(L)
      WZ_T = 0 
      WZ_1 = PDP              !!! this is just Layer L  being perturbed
      WZ_1 = WZ/PWAMNT(L)     !!! this is all layers above being perturbed, see notes book 47
      WZ_3 = 0
c     A_W = PWAMNT(L)/RWAMNT(L)
      A_W_T = 0.0
      A_W_1 = 1/RWAMNT(L)
      A_W_3 = 0.0
CjacX      WZREF = WZREF + PDP*RWAMNT(L)
CjacX      WZ = WZ + PDP*PWAMNT(L)
c      AZ_W = WZ/WZREF
      AZ_W_T = 0.0
      AZ_W_1 = WZ_1/WZREF      !!! this is just Layer L  being perturbed
      AZ_W_1 = AZ_W/PWAMNT(L)  !!! this is all layers above being perturbed, see notes book 47 !!!!! BUT KILLS O3 JACS so forget it
      AZ_W_3 = 0.0


CjacXC
CjacXC     Ozone terms
c      A_O = POAMNT(L)/ROAMNT(L)
      A_O_T = 0.0
      A_O_1 = 0.0
      A_O_3 = 1/ROAMNT(L)

CCCCC this is layer above stuff, so need to perturb ENTIRE layers 1--L, not just layer L
CjacX      XZREF = XZREF + ROAMNT(L)
CjacX      XZ = XZ + POAMNT(L)
c      XZ_O = XZ/XZREF
      XZ_O_T = 0.0
      XZ_O_1 = 0.0
      XZ_O_3 = 1/XZREF          !!! this is just Layer L  being perturbed
      XZ_O_3 = XZ_O/POAMNT(L)   !!! this is all layers above being perturbed, see notes book 47
CjacX      OZREF = OZREF + PDP*ROAMNT(L)
CjacX      OZ = OZ + PDP*POAMNT(L)
c      AZ_O = OZ/OZREF
      AZ_O_T = 0.0
      AZ_O_1 = 0.0
      AZ_O_3 = PDP/OZREF        !!! this is just Layer L  being perturbed
      AZ_O_3 = AZ_O/POAMNT(L)   !!! this is all layers above being perturbed, see notes book 47

CjacXC
CjacXC     Carbon monoxide terms
c      A_C = PCAMNT(L)/RCAMNT(L)
      A_C_T = 0
      A_C_1 = 0
      A_C_3 = 0
      A_C_5 = 1/RCAMNT(L)
CjacX      CZREF = CZREF + PDP*RCAMNT(L)
CjacX      CZ = CZ + PDP*PCAMNT(L)
c      AZ_C = CZ/CZREF
      AZ_C_T = 0.0
      AZ_C_1 = 0.0
      AZ_C_3 = 0.0
      AZ_C_5 = PDP              !!! this is just Layer L  being perturbed
      AZ_C_5 = AZ_C/PCAMNT(L)   !!! this is all layers above being perturbed, see notes book 47
CjacXC
CjacXC     Methane terms
c      A_M = PMAMNT(L)/RMAMNT(L)
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
c      AZ_M = MZ/MZREF
      AZ_M_T = MZ_T/MZREF
      AZ_M_1 = MZ_1/MZREF
      AZ_M_3 = MZ_3/MZREF
      AZ_M_6 = MZ_6/MZREF    !!! this is just Layer L  being perturbed
      AZ_M_6 = AZ_M/PMAMNT(L)   !!! this is all layers above being perturbed, see notes book 47
