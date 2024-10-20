C
C  MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C  AIRS
C
C  MICROWAVE FIRST-GUESS RETRIEVAL
C
C  Copyright (c) 2001 Massachusetts Institute of Technology
!ROUTINE NAME: OPAC2
!CALL INTERFACE:
      SUBROUTINE opac2(BFIELD2,CBTH2,LEVELPRO,SECANT,NLEV,IST,
     + PRES,TAIR,WVCD,O3CD,ALPHA,BETA,GAMMA,dALPHAdT,dBETAdV,
     + PAVG,tstd,plb,trcoef,maxcoef,nlstd)
C
!F77  LANGUAGE-  FORTRAN 77 
!ROUTINE HISTORY:
C  VERSION- 2.3  DATE- 4/5/93  PROGRAMMER- P.ROSENKRANZ
C           2.4  4/13/93 Sung-Yung Lee
C                   Common block deleted and other cosmetic changes.
C           2.5  10/21/93    P.ROSENKRANZ - REVISED MAG.FIELD
C                      CORRECTION FOR CH.14
C           2.6  3/15/94  P. Rosenkranz - fixed bug in o2 - beta section
C           2.7  4/5/94   P. Rosenkranz - put in slant-path adjustment
C           2.8  2/28/95  P. Rosenkranz - changed detection 
C                of coefficient types (to allow for MHS substitutes)
C           2.9  Sept.18,95 P.Rosenkranz - mod. to allow variable PRES 
C                 levels
C           2.10 Sept.29, 95 P.Rosenkranz - added zero pressure as anchor
C                     point for interp. of alpha
c  C Log for /home/syl/airs_level2/src/mit/opac.f[1.4]:
C   Initial version
C [1.1] Mon Jan 29 15:39:34 1996 syl@airs.jpl saved
C   First running version
C [1.2] Thu Jul 11 14:28:34 1996 lc@airs.jpl published
C   [Thu Jul 11 14:20:20 1996] Intention for change:
C   Renamed from opac2 to opac.
C [1.3] Wed Jul 24 16:55:42 1996 lc@airs.jpl saved
C   [Tue Jul 23 16:47:28 1996] Intention for change:
C   Convert to f90.
C [1.4] Wed Aug 21 13:01:38 1996 syl@airs.jpl saved
C   [Tue Aug 20 11:48:40 1996] Intention for change:
C   August 1996 version from P. Rosenkranz.
C   Aug. 6, 1996  P.Rosenkranz - added new code for H2O lines
C   Nov. 30, 1996 P.Rosenkranz - added code for dBETAdV to H2O chan.
C   Apr. 18, 1997 P.Rosenkranz - added code for dBETAdV to O2 chan.
c   Mar. 20, 1998 P.Rosenkranz - added code for dALPHAdT to O2 chan.
C   May 22, 2001  P.Rosenkranz - removed dead code; add LEVELPRO arg,
c                 chgd name
c   Oct.25, 2001 increased dimension of local arrays
c   Aug.27, 2003 don't calculate opacity below surface
c   May 24, 2005 P. Rosenkranz - lower bound for dry-air pressure
C   Mar. 30, 2006 P.Rosenkranz - added O3 code in H2O section
C
!ABSTRACT: COMPUTES ALPHA,BETA,GAMMA VECTORS
C          CORRESPONDING TO TAIR PROFILE, FOR AMSU/MHS.
C          TRANSMITTANCE BETWEEN LEVELS WILL BE COMPUTED AS
C          EXP(-(ALPHA+BETA*WVCD+GAMMA*WLCD)).
C
!ARGUMENTS:
C    SPECIFICATIONS-
C
      IMPLICIT NONE
      LOGICAL LEVELPRO
      INTEGER NLEV
      INTEGER IST
      INTEGER MAXCOEF
      INTEGER NLSTD
      REAL TAIR(NLEV)
      REAL O3CD(NLEV)
      REAL BFIELD2
      REAL CBTH2
      REAL SECANT
      REAL WVCD(NLEV-IST)
      REAL PRES(NLEV)
      REAL ALPHA(NLEV)
      REAL BETA(NLEV-IST)
      REAL GAMMA(NLEV-IST)
      real dBETAdV(NLEV-IST)
      real dALPHAdT(NLEV)
      real PAVG(NLSTD)
      real TSTD(NLSTD)
      REAL PLB(NLSTD)
      real TRCOEF(MAXCOEF,*)
C
C  NAME       I/O    UNITS     DESCRIPTION
C
C  BFIELD2     I    GAUSS**2   SQUARE OF MAGNETIC FIELD INTENSITY AT
C                              APPROX. 50 KM ALTITUDE
C  CBTH2       I               SQUARE OF COSINE OF ANGLE BETWEEN
C                                PROPAGATION DIRECTION AND MAG. FIELD
C                     *BOTH ABOVE USED ONLY FOR HIGH CHAN. LIKE AMSU-14*
C  LEVELPRO    I               FLAG FOR FORM OF TAIR PROFILE
C  SECANT      I               SECANT OF EARTH INCIDENCE ANGLE 
C                               (FROM VERTICAL)
C  NLEV        I               NUMBER OF DEFINED ELEMENTS IN PRES AND
C                                     TAIR (<=MAXLEV)
C  IST         I               DIFFERENCE IN LENGTH OF ALPHA VS
C                              BETA AND GAMMA VECTORS (NON-NEGATIVE)
C  PRES        I     MB        PRESSURE AT LOWER BOUNDARY LEVEL OF EACH
C                            LAYER. VALUES MUST INCREASE MONOTONICALLY.
C  TAIR        I    KELVIN     ATMOSPHERIC TEMPERATURE, 
C                              AT PRES LEVELS IF LEVELPRO=.TRUE., ELSE
C                              MEAN BETWEEN LEVELS I AND I-1
C  WVCD        I    G/CM**2    H2O VAPOR MASS DENSITY BETWEEN
C                                 LEVELS I AND I-1, STARTING AT IST+1
C  O3CD        I    MOL/CM**2  OZONE NUMBER DENSITY BETWEEN
C                               LEVELS I AND I-1.
C  ALPHA       O               O2 + O3 OPACITY ALONG THE SLANT-PATH BETWEEN 
C                                 LEVELS I AND I-1
C  BETA        O    CM**2/G    H2O VAPOR CROSS SECTIONS ON SLANT-PATH
C  GAMMA       O    CM**2/G    H2O LIQUID CROSS SECTIONS ON SLANT-PATH
C                              *PRES(IST+1) IS THE LEVEL FOR WHICH 
C                              CALCULATION OF BETA AND GAMMA BEGINS*
C  dALPHAdT    O    1/K        DERIVATIVE OF ALPHA(I) W.R.TO TAIR(I)
C                                (COMPUTED ONLY FOR O2 CHANNELS).
C  dBETAdV     O    CM**4/G**2 DERIVATIVE OF BETA(I) W.R.TO WVCD(I)
C  PAVG        I    MB         PRESSURE AT MIDPOINT OF LAYERS ASSOC
C                               WITH ROWS OF TRCOEF.
C  TSTD        I    K          STANDARD ATMOSPHERE TEMPERATURE AT PAVG.
C  PLB         I    MB         PRESSURE AT LOWER BOUNDARY OF LAYERS
C                              ASSOC WITH ROWS OF TRCOEF.
C  TRCOEF      I               TRANSMITTANCE CONSTANTS FOR THE SELECTED
C                               CHANNEL.
C  MAXCOEF     I               FIRST DIMENSION OF TRCOEF ARRAY.
C  NLSTD       I              LENGTH OF PAVG, PLB AND TSTD; NO. LAYERS
C                         FOR TRCOEF (32 =< NLSTD =< MAXLEV)
C
!ROUTINES CALLED: VLINT
!PARENT: MWTRAN
!RETURN VALUES:
!FILES ACCESSED:
!DESCRIPTION: P. W. Rosenkranz, IEEE Trans. Geosci. Rem. Sens.
c  v.33, pp.1135-1140 (1995); v. 41, pp. 362-368 (2003).
!KNOWN BUGS AND LIMITATIONS: present version computes dALPHAdT only for
c  O2 channels.
!END HEADER*************************************************************
C
C  LOCAL VARIABLES
      INTEGER NLEVH,I,IWV,II,INDT,I1,MAXLEV,MAXTRLEV
      parameter (MAXLEV=120, MAXTRLEV=100)
      REAL DELTA,TEMP(0:MAXTRLEV),PLB1(0:MAXTRLEV),TEMP2(0:MAXLEV)
      REAL PRD,TAV,PAV,PRD1,DI,SECE,TH1,A,PV,PD
      LOGICAL CH14,STRATO
C  New with August 1996 version
      real S,FREQ_LINE,FREQ_DISP,BETA_LOCAL,BETA_CONT,
     &     BANDS,WIDTH
      integer IBAND,NL
      real H2OM /2.991e-23/ !MASS OF H2O MOLECULE
      REAL dBETAdS,dBETAdF,dWdS,dWdF,dPdV,DENOM,FACT,dBETAdW,dAdQ
      REAL SO3,EO3,PADJ
c**********************************************************************
C
      NLEVH = NLEV-IST
C
C   IDENTIFY COEFFICIENTS FOR O2-LINE CHANNEL:
      IF(TRCOEF(1,26).NE.0.) GOTO 30 
C***********************************************************
C  COMPUTE ALPHA,BETA,GAMMA FOR A WINDOW CHANNEL (CONTINUUM ONLY)
C  OR FOR A CHANNEL WITH UP TO SIX PASSBANDS NEAR A SINGLE H2O LINE
C  (PLUS CONTINUUM).
C
      dWdS = 0.
      dWdF = 0.
      FREQ_LINE = TRCOEF(3,26)
      PRD = 0.
      DO 10 I=1,NLEV
      TAV = TAIR(I)
      IF(I.GT.1) THEN
        IF(LEVELPRO) TAV = (TAIR(I)+TAIR(I-1))/2.
        PAV = (PRES(I)+PRES(I-1))/2.
      ELSE
        PAV = PRES(1)/1.45
      ENDIF
      PRD1 = PRD
      IWV = I -IST
      IF(IWV.GE.1 .AND. I.GT.1) THEN
       dPdV = PAV*1.573/(PRES(I)-PRES(I-1))
       PV = dPdV*WVCD(IWV)
      ELSE
       PV = 0.
       dPdV = 0.
      ENDIF
      PD = AMAX1( PAV-PV, 0. )
      PRD = AMAX1( PRES(I)-PV, 0. )
      DI = (TAV-173.)/6.
      I1 = DI
      I1 = MIN0(24,MAX0(1,I1))
      DI = DI - FLOAT(I1)
C     O2 opacity
      ALPHA(I) = (TRCOEF(1,I1)*(1.-DI) + DI*TRCOEF(1,I1+1))*
     & (PRD-PRD1)*PAV*SECANT
      dALPHAdT(I) = 0. !only to set a value
      IF(O3CD(I).LE.0.) GOTO 11
C     ozone lines
      SO3 = TRCOEF(8,I1)*(1.-DI) + DI*TRCOEF(8,I1+1)
      IF(SO3.LE.0.) GOTO 11
      PADJ = PAV*(300./TAV)**TRCOEF(8,26)
      EO3 = 1./(PADJ/TRCOEF(8,27) +1.)
      ALPHA(I) = ALPHA(I) + O3CD(I)*SO3*EO3*SECANT
11    CONTINUE
      IF(IWV.LT.1) GOTO 10
C     H2O continuum
      dBETAdF = TRCOEF(5,I1)*(1.-DI) + TRCOEF(5,I1+1)*DI
      dBETAdS = TRCOEF(6,I1)*(1.-DI) + TRCOEF(6,I1+1)*DI
      BETA_CONT = dBETAdF*PD + dBETAdS*PV
      BETA_LOCAL = 0.
      dBETAdW = 0.
      IF(FREQ_LINE.LE.0.) GOTO 14
c     local H2O line
      S = TRCOEF(2,I1)*(1.-DI) + DI*TRCOEF(2,I1+1)
      dWdF = TRCOEF(3,I1)*(1.-DI) + DI*TRCOEF(3,I1+1)
      dWdS = TRCOEF(4,I1)*(1.-DI) + DI*TRCOEF(4,I1+1)
      WIDTH = dWdF*PD + dWdS*PV
      BANDS = 0.
      DO IBAND=1,6
        FREQ_DISP = TRCOEF(3,26+IBAND)
        IF(FREQ_DISP.EQ.0.) GOTO 12
        BANDS = BANDS + 1.
        FACT = ((FREQ_LINE+FREQ_DISP)/FREQ_LINE)**2
        DENOM = FREQ_DISP*FREQ_DISP + WIDTH*WIDTH
        BETA_LOCAL = FACT*WIDTH/DENOM + BETA_LOCAL
        dBETAdW = dBETAdW + FACT*(1./DENOM -2.*(WIDTH/DENOM)**2)
      END DO
12    CONTINUE
      FACT = S/(BANDS*H2OM*3.14159)
      BETA_LOCAL = FACT*BETA_LOCAL
      dBETAdW = FACT*dBETAdW
14    CONTINUE
      BETA(IWV) = (BETA_LOCAL + BETA_CONT)*SECANT
      BETA(IWV) = AMAX1( BETA(IWV), 0.)
      dBETAdV(IWV) = (dBETAdS-dBETAdF + dBETAdW*(dWdS-dWdF))*dPdV*SECANT
      IF(I1.GE.10) THEN
       GAMMA(IWV) = (TRCOEF(7,I1)*(1.-DI) + DI*TRCOEF(7,I1+1))*SECANT
       GAMMA(IWV) = AMAX1( GAMMA(IWV), 0.)
      ELSE
       GAMMA(IWV) = 0.
      ENDIF
10    CONTINUE
      RETURN
C*********************************************************
30    CONTINUE
C  COMPUTE ALPHA FOR AN O2-BAND CHANNEL
C  interpolate TAIR profile to std p and store in TEMP
      IF(LEVELPRO) THEN
        CALL VLINT(PRES,TAIR,NLEV,PAVG,TEMP(1),NLSTD)
      ELSE
        TEMP2(1) = PRES(1)/1.45
        DO I=2,NLEV
          TEMP2(I) = (PRES(I)+PRES(I-1))/2.
        END DO
        CALL VLINT(TEMP2(1),TAIR,NLEV,PAVG,TEMP(1),NLSTD)
      ENDIF
C
c  correction for angle dependence of transmittance on o2 chann.
C  For O2-band channels, there are 25 H2O coef in position 5.
C  If there is a 26th coef, this is a correction factor for the
C  angle dependence of O2.
      DELTA = TRCOEF(5,26)
      IF(DELTA.LE.0.) DELTA = 1.
      sece = 1. + (secant - 1.) * delta
C
C  Usually coef 7 is used for liquid water and first 9 of these
C  are zero; but if nonzero, then they are mag. field coeff.
      CH14 = TRCOEF(7,9).NE.0.
C
C  Compute integrated opacity from trcoef and store in TEMP;
c   store derivative of opacity in TEMP2.
      TEMP(0) = 0.
      TEMP2(0) = 0.
      PLB1(0) = 0.
      DO 32 I=1,NLSTD
      TH1 = TSTD(I)/TEMP(I) - 1.
      dAdQ = TRCOEF(2,I) + TH1*(TRCOEF(3,I) + TH1*TRCOEF(4,I))
      A = TRCOEF(1,I) + TH1*dAdQ
C
C  COMPUTE MAGNETIC CORRECTION TO ALPHA FOR AMSU CHAN 14 OR SIMILAR CH.
      IF(CH14 .AND. PAVG(I).LE.10.) THEN
       FACT = 1. + BFIELD2*(TRCOEF(6,I)+CBTH2*TRCOEF(7,I))/TRCOEF(1,I)
       A = A*FACT
       dAdQ = dAdQ*FACT
      ENDIF
C
      TEMP2(I) = TEMP2(I-1) - dAdQ*SECE*TSTD(I)/TEMP(I)**2
      TEMP(I) = TEMP(I-1) + AMAX1(A,0.)*SECE !reuse array for opacity
      PLB1(I) = PLB(I)
      NL = I
      IF(PLB1(NL).GT.PRES(NLEV)) GOTO 33
32    CONTINUE
C
33    CONTINUE
C  Interpolate opacity back to PRES levels
      CALL VLINT(PLB1,TEMP,NL+1,PRES,ALPHA,NLEV)
      CALL VLINT(PLB1,TEMP2,NL+1,PRES,dALPHAdT,NLEV)
      STRATO = ALPHA(IST+1).GT.1.
C
C  DIFFERENCE TO GET LAYER OPACITIES
      DO 34 I=2,NLEV
      II = NLEV+2-I
      ALPHA(II) = ALPHA(II) - ALPHA(II-1)
      dALPHAdT(II) = dALPHAdT(II) - dALPHAdT(II-1)
34    CONTINUE
C
      IF(STRATO .OR. CH14) GOTO 40
C  COMPUTE BETA AND GAMMA FOR O2 CHANNEL (NOT STRATOSPHERIC)
      DO 38 I=1,NLEVH
      INDT = IST +I
      TAV = TAIR(INDT)
      IF(INDT.GT.1) THEN
        IF(LEVELPRO) TAV = (TAIR(INDT)+TAIR(INDT-1))/2.
        PAV = (PRES(INDT)+PRES(INDT-1))/2.
        dPdV = PAV*1.573/(PRES(INDT)-PRES(INDT-1))
        PV = dPdV*WVCD(I)
      ELSE
        PAV = PRES(1)/1.45
        PV = 0.
        dPdV = 0.
      ENDIF
      PD = AMAX1( PAV-PV, 0. )
      DI = (TAV-173.)/6.
      I1 = DI
      I1 = MIN0(24,MAX0(1,I1))
      DI = DI - FLOAT(I1)
      dBETAdF = TRCOEF(5,I1)*(1.-DI) + TRCOEF(5,I1+1)*DI
      dBETAdS = TRCOEF(6,I1)*(1.-DI) + TRCOEF(6,I1+1)*DI
      BETA(I) = (dBETAdF*PD + dBETAdS*PV)*SECANT
      BETA(I) = AMAX1( BETA(I), 0.)
      dBETAdV(I) = (dBETAdS-dBETAdF)*dPdV*SECANT
      IF(I1.GE.10) THEN
       GAMMA(I) = (TRCOEF(7,I1)*(1.-DI) + DI*TRCOEF(7,I1+1))*SECANT
       GAMMA(I) = AMAX1( GAMMA(I), 0.)
      ELSE
       GAMMA(I) = 0.
      ENDIF
38    CONTINUE
      RETURN
C
C  SET BETA,GAMMA TO ZERO FOR STRATOSPHERIC CHANNELS
40    DO 42 I=1,NLEVH
      BETA(I) = 0.
      dBETAdV(I) = 0.
42    GAMMA(I) = 0.
      RETURN
      END
