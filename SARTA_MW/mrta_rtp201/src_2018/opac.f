C
C  MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C  AIRS
C
C  MICROWAVE FIRST-GUESS RETRIEVAL
C
!ROUTINE NAME: OPAC
!CALL INTERFACE:
      SUBROUTINE opac(BFIELD2,CBTH2,CHANNO,SECANT,NLEV,IST,
     + PRES,TAIR,WVCD,ALPHA,BETA,GAMMA,dALPHAdT,dBETAdV,
     + PAVG,tstd,plb,trcoef,maxcoef,nlstd, rtatype)
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
      INTEGER*4 CHANNO
      INTEGER*4 NLEV
      INTEGER*4 IST
      INTEGER*4 MAXCOEF
      INTEGER*4 NLSTD
      REAL*4   TAIR(NLEV)
      REAL*4   BFIELD2
      REAL*4   CBTH2
      REAL*4   SECANT
      REAL*4   WVCD(NLEV-IST)
      REAL*4   PRES(NLEV)
      REAL*4   ALPHA(NLEV)
      REAL*4   BETA(NLEV-IST)
      REAL*4   GAMMA(NLEV-IST)
      real*4   dBETAdV(NLEV-IST)
      real*4   dALPHAdT(NLEV)
      real*4   PAVG(NLSTD)
      real*4   TSTD(NLSTD)
      REAL*4   PLB(NLSTD)
      real*4   TRCOEF(MAXCOEF,NLSTD)
      integer*4 rtatype ! MIT type  1=WIN,2=WAT,3=O2,4=WIN(old)
C
C  NAME       I/O    UNITS     DESCRIPTION
C
C  BFIELD2     I    GAUSS**2   SQUARE OF MAGNETIC FIELD INTENSITY AT
C                              APPROX. 50 KM ALTITUDE
C  CBTH2       I               SQUARE OF COSINE OF ANGLE BETWEEN
C                                PROPAGATION DIRECTION AND MAG. FIELD
C                                *BOTH ABOVE USED ONLY FOR CH.14*
C  CHANNO      I               INDEX FOR AMSU-A CHANNEL (1-15) 
C                              (USED ONLY IN O2-CHANNEL SECTION)
C  SECANT      I               SECANT OF EARTH INCIDENCE ANGLE 
C                               (FROM VERTICAL)
C  NLEV        I               NUMBER OF DEFINED ELEMENTS IN PRES AND
C                                     TAIR
C  IST         I               DIFFERENCE IN LENGTH OF ALPHA VS
C                              BETA AND GAMMA VECTORS (NON-NEGATIVE)
C  PRES        I     MB        PRESSURE AT LOWER BOUNDARY LEVEL OF EACH
C                            LAYER. VALUES MUST INCREASE MONOTONICALLY.
C  TAIR        I    KELVIN     ATMOSPHERIC TEMPERATURE AT PRES LEVELS
C  WVCD        I    G/CM**2    H2O VAPOR MASS DENSITY BETWEEN
C                                    LEVELS I AND I-1, STARTING AT IST+1
C  ALPHA       O               O2 OPACITY ALONG THE SLANT-PATH BETWEEN 
C                                    LEVELS I AND I-1
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
C                         FOR TRCOEF (MIN 32, MAX = DIMEN OF TEMP ARRAY)
C
!ROUTINES CALLED: VLINT
!PARENT: general purpose or AMSUTAU
!RETURN VALUES:
!FILES ACCESSED:
!DESCRIPTION: Original version: P. W. Rosenkranz, IEEE Trans. Geosci.
c  Rem. Sens. v.33, pp.1135-1140 (1995);
c revision: IGARSS'98 Digest.
!KNOWN BUGS AND LIMITATIONS: present version computes dALPHAdT only for
c  O2 channels.
!END HEADER*************************************************************
C
C  LOCAL VARIABLES
      REAL*4   DELTA(15),TEMP(0:100),PLB1(0:100),TEMP2(0:100)
      REAL*4   PRD,TAV,PAV,PRD1,DI,SECE,TH1,A,PV,PD,B,
     & DELTANEW,DELTACHAN
      LOGICAL CH14,STRATO
      INTEGER*4 NLEVH,I,IWV,II,INDT,I1
C  New with August 1996 version
      real*4   S,FREQ_LINE,FREQ_DISP,BETA_LOCAL,BETA_CONT,
     &     BANDS,WIDTH
      integer*4 IBAND
      real*4    H2OM
      REAL*4   dBETAdS,dBETAdF,dWdS,dWdF,dPdV,DENOM,FACT,dBETAdW,dAdQ

      data      H2OM /2.991E-23/ ! MASS OF H2O MOLECULE
C  Correction factors for AMSU-A:
      data delta/3*1., .9947, .9872, .9772, .9694,
     & .9750, .9756, .9150, .8928, .8966, .8826, .9108, 1./
C
      NLEVH = NLEV-IST
C
C   IDENTIFY COEFFICIENTS FOR O2-LINE CHANNEL:
c      IF(TRCOEF(1,26).NE.0.) GOTO 30 
C   IDENTIFY (OLD-STYLE) COEFFICIENTS FOR H2O-LINE CHANNEL:
c     IF(TRCOEF(2,26).NE.0.) GOTO 60

c      print *,secant

      goto(1000,1000,3000,4000) rtatype

      print 5010, CHANNO, rtatype
      i = 8
clml      call softexit('OPAC    ', i)
 
C  ********************************************************************
C  COMPUTE ALPHA,BETA,GAMMA FOR A WINDOW CHANNEL (CONTINUUM ONLY)
C  OR FOR A CHANNEL WITH UP TO SIX PASSBANDS NEAR A SINGLE H2O LINE
C  (PLUS CONTINUUM).
C  ********************************************************************

 1000 dWdS = 0.
      dWdF = 0.
      FREQ_LINE = TRCOEF(3,26)
      PRD = 0.
      DO 10 I=1,NLEV
        IF(I.GT.1) THEN
          TAV = (TAIR(I)+TAIR(I-1))/2.
          PAV = (PRES(I)+PRES(I-1))/2.
        ELSE
          TAV = TAIR(1)
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
        PD = PAV - PV
        PRD = PRES(I) - PV
        DI = (TAV-173.)/6.
        I1 = DI
        I1 = MIN0(24,MAX0(1,I1))
        DI = DI - FLOAT(I1)
        ALPHA(I) = (TRCOEF(1,I1)*(1.-DI) + DI*TRCOEF(1,I1+1))*
     &             (PRD-PRD1)*PAV*SECANT
        dALPHAdT(I) = 0. !only to set a value
        IF(IWV.LT.1) GOTO 10
        dBETAdF = TRCOEF(5,I1)*(1.-DI) + TRCOEF(5,I1+1)*DI
        dBETAdS = TRCOEF(6,I1)*(1.-DI) + TRCOEF(6,I1+1)*DI
        BETA_CONT = dBETAdF*PD + dBETAdS*PV
        BETA_LOCAL = 0.
        dBETAdW = 0.
        IF(FREQ_LINE.LE.0.) GOTO 14
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

   12   CONTINUE
        FACT = S/(BANDS*H2OM*3.14159)
        BETA_LOCAL = FACT*BETA_LOCAL
        dBETAdW = FACT*dBETAdW

   14   CONTINUE

        BETA(IWV) = (BETA_LOCAL + BETA_CONT)*SECANT
        dBETAdV(IWV) = (dBETAdS-dBETAdF + dBETAdW*(dWdS-dWdF)) *
     &                 dPdV*SECANT
        IF(I1.GE.10) THEN
         GAMMA(IWV) = (TRCOEF(7,I1)*(1.-DI) + DI*TRCOEF(7,I1+1))*SECANT
        ELSE
         GAMMA(IWV) = 0.
        ENDIF
   10 CONTINUE  ! do i = 1, NLEV
      goto 5000

C*********************************************************
C  COMPUTE ALPHA FOR AN O2-BAND CHANNEL
C
c  correction for angle dependence of transmittance on o2 chann.
C  For O2-band channels, there are 25 H2O coef in position 5.
C  If there is a 26th coef, this is a correction factor for the
C  angle dependence of O2.
C  ********************************************************************
 3000 CONTINUE

C     interpolate TAIR profile to std p and store in TEMP
      CALL VLINT(PRES,TAIR,NLEV,PAVG,TEMP(1),NLSTD)

      DELTANEW = TRCOEF(5,26)
      IF(DELTANEW.GT.0.) THEN
        DELTACHAN = DELTANEW
      ELSE
        if(CHANNO.gt.15) then
          print 3010
          i = 7
clml          call softexit('OPAC    ', i)
        endif
        DELTACHAN = DELTA(CHANNO)
      ENDIF
      sece = 1.0 + (secant - 1.0) * deltachan
c      print *,sece,secant
C
C  Usually coef 7 is used for liquid water and first 9 of these
C  are zero; but if nonzero, then they are mag. field coeff.
      CH14 = TRCOEF(7,9).NE.0.
C
C  Compute integrated opacity from trcoef and store in TEMP
      TEMP(0) = 0.0
      TEMP2(0) = 0.0
      PLB1(0) = 0.0
      DO I = 1,NLSTD
        TH1 = TSTD(I)/TEMP(I) - 1.
        dAdQ = TRCOEF(2,I) + TH1*(TRCOEF(3,I) + TH1*TRCOEF(4,I))
        A = TRCOEF(1,I) + TH1*dAdQ
C
C  COMPUTE MAGNETIC CORRECTION TO ALPHA FOR CHAN 14
        IF(CH14 .AND. PAVG(I).LE.10.) THEN
          FACT = 1.0 +
     $           BFIELD2*(TRCOEF(6,I)+CBTH2*TRCOEF(7,I))/TRCOEF(1,I)
          A = A*FACT
          dAdQ = dAdQ*FACT
        ENDIF
C
        TEMP2(I) = TEMP2(I-1) - dAdQ*SECE*TSTD(I)/TEMP(I)**2
        TEMP(I) = TEMP(I-1) + AMAX1(A,0.)*SECE
        PLB1(I) = PLB(I)
      enddo
C
C  Interpolate opacity back to PRES levels
      CALL VLINT(PLB1,TEMP,NLSTD+1,PRES,ALPHA,NLEV)
      CALL VLINT(PLB1,TEMP2,NLSTD+1,PRES,dALPHAdT,NLEV)
      STRATO = ALPHA(IST+1).GT.1.
C
C  DIFFERENCE TO GET LAYER OPACITIES
      DO I = 2,NLEV
        II = NLEV+2-I
        ALPHA(II) = ALPHA(II) - ALPHA(II-1)
        dALPHAdT(II) = dALPHAdT(II) - dALPHAdT(II-1)
      enddo
C
      IF(STRATO .OR. CH14) then
C       ********************************************************************
C       SET BETA,GAMMA TO ZERO FOR STRATOSPHERIC CHANNELS
C       ********************************************************************
        DO I=1,NLEVH
          BETA(I) = 0.0
          dBETAdV(I) = 0.0
          GAMMA(I) = 0.0
        enddo

      else

C       ********************************************************************
C       COMPUTE BETA AND GAMMA FOR O2 CHANNEL (NOT STRATOSPHERIC)
C       ********************************************************************
        DO I = 1,NLEVH
          INDT = IST +I
          IF(INDT.GT.1) THEN
            TAV = (TAIR(INDT)+TAIR(INDT-1))/2.0
            PAV = (PRES(INDT)+PRES(INDT-1))/2.0
            dPdV = PAV*1.573/(PRES(INDT)-PRES(INDT-1))
            PV = dPdV*WVCD(I)
          ELSE
            TAV = TAIR(1)
            PAV = PRES(1)/1.45
            PV = 0.0
            dPdV = 0.0
          ENDIF
          PD = PAV - PV
          DI = (TAV-173.0)/6.0
          I1 = DI
          I1 = MIN0(24,MAX0(1,I1))
          DI = DI - FLOAT(I1)
          dBETAdF = TRCOEF(5,I1)*(1.-DI) + TRCOEF(5,I1+1)*DI
          dBETAdS = TRCOEF(6,I1)*(1.-DI) + TRCOEF(6,I1+1)*DI
          BETA(I) = (dBETAdF*PD + dBETAdS*PV)*SECANT
          dBETAdV(I) = (dBETAdS-dBETAdF)*dPdV*SECANT
          IF(I1.GE.10) THEN
            GAMMA(I) = (TRCOEF(7,I1)*(1.-DI) + DI*TRCOEF(7,I1+1))*SECANT
          ELSE
           GAMMA(I) = 0.
          ENDIF
        enddo
      endif
      goto 5000

C**********************************************************
C  CODE TO USE PRE-1996 H2O COEFFICIENTS (OBSOLETE)-
C
C  FOREIGN-GAS BROADENED CONTRIBUTION TO BETA FOR AN H2O CHANNEL
C**********************************************************
 4000 CONTINUE
      CALL VLINT(PRES,TAIR,NLEV,PAVG,TEMP(1),NLSTD)
      DO I=1,NLSTD
        TH1 = TSTD(I)/TEMP(I) -1.
        B = TRCOEF(2,I) + TH1*(TRCOEF(3,I)
     &    + TH1*(TRCOEF(4,I) + TH1*TRCOEF(5,I)))
        TEMP(I) = AMAX1(B,0.)
      enddo
C     Interpolate beta back to PRES levels
      CALL VLINT(PLB,TEMP(1),NLSTD,PRES(IST+1),BETA,NLEVH)
C  
C     COMPUTE ALPHA,BETA,GAMMA FOR H2O CHANNEL
      PRD = 0.
      DO 70 I=1,NLEV
        IF(I.GT.1) THEN
          TAV = (TAIR(I)+TAIR(I-1))/2.
          PAV = (PRES(I)+PRES(I-1))/2.
        ELSE
          TAV = TAIR(1)
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
        PRD = PRES(I) - PV
        DI = (TAV-173.)/6.
        I1 = DI
        I1 = MIN0(24,MAX0(1,I1))
        DI = DI - FLOAT(I1)
        ALPHA(I) = (TRCOEF(1,I1)*(1.-DI) + DI*TRCOEF(1,I1+1))*
     &             (PRD-PRD1)*PAV*SECANT
        dALPHAdT(I) = 0. !only to set a value
        IF(IWV.LT.1) GOTO 70
C       SELF-BROADENED CONTINUUM TERM:
        dBETAdS = TRCOEF(6,I1)*(1.-DI) + DI*TRCOEF(6,I1+1)
C       ADD TO FOREIGN-GAS TERM AND MULT BY SLANT
        BETA(IWV) = (BETA(IWV)+dBETAdS*PV)*SECANT
        dBETAdV(IWV) = dBETAdS*dPdV*SECANT
        IF(I1.GE.10) THEN
         GAMMA(IWV) = (TRCOEF(7,I1)*(1.-DI) + DI*TRCOEF(7,I1+1))*SECANT
        ELSE
         GAMMA(IWV) = 0.
        ENDIF
   70 CONTINUE

 5000 RETURN
 3010 format('OPAC: exceeded dimension of delta(15)')
 5010 format('OPAC: invalid rtatype(',i3,') =',i3)
      END
