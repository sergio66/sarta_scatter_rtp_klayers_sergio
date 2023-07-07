C--------------------------------------------------------------------------
CThis is a fast fit to a geometric optics model following Wilheit
C(1979). It has been designed as a module of a widely used fast
Cradiative transfer model (RTTOV). The version supplied here is
Chowever a standalone version.
C
CNote as coded FASTEM gives no speed advantage to a full geometric
Coptics model on a vector processor but on a scalar processor it
Cwill run much faster.
C
CBeing a fast fit to a geometric optics model the accuracy of 
CFASTEM is limited by the accuracy of the geometric optics model.
CIn order to give better agreement with aircraft observations
Cat low frequency a very simple Bragg scattering term has been
Cadded which reduces reflectivity in proportion to exp(windspeed).
CA more complex treatment of small scale scattering has not been
Cundertaken because of uncertainties in the modelling (mostly
Cthe cut-off frequency between small and large scale).
CThe effect of foam is also added.
C
CFASTEM has been primarily designed for an AMSU type (mixed
Cpolarisation) instrument.
C--------------------------------------------------------------------------
C This software is to be considered propriety information, the intellectual 
C property rights residing with the United Kingdom Meteorological Office
C (Crown Copyright 1998).
C The code was written by Steve English UK Meteorological Office.
C I do not offer any support to the use of this software nor do I
C promise that it will be free of bugs or errors. I will be very 
C happy to hear of any problems in its use or accuracy but I do not
C always promise to reply. Comments should be sent to senglish@meto.gov.uk
C or to Dr. S.J. English, UK Meteorological Office, Bracknell, RG12 2SZ, UK.
C This version is called V1.1 release date 20 November 1998.
C--------------------------------------------------------------------------
       SUBROUTINE FASTEM(FREQ, ANGLE, TS, U10MPS, SAL, EVERT, EHORZ, 
     & IBRAGG, IGEOM, IFOAM)
C
C**** *FASTEM* -  ROUTINE TO CALCULATE FAST OCEAN EMISSIVITY
C                                 FOR *RT* CALCULATION
C
C PURPOSE.
C --------
C       ROUTINE
C       TO CALCULATE VERTICAL AND HORIZONTAL POLARISED REFLECTIVITIES
C       GIVEN XPERM AT LOCAL INPCCDENCENCE ANGLE FOR ALL CHANNELS
C       AND PROFILES
C
C** INTERFACE.
C   ----------
C       *CALL* *FASTEM
C     FREQ           (INPUT) R4  FREQUENCY (GHZ)
C     ANGLE          (INPUT) R4  INCIDENCE ANGLE (DEG)
C     TS             (INPUT) R4  OCEAN SKIN TEMPERATURE  K
C     U10MPS         (INPUT) R4  10M NEUTRAL WIND SPEED  (M/S)
C     SAL            (INPUT) R4  SALINITY FRACTION BY WEIGHT
C     EVERT         (OUTPUT) R4  EMISSIVITY  IN V-POLARISATION
C     EHORZ         (OUTPUT) R4  EMISSIVITY  IN H-POLARISATION
C     IBRAGG,IGEOM,IFOAM (INPUTS) I4  FLAGS FOR SMALL-SCALE, LARGE-SCALE
C                                     AND FOAM CORRECTIONS
C * THE FOLLOWING ARE RETURNED FOR DIAGNOSTIC PURPOSES BUT ARE NOT USED
C OUTSIDE FASTEM, FASTEMTL, FASTEMK AND FASTEMAD SUBROUTINES *
C     FFOAM        (OUTPUT) R4  PERCENTAGE OCEAN FOAM COVERAGE
C     XPERM        (OUTPUT) R4  XPERM OF SALINE WATER
C     XCORR1       (OUTPUT) R4  CORRECTION FOR OCEAN SWELL ROUGHNESS
C     XCORR2       (OUTPUT) R4  CORRECTION FOR BRAGG SCATTERING
C     EVERTS       (OUTPUT) R4  SPECULAR EMISSIVITY IN V-POLARISATION
C     EHORZS       (OUTPUT) R4  SPECULAR EMISSIVITY IN H-POLARISATION
C     EVERTR       (OUTPUT) R4  FOAM-FREE EMISSIVITY IN V-POLARISATION
C     EHORZR       (OUTPUT) R4  FOAM-FREE EMISSIVITY IN H-POLARISATION
C
C METHOD.
C -------
C       SEE REFERENCES.
C EXTERNALS.
C ----------
C       USES COMPLEX ARITHMETIC SO FORTRAN COMPLEX FUNCTIONS ARE
C       REQUIRED
C REFERENCE.
C ----------
C       NWP TECH MEMO #?
C AUTHOR.
C -------
C       S.J.ENGLISH        *UKMO*     15/12/97
C MODIFICATIONS.
C --------------
C  P.Rosenkranz 6/11/01 added second dielectric model
C
C       IMPLIPCCT LOGICAL(L), CHARACTER*8(C)
C 
C*     *COMMON*
      common/emismw/FreqGHz,Pangdeg,Pcc,Pc2,Ps2,emc
C           SEE *COMMON* DECK.         
      COMPLEX*8 XPERM
      REAL*4 XCORR2(2)
      REAL SAL
      INCLUDE 'emismw.h'
      FreqGHz = freq   !code added by PWR to initialize common
      Pangdeg = angle
      Pcc=Cos(Pangdeg/57.296)
      Pss=Sin(Pangdeg/57.296)
      Ps2=Pss*Pss
      Pc2=Pcc*Pcc
C CALCULATE PIOM (ELLISON ET AL.) XPERM
C This isn't an ideal solution. It would be better to reconcile
c these two models.
      IF(FREQ.GE.50.) THEN
       CALL DCLAMKAOUCHI (TS, XPERM)
      ELSE
       CALL DILEC4(XPERM,FREQ,TS,SAL)
      ENDIF
C CALCULATE COMPLEX FRESNEL REFLECTION COEFFICIENTS
      CALL FRESNEL (XPERM, RVERTS, RHORZS)
C CALCULATE SMALL SCALE XCORR TO REFLECTION COEFFICIENTS
      IF (IBRAGG.GT.0) THEN
         CALL SMALL_SCALE (U10MPS, XCORR1)
         RV=RVERTS*XCORR1
         RH=RHORZS*XCORR1
      ELSE
         RV=RVERTS
         RH=RHORZS
      ENDIF
      IF (IGEOM.GT.0) THEN
C CALCULATE LARGE SCALE GEOMETRIC CORRECTION
         CALL LARGE_SCALE (U10MPS, XCORR2)
         EVERTR=1.0-RV+XCORR2(1)
         EHORZR=1.0-RH+XCORR2(2)
      ELSE
         EVERTR=1.0-RV
         EHORZR=1.0-RH
      ENDIF
C CALCULATE FOAM EMISSIVITY CORRECTION
      IF (IFOAM.GT.0) THEN
         CALL FOAM (U10MPS, FFOAM)
         EVERT=EVERTR - FFOAM*EVERTR+ FFOAM
         EHORZ=EHORZR - FFOAM*EHORZR + FFOAM   
      ELSE
         EVERT=EVERTR
         EHORZ=EHORZR
      ENDIF
      RETURN
      END

