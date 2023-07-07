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
        SUBROUTINE SMALL_SCALE(U10MPS, XCORR1)
C**** *SMALL_SCALE* -  MODIFIES FRESNEL
C       REFLECTION COEFFICIENTS FOR SMALL SCALE RIPPLES.
C PURPOSE.
C --------
C       TO CALCULATE A CORRECTION TO THE FRESNEL REFLECTION COEFFICIENTS
C       ALLOWING FOR THE PRESENCE OF SMALL RIPPLES.
C
C** INTERFACE.
C   ----------
C       *CALL* *SMALL_SCALE_CORRECTION.F
C       FREQGHZ      (INPUT) R4  FREQUENCY OF OBSERVATION (GHZ)
C       U10MPS       (INPUT) R4  10M NEUTRAL WIND SPEED (M/S)
C       PC2         (INPUT) R4  COSINE^2 LOCAL ZENITH ANGLE
C       XCORR1     (OUTPUT)  R4  REFLECTIVITY CORRECTION   
C
C       ALSO TAKES INPUT FROM *COMMON/(VARIABLES)/ *
C       AND *COMMON/CONSTANTS/ *.
C       AND PUTS OUTPUT IN *COMMON/(OUTPUT/ *.
C METHOD.
C -------
C       SEE REFERENCES.
C EXTERNALS.
C ----------
C       NONE
C REFERENCE.
C ----------
C       [1] ENGLISH S.J., 1997: A FAST EMISSIVITY MODEL FOR
C              ATOVS. NWP TECHNICAL REPORT (BEING WRITTEN).
C
C       [2] WU, S. T., AND FUNG A.K., 1972:  A NONCOHERENT MODEL
C	       FOR MICROWAVE EMISSION AND BACKSCATTERING FROM THE SEA
C              SURFACE. J.G.R., 77, 5917-5929.
C
C       [3] GUILLOU C., 1996: ETUDE DU TRANSFERT RADIATIF DANS
C              L'ATMOSPHERE EN MICRO-ONDES, A PARTIR D'OBSERVATIONS
C              RADIOMETRIQUES AEROPORTEES. THESE, UNIVERSITE DE PARIS.
C AUTHOR.
C -------
C       S.J.ENGLISH        *UKMO*     6/11/97
C MODIFICATIONS.
C --------------
C
C       IMPLIPCCT LOGICAL(L), CHARACTER*8(C)
C 
C*     *COMMON*        SEE *COMMON* DECK. - 'emismw.h'
      common/emismw/FreqGHz,Pangdeg,Pcc,Pc2,Ps2,emc(59)
        
      IF (FREQGHZ.GT.0.1) THEN
         XCORR1=EXP(EMC(21)*U10MPS*PC2/(FREQGHZ*FREQGHZ))
      ELSE
         XCORR1=1.0
      ENDIF
      RETURN
      END

