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
        SUBROUTINE LARGE_SCALE(U10MPS, XCORR2)
C**** *LARGE_SCALE_CORRECTION* -  MODIFIES FRESNEL
C       REFLECTION COEFFICIENTS FOR GEOMETRIC REFLECTIONS
C PURPOSE.
C --------
C       ROUTINE: TO CALCULATE A CORRECTION TO THE FRESNEL REFLECTION COEFFICIENTS
C       ALLOWING FOR THE PRESENCE OF LARGE SCALE ROUGHNESS      
C
C** INTERFACE.
C   ----------
C       *CALL* *LARGE_SCALE_CORRECTION.F
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
C         NONE
C REFERENCE.
C ----------
C        [1] ENGLISH S.J., 1997: A FAST EMISSIVITY MODEL FOR
C              ATOVS. NWP TECHNICAL REPORT (BEING WRITTEN).
C AUTHOR.
C -------
C       S.J.ENGLISH        *UKMO*     3/11/97
C MODIFICATIONS.
C --------------
C
C       IMPLIPCCT LOGICAL(L), CHARACTER*8(C)
C C*     *COMMON*
      common/emismw/FreqGHz,Pangdeg,Pcc,Pc2,Ps2,emc(59)
C           SEE *COMMON* DECK.- 'emismw.h'
      REAL*4 ZC(6), XCORR2(2)
C EFFECTIVE SPECULAR EMISSIVITY CALCULATED
C
C CALCULATE FREQUENCY, WINDSPEED AND VIEW ANGLE SQUARED FOR
C POLYNOMIALS
      FREQGHZ2=FREQGHZ*FREQGHZ
      U10MPS2=U10MPS*U10MPS
      SEC=1.0/PCC
      SEC2=SEC*SEC
      USEC=U10MPS*SEC
C JP=1 => V POL, JP=2 => H POL 
C JP: TWO POLARISTIONS
      DO 250 JP=1,2	
C JC: SIX COEFFICIENTS (CONSTANT, U, U^2, SEC, SEC^2, U*SEC)	
         DO 260 JC=1,6
C SELECT ELEMENT FROM DATA ARRAY EMC ELEMENTS 24-59 FOR THIS MODEL
            IEMC=24+(JC-1)*3+(JP-1)*18
C COEFFICIENTS 1-5 FOR THIS POLARISATION STORED IN ZC
            ZC(JC)=EMC(IEMC)+EMC(IEMC+1)*FREQGHZ
     &                +EMC(IEMC+2)*FREQGHZ2
 260     CONTINUE
C CALCULATE CORRECTION FOR THIS POLARISATION
         XCORR2(JP)=ZC(1)+ZC(2)*SEC+ZC(3)*SEC2
     &                            +ZC(4)*U10MPS+ZC(5)*U10MPS2+ZC(6)*USEC
         XCORR2(JP)=XCORR2(JP)/100.0
 250     CONTINUE
      RETURN
      END

