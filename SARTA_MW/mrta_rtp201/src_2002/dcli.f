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
        SUBROUTINE DCLAMKAOUCHI(TS, XPERM)

C**** *DCLAMKAOUCHI* -  DIRECT MODEL OF XPERM CALCULATION FOR *RT* CALCULATION
C
C PURPOSE.
C --------
C       ROUTINE
C       TO CALCULATE XPERM OF SALINE WATER BASED ON
C       PIOM MODEL.
C
C** INTERFACE.
C   ----------
C       *CALL* *DCLAMKAOUCHI (TS,FREQGHZ)
C       FREQGHZ      (INPUT) R4  FREQUENCY OF OBSERVATION (GHZ)
C       TS       (INPUT) R4  SEA SKIN TEMPERATURE INCREMENT (K)
C	XPERM (OUTPUT)  C8  XPERM INCREMENT OF SEA WATER
C
C       ALSO TAKES INPUT FROM *COMMON/(VARIABLES)/ *
C       AND *COMMON/CONSTANTS/ *.
C       AND PUTS OUTPUT IN *COMMON/(OUTPUT/ *.
C METHOD.
C -------
C       SEE REFERENCE.
C EXTERNALS.
C ----------
C       USES COMPLEX ARITHMETIC SO FORTRAN FUNCTIONS CMPLX, CSQRT
C       AND CONJG ARE REQUIRED.
C REFERENCES.
C ----------
C       [1] LAMKAOUCHI K., BALANA A. AND ELLISON W.J., 1997:
C       NEW XPERM DATA FOR SEA WATER (30-100 GHZ).
C       DRAFT ESA REPORT.
C       
C       [2] NWP TECH MEMO #?
C AUTHOR.
C -------
C       S.J.ENGLISH        *UKMO*     31/10/97
C MODIFICATIONS.
C --------------
C       NONE
C
C       IMPLIPCCT LOGICAL(L), CHARACTER*8(C)
C 
C*     *COMMON*
      common/emismw/FreqGHz,Pangdeg,Pcc,Pc2,Ps2,emc(59)
C           SEE *COMMON* DECK. - 'emismw.h'
        COMPLEX*8 XPERM
C CONVERT FROM KELVIN TO CENTIGRATE AND DEFINE QUADRATIC AND
C CUBIC FUNCTIONS FOR LATER POLYNOMIALS
        TC=TS-273.15
        TC2=TC*TC
        TC3=TC2*TC
C CHECK INPUTS
         IF (TC.LT.-5.0.OR.TC.GT.100.0.OR.FREQGHZ.LT.10.0.OR.
     &        FREQGHZ.GT.500.0) THEN
                   WRITE (*, *) 'SEVERE WARNING FROM DCLAMKAOUCHI'
                   WRITE (*, *) FreqGHz,TS
            ELSE
            IF (FREQGHZ.LT.20.0.OR.FREQGHZ.GT.200.0) THEN
                  WRITE (*,*) 'WARNING FROM DCLAMKAOUCHI'
                  WRITE (*, *) FreqGHz,TS
            ENDIF
        ENDIF
C DEFINE TWO RELAXATION FREQUENPCCES, TAU1 AND TAU2
        TAU1=EMC(1)+EMC(2)*TC+EMC(3)*TC2
        TAU2=EMC(4)+EMC(5)*TC+EMC(6)*TC2+EMC(7)*TC3
C STATIC XPERM ESTATIC=DEL1+DEL2+EINF
        DEL1=EMC(8)+EMC(9)*TC+EMC(10)*TC2+EMC(11)*TC3
        DEL2=EMC(12)+EMC(13)*TC+EMC(14)*TC2+EMC(15)*TC3
        EINF=EMC(18)+EMC(19)*TC
C CALCULATE XPERM USING DOUBLE-DEBYE FORMULA
        FEN=2.0*EMC(20)*FREQGHZ*0.001
        FEN2=FEN**2.0
        DEN1=1.0+FEN2*TAU1*TAU1
        DEN2=1.0+FEN2*TAU2*TAU2
        PERM_REAL1=DEL1/DEN1
        PERM_REAL2=DEL2/DEN2
        PERM_IMAG1=DEL1*FEN*TAU1/DEN1
        PERM_IMAG2=DEL2*FEN*TAU2/DEN2
        PERM_REAL=PERM_REAL1+PERM_REAL2+EINF
        PERM_IMAG=PERM_IMAG1+PERM_IMAG2
        XPERM=CMPLX(PERM_REAL,PERM_IMAG)
C RETURN VALUE OF XPERM TO MAIN PROGRAMME.
        RETURN
        END

