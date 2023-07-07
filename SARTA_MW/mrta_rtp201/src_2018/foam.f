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
        SUBROUTINE FOAM (U10MPS, FFOAM)
C**** *FOAM* -  MODEL OF FOAM COVERAGE FOR *RT* CALCULATION
C
C PURPOSE.
C --------
C       ROUTINE
C       TO CALCULATE FOAM COVERAGE AT A GIVEN NEUTRAL
C       STABILITY WINDSPEED. IN PRACTISE THE 10M WINDSPEED
C       IS USUALLY USED WITHOUT A STABILITY CORRECTION.
C
C** INTERFACE.
C   ----------
C       *CALL* *FOAM_CORRECTION
C       U10MPS       (INPUT)   R4  10M NEUTRAL WIND SPEED (M/S)
C       FFOAM        (OUTPUT)  R4  FOAM COVERAGE FRACTIONC
C       ALSO TAKES INPUT FROM *COMMON/(VARIABLES)/ *
C       AND *COMMON/CONSTANTS/ *.
C       AND PUTS OUTPUT IN *COMMON/(OUTPUT/ *.
C METHOD.
C -------
C       SEE REFERENCES.
C EXTERNALS.
C ----------
C       NONE.
C REFERENCE.
C ----------
C       NWP TECH MEMO #?
C
C       [1] MONAHAN E.C. AND O'MUIRCHEARTAIGH I.G., 1986: WHITECAPS
C             AND THE PASSIVE REMOTE SENSING OF THE OCEAN'S SURFACE., 
C             INT. J. OF REMOTE SENSING, 7 (NO 5), 627-642.
C
C       [2] WU S.T., 1979: OCEANIC WHITECAPS AND SEA STATE., J. PHYS.
C              OCEANOGRAPHY, 9, 1064-1078.
C
C AUTHOR.
C -------
C       S.J.ENGLISH        *UKMO*     6/11/97
C MODIFICATIONS.
C --------------
C
C       IMPLIPCCT LOGICAL(L), CHARACTER*8(C)
C 
C*     *COMMON*      SEE *COMMON* DECK.- 'emismw.h'
      common/emismw/FreqGHz,Pangdeg,Pcc,Pc2,Ps2,emc(59)
       
       FFOAM=EMC(22)*(U10MPS**EMC(23))

       RETURN
       END

