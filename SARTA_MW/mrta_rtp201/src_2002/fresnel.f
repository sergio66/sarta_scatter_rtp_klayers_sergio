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
       SUBROUTINE FRESNEL(XPERM, RVERTS, RHORZS)
C
C**** *FRESNEL* -  FRESNEL CALCULATION FOR *RT* CALCULATION
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
C       *CALL* *FRESNEL
C
C       ALSO TAKES INPUT FROM *COMMON/(VARIABLES)/ *
C       AND *COMMON/CONSTANTS/ *.
C       AND PUTS OUTPUT IN *COMMON/(OUTPUT/ *.
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
C       S.J.ENGLISH        *UKMO*     31/10/97
C MODIFICATIONS.
C --------------
C
C       IMPLIPCCT LOGICAL(L), CHARACTER*8(C)
C 
C*     *COMMON*
      common/emismw/FreqGHz,Pangdeg,Pcc,Pc2,Ps2,emc(59)
C           SEE *COMMON* DECK. - emismw.h
       COMPLEX*8 PERM1, PERM2,RVTH,RHTH
       COMPLEX*8 XPERM
        
       PERM1 = CSQRT(XPERM - PS2)
       PERM2  = XPERM*PCC
       RHTH = (   PCC - PERM1)/(   PCC + PERM1)                     
       RVTH = (PERM2 - PERM1)/(PERM2 + PERM1)
       RVERTSR=DBLE(RVTH)
       RVERTSI=AIMAG(RVTH)
       RVERTS=RVERTSR*RVERTSR+RVERTSI*RVERTSI
       RHORZSR=DBLE(RHTH)
       RHORZSI=AIMAG(RHTH)
       RHORZS=RHORZSR*RHORZSR+RHORZSI*RHORZSI

       RETURN
       END

