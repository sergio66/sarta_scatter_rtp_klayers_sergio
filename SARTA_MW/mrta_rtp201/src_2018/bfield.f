C  %Z%SCCS: %G%  FILE: %M% (%I%)
C
C  MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C  AIRS
C
C  MICROWAVE FIRST-GUESS RETRIEVAL
C
!ROUTINE NAME: BFIELD
!CALL INTERFACE:
      SUBROUTINE BFIELD(ALAT,ALON,R,BX,BY,BZ)
C
!F77  LANGUAGE-  FORTRAN 77
!ROUTINE HISTORY:
C  VERSION- 1.1 (3RD ORDER)  DATE- 4/20/93    PROGRAMMER- P.ROSENKRANZ
C  VERSION- 1.2 (5TH ORDER)  DATE- 12/21/94    PROGRAMMER- P.ROSENKRANZ
C
c  C Log for /moria/lc/src/mit/bfield.F[1.0]:
C   Initial version
C  Apr. 24, 1997 PWR mod. for Proto 5.1
C
!ABSTRACT: COMPUTE MAGNETIC FIELD COMPONENTS, IN NANOTESLA.
C
!ARGUMENTS:
C
C  NAME     TYPE    I/O   UNITS    DESCRIPTION
C
C  ALAT     R*4      I    DEG N    GEOCENTRIC LATITUDE
C  ALON     R*4      I    DEG E    GEOCENTRIC LONGITUDE
C  R        R*4      I    KM       DISTANCE FROM CENTER OF EARTH
C  BX       R*4      O    NT       NORTHWARD COMPONENT OF FIELD
C  BY       R*4      O    NT       EASTWARD COMPONENT OF FIELD
C  BZ       R*4      O    NT       DOWNWARD COMPONENT OF FIELD
C
!ROUTINE CALLED: SHVAL3
!PARENT: PROFIL3 or AMSUTAU
!RETURN VALUES:
!FILES ACCESSED: none
!DESCRIPTION:
!KNOWN BUGS AND LIMITATIONS: ADEQUATE FOR APPLICATIONS NOT REQUIRING 
c HIGH ACCURACY OR TRACKING OF THE SECULAR VARIATION.
!END HEADER*************************************************************
C
      IMPLICIT NONE
      INTEGER NDIM
      PARAMETER (NDIM=35)
      REAL GH(NDIM),EXT(3),ALAT,ALON,R,BX,BY,BZ
      INTEGER IEXT,NMAX
C  G(M,N), H(M,N) (EPOCH 1995.)
      DATA GH/  -29685.,  -1798.,   5330.,  -2200.,
     &  3070.,  -2357.,   1693.,   -449.,
     & 1331.,  -2274.,   -265.,   1246.,
     &  301.,    777.,   -402.,    941.,
     & 785.,    261.,    289.,   -231.,
     & -420.,    103.,    114.,   -306.,
     & -208.,    352.,     47.,    236.,
     & 156.,   -126.,   -152.,   -166.,
     & -61.,    -26.,    100./
      IEXT = 0
      NMAX = 5
      CALL SHVAL3(ALAT,ALON,R,NMAX,GH,IEXT,EXT,BX,BY,BZ)
      RETURN
      END
