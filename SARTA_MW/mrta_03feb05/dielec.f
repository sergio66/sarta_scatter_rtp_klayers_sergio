C
C  MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C  AIRS
C
C  MICROWAVE FIRST-GUESS RETRIEVAL
C
C
C------------------------------------------------------------------------------
C
C 1. DISCLAIMER
C    THE SOFTWARE AND/OR RELATED MATERIALS ARE PROVIDED "AS-IS" WITHOUT
C    WARRANTY OF ANY KIND INCLUDING ANY WARRANTIES OF PERFORMANCE OR
C    MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR PURPOSE (AS SET FORTH
C    IN UCC 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE LICENSED
C    PRODUCT, HOWEVER USED.
C
C    IN NO EVENT SHALL CALTECH/JPL BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
C    INCLUDING BUT NOT LIMITED TO INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY
C    KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C    REGARDLESS OF WHETHER CIT/JPL SHALL BE ADVISED, HAVE REASON TO KNOW,
C    OR IN FACT SHALL KNOW OF THE POSSIBILITY.
C
C    USER BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE SOFTWARE
C    AND/OR RELATED MATERIALS.
C
C 2. The user shall comply with all applicable U.S. export control laws and
C    regulations.  To the extent that the software is subject to U.S. export
C    control laws and regulations, the user has the responsibility to obtain
C    export licenses or other export authority as may be required before
C    exporting such information to foreign countries or providing access to
C    foreign persons.
C
C------------------------------------------------------------------------------
C
!ROUTINE NAME: dielec
!CALL INTERFACE:
      SUBROUTINE DIELEC (KAPPA,F,TK)
c***********************************************************************
C
!F77 LANGUAGE-  FORTRAN 77
C
!ROUTINE HISTORY:
c   Dec. 10, 2003  P. Rosenkranz, for version 4 of MW RTA
C
!ABSTRACT:  COMPUTES THE COMPLEX DIELECTRIC CONSTANT, KAPPA, FOR SYNTHETIC 
C           SEAWATER WITH SALINITY = .035
c
      IMPLICIT NONE
! ARGUMENTS:
C    SPECIFICATIONS:
      COMPLEX KAPPA
      REAL F
      REAL TK
c
C  NAME       I/O    UNITS     DESCRIPTION
c
C  KAPPA       O     none      dielectric constant
c  F           I     GHZ       frequency (validated range 3-105 GHz)
c  TK          I     KELVIN    temperature (validated range 271-303 K)
c
c
!ROUTINES CALLED: none
!PARENT: surface.F or fastem_mod.f
!RETURN VALUES:
!FILES ACCESSED:
!DESCRIPTION: 
c   W. J. Ellison et al., J.Geophys. Res. v.108, p.4663 (2003).
c   doi:10.1029/2002JD003213,2003
c
!KNOWN BUGS AND LIMITATIONS: salinity value is fixed
!END HEADER*************************************************************
C   LOCAL VARIABLES:
      REAL  T,R1,R2,C,EINF,DELTA1,DELTA2,TAU1,TAU2,SIGMA
C
      T = TK -273.15
      SIGMA = 2.906 + .09437*T
      EINF = 5.144 ! average value over -2  to 30 C
C  TAU IS IN PS
      TAU1 = 17.535 + (-.61767 + .0089481*T)*T
      TAU2 = 3.1842 + (.019189 + (-.010873 +.00025818*T)*T)*T
      DELTA1 = 68.396 + (-.40643 + (.022832 - .00053061*T)*T)*T
      DELTA2 = 4.7629 + (.1541 + (-.033717 + .00084428*T)*T)*T
      R1 = 6.2832E-3*TAU1*F
      R2 = 6.2832E-3*TAU2*F
      C = SIGMA/(5.5631E-2*F)
      KAPPA = EINF + DELTA1/CMPLX(1.,R1) + DELTA2/CMPLX(1.,R2) - 
     & CMPLX(0.,C)
      RETURN
      END
