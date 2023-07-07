      SUBROUTINE DILEC4 (KAPPA,F,TK,SL)
C   COMPUTES THE COMPLEX DIELECTRIC CONSTANT, KAPPA, FOR SEAWATER
      IMPLICIT NONE
      COMPLEX KAPPA
C   INPUT ARGUMENTS:
      REAL F  ! FREQUENCY IN GHZ
      REAL TK ! KELVIN TEMPERATURE
      REAL SL ! SALINITY FRACTION BY WEIGHT 
C               (FOR SEAWATER, AVERAGE SL = .035)
C   11/18/98 P.Rosenkranz 
c
C   REFERENCE FOR EQUATIONS:
c   C. Guillou et al., Radio Science v.33, pp.649-667 (1998).
c
C   LOCAL VARIABLES:
      REAL  S,T,R,C,EINF,ESTAT,TAU,SIGMA
C
      S = SL*1000.
      T = TK -273.15
      SIGMA = .08637 + (.0306 - 4.121E-4*T)*T
     & + S*(.07745 + (.001687 + 1.937E-5*T)*T)
      ESTAT = 81.82 + (-6.05E-2 + (-3.166E-2 + (3.109E-3 
     & + (-1.179E-4 + 1.483E-6*T)*T)*T)*T)*T 
     & - S*(.1254 + (9.403E-3 + (-9.555E-4 +(9.088E-5
     & + (-3.601E-6 + 4.713E-8*T)*T)*T)*T)*T)
      EINF = 6.458 + (-.04203 + (-.006588 + (6.492E-4
     & + (-1.2328E-5 + 5.043E-8*T)*T)*T)*T)*T
C  TAU IS IN PS
      TAU = 17.303 + (-.6665 + (5.148E-3 + (1.2145E-3
     & + (-5.032E-5 + 5.827E-7*T)*T)*T)*T)*T
     & + S*(-6.272E-3 + (2.357E-4 + (5.075E-4 + (-6.398E-5
     & + (2.463E-6 - 3.066E-8*T)*T)*T)*T)*T)
      R = 6.2832E-3*TAU*F
      C = SIGMA/(5.5633E-2*F)
      KAPPA = EINF + (ESTAT - EINF)/CMPLX(1.,R) - CMPLX(0.,C)
      RETURN
      END
