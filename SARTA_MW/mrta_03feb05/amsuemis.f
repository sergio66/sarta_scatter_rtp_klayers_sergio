      SUBROUTINE AMSUEMIS(FREQ,SECANT,SCANANG,QV,TSURF,WIND,
     & EMISOCEAN,SECANTRATIO)
!F77 LANGUAGE-  FORTRAN 77
!ROUTINE HISTORY:
c 2/2/05 initial version
c
!ABSTRACT: Computes ocean-surface emissivity for AMSU or a similar instrument
!ARGUMENTS:
C    SPECIFICATIONS:
      implicit none
C
C     TYPE  NAME          I/O  UNITS     DESCRIPTION
      REAL  FREQ        !  I   GHz       FREQUENCY
      REAL  SECANT      !  I             SECANT OF EARTH INCIDENCE ANGLE 
      REAL  SCANANG     !  I   degrees   INSTRUMENT SCAN ANGLE
      LOGICAL QV        !  I             .true. for quasi-vertical polarization;
c                                        i.e., near-vert pol. for small scan angles.
c                                        .false. for quasi-horizontal pol.
      REAL  TSURF       !  I   KELVIN    SURFACE TEMPERATURE
      REAL  WIND        !  I   m/sec     NEAR-SURFACE WIND SPEED
      REAL  EMISOCEAN   !  O             OCEAN-SURFACE EMISSIVITY
      REAL  SECANTRATIO !  O             RATIO OF DOWNWELLING PATHLENGTH
C                                         TO DIRECT PATHLENGTH
!ROUTINES CALLED: DIELEC
!PARENT: M_RTA
!DESCRIPTION: see Rosenkranz and Barnet, J. Geophys. Res. (2005)
!KNOWN BUGS AND LIMITATIONS: The surface roughness effect was derived 
c from AMSU data in a 705-km orbit and didn't distinguish between 
c dependence on incidence angle and dependence on polarization.
!END HEADER*************************************************************
C  Local variables
      REAL RV,RH,S2ALPHA,CALPHA,S2SCAN,CSCAN
      REAL EMISSP,S,FREQ1,F1,T0,T1,TINF,DELTA,DELTASECANT
      COMPLEX KAPPA,Z,RV1,RH1
      REAL SCOEF /.0061/
      REAL T0COEF(2) /0.1199, 0.1844/
      REAL TICOEF(2) /0.2795, 0.3863/
c***********************************************************************
C  compute specular ocean emissivity 
      s2scan = sin(SCANANG/57.296)**2
      CALPHA = 1./SECANT ! cosine of incidence angle
      S2ALPHA = 1. - CALPHA*CALPHA
      CALL DIELEC(KAPPA,FREQ,TSURF)
      Z=CSQRT(KAPPA-S2ALPHA)
      RV1=(KAPPA*CALPHA - Z)/(KAPPA*CALPHA + Z)
      RH1=(CALPHA - Z)/(CALPHA + Z)
      RV = REAL(RV1)**2 +AIMAG(RV1)**2
      RH = REAL(RH1)**2 +AIMAG(RH1)**2
      IF(QV) THEN
!      fraction of horizontal pol.= sin(scanangle)^2
       EMISSP = 1. - RH*S2SCAN - RV*(1.-S2SCAN)
      ELSE
!      fraction of vertical pol.= sin(scanangle)^2
       EMISSP = 1. - RV*S2SCAN - RH*(1.-S2SCAN)
      ENDIF
c
C  increased surface brightness due to roughness
      S = 1.5
      FREQ1 = 50.
      F1 = 1./( 1.+(FREQ1/FREQ)**S )
      CSCAN = COS(SCANANG/57.296)
      T0 = WIND*( T0COEF(1) + CSCAN*T0COEF(2) )
      TINF = WIND*( TICOEF(1) + CSCAN*TICOEF(2) )
      T1 = TINF - T0
      DELTA = T0 + T1*F1
      EMISOCEAN = EMISSP + DELTA/TSURF
C
      DELTASECANT = .02 + WIND*SCOEF
      SECANTRATIO = (SECANT+DELTASECANT)/SECANT
      RETURN
      END
