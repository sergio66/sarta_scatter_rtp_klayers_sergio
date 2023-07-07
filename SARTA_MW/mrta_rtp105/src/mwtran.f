C  MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C  AIRS
C
C  MICROWAVE TRANSMITTANCE
C
C  Copyright (c) 2001 Massachusetts Institute of Technology
!ROUTINE NAME: mwtran
!CALL INTERFACE:
       subroutine mwtran(secant,secratio,tran,tranr,numlev,
     &           PSurf,pres,tair,levelpro,h2ovcd,h2olcd,
     &           alat,alon,nchan)
c***********************************************************************
C
!F77 LANGUAGE-  FORTRAN 77
!ROUTINE HISTORY:
c June 4, 2001  programmer: P. Rosenkranz 
C Oct. 24, 2001   increased MAXLEV
c           
!ABSTRACT: COMPUTES LAYER TRANSMITTANCE VECTOR CORRESPONDING
C             TO TAIR AND H2O VAPOR AND LIQUID PROFILES
c             This is a more generally useful alternative to amsutau.
C
!ARGUMENTS:
C    SPECIFICATIONS:
C
      implicit none
      integer numlev, nchan
      logical levelpro
      real secant,secratio,tran(numlev),tranr(numlev)
      real tair(numlev)
      real h2ovcd(numlev),h2olcd(numlev),alat,alon,pres(numlev)
C
C  NAME       I/O    UNITS     DESCRIPTION
C
C  SECANT      I               SECANT OF EARTH INCIDENCE ANGLE 
c                               (FROM VERTICAL).
C  SECRATIO    I               RATIO OF REFLECTED PATH LENGTH TO DIRECTLY
C                               OBSERVED PATH LENGTH.
C  TRAN        O               LAYER TRANSMITTANCES PRES(I-1) TO PRES(I)
C                              ALONG THE SLANT PATH OF SECANT.
C  TRANR       O               LAYER TRANSMITTANCES PRES(I-1) TO PRES(I)
C                              ALONG THE SLANT PATH OF SECANT*SECRATIO.
C  NUMLEV      I               NUMBER OF LEVELS INCLUDING THE SURFACE,
C                                     IN PRES AND TAIR (<=MAXLEV).
C  psurf       I    MB         PRESSURE AT SURFACE
C  PRES        I    MB         PRESSURE AT LOWER BOUNDARY LEVEL OF EACH
C                                LAYER. LAST VALUE IS SURFACE PRESSURE.
C  TAIR        I    KELVIN     ATMOSPHERIC TEMPERATURE AT PRES LEVELS (IF
C                               LEVELPRO) OR MEAN TEMPERATURE OF LAYERS
C                               PRES(I-1) TO PRES(I) (IF LEVELPRO=.FALSE.)
C  LEVELPRO    I               FLAG SPECIFYING FORM OF TEMPERATURE PROFILE
C  H2OVCD      I    MOL/CM**2  H2O VAPOR NUMBER DENSITY BETWEEN
C                               LEVELS I AND I-1.
C  H2OLCD      I    MOL/CM**2  H2O LIQUID NUMBER DENSITY BETWEEN
C                               LEVELS I AND I-1.
C  ALAT        I     DEG N     LATITUDE OF OBSERVATION
C  ALON        I     DEG E     LONGITUDE OF OBSERVATION
C                   ABOVE 2 USED ONLY FOR AMSU CH.14 OR SIMILAR MAG. FIELD
C                   DEPENDENCE.  IF ABS(ALAT)>90, TRANSMITTANCE
C                   IS COMPUTED FOR ZERO FIELD CASE.
C  NCHAN       I    INDEX FOR CHANNEL COEFFICIENTS IN THE COMMON BLOCK.
C
!ROUTINES CALLED: OPAC2, BFIELD
!PARENT: general purpose
!RETURN VALUES:
!FILES ACCESSED:
!DESCRIPTION: Original version: P. W. Rosenkranz, IEEE Trans. Geosci.
c  Rem. Sens. v.33, pp.1135-1140 (1995);
c revision: IGARSS'98 Digest.
!KNOWN BUGS AND LIMITATIONS:
!END HEADER*************************************************************
c      INCLUDE 'mwtran.inc'
      INCLUDE 'amsurta.com'
C  Local variables
      integer  maxlev
      parameter (MAXLEV=120)
      real    h2om
      parameter (h2om = 2.991e-23)
c
      integer ist, nlevh, iwv, i, j
      real    bx, by, bz, bfield2, cbth2, op
      real alpha(MAXLEV),beta(MAXLEV),gamma(MAXLEV),dbetadv(maxlev)
      real dalphadt(MAXLEV)
      real*4 wvcd(maxlev)       ! Water vapor column density in g/cm**2
      real*4 lqcd(MAXLEV)   ! Liquid Water column density in gr/cm**2
c
c for fix at surface
      real*4   ratio,psurf
c
c***********************************************************************
      if(nchan.gt.nchread .or. numlev.gt.MAXLEV) then
       do i=1,numlev
       tran(i) = 0.
       tranr(i) = 0.
       end do
      return
      endif
c
      rtatype(nchan) = 1                          ! assume window or new H2O
      if(TRCOEF(1,26,nchan).ne.0.0) then
         rtatype(nchan) = 3  ! O2 channel
      endif 

      ist = 0
      nlevh = numlev-ist
      do 10 i=1,nlevh
      j = i + ist
      lqcd(i) = h2olcd(j)*h2om
10    wvcd(i) = h2ovcd(j)*h2om
C
C  GET MAGNETIC FIELD VALUES: only if magnetic field correction is needed.
C
      IF(TRCOEF(7,9,NCHAN).NE.0. .and. ABS(ALAT).LE.90.) THEN
       CALL BFIELD(ALAT,ALON,6413.,BX,BY,BZ)
       BFIELD2 = (BX*BX + BY*BY + BZ*BZ)*1.E-10
      ELSE
       BFIELD2 = 0.
      ENDIF
      CBTH2 = .5 !not critical for AMSU
c
c    The new version of opac computes the derivatives dalphadt and 
c    dbetadv, but they are not used here.
C
c       call opac(bfield2,cbth2,nchan,secant,numlev,ist,
c     &     pres, tair, wvcd, alpha, beta, gamma,dalphadt,dbetadv,
c     &     pstd, tstd, plb, trcoef(1,1,nchan),
c     &     MAXRTACOEF, ntrlev, rtatype(nchan))
c
       call opac2(bfield2,cbth2,levelpro,secant,numlev,ist,
     &     pres,tair,wvcd,alpha,beta,gamma,dalphadt,dbetadv,
     &     pstd,tstd,plb,trcoef(1,1,nchan),maxcoef,ntrlev)
c
c Phil's original code ignores the thin layer at the surface.
c  So added the next little bit...
c
      ratio = (Psurf-pres(numlev-1))
     &                  / (pres(numlev)-pres(numlev-1))
      alpha(numlev) = alpha(numlev) * ratio
      beta(numlev)  = beta(numlev) * ratio
      gamma(numlev) = gamma(numlev) * ratio
c Phil's code
      do i = 1, numlev
         op = alpha(i)
         if(i.gt.ist) then
           iwv    = i-ist
           op = op + wvcd(iwv)*beta(iwv) + h2olcd(i)*h2om*gamma(iwv)
         end if
         tran(i) = DBLE(exp(-op))
         tranr(i) = DBLE(exp(-op*secratio))
      end do

      return
      end
