C  MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C  AIRS
C
C  MICROWAVE TRANSMITTANCE
C
C  Copyright (c) 2001 Massachusetts Institute of Technology
!ROUTINE NAME: mwtran
!CALL INTERFACE:
      subroutine mwtran(secant,secratio,tran,tranr,
     &                   numlev,pres,tair,levelpro,h2ovcd,h2olcd,o3cd,
     &                   alat,alon,nchan)
c***********************************************************************
C
!F77 LANGUAGE-  FORTRAN 77
!ROUTINE HISTORY:
c June 4, 2001  programmer: P. Rosenkranz 
C Oct. 24, 2001   increased MAXLEV
c Aug. 27, 2003   Lambertian surface reflection option
c Mar. 29, 2006   add O3 argument
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
      real o3cd(numlev)
C
C  NAME       I/O    UNITS     DESCRIPTION
C
C  SECANT      I               SECANT OF EARTH INCIDENCE ANGLE 
c                               (FROM VERTICAL).
C  SECRATIO    I               RATIO OF REFLECTED PATH LENGTH TO DIRECTLY
C                               OBSERVED PATH LENGTH; OR SET TO ZERO TO
C                               USE A LAMBERTIAN APPROX.
C  TRAN        O               LAYER TRANSMITTANCES PRES(I-1) TO PRES(I)
C                              ALONG THE SLANT PATH OF SECANT.
C  TRANR       O               LAYER TRANSMITTANCES PRES(I-1) TO PRES(I)
C                              ALONG THE SLANT PATH OF SECANT*SECRATIO.
C  NUMLEV      I               NUMBER OF LEVELS IN PRES AND TAIR (<=MAXLEV).
C  PRES        I    MB         PRESSURE AT LOWER BOUNDARY LEVEL OF EACH
C                                LAYER. LAST VALUE IS SURFACE PRESSURE.
C  TAIR        I    KELVIN     ATMOSPHERIC TEMPERATURE AT PRES LEVELS (IF
C                               LEVELPRO=.TRUE.) OR MEAN TEMPERATURE OF LAYERS
C                               PRES(I-1) TO PRES(I) (IF LEVELPRO=.FALSE.)
C  LEVELPRO    I               FLAG SPECIFYING FORM OF TAIR PROFILE
C  H2OVCD      I    MOL/CM**2  H2O VAPOR NUMBER DENSITY BETWEEN
C                               LEVELS I AND I-1.
C  H2OLCD      I    MOL/CM**2  H2O LIQUID NUMBER DENSITY BETWEEN
C                               LEVELS I AND I-1.
C  O3CD        I    MOL/CM**2  OZONE NUMBER DENSITY BETWEEN
C                               LEVELS I AND I-1.
C  ALAT        I     DEG N     LATITUDE OF OBSERVATION
C  ALON        I     DEG E     LONGITUDE OF OBSERVATION
C                   ALAT,ALON USED ONLY FOR AMSU CH.14 OR SIMILAR MAG. FIELD
C                   DEPENDENCE.  IF ABS(ALAT)>90, TRANSMITTANCE
C                   IS COMPUTED FOR ZERO FIELD CASE.
C  NCHAN       I    INDEX FOR CHANNEL COEFFICIENTS IN THE COMMON BLOCK.
C
!ROUTINES CALLED: OPAC2, BFIELD
!PARENT: general purpose
!RETURN VALUES:
!FILES ACCESSED:
!DESCRIPTION: P. W. Rosenkranz, IEEE Trans. Geosci. Rem. Sens.
c  v.33, pp.1135-1140 (1995); v. 41, pp. 362-368 (2003).
!KNOWN BUGS AND LIMITATIONS:
!END HEADER*************************************************************
      INCLUDE 'mwtran.inc'
C  Local variables
      integer ist, nlevh, iwv, i, j, maxlev
      real    h2om, bx, by, bz, bfield2, cbth2, op
      parameter (MAXLEV=120)
      real alpha(MAXLEV),beta(MAXLEV),gamma(MAXLEV),dbetadv(maxlev)
      real dalphadt(MAXLEV)
      real wvcd(maxlev)       ! Water vapor column density in g/cm**2
      real totop1,secref,op1(MAXLEV)
c***********************************************************************
      if(nchan.gt.nchread .or. numlev.gt.MAXLEV) then
       do i=1,numlev
       tran(i) = 0.
       tranr(i) = 0.
       end do
      return
      endif
c
      h2om = 2.991e-23
      ist = 0
      nlevh = numlev-ist
      do 10 i=1,nlevh
      j = i + ist
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
      call opac2(bfield2,cbth2,levelpro,secant,numlev,ist,
     & pres,tair,wvcd,o3cd,alpha,beta,gamma,dalphadt,dbetadv,pstd,
     & tstd,plb,trcoef(1,1,nchan),maxcoef,ntrlev)
      TOTOP1 = 0.
      do i = 1, numlev
         op = alpha(i)
         if(i.gt.ist) then
           iwv    = i-ist
           op = op + wvcd(iwv)*beta(iwv) + h2olcd(i)*h2om*gamma(iwv)
         end if
         tran(i) = exp(-op)
         OP1(I) = OP/SECANT  !nadir opacity
         TOTOP1 = TOTOP1 + OP1(I)
      END DO
      IF(SECRATIO.LE.0.) THEN ! approximate Lambertian
        SECREF = 1.55 - .16*ALOG(TOTOP1 + .06)
        SECREF = AMAX1(SECREF, 1.)
      ELSE ! use quasi-specular model
        SECREF = SECANT*SECRATIO
      ENDIF
      do i = 1, numlev
        TRANR(I) = EXP(-OP1(I)*SECREF)
      end do
      return
      end
