C
C  MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C
C  AIRS
C
C  MICROWAVE FIRST-GUESS RETRIEVAL
C
!ROUTINE NAME: amsutau
!CALL INTERFACE:
      subroutine amsutau(nchan, secant, hsbflg,
     &      numlev, Psurf, Pobs, Tair, h2ovcd, h2olcd,
     &      alat, alon, taufix, tauh2o, tauliq,
     &      taufix_r, tauh2o_r, tauliq_r)
      implicit none
      include 'amsurta.com'
      include 'paramet.com'

c***********************************************************************
C
!F77 LANGUAGE-  FORTRAN 77
!ROUTINE HISTORY:
C  VERSION- 2.0  DATE- 2/22/93  PROGRAMMER- P.ROSENKRANZ
C  VERSION- 2.1  DATE- 4/06/93  PROGRAMMER- Sung-Yung Lee
C                                   Commmon block deleted in opac.f,
C                                   and other cosmetic changes.
C           2.2        10/21/93 P.ROSENKRANZ -
C                                changed tau_ozo to tau_liq
C                                changed dimensions of local arrays
C           2.3        3/16/94 P. Rosenkranz - error occuring
c                      when ist=0; also fixed underflow by using double
c                      precision.
c           2.4        4/5/94  P. Rosenkranz - correction for
c                        angles away from vertical is now done in opac.
c           2.5        2/27/95 P.Rosenkranz - changed call to opac.
c           2.6        Sept.19,95 P.Rosenkranz- new version of OPAC
c                     removes IST restriction and allows any PRES array
c  C Log for /home/jsg/airs_level2_f77/src/mit/amsutau.F[1.2]:
C   Initial version
C [1.1] Mon Jan 29 15:39:28 1996 syl@airs.jpl saved
C   First running version
C [1.2] Thu Jul 11 14:44:18 1996 lc@airs.jpl published
C   [Thu Jul 11 14:19:06 1996] Intention for change:
C   Installed new MW rapid algorithm.
C   Added include files for transmission coefficient data.
C   Changed argument list in the subroutine statement.
C   Removed local parameter variable MAXLEV.
C   Renamed array variable PRES to POBS.
C   Call opac instead opac2.
C Apr. 18, 1997 pwr - changed call to OPAC for Proto 5.1 version.
c June 24, 1997 pwr - get tstd from common trh.
c Sept 10, 1997 syl - take care of thin layer at the surface
c Mar. 18, 1998 pwr - new call for OPAC
c
!ABSTRACT: COMPUTES TRANSMITTANCE VECTOR CORRESPONDING
C             TO TAIR AND H2O VAPOR AND LIQUID PROFILES
C             FOR AMSU/MHS. IF ABS(ALAT)>90, TRANSMITTANCE
C             IS COMPUTED FOR ZERO FIELD CASE.
C
!ARGUMENTS:
C    SPECIFICATIONS:
C
C
C  NAME       I/O    UNITS     DESCRIPTION
C
C  NCHAN       I               INDEX FOR AMSU/MHS CHANNEL (1-20)
C                                1-15 for AMSU and 16-20 for MHS
C  SECANT      I               SECANT OF EARTH INCIDENCE ANGLE
c                               (FROM VERTICAL)
c  taufix      R               transmittance to space for fixed gases
c  tauh2o      R               transmittance to space for water vapor
c  tauliq      R               transmittance to space for liquid water
C  NUMLEV      I               NUMBER OF LEVELS ABOVE THE SURFACE (IN
C                                     pobs AND TAIR)
C  psurf       I    MB         PRESSURE AT SURFACE
C  pobs        I    MB         PRESSURE AT LOWER BOUNDARY LEVEL OF EACH
C                                 LAYER
C  TAIR        I    KELVIN     ATMOSPHERIC TEMPERATURE AT PRES LEVELS
C  H2OVCD      I    MOL/CM**2  H2O VAPOR NUMBER DENSITY BETWEEN
C                               LEVELS I AND I-1
C  H2OLCD      I    MOL/CM**2  H2O LIQUID NUMBER DENSITY BETWEEN
C                               LEVELS I AND I-1
C  ALAT        I     DEG N     LATITUDE OF OBSERVATION
C  ALON        I     DEG E     LONGITUDE OF OBSERVATION
c  Pstd_a     COM              pressure level within layer
c  Tstd_a     COM              standard temperature profile at Pstd(L)
c  trcoef     COM              rapid transmittance coefficients for channel NCHAN
C
!ROUTINES CALLED: OPAC,BFIELD
!PARENT: general purpose
!RETURN VALUES:
!FILES ACCESSED:
!DESCRIPTION: Original version: P. W. Rosenkranz, IEEE Trans. Geosci.
c  Rem. Sens. v.33, pp.1135-1140 (1995);
c Revision: IGARSS'98 Digest.
!KNOWN BUGS AND LIMITATIONS:
!END HEADER*************************************************************
C
C  input variables

      real*4    secant,Tair(*)
      real*4    h2ovcd(*),h2olcd(*),alat,alon,psurf,pobs(*)
      integer*4 nchan,numlev

C  output variables
      real*8    taufix(*), tauh2o(*), tauliq(*)
      real*8    taufix_r(*), tauh2o_r(*), tauliq_r(*)

C  Local variables
      real*4    h2om
      parameter (h2om=2.991e-23)    ! watmol/avogad -This is the value

      integer*4 ist, nlevh, iwv, i, j
      logical*4 hsbflg
      real*4    bx, by, bz, bfield2, cbth2
c      used in computation of the mw transmittance coefficients.
      real*4   alpha(MAXLEV),beta(MAXLEV)
      real*4   gamma(MAXLEV),dbetadv(MAXLEV)
      real*4   dalphadt(MAXLEV)
      real*4   wvcd(MAXLEV)       ! Water vapor column density in gr/cm**2
      real*4   lqcd(MAXLEV)       ! Liquid Water column density in gr/cm**2
      real*8   trnliq, trnfix, trnh2o ! summation variables
      real*4   ratio
      real*4   dfact   ! downwelling optical depth parameter
c
c      logical levelpro
C
c***********************************************************************

      if(nchan.gt.MAXRTACHAN) then
        print *, 'channel index out of bounds', nchan, maxchan
        stop
      endif

      if(numlev.gt.MAXLEV) then
        print *, '# levels index is out of bounds'
        print *, 'numlev=', numlev, '  MAXLEV=',MAXLEV
        stop
      endif

      rtatype(nchan) = 1                          ! assume window or new H2O
      if(TRCOEF(1,26,nchan).ne.0.0) then
         rtatype(nchan) = 3  ! O2 channel
      endif 

      ist = 0
      nlevh = numlev-ist
      do i = 1,nlevh
        j = i + ist
        wvcd(i) = h2ovcd(j)*h2om  ! convert molecules/cm^2 to  gram/cm^2
        lqcd(i) = h2olcd(j)*h2om
      enddo
C
C  GET MAGNETIC FIELD VALUES:
C          only channel 14 (AMSU) needs magnetic field correction.
C          which is channel 15 (ATMS)
C
      IF(TRCOEF(7,9,NCHAN).NE.0. .and. ABS(ALAT).LE.90.) THEN
       CALL BFIELD(ALAT,ALON,6413.,BX,BY,BZ)
       BFIELD2 = (BX*BX + BY*BY + BZ*BZ)*1.E-10
      ELSE
       BFIELD2 = 0.
      ENDIF
c
c      IF(magchl(nchan) .and. ABS(ALAT).LE.90.) THEN
cC                  What is 6413.? Supposedly distance, in km,
cC                  from the center of earth.
c        CALL BFIELD(ALAT,ALON,6413.,BX,BY,BZ)
c        BFIELD2 = (BX*BX + BY*BY + BZ*BZ)*1.E-10
c      ELSE
c        BFIELD2 = 0.0
c      ENDIF
c
      CBTH2 = 0.5
c
c    MW channels (AMSU and MHS) have same numlev, Plba and Pstda.
c    The new version of opac computes the derivatives dalphadt and 
c    dbetadv, but they are not used here.
cc
c      if(hsbflg) then
c        call opac(bfield2,cbth2,nchan,secant,numlev,ist,
c     &     Pobs, Tair, wvcd, alpha, beta, gamma,dalphadt,dbetadv,
c     &     Pstd_b, Tstd_b, Plb_b, trcoef(1,1,nchan),
c     &     MAXRTACOEF, ntrlev_b, rtatype(nchan))
c      else
c        call opac(bfield2,cbth2,nchan,secant,numlev,ist,
c     &     Pobs, Tair, wvcd, alpha, beta, gamma,dalphadt,dbetadv,
c     &     Pstd_a, Tstd_a, Plb_a, trcoef(1,1,nchan),
c     &     MAXRTACOEF, ntrlev_a, rtatype(nchan))
c      endif
cc
c      if(hsbflg) then
        call opac(bfield2,cbth2,nchan,secant,numlev,ist,
     &     Pobs, tair, wvcd, alpha, beta, gamma,dalphadt,dbetadv,
     &     pstd, tstd, plb, trcoef(1,1,nchan),
     &     MAXRTACOEF, ntrlev, rtatype(nchan))
c      else
c        call opac(bfield2,cbth2,nchan,secant,numlev,ist,
c     &     Pobs, tair, wvcd, alpha, beta, gamma,dalphadt,dbetadv,
c     &     pstd, tstd, plb, trcoef(1,1,nchan),
c     &     MAXRTACOEF, ntrlev, rtatype(nchan))
c      endif
!
!     Correct tau for the thin layer at the surface.
!     Sung-Yung Lee / Sept 10, 1997
!     This correction assumes that the input h2o column densities are 
!     defined as though the profile extends below the surface, hence 
!     the mult by the same ratio throughout.
!
      ratio = (Psurf-Pobs(numlev-1))
     &                  / (pobs(numlev)-pobs(numlev-1))
      alpha(numlev) = alpha(numlev) * ratio
      beta(numlev)  = beta(numlev) * ratio
      gamma(numlev) = gamma(numlev) * ratio

      trnfix = 1.d0
      trnh2o = 1.d0
      trnliq = 1.d0
      do i = 1, numlev
        trnfix = trnfix * DBLE(exp(-alpha(i)))
        if(i.gt.ist) then
          iwv    = i-ist
          trnh2o = trnh2o * DBLE(exp (-wvcd(iwv)*beta(iwv)))
          trnliq = trnliq * DBLE(exp (-lqcd(iwv)*gamma(iwv)))
        end if
        taufix(i) = trnfix
        tauh2o(i) = trnh2o
        tauliq(i) = trnliq
      end do

      dfact = TRCOEF(5,27,nchan)    ! kappa_down/kappa
c
      if(dfact.eq.0.0) dfact = 1.0  ! old files
c
      if(dfact.eq.1.0) then
        do i = 1, numlev
          taufix_r(i) = taufix(i)
          tauh2o_r(i) = tauh2o(i)
          tauliq_r(i) = tauliq(i)
        enddo
      else
        trnfix = 1.d0
        trnh2o = 1.d0
        trnliq = 1.d0
        do i = 1, numlev
          trnfix = trnfix * DBLE(exp(-dfact*alpha(i)))
          if(i.gt.ist) then
            iwv    = i-ist
            trnh2o = trnh2o * DBLE(exp (-wvcd(iwv)*dfact*beta(iwv)))
            trnliq = trnliq * DBLE(exp (-lqcd(iwv)*dfact*gamma(iwv)))
          end if
          taufix_r(i) = trnfix
          tauh2o_r(i) = trnh2o
          tauliq_r(i) = trnliq
        enddo
      endif
 
      return
      end
