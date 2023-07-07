C=======================================================================
C=======================================================================
C
C        GODDARD SPACE FLIGHT CENTER
C
C        AIRS
C
C        SUSSKIND RETRIEVALS
C
!F77====================================================================

!ROUTINE NAME:  AIRSRAD1
C
!ABSTRACT:
C
c     AIRSRAD1() computes the clear radiance.  This routine is
c        used to compute the forward problem for cloud cleared
C        or CLEAR component of the cloudy radiance.
C
C
!CALL PROTOCOL:
C
C      SUBROUTINE MW_BT(iphys, numlev, Pobs, tau_up,
C     $    Tlayer, Tsurf, emissiv, tau_r, bigbang, rho_therm,
C     $    raddn, radcof)
C
!INPUT PARAMETERS:
C
C   Name        Type     DESCRIPTION
C   ----        ----     -----------
C   iphys       integer  level of physics support
C   numlev      integer  # of levels
C   Pobs(*)     real     layer boundary pressure levels
C   tau_up(*)   real     layer to space transmittance between surface
C                        and instrument at the computed secant angle
C   Tlayer(*)     real     T(Pavg(L)), layer temperature, K  (NOT LEVEL)
C   Tsurf       real     surface skin temperature
C   emissiv     real     surface emissivity (dimensionless)
c   bigbang     real     cosmic background radiation for this channel
c   tau_r(*)    real     reflective upwelling transmissivity
C   rho_therm   real     surface thermal bi-directional reflectance (per steradian)
C   raddn       real     component of downwelling radiance at spacecraft
C                        for uW: raddn = downwelling integral
!OUTPUT PARAMETERS:
C
C   Name        Type     DESCRIPTION
C   ----        ----     -----------
C   radcof      real     total radiance, mW/M^2/steradian/cm^-1
C
!INPUT/OUTPUT PARAMETERS: none
C
C
!RETURN VALUES: none
C
!PARENT(S):
C
C       any routine computing radiance
C
!ROUTINES CALLED:
C
!FILES ACCESSED:  None
C
!COMMON BLOCKS:
C
C      paramet.com --> src_jpl/paramet.com  to get PLCON1, PLCON2, bigbang
!DESCRIPTION:
C
c     It is also true that Pobs(numlev-1) < Psurf <= Pobs(numlev)
c     which is ensured elsewhere.
c
c     There has been much confusion regarding the meaning of the
c     reflectance terms.  The literature and handbooks do not
c     help, since definitions vary with regard to the use 1/pi
c     with this term.  We have defined the reflectance to have
c     the units of steradian^-1, which would have a value of
c     (1-emissvity)/pi for a Lambertion (diffuse) surface or
c     cloud.
c
C
c     iphys
c        0  same as 1 except getemis1() has variable MW emissivity
c        1  64 level, Psurf=1000, clouds are transmissive
c                     (cldy_rad & airsrad2())
c        2  66 level, variable Psurf, clouds are opaque
c                and surface reflectance divided by PI
c        3  same as 2, except surface reflectance bug fixed
c        4  same as 2 (12/95 UMBC RTA, matches JPL code)
c        5  MW downwelling includes cosmic background radiation
C     
!ALGORITHM REFERENCES:  None
C
!KNOWN BUGS AND LIMITATIONS:  None
C
!ROUTINE HISTORY:
C   Date      Programmer      Description
C    4/98     C.Barnet        iphys=4 upgrades
C    8/98     C.Barnet        iphys=5 upgrades
C    9/99     C.Barnet        imbedded planck function
C    3/00     C.Barnet        separated IR from MW and renamed f/ AIRSRAD1()
C
!END====================================================================
C
C
C=======================================================================


      SUBROUTINE MW_BT(iphys, numlev, Pobs, tau_up,
     $    Tlayer, Tsurf, emissiv,
     $    tau_r, bigbang, rho_therm, raddn, radcof)
      implicit none

c  paramet.com: PLCON1, PLCON2, tauRTthresh
c#include "../src_jpl/paramet.com"
      include 'paramet.com'

c  Input variables.
c  ----------------
      integer*4    iphys
      integer*4    numlev
      real*4       Pobs(*)
      real*4       tau_up(*)  ! upwelling transmittance (f/MW_RTA())
      real*4       Tlayer(*)  ! LAYER temperature
      real*4       Tsurf
      real*4       emissiv
      real*4       tau_r(*)   ! upwelling reflective transmittance (f/MW_RTA())
      real*4       rho_therm  ! thermal reflectance
      real*4       bigbang

c  Output variables
c  ----------------
      real*4       radcof
      real*4       raddn
c
c  Local variables.
c  -----------------
      integer*4    L
      real*4       tauddn, taudup, taubot, taubot_r

c     variables for iphys=2 modification
c     ----------------------------------

      taubot = tau_up(numlev)
      raddn = 0.0    ! BUG fixed 7/31/98  ######

      radcof = ( 1.0 - tau_up(1) ) * Tlayer(1)
      do L = 2, numlev
         radcof  = radcof + ( tau_up(l-1) - tau_up(L) ) * Tlayer(L)
      enddo

      if ( taubot .gt. tauRTthresh ) then

c       ------------------------------------
c       add radiances emitted by the surface
c       ------------------------------------

        radcof = radcof + emissiv * taubot * Tsurf

c       ----------------------------------------------
c       Add downward radiance reflected by the surface
c       ----------------------------------------------

        taubot_r = tau_r(numlev)
        tauddn = taubot_r

        raddn = 0.0
        do L = 1, numlev
          taudup  = tauddn
          tauddn  = taubot_r / tau_r(L)
          raddn   = raddn + (tauddn-taudup)*Tlayer(L)
        enddo

c       add cosmic background
c       ---------------------

        if(iphys.ge.5) then
           raddn = raddn + bigbang * tau_r(numlev)  ! incoming transmittance
        endif

        radcof   = radcof + raddn * rho_therm * taubot  ! outgoing trans.

      endif  ! tau.gt.tauRTthresh

      return
      end

