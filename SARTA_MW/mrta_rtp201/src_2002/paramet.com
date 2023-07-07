c     items in this common block are used by
c     routines in
c       ~/umbc9803/
c       ~/src_gsfc/
c       ~/simrad/
c       ~/src/
c
      integer*4 MAXCLD
      integer*4 MAXEMIS
      integer*4 MAXDESCR
      integer*4 MAXLEV
      integer*4 MAXFOV
      integer*4 MAXSPARE

      parameter (MAXCLD   =   8 )
      parameter (MAXEMIS  = 100 )
      parameter (MAXDESCR =  50 )
      parameter (MAXLEV   = 100 )
      parameter (MAXFOV   = MAXCLD + 1 )
      parameter (MAXSPARE = 200 )

c     fundamental constants
c     ---------------------
c     PT99 = American Institute of Physics 1999 publication
c            Cohen, E.R, B.N. Taylor.  "The Fundamental Physical
c            Constants", Physics Today Buyers Guide,  Aug. 1999
c     ESAA = Explanatory Supplement to the Astronomical Almanac, 1992
c     CRC  = CRC Chemistry and Physics Handbook, 1994
c     NIST = http://physics.nist.gov/cuu/Constants/index.html & ref/ contained 3/16/00

      real*4    pi, deg2rad, rad2deg
      parameter (pi = 3.14159265359)
      parameter (deg2rad = pi/180.0)
      parameter (rad2deg = 180.0/pi)

      real*4    Navog      ! NIST Avogadro's number (molecules/mole)
      real*4    BOLTZMNS   ! NIST Boltzmann's constant (erg/K), k
      real*4    PLANCKS    ! NIST Planck's constant (erg*s), h
      real*4    CLIGHT     ! NIST speed of light (cm/s), c
      real*4    Nloschmidt ! NIST Loschmidt's molecules/cm^3 per atm @ STP
c
      parameter (Navog=6.02214199E+23)       
      parameter ( BOLTZMNS = 1.3806503E-16 )
      parameter ( PLANCKS  = 6.62606876E-27 )
      parameter ( CLIGHT   = 2.99792458E+10 )
      parameter (Nloschmidt= 2.6867775E+19 )
c
c     JPL values in parameters.inc
c--   parameter (Navog=6.0221367E+23)
c--   parameter (BOLTZMNS=1.380658E-16)
c--   parameter (Nloschmidt=2.686763E+19)
c--   parameter (CLIGHT= 2.99792458E+10)
c--   parameter (PLANCKS=6.6260755E-27)
c--   parameter (CDair=2.152e25)         !        = Ps*AVOGAD/(GRAV*MWT)

c     chemical constants
c     ------------------
      real*4     P_std
      real*4     T_std
      real*4     MW_d, MW_w, MW_O3, MW_CH4, MW_CO, MW_CO2
      parameter (P_std = 1013.250)        ! PT99, standard pressure, milli-bar
      parameter (T_std =  273.150)        ! standard temperature, K
      parameter (MW_w  = 18.0151)         ! gm/mole  water
      parameter (MW_d  = 28.9644)         ! gm/mole  dry air
      parameter (MW_O3 = 47.9982)         ! CRC, gm/mole  ozone
      parameter (MW_CH4= 16.04303)        ! gm/mole  methane
      parameter (MW_CO = 28.0104)         ! gm/mole  carbon monoxide
      parameter (MW_CO2= 44.00995)        ! gm/mole  carbon dioxide

      real*4     eps
      parameter  (eps = MW_w/MW_d)        ! MW_W/MW_D

c     Earth constants
c     ---------------
      real*4     Rearth  ! <R> = (Req*Req*Rpl)^(1/3) (volume of oblate sphere)
      parameter  (Rearth = 6371.004)  ! ESAA values of Req, Rpl
      real*4     G_std
      parameter  (G_std =  980.665)       ! PT99, standard gravity, cm/s

c     for radiative transfer
c     ----------------------

      real SEA_WATER_FREEZE               ! freezing point for seawater
      parameter (SEA_WATER_FREEZE=271.3)

      real*4      tauRTthresh
      parameter (tauRTthresh=0.00001)  ! changed from 0.0001 11/1/99 (per S.Hannon)
      real*4    ASINLIM
      parameter (ASINLIM=0.9999985)  ! ASIN(0.9999985) = 89.9 degrees


c     derived constants
c     -----------------
      real*4     Runiv
      real*4     Rgas
      real*4     CDair
      real*4     PLCON1, PLCON2
c
      parameter (Runiv=Navog*BOLTZMNS)    ! universal gas constant, erg/mole/K
      parameter (Rgas= Runiv/MW_d)        ! erg/gm/K, R of dry air - R*/MW_D
c
c     column density of air @ P_std, molecules/cm^2
      parameter  (CDair=1000.0*P_std*Navog/MW_d/G_std) !  = 2.14824e+25
c
C
C             1000*P_std*Navog
C     CDair = ----------------    Cd(L) = CDair*(P(L-1)-P(L))/P_std
C               MW_D*G_std

c     for imbedded Planck functions
c     -----------------------------

      parameter ( PLCON1   = 2.0*PLANCKS*CLIGHT*CLIGHT)
      parameter ( PLCON2   = PLANCKS*CLIGHT/BOLTZMNS )

