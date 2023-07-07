  
        integer*4 MAXRTACOEF,MAXRTALEV,MAXRTACHAN
        parameter (MAXRTACOEF=7,MAXRTALEV=100,MAXRTACHAN=50)
        integer*4 MAXUSERLEV
        parameter (MAXUSERLEV=100)
        integer*4 maxchan
        logical*4 magchl(MAXRTACHAN)   ! magnetic channel flag
        integer*4 rtatype(MAXRTACHAN)  ! MIT type  1=WIN,2=WAT,3=O2,4=WIN(old)
c
        integer*4 lastIRchl  ! last IR channel number
        integer*4 lastAMchl  ! last AMSU channel number
        real*4    TRCOEF(MAXRTACOEF,MAXRTALEV,MAXRTACHAN)

        integer*4 ntrlev_a, ncoef_a
        real*4 pstd_a(MAXRTALEV), Tstd_a(MAXRTALEV), Plb_a(MAXRTALEV)
        real*4 wcdstd_a(MAXRTALEV), ocdstd_a(MAXRTALEV)

        integer*4 ntrlev_b, ncoef_b
        real*4 pstd_b(MAXRTALEV), Tstd_b(MAXRTALEV), Plb_b(MAXRTALEV)
        real*4 wcdstd_b(MAXRTALEV), ocdstd_b(MAXRTALEV)

        real*4 alat_amsu, alon_amsu, secant_amsu, secant_mhs
c
c From mwtran.inc
        INTEGER MAXTRLEV, MAXCOEF
c, MAXCHAN
        PARAMETER (MAXTRLEV=100,MAXCOEF=7
c,MAXCHAN=40)
        INTEGER NCHREAD     !number of channels read
        INTEGER NTRLEV      !NUMBER OF LAYERS FOR TRANSMITTANCE
        REAL PSTD(MAXRTALEV) !PRESSURE (MB) AT MIDPOINT OF LAYERS 
        REAL TSTD(MAXRTALEV) !STANDARD ATMOSPHERE TEMPERATURE (K) AT PSTD
        REAL PLB(MAXRTALEV)  !LOWER BOUNDARY PRESSURES (MB) OF LAYERS

      common/amsurta/ NCHREAD,NTRLEV,PSTD,TSTD,PLB,trcoef,
     &      ntrlev_a,ncoef_a,Pstd_a,Plb_a,Tstd_a,wcdstd_a,ocdstd_a,
     &      ntrlev_b,ncoef_b,Pstd_b,Plb_b,Tstd_b,wcdstd_b,ocdstd_b,
     &      maxchan, magchl, rtatype, lastIRchl, lastAMchl,
     &      alat_amsu, alon_amsu, secant_amsu, secant_mhs

c     MAXRTALEV = maximum number of "levels" in MIT RTA
c       "levels" refers to the fact that the MIT coefficients
c       were formatted in a manner to be consistent with the AIRS
c       RTA.  Thus, they are not really levels in trcoef(), but
c       they are levels in Pstd(),Tstd(), etc.
c
c     MAXUSERLEV = maximum number of levels that the user can
c       pass into the amsutau() routine (for local variables)
c
c     Plb_x()   layer boundaries
c     Pstd_x()  mid-layer -- level calculations
c     Tstd_x(L) temperature in K at Pstd(L)
c     wcdstd_x(L) is the water vapor column density in molecules/cm^2
c              from Pstd(L-1) to Pstd(L).
c              For L=1, Pstd(L-1) is assumed to be from top of atmosphere
c              to Pstd(1)
c     trcoef(7,66,50)  7 coef x 66 levels (or 25 temperatures) x channels
c              + other stuff
c  after 11/2001 the MIT RTA'S have 100 possible levels
c     trcoefa(7,100,50)  7 coef x 100 levels (or 25 temperatures) x channels
c              + other stuff

