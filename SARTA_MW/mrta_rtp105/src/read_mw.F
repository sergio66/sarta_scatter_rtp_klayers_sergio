
      subroutine read_mw(readmode, iounit, trform, TCB, fcenter,
     $           ncoefx, nchl, fx)
      implicit none
      logical*4 trform
      real*4    fx
      real*4    Tcb     ! cosmic background temperature
      real*4    fcenter
      integer*4 readmode
      integer*4 iounit 
      integer*4 ncoefx  ! NOT USED
      integer*4 nchl

#include "amsurta.com"
#include "../src_jpl/paramet.com"

      integer*4 i, n, L, iprt, iret, nlev, ncoef

c     stuff for printout of RTA COEF's
c     --------------------------------
      logical*4 debug
      integer*4 Nfreq
      character*10 rtaname(4)
      character*10 TRLINE(8)
      parameter (debug=.false.)
      data rtaname/'  WINDOW  ','   WATER  ','  OXYGEN  ','old-WATER '/

      ncoef = ncoefx

      if(readmode.eq.2) then
        nlev = ntrlev_a
        ncoef = ncoef_a
      elseif(readmode.eq.3) then
        nlev = ntrlev_b
        ncoef = ncoef_b
      else
        i = 4
        call softexit('READ_MW ',i)
      endif

      iprt = 0

      n = nchl - lastIRchl   ! ORIGINAL AMSU or HSB+AMSU ID number
      call readtram(iounit, n, fx, nlev, ncoef, iret)

      if(iret.ne.0) fx = -1.0  ! freq set negative upon error

      magchl(n) = .false.
      if(TRCOEF(1,26,n).ne.0.0 .and. TRCOEF(7,9,n).ne.0.0) then
        magchl(n) = .TRUE.
        print 4000, n
      endif

      rtatype(n) = 1                          ! assume window or new H2O
      if(TRCOEF(1,26,n).ne.0.0) then
         rtatype(n) = 3  ! O2 channel
      else
         if(TRCOEF(2,26,n).ne.0.0) rtatype(n) = 4  ! old H2o type
      endif
      Nfreq = 0
      if(rtatype(n).eq.1) then
        do i = 1, 6
          if(TRCOEF(3,26+i,n).ne.0.0) Nfreq = Nfreq + 1
        enddo
      endif
      if(Nfreq.gt.0) rtatype(n) = 2   ! H2O, NOT window

      if(rtatype(n).eq.4.or.TRCOEF(5,28,n).eq.0.0) then
         TCB = 2.73
         fcenter = fx
      else
        fcenter = TRCOEF(5,28,n)
        TCB = TRCOEF(5,29,n)
      endif

      if(debug.and.readmode.ge.2) then   ! AMSU's

        print 1010, CHAR(12), rtaname(rtatype(n)), n, fx  ! unknown type
        do L = 1, 66
          if(L.lt.36.or.L.gt.64) then
            do i = 1, 7
              if(TRCOEF(i,L,n).eq.0.0) then
                TRLINE(i) = '      zero'
              else
                write(TRLINE(i),'(e10.3)') TRCOEF(i,L,n)
              endif
            enddo
            print 1020, L, (TRLINE(i), i=1,7)
          elseif(L.eq.55) then
            print 1030
          endif
        enddo
      endif

      return
  100 format(/,3i6,f12.8)
  110 format(i3,7e12.5)
 1010 format(a1,/,'Transmittance coefficients for ',a9,' channel',
     1   ' f(',i3,') =',f9.4)
 1020 format(i4,7(1x,a10))
 1030 format(' ...',7('        ...'))
 4000 format('READ_MW: mag. field correction applied to chl=',i3)
      end

