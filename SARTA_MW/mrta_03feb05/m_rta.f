c
c Stand-alone Microwave RTA
c
c this reads the RTP files specified by fin
c and writes  an RTP file specified by fout
c
c  P. Rosenkranz, June 27, 2001 for RTP vers. 0.99
c                 Sept.27, 2001 RTP vers. 1.01
c                 Oct. 18, 2001 check raob top P
C                 Oct. 24, 2001 added usestb control
c                 Mar. 20, 2002 multiple H2O units
c                 Mar. 28, 2002 mod for RTP vers. 1.05
c                 Aug. 20, 2002 interp missing data
C                 Dec. 12, 2002 code for gunit=1
c                 Jan. 2,  2003 fixed bug in pressure for layer files
c                 June 3,  2003 code for gunit=21
c                 June 18, 2003 flags added to namelist
C                 Sept.10, 2003 optional Lambertian approx.
c                 Mar. 22, 2004 use retrieved secant ratio along with stb
c                 Feb. 2, 2005  replaced fastem with amsuemis
c
c subroutines called:
c  rtpopen, rtpwrite, rtpclose, getcoef, amsuemis, mwtran, tb11
c
c Known bugs and limitations:  roughness effects in amsuemis are derived 
c  from a 705 km orbit.
c 
c Required input fields:
c  prof.txover, prof.plevs, prof.ptemp, prof.gamnt, prof.landfrac
c  prof.mwstb, prof.udef(9) (if usestb > 0)
c
        program m_rta
        implicit none
c
c RTP declarations
c
        include 'rtpdefs.f'
        integer rtpopen, rtpread, rtpwrite, rtpclose
        record /RTPHEAD/ head
        record /RTPPROF/ prof
        record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
c 
c other variables
c
        integer i, j
        integer status
        integer chan1, chan2
        character*64 mode
c
c namelist inputs
c
        character*64 fin      !  input file name
        character*64 fout     !  output file name
        character*64 tramsua  !  transmittance coeff for AMSU-A
        character*64 trhsb    !  transmittance coeff for AMSU-B/HSB
        integer iprt          !  flag = 1 for detailed printout
        integer usestb        !  0 never use stb
                              !  1 use stb if exists and landfrac > 0.1
                              !  2 always use stb (if it exists)
c
c declarations for the MW calculation
c
        integer NMWEM
        parameter (NMWEM=7)
        real h2om, ratio, grav
        parameter (h2om=2.991e-23, grav=.980665)
        logical levelpro,quasiv
        integer nlayer, ior, iounit, iret, k, ii, nlayer1
        real p(MAXLEV), t(MAXLEV), tran(MAXLEV), tranr(MAXLEV)
        real h2ovcd(MAXLEV), h2olcd(MAXLEV), mmr, secant, eps
        real h2o(MAXLEV),ti(MAXLEV),tsurf,wind,fracland,emisocean
        real secrat ! assumed ratio of reflected to direct 
c                     path length over water surface
        real td !atmospheric emission from direct path
        real tr !atmospheric emission from reflected path
        real e !total atmospheric transmittance on direct path 
        real er !total atmospheric transmittance on reflected path 
        real freqs(NMWEM)  ! frequencies for surface emissivity
        real emismw(NMWEM) ! microwave emissivities at nominal freqs.
        real tcb(NMWEM) ! equiv. cosmic background temps at nominal freqs,
c                    to linearize Planck's function
        integer emch(MWMAXCHAN)   ! index into emismw list from MW
c                                 channel num.
        real freq(MWMAXCHAN)      ! channel center freq (GHz)
        namelist /m_rta_fn/ tramsua,trhsb,fin,fout,iprt,usestb
c 
        DATA freqs/23.8, 31.4, 50.3, 52.8, 89., 150., 183.3/
        DATA TCB/2.73, 2.73, 2.91, 2.91, 3.27, 4.16, 4.77/
        data emch/1, 2, 3, 11*4, 2*5, 6, 3*7/
        iounit = 1
        quasiv = .true. ! AMSU polarization
c
c zero out liquid water
        do k=1,MAXLEV
          h2olcd(k) = 0.
        end do
c
c get filenames
        open(iounit,file='m_rta.nl',status='old',form='formatted')
        read(iounit,m_rta_fn)
        close(iounit)
c
c read MW transmittance coefficients
c
        open(iounit,file=tramsua,status='old',form='formatted')
        call getcoef(iounit, iret, freq )
        print *, 'AMSU coeff file status =', iret
        close(iounit)
        open(iounit,file=trhsb,status='old',form='formatted')
        call getcoef(iounit, iret, freq(16) )
        print *, 'HSB coeff file status =', iret
        close(iounit)
c
c open the input file for reading
c
        mode = 'r'
        status = rtpopen(fin, mode, head, hatt, patt, chan1)
        print *, 'read open status = ', status
        print *, 'ptype=', head.ptype
        print *, 'pfields=', head.pfields
        print *, 'num. IR chan =',head.nchan
        print *, 'H2O unit =',head.gunit(1)
        levelpro = head.ptype .eq. LEVPRO
c
c define conversion factor for H2O units to mmr (g/g)
c
        if(head.gunit(1).eq.1) then
          ratio = 0.
        else if(head.gunit(1).eq.10) then
          ratio = .6221e-6
        else if(head.gunit(1).eq.20) then
          ratio = 1.e-3
        else if(head.gunit(1).eq.21) then
          ratio = 1.
        else
          print *, 'unacceptable H2O units'
          goto 22
        endif
c
c add an attribute
c
        i = 1
        do while (hatt(i).fname(1:1) .ne. char(0) .and. i .lt. MAXNATTR)
          print *, i, hatt(i).fname
          print *, i, hatt(i).aname
          print *, i, hatt(i).atext
          i = i + 1
        end do
        hatt(i).fname = 'header'//char(0)
        hatt(i).aname = 'm_rta'//char(0)
        hatt(i).atext = 'vers. 1.05 MW TB calc.'//char(0)
        hatt(i+1).fname = char(0)
c
c update the header data
c
        head.pfields = ior(head.pfields, MWCALCBIT)
        head.mwmemis = NMWEM
c Note: head.mwmemis is obsolete but still necessary, rtp 1.01 bug.
C       Assign channel frequencies; added by Scott H.
        do i = 1,MWMAXCHAN
           head.mwfchan(i)=freq(i)
        enddo
c
c create the output file
c
        mode = 'c'
        status = rtpopen(fout, mode, head, hatt, patt, chan2)
        print *, 'create open status = ', status
c
c loop on profiles, until EOF
c
        i = 0
10      continue !head of loop
          status = rtpread(chan1, prof)
          if (status .eq. -1) goto 22
c
c check raob quality
        if(prof.txover.gt.500.) goto 10 ! check top pressure of raob
c   check bottom pressure
        if(amax1(prof.plevs(1),prof.plevs(prof.nlevs)).lt.800.) goto 10
c
c
c put profiles in vectors for use by mwtran - note that for
c mwtran, the first layer extends from the top of the atmosphere
c to the first pressure level, which differs from the standard RTP definition.
c This won't make much difference as long as the profile extends
c beyond the weighting function.
c
        if(levelpro) then 
            nlayer = prof.nlevs
        else
            nlayer = prof.nlevs - 1
        endif
        if(prof.plevs(2) .lt. prof.plevs(1)) then ! reverse the order
           do j=1,nlayer
              ii = nlayer + 1 - j
              p(j) = prof.plevs(ii)
              ti(j) = prof.ptemp(ii)
              if(ti(j).lt.100. .or. ti(j).gt.400.) ti(j) = 0.
              h2o(j) = amax1(prof.gamnt(ii,1) , 0.)
              if(ratio*h2o(j).ge.1.) h2o(j) = 0.
           end do
        else
           do j=1,nlayer
              ti(j) = prof.ptemp(j)
              if(ti(j).lt.100. .or. ti(j).gt.400.) ti(j) = 0.
              h2o(j) = amax1(prof.gamnt(j,1) , 0.)
              if(ratio*h2o(j).ge. 1.) h2o(j) = 0.
              if(levelpro) then
                p(j) = prof.plevs(j)
              else
                p(j) = prof.plevs(j+1)
              endif
           end do
        endif
c
c  fill in isolated missing values by interpolation
c
        if(ti(1).le.0. .or. ti(nlayer).le.0.) goto 10 !reject
        nlayer1 = nlayer -1
        do j=2,nlayer1
         if(ti(j).le.0.) then
           if(ti(j+1).le.0.) goto 10 !reject 2 adj. missing 
           ti(j) = ( (p(j+1)-p(j))*ti(j-1) + (p(j)-p(j-1))*ti(j+1) )/
     &       (p(j+1)-p(j-1))
         endif
         if(h2o(j).le.0. .and. h2o(j-1).gt.0. .and. h2o(j+1).gt.0.)
     &    h2o(j) = ( (p(j+1)-p(j))*h2o(j-1) + (p(j)-p(j-1))*h2o(j+1))/
     &       (p(j+1)-p(j-1))
        end do
c
c  compute layer values
c
        t(1) = ti(1)
        h2ovcd(1) = 0.
        do k=2,nlayer
         if(levelpro) then
           t(k) = (ti(k)+ti(k-1))/2.
            mmr = ratio*(h2o(k)+h2o(k-1))/2.
         else
            t(k) = ti(k)
            mmr = ratio*h2o(k)
         endif
         if(ratio.le.0.) then
            h2ovcd(k) = h2o(k)
         else
            h2ovcd(k) = mmr*(p(k)-p(k-1))/((1.+mmr)*h2om*grav)
         endif
        end do
c
        i = i + 1
        if(iprt.gt.0) then
          print *, 'num.levs & layers', prof.nlevs, nlayer, i
          do k=1,nlayer
            write(*,1) prof.plevs(k),p(k),prof.ptemp(k),t(k),
     &      prof.gamnt(k,1),h2ovcd(k)
1           format(1x,4f8.2,2e12.3)
          end do
        endif
c
      if(prof.stemp.gt.100. .and. prof.stemp.lt.400.) then
         tsurf = prof.stemp
      else
         tsurf = ti(nlayer)
      endif
c
      wind = amax1( prof.wspeed, 0. )
      fracland = prof.landfrac
c
c compute emissivities at the 7 nominal frequencies
c
      secant = 1./cos(prof.mwaszen/57.296)
      if(prof.mwnstb.eq.NMWEM .and. 
     &  ((fracland.gt. 0.1 .and. usestb.eq.1) .or. usestb.ge.2) ) then
c        use observed surface brightness 
         do j=1,NMWEM
           emismw(j) = amin1(prof.mwstb(j)/tsurf, 1.)
         end do
         secrat = prof.udef(9) ! must contain the retrieved value of secratio
c                               when observed surface brightness is used.
c
      elseif(tsurf.gt.273.2) then 
c      compute ocean-surface emissivities from model
       do j=1,NMWEM
        call amsuemis(freqs(j),secant,prof.mwasang,quasiv,tsurf,wind,
     &   emisocean,secrat)
        emismw(j) = (1.-fracland)*emisocean + fracland
       end do
       if(fracland.gt.0.5) then
           secrat = 0. ! flag for Lambertian calc.
       endif
c
      else !use default values
        do j=1,NMWEM
            emismw(j) = 1.
        end do
        secrat = 0. ! flag for Lambertian calc.
        print *, i,' WARNING: using unity for MW emissivities'
      endif
c
c store emissivities
        do j=1,NMWEM
          prof.mwemis(j) = emismw(j)
          prof.mwefreq(j) = freqs(j)
        end do
        prof.mwnemis = NMWEM
c
c add MW radiance data (profiles are always in layer mode at this point)
c
        do j = 1, MWMAXCHAN
          if(j.eq.16) secant = 1./cos(prof.mwbszen/57.296)
          call mwtran(secant,secrat,tran,tranr,nlayer,
     &       p,t,.false.,h2ovcd,h2olcd,prof.plat,prof.plon,j)
          call tb11(td,tr,e,er,nlayer,t,tran,tranr)
          eps = emismw(emch(j))
          prof.mwcalc(j) = td + e*( eps*tsurf +
     &     (1.-eps)*(tr+er*tcb(emch(j))) )
        end do
c
c
c write out the updated profile
c
        status = rtpwrite(chan2, prof)
        if(status.ne.0) print *, 'write status = ', status, i

        goto 10

 22     continue

        status = rtpclose(chan1)
        print *, 'read close status = ', status

        status = rtpclose(chan2)
        print *, 'write close status = ', status, 
     &   '  num.recs = ',i

        stop
        end

