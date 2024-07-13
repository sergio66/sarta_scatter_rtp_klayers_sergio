c
c Stand-alone Microwave RTA
c
c this reads the RTP file specified by fin
c and writes  an RTP file specified by fout
c
c  P. Rosenkranz, June 27, 2001 for RTP vers. 0.99
c                 Sept.27, 2001 RTP vers. 1.01
c                 Oct. 18, 2001 check raob top P
C                 Oct. 24, 2001 added usestb control
c                 Mar. 20, 2002 multiple H2O units
c                 Mar. 28, 2002 mod for RTP vers. 1.05
c L. LeBlanc      Apr. 18, 2002 added capability to accept H2O in mol/cm**2
c                 Apr. 25, 2002 added comments
c                 May   2, 2002 modified to call C. Barnet's code
c------------------------------------------------------------------
c
c subroutines called:
c  rtpopen, rtpwrite, rtpclose, getcoef, fastem, amsutau, mw_bt
c 
       program m_rta
       implicit none
c####
c RTP declarations
       include '/home/sergio/git/rtp/rtpV221/include/rtpdefs.f'
       integer rtpopen, rtpread, rtpwrite, rtpclose
       record /RTPHEAD/ head
       record /RTPPROF/ prof
       record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
c####
c other variables
       integer i, j, k, n  ! for profile, hinge point, altitude, channel
       integer rtpstat
       integer chan1, chan2
       character*64 mode
c####
c namelist inputs
       character*70 fin     ! input file name
       character*70 fout    ! output file name
       character*80 buf
       character*80 val
       character*80 var

       character*90 tramsua ! transmittance coeff for AMSU-A
       character*90 tramsub ! transmittance coeff for AMSU-B       
       character*90 trhsb   ! transmittance coeff for AMSU-B/HSB
       character*90 tratms  ! transmittance coeff for ATMS

c123456789012345678901234567890123456789012345678901234567890
c0        1         2         3         4         5
c/asl/packages/mrta_rtp105/src/tr_100a.eos
c/home/sergio/SARTA_MW/mrta_rtp201/src_2018/tr_amsua.dat

c set these to remove need for input file
       integer iprt         ! flag = 1 for detailed printout
       parameter (iprt=0)
       integer usestb       ! 0 never use stb, 1 use with landfrac, 2 always use
       parameter (usestb=0)  !no L2 AIRS data
       integer iInstr
       parameter (iInstr=1)  ! AMSU A
c       parameter (iInstr=2)  ! AMSU B
c       parameter (iInstr=3)  ! ATMS
c       parameter (iInstr=4)  ! MHS
cc       parameter (tramsua='/asl/packages/mrta_rtp105/src/tr_100a.eos')
cc       parameter (trhsb='/asl/packages/mrta_rtp105/src/tr_100b.eos')
       parameter (tramsua='/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/SARTA_MW/mrta_rtp201/src_2018/tr_amsua.dat')
       parameter (tramsub='/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/SARTA_MW/mrta_rtp201/src_2018/tr_amsub.dat')
       parameter (tratms='/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/SARTA_MW/mrta_rtp201/src_2018/tr_atms.dat')
       parameter (trhsb='/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/SARTA_MW/mrta_rtp201/src_2018/tr_mhs.dat')                     

c####
c declarations for the MW calculation
c set up the defined parameters first
       integer NMWEM
       parameter (NMWEM=20)  ! number of hinge points
       real h2om,o3m
       parameter (h2om=2.991e-23) ! h2o molec/avogadro number
       parameter (o3m=9.9691e-23) ! o3  molec/avogadro number       
       real grav
       parameter (grav=.980665)  ! gravity
       real degtorad
       parameter (degtorad=57.29578) ! 180/pi to convert from degrees to radians
       integer iphys
       parameter (iphys=5)  ! Barnet's phsyics level parameter
c
       logical levelpro,airspro
       integer nlayer, ior, iounit, iret, levk, kindex
c
c MAXLEV(=120) and MAXCHAN(=20) are set up in rtpdef.f
c
       real p(MAXLEV) ! pressure profile
       real t(MAXLEV) ! temperature profile
       real h2ovcd(MAXLEV) ! water vapour profile in mol/cm^2
       real h2olcd(MAXLEV) ! liquid water profile in mol/cm^2
       real o3cd(MAXLEV)   ! ozone profile in mol/cm^2       
       real mmr, o3mmr, secant, eps,epsv,epsh
       real o3(MAXLEV),h2o(MAXLEV),ti(MAXLEV),tsurf,wind,fracland,fv,emisocean
       real zenang, ratio
       real secratio(MAXCHAN) ! assumed ratio of reflected to direct 
c                                       path length
       real rcalc(MAXCHAN)
       real freqs(NMWEM)  ! frequencies for surface emissivity
       real emismw(NMWEM) ! microwave emissivities at nominal freqs.
       real tcb(NMWEM) ! equiv. cosmic background temps at nominal freqs,
c                   to linearize Planck's function
       integer emch(MAXCHAN) ! index into emismw list from MW channel num.
       real freq(MAXCHAN)      ! channel center freq (GHz)
c       namelist /m_rta_fn/ tramsua,tramsub,tratms,trhsb,fin,fout,iprt,usestb,codev
c####
c Variables needed for Rosencranz's subroutines
       real tran(MAXLEV)
       real tranr(MAXLEV)
       real td !atmospheric emission from direct path
       real tr !atmospheric emission from reflected path
       real e !total atmospheric transmittance on direct path 
       real er !total atmospheric transmittance on reflected path

       integer iNumChan
       
c####
c Variables needed for Barnet's subroutines
       integer*4 numlev
       integer*4 lsurface ! function
c####
c should pull freqs out of here
c AMSU-A see tramsua 23.8 31.4 50.3 52.799 53.595 54.4 54.941 55.499 57.290 57.290 57.290 57.290 57.290 88.997
c AMSU-B see tramsub 89.00 150.00 183.3 183.3 183.3
c ATMS   see tratms  23.8 31.4 50.3 51.76 52.8 53.594 54.4 54.940 55.5 57.290 57.290 57.290 57.290 57.290 88.200 165.500 183.310 183.310 183.310
c MHS    see trahsb 89.0 157.0 183.310 183.310 190.310

       IF (iInstr .EQ. 1) THEN
         iNumChan = 13
         DATA freqs(1:13) /23.8, 31.4, 50.3, 52.8, 53.6, 54.4, 54.9, 55.5, 57.3, 57.31, 57.32,57.33, 89.0/
       ELSEIF (iInstr .EQ. 2) THEN
         iNumChan = 5
         DATA freqs(1:5) /89.00, 150.00, 183.31, 183.32, 183.33/
       ELSEIF (iInstr .EQ. 3) THEN
         iNumChan = 19
         DATA freqs(1:19) /23.8, 31.4, 50.3, 51.76, 52.8, 53.594, 54.4, 54.940, 55.5, 57.290, 57.291, 57.292,
     $	                   57.293, 57.294, 88.200, 165.500, 183.310, 183.310, 183.310/
       ELSEIF (iInstr .EQ. 4) THEN
         iNumChan = 5
         DATA freqs(1:5) /89.0, 157.0, 183.310, 183.311, 190.310/	 
       ENDIF
       
       DATA TCB(1:7)/2.73, 2.73, 2.91, 2.91, 3.27, 4.16, 4.77/
       DATA TCB(8:20) /13*2.73/       
       DATA emch(1:20)/1, 2, 3, 11*4, 2*5, 6, 3*7/
       DATA secratio(1:20) / 3*1.15, 11*1.10, 6*1.15 /
c==========================================================================
c
c zero out liquid water
       do k=1,MAXLEV
         h2olcd(k) = 0.
       end do
c
c get filenames
c       open(iounit,file='m_rta.nl',status='old',form='formatted')
c       read(iounit,m_rta_fn)
c       close(iounit)
c
c read MW transmittance coefficients
c
       iounit = 11
       if (iInstr .EQ. 1) THEN
         open(iounit,file=tramsua,status='old',form='formatted')
         call getcoef(iounit, iret, freq )
c         print *, 'AMSU-A coeff file status =', iret
         close(iounit)

       elseif (iInstr .EQ. 2) THEN
         open(iounit,file=tramsub,status='old',form='formatted')
         call getcoef(iounit, iret, freq )
c         print *, 'AMSU-B coeff file status =', iret
         close(iounit)

       elseif (iInstr .EQ. 3) THEN
         open(iounit,file=tratms,status='old',form='formatted')
         call getcoef(iounit, iret, freq )
c         print *, 'ATMS coeff file status =', iret
         close(iounit)

       elseif (iInstr .EQ. 4) THEN
         open(iounit,file=trhsb,status='old',form='formatted')
         call getcoef(iounit, iret, freq(16) )
c         print *, 'HSB coeff file status =', iret
         close(iounit)
       end if
       
c
c get the filenames from the command line
c
       do i = 1, 2
         call GETARG(i, buf)
c find the "=" character in the command-line argument string
          j=index(buf, '=')
         if (j .ne. 0) then
            var = buf(1:j-1)
            val = buf(j+1:len(buf))
            if (var(1:3) .eq. 'fin') then
               fin=val
            elseif (var(1:4) .eq. 'fout') then
               fout=val
            else
               write(6,1020) var(1:6)
 1020          format('unknown command-line argument: ',a6)
               stop
            endif
         endif
       enddo  ! end of loop over command-line arguments
c
c open the input file for reading
c
       mode = 'r'
       rtpstat = rtpopen(fin, mode, head, hatt, patt, chan1)
c       print *, 'read open status = ', rtpstat
c       print *, 'ptype=', head.ptype
c       print *, 'pfields=', head.pfields
c       print *, 'num. IR chan =',head.nchan
c       print *, 'H2O unit =',head.gunit(1)
       
       levelpro = head.ptype .eq. LEVPRO
c
c levelpro is a check because temps and pressures are converted
c  to layers (for Phil's code). But AIRS has temps and pressures in
c  levels, even though gas profiles are in layers. So if head.ptype
c  is AIRSLAY, levelpro is TRUE for working out temp & press but
c  NOT for determining nlayer!! Don't blame me, I discovered this
c  the hard way, too.
c
       airspro = head.ptype .eq. AIRSLAY
c
c add an attribute
c
c       print *,'maxnattr ',MAXNATTR
       i = 1
c       do while (hatt(i).fname(1:1) .ne. char(0) .and. i .lt. MAXNATTR)
       do while (i .lt. MAXNATTR)
c         print *, i, hatt(i).fname
c         print *, i, hatt(i).aname
c         print *, i, hatt(i).atext
         i = i + 1
       end do
       i = i-1
       hatt(i).fname = 'header'//char(0)
       hatt(i).aname = 'm_rta'//char(0)
       hatt(i).atext = 'vers. 1.05 MW TB calc.'//char(0)
       hatt(i+1).fname = char(0)
c
c update the header data
c
       head.pfields = ior(head.pfields, IRCALCBIT)
C Assign channel frequencies; added by Scott H.
       do i = 1,iNumCHAN
         head.vchan(i)=freq(i)
         if (iprt .eq. 0) print *,'iInstr,iNumChan = ',iInstr,iNumChan,i,freq(i)
       enddo       
c
c create the output file
c
       mode = 'c'
       rtpstat = rtpopen(fout, mode, head, hatt, patt, chan2)
c       print *, 'create open status = ', rtpstat
c
c loop on profiles, until EOF
c
       i = 0
10     continue !head of loop
       rtpstat = rtpread(chan1, prof)
c	  
       if (rtpstat .eq. -1) goto 22
c
c check raob quality
       if(prof.txover.gt.500.) goto 10 ! check top pressure of raob
c       
       if((head.gunit(1).ne.1) .and. (head.gunit(1).ne.10) .and.
     &  (head.gunit(1).ne.20) .and. (head.gunit(1).ne.21)) then
         print *,'water unit must be 1,10,20 or 21 - ',head.gunit(1)
	 goto 22
       endif
c check to make sure angles are not equal to -9999; ditto for landfrac
       if((prof.satzen.eq.-9999).or.(prof.scanang.eq.-9999) ) then
	 print *,'invalid AMSU angles',prof.satzen,prof.scanang
	 goto 22
       endif
c       if((prof.mwbszen.eq.-9999).or.(prof.mwbsang.eq.-9999) ) then
c	 print *,'invalid HSB angles'
c	 goto 22
c       endif
       if(prof.landfrac.eq.-9999) then
	 print *,'invalid landfrac'
	 goto 22
       endif
c
c put profiles in vectors for use by mwtran - note that for
c mwtran, the first layer extends from the top of the atmosphere
c to the first pressure level.
c
c Further note that, with RTP, prof.nlevs = 101 and prof.nlays=100;
c There are 100 AIRSlayers in the profile.
c
c using ncep/ecmwf model data results in altitudes with no data in
c them, because of the way matlab requires its structures
       if(levelpro) then 
         nlayer = prof.nlevs
         if(prof.ptemp(nlayer) .eq. 0) then
	   nlayer = nlayer-1 
	 endif
       else
         nlayer = prof.nlevs - 1
         if(prof.ptemp(nlayer-1) .eq. 0) then
	   nlayer = nlayer-2 
	 endif
       endif
c
c LML 4/16/02
c Modified code to now assume that water profile data are given in
c mol/cm^2 and only modify the input if the data aren't given in
c the units required by mwtran
c       
       if(prof.plevs(2) .lt. prof.plevs(1)) then ! reverse the order
         do k=1,nlayer
           levk = nlayer + 1 - k
           p(k) = prof.plevs(levk)
           ti(k) = prof.ptemp(levk)
           if(ti(k).lt.100. .or. ti(k).gt.400.) goto 10
           h2o(k) = amax1(prof.gamnt(levk,1) , 0.)
           o3(k)  = amax1(prof.gamnt(levk,3) , 0.)	   
         end do
       else
         do k=1,nlayer
           ti(k) = prof.ptemp(k)
           if(ti(k).lt.100. .or. ti(k).gt.400.) goto 10
           h2o(k) = amax1(prof.gamnt(k,1) , 0.)
           o3(k)  = amax1(prof.gamnt(k,3) , 0.)	   
           p(k) = prof.plevs(k)
         end do
       endif
c
c define conversion factor for H2O units to mmr (g/g)
       ratio = 0.
       if(head.gunit(1).eq.1) ratio = 1           !mol/cm^2
       if(head.gunit(1).eq.10) ratio = .6221e-6   !ppmv (ratio H2O to air)
       if(head.gunit(1).eq.20) ratio = 1.e-3      !g/kg
       if(head.gunit(1).eq.21) ratio = 1          !g/g
c
c If head.gunit(1)=1, then everything is already in mol/cm**2
c and we don't need to calculate mmr and all that
c       
       t(1) = ti(1)
       if(levelpro) then
         h2ovcd(1) = 0.
         o3cd(1) = 0.	 
       else 
	 h2ovcd(1) = ratio*h2o(1)
	 o3cd(1)   = ratio*o3(1)	 
       endif
       do k=2,nlayer
         if(levelpro) then
           t(k) = ti(k) 
           mmr = ratio*(h2o(k)+h2o(k-1))/2.
           o3mmr = ratio*(o3(k)+o3(k-1))/2.           	   
         else 
	   if(airspro) then  !i.e., actually levels!
	     t(k) = ti(k) 
             mmr = ratio*h2o(k)
             o3mmr = ratio*o3(k)	     
           else              !really layers
	     t(k) = ti(k) 
             mmr = ratio*h2o(k)
             o3mmr = ratio*o3(k)	     
           endif
         endif
         if(head.gunit(1).eq.1) then
	   h2ovcd(k) = mmr
	   o3cd(k) = o3mmr 	   
	 else
	   h2ovcd(k) = mmr*(p(k)-p(k-1))/((1.+mmr)*h2om*grav)
	   o3cd(k) = o3mmr*(p(k)-p(k-1))/((1.+o3mmr)*o3m*grav)	   
	 endif
       end do
c
       i = i + 1
c	
       if(iprt.gt.0) then
c          print *, 'num.levs & layers', prof.nlevs, nlayer, i
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
       if(prof.nemis.eq.NMWEM .and. 
     &  ((fracland.gt. 0.1 .and. usestb.eq.1) .or. usestb.ge.2) ) then
c use observed surface brightness 
         do j=1,NMWEM
           emismw(j) = amin1(prof.robs1(j)/tsurf, 1.)
         end do
       elseif(tsurf.gt.273.15) then !compute ocean-surface emissivities from model
         if (iprt .gt. 0) print *, i,' computing OCEAN emiss tsurf,wind=',tsurf,wind
         zenang = abs(prof.satzen)
         fv = cos(prof.scanang/degtorad)**2 ! polarization rotation
         do j=1,NMWEM
           if(j.eq.6) then  !hinge points 6 and 7 come from AMSU-B/HSB
             zenang = abs(prof.satzen)
             fv = cos(prof.scanang/degtorad)**2 
           endif
           call fastem(freqs(j),zenang,tsurf,wind,.035,epsv,epsh,1,1,1)
           emisocean = fv*epsv + (1.-fv)*epsh
           emismw(j) = (1.-fracland)*emisocean + fracland
           if (i .EQ. 2000 .and. iprt .gt. 0) write(*,'(A,I3,8(F8.4,1X))') 'emiss cal',j,prof.satzen,
     $          prof.scanang,fv,epsv,epsh,emisocean,fracland,emismw(j)
         end do
c
      else !use default values
        do j=1,NMWEM
            emismw(j) = 1.
        end do
        if (iprt .gt. 0) print *, i,' WARNING: using unity for MW emissivities'
      endif
c
c store emissivities
        do j=1,NMWEM
          prof.emis(j) = emismw(j)
          prof.efreq(j) = freqs(j)
        end do
        prof.nemis = NMWEM
C
c numlev gives you the layer number of the surface
c
        numlev = lsurface(nlayer, p, prof.spres, 100., 1100. )

        secant = 1./cos(prof.satzen/degtorad)
c
        if (i .eq. 2000 .and. iprt .gt. 0) print *,'prof 2000',secant,numlev

        do n = 1, iNumChan
c use the AMSU-B/HSB angles for 16-20
          if(n.eq.16) secant = 1./cos(prof.satzen/degtorad)
          eps = emismw(emch(n))
c
          if(airspro) then
            call mwtran(secant,secratio(n),tran,tranr,numlev,
     &        p,t,airspro,h2ovcd,h2olcd,o3cd,prof.plat,prof.plon,n)
          else
            call mwtran(secant,secratio(n),tran,tranr,numlev,
     &        p,t,levelpro,h2ovcd,h2olcd,o3cd,prof.plat,prof.plon,n)
          endif
c
  	  call tb11(td,tr,e,er,numlev,t,tran,tranr)
c
          rcalc(n) = td + e*( eps*tsurf + (1.-eps)*(tr+er*tcb(emch(n))) )
	  prof.rcalc(n) = rcalc(n)
          if (i .eq. 2000 .and. iprt .gt. 0) then 
            write(*,'(A,I3,1X,8(F8.4,1X))') 'prof 2000 : ',n,tsurf,eps,
     $            td,e,tr,er,tcb(emch(n)),rcalc(n)	  
         end if
        end do
c
c write out the updated profile
c
c        rtpstat = rtpwrite(chan2, prof)
        CALL wrtrtp(i,chan2,iNumChan,rcalc,prof)
c        if(rtpstat.ne.0) print *, 'write status = ', rtpstat, i
c
        goto 10
c
 22     continue
c
        rtpstat = rtpclose(chan1)
        print *, 'read close status = ', rtpstat
c
        rtpstat = rtpclose(chan2)
        print *, 'write close status = ', rtpstat, '  num.recs = ',i
c
        stop
        end
