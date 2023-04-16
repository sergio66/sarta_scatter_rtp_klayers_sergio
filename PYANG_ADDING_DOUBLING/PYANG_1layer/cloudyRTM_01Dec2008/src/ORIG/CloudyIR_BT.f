c 	subroutine CloudyIR_BT(tsurf,temper,ttrans,ncld,nlev,IWindex,
c     &      wvnum,De,Tauvis,obszenang,emis,R,T, TopBTcld,radtop,istat)
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c * Calculate the BRIGHTNESS TEMPERATURE at the top of atmosphere with clouds

c  INPUTS ...
c      tsurf:   surface temperature (K)
c     temper:   temperature profile (K)
c     ttrans:   clear level-to-space transmittance
c       ncld:   level number at which the clouds are located, 2 < ncld < nlev
c       nlev:   level number at the bottom of the atmosphere (1 is the top)
c    IWindex:   1 = ice clouds, 2 = water clouds
c      wvnum:   wavenumber within the range 500-2500 cm^-1
c         De:   effective size of cloud particles in micrometers
c     Tauvis:   optical thickness of clouds in the visible
c  obszenang:   observing zenith angle in degrees
c	 emis:   surface emissivity at wvnum
c

c  OUTPUT ...
c   TopBTcld:   brightness temperature (K) at top of atmosphere 
c          R:   cloud reflectivity
c          T:   cloud transmissivity
c

c  SUBROUTINES and FUNCTIONS used

c     INTERPOLATION_RT: Get the cloud reflectance and transmitance functions from
c               the pre-calculated database  
c     radwnplank:       Calculate the Planck radiance from wavenumber and temperature
c     wnbrit:           Convert radiance to brightness temperature
C	Program history
C	This program is based on the one-layer cloudy radiative transfer model 
C	Developed by Wei et al. 2004 (IEEE Transactions on Geoscience and 
C	Remote Sensing)
C 	2006: Re-write IR cloudy radiative transfer equation by Jun Li 
C	Jun Li (Jun.Li@ssec.wisc.edu)
C	2007: Coding by Hal Woolf
C	Hal Woolf (Hal.Woolf@ssec.wisc.edu)
C 	Contact Jun Li for any questions regarding this science codes 
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 	subroutine CloudyIR_BT(tsurf,temper,ttrans,ncld,nlev,IWindex,
     &      wvnum,De,Tauvis,obszenang,emis,R,T, TopBTcld,radtop,istat)


 	parameter (mlev=101)
	real*8 wvx
	integer istat
	dimension temper(*),ttrans(*)
	dimension plapro(mlev),taustar(mlev),tautstar(mlev)
	real*8 wnplan

c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c	write(*,'(''CLD TSURF, PSURF, NLEV:'',2f10.3,i7)') 
c     &    tsurf,psurf,nlev

c * Get the cloud reflectance and transmittance

	if(Tauvis.lt.0.04) then
c ......... optical thickness of clouds is small enough to ignore
	   R=0.0
	   T=1.0
	else	  
	   wvlen=10000./wvnum   ! wavenumber in cm^-1 --> wavelength in microns
	   call INTERPOLATION_RT(De,Tauvis,wvlen,obszenang,IWindex,R,T,istat) 
	   if (istat .gt. 0) then
	     write (*,*) 'Interpolation_RT status=',istat
	   endif
	endif
	
	if (istat .eq. 0) then
	
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	wvx = wvnum
	 
c * Define surface reflectivity

	rsurf=1.0-emis

c * Transform the temperature profile to a Planck radiance profile at the specified
c	wavenumber ... to avoid doing it over and over!

	do l=1,nlev
	   plapro(l)=wnplan(wvx,temper(l))
	enddo

c * Calculate some auxiliary transmittance quantities, making sure to guard against
c	dividing by zero!

	taus2=ttrans(nlev)**2
	tauc2=ttrans(ncld)**2

	do l=1,nlev
	   temp=ttrans(l)
	   if(temp.gt.0.001) then
	      taustar(l)=taus2/temp
	      tautstar(l)=tauc2/temp
	   else
	      taustar(l)=taustar(l-1)
	      tautstar(l)=tautstar(l-1)
	   endif
	enddo


c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c * Ready to calculate the many contributions to the radiance ...

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c + surface  

	radsurf = wnplan(wvx,tsurf)*emis*ttrans(nlev)*T
     
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c + below the cloud

	radba=0.
	radbb=0.
	b1=plapro(ncld)
	t1=ttrans(ncld)
	ts1=taustar(ncld)
	do l=ncld+1,nlev
	   b2=plapro(l)
	   bbar=0.5*(b1+b2)
	   t2=ttrans(l)
	   ts2=taustar(l)
	   radba=radba+bbar*(t1-t2)
	   radbb=radbb+bbar*(ts2-ts1)   
	   b1=b2
	   t1=t2
	   ts1=ts2
	enddo

	radbelow = radba*T + radbb*T*rsurf
      
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c + the cloud

	if(T.ne.1.0) then
	   btc=plapro(ncld)*ttrans(ncld)
	  if(tauc2.ge.0.0001) then
	   radcloud=(1.0-T-R)*(1.0+rsurf*T*(taus2/tauc2))*btc
	  else
	   radcloud=(1.0-T-R)*btc
	  endif
	else
	   radcloud=0.0
	endif

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c + above the cloud

	radba=0.
	radbb=0.
	radbc=0.
	b1=plapro(1)
	t1=ttrans(1)
	ts1=taustar(1)
	tts1=tautstar(1)
	do l=2,ncld
	   b2=plapro(l)
	   bbar=0.5*(b1+b2)
	   t2=ttrans(l)
	   ts2=taustar(l)
	   tts2=tautstar(l)
	   radba=radba+bbar*(t1-t2)
	   radbb=radbb+bbar*(ts2-ts1)    
	   radbc=radbb+bbar*(tts2-tts1)  
	   b1=b2
	   t1=t2
	   ts1=ts2
	   tts1=tts2
	enddo

	radabove = radba + radbb*rsurf*T*T + radbc*R
        
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c + total top-of-atmosphere cloudy radiance

	radtop = radsurf + radbelow + radcloud + radabove

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c + top-of-atmosphere cloudy brightness temperature

	wvx = wvnum
	TopBTcld=wnbrit(wvx,radtop)
        end if !if istat==0
	return
	end
