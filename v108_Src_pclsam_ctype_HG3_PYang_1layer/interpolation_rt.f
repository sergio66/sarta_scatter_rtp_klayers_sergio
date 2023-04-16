c==========================================================================================
c The subroutine Interpolation_RT is used to get the Transmissive and Reflective funtions from
c the fitting and interpolation based on the pre-calculated databas
c INPUT: 
c      De: Effective size of cloud particles (Cirrus :10-155micron, Water clouds:2-100 u )
c      Tau: Visible optical thickness       (Cirrus: 0.04-100.0,      Water clouds:0.06-150)
c      Wavelength_in: wavelength in u    (3-20.0micron)
c      Obszenang:     Satellite Observing Angle (0-80 Degree, downward:0)
c      IWindex:       1: Ice clouds, 2: Water clouds, other value are illegal.
c OUTPUT:
c      R: Albedo function  at above condition        
c      T: Transmissive Function (including direct and difuse transmissive function)
c  Version: Oct,2002, by Heli Wei   hlwei@ariel.met.tamu.edu
c
C E. Weisz, June 2008: added istat
C S. Machado, January 2008 Changed Nwl to 251, made everything implicit none

c========================================================================================
    
       SUBROUTINE INTERPOLATION_RT(De,Tau,WAVelength_IN, 
     & Obszenang,IWindex,R,T,istat) 

       IMPLICIT NONE
       INCLUDE 'py_include.param'

       REAL De,Tau,WAVelength_IN,Obszenang,R,T
       integer istat,IWindex       

c      USE COMMON_SPACE
        real Wave_length(NWL)
        real Taucirr(Ntcld), Tauwater(Ntcld)
        REAL FitcoeRc(Nwl,Ntcld,MM*MM), FitcoeTc(Nwl,Ntcld,MM*MM)
        REAL FitcoeRw(Nwl,Ntcld,MM*MM), FitcoeTw(Nwl,Ntcld,MM*MM)
        common /py_fittedcoe/Wave_length, Taucirr, Tauwater,
     &                  FitcoeRc, FitcoeTc, FitcoeRw, FitcoeTw

       real tem(MM*MM)   
       real U, Ref1, Ref2, Rwl1, Rwl2, sumfit
       integer k
       integer Nwl1, Nwl2, Nt1, Nt2

c       print *,'begI', De,Tau,WAVelength_IN, Obszenang, IWindex, R, T

       istat=0
       if ((Wavelength_IN.lt.3.0) .or. (Wavelength_IN.gt.20.0) ) then
         write(*,*) 'Wavelength is out of range!!!!!!'  
	 istat=1       
       end if
 
       if ((obszenang.lt.0.0) .or. (obszenang.gt.80.0) ) then
         write(*,*) 'Satellite observing Zenith angle is out of range 
     & (0-80Degree)' 
        istat=2         
       end if
 
       if ((Tau.lt.0.04) .and. (IWindex.eq.1))  then
          tau=0.04
        end if
      if ((Tau.gt.100.0) .and. (IWindex.eq.1))  then
          tau=100.0
        end if

       if (((De.lt.10) .or. (De.gt.155.0)).and. (IWindex.eq.1))  then
       write(*,*) 'Effective size of cirrus clouds is out of the 
     & range (10-155 micron)'
       istat=3         
       end if


      if ((Tau.lt.0.06) .and. (IWindex.eq.2))  then
          tau=0.06
      end if

       if ((Tau.gt.150.0) .and. (IWindex.eq.2))  then
          tau=150.0
        end if

       if (((De.lt.10) .or. (De.gt.155.0)).and. (IWindex.eq.1))  then
         
       write(*,*) 'Effective size of cirrus clouds is out of the 
     & range (10-155 micron)'   
       istat=4      
       end if

      
       if (((De.lt.2) .or. (De.gt.100.0)).and. (IWindex.eq.2))  then
       write(*,*) 'Effective size of water clouds is out of the 
     & range (2-100 micron)' 
       istat=5        
       end if

       if ((IWindex.NE.1) .and. (IWindex.Ne.2))  then
         print *,'Please check IWindex',IWindex
         write(*,*) 'Only water clouds or Cirrus clouds are permitted' 
	 istat=6        
       end if
       
       if (istat .eq. 0) then

         U=cos(3.141593*obszenang/180.)
	 
	Nwl1=(3000.-10000./wavelength_IN)/10.  
c Because Wavenumeber increase from 500 to 3000 in the step of 10 cm-1

        if (Nwl1.lt.1 ) then
	  Nwl1=1
	end if
	  
        if (Nwl1.gt.Nwl-1 ) then
	  Nwl1=Nwl-1
	end if
	
        Nwl2=Nwl1+1

c==================================================================
cfor Cirrus clouds 
      if (IWindex.eq.1) then

         Nt1=sqrt(Tau/0.04)  !Because Tauvis=KK**2*0.04,kk=1,50 
 
        if (Nt1.lt.1 ) then
	  Nt1=1
	end if
	  
        if (Nt1.gt.Ntcld-1 ) then
	  Nt1=Ntcld-1
	end if
 
 	 Nt2=Nt1+1

c        print *,Nwl1,Nwl2,Tau,Nt1,Nt2,De,U

         do k=1,MM*MM
           tem(k)=fitcoeRc(Nwl1,Nt1,k)
c           print *,k,tem(k)
         end do
         Ref1=sumfit(tem,De,U) !!took out MM as first param

         do k=1,MM*MM
           tem(k)=fitcoeRc(Nwl1,Nt2,k)
c           print *,k,tem(k)
         end do
         Ref2=sumfit(tem,De,U)
         Rwl1=Ref2+(Ref1-Ref2)*(Tau-Taucirr(Nt2))/
     & (Taucirr(Nt1)-Taucirr(Nt2))

c        print *,'ref : ',Taucirr(Nt1),Taucirr(Nt2),Ref1,Ref2,Rwl1

        if (Rwl1.lt.0.0) then 
           Rwl1=0.0
         end if

         do k=1,MM*MM
         tem(k)=fitcoeRc(Nwl2,Nt1,k)
         end do
         Ref1=sumfit(tem,De,U) 
         do k=1,MM*MM
         tem(k)=fitcoeRc(Nwl2,Nt2,k)
         end do
         Ref2=sumfit(tem,De,U)
         Rwl2=Ref2+(Ref1-Ref2)*(Tau-Taucirr(Nt2))/
     & (Taucirr(Nt1)-Taucirr(Nt2))
 
c        print *,Ref1,Ref2,Rwl2

         if (Rwl2.lt.0.0) then 
            Rwl2=0.0
     	 end if

        R=Rwl2+(Rwl1-Rwl2)*(Wavelength_in-wave_length(Nwl2))/
     & (wave_length(Nwl1)-wave_length(Nwl2))


         do k=1,MM*MM
         tem(k)=fitcoeTc(Nwl1,Nt1,k)
         end do
         Ref1=sumfit(tem,De,U) 
         do k=1,MM*MM
         tem(k)=fitcoeTc(Nwl1,Nt2,k)
         end do
         Ref2=sumfit(tem,De,U)
         Rwl1=Ref2+(Ref1-Ref2)*(Tau-Taucirr(Nt2))/
     & (Taucirr(Nt1)-Taucirr(Nt2))
c        print *,Ref1,Ref2,Rwl1

         do k=1,MM*MM
         tem(k)=fitcoeTc(Nwl2,Nt1,k)
         end do
         Ref1=sumfit(tem,De,U) 
         do k=1,MM*MM
         tem(k)=fitcoeTc(Nwl2,Nt2,k)
         end do
         Ref2=sumfit(tem,De,U)
         Rwl2=Ref2+(Ref1-Ref2)*(Tau-Taucirr(Nt2))/
     & (Taucirr(Nt1)-Taucirr(Nt2))
c        print *,Ref1,Ref2,Rwl2
 
        T=Rwl2+(Rwl1-Rwl2)*(Wavelength_in-wave_length(Nwl2))
     & /(wave_length(Nwl1)-wave_length(Nwl2))

c         write(*,*)'cirrus',WAVelength_IN,R,T
       end if

c===============================================================
cEnd of water clouds albedo and transmissivity function     


        if (IWindex.eq.2) then
c==================================================================
cfor water clouds  
        Nt1=sqrt(Tau/0.06)  !Because Tauvis=KK**2*0.06,kk=1,50 
 
        if (Nt1.lt.1 ) then
	  Nt1=1
	end if
	  
        if (Nt1.gt.Ntcld-1 ) then
	  Nt1=Ntcld-1
	end if
 
 	 Nt2=Nt1+1

         do k=1,MM*MM
         tem(k)=fitcoeRw(Nwl1,Nt1,k)
         end do
	 Ref1=sumfit(tem,De,U) 

         do k=1,MM*MM
         tem(k)=fitcoeRw(Nwl1,Nt2,k)
         end do
	 Ref2=sumfit(tem,De,U)
	 
         Rwl1=Ref2+(Ref1-Ref2)*(Tau-Tauwater(Nt2))/
     & (Tauwater(Nt1)-Tauwater(Nt2))

        if (Rwl1.lt.0.0) then 
           Rwl1=0.0
         end if

         do k=1,MM*MM
         tem(k)=fitcoeRw(Nwl2,Nt1,k)
         end do
         Ref1=sumfit(tem,De,U) 
         do k=1,MM*MM
         tem(k)=fitcoeRw(Nwl2,Nt2,k)
         end do
         Ref2=sumfit(tem,De,U)
         Rwl2=Ref2+(Ref1-Ref2)*(Tau-Tauwater(Nt2))/
     & (Tauwater(Nt1)-Tauwater(Nt2))
 
         if (Rwl2.lt.0.0) then 
            Rwl2=0.0
     	 end if

        R=Rwl2+(Rwl1-Rwl2)*(Wavelength_in-wave_length(Nwl2))/
     & (wave_length(Nwl1)-wave_length(Nwl2))

         do k=1,MM*MM
         tem(k)=fitcoeTw(Nwl1,Nt1,k)
         end do
         Ref1=sumfit(tem,De,U) 
         do k=1,MM*MM
         tem(k)=fitcoeTw(Nwl1,Nt2,k)
         end do
         Ref2=sumfit(tem,De,U)
         Rwl1=Ref2+(Ref1-Ref2)*(Tau-Tauwater(Nt2))/
     & (Tauwater(Nt1)-Tauwater(Nt2))

         do k=1,MM*MM
         tem(k)=fitcoeTw(Nwl2,Nt1,k)
         end do
         Ref1=sumfit(tem,De,U) 
         do k=1,MM*MM
         tem(k)=fitcoeTw(Nwl2,Nt2,k)
         end do
         Ref2=sumfit(tem,De,U)
         Rwl2=Ref2+(Ref1-Ref2)*(Tau-Tauwater(Nt2))/
     & (Tauwater(Nt1)-Tauwater(Nt2))
 
        T=Rwl2+(Rwl1-Rwl2)*(Wavelength_in-wave_length(Nwl2))
     & /(wave_length(Nwl1)-wave_length(Nwl2))
c===============================================================
cEnd of water clouds albedo and transmissivity function 
c         write(*,*)'water',WAVelength_IN,R,T
       end if
        
c        write(*,*)De,Tau,WAVelength_IN,Obszenang,IWindex,R,T
c        read(*,*)
	 T=exp(T-6.0)

         if (R.lt.0.0) then
           R=0.0
         end if
         if (R.gt.1.0) then
           R=1.0
         end if

        if (T.lt.0.0) then
           T=0.0
         end if
         if (T.gt.1.0) then
           T=1.0
         end if
	 
	 
	 end if !istat ==0


        IF (istat .NE. 0) THEN 
          print *, 'input deff, tau550, freq, theta, IWindex = '
          print *, De,Tau,WAVelength_IN, Obszenang, IWindex
          stop
          END IF
c        print *, 'endI', De,Tau,WAVelength_IN, Obszenang, IWindex, R, T
c        stop

        RETURN
        END
