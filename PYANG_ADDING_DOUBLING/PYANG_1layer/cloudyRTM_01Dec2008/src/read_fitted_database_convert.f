c to compile and run, do f77 -f -N3 read_fitted_database_convert.f; a.out
c USE THIS as absoft did not like the unformatted statement when opening unit-23
c    /usr/bin/g77  -fcase-lower -ffixed-line-length-120 read_fitted_database_convert.f; a.out
c this takes in ping yang's database, gets rid of those comments which I am having trouble
c with in Matlab, and just dumps back to another binary file
c notice we do not add record info to unformatted files ie NOT using -N3 option

      call read_fitted_database
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c============================================================================
cThe code belowings are written by Heli Wei, Texaes A&M University.
c Nov. of 2002
c============================================================================

c....................................................................
c    This subroutine read in the pre-calculated fitted coefficeients for 
c    the transmissive function and albedo function of cirrus clouds 
c    and water clouds at various wavelengths, optical thickness.
c Wavelength: (200) 4.0-20.micron, step 10cm-1  
c Cirrus: Optical thickness (50):0.04-100.0  (0.04*K*K)
c Waterclouds: Optical thickness (50):0.06-150.0  (0.06*K*K)
c     Fitting range: Effective size of ice crystals from 10-154micron
c                    Effective size of waterdroplets from 2-100micron (diameter)           
c                    Observing zenith angle (0-80 Degree)   
c The data were produced using DISORT  (rtcoe.f90)
c.....................................................................
  

       SUBROUTINE READ_fitted_database       

c       USE COMMON_SPACE
        parameter (Nwl=201,Ntcld=50, MM=4)
        real Wave_length(NWL)
        real Taucirr(Ntcld), Tauwater(Ntcld)
        real FitcoeRc(Nwl,Ntcld,MM*MM), FitcoeTc(Nwl,Ntcld,MM*MM)
        real FitcoeRw(Nwl,Ntcld,MM*MM), FitcoeTw(Nwl,Ntcld,MM*MM)
        common /fittedcoe/Wave_length, Taucirr, Tauwater,
     &                  FitcoeRc, FitcoeTc, FitcoeRw, FitcoeTw

        character*50 cc  
  
c looks like PingYang stores wavelength
c                            tau, (fitcoeff(k),k=1,MM*MM)
c For the Reflective function of cirrus clouds
 50   FORMAT(A50)
      OPEN(3,FILE='../cloudcoeffs/fitedcoerc.dat')
      OPEN(UNIT=23,FILE='../cloudcoeffs/fitedcoerc_simple.dat',FORM='UNFORMATTED')
      read(3,*)cc
       Do i=1,NWL
	 read(3,*)Wave_length(i),cc
         write(23) Wave_length(i)
     	      do j=1,Ntcld
		      read(3,*) Taucirr(j),(fitcoeRc(i,j,k),k=1,MM*MM)
		      write(23) Taucirr(j)
                      write(23) (fitcoeRc(i,j,k),k=1,MM*MM)
                      print *,i,j,Wave_length(i),Taucirr(j),MM*MM
	      end do         
       
c       print *,'WL=',Wave_length(i),Taucirr
       End do
        close(3)
        close(23)

c For the albedo function of cirrus clouds
      OPEN(3,FILE='../cloudcoeffs/fitedcoetc.dat')
      read(3,*)cc
       Do i=1,NWL
	 read(3,*)Wave_length(i),cc
	      do j=1,Ntcld
		      read(3,*)Taucirr(j),(fitcoeTc(i,j,k),k=1,MM*MM)
	      end do
       End do
      close(3)
      
c For the Reflective function of water clouds
      OPEN(3,FILE='../cloudcoeffs/fitedcoerw.dat')
      read(3,*)cc
       Do i=1,NWL
	 read(3,*)Wave_length(i),cc
     	      do j=1,Ntcld
		      read(3,*)Tauwater(j),(fitcoeRw(i,j,k),k=1,MM*MM)
	      end do         
       
c       print *,'WL=',Wave_length(i),Taucirr
       End do
        close(3)

c For the Reflective function of water clouds
      OPEN(3,FILE='../cloudcoeffs/fitedcoetw.dat')
      read(3,*)cc
       Do i=1,NWL
	 read(3,*)Wave_length(i),cc
	      do j=1,Ntcld
		      read(3,*)Tauwater(j),(fitcoeTw(i,j,k),k=1,MM*MM)
	      end do
       End do
      close(3)

      RETURN
      END
