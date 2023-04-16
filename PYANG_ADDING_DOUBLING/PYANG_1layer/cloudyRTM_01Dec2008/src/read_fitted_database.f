c to compile and run, do f77 -f -N3 read_fitted_database.f; a.out
c USE THIS as absoft did not like the unformatted statement when opening unit-23
c    /usr/bin/g77  -fcase-lower -ffixed-line-length-120 read_fitted_database.f; a.out
c this takes in ping yang's database, gets rid of those comments which I am having trouble
c with in Matlab, and just dumps back to another binary file
c notice we do not add record info to unformatted files ie NOT using -N3 option

      call read_fitted_database
      end
c============================================================================
cThe code belowings are written by Heli Wei, Texaes A&M University.
c Nov. of 2002
c Modified by Sergio DeSouza-Machado so that dummy junk extends the database from
c   [4um = 2500 cm-1  --> 20 um = 500 cm-1] (2500-500)/10 = 200+1 = 201 pts  to 
c   [4um = 2500 cm-1  --> 20 um = 500 cm-1] (3000-500)/10 = 250+1 = 251 pts
c Hence Nwl has been changed from 201 to 251
c and we fill (1:50) of FitCoefXY with numbers from FitCoefXY(51,j,k) ie the 4 um coeffs
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

       IMPLICIT NONE
       INCLUDE 'py_include.param'

c       USE COMMON_SPACE
        real Wave_length(NWL)
        real Taucirr(Ntcld), Tauwater(Ntcld)
        real FitcoeRc(Nwl,Ntcld,MM*MM), FitcoeTc(Nwl,Ntcld,MM*MM)
        real FitcoeRw(Nwl,Ntcld,MM*MM), FitcoeTw(Nwl,Ntcld,MM*MM)
        common /fittedcoe/Wave_length, Taucirr, Tauwater,
     &                  FitcoeRc, FitcoeTc, FitcoeRw, FitcoeTw
        character*50 cc  
        integer iOffSet,i,j,k,iDebug
        real wavenumber

      !recall Ping Yang only supplied 4 to 20 um, so put zeros for 3 to 4 um
      !this should be the wavenumber spacing
      DO i=1,NWL
        wavenumber = 3000.0 -(i-1)*10.0
        Wave_length(i) = 10000/wavenumber
        IF (wavenumber .gt. 2501.0) iOffSet = i
        END DO

      iDebug = -1
      IF (iDebug .GT. 0) THEN
        print *,'3 to <4 um  = 1 to ',iOffSet, ' = ',(iOffset-1)+1,' points'
        print *,'=4 to 20 um = ', iOffSet+1,' to 251 = ',(251-(iOffSet+1)+1),' points'
        END IF

      iOffSet = iOffSet + 1
      IF (iDebug .GT. 0) print *,'PYs database starts at ',iOffSet

c For the Reflective function of cirrus clouds
      OPEN(3,FILE='../cloudcoeffs/fitedcoerc.dat')
      read(3,*)cc
      Do i=iOffSet,NWL
	read(3,*)Wave_length(i),cc
     	do j=1,Ntcld
	  read(3,*)Taucirr(j),(fitcoeRc(i,j,k),k=1,MM*MM)
	  end do         
        End do
      close(3)

c For the albedo function of cirrus clouds
      OPEN(3,FILE='../cloudcoeffs/fitedcoetc.dat')
      read(3,*)cc
       Do i=iOffSet,NWL
	 read(3,*)Wave_length(i),cc
         do j=1,Ntcld
	   read(3,*)Taucirr(j),(fitcoeTc(i,j,k),k=1,MM*MM)
	   end do
         End do
      close(3)
      
c For the Reflective function of water clouds
      OPEN(3,FILE='../cloudcoeffs/fitedcoerw.dat')
      read(3,*)cc
       Do i=iOffSet,NWL
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
      Do i=iOffSet,NWL
        read(3,*)Wave_length(i),cc
	do j=1,Ntcld
	  read(3,*)Tauwater(j),(fitcoeTw(i,j,k),k=1,MM*MM)
	  end do
         End do
      close(3)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c set 3 --> 4 um as zeros
      Do i=1,iOffSet-1
        do j=1,Ntcld
          do k=1,MM*MM
	    fitcoeRc(i,j,k) = 0.0
	    fitcoeTc(i,j,k) = 1.0
	    fitcoeRw(i,j,k) = 0.0
	    fitcoeTw(i,j,k) = 1.0
            end do
	  end do         
        END do
c set 3 --> 4 um as 4 um results
      Do i=1,iOffSet-1
        do j=1,Ntcld
          do k=1,MM*MM
	    fitcoeRc(i,j,k) = fitcoeRc(iOffSet,j,k)
	    fitcoeTc(i,j,k) = fitcoeTc(iOffSet,j,k)
	    fitcoeRw(i,j,k) = fitcoeRw(iOffSet,j,k)
	    fitcoeTw(i,j,k) = fitcoeTw(iOffSet,j,k)
            end do
	  end do         
        END do
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (iDebug .GT. 0) THEN
        Do i=1,NWL
          Do j=1,Ntcld
            print *,i,j,Wave_length(i),10000/Wave_length(i),
     $        Taucirr(j),fitcoeTc(i,j,1),fitcoeRc(i,j,1)
            END DO
          END DO

        Do i=1,NWL
          Do j=1,Ntcld
            print *,i,j,Wave_length(i),10000/Wave_length(i),
     $      Tauwater(j),fitcoeTw(i,j,1),fitcoeRw(i,j,1)
            END DO
          END DO
        END IF

      RETURN
      END
