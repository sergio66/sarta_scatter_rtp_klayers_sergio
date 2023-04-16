c============================================================================
c    To sum the 2-dimension polynomial fitting according the fitted coefficients
c INPUT:
c       MM:  The order of polynomial
c     Tem1:  The fitted coefficients (MM*MM)
c       De:  The effective size of clouds
c        U:  Cos(theta), where theta is the observing zenith angle.
c
c OUTPUT:
c       Sumfit:  The Albedo or transmissive function fitted according the De and U
c
c Sergio Machado Jan 2009 : Took out "MM" from parameter list 
c===========================================================================
	real Function sumfit(tem1,De,U)

        implicit none
        include 'py_include.param'

        real De,U,tem1(MM*MM)

        real tcoe(MM,MM),A(MM)
        integer k,j

	do  k=1,MM 
	 do j=1,MM 
	   tcoe(k,j)=tem1((k-1)*MM+j)
	  end do
	 end do

	 sumfit=0.0
	 do j=1,MM
           a(j)=0.0
	   do k=1,MM
	     a(j)=a(j)+tcoe(j,k)*U**(k-1)
	   end do
	  sumfit=sumfit+a(j)*De**(j-1)
          end do
          return
          end
