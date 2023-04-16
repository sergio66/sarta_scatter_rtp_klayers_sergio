C **
	function wnplan(vn,tem)
c * Temperature to Planck radiance ... monochromatic
c .... version of 01.05.00

	implicit real*8 (a-h,o-z)
	real tem
	real*8 wnplan
	

c * Boltzmann Constant
	parameter (b = 1.380658d-16)
c * Velocity of Light
	parameter (c = 2.99792458d+10)
c * Planck Constant
	parameter (h = 6.6260755d-27)

	parameter (c1 = 2.d0*h*c*c)
	parameter (c2 = h*c/b)
	bnt(x,y,z)=x/(dexp(y/z)-1.d0)
	


	f1=c1*vn**3
	f2=c2*vn
	t=tem
	
	rad=bnt(f1,f2,t)
	wnplan=rad

	return
	end
