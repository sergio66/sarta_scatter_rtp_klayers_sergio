	function wnbrit(vn,rad)
c * Radiance to brightness temperature ... monochromatic
c .... version of 01.05.00

	implicit real*8 (a-h,o-z)
	real rad,wnbrit

c * Boltzmann Constant
	parameter (b = 1.380658d-16)
c * Velocity of Light
	parameter (c = 2.99792458d+10)
c * Planck Constant
	parameter (h = 6.6260755d-27)

	parameter (c1 = 2.d0*h*c*c)
	parameter (c2 = h*c/b)
	tnb(x,y,z)=y/dlog(x/z+1.d0)

	f1=c1*vn**3
	f2=c2*vn
	r=rad
	tbb=tnb(f1,f2,r)
	wnbrit=tbb

	return
	end
