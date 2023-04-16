       real*8 function wnplan(vn,tem)

c Temperature to Planck radiance ... monochromatic
c .... version of 01.05.00

        include 'py_include.param'

	real*8 f1,f2,t,x,y,z
	real tem

c define this Planck fcn here	
	bnt(x,y,z)=x/(dexp(y/z)-1.d0)
c define this Planck fcn here	
	
	f1=c1*vn**3
	f2=c2*vn
	t=tem
	
	rad=bnt(f1,f2,t)
	wnplan=rad

	return
	end
