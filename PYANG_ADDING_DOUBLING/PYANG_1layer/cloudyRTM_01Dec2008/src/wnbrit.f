	real*8 function wnbrit(vn,rad)
c * Radiance to brightness temperature ... monochromatic
c .... version of 01.05.00

        implicit none
        include 'py_include.param'

	real*8 tnb,x,y,z,r
	real rad

c define this Planck fcn here	
	tnb(x,y,z)=y/dlog(x/z+1.d0)
c define this Planck fcn here	

	f1=c1*vn**3
	f2=c2*vn
	r=rad
	tbb=tnb(f1,f2,r)
	wnbrit=tbb

	return
	end
