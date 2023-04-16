       real function wnplan(vn,tem)

c Temperature to Planck radiance ... monochromatic
c .... version of 01.05.00

        include 'py_include.param'
        include 'incFTC.f'

	real f1,f2,t,x,y,z
	real tem

c define this Planck fcn here	
	bnt(x,y,z)=x/(exp(y/z)-1.0)
c define this Planck fcn here	
	
	f1=c1*vn**3    !!have a factor of 1000 between PY and Scott
	f2=c2*vn
	t=tem
	
	rad=bnt(f1,f2,t)
	wnplan=rad

	return
	end
