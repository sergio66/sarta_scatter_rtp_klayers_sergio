       subroutine planck(vn,tem,rad)

c Temperature to Planck radiance ... monochromatic
c .... version of 01.05.00

        implicit none
        include 'incFTC.f'
        include 'pyang2cld.param'

	real f1,f2,x,y,z,bnt
	real vn,tem,rad

c define this Planck fcn here	
	bnt(x,y,z)=x/(exp(y/z)-1.0)
c define this Planck fcn here	
	
	f1=c1*vn**3
	f2=c2*vn
	
	rad=bnt(f1,f2,tem)

	return
	end
