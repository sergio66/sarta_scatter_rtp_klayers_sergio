       subroutine planck(tem, C1V3, C2V, rad)

c Temperature to Planck radiance ... monochromatic
c .... version of 01.05.00

        implicit none
        include 'incFTC.f'
        include 'pyang2cld.param'

c	real f1,f2,x,y,z,bnt
	real vn,tem,rad, C1V3,C2V

c define this Planck fcn here	
c	bnt(x,y,z)=x/(exp(y/z)-1.0)
c define this Planck fcn here	
	
c	f1=c1*vn**3
c	f2=c2*vn
c	rad=bnt(f1,f2,tem)
c	rad=bnt(C1V3,C2V,tem)

        rad = C1V3/(exp(C2V/tem)-1.0)

	return
	end
