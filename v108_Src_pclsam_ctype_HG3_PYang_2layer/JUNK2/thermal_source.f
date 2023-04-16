c this calculates the emission layer contribution
c  call thermal_source(slopd,ssa,TRANSMI_C(i,j),B_TOP,B_BOT,B_UP(i,j),B_DN(i,j))
      SUBROUTINE thermal_source(slopd,ssa,tau,B_TOP,B_BOT,B_AVG,BUP,BDN)

        IMPLICIT NONE
        INCLUDE 'pyang2cld.param'

       REAL slopd,ssa,tau, B_TOP, B_BOT, B_AVG, BUP, BDN
  
c try weighted planck radiance avg
       bup = (1.0 - tau) * (B_TOP + B_BOT)/2
       bdn = (1.0 - tau) * (B_TOP + B_BOT)/2

c try avg temp radiance
       bup = (1.0 - tau) * B_AVG
       bdn = (1.0 - tau) * B_AVG

       RETURN
       END

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
