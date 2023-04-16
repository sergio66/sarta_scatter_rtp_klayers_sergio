c this calculates the emission layer contribution
c  call thermal_source(slopd,ssa,TRANSMI_C(i,j),B_TOP,B_BOT,B_UP(i,j),B_DN(i,j))
      SUBROUTINE thermal_source(slopd,ssa,tau,B_TOP,B_BOT,BUP,BDN)

        IMPLICIT NONE
        INCLUDE 'pyang2cld.param'

       REAL slopd,ssa,tau, B_TOP, B_BOT, BUP, BDN
  
       bup = (1.0 - tau) * (B_TOP + B_BOT)/2
       bdn = (1.0 - tau) * (B_TOP + B_BOT)/2

       RETURN
       END

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
