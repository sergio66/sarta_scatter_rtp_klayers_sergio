c this calculates the emission layer contribution
c  call thermal_source(slopd,ssa,TRANSMI_C(i,j),B_TOP,B_BOT,B_UP(i,j),B_DN(i,j))
      SUBROUTINE thermal_source(tau,B_AVG,BUP,BDN)

        IMPLICIT NONE
        INCLUDE 'pyang2cld.param'

       REAL tau, B_AVG, BUP, BDN
  
       bup = (1.0 - tau) * B_AVG
       bdn = bup

       RETURN
       END

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
