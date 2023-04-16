       REAL FUNCTION cumsum2prod(raaARR,IANG_INDEX,NLAYER,START,END)

        IMPLICIT NONE
        INCLUDE 'pyang2cld.param'

c-----------
c Arguments
c-----------
       INTEGER NLAYER,START,END, IANG_INDEX
       REAL raaARR( 0:NLEVEL_MAX - 1, NANG_MAX )

c------------------
c Local variables 
c------------------
       REAL MULTI,QIKEXP

       IF ( START .GT. END ) THEN
         MULTI = 0.0
       ELSE
         MULTI = raaARR(END,IANG_INDEX) - raaARR(START-1,IANG_INDEX)
       ENDIF

      CUMSUM2PROD = QIKEXP(-MULTI)

      RETURN
      END

