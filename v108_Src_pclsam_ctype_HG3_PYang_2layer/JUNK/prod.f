
       REAL FUNCTION prod(raaARR,IANG_INDEX,NCLMN,STARTCLMN,ENDCLMN)

        IMPLICIT NONE
        INCLUDE 'pyang2cld.param'

c-----------
c Arguments
c-----------
       INTEGER NCLMN,STARTCLMN,ENDCLMN, IANG_INDEX
       REAL raaARR( NLEVEL_MAX - 1, NANG_MAX )

c------------------
c Local variables 
c------------------
       REAL MULTI
       INTEGER i

       IF ( STARTCLMN .GT. ENDCLMN ) THEN
         MULTI = 1.0
       ELSE
         MULTI = 1.0
         DO i = STARTCLMN, ENDCLMN
           MULTI = MULTI * raaARR( i, IANG_INDEX )
         END DO
       ENDIF

      PROD = MULTI

      RETURN
      END

