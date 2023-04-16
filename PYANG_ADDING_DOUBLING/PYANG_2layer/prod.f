
       REAL FUNCTION prod( ARR, NCLMN, STARTCLMN, ENDCLMN )

        IMPLICIT NONE
        INCLUDE 'pyang2cld.param'

c-----------
c Arguments
c-----------
       INTEGER NCLMN,STARTCLMN,ENDCLMN
       REAL ARR(NCLMN)

c------------------
c Local variables 
c------------------
       REAL MULTI
       INTEGER i

       IF ( STARTCLMN .GT. ENDCLMN ) THEN
         PROD = 1
         RETURN
       ELSE
         MULTI = ONE
         DO i = STARTCLMN, ENDCLMN
           MULTI = MULTI * ARR( i )
         END DO
         PROD = MULTI
       ENDIF

      RETURN
      END

