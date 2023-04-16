
   REAL FUNCTION prod( ARR,       ! c Input
                       NCLMN,     ! c Input
                       STARTCLMN, ! c Input
                       ENDCLMN )  ! c Input

c-----------
c Arguments
c-----------
       REAL ARR(NCLMN),
       INTEGER NCLMN,STARTCLMN,ENDCLMN

c------------------
c Local variables 
c------------------
       REAL MULTI
       INTEGER i

       IF ( STARTCLMN > ENDCLMN ) THEN
         PROD = 1
         RETURN
       ELSE
         MULTI = ONE
         DO i = STARTCLMN, ENDCLMN
           MULTI = MULTI * ARR( i )
         END DO
         PROD = MULTI
       ENDIF

    END FUNCTION prod

