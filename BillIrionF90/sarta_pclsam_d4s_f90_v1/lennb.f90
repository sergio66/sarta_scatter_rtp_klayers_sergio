!*****************************************************************************  
      FUNCTION LENNB(STRING)
!*****************************************************************************  
      IMPLICIT REAL (A-H,O-Z)
      CHARACTER*(*) STRING
!
!  FIND LENGTH OF NON-BLANK LEFT SUB-STRING OF STRING:
!
      LENNB = LEN(STRING)
      IF(LENNB.LE.0) THEN
        RETURN
      ELSE IF (LENNB.GT.511) THEN
        LENNB = 511
      END IF
      NL = LENNB
      DO 10 I = NL,1,-1
      IF(ICHAR(STRING(I:I)).EQ.32) GOTO 10
        LENNB = I
        RETURN
   10 CONTINUE
      LENNB = 0
      RETURN
      END
