C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              CHKGUC
C
!F77====================================================================


!ROUTINE NAME: CHKGUC


!ABSTRACT:
C    Check the input gas amount units code number to see if it on the
C    list of allowed code numbers for use with KLAYERS.  Returns the
C    index of the units code number on the list if found, or -1 if
C    not found.


!CALL PROTOCOL:
C    CHKGUC(ICODE, IND)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   ICODE  gas units code number        none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   IND     index on units code list    none


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): OPNRTP


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    none


!COMMON BLOCKS:
C      Allowed gas ICODE from "cbgasu.f"
C      COMMON /COMGUC/ GUCS
C      INTEGER GUCS(MAXGUC)


!DESCRIPTION:
C    Loops over the list of allowed gas amount units code numbers
C    looking for the input ICODE.  If it finds it, it returns the
C    index of the code number on the list, otherwise it returns
C    a -1 value to indicate the gas was not found.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 21 Feb 2001 Scott Hannon      Created based on chkgas


!END====================================================================


C      =================================================================
       SUBROUTINE CHKGUC(ICODE, IND)
C      =================================================================


C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
       include 'incLAY.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input:
       INTEGER ICODE
C
C      Output:
       INTEGER  IND
C
C      Allowed gas units code numbers
       COMMON /COMGUC/ GUCS
       INTEGER GUCS(MAXGUC)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER I

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************
C
C      Loop over list of allowed code numbers until ICODE is found or
C      we reach the end of list
       I=1
       IND=-1
 10    IF (ICODE .EQ. GUCS(I)) THEN
          IND=I
       ELSEIF (I .LT. MAXGUC) THEN
          I=I  + 1
          GOTO 10
       ENDIF
C
       RETURN
       END
