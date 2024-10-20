C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              CHKGAS
C
!F77====================================================================


!ROUTINE NAME: CHKGAS


!ABSTRACT:
C    Check the input gas ID to see if it on the list of allowed
C    IDs for use with KLAYERS.  Returns the index of the gas
C    on the list if found, or -1 if not found.  Also returns the
C    molecular mass of the gas.


!CALL PROTOCOL:
C    CHKGAS(ID, IND, MASS)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   ID      gas ID                      none (HITRAN+XSEC)


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   IND     index on gas list           none
C    REAL      MASS    molecular mass              AMU


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): RDINFO


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    none


!COMMON BLOCKS:
C      Allowed gas IDs from "cbgids.f"
C      COMMON /COMGID/ GIDS, GMASS
C      INTEGER GIDS(NGIDS)
C      REAL GMASS(NGIDS)


!DESCRIPTION:
C    Loops over the list of allowed gas IDs looking for the input gas
C    ID.  If it finds it, it returns the index of the gas on the list,
C    otherwise it returns a -1 value to indicate the gas was not found.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C  2 Feb 2001 Scott Hannon      Created
C 27 Feb 2001 Scott Hannon      Added MASS


!END====================================================================


C      =================================================================
       SUBROUTINE CHKGAS(ID, IND, MASS)
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
       INTEGER ID
C
C      Output:
       INTEGER  IND
       REAL   MASS
C
C      Allowed gas IDs
       COMMON /COMGID/ GIDS, GMASS
       INTEGER GIDS(NGIDS)
       REAL GMASS(NGIDS)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER IG

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
C      Loop over list of allowed IDs until input ID is found or reach
C      end of list
       IG=1
       IND=-1
       MASS=0.0
 10    IF (ID .EQ. GIDS(IG)) THEN
          IND=IG
          MASS=GMASS(IG)
       ELSEIF (IG .LT. NGIDS) THEN
          IG=IG  + 1
          GOTO 10
       ENDIF
C
       RETURN
       END
