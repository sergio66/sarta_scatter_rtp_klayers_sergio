C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              CBGIDS
C
!F77====================================================================

!ROUTINE NAME: CBGIDS

!ABSTRACT:
C    Common Block for allowed gas IDs, and their molecular mass.

!CALL PROTOCOL: none

!INPUT PARAMETERS: none

!OUTPUT PARAMETERS: none

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S):
C    KLAYERS

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C    Assign values to COMGID

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS: none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C  1 Feb 2001 Scott Hannon      Created
C 27 Feb 2001 Scott Hannon      Added GMASS
C  3 Jan 2003 Scott Hannon      Correct GMASS of SO2 to ~64 (was ~32)
C 13 Jul 2010 Scott Hannon      Add gases 32-42 and 64-80
C 16 Jul 2010 Scott Hannon      Add xsec 81


!END====================================================================

C      =================================================================
       BLOCK DATA CBGIDS
C      =================================================================
C
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
C      none

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      EXECUTABLE CODE
C-----------------------------------------------------------------------
C      none; common block data
C
C      Edit as needed, but GIDS must match what is the AFGL model
C      atmospheres data file.
C***********************************************************************
C***********************************************************************
       COMMON /COMGID/ GIDS, GMASS
C
       INTEGER GIDS(NGIDS)  ! gas ID numbers
       REAL   GMASS(NGIDS)  ! approx molecular mass in AMU
C
C      Valid HITRAN gas IDs
       DATA GIDS( 1) /  1 /   ! H2O (water vapor)
       DATA GIDS( 2) /  2 /   ! CO2 (carbon dioxide)
       DATA GIDS( 3) /  3 /   ! N2 ("nitrogen")
       DATA GIDS( 4) /  4 /   ! O2 ("oxygen")
       DATA GIDS( 5) /  5 /   ! CO (carbon monoxide)
       DATA GIDS( 6) /  6 /   ! NO (nitric oxide)
C
C      Approximate molecular mass in AMU
C      Values from sum of atomic weights of elements, pg 7-9 & 7-10
C      of American Institute of Physics Handbook (2nd Edition), 1963
       DATA GMASS( 1) /  18.0153 /   ! H2O (water vapor)
       DATA GMASS( 2) /  44.0100 /   ! CO2 (carbon dioxide)
       DATA GMASS( 3) /  28.0134 /   ! N2 ("nitrogen")
       DATA GMASS( 4) /  31.9988 /   ! O2 ("oxygen")
       DATA GMASS( 5) /  28.0106 /   ! CO (carbon monoxide)
       DATA GMASS( 6) /  30.0061 /   ! NO (nitric oxide)
C
       END
