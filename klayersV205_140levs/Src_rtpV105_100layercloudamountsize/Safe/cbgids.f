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
       DATA GIDS( 3) /  3 /   ! O3 (ozone)
       DATA GIDS( 4) /  4 /   ! N2O (nitrous oxide)
       DATA GIDS( 5) /  5 /   ! CO (carbon monoxide)
       DATA GIDS( 6) /  6 /   ! CH4 (methane)
       DATA GIDS( 7) /  7 /   ! O2 ("oxygen")
       DATA GIDS( 8) /  8 /   ! NO (nitric oxide)
       DATA GIDS( 9) /  9 /   ! SO2 (sulfur dioxide)
       DATA GIDS(10) / 10 /   ! NO2 
       DATA GIDS(11) / 11 /   ! NH3 (ammonia)
       DATA GIDS(12) / 12 /   ! HNO3 (nitric acid)
       DATA GIDS(13) / 13 /   ! OH
       DATA GIDS(14) / 14 /   ! HF
       DATA GIDS(15) / 15 /   ! HCl (hydrogen chloride)
       DATA GIDS(16) / 16 /   ! HBr (hydrogen bromide)
       DATA GIDS(17) / 17 /   ! HI (hydrogen iodide)
       DATA GIDS(18) / 18 /   ! ClO
       DATA GIDS(19) / 19 /   ! OCS
       DATA GIDS(20) / 20 /   ! H2CO
       DATA GIDS(21) / 21 /   ! HOCl
       DATA GIDS(22) / 22 /   ! N2 ("nitrogen")
       DATA GIDS(23) / 23 /   ! HCN
       DATA GIDS(24) / 24 /   ! CH3Cl (methyl chloride)
       DATA GIDS(25) / 25 /   ! H2O2
       DATA GIDS(26) / 26 /   ! C2H2 (acetylene)
       DATA GIDS(27) / 27 /   ! C2H6 (ethane)
       DATA GIDS(28) / 28 /   ! PH3
       DATA GIDS(29) / 29 /   ! COF2
       DATA GIDS(30) / 30 /   ! SF6
       DATA GIDS(31) / 31 /   ! H2S (hydrogen sulfide)
C
C      Valid cross-section gas IDs
       DATA GIDS(32) / 51 /   ! CCl3F (CFC-11)
       DATA GIDS(33) / 52 /   ! CCl2F2 (CFC-12)
       DATA GIDS(34) / 53 /   ! CClF3 (CFC-13)
       DATA GIDS(35) / 54 /   ! CF4 (FCF-14)
       DATA GIDS(36) / 55 /   ! CHCl2F (CFC-21)
       DATA GIDS(37) / 56 /   ! CHClF2 (CFC-22)
       DATA GIDS(38) / 57 /   ! C2Cl3F3 (CFC-113)
       DATA GIDS(39) / 58 /   ! C2Cl2F4 (CFC-114)
       DATA GIDS(40) / 59 /   ! C2ClF5 (CFC-115)
       DATA GIDS(41) / 60 /   ! CCl4 (carbon tetrachloride)
       DATA GIDS(42) / 61 /   ! ClONO2 {aka ClNO3?}
       DATA GIDS(43) / 62 /   ! N2O5
       DATA GIDS(44) / 63 /   ! HNO4
C
C
C      Approximate molecular mass in AMU
C      Values from sum of atomic weights of elements, pg 7-9 & 7-10
C      of American Institute of Physics Handbook (2nd Edition), 1963
       DATA GMASS( 1) /  18.0153 /   ! H2O (water vapor)
       DATA GMASS( 2) /  44.0100 /   ! CO2 (carbon dioxide)
       DATA GMASS( 3) /  47.9982 /   ! O3 (ozone)
       DATA GMASS( 4) /  44.0128 /   ! N2O (nitrous oxide)
       DATA GMASS( 5) /  28.0106 /   ! CO (carbon monoxide)
       DATA GMASS( 6) /  16.0430 /   ! CH4 (methane)
       DATA GMASS( 7) /  31.9988 /   ! O2 ("oxygen")
       DATA GMASS( 8) /  30.0061 /   ! NO (nitric oxide)
       DATA GMASS( 9) /  64.0628 /   ! SO2 (sulfur dioxide)
       DATA GMASS(10) /  46.0055 /   ! NO2 
       DATA GMASS(11) /  17.0306 /   ! NH3 (ammonia)
       DATA GMASS(12) /  63.0129 /   ! HNO3 (nitric acid)
       DATA GMASS(13) /  17.0074 /   ! OH
       DATA GMASS(14) /  20.0064 /   ! HF
       DATA GMASS(15) /  36.4610 /   ! HCl (hydrogen chloride)
       DATA GMASS(16) /  80.9170 /   ! HBr (hydrogen bromide)
       DATA GMASS(17) / 127.9124 /   ! HI (hydrogen iodide)
       DATA GMASS(18) /  51.4524 /   ! ClO
       DATA GMASS(19) /  60.0746 /   ! OCS
       DATA GMASS(20) /  30.0265 /   ! H2CO
       DATA GMASS(21) /  52.4604 /   ! HOCl
       DATA GMASS(22) /  28.0134 /   ! N2 ("nitrogen")
       DATA GMASS(23) /  27.0258 /   ! HCN
       DATA GMASS(24) /  50.4881 /   ! CH3Cl (methyl chloride)
       DATA GMASS(25) /  34.0147 /   ! H2O2
       DATA GMASS(26) /  26.0382 /   ! C2H2 (acetylene)
       DATA GMASS(27) /  30.0701 /   ! C2H6 (ethane)
       DATA GMASS(28) /  33.9977 /   ! PH3
       DATA GMASS(29) /  66.0074 /   ! COF2
       DATA GMASS(30) / 146.0544 /   ! SF6
       DATA GMASS(31) /  34.0799 /   ! H2S (hydrogen sulfide)
C
       DATA GMASS(32) / 137.3685 /   ! CCl3F (CFC-11)
       DATA GMASS(33) / 120.9139 /   ! CCl2F2 (CFC-12)
       DATA GMASS(34) / 104.4594 /   ! CClF3 (CFC-13)
       DATA GMASS(35) /  88.0048 /   ! CF4 (FCF-14)
       DATA GMASS(36) / 102.9235 /   ! CHCl2F (CFC-21)
       DATA GMASS(37) /  86.4689 /   ! CHClF2 (CFC-22)
       DATA GMASS(38) / 187.3765 /   ! C2Cl3F3 (CFC-113)
       DATA GMASS(39) / 170.9219 /   ! C2Cl2F4 (CFC-114)
       DATA GMASS(40) / 154.4673 /   ! C2ClF5 (CFC-115)
       DATA GMASS(41) / 153.8231 /   ! CCl4 (carbon tetrachloride)
       DATA GMASS(42) /  97.4579 /   ! ClONO2 {aka ClNO3?}
       DATA GMASS(43) / 108.0104 /   ! N2O5
       DATA GMASS(44) /  79.0123 /   ! HNO4
C
       END
