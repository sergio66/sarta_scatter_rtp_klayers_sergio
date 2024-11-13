C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              CBPLEV_aeri999
C
!F77====================================================================

!ROUTINE NAME: CBPLEV

!ABSTRACT:
C    Common Block for PLEV
C    This version for AERI layers starting at 999 mb

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
C    Assign values to PLEV

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS: none

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 15 Nov 2000 Scott Hannon      PLEV moved from klayers.f to here
C 19 Jul 2001 Scott Hannon      Create version for AERI 999 mb
C 20 Nov 2001 Scott Hannon      Add NAMGRD

!END====================================================================

C      =================================================================
       BLOCK DATA CBPLEV
C      =================================================================
C
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE

C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
C      none

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
       INTEGER I  ! generic looping variable

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      EXECUTABLE CODE
C-----------------------------------------------------------------------
C      none; common block data
C
C      Edit as needed, but PLEV must match NBPLEV in incLAY.f
C      PLEV must be in mb or hPa (these are equivalent), and sorted
C      so that level 1 is the max pressure and NLEV is the min.
C***********************************************************************
C***********************************************************************
       COMMON /COMLEV/ NAMGRD, PLEV
C
       CHARACTER*40 NAMGRD
       REAL PLEV(101)
C
C      NAMGRD is the name of the layering grid
C      template    /'1234567890123456789012345678901234567890'/
       DATA NAMGRD /'AERI999 100 layer grid, 0.005 - 999 mb  '/
C
       DATA  ( PLEV(I), I = 101, 52, -1 )
     $    /   0.0050,    0.0122,    0.0274,    0.0564,    0.1078,
     $        0.1932,    0.3261,    0.5220,    0.7965,    1.1642,
     $        1.6371,    2.2237,    2.9278,    3.7485,    4.6800,
     $        5.7126,    6.8330,    8.1306,    9.6259,   11.3403,
     $       13.2965,   15.5182,   18.0298,   20.8567,   24.0248,
     $       27.5604,   31.4904,   35.8416,   40.6406,   45.9140,
     $       51.6879,   57.9876,   64.8379,   72.2623,   80.2833,
     $       88.9218,   98.1974,  108.1280,  118.7296,  130.0163,
     $      142.0004,  154.6917,  168.0981,  182.2253,  197.0764,
     $      212.6527,  228.9527,  245.9731,  263.7080,  282.1496/
       DATA  ( PLEV(I), I = 51, 1, -1 )
     $    / 301.2878,  321.1103,  341.6032,  362.7504,  384.5341,
     $      406.9350,  429.9320,  453.5027,  477.6235,  502.2693,
     $      527.4142,  553.0314,  579.0933,  605.5714,  632.4371,
     $      659.6609,  687.2136,  715.0652,  743.1862,  769.0486,
     $      792.7365,  814.3535,  834.0170,  851.8521,  867.9877,
     $      882.5530,  895.6743,  907.4739,  918.0683,  927.5671,
     $      936.0733,  943.6821,  950.4816,  956.5526,  961.9690,
     $      966.7982,  971.1011,  974.9331,  978.3441,  981.3792,
     $      984.0786,  986.4787,  988.6122,  990.5081,  992.1924,
     $      993.6886,  995.0173,  996.1972,  997.2447,  998.1746,
     $      999.0000/
       END