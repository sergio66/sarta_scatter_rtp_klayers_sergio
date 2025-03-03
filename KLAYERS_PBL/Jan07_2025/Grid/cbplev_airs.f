C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              incLAY
C
!F77====================================================================

!ROUTINE NAME: CBPLEV

!ABSTRACT:
C    Common Block for PLEV

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
       CHARACTER*50 NAMGRD
       REAL PLEV(101)
C
C      NAMGRD is the name of the layering grid
C      template    /'1234567890123456789012345678901234567890'/
       DATA NAMGRD /'AIRS 100 layer grid, 0.005 to 1100 mb   '/
C
       DATA  ( PLEV(I), I = 101, 52, -1 )
     $         /    0.0050,    0.0161,    0.0384,    0.0769,    0.1370,
     $              0.2244,    0.3454,    0.5064,    0.7140,    0.9753,
     $              1.2972,    1.6872,    2.1526,    2.7009,    3.3398,
     $              4.0770,    4.9204,    5.8776,    6.9567,    8.1655,
     $              9.5119,   11.0038,   12.6492,   14.4559,   16.4318,
     $             18.5847,   20.9224,   23.4526,   26.1829,   29.1210,
     $             32.2744,   35.6505,   39.2566,   43.1001,   47.1882,
     $             51.5278,   56.1260,   60.9895,   66.1253,   71.5398,
     $             77.2396,   83.2310,   89.5204,   96.1138,  103.0172,
     $            110.2366,  117.7775,  125.6456,  133.8462,  142.3848/
       DATA  ( PLEV(I), I = 51, 1, -1 )
     $         /  151.2664,  160.4959,  170.0784,  180.0183,  190.3203,
     $            200.9887,  212.0277,  223.4415,  235.2338,  247.4085,
     $            259.9691,  272.9191,  286.2617,  300.0000,  314.1369,
     $            328.6753,  343.6176,  358.9665,  374.7241,  390.8926,
     $            407.4738,  424.4698,  441.8819,  459.7118,  477.9607,
     $            496.6298,  515.7200,  535.2322,  555.1669,  575.5248,
     $            596.3062,  617.5112,  639.1398,  661.1920,  683.6673,
     $            706.5654,  729.8857,  753.6275,  777.7897,  802.3714,
     $            827.3713,  852.7880,  878.6201,  904.8659,  931.5236,
     $            958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170,
     $           1100.0000 /
       END
