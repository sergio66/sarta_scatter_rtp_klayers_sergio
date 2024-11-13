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
       CHARACTER*40 NAMGRD
       REAL PLEV(101)
C
C      NAMGRD is the name of the layering grid
C      template    /'1234567890123456789012345678901234567890'/
       DATA NAMGRD /'AIRS 100 layer grid, 0.00025 to 1100 mb '/
C
C include xcbplev_airs120km.f from ../Src_rtpV105/plevs_0_120km.m
       DATA  ( PLEV(I), I = 101, 52, -1 ) 
     $         /   0.000022,  0.000062,  0.000143,  0.000293,  0.000549,
     $             0.000964,  0.001607,  0.002566,  0.003954,  0.005915,
     $             0.008625,  0.012297,  0.017192,  0.023622,  0.031957,
     $             0.042633,  0.056161,  0.073136,  0.094245,  0.120280,
     $             0.152148,  0.190879,  0.237645,  0.293771,  0.360745,
     $             0.440241,  0.534126,  0.644485,  0.773633,  0.924136,
     $             1.098831,  1.300846,  1.533619,  1.800927,  2.106904,
     $             2.456068,  2.853345,  3.304101,  3.814165,  4.389859,
     $             5.038034,  5.766093,  6.582033,  7.494472,  8.512690,
     $             9.646663, 10.907103, 12.305495, 13.854142, 15.566205/ 
       DATA  ( PLEV(I), I = 51, 1, -1 ) 
     $         /    17.4557,   19.5378,   21.8283,   24.3444,   27.1042,
     $              30.1269,   33.4331,   37.0446,   40.9843,   45.2766,
     $              49.9475,   55.0243,   60.5357,   66.5125,   72.9866,
     $              79.9921,   87.5648,   95.7422,  104.5641,  114.0721,
     $             124.3099,  135.3236,  147.1615,  159.8742,  173.5149,
     $             188.1391,  203.8052,  220.5742,  238.5099,  257.6791,
     $             278.1515,  300.0000,  323.3007,  348.1330,  374.5798,
     $             402.7275,  432.6661,  464.4894,  498.2951,  534.1850,
     $             572.2648,  612.6447,  655.4391,  700.7670,  748.7519,
     $             799.5223,  853.2114,  909.9575,  969.9041, 1033.2002,
     $            1100.0000  /
       END 