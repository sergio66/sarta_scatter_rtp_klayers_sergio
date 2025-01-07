C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              CBPLEV_aeri961
C
!F77====================================================================

!ROUTINE NAME: CBPLEV

!ABSTRACT:
C    Common Block for PLEV
C    This version for AERI layers starting at 961 mb

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
C 19 Jul 2001 Scott Hannon      Create version for AERI 961 mb
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
       DATA NAMGRD /'AERI961 100 layer grid, 0.005 - 961 mb  '/
C
       DATA  ( PLEV(I), I = 101, 52, -1 )
     $    /   0.0050,    0.0122,    0.0272,    0.0559,    0.1068,
     $        0.1910,    0.3218,    0.5144,    0.7838,    1.1442,
     $        1.6073,    2.1810,    2.8691,    3.6705,    4.5794,
     $        5.5862,    6.6780,    7.9418,    9.3974,   11.0653,
     $       12.9675,   15.1267,   17.5666,   20.3115,   23.3863,
     $       26.8163,   30.6272,   34.8447,   39.4945,   44.6019,
     $       50.1918,   56.2887,   62.9160,   70.0962,   77.8506,
     $       86.1994,   95.1610,  104.7524,  114.9889,  125.8837,
     $      137.4483,  149.6921,  162.6223,  176.2440,  190.5603,
     $      205.5718,  221.2772,  237.6728,  254.7530,  272.5097/
       DATA  ( PLEV(I), I = 51, 1, -1 )
     $    / 290.9333,  310.0119,  329.7316,  350.0771,  371.0310,
     $      392.5747,  414.6877,  437.3485,  460.5341,  484.2207,
     $      508.3832,  532.9957,  558.0317,  583.4639,  609.2648,
     $      635.4061,  661.8596,  688.5967,  715.5891,  740.4108,
     $      763.1429,  783.8859,  802.7529,  819.8643,  835.3442,
     $      849.3168,  861.9036,  873.2220,  883.3839,  892.4947,
     $      900.6531,  907.9506,  914.4717,  920.2941,  925.4886,
     $      930.1199,  934.2464,  937.9212,  941.1923,  944.1028,
     $      946.6914,  948.9930,  951.0389,  952.8569,  954.4721,
     $      955.9068,  957.1810,  958.3124,  959.3168,  960.2086,
     $      961.0000/
       END
