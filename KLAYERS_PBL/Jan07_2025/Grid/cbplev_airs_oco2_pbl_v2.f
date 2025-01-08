C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              incLAY
C
CF77====================================================================

CROUTINE NAME: airs-oco2-pbl-plevs  (heritage: CBPLEV)

CABSTRACT:
C    Common Block for PLEV

CCALL PROTOCOL: none

CINPUT PARAMETERS: none

COUTPUT PARAMETERS: none

CINPUT/OUTPUT PARAMETERS: none

CRETURN VALUES: none

CPARENT(S):
C    KLAYERS

CROUTINES CALLED: none

CFILES ACCESSED: none

CCOMMON BLOCKS: none

CDESCRIPTION:
C    Assign values to PLEV

CALGORITHM REFERENCES: see DESCRIPTION

CKNOWN BUGS AND LIMITATIONS: none

CROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 15 Nov 2000 Scott Hannon      PLEV moved from klayers.f to here
C 20 Nov 2001 Scott Hannon      Add NAMGRD
C    Nov 2024 C.L.Hepplewhite   redifened levels for airs-oco2-pbl (JPL)

CEND====================================================================

C      =================================================================
       BLOCK DATA CBPLEV
C      =================================================================
C
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE

C
C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER I  ! generic looping variable

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
C      template    /'1234567890123456789012345678901234567890123456789'/
       DATA NAMGRD /'AIRS-OCO2 PBL 100 layer grid, 0.005 to 1044 mb '/
C
       DATA  ( PLEV(I), I = 101, 52, -1 )
     $   /  0.0050,     0.0215,     0.0610,     0.1370,     0.2591,
     $      0.4458,     0.7140,     1.0726,     1.5457,     2.1526,
     $      2.8990,     3.8148,     4.9204,     6.2173,     7.7409,
     $      9.5119,    11.5270,    13.8267,    16.4318,    19.3334,
     $     22.5769,    26.1829,    30.1363,    34.4876,    39.2566,
     $     44.4219,    50.0386,    56.1260,    62.6555,    69.6874,
     $     77.2396,    85.2768,    93.8637,   103.0172,   112.6950,
     $    122.9662,   133.8462,   145.2858,   157.3585,   170.0784,
     $    183.3888,   197.3677,   212.0277,   227.3051,   243.2818,
     $    259.9691,   277.2961,   295.3488,   314.1369,   333.5824 /
       DATA  ( PLEV(I), I = 51, 1, -1 )
     $   / 353.7755,   374.7241,   396.3433,   418.7271,   441.8819,
     $     465.7160,   490.3271,   515.7200,   541.7963,   568.6572,
     $     596.3062,   624.6382,   653.7580,   683.6673,   692.6160,
     $     701.6819,   710.8045,   719.9743,   729.2623,   738.5444,
     $     747.9356,   757.3921,   766.8860,   776.4989,   786.1133,
     $     795.8282,   805.6172,   815.4336,   825.3696,   835.3147,
     $     845.3513,   855.4710,   865.6077,   875.8645,   886.1383,
     $     896.4941,   906.9424,   917.3968,   927.9716,   938.5716,
     $     949.2438,   960.0179,   970.7868,   981.6766,   992.5999,
     $    1003.5851,  1014.6820,  1025.7620,  1036.9630,  1048.2863,
     $    1059.700 /
       END
