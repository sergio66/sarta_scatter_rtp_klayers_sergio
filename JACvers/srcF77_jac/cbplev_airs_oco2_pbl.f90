!=======================================================================
!
!              University of Maryland Baltimore County [UMBC]
!
!              AIRS
!
!              incLAY
!
!F90====================================================================

!ROUTINE NAME: airs-oco2-pbl-plevs  (heritage: CBPLEV)

!ABSTRACT:
!    Common Block for PLEV

!CALL PROTOCOL: none

!INPUT PARAMETERS: none

!OUTPUT PARAMETERS: none

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S):
!    KLAYERS

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
!    Assign values to PLEV

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS: none

!ROUTINE HISTORY:
!    Date     Programmer        Comments
!------------ ----------------- ----------------------------------------
! 15 Nov 2000 Scott Hannon      PLEV moved from klayers.f to here
! 20 Nov 2001 Scott Hannon      Add NAMGRD
!    Nov 2024 C.L.Hepplewhite   redifened levels for airs-oco2-pbl (JPL)

!END====================================================================

!-----------------------------------------------------------------------
!      EXECUTABLE CODE
!-----------------------------------------------------------------------
!      none; common block data
!
!      Edit as needed, but PLEV must match NBPLEV in incLAY.f
!      PLEV must be in mb or hPa (these are equivalent), and sorted
!      so that level 1 is the max pressure and NLEV is the min.
!***********************************************************************
!***********************************************************************
       COMMON /COMLEV/ NAMGRD, PLEV
!
       CHARACTER*50 NAMGRD
       REAL PLEV(101)
!
!      NAMGRD is the name of the layering grid
!      template    /'1234567890123456789012345678901234567890123456789'/
       DATA NAMGRD /'AIRS-OCO2 PBL 100 layer grid, 0.005 to 1044 mb '/
!
       DATA  ( PLEV(I), I = 101, 52, -1 )
      $   /   0.0050,     0.0161,     0.0384,     0.1370,     0.2244,
      $       0.5064,     0.7140,     1.2972,     1.6872,     2.7009,
      $       3.3398,     4.9204,     5.8776,     8.1655,     9.5119,
      $      12.6492,    14.4559,    18.5847,    20.9224,    26.1829,
      $      32.2744,    35.6505,    43.1001,    47.1882,    51.5278,
      $      60.9895,    66.1253,    71.5398,    77.2396,    83.2310,
      $      89.5204,    96.1138,   103.0172,   110.2366,   117.7775,
      $     125.6456,   133.8462,   142.3848,   151.2664,   160.4959,
      $     170.0784,   180.0183,   190.3203,   200.9887,   212.0277,
      $     223.4415,   235.2338,   247.4085,   259.9691,   272.9191/
       DATA  ( PLEV(I), I = 51, 1, -1 )
      $   / 286.2617,   300.0000,   314.1369,   328.6753,   343.6176,
      $     358.9665,   374.7241,   390.8926,   407.4738,   424.4698,
      $     441.8819,   459.7118,   477.9607,   496.6298,   515.7200,
      $     535.2322,   555.1669,   575.5248,   596.3062,   617.5112,
      $     639.1398,   661.1920,   683.6673,   706.5654,   729.8857,
      $     740.9087,   751.9317,   763.1198,   774.3380,   785.6910,
      $     797.1039,   808.6214,   820.2285,   831.9100,   843.7106,
      $     855.5557,   867.5492,   879.5574,   891.7430,   903.9286,
      $     916.2906,   928.6674,   941.1906,   953.7576,   966.4412,
      $     979.1977,   992.0411,  1004.9859,  1017.9882,  1031.1202,
      $    1044.2522 / 
       END
