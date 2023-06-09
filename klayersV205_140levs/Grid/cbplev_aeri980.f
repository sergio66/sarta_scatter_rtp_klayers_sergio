C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              CBPLEV_aeri980
C
!F77====================================================================

!ROUTINE NAME: CBPLEV

!ABSTRACT:
C    Common Block for PLEV
C    This version for AERI layers starting at 980 mb

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
C 19 Jul 2001 Scott Hannon      Create version for AERI 980 mb
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
       DATA NAMGRD /'AERI980 100 layer grid, 0.005 - 980 mb  '/
C
       DATA  ( PLEV(I), I = 101, 52, -1 )
     $    /   0.0050,    0.0122,    0.0273,    0.0561,    0.1073,
     $        0.1921,    0.3240,    0.5182,    0.7902,    1.1542,
     $        1.6223,    2.2025,    2.8986,    3.7096,    4.6299,
     $        5.6497,    6.7558,    8.0366,    9.5120,   11.2033,
     $       13.1326,   15.3231,   17.7989,   20.5849,   23.7064,
     $       27.1894,   31.0600,   35.3445,   40.0690,   45.2595,
     $       50.9416,   57.1401,   63.8790,   71.1815,   79.0694,
     $       87.5632,   96.6820,  106.4432,  116.8624,  127.9534,
     $      139.7279,  152.1956,  165.3641,  179.2387,  193.8226,
     $      209.1166,  225.1195,  241.8276,  259.2352,  277.3345/
       DATA  ( PLEV(I), I = 51, 1, -1 )
     $    / 296.1155,  315.5661,  335.6725,  356.4189,  377.7877,
     $      399.7600,  422.3150,  445.4307,  469.0838,  493.2499,
     $      517.9035,  543.0182,  568.5670,  594.5221,  620.8551,
     $      647.5375,  674.5403,  701.8345,  729.3908,  754.7327,
     $      777.9424,  799.1222,  818.3871,  835.8602,  851.6678,
     $      865.9365,  878.7904,  890.3493,  900.7272,  910.0319,
     $      918.3641,  925.8171,  932.4774,  938.4240,  943.7294,
     $      948.4595,  952.6741,  956.4275,  959.7685,  962.7412,
     $      965.3852,  967.7361,  969.8257,  971.6826,  973.3324,
     $      974.7978,  976.0992,  977.2548,  978.2808,  979.1916,
     $      980.0000/
       END
