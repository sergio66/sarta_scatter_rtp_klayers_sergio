C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              CBPLEV_910
C
!F77====================================================================

!ROUTINE NAME: CBPLEV

!ABSTRACT:
C    Common Block for PLEV
C    This version for layers 910:1.2:1030

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
C 04 May 2010 Scott Hannon      Created

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
C      so that level 1 is the max pressure and NLEV is the max.
C***********************************************************************
C***********************************************************************
       COMMON /COMLEV/ NAMGRD, PLEV
C
       CHARACTER*40 NAMGRD
       REAL PLEV(101)
C
C      NAMGRD is the name of the layering grid
C      template    /'1234567890123456789012345678901234567890'/
       DATA NAMGRD /'910: 100 layer grid 910:1.2:1030  '/
C
       DATA  ( PLEV(I), I = 101,52,-1 )
     $    / 910.0000,  911.2000,  912.4000,  913.6000,  914.8000,
     $      916.0000,  917.2000,  918.4000,  919.6000,  920.8000,
     $      922.0000,  923.2000,  924.4000,  925.6000,  926.8000,
     $      928.0000,  929.2000,  930.4000,  931.6000,  932.8000,
     $      934.0000,  935.2000,  936.4000,  937.6000,  938.8000,
     $      940.0000,  941.2000,  942.4000,  943.6000,  944.8000,
     $      946.0000,  947.2000,  948.4000,  949.6000,  950.8000,
     $      952.0000,  953.2000,  954.4000,  955.6000,  956.8000,
     $      958.0000,  959.2000,  960.4000,  961.6000,  962.8000,
     $      964.0000,  965.2000,  966.4000,  967.6000,  968.8000/
       DATA  ( PLEV(I), I = 51,1,-1 )
     $    / 970.0000,  971.2000,  972.4000,  973.6000,  974.8000,
     $      976.0000,  977.2000,  978.4000,  979.6000,  980.8000,
     $      982.0000,  983.2000,  984.4000,  985.6000,  986.8000,
     $      988.0000,  989.2000,  990.4000,  991.6000,  992.8000,
     $      994.0000,  995.2000,  996.4000,  997.6000,  998.8000,
     $     1000.0000, 1001.2000, 1002.4000, 1003.6000, 1004.8000,
     $     1006.0000, 1007.2000, 1008.4000, 1009.6000, 1010.8000,
     $     1012.0000, 1013.2000, 1014.4000, 1015.6000, 1016.8000,
     $     1018.0000, 1019.2000, 1020.4000, 1021.6000, 1022.8000,
     $     1024.0000, 1025.2000, 1026.4000, 1027.6000, 1028.8000,
     $     1030.0000/
C
       END