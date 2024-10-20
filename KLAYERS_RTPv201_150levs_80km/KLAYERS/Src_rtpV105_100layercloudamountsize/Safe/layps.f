C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              LAYPS
C
!F77====================================================================

!ROUTINE NAME: LAYPS

!ABSTRACT:
C    Calculate mean layer cloud particle size.

!CALL PROTOCOL:
C    LAYPS(NPS, INDPS, LAYBOT, NSUB, PSUB, TSUB, MRSUB, DZSUB, ALAY)

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   NPS     # of particle size profs    none
C    INT arr   INDPS   indices of p.size profs     none
C    INTEGER   LAYBOT  bottom layer number         none
C    INT arr   NSUB    # of sub-levels each layer  none
C    REAL arr  PSUB    avg pres of each sub-layer  mb
C    REAL arr  TSUB    avg air temp of sub-layers  K
C    REAL arr  MRSUB   avg mix ratio each layer    ppmv
C    REAL arr  DZSUB   thickness of each sublayer  m


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  ALAY    avg gas amount, AIRS layer  k.mol/cm^2


!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S): KLAYERS

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C    Integrates (sums) the sub-layers for each layer using averaged
C    sub-layer values for temperature, mixing ratio, and pressure.
C    An average layer particle size is computed using the "air amount"
C    in the sub-layers as a weight.

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS:
C    No error checking of the input profile 

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 31 Aug 2007 Scott Hannon      created based on integ.f


!END====================================================================


C      =================================================================
       SUBROUTINE LAYPS(NPS, INDPS, LAYBOT, NSUB,
     $    PSUB, TSUB, MRSUB, DZSUB, ALAY)
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
C      Input parameters:
       INTEGER    NPS
       INTEGER  INDPS(MXCLD)
       INTEGER LAYBOT
       INTEGER   NSUB(MYNLAY)
       REAL   PSUB(NSUBLV-1)
       REAL   TSUB(NSUBLV-1)
       REAL  MRSUB(NSUBLV-1,MXGAS)
       REAL  DZSUB(NSUBLV-1)
C
C      Output parameters:
       REAL   ALAY(MYNLAY,MXGAS)

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       REAL   ASUB
       REAL   ASUM
       REAL  PSSUM(MXCLD)
       INTEGER      I
       INTEGER      J
       INTEGER      K
       INTEGER      L
       INTEGER      M

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************
C
C      -----------------------------
C      Loop over the 100 AIRS layers
C      -----------------------------
       IF (NPS .GT. 0) THEN
          J=0
          DO L=LAYBOT,MYNLAY
C
             ASUM=0.0
             DO M=1,NPS
                PSSUM(M)=0.0
             ENDDO

C            ------------------------
C            Loop over the sub-layers
C            ------------------------
             DO I=1,NSUB(L)

C               sublay index
                K=I+J

C               sublay gas amount (excluding constants)
                ASUB=DZSUB(K)*PSUB(K)/TSUB(K)

C               increment layer sums
                ASUM=ASUM + ASUB
                DO M=1,NPS
                   PSSUM(M)=PSSUM(M) + MRSUB(K,INDPS(M))*ASUB
                ENDDO
C
             ENDDO
C
C            layer mean particle size
             DO M=1,NPS
                ALAY(L,INDPS(M))=PSSUM(M)/ASUM
             ENDDO
C
             J=J+NSUB(L)
          ENDDO
C
       ENDIF
C
       RETURN
       END
