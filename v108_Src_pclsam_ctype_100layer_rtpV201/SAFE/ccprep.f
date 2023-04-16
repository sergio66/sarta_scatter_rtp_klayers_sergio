C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CCPREP
C
!F77====================================================================


!ROUTINE NAME:
C    CCPREP


!ABSTRACT:
C    Prepare 100 layer complex cloud for RTA calculation.


!CALL PROTOCOL:
C    CCPREP( NCHAN, LBOT, INDMIE, MIENPS, CNGWAT, CPSIZE, 
C       MIEPS, MIEABS, MIEEXT, MIEASY,
C       LBOTC, LTOPC, G_ASYM, NEXTB, NSCAB )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   NCHAN   number of channels          none
C    INTEGER   LBOT    bottom layer                none
C    INTEGER   INDMIE  index into MIE arrays       none
C    INT arr   MIENPS  # of particle sizes         none
C    REAL      CNGWAT  cloud non-gas water         g/m^2
C    REAL      CPSIZE  cloud particle size         um
C    REAL arr  MIEPS   Mie table particle sizes    um
C    REAL arr  MIEABS  Mie table absorption data   m^2/g
C    REAL arr  MIEEXT  Mie table extinction data   m^2/g
C    REAL arr  MIEASY  Mie table asymmetry data    ?


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER  LBOTC    layer containing cloud bot  none
C    INTEGER  LTOPC    layer containing cloud top  none
C    REAL arr G_ASYM   "g" asymmetry at cpsize     ?
C    REAL arr NEXTB    extinction at cpsize        m^2/g
C    REAL arr NSCAB    scattering at cpsize        m^2/g


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    SARTA


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    Calculates the transmission thru a cloud



!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date        Programmer     Comments
C    ----------- -------------- ----------------------------------------
C    21 Feb 2007 Scott Hannon   created (based on ccprep.f for slab cld)


!END====================================================================

C      =================================================================
       SUBROUTINE CCPREP( NCHAN, LBOT, INDMIE, MIENPS, CNGWAT, CPSIZE, 
     $    MIEPS, MIEABS, MIEEXT, MIEASY,
     $    LBOTC, LTOPC, G_ASYM, NEXTB, NSCAB )
C      =================================================================

C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
       include 'incFTC.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      QIKEXP


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input
       INTEGER  NCHAN              ! number of channels
       INTEGER   LBOT              ! bottom layer
       INTEGER INDMIE              ! index of CTYPE in MIE arrays
       INTEGER MIENPS(NMIETY)      ! # of particle sizes
       REAL CNGWAT(MAXLAY)         ! cloud non-gas water/whatever
       REAL CPSIZE                 ! cloud particle size
       REAL  MIEPS(MXMIEA,NMIETY)  ! particle size
       REAL MIEABS(MXCHAN,MXMIEA,NMIETY) ! scattering absorption
       REAL MIEEXT(MXCHAN,MXMIEA,NMIETY) ! scattering extinction
       REAL MIEASY(MXCHAN,MXMIEA,NMIETY) ! scattering asymmetry

C      Output
       INTEGER  LBOTC         ! layer containing cloud bottom
       INTEGER  LTOPC         ! layer containing cloud top
       REAL G_ASYM(MXCHAN)    ! "g" asymmetry
       REAL  NEXTB(MXCHAN)    ! nadir extinction
       REAL  NSCAB(MXCHAN)    ! nadir scattering


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      I         ! looping variable for channel
       INTEGER    IHI         ! high index
       INTEGER    ILO         ! low index
       INTEGER      L         ! looping variable for layer
       INTEGER     LR         ! reversed layer index
       INTEGER    NPS         ! # of particle sizes for this CTYPE
       REAL   ABSB            ! interpolated absorption
       REAL      X            ! generic real variable


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
C      Determine top and bottom layers of cloud
       LTOPC=0
       LBOTC=0
       DO L=1,LBOT
          LR = LBOT + 1 - L
          IF (LTOPC .EQ. 0 .AND. CNGWAT(L)  .GT. CNGMIN) LTOPC=L
          IF (LBOTC .EQ. 0 .AND. CNGWAT(LR) .GT. CNGMIN) LBOTC=LR
       ENDDO
C

C      Interpolate Mie tables for particle size
C      Note: extinction = scattering + absorption, so sca=ext - abs
C
C      Number of particle sizes for current CTYPE
       NPS=MIENPS(INDMIE)
C
C      Minimum particle size
       IF (CPSIZE .LE. MIEPS(1,INDMIE)) THEN
          DO I=1,NCHAN
             NEXTB(I)=MIEEXT(I,1,INDMIE)
             NSCAB(I)=MIEEXT(I,1,INDMIE) - MIEABS(I,1,INDMIE)
             G_ASYM(I)=MIEASY(I,1,INDMIE)
          ENDDO
C
C      Maximum particle size
       ELSEIF (CPSIZE .GE. MIEPS(NPS,INDMIE)) THEN
          DO I=1,NCHAN
             NEXTB(I)=MIEEXT(I,NPS,INDMIE)
             NSCAB(I)=MIEEXT(I,NPS,INDMIE) - MIEABS(I,NPS,INDMIE)
             G_ASYM(I)=MIEASY(I,NPS,INDMIE)
          ENDDO
C
C      Intermediate particle size
       ELSE
          IHI=1
 10       IF (MIEPS(IHI,INDMIE) .LT. CPSIZE .AND. IHI .LT. NPS) THEN
             IHI=IHI + 1
             GOTO 10
          ENDIF
          ILO=IHI - 1
C
          X=( LOG(CPSIZE) - LOG(MIEPS(ILO,INDMIE)) ) /
     $      ( LOG(MIEPS(IHI,INDMIE)) - LOG(MIEPS(ILO,INDMIE)) )
C
          DO I=1,NCHAN
             NEXTB(I) =MIEEXT(I,ILO,INDMIE) +
     $          X*(MIEEXT(I,IHI,INDMIE) - MIEEXT(I,ILO,INDMIE))
             ABSB     =MIEABS(I,ILO,INDMIE) +
     $          X*(MIEABS(I,IHI,INDMIE) - MIEABS(I,ILO,INDMIE))
             NSCAB(I) =NEXTB(I) - ABSB
             G_ASYM(I)=MIEASY(I,ILO,INDMIE) +
     $          X*(MIEASY(I,IHI,INDMIE) - MIEASY(I,ILO,INDMIE))
          ENDDO
       ENDIF
C

       RETURN
       END
