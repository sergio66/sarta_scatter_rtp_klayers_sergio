C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County [UMBC]
C
C    AIRS
C
C    CCPREP
C
c copied from home/sergio/SARTA_CLOUDY/v107_Src_pclsam_ctype_100layer_rtpV201/SAFECOPY/ccprep.f 
!F77====================================================================


!ROUTINE NAME:
C    CCPREP


!ABSTRACT:
C    Prepare 100 layer complex cloud for RTA calculation.


!CALL PROTOCOL:
C    CCPREP( NCHAN, LBOT, INDMIE, MIENPS, CNGWAT, CPSIZE, 
C       MIEPS, MIEABS, MIEEXT, MIEASY,
C       LBOTC, LTOPC, G_ASYM, NEXTB, NSCAB, xG_ASYM, xNEXTB, xNSCAB,
C       TEMP, CTYPE)


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
C new
C    REAL arr  TEMP    Layer temperature           K
C    INTEGER CTYPE         ! cloud1 type code number

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
       SUBROUTINE CCPREP_PROFILE( NCHAN, LBOT, INDMIE, MIENPS, 
     $    CNGWAT, CPSIZE, MIEPS, MIEABS, MIEEXT, MIEASY,
     $    LBOTC, LTOPC, 
     $    G_ASYM, NEXTB, NSCAB, 
     $    xG_ASYM, xNEXTB, xNSCAB, 
     $    TEMP,  CTYPE)
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
c new
       REAL TEMP(MAXLAY)           ! layer temperatures
       INTEGER CTYPE               ! cloud1 type code number

C      Output
       INTEGER  LBOTC         ! layer containing cloud bottom
       INTEGER  LTOPC         ! layer containing cloud top
       REAL G_ASYM(MXCHAN)    ! "g" asymmetry       
       REAL  NEXTB(MXCHAN)    ! nadir extinction    
       REAL  NSCAB(MXCHAN)    ! nadir scattering    
       REAL xG_ASYM(MAXLAY,MXCHAN)    ! "g" asymmetry       
       REAL  xNEXTB(MAXLAY,MXCHAN)    ! nadir extinction    
       REAL  xNSCAB(MAXLAY,MXCHAN)    ! nadir scattering    

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
       REAL      XCPSIZE      ! generic size
       REAL      TCLD         ! cld layer temp
       REAL  xC0,xC1,xC2,xC3  ! Liou/Oh ice particle size

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
       IHI = -9999
       ILO = -9999
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

ccc uncomment for testing
c       print *,'LTOPC,LBOTC,CTYPE,INDMIE = ',LTOPC,LBOTC,CTYPE,INDMIE
c       print *,CPSIZE,NPS,MIEPS(1,INDMIE),MIEPS(NPS,INDMIE)
c       print *,ILO,IHI
c       print *,(MIEPS(I,INDMIE),I=1,NPS)
c       IF (ILO .LT. 0) THEN
c         print *,'ext,abs,asym = ',MIEEXT(903,5,INDMIE),
c     $         MIEABS(903,5,INDMIE),MIEASY(903,5,INDMIE)
c       ELSE
c         print *,X,MIEPS(ILO,INDMIE),CPSIZE,MIEPS(IHI,INDMIE),ILO,IHI
c         print *,'LO ext,abs,asym = ',ILO,MIEEXT(903,ILO,INDMIE),
c     $         MIEABS(903,ILO,INDMIE),MIEASY(903,ILO,INDMIE)
c         print *,'HI ext,abs,asym = ',IHI,MIEEXT(903,IHI,INDMIE),
c     $         MIEABS(903,IHI,INDMIE),MIEASY(903,IHI,INDMIE)
c      END IF
ccc
C

C      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c      now do the looping over layer temps
       IF (CTYPE .GE. 300) THEN
         !! this is aerosol, assume particle size indpt of layer temp
         DO I = 1,NCHAN
           DO L = 1,LBOT
             xNEXTB(L,I) = NEXTB(I)
             xNSCAB(L,I) = NSCAB(I)
             xG_ASYM(L,I) = G_ASYM(I)
c             IF (I .EQ. 903) THEN
c               PRINT *,L,I,xNEXTB(L,I),xNSCAB(L,I),xG_ASYM(L,I)
c             END IF
           END DO
         END DO

       ELSEIF ((CTYPE .GE. 100) .AND. (CTYPE .LT. 300)) THEN
         !! this is ice,    particle size depend on layer temp
         !! this is water,  particle size depend on layer temp

         !! ice coeffs from coefficents from S-C Ou, K-N. Liou, Atmospheric Research
         !! 35(1995):127-138
         xC0 = 326.3
         xC1 = 12.42
         xC2 = 0.197
         xC3 = 0.0012

         DO I = 1,NCHAN
           DO L = 1,LBOT
             xNEXTB(L,I) = NEXTB(I)
             xNSCAB(L,I) = NSCAB(I)
             xG_ASYM(L,I) = G_ASYM(I)
           END DO
         END DO

         DO L = 1,LBOT
           XCPSIZE = CPSIZE
           IF ((CTYPE .GE. 100) .AND. (CTYPE .LT. 200)) THEN
             XCPSIZE = 20.0      !! water
           ELSEIF ((CTYPE .GE. 200) .AND. (CTYPE .LT. 300)) THEN
             TCLD = TEMP(L) - 273.16
             XCPSIZE = 60.0      !! ice
             IF (TCLD .LT. -60) TCLD = -60
             IF (TCLD .GT. -20) TCLD = -20
             XCPSIZE = xC0 + xC1 * TCLD + xC2 * (TCLD ** 2) + 
     $                 xC3 * TCLD * (TCLD**2)
           ENDIF 
c           print *,CTYPE,L,TEMP(L) - 273.16,XCPSIZE
           IF (XCPSIZE .LE. MIEPS(1,INDMIE)) THEN
             ! Minimum particle size
             DO I=1,NCHAN
               xNEXTB(L,I)=MIEEXT(I,1,INDMIE)
               xNSCAB(L,I)=MIEEXT(I,1,INDMIE) - MIEABS(I,1,INDMIE)
               xG_ASYM(L,I)=MIEASY(I,1,INDMIE)
             ENDDO
           ELSEIF (XCPSIZE .GE. MIEPS(NPS,INDMIE)) THEN
             ! Maximum particle size
             DO I=1,NCHAN
               xNEXTB(L,I)=MIEEXT(I,NPS,INDMIE)
               xNSCAB(L,I)=MIEEXT(I,NPS,INDMIE) - MIEABS(I,NPS,INDMIE)
               xG_ASYM(L,I)=MIEASY(I,NPS,INDMIE)
             ENDDO
           ELSE
             ! Intermediate particle size
             IHI=1
 20          IF (MIEPS(IHI,INDMIE) .LT. XCPSIZE .AND. IHI .LT. NPS) THEN
               IHI=IHI + 1
               GOTO 20
             ENDIF
             ILO=IHI - 1
             X=( LOG(XCPSIZE) - LOG(MIEPS(ILO,INDMIE)) ) /
     $          ( LOG(MIEPS(IHI,INDMIE)) - LOG(MIEPS(ILO,INDMIE)) )
             DO I=1,NCHAN
               xNEXTB(L,I) =MIEEXT(I,ILO,INDMIE) +
     $           X*(MIEEXT(I,IHI,INDMIE) - MIEEXT(I,ILO,INDMIE))
               ABSB     =MIEABS(I,ILO,INDMIE) +
     $            X*(MIEABS(I,IHI,INDMIE) - MIEABS(I,ILO,INDMIE))
               xNSCAB(L,I) =xNEXTB(L,I) - ABSB
               xG_ASYM(L,I)=MIEASY(I,ILO,INDMIE) +
     $            X*(MIEASY(I,IHI,INDMIE) - MIEASY(I,ILO,INDMIE))
c               IF (I .EQ. 903) THEN
c                 PRINT *,L,I,xNEXTB(L,I),xNSCAB(L,I),xG_ASYM(L,I)
c               END IF
             ENDDO   !! loop over chans
           END IF    !! if beginning, end or mod table
         ENDDO       !! loop over layers
       ENDIF         !! if INDMI = 101,201 or 301
C

C      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       RETURN
       END
