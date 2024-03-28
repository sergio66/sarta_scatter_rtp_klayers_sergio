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
C    Prepapre lookup table etc for a complex cloud calculation.


!CALL PROTOCOL:
C    CCPREP( NCHAN, LBOT, INDMIE, MIENPS,
C       CNGWAT, CPSIZE, CPRTOP, CPRBOT, PLEV, TEMP, SECANG, SECSUN,
C       MIEPS, MIEABS, MIEEXT, MIEASY, LCBOT, LCTOP, CLEARB, CLEART,
C       TCBOT, TCTOP, MASEC, MSSEC, CFRCL, G_ASYM, NEXTOD, NSCAOD, 
C       POLYNOM_BACKSCAT(4), ISCALING, CTYPE,
C       DOJAC, JACA_G_ASY, JACA_NEXTO, JACA_NSCAO, JACA_FINAL, 
C              JACS_G_ASY, JACS_NEXTO, JACS_NSCAO, JACS_FINAL,
C              JACTOP_CFRCL,JACBOT_CFRCL,JACTOP_CFRCL_v,JACBOT_CFRCL_v)

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   NCHAN   number of channels          none
C    INTEGER   LBOT    bottom layer                none
C    INTEGER   INDMIE  index into MIE arrays       none
C    INT arr   MIENPS  # of particle sizes         none
C    REAL      CNGWAT  cloud non-gas water         g/m^2
C    REAL      CPSIZE  cloud particle size         um
C    REAL      CPRTOP  cloud top pressure          mb
C    REAL      CPRBOT  cloud bottom pressure       mb
C    REAL arr  PLEV    layer pres boundary levels  mb
C    REAL arr  TEMP    average layer temperature   K
C    REAL arr  SECANG  path secant angles          none
C    REAL arr  SECSUN  sun path secant angles      none
C    REAL arr  MIEPS   Mie table particle sizes    um
C    REAL arr  MIEABS  Mie table absorption data   m^2/g
C    REAL arr  MIEEXT  Mie table extinction data   ?
C    REAL arr  MIEASY  Mie table asymmetry data    ?
C    INTEGER    CTYPE  cloud type


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER  LCBOT    layer containing cloud bottom
C    INTEGER  LCTOP    layer containing cloud top
C    REAL     CLEARB   frac of layer at bottom of cloud clear
C    REAL     CLEART   frac of layer at top of cloud clear
C    REAL     TCBOT    temperature at cloud bottom
C    REAL     TCTOP    temperature at cloud top
C    REAL     MASEC    mean cloud view angle secant
C    REAL     MSSEC    mean cloud sun-only angle secant
C    REAL arr CFRCL    fraction of cloud in layer
C    REAL arr G_ASYM   "g" asymmetry
C    REAL arr NEXTOD   nadir extinction optical depth
C    REAL arr NSCAOD   nadir scattering optical depth
C    REAL arr POLYNOM_BACKSCAT back scatter coeffs

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
C    23 Jan 2004 Scott Hannon   Created from a re-write of calcc1 to
C                                  output results and not call calcc2.
C    31 Mar 2006 Scott Hannon   Revised for flexible CTYPE; add INDMIE
C                               and MIENPS.
C    26 Apr 2006 Scott Hannon   Add LBLACK argument and "if" block.
C    14 Nov 2007 Scott Hannon   Remove LBLACK

!END====================================================================

C      =================================================================
       SUBROUTINE CCPREP( NCHAN, LBOT, INDMIE, MIENPS,
     $    CNGWAT, CPSIZE, CPRTOP, CPRBOT, PLEV, TEMP, SECANG, SECSUN,
     $    MIEPS, MIEABS, MIEEXT, MIEASY,
     $    LCBOT, LCTOP, CLEARB, CLEART, TCBOT, TCTOP, MASEC, MSSEC,
     $    CFRCL, G_ASYM, NEXTOD, NSCAOD, POLYNOM_BACKSCAT, ISCALING, CTYPE,
     $    DOJAC, JACA_G_ASYM, JACA_NEXTOD, JACA_NSCAOD, JACA_FINAL, 
     $           JACS_G_ASYM, JACS_NEXTOD, JACS_NSCAOD, JACS_FINAL,
     $           JACTOP_CFRCL,JACBOT_CFRCL,JACTOP_CFRCL_v,JACBOT_CFRCL_v)

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
       INTEGER CTYPE               ! 101-199 = W, 201-299 = I, 301-399 = A
       REAL CNGWAT                 ! cloud non-gas water
       REAL CPSIZE                 ! cloud particle size
       REAL CPRTOP                 ! cloud top pressure
       REAL CPRBOT                 ! cloud bottom pressure
       REAL   PLEV(MAXLAY+1)       ! pressure levels
       REAL   TEMP(MAXLAY)         ! temperature
       REAL SECANG(MAXLAY)         ! secant of view path
       REAL SECSUN(MAXLAY)         ! secant of total sun path
       REAL  MIEPS(MXMIEA,NMIETY)  ! particle size
       REAL MIEABS(MXCHAN,MXMIEA,NMIETY) ! scattering absorption
       REAL MIEEXT(MXCHAN,MXMIEA,NMIETY) ! scattering extinction
       REAL MIEASY(MXCHAN,MXMIEA,NMIETY) ! scattering asymmetry
       LOGICAL DOJAC

C      Output
       INTEGER  LCBOT            ! layer containing cloud bottom
       INTEGER  LCTOP            ! layer containing cloud top
       REAL CLEARB               ! frac of layer at bottom of cloud clear
       REAL CLEART               ! frac of layer at top of cloud clear
       REAL  TCBOT               ! temperature at cloud bottom
       REAL  TCTOP               ! temperature at cloud top
       REAL  MASEC               ! mean cloud view angle secant
       REAL  MSSEC               ! mean cloud sun-only angle secant
       REAL  CFRCL(MAXLAY)       ! fraction of cloud in layer
       REAL G_ASYM(MXCHAN)       ! "g" asymmetry
       REAL NEXTOD(MXCHAN)       ! nadir extinction optical depth
       REAL NSCAOD(MXCHAN)       ! nadir scattering optical depth
       REAL POLYNOM_BACKSCAT(4)  ! backscattering coefs : similarity 1999, Chou 1999 or Maestri/Martinazzi JSQRT 2021
       INTEGER ISCALING          ! -1 for similarity, -2 for chou, -3 for Maestri/Martinazzo
       REAL JACA_G_ASYM(MXCHAN),JACA_NEXTOD(MXCHAN),JACA_NSCAOD(MXCHAN),JACA_FINAL(MXCHAN)  !! amount jacs
       REAL JACS_G_ASYM(MXCHAN),JACS_NEXTOD(MXCHAN),JACS_NSCAOD(MXCHAN),JACS_FINAL(MXCHAN)  !! sze jacs
       REAL JACTOP_CFRCL(MAXLAY),JACBOT_CFRCL(MAXLAY)    ! derivatives of fraction of cloud in layer
       REAL JACTOP_CFRCL_v(MAXLAY,MXCHAN),JACBOT_CFRCL_v(MAXLAY,MXCHAN)    ! spectral derivatives of fraction of cloud in layer
C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      I         ! looping variable for channel
       INTEGER    IHI         ! high index
       INTEGER    ILO         ! low index
       INTEGER      L         ! looping variable for layer
       INTEGER     LR         ! reversed layer index
       INTEGER    NPS         ! # of particle sizes for this CTYPE
       REAL  ABSOD            ! interpolated absorption optical depth
       REAL   PAVG            ! layer average pressure
       REAL  PAVG2            ! adjacent layer average pressure
       REAL      X            ! generic junk real variable
       CHARACTER*50 caStr     ! junk

c to help with jacs
       REAL  ASYM             ! interpolated asymmetry
       REAL  JACABSOD,DX      ! temporary variables for jacobians


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
C which output variable depends on cprtop/cprbot?
C   LCTOP/LCBOT           -- used by docloudyTwoSlab_RT.f for black/TwoSlab clouds
C   CFRCL(MAXLAY)         -- used by docloudyTwoSlab_RT.f for TwoSLab clouds
C   CLEART/CLEARB         -- used by docloudyTwoSlab_RT.f for black clouds 
C   TCTOP only on CPRTOP  -- used by docloudyTwoSlab_RT.f for black clouds 
C
C***********************************************************************
C***********************************************************************

C      --------------------------------
C      Find top and bottom cloud layers
C      --------------------------------
       DO L=1,LBOT
          CFRCL(L)=0.0 ! initialize to zero
C          LR = MAXLAY + 1 - L ! replaced by line below 21 July 2003
          LR = LBOT + 1 - L
          IF (PLEV(L) .LE. CPRTOP) LCTOP=L
          IF (PLEV(LR+1) .GE. CPRBOT) LCBOT=LR
       ENDDO
C
C      Calc fraction of layer at top & bottom of cloud that is
C      clear (of this cloud; there may be another cloud there).
       CLEART=(CPRTOP   - PLEV(LCTOP))/(PLEV(LCTOP+1) - PLEV(LCTOP))
       CLEARB=(PLEV(LCBOT+1) - CPRBOT)/(PLEV(LCBOT+1) - PLEV(LCBOT))
C

C      --------------------------
C      Calc cloud top temperature
C      --------------------------
       L=LCTOP
       PAVG=(PLEV(L+1) - PLEV(L))/LOG( PLEV(L+1)/PLEV(L) )
       IF (PAVG .GT. CPRTOP .OR. L .EQ. LBOT) THEN
          PAVG2=(PLEV(L) - PLEV(L-1))/LOG( PLEV(L)/PLEV(L-1) )
          TCTOP=TEMP(L) + LOG(CPRTOP/PAVG)*
     $       (TEMP(L-1) - TEMP(L))/LOG( PAVG2/PAVG )
       ELSE
          PAVG2=(PLEV(L+2) - PLEV(L+1))/LOG( PLEV(L+2)/PLEV(L+1) )
          TCTOP=TEMP(L) + LOG(CPRTOP/PAVG)*
     $       (TEMP(L+1) - TEMP(L))/LOG( PAVG2/PAVG )
       ENDIF
C
ccccccccc
c      print *, 'top PLEV(L-1)=', PLEV(L-1)
c      print *, 'top PLEV(L  )=', PLEV(L)
c      print *, 'top PLEV(L+1)=', PLEV(L+1)
c      print *, 'top PLEV(L+2)=', PLEV(L+2)
c      print *, 'top pavg=', PAVG
c      print *, 'top pavg2=', PAVG2
c      print *, 'top TEMP(L-1)=',TEMP(L-1)
c      print *, 'top TEMP(L  )=',TEMP(L)
c      print *, 'top TEMP(L+1)=',TEMP(L+1)
ccccccccc
C
C      -----------------------------
C      Calc cloud bottom temperature
C      -----------------------------
       L=LCBOT
       PAVG=(PLEV(L+1) - PLEV(L))/LOG( PLEV(L+1)/PLEV(L) )
       IF (PAVG .GT. CPRBOT .OR. L .EQ. LBOT) THEN
          PAVG2=(PLEV(L) - PLEV(L-1))/LOG( PLEV(L)/PLEV(L-1) )
          TCBOT=TEMP(L) + LOG(CPRBOT/PAVG)*
     $       (TEMP(L-1) - TEMP(L))/LOG( PAVG2/PAVG )
       ELSE
          PAVG2=(PLEV(L+2) - PLEV(L+1))/LOG( PLEV(L+2)/PLEV(L+1) )
          TCBOT=TEMP(L) + LOG(CPRBOT/PAVG)*
     $       (TEMP(L+1) - TEMP(L))/LOG( PAVG2/PAVG )
       ENDIF
C
ccccccccc
c      print *, 'bot PLEV(L-1)=', PLEV(L-1)
c      print *, 'bot PLEV(L  )=', PLEV(L)
c      print *, 'bot PLEV(L+1)=', PLEV(L+1)
c      print *, 'bot PLEV(L+2)=', PLEV(L+2)
c      print *, 'bot pavg=', PAVG
c      print *, 'bot pavg2=', PAVG2
c      print *, 'bot TEMP(L-1)=',TEMP(L-1)
c      print *, 'bot TEMP(L  )=',TEMP(L)
c      print *, 'bot TEMP(L+1)=',TEMP(L+1)
ccccccccc

ccccccccc
c      print *, 'lctop, plev(lctop) = ', LCTOP, PLEV(LCTOP)
c      print *, 'lcbot, plev(lcbot+1) = ', LCBOT, PLEV(LCBOT+1)
c      print *, 'tctop=', TCTOP
c      print *, 'tcbot=', TCBOT
c      print *, 'cleart=', CLEART
c      print *, 'clearb=', CLEARB
ccccccccc
       JACTOP_CFRCL = 0
C      -----------------------------------------------------------------
C      Calc mean secant angles thru cloud and fraction of cloud in layer
C      -----------------------------------------------------------------
       IF (LCTOP .EQ. LCBOT) THEN
          MASEC=SECANG(LCTOP)
          MSSEC=SECSUN(LCTOP)
          CFRCL(LCTOP)=1.0
          IF (DOJAC) THEN
            L = LCTOP
            JACTOP_CFRCL(L) = 0;
            JACBOT_CFRCL(L) = 0;
c            write(*,'(I5,6(F12.5))') L,CFRCL(L),CPRBOT,CPRTOP,PLEV(L),JACTOP_CFRCL(L),JACBOT_CFRCL(L)
          END IF
       ELSE
C         top & bottom layers
          MASEC=SECANG(LCTOP)*(1-CLEART) + SECANG(LCBOT)*(1-CLEARB)
          MSSEC=SECSUN(LCTOP)*(1-CLEART) + SECSUN(LCBOT)*(1-CLEARB)
          X=(1-CLEART) + (1-CLEARB)
          CFRCL(LCTOP)=(PLEV(LCTOP+1)-CPRTOP)/(CPRBOT-CPRTOP)
          CFRCL(LCBOT)=(CPRBOT-PLEV(LCBOT))/(CPRBOT-CPRTOP)
          IF (DOJAC) THEN
            L = LCTOP
            JACTOP_CFRCL(L) = 1/(CPRBOT-CPRTOP)/(CPRBOT-CPRTOP) * ( (CPRBOT-CPRTOP)*(-1) - (PLEV(LCTOP+1)-CPRTOP)*(-1) )
            JACBOT_CFRCL(L) = 1/(CPRBOT-CPRTOP)/(CPRBOT-CPRTOP) * ( (CPRBOT-CPRTOP)*(0)  - (PLEV(LCTOP+1)-CPRTOP)*(+1) )
c            write(*,'(I5,6(F12.5))') L,CFRCL(L),CPRBOT,CPRTOP,PLEV(L),JACTOP_CFRCL(L),JACBOT_CFRCL(L)

            L = LCBOT
            JACTOP_CFRCL(L) = 1/(CPRBOT-CPRTOP)/(CPRBOT-CPRTOP) * ( (CPRBOT-CPRTOP)*(0)  - (CPRBOT-PLEV(LCBOT))*(-1) )
            JACBOT_CFRCL(L) = 1/(CPRBOT-CPRTOP)/(CPRBOT-CPRTOP) * ( (CPRBOT-CPRTOP)*(1)  - (CPRBOT-PLEV(LCBOT))*(+1) )
c            write(*,'(I5,6(F12.5))') L,CFRCL(L),CPRBOT,CPRTOP,PLEV(L),JACTOP_CFRCL(L),JACBOT_CFRCL(L)
          END IF

C         other layers
          DO L=LCTOP+1,LCBOT-1
             MASEC=MASEC + SECANG(L)
             MSSEC=MSSEC + SECSUN(L)
             X=X + 1
             CFRCL(L)=(PLEV(L+1)-PLEV(L))/(CPRBOT-CPRTOP)
             IF (DOJAC) THEN
               JACTOP_CFRCL(L) = -(PLEV(L+1)-PLEV(L))/(CPRBOT-CPRTOP)/(CPRBOT-CPRTOP)*(-1)
               JACBOT_CFRCL(L) = -(PLEV(L+1)-PLEV(L))/(CPRBOT-CPRTOP)/(CPRBOT-CPRTOP)*(+1)
               ! see test_cld_jacs_rads_cldfracsL.m
c               write(*,'(I5,6(F12.5))') L,CFRCL(L),CPRBOT,CPRTOP,PLEV(L),JACTOP_CFRCL(L),JACBOT_CFRCL(L)
             END IF
          ENDDO
C         Divide secant sum by weight sum to get mean secant
          MASEC=MASEC/X
          MSSEC=MSSEC/X
       ENDIF
       MSSEC=MSSEC - MASEC ! convert total secant to sun-only secant
C

C backscatter coeffs, see Maestr/Martinazzo Journal of Quantitative Spectroscopy & Radiative Transfer 271 (2021) 107739
c Assessment of the accuracy of scaling methods for radiance simulations at far and mid infrared wavelengths
c Michele Martinazzoa, Davide Magurnoa, William Cossicha, Carmine Seriob, Guido Masiellob, Tiziano Maestri
c       ISCALING = 1     !!! similarity, been using this for years
c       ISCALING = 2     !!! chou
c       ISCALING = 3     !!! Maestri/Martinazzo
       ISCALING = ABS(FLOOR(rXTang))
cc       write(*,'(A,F8.3,A,I3)') 'rXtang in incFTC = ',rXTang,' so ISCALING = ',ISCALING
       IF ((ISCALING .EQ. 0) .AND. (rXTang .LE. 0.0)) THEN
         POLYNOM_BACKSCAT = (/0.5000, 0.5000, 0.0000, 0.0000/)        !!!! similarity scaling , eqn 10
         caStr = ' : SIMILARITY scaling'
       ELSEIF ((ISCALING .EQ. 0) .AND. (rXTang .GT. 0.0)) THEN
         POLYNOM_BACKSCAT = (/0.5000, 0.5000, 0.0000, 0.0000/)        !!!! similarity scaling , eqn 10
         caStr = ' : SIMILARITY scaling w/ Tang Correction'
       ELSEIF ((ISCALING .EQ. 1) .AND. (rXTang .LE. 0.0)) THEN
         POLYNOM_BACKSCAT = (/0.5000, 0.5000, 0.0000, 0.0000/)        !!!! similarity scaling , eqn 10, default till July 2023
         caStr = ' : SIMILARITY scaling'
       ELSEIF ((ISCALING .EQ. 1) .AND. (rXTang .GT. 0.0)) THEN
         POLYNOM_BACKSCAT = (/0.5000, 0.5000, 0.0000, 0.0000/)        !!!! similarity scaling , eqn 10, default till July 2023
         caStr = ' : SIMILARITY scaling w/ Tang Correction'
       ELSEIF ((ISCALING .EQ. 2).AND. (rXTang .LE. 0.0))  THEN
         POLYNOM_BACKSCAT = (/0.5000, 0.3738, 0.0076, 0.1186/)        !!!! chou scaling, Table 3
         caStr = '  : CHOU scaling'
       ELSEIF ((ISCALING .EQ. 2).AND. (rXTang .GT. 0.0))  THEN
         POLYNOM_BACKSCAT = (/0.5000, 0.3738, 0.0076, 0.1186/)        !!!! chou scaling, Table 3
         caStr = '  : CHOU scaling w/ Tang Correction'
       ELSEIF ((ISCALING .EQ. 3) .AND. ((CTYPE .GE. 101) .AND. (CTYPE .LE. 199))) THEN
         POLYNOM_BACKSCAT = (/0.5000, 0.2884, 0.5545,-0.3429/)        !!!! uniboWAT scaling, Table 3
         caStr = ' : Meaestri/Martinazzo water scaling, NO TANG'
       ELSEIF ((ISCALING .EQ. 3) .AND. ((CTYPE .GE. 201) .AND. (CTYPE .LE. 299))) THEN
         POLYNOM_BACKSCAT = (/0.5000, 0.4452,-0.3189, 0.3737/)        !!!! uniboICE scaling, Table 3
         caStr = ' : Meaestri/Martinazzo ice scaling, NO TANG'
       ELSE
         print *,'unknown rXtang!!!! ',rXTang
         STOP
       END IF

c       WRITE(*,'(A,F8.3,I3,A)') 'ccprep_slab.f : rXtang,ISCALING = ',rXtang,ISCALING,caStr
c       WRITE(*,'(A,4(F8.3))')   '  backscatter polynomial = ',POLYNOM_BACKSCAT
      
C      --------------------------------------------------
C      Interpolate tables for particle size and scale for
C      nadir total cloud water
C      -----------------------
C      Note: extinction = scattering + absorption, so sca = ext - abs
C
C      Number of particle sizes for current CTYPE
       NPS=MIENPS(INDMIE)
c
C      *************************
C      Minimum particle size
       IF (CPSIZE .LE. MIEPS(1,INDMIE)) THEN
          !!! need these 3 lines for jacos
          ILO = 1
          IHI = 2
          X=( LOG(CPSIZE) - LOG(MIEPS(ILO,INDMIE)) ) /
     $      ( LOG(MIEPS(IHI,INDMIE)) - LOG(MIEPS(ILO,INDMIE)) )
c          print *,'MIN PARTICLE SIZE',CPSIZE,MIEPS(1,INDMIE),ILO,IHI,X
          DO I=1,NCHAN
             NEXTOD(I)=CNGWAT*MIEEXT(I,1,INDMIE)
             NSCAOD(I)=CNGWAT*(MIEEXT(I,1,INDMIE) - MIEABS(I,1,INDMIE))
             G_ASYM(I)=MIEASY(I,1,INDMIE)
          ENDDO
          IF (DOJAC) THEN
            DO I=1,NCHAN
              JACA_NEXTOD(I)=MIEEXT(I,1,INDMIE)
              JACA_NSCAOD(I)=(MIEEXT(I,1,INDMIE) - MIEABS(I,1,INDMIE))
              JACA_G_ASYM(I)=0

              JACA_FINAL(I) = NEXTOD(I) - NSCAOD(I)*(1.0+G_ASYM(I))/2.0
              DO L = LCTOP,LCBOT
                JACTOP_CFRCL_v(L,I) = JACA_FINAL(I) * JACTOP_CFRCL(L)
                JACBOT_CFRCL_v(L,I) = JACA_FINAL(I) * JACBOT_CFRCL(L)
              END DO
              JACA_FINAL(I) = 1/(CNGWAT+1e-16)*JACA_FINAL(I)
                
               DX = 1 / ( LOG(MIEPS(IHI,INDMIE)) - LOG(MIEPS(ILO,INDMIE)) )
               DX = DX / (CPSIZE+1e-16)   !!!! since we have log(sze) ==> d/dsze = 1/sze
               JACS_NEXTOD(I)=CNGWAT*( 0*MIEEXT(I,ILO,INDMIE) +
     $            DX*(MIEEXT(I,IHI,INDMIE) - MIEEXT(I,ILO,INDMIE)) )
               JACABSOD    =CNGWAT*( 0*MIEABS(I,ILO,INDMIE) +
     $            DX*(MIEABS(I,IHI,INDMIE) - MIEABS(I,ILO,INDMIE)) )
               JACS_NSCAOD(I)=JACS_NEXTOD(I) - JACABSOD
               JACS_G_ASYM(I)=0*MIEASY(I,ILO,INDMIE) +
     $                  DX*(MIEASY(I,IHI,INDMIE) - MIEASY(I,ILO,INDMIE))
               JACS_FINAL(I) = JACS_NEXTOD(I) - 0.5*JACS_NSCAOD(I)*(1+G_ASYM(I)) - 0.5*NSCAOD(I)*JACA_G_ASYM(I)
            ENDDO
          END IF

C      *************************
C      Maximum particle size
       ELSEIF (CPSIZE .GE. MIEPS(NPS,INDMIE)) THEN
          !!! need these 3 lines for jacos
          IHI = NPS
          ILO = IHI - 1
          X=( LOG(CPSIZE) - LOG(MIEPS(ILO,INDMIE)) ) /
     $      ( LOG(MIEPS(IHI,INDMIE)) - LOG(MIEPS(ILO,INDMIE)) )
c          print *,'MAX PARTICLE SIZE',CPSIZE,MIEPS(NPS,INDMIE),ILO,IHI,X
          DO I=1,NCHAN
             NEXTOD(I)=CNGWAT*MIEEXT(I,NPS,INDMIE)
             NSCAOD(I)=CNGWAT*(MIEEXT(I,NPS,INDMIE) - MIEABS(I,NPS,INDMIE))
             G_ASYM(I)=MIEASY(I,NPS,INDMIE)
          ENDDO
          IF (DOJAC) THEN
            DO I=1,NCHAN
              JACA_NEXTOD(I)=MIEEXT(I,NPS,INDMIE)
              JACA_NSCAOD(I)=(MIEEXT(I,NPS,INDMIE) - MIEABS(I,NPS,INDMIE))
              JACA_G_ASYM(I)=0

              JACA_FINAL(I) = NEXTOD(I) - NSCAOD(I)*(1.0+G_ASYM(I))/2.0
              DO L = LCTOP,LCBOT
                JACTOP_CFRCL_v(L,I) = JACA_FINAL(I) * JACTOP_CFRCL(L)
                JACBOT_CFRCL_v(L,I) = JACA_FINAL(I) * JACBOT_CFRCL(L)
              END DO
              JACA_FINAL(I) = 1/(CNGWAT+1e-16)*JACA_FINAL(I)

               DX = 1 / ( LOG(MIEPS(IHI,INDMIE)) - LOG(MIEPS(ILO,INDMIE)) )
               DX = DX / (CPSIZE+1e-16)   !!!! since we have log(sze) ==> d/dsze = 1/sze
               JACS_NEXTOD(I)=CNGWAT*( 0*MIEEXT(I,ILO,INDMIE) +
     $            DX*(MIEEXT(I,IHI,INDMIE) - MIEEXT(I,ILO,INDMIE)) )
               JACABSOD    =CNGWAT*( 0*MIEABS(I,ILO,INDMIE) +
     $            DX*(MIEABS(I,IHI,INDMIE) - MIEABS(I,ILO,INDMIE)) )
               JACS_NSCAOD(I)=JACS_NEXTOD(I) - JACABSOD
               JACS_G_ASYM(I)=0*MIEASY(I,ILO,INDMIE) +
     $                  DX*(MIEASY(I,IHI,INDMIE) - MIEASY(I,ILO,INDMIE))
               JACS_FINAL(I) = JACS_NEXTOD(I) - 0.5*JACS_NSCAOD(I)*(1+G_ASYM(I)) - 0.5*NSCAOD(I)*JACA_G_ASYM(I)
            ENDDO
          END IF

C      *************************
C      Intermediate particle size
       ELSE
c          print *,'YAY PARTICLE SIZE'
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
             NEXTOD(I)=CNGWAT*( MIEEXT(I,ILO,INDMIE) +
     $          X*(MIEEXT(I,IHI,INDMIE) - MIEEXT(I,ILO,INDMIE)) )
             ABSOD    =CNGWAT*( MIEABS(I,ILO,INDMIE) +
     $          X*(MIEABS(I,IHI,INDMIE) - MIEABS(I,ILO,INDMIE)) )
             NSCAOD(I)=NEXTOD(I) - ABSOD
             G_ASYM(I)=MIEASY(I,ILO,INDMIE) +
     $          X*(MIEASY(I,IHI,INDMIE) - MIEASY(I,ILO,INDMIE))
          ENDDO

          IF (DOJAC) THEN
c           from eg calrad1
C           Optical depth of cloud1 including scattering adjustment
c           K1=NEXTO1(I) - NSCAO1(I)*(1.0+G_ASY1(I))/2.0
cccc        so
c           d(K1) = d(NEXTO1(I)) - 0.5*(d(NSCAO1(I))*(1+G_ASY1(I)) + NSCAO1(I)*d(G_ASY1(I)))
            DO I=1,NCHAN
               ABSOD  = CNGWAT*( MIEABS(I,ILO,INDMIE) +
     $          X*(MIEABS(I,IHI,INDMIE) - MIEABS(I,ILO,INDMIE)) )
               ASYM = MIEASY(I,ILO,INDMIE) +
     $          X*(MIEASY(I,IHI,INDMIE) - MIEASY(I,ILO,INDMIE))

               JACA_NEXTOD(I)=( MIEEXT(I,ILO,INDMIE) +
     $            X*(MIEEXT(I,IHI,INDMIE) - MIEEXT(I,ILO,INDMIE)) )
               JACABSOD    =( MIEABS(I,ILO,INDMIE) +
     $            X*(MIEABS(I,IHI,INDMIE) - MIEABS(I,ILO,INDMIE)) )
               JACA_NSCAOD(I)=JACA_NEXTOD(I) - JACABSOD
               JACA_G_ASYM(I)=0
c               JACA_FINAL(I) = JACA_NEXTOD(I) - 0.5*JACA_NSCAOD(I)*(1+G_ASYM(I))    !!!!!  - 0.5*NSCAOD*JACA_G_ASYM(I) == 0

               !!! this is simple and correct form of derivative, except make sure you use (CNGWAT+1e-16) in denominator
               JACA_FINAL(I) = NEXTOD(I) - NSCAOD(I)*(1.0+G_ASYM(I))/2.0
               DO L = LCTOP,LCBOT
                 JACTOP_CFRCL_v(L,I) = JACA_FINAL(I) * JACTOP_CFRCL(L)
                 JACBOT_CFRCL_v(L,I) = JACA_FINAL(I) * JACBOT_CFRCL(L)
               END DO
               JACA_FINAL(I) = 1/(CNGWAT+1e-16)*JACA_FINAL(I)

               DX = 1 / ( LOG(MIEPS(IHI,INDMIE)) - LOG(MIEPS(ILO,INDMIE)) )
               DX = DX / (CPSIZE+1e-16)   !!!! since we have log(sze) ==> d/dsze = 1/sze
               JACS_NEXTOD(I)=CNGWAT*( 0*MIEEXT(I,ILO,INDMIE) +
     $            DX*(MIEEXT(I,IHI,INDMIE) - MIEEXT(I,ILO,INDMIE)) )
               JACABSOD    =CNGWAT*( 0*MIEABS(I,ILO,INDMIE) +
     $            DX*(MIEABS(I,IHI,INDMIE) - MIEABS(I,ILO,INDMIE)) )
               JACS_NSCAOD(I)=JACS_NEXTOD(I) - JACABSOD
               JACS_G_ASYM(I)=0*MIEASY(I,ILO,INDMIE) +
     $                  DX*(MIEASY(I,IHI,INDMIE) - MIEASY(I,ILO,INDMIE))
               JACS_FINAL(I) = JACS_NEXTOD(I) - 0.5*JACS_NSCAOD(I)*(1+G_ASYM(I)) - 0.5*NSCAOD(I)*JACA_G_ASYM(I)

            END DO
          END IF           
       ENDIF

       IF (DOJAC) THEN
         !! remember SARTA uncompresses ODs at angle theta, so automatically doing 1/mu dtau_atmosphere_nadir/dX
         !! so turn these cloud nadir jacs into cloud view angle jacs using masec = 1/cos(theta)
         JACA_FINAL(1:NCHAN) = MASEC * JACA_FINAL(1:NCHAN)
         JACS_FINAL(1:NCHAN) = MASEC * JACS_FINAL(1:NCHAN)
         JACTOP_CFRCL = MASEC * JACTOP_CFRCL
         JACBOT_CFRCL = MASEC * JACBOT_CFRCL
         JACTOP_CFRCL_v(:,1:NCHAN) = MASEC * JACTOP_CFRCL_v(:,1:NCHAN)
         JACBOT_CFRCL_v(:,1:NCHAN) = MASEC * JACBOT_CFRCL_v(:,1:NCHAN)
       END IF
C
c see test_cld_jacs.m
c       I = 1
c       write (*,'(A,I4,12(F12.4))') 'cld jacs',I,NEXTOD(I),NSCAOD(I),G_ASYM(I),
c     $ NEXTOD(I) - NSCAOD(I)*(1.0+G_ASYM(I))/2.0,
c     $ JACA_NEXTOD(I),JACA_NSCAOD(I),JACA_G_ASYM(I),JACA_FINAL(I),
c     $ JACS_NEXTOD(I),JACS_NSCAOD(I),JACS_G_ASYM(I),JACS_FINAL(I)

       RETURN
       END
