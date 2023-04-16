C FILE NAME: hffp_init_dec.f
C
C=======================================================================
C
C    Joint Center for Earth Systems Technology (JCET)
C    University of Maryland Baltimore County   (UMBC)
C
C    HIRS Fast Forward Program
C
C    Version number see HFFPVS in hffp_glob_dec.f
C
C                                   Tobias Wehr       wehr@umbc.edu
C                                   L.Larrabee Strow  strow@umbc.edu
C
C=======================================================================
C

C     PARAMETER DECLARATIONS

C     INITIALIZATION PARAMETERS
      CHARACTER*80   CTLFIL,            ! NAME OF CONTROL FILE
     &               CTLDMY,            ! INITIAL CONTROL FILE NAME DUMMY
     &               ATMFLN(MAXJOB)     ! NAME OF ATMOSPHERE FILE
      INTEGER        SATNUM,            ! SATELLITE NUMBER (7,9,10,11,12,...)
     &               NOJOBS,            ! NUMBER OF ATMOSPHERES TO BE PROCESSED
     &               TOTIRF,            ! TOTAL NUMBER OF IRFs
     &               ID_IRF(MAXIRF),    ! IDs OF IRFs 
     &               PREDID(MAXIRF,MAXMOL,MAXPRD+1),! predictor identifiers (IDs)
                                        ! for every channel (if used or not!) 
                                        ! and every molecule (F,W,O):
                                        ! 3rd-index vector is:
                                        !   number of predictors, [predictor IDs]
     &               NCOEFF(MAXIRF,MAXMOL) ! number of coefficients
      REAL*8         COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER), ! coefficients
     &               T_TAUF(MAXIRF),  ! TUNING FACTORS OF TAU(FIXED)
     &               T_TAUW(MAXIRF),  ! TUNING FACTORS OF TAU(H2O)
     &               T_TAUO(MAXIRF),  ! TUNING FACTORS OF TAU(O3)
     &               T_OFFS(MAXIRF),  ! OFFSET TUNING OF EITHER THE RADIANCES
                                      !   OR THE BRIGHTNESS TEMPERATURES
     &               T_SEC(MAXIRF),   ! OFFSET TO ZENITH SECANT, CHANNEL DEPENDEND
     &               SURF_T(MAXJOB),  ! SURFACE TEMPERATURES
     &               SURF_P(MAXJOB),  ! SURFACE PRESSURES
     &               SECANT(MAXJOB),  ! OBSERVATION SECANT
     &               SUNSEC(MAXJOB),  ! SUN SECANT
     &               SUREMI(MAXJOB),  ! SURFACE EMISSIVITY
     &               SUNSOL(MAXJOB),  ! SOLID ANGLE OF THE SUN
     &               SOLRAD(MAXIRF),  ! CONVOLVED SOLAR RADIANCES
     &               ATEMP(NLAYER),   ! ATMOSPHERE TEMPERATURES
     &               AFIXED(NLAYER),  ! ATMOSPHERE FIXED AMOUNT
     &               AWATER(NLAYER),  ! ATMOSPHERE WATER AMOUNT
     &               AOZONE(NLAYER)   ! ATMOSPHERE OZONE AMOUNT
      INTEGER        WRTOUT           ! WRITE OUTPUT TO FILE (1=YES)

      DATA CTLDMY /'sample.ctl'/
      





