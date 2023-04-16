C FILE NAME: hffp_kernel_dec.f
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
C     This file contains the INPUT/OUTPUT data for SUBROUTINE KERNEL
C     which is defined in hffp_util.f
C
C=======================================================================



C     KERNEL-INPUT:
      INTEGER        IJOB,           ! job number (only for error messages)
     &               SATNUM,         ! satellite number
     &               TOTIRF,         ! total number of filter channels (IRFs)
     &               ID_IRF(MAXIRF), ! IDs of IRFs
     &               PREDID(MAXIRF,MAXMOL,MAXPRD+1),! predictor IDs   
     &               NCOEFF(MAXIRF,MAXMOL),         ! number of coefficients
     &               VERBOS,                        ! verbose
     &               JACTYP,         ! Jacobians type: 1: radiance Jacobians
                                     !                 2: T_Bright Jacobians
     &               CLCTBR,         ! calc T_bright (0=no, 1=yes)
     &               OT_RTB,         ! offset tuning of      radiances (1)
                                     !      or brightness temperatures (2)
     &               T_ADJ           ! 1: the transmittance tuning parameters
                                     !    can be changed by SUBROUTINE C_K, 
                                     !    if found to be unphysical
                                     ! 0: the transmittance tuning parameters
                                     !    will not be changed by SUBROUTINE C_K, 
                                     !    even if unphysical
      REAL*8         COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER), ! coefficients
     &               T_TAUF(MAXIRF),    ! TUNING FACTORS OF TAU(FIXED)
     &               T_TAUW(MAXIRF),    ! TUNING FACTORS OF TAU(H2O)
     &               T_TAUO(MAXIRF),    ! TUNING FACTORS OF TAU(O3)
     &               T_OFFS(MAXIRF),    ! RAD OR T_BRIGHT TUNING OFFSET
     &               T_SEC(MAXIRF),     ! OFFSET TO SECANT, CHANNEL DEPENDEND 
                                        !         (secant "tuning" parameter)
     &               SURF_T,            ! SURFACE TEMPERATURE
     &               SURF_P,            ! SURFACE PRESSURE
     &               SECANT,            ! OBSERVATION SECANT
     &               SUNSEC,            ! SUN SECANT
     &               SUREMI,            ! SURFACE EMISSIVITY
     &               SUNSOL,            ! SOLID ANGLE OF THE SUN
     &               SOLRAD(MAXIRF),    ! CONVOLVED SOLAR RADIANCES
     &               ATEMP(NLAYER),     ! ATMOSPHERE TEMPERATURES
     &               AFIXED(NLAYER),    ! ATMOSPHERE FIXED AMOUNT
     &               AWATER(NLAYER),    ! ATMOSPHERE WATER AMOUNT
     &               AOZONE(NLAYER)     ! ATMOSPHERE OZONE AMOUNT
      INTEGER        CLCJAC(3)          ! CALCULATE JACOBIANS FOR T,W,O:
      ! 0=NO,1=YES(ANALYTICAL J.),2=YES(FINITE DIFFERENCES), 3=YES(both)
      ! NOTE: finite differences Jacobians are NOT calculated by this
      !       routine. This routine calculates the analytical Jacobians
      !       if CLCJAC is 1 or 3.

C     KERNEL-OUTPUT
      REAL*8         CW_M(MAXIRF),      ! CHANNEL CENTER WAVE NUMBERS as used 
     &               RADTOT(MAXIRF),    ! RADIANCES TOTAL (UPWELLING+REFLECED)
     &               TBMEAN(MAXIRF),    ! BRIGTH. TEMP FOR MEAN CENTER WAVENO.
     &               JACOBT(MAXIRF,NLAYER),        ! JACOBIANS dY/dT
     &               JACOBW(MAXIRF,NLAYER),        ! JACOBIANS dY/dW
     &               JACOBO(MAXIRF,NLAYER)         ! JACOBIANS dY/dO
      INTEGER        NOTECT             ! NOTE IF T_TAUX CHANGED IN C_TAU
