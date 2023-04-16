C FILE NAME: hffp_glob_dec.f
C
C=======================================================================
C
C    Joint Center for Earth Systems Technology (JCET)
C    University of Maryland Baltimore County   (UMBC)
C
C    HIRS Fast Forward Program
C
C    Version number see HFFPVS (in this file)
C
C                                   Tobias Wehr       wehr@umbc.edu
C                                   L.Larrabee Strow  strow@umbc.edu
C
C=======================================================================

C     ===========================================================
C     VERSION INFORMATION
C     ===========================================================
      CHARACTER*40 HFFPVS
      DATA HFFPVS /'HFFP VERSION 1.04 (UMBC, Sep. 1998)'/
      ! 1.04 includes NOAA-15

C     ===========================================================
C     YOU MIGHT WANT TO EDIT THIS
C     ===========================================================

      INTEGER MAXJOB         ! MAXIMAL NUMBER OF ATMOSPHERE PROFILE FILE NAMES
      PARAMETER(MAXJOB=50)   ! YOU CAN SET THIS TO 1 FOR DAO RETRIEVAL SYSTEM

C     ===========================================================
C     DO NOT EDIT BELOW
C     ===========================================================

      INTEGER MATLAB         ! VERBOSE PARAMETER FOR SOME OUTPUT TO SCREEN
      PARAMETER(MATLAB=0)    ! IN MATLAB FORMAT (1=ON, 0=OFF)

      INTEGER NLAYER         ! NUMBER OF ATMOSPHERIC LAYER
      PARAMETER(NLAYER=100)

      INTEGER MAXMOL         ! MAXIMAL NUMBER OF MOLECULES
      PARAMETER(MAXMOL=3)    ! 3 IN THE ORDER OF: FIXED, WATER, OZONE
                             ! This is also used for the Jacobian-variables,
                             ! where the order of the index would be de-
                             ! rivative with respect to temperatures, water,
                             ! ozone. (I.e., temperature replaces fixed.)

      INTEGER MAX_QP         ! MAXIMAL NUMBER OF Q-PROFILES
      PARAMETER(MAX_QP=11)

      INTEGER MAXPRD         ! MAXIMAL NUMBER OF PREDICTORS PER IRF PER MOL
      PARAMETER(MAXPRD=33)

      INTEGER MAXIRF         ! MAXIMAL NUMBER OF IRFs (FILTER CHANNELS)
      PARAMETER(MAXIRF=19)

      INTEGER SNIIRF         ! INDEX OF FIRST SOLAR CHANNEL: FOR ALL 
      PARAMETER(SNIIRF=13)   ! CHANNELS WITH CHANNEL ID .GE. SNIIRF, THE
                             ! REFLECTED SOLAR MIGHT BE CALCULATED

      REAL*8 MAXMSV          ! MAXIMAL MODELED SECANT VALUE in the 
      PARAMETER(MAXMSV=8.95) ! coefficient fits

      REAL*8 MAXSNS
      PARAMETER(MAXSNS=100.0)! MAXIMAL SUN SECANT VALUE, when sun secant
                             ! is larger than MAXSNS or negative, it will
                             ! be set to MAXSNS. The value of MAXSNS is
                             ! quite arbitrary. It must be "large" compared
                             ! to MAXMSV

      INTEGER MILAYT         ! MAXIMUM NUMBER OF LAYER-ABOVE CONSIDERATION
      PARAMETER(MILAYT=30)   ! for calculating the analytic T-Jacobians
                             ! (NLAYER is the maximum useful and best value)
                             ! proposed: 30

      INTEGER MILAYW         ! MAXIMUM NUMBER OF LAYER-ABOVE CONSIDERATION
      PARAMETER(MILAYW=80)   ! for calculating the analytic W-Jacobians
                             ! (NLAYER is the maximum useful and best value)
                             ! proposed: 80

      INTEGER MILAYO         ! MAXIMUM NUMBER OF LAYER-ABOVE CONSIDERATION
      PARAMETER(MILAYO=70)   ! for calculating the analytic O-Jacobians
                             ! (NLAYER is the maximum useful and best value)
                             ! proposed: 70

