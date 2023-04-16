C FILE NAME: hffp_aux_dec.f
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

C     ===============================================================
C     EDIT PATH NAMES
C     ATTENTION: don't forget the / AND a blank after the path name !
C     ===============================================================

C     PATH FOR CENTER WAVENUMBERS FILES AND SOLAR RADIANCE FILES
      CHARACTER*80 CFFPTH
C     PARAMETER (CFFPTH='/salsify/users/wehr/HIRS/filter/ ')    !salsify
      PARAMETER (CFFPTH='./INPUT/ ')               !squirrel,muir,salsify
C     PARAMETER (CFFPTH='/home/wehr/HIRS/JCET/INPUT/ ')         !molotov

      CHARACTER*80 PCPATH    ! PATH NAME FOR PREDICTOR AND COEFF. PATH
C     PARAMETER(PCPATH='/salsify/users/wehr/HIRS/fmodel/ ')     !salsify
      PARAMETER(PCPATH='./INPUT/ ')                !squirrel,muir,salsify
c     PARAMETER(PCPATH='/home/wehr/HIRS/JCET/INPUT/ ')          !molotov

      ! The following line is only used by the interface program, but
      ! not by KERNEL
      CHARACTER*80 OUTPTH ! OUTPUT PATH NAME
c     PARAMETER(OUTPTH='/salsify/data/Wehr/hirs/fast_output/ ') !salsify
      PARAMETER(OUTPTH='./OUTPUT/ ')               !squirrel,muir,salsify

C     =============================================================
C     DO NOT CHANGE PARAMETERS BELOW
C     =============================================================

C     NAMES OF CENTER WAVENUMBERS FILES
      ! NOTE: these files contain wavenumbers in units [cm-1]
      CHARACTER*80 CFFILE(7)
      DATA CFFILE/'noaa7.centerfrq','noaa9.centerfrq',
     &            'noaa10.centerfrq','noaa11.centerfrq',
     &            'noaa12.centerfrq','noaa14.centerfrq',
     &            'noaa15.centerfrq'/
      
C     NAMES OF PREDICTOR ID FILES
      CHARACTER*11 PIDFLS(7)
      DATA PIDFLS/'PRED_7.dat','PRED_9.dat','PRED_10.dat',
     &            'PRED_11.dat','PRED_12.dat','PRED_14.dat',
     &            'PRED_15.dat'/

C     NAMES OF COEFFICIENT FILES
      CHARACTER*11 COFFLS(7)
      DATA COFFLS/'COEF_7.dat','COEF_9.dat','COEF_10.dat',
     &            'COEF_11.dat','COEF_12.dat','COEF_14.dat',
     &            'COEF_15.dat'/

C     NAMES OF SOLAR RADIANCES FILES
      CHARACTER*13  SOLRDF(7)
      DATA SOLRDF/'SOLRAD_7.dat','SOLRAD_9.dat','SOLRAD_10.dat',
     &            'SOLRAD_11.dat','SOLRAD_12.dat','SOLRAD_14.dat',
     &            'SOLRAD_15.dat'/

C     CORRESPONDING NUMBERS:
C     POSSIBLE NUMBERS IN CONTROL FILE, THE INDECEES OF NOAANR
C     CORRESPOND TO THE INDECEES OF PIDFLS AND COFFLS
      INTEGER NOAANR(7)
      DATA NOAANR/7,9,10,11,12,14,15/

