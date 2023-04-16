C FILE NAME: hffp_tau.f
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
C ====================================================================
C     Routines to calculate transmittances
C ====================================================================

      SUBROUTINE C_TAU(IJOB,SATNUM,TOTIRF,ID_IRF,PREDID,
     &                  NCOEFF,COEFF,SECANT,SUNSEC,
     &                  T_TAUF,T_TAUW,T_TAUO,T_SEC,TAU,TAU_RS,
     &                  ATEMP,AFIXED,AWATER,AOZONE,
     &                  RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEG,K_NEGS,
     &                  VERBOS,T_ADJ,NOTECT)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     calculates the total layer transmittances TAU for the upwelling
C     thermal radiation and TAU_RS for the reflected solar
C
C     INPUT:
C       IJOB           job number (= atmosphere number) from the main
C                      program loop. This parameter is only used in
C                      an warning message for the reflected solar.
C       SATNUM         satellite number (e.g. 7 for NOAA7)
C       TOTIRF         total number of IRFs (channels)
C       ID_IRF         IDs of IRFs
C       PREDID         predictor identifiers
C       NCOEFF         number of coefficients
C       COEFF          coefficients
C       SECANT         secant (of observation)
C       SUNSEC         secant of sun
C       ATEMP          atmosphere temperature profile
C       AFIXED         atmosphere fixed gases profile
C       AWATER         atmosphere water profile
C       AOZONE         atmosphere ozone profile
C       RTEMP          reference atmosphere temperature profile
C       RFIXED         reference atmosphere fixed gases profile
C       RWATER         reference atmosphere water profile
C       ROZONE         reference atmosphere ozone profile
C       PRES           pressure profile
C       VERBOS         screen output on (1) or off (0)
C       T_ADJ           1: the transmittance tuning parameters can be changed
C                          by SUBROUTINE C_K, if found to be unphysical
C                       0: the transmittance tuning parameters will not be 
C                          changed by SUBROUTINE C_K, even if unphysical
C       T_SEC          channel dependend zenith secant offset
C
C     INPUT/OUTPUT
C       T_TAUF         tuning parameters for tau(fixed)
C       T_TAUW         tuning parameters for tau(water)
C       T_TAUO         tuning parameters for tau(ozone)
C       NOTECT         1 if T_TAUX has been changed, 0 if not
C
C     OUTPUT:
C       TAU            transmittances for upwelling thermal
C       TAU_RS         transmittances for reflected solar
C       K_NEG          indicator if the absorption coefficient K is 
C                      negative for IIRF,IMOL,ILAY. This is an information
C                      that is needed for the Jacobian calculations.
C                      0: K is non-negative, 1: K is negative and disregarded
C       K_NEGS         (same for solar)
C 
C     CALLING ROUTINE:
C       KERNEL
C
C     COMMENTS:
C       for k > MAX_K (locally defined constant) the corresponding
C       TAU and TAU_RS will be set to zero
C
C     HISTORY:
C       written 10/17/97 by Tobias Wehr
C       solar   11/13/97 by Tobias Wehr
C --------------------------------------------------------------------
      IMPLICIT NONE
      include "hffp_glob_dec.f"

C     SUBROUTINE ARGUMENT PARAMETER
      INTEGER        SATNUM,TOTIRF,ID_IRF(MAXIRF),IJOB,T_ADJ
      REAL*8         T_TAUF(MAXIRF),T_TAUW(MAXIRF),T_TAUO(MAXIRF),
     &               T_SEC(MAXIRF)
      REAL*8         TAU(MAXIRF,NLAYER),     
     &               TAU_RS(MAXIRF,NLAYER),  
     &               SECANT,SUNSEC
      INTEGER        K_NEG(MAXIRF,MAXMOL,NLAYER),
     &               K_NEGS(MAXIRF,MAXMOL,NLAYER)
C     atmospheric profile (for which tau has to be calculated)
      REAL*8         ATEMP(NLAYER),AFIXED(NLAYER),
     &               AWATER(NLAYER),AOZONE(NLAYER)
C     reference profile (of fast algorithm development)
      REAL*8         RTEMP(NLAYER),RFIXED(NLAYER),
     &               RWATER(NLAYER),ROZONE(NLAYER),
     &               PRES(NLAYER)
      INTEGER        PREDID(MAXIRF,MAXMOL,MAXPRD+1) ! predictor IDs
                     ! for every channel (if used or not!) and every molecule (F,W,O):
                     ! 3rd-index vector is:
                     !   [number of predictors, [predictor IDs]]

      INTEGER        NCOEFF(MAXIRF,MAXMOL)
      REAL*8         COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER)
      INTEGER        VERBOS
C     marker, if tuning paramters T_TAUX have been changed in C_K
      INTEGER        NOTECT ! 0=no, 1=yes

C     LOCAL PARAMETER

      ! Q-profiles are profiles used for the predictor calculations. 
      ! They are composed from temperature, pressure, water- and 
      ! ozone-profiles from both the actual and the reference 
      ! atmosphere profile. The Q-profiles contains layer-above-
      ! information.
      REAL*8         Q_PROF(MAX_QP,NLAYER)
      ! the derivatives of the Q-profiles with respect to
      ! temperature, water, ozone. The first index of DQDX is
      ! 1 for derivative after T, 2 for W, 3 for O

C     QPINIT(j)=1 if Q-profile has been already initialized
      ! this variable shall avoid unnecessary calculation of a 
      ! Q-profile. A Q-profile will only be calculated if is 
      ! will be used. If it has already being used before, this 
      ! flag is set to 1 and it will not be calculated anymore.
      INTEGER        QPINIT(MAX_QP) 

C     k-profiles
      REAL*8         K_PROF(NLAYER),   ! thermal upwelling abs.coeff.
     &               K_RSPF(NLAYER)    ! reflected solar abs.coeff.
      INTEGER        RSWARN            ! warning already given (=1)
      INTEGER        IIRF,IIRF_A,ILAY,IMOL
      INTEGER        IHELP
C     effective secant (sun + satellite) for refl. solar
      REAL*8         EFFSEC
c     maximum value for abs.coeff. k, for k>MAX_K we set tau=0
      REAL*8         MAX_K
      PARAMETER(MAX_K=30.0)
C     SCOUNT counts how many times "BLOCK K" has been executed:
      INTEGER        SCOUNT
C     marker, if tuning paramters T_TAUX have been changed in C_K
      INTEGER        CHTTAU ! 0: no, >0: yes
C     corrected secant (SECANT + T_SEC) for a respective channel
      REAL*8         C_SEC

      DATA           SCOUNT/0/
      SAVE           SCOUNT

C     ============================================================
C     !                        BEGIN                             !
C     ============================================================

C     INIT RSWARN
      ! if solar reflected will not be calculated due to too large
      ! effective angle, this parameter will be set to 1 and refl.
      ! solar will not be considered for this job anymore
      RSWARN=0
C     INITE NOTECT
      NOTECT=0

C     RESET QPINIT TO ZERO
      DO 10 IHELP=1,MAX_QP ! do not use a "DATA" statement for this, 
         QPINIT(IHELP)=0   ! since it might not reset QPINIT every 
 10   CONTINUE             ! time when C_TAU is called (hffp_main
                           ! might be called multiple times!!!)

C     check consistency of PREDID and NCOEFF
      IF (IJOB .EQ. 1) THEN
         DO 100 IIRF=1,MAXIRF
            DO 200 IMOL=1,MAXMOL
               IF (PREDID(IIRF,IMOL,1) .NE. NCOEFF(IIRF,IMOL)) THEN
 201              FORMAT('FATAL ERROR IN C_TAU: ','PREDID(',I2,',',
     &                   I1,',1) .NE. NCOEFF(',I2,',',I1,')')
 202              FORMAT('PREDID(IIRF,IMOL,1)= ',I3,', ',
     &                 'NCOEFF(IIRF,IMOL)=   ',I3)
 203              FORMAT('IIRF= ',I2,', IMOL= ',I2)
 204              FORMAT('DETAILS:')
                  WRITE (*,201) IIRF,IMOL,IIRF,IMOL
                  WRITE (*,204)
                  WRITE (*,203) IIRF,IMOL
                  WRITE (*,202) PREDID(IIRF,IMOL,1),NCOEFF(IIRF,IMOL)
                  WRITE (*,*) '>>> PROGRAM ABORTED <<<'
                  STOP
               ENDIF
 200        CONTINUE
 100     CONTINUE
      ENDIF

 500  CONTINUE ! THIS IS A GOTO-LOOP DESTINATION
      CHTTAU=0
C     ---------------------------------------------------
C     ---------- BEGIN "BLOCK K" ------------------------
C     ---------------------------------------------------
C     loop desired channels (IRFs)
      DO 1000 IIRF=1,TOTIRF
         ! reset K_PROF for each channel
         Do 1001 ILAY=1,NLAYER
            K_PROF(ILAY)=0.0
            K_RSPF(ILAY)=0.0
 1001    CONTINUE
         IIRF_A=ID_IRF(IIRF)    ! IIRF_A contains now the channel ID number

         ! calculate corrected ("tuned") secant for the respective channel
         CALL ADDSEC(SECANT,T_SEC(IIRF_A),C_SEC)

         CALL C_K(IIRF_A,PREDID,NCOEFF,COEFF,
     &        QPINIT,Q_PROF,K_PROF,C_SEC,
     &        ATEMP,AFIXED,AWATER,AOZONE,
     &        RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &        T_TAUF,T_TAUW,T_TAUO,K_NEG,VERBOS,
     &        SCOUNT,T_ADJ,CHTTAU)

c        write (*,*) 'channel ',iirf
c        do 9999 ilay=1,nlayer
c           write (*,*) ilay, k_prof(ilay)
c9999    continue

C        REFLECTED SOLAR
         IF ((IIRF_A .GE. SNIIRF) .AND. (RSWARN .EQ. 0)) THEN
            ! calculate EFFSEC first as the effective secant for sun
            CALL ADDSEC(C_SEC,SUNSEC,EFFSEC)

            IF (EFFSEC .GT. MAXMSV) THEN
               RSWARN=1
               IF ((VERBOS .EQ. 1) .OR. (MATLAB .EQ. 1)) THEN
                  WRITE (*,*) 'WARNING: SUN ANGLE LARGER THAN MODELED'
                  WRITE (*,*) '  refleced solar will not be calculated'
                  WRITE (*,1101) IJOB
 1101             FORMAT('   for atmosphere no. ',I4)
               ENDIF
            ELSE
               ! CALCULATE ABSORPTION COEFFICIENTS FOR REFLECTED SOLAR
               CALL C_K(IIRF_A,PREDID,NCOEFF,COEFF,
     &              QPINIT,Q_PROF,K_RSPF,EFFSEC,
     &              ATEMP,AFIXED,AWATER,AOZONE,
     &              RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &              T_TAUF,T_TAUW,T_TAUO,K_NEGS,VERBOS,
     &              SCOUNT,T_ADJ,CHTTAU)
            ENDIF
         ENDIF

C        ! ===============================================
C        WRITE (*,*) 'AbsCoeffF(',IIRF_A,',:)=[ ...'
C        DO 1999 IHELP=1,NLAYER-1
C        WRITE (*,1998) K_PROF(IHELP)
C 1998   FORMAT(E,',...')
C 1999   CONTINUE
C        WRITE (*,*) K_PROF(NLAYER),'];'
C        ! ===============================================

         ! calculate tau; set tau=0 for k>MAX_K
         DO 2001 ILAY=1,NLAYER
            IF (K_PROF(ILAY) .GT. MAX_K) THEN
               TAU(IIRF_A,ILAY)=0.0
            ELSE
               TAU(IIRF_A,ILAY)=DEXP(-1*K_PROF(ILAY))
            ENDIF
 2001    CONTINUE
         !WRITE (*,*) TAU(IIRF_A,4) ?????
         ! calculate tau for reflected solar
         ! set tau-solar to zero for non-solar channels
         DO 2002 ILAY=1,NLAYER
            IF (IIRF_A .GE. SNIIRF) THEN
               IF (K_RSPF(ILAY) .GT. MAX_K) THEN
                  TAU_RS(IIRF_A,ILAY)=0.0
               ELSE
                  TAU_RS(IIRF_A,ILAY)=DEXP(-1*K_RSPF(ILAY))
               ENDIF
            ELSE
               TAU_RS(IIRF_A,ILAY)=0.0
            ENDIF
 2002    CONTINUE
         
 1000 CONTINUE
C     ---------------------------------------------------
C     ---------- END "BLOCK K" --------------------------
C     ---------------------------------------------------
      SCOUNT=SCOUNT+1

      IF ((CHTTAU .GT. 0) .AND. (SCOUNT .EQ. 1)) THEN
         IF (VERBOS .EQ. 1) THEN
            WRITE (*,*) 'NOTE: CALCULATION OF TRANSMITTANCES WILL BE'
            WRITE (*,*) '      REPEATED BECAUSE OF AUTOMATICALLY'
            WRITE (*,*) '      CHANGED TUNING PARAMETERS'
         ENDIF
         NOTECT=1
         GOTO 500
      ENDIF

 
      RETURN
      END

C ====================================================================

      SUBROUTINE C_K(IIRF,PREDID,NCOEFF,COEFF,
     &                    QPINIT,Q_PROF,K_PROF,SECANT,
     &                    ATEMP,AFIXED,AWATER,AOZONE,
     &                    RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &                    T_TAUF,T_TAUW,T_TAUO,K_NEG,VERBOS,
     &                    SCOUNT,T_ADJ,CHTTAU)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     Calculate tuned absoption coefficients for IIRF.
C     This routine might also change the transmittance tuning paramters,
C     so that the tuned absorption coefficients of a certain gas (F,W,O)
C     cannot be smaller than zero (in every layer).
C
C     INPUT:
C       IIRF            IRF counter = actual IRF (channel, or channel ID)
C       PREDID          IDs of the predictors to be used (see COMMENT)
C       NCOEFF          number of coefficients to be used (see COMMENT)
C       COEFF           coefficients to be used (see COMMENT)
C       SECANT          observation secant
C       ATEMP,AFIXED,AWATER,AOZONE       atmosphere profile          
C       RTEMP,RFIXED,RWATER,ROZONE,PRES  reference atmosphere profile
C       VERBOS          screen output on/off (1=on, 0=off)
C       SCOUNT          counter how many times C_TAU has been called yet
C       T_ADJ           1: the transmittance tuning parameters can be
C                          changed by this routine, if found to be unphysical
C                       0: the transmittance tuning parameters will not be
C                          changed by this routine, even if unphysical
C
C     INPUT AND OUTPUT:
C       QPINIT          Q-profiles initialization flags (see COMMENT)
C       Q_PROF          Q-profiles (see COMMENT)
C       T_TAUF,         tuning parameters might be changed by this 
C       T_TAUW,         routine so that the tuned absorption coefficient
C       T_TAUO          of a certain gas (F,W,O) is never negative
C
C     OUTPUT:
C       K_PROF          k-profile
C       K_NEG           indicator if the absorption coefficient K is 
C                       negative for IIRF,IMOL,ILAY. This is an information
C                       that is needed for the Jacobian calculations.
C                       0: K is non-negative, 1: K is negative and disregarded
C       CHTTAU          INCREMENTED BY 0: T_TAUX has not been changed by C_K
C                       INCREMENTED BY 1: T_TAUX has been changed by C_K
C
C     CALLING ROUNTINE
C       C_TAU
C
C     COMMENTS:
C       1. PREDID,NCOEFF,COEFF contain the information for ALL channels,
C          no matter whether all channels are used or not.
C       2. The Q-profiles are the profiles which are needed to calculate
C          Q (for example the "layer-above-profiles"); Q-profiles are NOT
C          the Q-matrix!
C       3. there is no need to store the calculated predictors Q explicitly
C          in a variable; instead k will be calculated immediately as k=Qc 
C          with Q=predictor, c=coefficients
C       4. the k-profile will be calculated successively for each molecule
C          (K_MOL) and added afterwards to the total k-profile (K_PROF).
C          K_MOL(layer) will only be added to K_PROF(layer) if > 0.
C
C     HISTORY:
C       written 9/15/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE CALLING VARIABLES
      INTEGER     IIRF,
     &            NCOEFF(MAXIRF,MAXMOL),QPINIT(MAX_QP),
     &            PREDID(MAXIRF,MAXMOL,MAXPRD+1),
     &            K_NEG(MAXIRF,MAXMOL,NLAYER),VERBOS,
     &            SCOUNT,CHTTAU,T_ADJ
      REAL*8      COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER),
     &            Q_PROF(MAX_QP,NLAYER),K_PROF(NLAYER),SECANT
      REAL*8      ATEMP(NLAYER),AFIXED(NLAYER),AWATER(NLAYER),
     &            AOZONE(NLAYER),RTEMP(NLAYER),RFIXED(NLAYER),
     &            RWATER(NLAYER),ROZONE(NLAYER),PRES(NLAYER)
      REAL*8      T_TAUF(MAXIRF),T_TAUW(MAXIRF),T_TAUO(MAXIRF) ! tuning parameters

C     LOCAL VARIABLES
      INTEGER     IMOL,           ! molecule counter
     &            ILAY,           ! layer counter
     &            IPRED,          ! predictor counter
     &            PID             ! predictor ID
      REAL*8      KADD(NLAYER),   ! K_PROF increment. NOTE: KADD could 
                                  ! as well be a scalar rather than a
                                  ! vector. It has been choosen a vector only
                                  ! for debugging/monitoring purposes!
     &            K_MOL(NLAYER),  ! K-profile of a certain molecule
     &            TUNING,         ! tuning parameter
     &            ADDTOK          ! = k(molecule) + tuning parameter

C     VARIABLES ONLY USED FOR DEBUGGING
!     INTEGER     I

C     ===========================================================
C     BEGIN
C     ===========================================================

C     ........ CALCULATE PREDICTORS: ...........
C     initialize k-profile with zeros
      DO 10 ILAY=1,NLAYER
         K_PROF(ILAY)=0.0
 10   CONTINUE

C     loop molecules
      DO 700 IMOL=1,MAXMOL
         ! reset K_MOL
         DO 11 ILAY=1,NLAYER
            K_MOL(ILAY)=0.0
 11      CONTINUE
         
C        loop all required predictors
         DO 600 IPRED=1,PREDID(IIRF,IMOL,1) ! PREDID(IIRF,IMOL,1) is the total
                                            ! number of predictors for (IIRF,IMOL)
            ! PID is introduced for programmer's convenience:
            PID = PREDID(IIRF,IMOL,1+IPRED) 

c           !==========================================
c           IF ((IMOL .EQ. 2) .OR. (IMOL .EQ. 3)) THEN
c              WRITE (*,678) PID,IMOL
c           ENDIF
c678        FORMAT(' C_K:    PID,IMOL=',I3,I3) 
c           !==========================================

C ................................................................
C        PREDICTORS FOR FIXED
C ................................................................
C        FIXED 1
C ................................................................
            IF ((IMOL .EQ. 1) .AND. (PID .EQ. 1)) THEN
               ! initialize Q-profile (the subroutine I_Q will
               ! check if already initialized or not)
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 101 ILAY=1,NLAYER
                  ! sum up k
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT * DSQRT(Q_PROF(1,ILAY))
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 101           CONTINUE
C ................................................................
C        FIXED 2
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 2)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 102 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)/Q_PROF(1,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 102           CONTINUE
C ................................................................
C        FIXED 3
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 3)) THEN
               DO 103 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*SECANT
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 103           CONTINUE
C ................................................................
C        FIXED 4
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 4)) THEN
               DO 104 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*(SECANT**2)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 104           CONTINUE
C ................................................................
C        FIXED 5
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 5)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 105 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(1,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 105           CONTINUE
C ................................................................
C        FIXED 6
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 6)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(2,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 106 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(2,ILAY)/Q_PROF(1,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 106           CONTINUE
C ................................................................
C        FIXED 7
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 7)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 107 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(1,ILAY)**3
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 107           CONTINUE
C ................................................................
C        FIXED 8
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 8)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 108 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                       Q_PROF(1,ILAY)**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 108           CONTINUE
C ................................................................
C        FIXED 9
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 9)) THEN
               CALL I_Q(2,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 109 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                       SECANT*Q_PROF(2,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 109           CONTINUE
C ................................................................
C        FIXED 10
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 10)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 110 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(1,ILAY)**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 110           CONTINUE
C ................................................................
C        FIXED 11
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 11)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 111 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                       SECANT*Q_PROF(1,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 111           CONTINUE
C ................................................................
C        FIXED 12
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 12)) THEN
               DO 112 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*SECANT**(1.5)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 112           CONTINUE
C ................................................................
C        FIXED 13
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 13)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 113 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 Q_PROF(1,ILAY)*SECANT**(1.5) 
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 113           CONTINUE
C ................................................................
C        FIXED 14
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 14)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 114 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (Q_PROF(1,ILAY)*SECANT)**(1.5)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 114           CONTINUE
C ................................................................
C     FIXED 15
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 15)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 115 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 Q_PROF(1,ILAY)*SECANT**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 115           CONTINUE
C ................................................................
C        FIXED 16
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .EQ. 16)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 116 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (Q_PROF(1,ILAY)*SECANT)**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 116           CONTINUE
C ................................................................
C        FIXED - UNDEFINED
C ................................................................
            ELSE IF ((IMOL .EQ. 1) .AND. (PID .GT. 16)) THEN
               WRITE (*,*) 'FATAL ERROR IN C_K:'
               WRITE (*,*) '((IMOL .EQ. 1) .AND. (PID .GT. 16))'
               WRITE (*,*) 'PROGRAM ABORTED'
               STOP
C ................................................................
C        PREDICTORS FOR WATER
C ................................................................
C        WATER 1
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 1)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 201 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 Q_PROF(4,ILAY)*SECANT
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 201           CONTINUE
C ................................................................
C        WATER 2
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 2)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 202 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DSQRT(Q_PROF(4,ILAY)*SECANT)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 202           CONTINUE
C ................................................................
C        WATER 3
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 3)) THEN
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 203 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DSQRT(Q_PROF(5,ILAY)*SECANT)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 203           CONTINUE
C ................................................................
C     WATER 4
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 4)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 204 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (Q_PROF(4,ILAY)*SECANT)**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 204           CONTINUE
C ................................................................
C        WATER 5
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 5)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(6,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 205 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 Q_PROF(6,ILAY)*Q_PROF(4,ILAY)*SECANT
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 205           CONTINUE
C ................................................................
C        WATER 6
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 6)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 206 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (Q_PROF(4,ILAY)*SECANT)**3
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 206           CONTINUE
C ................................................................
C        WATER 7
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 7)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 207 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*(Q_PROF(4,ILAY)**2)/Q_PROF(5,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 207           CONTINUE
C ................................................................
C        WATER 8
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 8)) THEN
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 208 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(5,ILAY))**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 208           CONTINUE
C ................................................................
C        WATER 9
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 9)) THEN
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 209 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(4,ILAY)*Q_PROF(3,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 209           CONTINUE
C ................................................................
C        WATER 10
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 10)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 210 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(4,ILAY))**0.25
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 210           CONTINUE
C ................................................................
C        WATER 11
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 11)) THEN
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 211 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(5,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 211           CONTINUE
C ................................................................
C        WATER 12
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 12)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(6,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 212 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(4,ILAY)*(Q_PROF(6,ILAY)**2)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 212           CONTINUE
C ................................................................
C        WATER 13
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 13)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 213 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 Q_PROF(4,ILAY)/Q_PROF(1,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 213           CONTINUE
C ................................................................
C        WATER 14
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 14)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 214 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)/
     &                 Q_PROF(1,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 214           CONTINUE
C ................................................................
C        WATER 15
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 15)) THEN
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 215 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(3,ILAY)*
     &                 DSQRT(SECANT*Q_PROF(4,ILAY))
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 215           CONTINUE
C ................................................................
C        WATER 16
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 16)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 216 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*SECANT*
     &                 (Q_PROF(4,ILAY)**1.5)/Q_PROF(5,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 216           CONTINUE
C ................................................................
C        WATER 17
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 17)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 217 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*(SECANT**2)*
     &                 (Q_PROF(4,ILAY)**3)/Q_PROF(5,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 217           CONTINUE
C ................................................................
C        WATER 18
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 18)) THEN
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 218 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*(Q_PROF(4,ILAY)**2)*
     &                 Q_PROF(3,ILAY)/Q_PROF(5,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 218           CONTINUE
C ................................................................
C        WATER 19
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 19)) THEN
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 219 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(5,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 219           CONTINUE
C ................................................................
C        WATER 20
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 20)) THEN
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 220 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(3,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 220           CONTINUE
C ................................................................
C        WATER 21
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 21)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 221 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(1,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 221           CONTINUE
C ................................................................
C        WATER 22
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 22)) THEN
               DO 222 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*SECANT
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 222           CONTINUE
C ................................................................
C        WATER 23
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 23)) THEN
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 223 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(4,ILAY)*
     &                 Q_PROF(3,ILAY)*ABS(Q_PROF(3,ILAY))
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 223           CONTINUE
C ................................................................
C        WATER 24
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 24)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 224 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(4,ILAY))**1.25
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 224           CONTINUE
C ................................................................
C        WATER 25
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 25)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 225 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(4,ILAY))**1.5
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 225           CONTINUE
C ................................................................
C        WATER 26
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 26)) THEN
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 226 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(3,ILAY)*
     &                 (SECANT*Q_PROF(4,ILAY))**1.5
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 226           CONTINUE
C ................................................................
C        WATER 27
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 27)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 227 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(4,ILAY))**2.5
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 227           CONTINUE
C ................................................................
C        WATER 28
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 28)) THEN
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 228 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(3,ILAY)*
     &                 (SECANT*Q_PROF(4,ILAY))**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 228           CONTINUE
C ................................................................
C        WATER 29
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 29)) THEN
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 229 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(4,ILAY)/Q_PROF(5,ILAY))**1.5
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 229           CONTINUE
C ................................................................
C        WATER 30
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 30)) THEN
               CALL I_Q(5,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 230 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(4,ILAY)/Q_PROF(5,ILAY))**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 230           CONTINUE
C ................................................................
C        WATER 31
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 31)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 231 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DLOG(1.1*SECANT)*Q_PROF(4,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 231           CONTINUE
C ................................................................
C        WATER 32
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 32)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 232 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DLOG(1.1*SECANT*Q_PROF(4,ILAY))
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 232           CONTINUE
C ................................................................
C        WATER 33
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .EQ. 33)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 233 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DEXP(SECANT/2.0)*Q_PROF(4,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 233           CONTINUE
C ................................................................
C        WATER - UNDEFINED
C ................................................................
            ELSE IF ((IMOL .EQ. 2) .AND. (PID .GT. 33)) THEN
               WRITE (*,*) 'FATAL ERROR IN C_K:'
               WRITE (*,*) '((IMOL .EQ. 2) .AND. (PID .GT. 33))'
               WRITE (*,*) 'PROGRAM ABORTED'
               STOP
C ................................................................
C        PREDICTORS FOR OZONE
C ................................................................
C        OZONE 1
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 1)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 301 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DSQRT(SECANT*Q_PROF(7,ILAY))
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 301           CONTINUE
C ................................................................
C        OZONE 2
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 2)) THEN
               CALL I_Q(1,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 302 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 Q_PROF(7,ILAY)/Q_PROF(1,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 302           CONTINUE
C ................................................................
C        OZONE 3
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 3)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 303 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(7,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 303           CONTINUE
C ................................................................
C        OZONE 4
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 4)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(10,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 304 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DSQRT(SECANT*Q_PROF(7,ILAY))*
     &                 Q_PROF(7,ILAY)/Q_PROF(10,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 304           CONTINUE
C ................................................................
C        OZONE 5
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 5)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 305 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 Q_PROF(7,ILAY)*Q_PROF(3,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 305           CONTINUE
C ................................................................
C        OZONE 6
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 6)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(11,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 306 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(7,ILAY)*Q_PROF(11,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 306           CONTINUE
C ................................................................
C        OZONE 7
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 7)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 307 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(7,ILAY))**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 307           CONTINUE
C ................................................................
C        OZONE 8
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 8)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(8,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 308 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*SECANT* 
     &                 Q_PROF(7,ILAY)*DSQRT(SECANT*Q_PROF(8,ILAY))
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 308           CONTINUE
C ................................................................
C        OZONE 9
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 9)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 309 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(3,ILAY)*
     &                 DSQRT(SECANT*Q_PROF(7,ILAY))
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 309           CONTINUE
C ................................................................
C        OZONE 10
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 10)) THEN
               CALL I_Q(8,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 310 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(8,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 310           CONTINUE
C ................................................................
C        OZONE 11
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 11)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(9,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 311 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(9,ILAY)*
     &                 DSQRT(SECANT*Q_PROF(7,ILAY))
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 311           CONTINUE
C ................................................................
C        OZONE 12
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 12)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(9,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 312 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(7,ILAY)*Q_PROF(9,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 312           CONTINUE
C ................................................................
C        OZONE 13
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 13)) THEN
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 313 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(4,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 313           CONTINUE
C ................................................................
C        OZONE 14
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 14)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(10,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 314 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*(Q_PROF(7,ILAY)**2)/Q_PROF(10,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 314           CONTINUE
C ................................................................
C        OZONE 15
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 15)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(4,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 315 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT**2)*Q_PROF(7,ILAY)*Q_PROF(4,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 315           CONTINUE
C ................................................................
C        OZONE 16
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 16)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 316 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 SECANT*Q_PROF(7,ILAY)*Q_PROF(3,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 316           CONTINUE
C ................................................................
C        OZONE 17
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 17)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 317 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(7,ILAY))**1.25
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 317           CONTINUE
C ................................................................
C        OZONE 18
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 18)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 318 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(7,ILAY))**1.5
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 318           CONTINUE
C ................................................................
C        OZONE 19
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 19)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 319 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(3,ILAY)*
     &                 (SECANT*Q_PROF(7,ILAY))**1.5
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 319           CONTINUE
C ................................................................
C        OZONE 20
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 20)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 320 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 (SECANT*Q_PROF(7,ILAY))**2.5
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 320           CONTINUE
C ................................................................
C        OZONE 21
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 21)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               CALL I_Q(3,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 321 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*Q_PROF(3,ILAY)*
     &                 (SECANT*Q_PROF(7,ILAY))**2
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 321           CONTINUE
C ................................................................
C        OZONE 22
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 22)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 322 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DLOG(1.1*SECANT) *Q_PROF(7,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 322           CONTINUE
C ................................................................
C        OZONE 23
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 23)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 323 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DLOG(1.1*SECANT*Q_PROF(7,ILAY))
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 323           CONTINUE
C ................................................................
C        OZONE 24
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .EQ. 24)) THEN
               CALL I_Q(7,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,
     &              AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
               DO 324 ILAY=1,NLAYER
                  KADD(ILAY)=COEFF(IIRF,IMOL,IPRED,ILAY)*
     &                 DEXP(SECANT/2.0)*Q_PROF(7,ILAY)
                  K_MOL(ILAY)=K_MOL(ILAY)+KADD(ILAY)
 324           CONTINUE
C ................................................................
C        OZONE - UNDEFINED
C ................................................................
            ELSE IF ((IMOL .EQ. 3) .AND. (PID .GT. 24)) THEN
               WRITE (*,*) 'FATAL ERROR IN C_K:'
               WRITE (*,*) '((IMOL .EQ. 3) .AND. (PID .GT. 24))'
               WRITE (*,*) 'PROGRAM ABORTED'
               STOP
C ................................................................
C        ALL - UNDEFINED
C ................................................................
            ELSE
               WRITE (*,*) 'FATAL ERROR IN C_K:'
               WRITE (*,*) 'INVALID MOLECULE OR PREDICTOR ID'
               WRITE (*,*) 'PROGRAM ABORTED'
               STOP
            ENDIF
           
 600     CONTINUE
C ................................................................
C ................................................................

         ! TUNING:
         IF (IMOL .EQ. 1) THEN
C           TUNING=DLOG(T_TAUF(IIRF)) ! if tuning is transm * tuning
            TUNING=T_TAUF(IIRF)       ! if tuning is transm ** tuning
         ELSE IF (IMOL .EQ. 2) THEN
C           TUNING=DLOG(T_TAUW(IIRF))
            TUNING=T_TAUW(IIRF)
         ELSE IF (IMOL .EQ. 3) THEN
C           TUNING=DLOG(T_TAUO(IIRF))
            TUNING=T_TAUO(IIRF)
         ELSE
            WRITE (*,*) 'FATAL ERROR IN C_K DURING TUNING:'
            WRITE (*,*) 'IMOL IS NEITHER 1,2,3'
            WRITE (*,*) 'PROGRAM ABORTED'
            STOP
         ENDIF

         ! now add the k-profile of the actual molecule to the total k
         ! DO NOT add anything if K_MOL(ILAY) is negative
         ! also add TUNING and set K_NEG (that is needed for Jacobians)
         DO 800 ILAY=1,NLAYER
            IF (K_MOL(ILAY) .GT. 0.0) THEN
               
C              ADDTOK=K_MOL(ILAY)+TUNING ! if tuning is transm * tuning
               ADDTOK=K_MOL(ILAY)*TUNING ! if tuning is transm ** tuning

C     ---------------------------------------------------------------------
C              !!! THE FOLLOWING LOOP ONLY IF TUNING IS TRANSM * TUNING !!!
C              IF (T_ADJ .EQ. 1) THEN
C                 ! ----------------------------------------------------
C                 ! CHECK / CHANGE TUNING PARAMETERS
C                 ! ----------------------------------------------------
C                 ! This if-loop prevents from adding a negative number
C                 ! to the total absorption coefficient (K_PROF) for a
C                 ! too small tuning parameter. If the tuning parameter is
C                 ! found to be too small, it will be adjusted to the 
C                 ! smallest possible number.
C                 ! This if-loop will be entered only if this is the first
C                 ! run of C_TAU. It is assumed, that a 2nd, 3rd... call of 
C                 ! C_TAU is for finite differences Jacobians, where the
C                 ! tuning parameters must not be changed anymore.
C                 ! NOTE: using this loop would cause that the tuning
C                 ! parameters would usually be changed if they are much
C                 ! smaller than 1.
C                 IF ((ADDTOK .LT. 0.0) .AND. (SCOUNT .EQ. 0)) THEN
C                    TUNING = -K_MOL(ILAY)
C                    IF (IMOL .EQ. 1) THEN
C                       T_TAUF(IIRF)=DEXP(TUNING)
C                       IF (VERBOS .EQ. 1) 
C    &                       WRITE (*,801) IIRF,T_TAUF(IIRF),ILAY
C                    ELSE IF (IMOL .EQ. 2) THEN
C                       T_TAUW(IIRF)=DEXP(TUNING)
C                       IF (VERBOS .EQ. 1) 
C    &                       WRITE (*,802) IIRF,T_TAUW(IIRF),ILAY
C                    ELSE IF (IMOL .EQ. 3) THEN
C                       T_TAUO(IIRF)=DEXP(TUNING)
C                       IF (VERBOS .EQ. 1) 
C    &                       WRITE (*,803) IIRF,T_TAUO(IIRF),ILAY
C                    ENDIF
C                    ADDTOK=0.0
C                    CHTTAU=CHTTAU+1
C                 ENDIF
C              ! ---- END OF CHECK / CHANGE TUNING PARAMETERS -------
C              ENDIF ! of: IF (T_ADJ .EQ. 1) THEN
C     ---------------------------------------------------------------------

               K_PROF(ILAY)=K_PROF(ILAY)+ADDTOK
               K_NEG(IIRF,IMOL,ILAY)=0
            ELSE
               K_NEG(IIRF,IMOL,ILAY)=1
            ENDIF
 800     CONTINUE

 700  CONTINUE

 801  FORMAT(' WARNING: ADJUST TUNING OF TAU(FIXED) FOR',
     &       ' CHANNEL ',I2,' TO ',F10.5,' (LAYER ',I3,')')
 802  FORMAT(' WARNING: ADJUST TUNING OF TAU(WATER) FOR',
     &       ' CHANNEL ',I2,' TO ',F10.5,' (LAYER ',I3,')')
 803  FORMAT(' WARNING: ADJUST TUNING OF TAU(OZONE) FOR',
     &       ' CHANNEL ',I2,' TO ',F10.5,' (LAYER ',I3,')')
      RETURN
      END

C     END OF SUBROUTINE C_K



C ====================================================================

      SUBROUTINE I_Q(QPID,QPINIT,Q_PROF,ATEMP,AFIXED,AWATER,AOZONE,
     &               RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     initialize Q-profile number QPID
C
C     INPUT:
C       QPID            ID of Q-profile which shall be initialized
C       ATEMP,AFIXED,AWATER,AOZONE       atmosphere profile          
C       RTEMP,RFIXED,RWATER,ROZONE,PRES  reference atmosphere profile
C       VERBOS          screen output on/off (1=on,0=off)
C
C     INPUT AND OUTPUT:
C       QPINIT          Q-profiles initialization flags (see COMMENT)
C       Q_PROF          Q-profiles (see COMMENT)
C
C     CALLING ROUNTINE
C       C_K
C
C     COMMENTS:
C       The Q-profiles are the profiles which are needed to calculate
C       Q (for example the "layer-above-profiles"); Q-profiles are NOT
C       the Q-matrix!
C
C     HISTORY:
C       written 9/2/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE CALLING VARIABLES
      INTEGER     QPID,QPINIT(MAX_QP),VERBOS
      REAL*8      Q_PROF(MAX_QP,NLAYER),
     &            ATEMP(NLAYER),AFIXED(NLAYER),AWATER(NLAYER),
     &            AOZONE(NLAYER),RTEMP(NLAYER),RFIXED(NLAYER),
     &            RWATER(NLAYER),ROZONE(NLAYER),PRES(NLAYER)
C     LOCAL VARIABLES
      INTEGER     ILAY
      REAL*8      RDMMY1,RDMMY2,RDMMY3

C     check if initialization flag is already set
      IF (QPINIT(QPID) .EQ. 1) THEN
         RETURN
      ELSE
         IF (VERBOS .EQ. 1) THEN
            WRITE (*,20) QPID
 20         FORMAT(' initializing Q-profile ',I2)
         ENDIF
      ENDIF
      

C     calculate Q-profiles:
C     NOTE: Q-profiles are not initialized by any other subroutine!

C............................................................................
C     Q-PROFILE 1
C............................................................................
      IF (QPID .EQ. 1) THEN
         DO 1010 ILAY=1,NLAYER
            Q_PROF(1,ILAY)=ATEMP(ILAY) / RTEMP(ILAY)
 1010    CONTINUE
         QPINIT(1)=1
C............................................................................
C     Q-PROFILE 2
C............................................................................
      ELSE IF (QPID .EQ. 2) THEN
         ! Q-profile 2 needs Q-profile 1
         IF (QPINIT(1) .NE. 1) THEN
            DO 1021 ILAY=1,NLAYER
               Q_PROF(1,ILAY)=ATEMP(ILAY) / RTEMP(ILAY)
 1021       CONTINUE
            QPINIT(1)=1
         ENDIF
         ! now do Q-profile 2
         Q_PROF(2,NLAYER)=0.0
         DO 1022 ILAY=(NLAYER-1),1,-1
            Q_PROF(2,ILAY)=Q_PROF(2,ILAY+1) +
     &                     PRES(ILAY) * (PRES(ILAY)-PRES(ILAY+1)) *
     &                     Q_PROF(1,ILAY)
 1022    CONTINUE
         QPINIT(2)=1
C............................................................................
C     Q-PROFILE 3
C............................................................................
      ELSE IF (QPID .EQ. 3) THEN
         DO 1030 ILAY=1,NLAYER
            Q_PROF(3,ILAY)=ATEMP(ILAY) - RTEMP(ILAY)
 1030    CONTINUE
         QPINIT(3)=1
C............................................................................
C     Q-PROFILE 4
C............................................................................
      ELSE IF (QPID .EQ. 4) THEN
         DO 1040 ILAY=1,NLAYER
            Q_PROF(4,ILAY)=AWATER(ILAY) / RWATER(ILAY)
 1040    CONTINUE
         QPINIT(4)=1
C............................................................................
C     Q-PROFILE 5
C............................................................................
      ELSE IF (QPID .EQ. 5) THEN
         Q_PROF(5,NLAYER)=AWATER(NLAYER) / RWATER(NLAYER)
         RDMMY1=0.0
         RDMMY2=0.0
         DO 1050 ILAY=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY)*(PRES(ILAY)-PRES(ILAY+1))
            RDMMY1=RDMMY1+RDMMY3*AWATER(ILAY)
            RDMMY2=RDMMY2+RDMMY3*RWATER(ILAY)
            Q_PROF(5,ILAY)=RDMMY1/RDMMY2
 1050    CONTINUE
         QPINIT(5)=1
C............................................................................
C     Q-PROFILE 6
C............................................................................
      ELSE IF (QPID .EQ. 6) THEN
         Q_PROF(6,NLAYER)=(AWATER(NLAYER)*ATEMP(NLAYER)) / 
     &                    (RWATER(NLAYER)*RTEMP(NLAYER))
         RDMMY1=0.0
         RDMMY2=0.0
         DO 1060 ILAY=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY)*(PRES(ILAY)-PRES(ILAY+1))
            RDMMY1=RDMMY1+RDMMY3*AWATER(ILAY)*ATEMP(ILAY)
            RDMMY2=RDMMY2+RDMMY3*RWATER(ILAY)*RTEMP(ILAY)
            Q_PROF(6,ILAY)=RDMMY1/RDMMY2
 1060    CONTINUE
         QPINIT(6)=1
C............................................................................
C     Q-PROFILE 7
C............................................................................
      ELSE IF (QPID .EQ. 7) THEN
         !write (*,*) AOZONE ! ********
         DO 1070 ILAY=1,NLAYER
            !write (*,*) '1: ',ilay,aozone(ilay),rozone(ilay) ! *******
            Q_PROF(7,ILAY)=AOZONE(ILAY) / ROZONE(ILAY)
 1070    CONTINUE
         QPINIT(7)=1
C............................................................................
C     Q-PROFILE 8
C............................................................................
      ELSE IF (QPID .EQ. 8) THEN
         ! Q-profile 8 needs Q-profile 7
         IF (QPINIT(7) .NE. 1) THEN
            DO 1081 ILAY=1,NLAYER
               !write (*,*) '2: ',ilay,aozone(ilay),rozone(ilay) ! *******
               Q_PROF(7,ILAY)=AOZONE(ILAY) / ROZONE(ILAY)
 1081       CONTINUE
            QPINIT(7)=1
         ENDIF 
         ! now do Q-profile 8
         Q_PROF(8,NLAYER)=0.0
         DO 1082 ILAY=(NLAYER-1),1,-1
            Q_PROF(8,ILAY)=Q_PROF(8,ILAY+1) +
     &                     PRES(ILAY) * (PRES(ILAY)-PRES(ILAY+1)) *
     &                     Q_PROF(7,ILAY+1)
 1082    CONTINUE
         QPINIT(8)=1
C............................................................................
C     Q-PROFILE 9
C............................................................................
      ELSE IF (QPID .EQ. 9) THEN
         ! Q-profile 9 needs Q-profile 1 and Q-profile 7
         IF (QPINIT(1) .NE. 1) THEN
            DO 1091 ILAY=1,NLAYER
               Q_PROF(1,ILAY)=ATEMP(ILAY) / RTEMP(ILAY)
 1091       CONTINUE
            QPINIT(1)=1
         ENDIF
         IF (QPINIT(7) .NE. 1) THEN
            DO 1092 ILAY=1,NLAYER
               !write (*,*) '3: ',ilay,aozone(ilay),rozone(ilay) ! *******
               Q_PROF(7,ILAY)=AOZONE(ILAY) / ROZONE(ILAY)
 1092       CONTINUE
            QPINIT(7)=1
         ENDIF 
         ! now do Q-profile 9
         Q_PROF(9,NLAYER)=0.0
         DO 1093 ILAY=(NLAYER-1),1,-1
            Q_PROF(9,ILAY)=Q_PROF(9,ILAY+1) +
     &                     PRES(ILAY) * (PRES(ILAY)-PRES(ILAY+1)) *
     &                     Q_PROF(7,ILAY+1) * Q_PROF(1,ILAY+1)
 1093    CONTINUE
         QPINIT(9)=1
C............................................................................
C     Q-PROFILE 10
C............................................................................
      ELSE IF (QPID .EQ. 10) THEN
         Q_PROF(10,NLAYER)=AOZONE(NLAYER)/ROZONE(NLAYER)
         RDMMY1=0.0
         RDMMY2=0.0
         DO 1100 ILAY=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY)*(PRES(ILAY)-PRES(ILAY+1))
            RDMMY1=RDMMY1+RDMMY3*AOZONE(ILAY)
            RDMMY2=RDMMY2+RDMMY3*ROZONE(ILAY)
            Q_PROF(10,ILAY)=RDMMY1/RDMMY2
 1100    CONTINUE
         QPINIT(10)=1
C............................................................................
C     Q-PROFILE 11
C............................................................................
      ELSE IF (QPID .EQ. 11) THEN
         Q_PROF(11,NLAYER)=(AOZONE(NLAYER)*ATEMP(NLAYER)) /
     &                     (ROZONE(NLAYER)*RTEMP(NLAYER))
         RDMMY1=0.0
         RDMMY2=0.0
         DO 1110 ILAY=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY)*(PRES(ILAY)-PRES(ILAY+1))
            RDMMY1=RDMMY1+RDMMY3*AOZONE(ILAY)*ATEMP(ILAY)
            RDMMY2=RDMMY2+RDMMY3*ROZONE(ILAY)*RTEMP(ILAY)
            Q_PROF(11,ILAY)=RDMMY1/RDMMY2
 1110    CONTINUE
         QPINIT(11)=1
C............................................................................
      ELSE
         WRITE (*,*) 'FATAL ERROR IN I_Q: NO SUCH Q-PROFILE'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF

      RETURN
      END
C     END OF SUBROUTINE I_Q

C ====================================================================


      SUBROUTINE C_DTAU(CLCJAC,SATNUM,TOTIRF,ID_IRF,PREDID,
     &                  NCOEFF,COEFF,SECANT,T_SEC,SUNSEC,
     &                  TAU,TAU_RS,
     &                  ATEMP,AFIXED,AWATER,AOZONE,
     &                  RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &                  DTAUT,DTAUW,DTAUO,DTAUTS,DTAUWS,DTAUOS,
     &                  K_NEG,K_NEGS,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     calculates the derivatives of TAU (layer transmittances)
C     for the upwelling thermal radiation: dTAU/dT, dTAU/dW, dTAU/dO,
C     and the derivatives of TAU_RS for the reflected solar
C     NOTE: the output is derivatives of layer transmittances and
C           not (!) layer-to-space transmittances !
C
C     INPUT:
C       SATNUM         satellite number (e.g. 7 for NOAA7)
C       TOTIRF         total number of IRFs (channels)
C       ID_IRF         IDs of IRFs
C       PREDID         predictor identifiers
C       NCOEFF         number of coefficients
C       COEFF          coefficients
C       SECANT         secant (of observation)
C       T_SEC          channel dependend secant offset (tuning parameter)
C       SUNSEC         secant of sun
C       ATEMP          atmosphere temperature profile
C       AFIXED         atmosphere fixed gases profile
C       AWATER         atmosphere water profile
C       AOZONE         atmosphere ozone profile
C       RTEMP          reference atmosphere temperature profile
C       RFIXED         reference atmosphere fixed gases profile
C       RWATER         reference atmosphere water profile
C       ROZONE         reference atmosphere ozone profile
C       PRES           pressure profile
C       TAU            transmittances for upwelling thermal
C       TAU_RS         transmittances for reflected solar
C       K_NEG          indecees if K has been disregarded in SUBROUTINE C_K
C                      (=1) because it was negative, or used (=0)
C       K_NEGS         (same for solar)
C       VERBOS         screen output on/off (1=on,0=off)
C
C     OUTPUT:
C       DTAUT          derivatives of transmittances for upwelling thermal
C       DTAUW            for temperature, water, ozone, respectively 
C       DTAUO
C       DTAUTS         derivatives of transmittances for reflected solar
C       DTAUWS
C       DTAUOS
C
C     CALLING ROUTINE:
C       KERNEL
C
C     COMMENTS:
C       for k > MAX_K (locally defined constant) the corresponding
C       TAU and TAU_RS will be set to zero
C
C     HISTORY:
C       written 12/03/97 by Tobias Wehr
C --------------------------------------------------------------------
      IMPLICIT NONE
      include "hffp_glob_dec.f"

C     SUBROUTINE ARGUMENT PARAMETER (INPUT)
      INTEGER        CLCJAC(3), ! calc. Jacobians for T,W,O
     &               SATNUM,TOTIRF,ID_IRF(MAXIRF),
     &               PREDID(MAXIRF,MAXMOL,MAXPRD+1),
     &               NCOEFF(MAXIRF,MAXMOL),
     &               K_NEG(MAXIRF,MAXMOL,NLAYER),
     &               K_NEGS(MAXIRF,MAXMOL,NLAYER),
     &               VERBOS
      REAL*8         TAU(MAXIRF,NLAYER),     
     &               TAU_RS(MAXIRF,NLAYER),  
     &               SECANT,T_SEC(MAXIRF),SUNSEC,
     &               ATEMP(NLAYER),AFIXED(NLAYER),
     &               AWATER(NLAYER),AOZONE(NLAYER),
     &               RTEMP(NLAYER),RFIXED(NLAYER),
     &               RWATER(NLAYER),ROZONE(NLAYER),
     &               PRES(NLAYER),
     &               COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER)
C     SUBROUTINE ARGUMENT PARAMETER (OUTPUT)
      REAL*8         DTAUT(MAXIRF,NLAYER,NLAYER),    ! dtau/dT, thermal upw.
     &               DTAUW(MAXIRF,NLAYER,NLAYER),    ! dtau/dW, thermal upw.
     &               DTAUO(MAXIRF,NLAYER,NLAYER),    ! dtau/dO, thermal upw.
     &               DTAUTS(MAXIRF,NLAYER,NLAYER),   ! dtau/dT, refl. solar
     &               DTAUWS(MAXIRF,NLAYER,NLAYER),   ! dtau/dW, refl. solar
     &               DTAUOS(MAXIRF,NLAYER,NLAYER)    ! dtau/dO, refl. solar
C     LOCAL PARAMETER

      ! Q-profiles are profiles used for the predictor calculations. 
      ! They are composed from temperature, pressure, water- and 
      ! ozone-profiles from both the actual and the reference 
      ! atmosphere profile. The Q-profiles contains layer-above-
      ! information.
C      REAL*8         Q_PROF(MAX_QP,NLAYER)
      ! the derivatives of the Q-profiles with respect to
      ! temperature, water, ozone.
      REAL*8         DQDT(MAX_QP,NLAYER,NLAYER)
      REAL*8         DQDW(MAX_QP,NLAYER,NLAYER)
      REAL*8         DQDO(MAX_QP,NLAYER,NLAYER)


C     DQDX_I(j)=1 if DQDX-profile has been already initialized
      ! this variable shall avoid unnecessary calculation of a 
      ! Q-profile. A Q-profile will only be calculated if is 
      ! will be used. If it has already being used before, this 
      ! flag is set to 1 and it will not be calculated anymore.
      INTEGER        DQDT_I(MAX_QP),DQDW_I(MAX_QP),DQDO_I(MAX_QP) 

C     k-profiles
      REAL*8         DKDT(NLAYER,NLAYER)   ! d(abscoeff upwelling)/dT
      REAL*8         DKDW(NLAYER,NLAYER)   ! d(abscoeff upwelling)/dW
      REAL*8         DKDO(NLAYER,NLAYER)   ! d(abscoeff upwelling)/dO
      REAL*8         DKDT_S(NLAYER,NLAYER) ! d(abscoeff solar)/dT
      REAL*8         DKDW_S(NLAYER,NLAYER) ! d(abscoeff solar)/dW
      REAL*8         DKDO_S(NLAYER,NLAYER) ! d(abscoeff solar)/dO

      INTEGER        RSWARN            ! warning already given (=1)
      INTEGER        IIRF,IIRF_A,ILAY1,ILAY2
      INTEGER        IHELP
C     corrected ("tuned") secant (SECANT + T_SEC)
      REAL*8         C_SEC
C     effective secant (sun + satellite) for refl. solar
      REAL*8         EFFSEC
C     maximum value for abs.coeff. k, for k>MAX_K we set tau=0
      REAL*8         MAX_K
      PARAMETER(MAX_K=30.0)
C     INITIALIZATIONS
      

C     ============================================================
C     !                        BEGIN                             !
C     ============================================================

C     INIT RSWARN
      ! if solar reflected will not be calculated due to too large
      ! effective angle, this parameter will be set to 1 and refl.
      ! solar will not be considered for this job anymore
      RSWARN=0

C     RESET QPINIT TO ZERO
      IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
         DO 10 IHELP=1,MAX_QP   ! do not use a "DATA" statement for this
            DQDT_I(IHELP)=0     ! (see C_TAU for reason)
            DO 11 ILAY1=1,NLAYER
               DO 12 ILAY2=1,NLAYER
                  DQDT(IHELP,ILAY1,ILAY2)=0.0
 12            CONTINUE
 11         CONTINUE
 10      CONTINUE
      ENDIF
      IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
         DO 20 IHELP=1,MAX_QP   ! do not use a "DATA" statement for this
            DQDW_I(IHELP)=0     ! (see C_TAU for reason)
            DO 21 ILAY1=1,NLAYER
               DO 22 ILAY2=1,NLAYER
                  DQDW(IHELP,ILAY1,ILAY2)=0.0
 22            CONTINUE
 21         CONTINUE
 20      CONTINUE
      ENDIF
      IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
         DO 30 IHELP=1,MAX_QP   ! do not use a "DATA" statement for this
            DQDO_I(IHELP)=0     ! (see C_TAU for reason)
            DO 31 ILAY1=1,NLAYER
               DO 32 ILAY2=1,NLAYER
                  DQDO(IHELP,ILAY1,ILAY2)=0.0
 32            CONTINUE
 31         CONTINUE
 30      CONTINUE
      ENDIF
C     loop desired channels (IRFs)
      DO 1000 IIRF=1,TOTIRF
         ! reset K_PROF for each channel
         ! ATTENTION: loop both ILAY1 and ILAY2 from 1 to NLAYER, because
         !   it will be used to calculate DTAUX (see "ATTENTION" below)
         IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
            DO 1012 ILAY1=1,NLAYER
               DO 1013 ILAY2=1,NLAYER
                  DKDT(ILAY1,ILAY2)=0.0
                  DKDT_S(ILAY1,ILAY2)=0.0
 1013          CONTINUE
 1012       CONTINUE
         ENDIF
         IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
            DO 1022 ILAY1=1,NLAYER
               DO 1023 ILAY2=1,NLAYER
                  DKDW(ILAY1,ILAY2)=0.0
                  DKDW_S(ILAY1,ILAY2)=0.0
 1023          CONTINUE
 1022       CONTINUE
         ENDIF
         IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
            DO 1032 ILAY1=1,NLAYER
               DO 1033 ILAY2=1,NLAYER
                  DKDO(ILAY1,ILAY2)=0.0
                  DKDO_S(ILAY1,ILAY2)=0.0
 1033          CONTINUE
 1032       CONTINUE
         ENDIF
         IIRF_A=ID_IRF(IIRF)    ! IIRF_A contains now the channel ID number
         CALL ADDSEC(SECANT,T_SEC(IIRF_A),C_SEC)

         IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
            CALL C_DKDT(IIRF_A,PREDID,NCOEFF,COEFF,
     &           DQDT_I,DQDT,DKDT,C_SEC,
     &           ATEMP,AFIXED,AWATER,AOZONE,
     &           RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEG,VERBOS)
         ENDIF
         IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
            CALL C_DKDW(IIRF_A,PREDID,NCOEFF,COEFF,
     &           DQDW_I,DQDW,DKDW,C_SEC,
     &           ATEMP,AFIXED,AWATER,AOZONE,
     &           RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEG,VERBOS)
         ENDIF
         IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
            CALL C_DKDO(IIRF_A,PREDID,NCOEFF,COEFF,
     &           DQDO_I,DQDO,DKDO,C_SEC,
     &           ATEMP,AFIXED,AWATER,AOZONE,
     &           RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEG,VERBOS)
         ENDIF
C        REFLECTED SOLAR
         IF ((IIRF_A .GE. SNIIRF) .AND. (RSWARN .EQ. 0)) THEN
            ! calculate EFFSEC first as the effective secant
            CALL ADDSEC(C_SEC,SUNSEC,EFFSEC)
            !write (*,*) 'EFFSEC=',EFFSEC
            IF (EFFSEC .GT. MAXMSV) THEN
               RSWARN=1
            ELSE
               ! CALCULATE ABSORPTION COEFFICIENTS FOR REFLECTED SOLAR
               IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
                  CALL C_DKDT(IIRF_A,PREDID,NCOEFF,COEFF,
     &                 DQDT_I,DQDT,DKDT_S,EFFSEC,
     &                 ATEMP,AFIXED,AWATER,AOZONE,
     &                 RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEGS,VERBOS)
               ENDIF
               IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
                  CALL C_DKDW(IIRF_A,PREDID,NCOEFF,COEFF,
     &                 DQDW_I,DQDW,DKDW_S,EFFSEC,
     &                 ATEMP,AFIXED,AWATER,AOZONE,
     &                 RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEGS,VERBOS)
               ENDIF
               IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
                  CALL C_DKDO(IIRF_A,PREDID,NCOEFF,COEFF,
     &                 DQDO_I,DQDO,DKDO_S,EFFSEC,
     &                 ATEMP,AFIXED,AWATER,AOZONE,
     &                 RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEGS,VERBOS)
               ENDIF
            ENDIF
         ENDIF

C        calculated d(tau)/dX   X=T,W,O
C        dtau(i)/dX(j) = -tau(i) * dk(i)/dX(j)
         ! ATTENTION: loop through all NLAYER layers of both indecees
         !            ILAY1 and ILAY2, since DTAUX and DTAUXS has not
         !            been initialized before!
         ! thermal upwelling
         IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
            DO 2011 ILAY1=1,NLAYER
               DO 2012 ILAY2=1,NLAYER
                  DTAUT(IIRF_A,ILAY1,ILAY2)=
     &                 (-1.0)*TAU(IIRF_A,ILAY1)*DKDT(ILAY1,ILAY2)
 2012          CONTINUE
 2011       CONTINUE
         ENDIF
         IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
            DO 2021 ILAY1=1,NLAYER
               DO 2022 ILAY2=1,NLAYER
                  DTAUW(IIRF_A,ILAY1,ILAY2)=
     &                 (-1.0)*TAU(IIRF_A,ILAY1)*DKDW(ILAY1,ILAY2)
 2022          CONTINUE
 2021       CONTINUE
         ENDIF
         IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
            DO 2031 ILAY1=1,NLAYER
               DO 2032 ILAY2=1,NLAYER
                  DTAUO(IIRF_A,ILAY1,ILAY2)=
     &                 (-1.0)*TAU(IIRF_A,ILAY1)*DKDO(ILAY1,ILAY2)
 2032          CONTINUE
 2031       CONTINUE
         ENDIF
         ! reflected solar
         IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
            IF (IIRF_A .GE. SNIIRF) THEN
               DO 3013 ILAY1=1,NLAYER
                  DO 3014 ILAY2=1,NLAYER
                     DTAUTS(IIRF_A,ILAY1,ILAY2)=
     &                 (-1.0)*TAU_RS(IIRF_A,ILAY1)*DKDT_S(ILAY1,ILAY2)
 3014             CONTINUE
 3013          CONTINUE
            ELSE
               DO 3015 ILAY1=1,NLAYER
                  DO 3016 ILAY2=1,NLAYER
                     DTAUTS(IIRF_A,ILAY1,ILAY2)=0.0
 3016             CONTINUE
 3015          CONTINUE
            ENDIF
         ENDIF
         IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
            IF (IIRF_A .GE. SNIIRF) THEN
               DO 3023 ILAY1=1,NLAYER
                  DO 3024 ILAY2=1,NLAYER
                     DTAUWS(IIRF_A,ILAY1,ILAY2)=
     &                 (-1.0)*TAU_RS(IIRF_A,ILAY1)*DKDW_S(ILAY1,ILAY2)
 3024             CONTINUE
 3023          CONTINUE
            ELSE
               DO 3025 ILAY1=1,NLAYER
                  DO 3026 ILAY2=1,NLAYER
                     DTAUWS(IIRF_A,ILAY1,ILAY2)=0.0
 3026             CONTINUE
 3025          CONTINUE
            ENDIF
         ENDIF
         IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
            IF (IIRF_A .GE. SNIIRF) THEN
               DO 3033 ILAY1=1,NLAYER
                  DO 3034 ILAY2=1,NLAYER
                     DTAUOS(IIRF_A,ILAY1,ILAY2)=
     &                 (-1.0)*TAU_RS(IIRF_A,ILAY1)*DKDO_S(ILAY1,ILAY2)
 3034             CONTINUE
 3033          CONTINUE
            ELSE
               DO 3035 ILAY1=1,NLAYER
                  DO 3036 ILAY2=1,NLAYER
                     DTAUOS(IIRF_A,ILAY1,ILAY2)=0.0
 3036             CONTINUE
 3035          CONTINUE
            ENDIF
         ENDIF

 1000 CONTINUE ! end of DO 1000 IIRF=1,TOTIRF
      RETURN
      END
C     END OF C_DTAU
C ====================================================================


      SUBROUTINE C_DTAU_FD(CLCJAC,SATNUM,TOTIRF,ID_IRF,PREDID,
     &                  NCOEFF,COEFF,SECANT,T_SEC,SUNSEC,
     &                  ATEMP,AFIXED,AWATER,AOZONE,
     &                  RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &                  DTAUT,DTAUW,DTAUO,DTAUTS,DTAUWS,DTAUOS,
     &                  VERBOS,IJOB,T_TAUF,T_TAUW,T_TAUO)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     calculates the derivatives of TAU (layer transmittances)
C     for the upwelling thermal radiation: dTAU/dT, dTAU/dW, dTAU/dO,
C     and the derivatives of TAU_RS for the reflected solar
C     NOTE: the output is derivatives of layer transmittances and
C           not (!) layer-to-space transmittances !
C     THIS SUBROUTINE CALCULATES DERIVATIVES OF TAU USING THE
C     FINITE DIFFERENCES METHODE
C
C     INPUT:
C       SATNUM         satellite number (e.g. 7 for NOAA7)
C       TOTIRF         total number of IRFs (channels)
C       ID_IRF         IDs of IRFs
C       SECANT         secant (of observation)
C       T_SEC          channel dependend secant offset (tuning parameter)
C       SUNSEC         secant of sun
C       ATEMP          atmosphere temperature profile
C       AFIXED         atmosphere fixed gases profile
C       AWATER         atmosphere water profile
C       AOZONE         atmosphere ozone profile
C       RTEMP          reference atmosphere temperature profile
C       RFIXED         reference atmosphere fixed gases profile
C       RWATER         reference atmosphere water profile
C       ROZONE         reference atmosphere ozone profile
C       PRES           pressure profile
C       VERBOS         screen output on/off (1=on,0=off)
C       T_TAUF         tuning factors of tau(fixed)
C       T_TAUW         tuning factors of tau(water)
C       T_TAUO         tuning factors of tau(ozone)
C
C     OUTPUT:
C       DTAUT          derivatives of transmittances for upwelling thermal
C       DTAUW            for temperature, water, ozone, respectively 
C       DTAUO
C       DTAUTS         derivatives of transmittances for reflected solar
C       DTAUWS
C       DTAUOS
C
C     CALLING ROUTINE:
C       KERNEL
C
C     COMMENTS:
C       for k > MAX_K (locally defined constant) the corresponding
C       TAU and TAU_RS will be set to zero
C
C     HISTORY:
C       written 12/03/97 by Tobias Wehr
C --------------------------------------------------------------------
      IMPLICIT NONE
      include "hffp_glob_dec.f"

C     SUBROUTINE ARGUMENT PARAMETER (INPUT)
      INTEGER        CLCJAC(3), ! calc. Jacobians for T,W,O
     &               SATNUM,TOTIRF,ID_IRF(MAXIRF),
     &               PREDID(MAXIRF,MAXMOL,MAXPRD+1),
     &               NCOEFF(MAXIRF,MAXMOL),
     &               VERBOS,IJOB
      REAL*8         SECANT,T_SEC(MAXIRF),SUNSEC,
     &               ATEMP(NLAYER),AFIXED(NLAYER),
     &               AWATER(NLAYER),AOZONE(NLAYER),
     &               RTEMP(NLAYER),RFIXED(NLAYER),
     &               RWATER(NLAYER),ROZONE(NLAYER),
     &               PRES(NLAYER),
     &               COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER),
     &               T_TAUF(MAXIRF),
     &               T_TAUW(MAXIRF),
     &               T_TAUO(MAXIRF)
C     SUBROUTINE ARGUMENT PARAMETER (OUTPUT)
      REAL*8         DTAUT(MAXIRF,NLAYER,NLAYER),    ! dtau/dT, thermal upw.
     &               DTAUW(MAXIRF,NLAYER,NLAYER),    ! dtau/dW, thermal upw.
     &               DTAUO(MAXIRF,NLAYER,NLAYER),    ! dtau/dO, thermal upw.
     &               DTAUTS(MAXIRF,NLAYER,NLAYER),   ! dtau/dT, refl. solar
     &               DTAUWS(MAXIRF,NLAYER,NLAYER),   ! dtau/dW, refl. solar
     &               DTAUOS(MAXIRF,NLAYER,NLAYER)    ! dtau/dO, refl. solar
C     LOCAL PARAMETER
      INTEGER        IIRF,IIRF_A,ILAY1,ILAY2,ILAY3
      INTEGER        IHELP,IDUMMY
C     effective secant (sun + satellite) for refl. solar
      REAL*8         EFFSEC

      REAL*8         TAU(MAXIRF,NLAYER),    ! undist. tau
     &               TAU_RS(MAXIRF,NLAYER), ! undist. tau refl. solar
     &               TAUD(MAXIRF,NLAYER),   ! disturbed tau
     &               TAURSD(MAXIRF,NLAYER), ! disturbed tau refl. solar
     &               DMARR1(MAXIRF,MAXMOL,NLAYER),  ! dummy
     &               DMARR2(MAXIRF,MAXMOL,NLAYER),  ! dummy
     &               DATM(NLAYER),          ! disturbed atmosphere
     &               DELTA,DELTAT,DELTAW,DELTAO ! deltas
      DATA DELTAT/1.0/  ! 1 K
      DATA DELTAW/1.0/  ! 1 %
      DATA DELTAO/1.0/  ! 1 %

C     BEGIN
      WRITE (*,*) 'CALCULATE D(TAU)/DX WITH FINITE DIFFERENCES'

      ! CALCULATE "UNDISTURBED" TAU
      CALL C_TAU(IJOB,SATNUM,TOTIRF,ID_IRF,PREDID,
     &           NCOEFF,COEFF,SECANT,SUNSEC,
     &           T_TAUF,T_TAUW,T_TAUO,T_SEC,TAU,TAU_RS,
     &           ATEMP,AFIXED,AWATER,AOZONE,
     &           RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &           DMARR1,DMARR2,0,0,IDUMMY)

      ! calculate d(tau)/dT
      IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
         WRITE (*,*) '   CALCULATE D(TAU)/DT'
         DO 1100 ILAY1=1,NLAYER
            DO 1200 ILAY2=1,NLAYER
               DO 1201 ILAY3=1,NLAYER
                  DATM(ILAY3)=ATEMP(ILAY3)
 1201          CONTINUE
               DELTA=DELTAT
               DATM(ILAY2)=ATEMP(ILAY2)+DELTA
               CALL C_TAU(IJOB,SATNUM,TOTIRF,ID_IRF,PREDID,
     &              NCOEFF,COEFF,SECANT,SUNSEC,
     &              T_TAUF,T_TAUW,T_TAUO,T_SEC,TAUD,TAURSD,
     &              DATM,AFIXED,AWATER,AOZONE,
     &              RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &              DMARR1,DMARR2,0,0,IDUMMY)
               DO 1300 IIRF=1,TOTIRF   
                  IIRF_A=ID_IRF(IIRF)
                  DTAUT(IIRF_A,ILAY1,ILAY2)=
     &                 (TAUD(IIRF_A,ILAY1)-TAU(IIRF_A,ILAY1))/DELTA
                  DTAUTS(IIRF_A,ILAY1,ILAY2)=
     &                 (TAURSD(IIRF_A,ILAY1)-TAU_RS(IIRF_A,ILAY1))/DELTA
 1300          CONTINUE
 1200       CONTINUE
 1100    CONTINUE
      ENDIF
      ! calculate d(tau)/dW
      IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
         WRITE (*,*) '   CALCULATE D(TAU)/DW'
         DO 2100 ILAY1=1,NLAYER
            DO 2200 ILAY2=1,NLAYER
               DO 2201 ILAY3=1,NLAYER
                  DATM(ILAY3)=AWATER(ILAY3)
 2201          CONTINUE
               DELTA=AWATER(ILAY2)*DELTAW/100.0
               DATM(ILAY2)=AWATER(ILAY2)+DELTA
               CALL C_TAU(IJOB,SATNUM,TOTIRF,ID_IRF,PREDID,
     &              NCOEFF,COEFF,SECANT,SUNSEC,
     &              T_TAUF,T_TAUW,T_TAUO,TAUD,T_SEC,TAURSD,
     &              ATEMP,AFIXED,DATM,AOZONE,
     &              RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &              DMARR1,DMARR2,0,0,IDUMMY)
               DO 2300 IIRF=1,TOTIRF   
                  IIRF_A=ID_IRF(IIRF)
                  DTAUW(IIRF_A,ILAY1,ILAY2)=
     &                 (TAUD(IIRF_A,ILAY1)-TAU(IIRF_A,ILAY1))/DELTA
                  DTAUWS(IIRF_A,ILAY1,ILAY2)=
     &                 (TAURSD(IIRF_A,ILAY1)-TAU_RS(IIRF_A,ILAY1))/DELTA
 2300          CONTINUE
 2200       CONTINUE
 2100    CONTINUE
      ENDIF
      ! calculate d(tau)/dO
      IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
         WRITE (*,*) '   CALCULATE D(TAU)/DO'
         DO 3100 ILAY1=1,NLAYER
            DO 3200 ILAY2=1,NLAYER
               DO 3201 ILAY3=1,NLAYER
                  DATM(ILAY3)=AOZONE(ILAY3)
 3201          CONTINUE
               DELTA=AOZONE(ILAY2)*DELTAO/100.0
               DATM(ILAY2)=AOZONE(ILAY2)+DELTA
               CALL C_TAU(IJOB,SATNUM,TOTIRF,ID_IRF,PREDID,
     &              NCOEFF,COEFF,SECANT,SUNSEC,
     &              T_TAUF,T_TAUW,T_TAUO,TAUD,T_SEC,TAURSD,
     &              ATEMP,AFIXED,AWATER,DATM,
     &              RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &              DMARR1,DMARR2,0,0,IDUMMY)
               DO 3300 IIRF=1,TOTIRF   
                  IIRF_A=ID_IRF(IIRF)
                  DTAUO(IIRF_A,ILAY1,ILAY2)=
     &                 (TAUD(IIRF_A,ILAY1)-TAU(IIRF_A,ILAY1))/DELTA
                  DTAUOS(IIRF_A,ILAY1,ILAY2)=
     &                 (TAURSD(IIRF_A,ILAY1)-TAU_RS(IIRF_A,ILAY1))/DELTA
 3300          CONTINUE
 3200       CONTINUE
 3100    CONTINUE
      ENDIF

      WRITE (*,*) 'CALCULATE D(TAU)/DX WITH FINITE DIFFERENCES - DONE'
      RETURN
      END
C     END OF C_DTAU_FD

C ====================================================================

      SUBROUTINE C_DKDT(IIRF,PREDID,NCOEFF,COEFF,
     &                    DQDT_I,DQDT,DKDT,SECANT,
     &                    ATEMP,AFIXED,AWATER,AOZONE,
     &                    RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEG,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     calculate derivatives dk/dT for IIRF
C     (k=absorption koefficients, T=temperature)
C
C     INPUT:
C       IIRF            IRF counter = actual IRF (channel, or channel ID)
C       PREDID          IDs of the predictors to be used (see COMMENT)
C       NCOEFF          number of coefficients to be used (see COMMENT)
C       COEFF           coefficients to be used (see COMMENT)
C       SECANT          observation secant
C       ATEMP,AFIXED,AWATER,AOZONE       atmosphere profile          
C       RTEMP,RFIXED,RWATER,ROZONE,PRES  reference atmosphere profile
C       K_NEG           indecees if K has been disregarded in SUBROUTINE C_K
C                       (=1) because it was negative, or used (=0)
C       VERBOS          screen output on/off (1=on,0=off)
C
C     INPUT AND OUTPUT:
C       DQDT_I          Q-profiles initialization flags
C       DQDT            derivatives of Q-profiles with respect to T
C
C     OUTPUT:
C       DKDT            profile of derivatives dk/dT
C
C     CALLING ROUNTINE
C       C_DTAU
C
C     COMMENTS:
C       1. PREDID,NCOEFF,COEFF contain the information for ALL channels,
C          no matter whether all channels are used or not.
C
C     HISTORY:
C       written 12/03/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE CALLING VARIABLES
      INTEGER IIRF,
     &        NCOEFF(MAXIRF,MAXMOL),
     &        DQDT_I(MAX_QP),
     &        PREDID(MAXIRF,MAXMOL,MAXPRD+1),
     &        K_NEG(MAXIRF,MAXMOL,NLAYER),
     &        VERBOS
      REAL*8  COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER),
     &        DKDT(NLAYER,NLAYER),
     &        DQDT(MAX_QP,NLAYER,NLAYER),
     &        SECANT,
     &        ATEMP(NLAYER),AFIXED(NLAYER),AWATER(NLAYER),
     &        AOZONE(NLAYER),RTEMP(NLAYER),RFIXED(NLAYER),
     &        RWATER(NLAYER),ROZONE(NLAYER),PRES(NLAYER)

C     LOCAL VARIABLES
      INTEGER IMOL,                  ! molecule counter
     &        ILAY1,ILAY2,           ! layer counters
     &        IPRED,                 ! predictor counter
     &        PID,                   ! predictor ID
     &        RCODE                  ! subroutine identifier, used in MISSIT
      REAL*8  DKADD,                 ! dk/dT-profile increment
     &        DKDT_M(NLAYER,NLAYER), ! dk/dT-profile of a certain molecule
     &        TR(NLAYER),            ! profile ATEMP/RTEMP
     &        WR(NLAYER),            ! profile AWATER/RWATER
     &        TD(NLAYER),            ! profile ATEMP-RTEMP
     &        OR(NLAYER),            ! profile AOZONE/ROZONE
     &        TWZ2(NLAYER),          ! profile T_{Wz2} (see latex doc) 
     &        RDMMY1,RDMMY2,RDMMY3

      PARAMETER(RCODE=1) ! SUBROUTINE C_DKDT=1, C_DKDW=2, C_DKDO=3

C     ===========================================================
C     BEGIN
C     ===========================================================

C     ........ CALCULATE PREDICTORS: ...........
C     initialize dk/dT-profiles with zeros
      DO 101 ILAY1=1,NLAYER
         DO 102 ILAY2=1,NLAYER
            DKDT(ILAY1,ILAY2)=0.0
 102     CONTINUE
 101  CONTINUE

C     calculate profiles TR,WR,TD,TWZ2,OR
      RDMMY1=0.0
      RDMMY2=0.0
      DO 103 ILAY1=NLAYER,1,-1
         TR(ILAY1)=ATEMP(ILAY1)/RTEMP(ILAY1)
         WR(ILAY1)=AWATER(ILAY1)/RWATER(ILAY1)
         TD(ILAY1)=ATEMP(ILAY1)-RTEMP(ILAY1)
         OR(ILAY1)=AOZONE(ILAY1)/ROZONE(ILAY1)
         IF (ILAY1 .EQ. NLAYER) THEN
            TWZ2(ILAY1)=(AWATER(ILAY1)*ATEMP(ILAY1))/
     &           (RWATER(ILAY1)*RTEMP(ILAY1))
         ELSE
            RDMMY3=PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
            RDMMY1=RDMMY1+RDMMY3*AWATER(ILAY1)*ATEMP(ILAY1)
            RDMMY2=RDMMY2+RDMMY3*RWATER(ILAY1)*RTEMP(ILAY1)
            TWZ2(ILAY1)=RDMMY1/RDMMY2
         ENDIF
 103  CONTINUE

      DO 700 IMOL=1,MAXMOL
         ! reset DKDT_M
         DO 701 ILAY1=1,NLAYER
            DO 702 ILAY2=ILAY1,NLAYER
               DKDT_M(ILAY1,ILAY2)=0.0
 702        CONTINUE
 701     CONTINUE
         
C        loop all required predictors
         DO 600 IPRED=1,PREDID(IIRF,IMOL,1) ! PREDID(IIRF,IMOL,1) is the total
                                            ! number of predictors for (IIRF,IMOL)
            ! PID is introduced for programmer's convenience:
            PID = PREDID(IIRF,IMOL,1+IPRED) 

            IF (IMOL .EQ. 1) THEN
C ................................................................
C        PREDICTORS FOR FIXED
C ................................................................
C        FIXED 1
C ................................................................
               IF (PID .EQ. 1) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1011 ILAY1=1,NLAYER
                     DO 1012 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                       0.5*SECANT*DQDT(1,ILAY1,ILAY2)/
     &                       DSQRT(TR(ILAY1))
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1011
 1012                CONTINUE
 1011             CONTINUE
C................................................................
C        FIXED 2
C ................................................................
               ELSEIF (PID .EQ. 2) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1021 ILAY1=1,NLAYER
                     DO 1022 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                       (-1.0)*DQDT(1,ILAY1,ILAY2)/
     &                       (TR(ILAY1)**2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1021
 1022                CONTINUE
 1021             CONTINUE
C ................................................................
C        FIXED 5
C ................................................................
               ELSEIF (PID .EQ. 5) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1051 ILAY1=1,NLAYER
                     DO 1052 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                       DQDT(1,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1051
 1052                CONTINUE
 1051             CONTINUE
C ................................................................
C        FIXED 6
C ................................................................
               ELSEIF (PID .EQ. 6) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  CALL I_DQDT(2,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  RDMMY1=0
                  DO 1061 ILAY1=NLAYER,1,-1
                     IF (ILAY1 .EQ. NLAYER) THEN
                        DO 1062 ILAY2=ILAY1,NLAYER
                           DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                          SECANT*DQDT(2,ILAY1,ILAY2)/TR(ILAY1)
                           DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                           IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1061
 1062                   CONTINUE
                     ELSE
                        ! RDMMY1 is elements of Tz profile
                        RDMMY1=RDMMY1 + TR(ILAY1)*
     &                         PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
                        DO 1063 ILAY2=ILAY1,NLAYER
                           DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                          SECANT*( DQDT(2,ILAY1,ILAY2)/TR(ILAY1)-
     &                          DQDT(1,ILAY1,ILAY2)*RDMMY1/TR(ILAY1)**2)
                           DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                           IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1061
 1063                   CONTINUE
                     ENDIF
 1061             CONTINUE
C ................................................................
C        FIXED 7
C ................................................................
               ELSEIF (PID .EQ. 7) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1071 ILAY1=1,NLAYER
                     DO 1072 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*3.0*SECANT*
     &                        DQDT(1,ILAY1,ILAY2)*TR(ILAY1)**2
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1071
 1072                CONTINUE
 1071             CONTINUE
C ................................................................
C        FIXED 8
C ................................................................
               ELSEIF (PID .EQ. 8) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1081 ILAY1=1,NLAYER
                     DO 1082 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*2.0*
     &                        DQDT(1,ILAY1,ILAY2)*TR(ILAY1)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1081
 1082                CONTINUE
 1081             CONTINUE
C ................................................................
C        FIXED 9
C ................................................................
               ELSEIF (PID .EQ. 9) THEN
                  CALL I_DQDT(2,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1091 ILAY1=1,NLAYER
                     DO 1092 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*SECANT*
     &                        DQDT(2,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1091
 1092                CONTINUE
 1091             CONTINUE
                  
C ................................................................
C        FIXED 10
C ................................................................
               ELSEIF (PID .EQ. 10) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1101 ILAY1=1,NLAYER
                     DO 1102 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*2.0*SECANT*
     &                        DQDT(1,ILAY1,ILAY2)*TR(ILAY1)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1101
 1102                CONTINUE
 1101             CONTINUE
C ................................................................
C        FIXED 11
C ................................................................
               ELSEIF (PID .EQ. 11) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1111 ILAY1=1,NLAYER
                     DO 1112 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*SECANT*
     &                        DQDT(1,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1111
 1112                CONTINUE
 1111             CONTINUE
C ................................................................
C        FIXED 13
C ................................................................
               ELSEIF (PID .EQ. 13) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1131 ILAY1=1,NLAYER
                     DO 1132 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        (SECANT**1.5)*DQDT(1,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1131
 1132                CONTINUE
 1131             CONTINUE
C ................................................................
C        FIXED 14
C ................................................................
               ELSEIF (PID .EQ. 14) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1141 ILAY1=1,NLAYER
                     DO 1142 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*1.5*
     &                        (SECANT**1.5)*DSQRT(TR(ILAY1))*
     &                        DQDT(1,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1141
 1142                CONTINUE
 1141             CONTINUE
C ................................................................
C        FIXED 15
C ................................................................
               ELSEIF (PID .EQ. 15) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1151 ILAY1=1,NLAYER
                     DO 1152 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*(SECANT**2)*
     &                        DQDT(1,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1151
 1152                CONTINUE
 1151             CONTINUE
C ................................................................
C        FIXED 16
C ................................................................
               ELSEIF (PID .EQ. 16) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 1161 ILAY1=1,NLAYER
                     DO 1162 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*2.0*
     &                        (SECANT**2)*TR(ILAY1)*DQDT(1,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1161
 1162                CONTINUE
 1161             CONTINUE
C ................................................................
C        FIXED - UNDEFINED
C ................................................................
               ELSEIF (PID .GT. 16) THEN
                  WRITE (*,*) 'FATAL ERROR IN C_DKDT:'
                  WRITE (*,*) '((IMOL .EQ. 1) .AND. (PID .GT. 16))'
                  WRITE (*,*) 'PROGRAM ABORTED'
                  STOP
               ENDIF
C ................................................................
            ELSEIF (IMOL .EQ. 2) THEN
C ................................................................
C        PREDICTORS FOR WATER
C ................................................................
C        WATER 5
C ................................................................
               IF (PID .EQ. 5) THEN
                  CALL I_DQDT(6,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2051 ILAY1=1,NLAYER
                     DO 2052 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*SECANT*
     &                        WR(ILAY1)*DQDT(6,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2051
 2052                CONTINUE
 2051             CONTINUE
C ................................................................
C        WATER 9
C ................................................................
               ELSEIF (PID .EQ. 9) THEN
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2091 ILAY1=1,NLAYER
                     DO 2092 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*SECANT*
     &                        WR(ILAY1)*DQDT(3,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2091
 2092                CONTINUE
 2091             CONTINUE
C ................................................................
C        WATER 12
C ................................................................
               ELSEIF (PID .EQ. 12) THEN
                  CALL I_DQDT(6,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2121 ILAY1=1,NLAYER
                     DO 2122 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*2.0*SECANT*
     &                       WR(ILAY1)*TWZ2(ILAY1)*DQDT(6,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2121
 2122                CONTINUE
 2121             CONTINUE
C ................................................................
C        WATER 13
C ................................................................
               ELSEIF (PID .EQ. 13) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2131 ILAY1=1,NLAYER
                     DO 2132 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        (-1.0)*WR(ILAY1)*DQDT(1,ILAY1,ILAY2)/
     &                        (TR(ILAY1)**2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2131
 2132                CONTINUE
 2131             CONTINUE
C ................................................................
C        WATER 14
C ................................................................
               ELSEIF (PID .EQ. 14) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2141 ILAY1=1,NLAYER
                     DO 2142 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        (-1.0)*DQDT(1,ILAY1,ILAY2)/(TR(ILAY1)**2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2141
 2142                CONTINUE
 2141             CONTINUE
C ................................................................
C        WATER 15
C ................................................................
               ELSEIF (PID .EQ. 15) THEN
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2151 ILAY1=1,NLAYER
                     DO 2152 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                       DSQRT(SECANT*WR(ILAY1))*DQDT(3,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2151
 2152                CONTINUE
 2151             CONTINUE
C ................................................................
C        WATER 20
C ................................................................
               ELSEIF (PID .EQ. 20) THEN
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2201 ILAY1=1,NLAYER
                     DO 2202 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        DQDT(3,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2201
 2202                CONTINUE
 2201             CONTINUE
C ................................................................
C        WATER 21
C ................................................................
               ELSEIF (PID .EQ. 21) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2211 ILAY1=1,NLAYER
                     DO 2212 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        DQDT(1,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2211
 2212                CONTINUE
 2211             CONTINUE
C ................................................................
C        WATER 23
C ................................................................
               ELSEIF (PID .EQ. 23) THEN
!                 CALL MISSIT(RCODE,IMOL,PID)
!                 new: coded 9/24/98, T. Wehr
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2231 ILAY1=1,NLAYER
                     DO 2232 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        SECANT*WR(ILAY1)*DQDT(3,ILAY1,ILAY2)*
     &                        ABS(TD(ILAY1))
                        IF (ILAY1 .EQ. ILAY2) THEN
                           IF (TD(ILAY2) .GT. 0.0) THEN
                              DKADD=DKADD+COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                              SECANT*WR(ILAY1)*TD(ILAY1)
                           ELSEIF (TD(ILAY2) .LT. 0.0) THEN
                              DKADD=DKADD-COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                              SECANT*WR(ILAY1)*TD(ILAY1)
                           ENDIF
                        ENDIF
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2231
 2232                CONTINUE
 2231             CONTINUE
                  



C ................................................................
C        WATER 26
C ................................................................
               ELSEIF (PID .EQ. 26) THEN
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2261 ILAY1=1,NLAYER
                     DO 2262 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        ((SECANT*WR(ILAY1))**1.5)*
     &                        DQDT(3,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2261
 2262                CONTINUE
 2261             CONTINUE
C ................................................................
C        WATER 28
C ................................................................
               ELSEIF (PID .EQ. 28) THEN
                  CALL MISSIT(RCODE,IMOL,PID)
C ................................................................
C        WATER - UNDEFINED
C ................................................................
               ELSEIF (PID .GT. 33) THEN
                  WRITE (*,*) 'FATAL ERROR IN C_DKDT:'
                  WRITE (*,*) '((IMOL .EQ. 2) .AND. (PID .GT. 32))'
                  WRITE (*,*) 'PROGRAM ABORTED'
                  STOP
               ENDIF
C ................................................................
            ELSEIF (IMOL .EQ. 3) THEN
C ................................................................
C        PREDICTORS FOR OZONE
C ................................................................
C        OZONE 2
C ................................................................
               IF (PID .EQ. 2) THEN
                  CALL I_DQDT(1,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3021 ILAY1=1,NLAYER
                     DO 3022 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*(-1.0)*
     &                        OR(ILAY1)*DQDT(1,ILAY1,ILAY2)/
     &                        (TR(ILAY1)**2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 3021
 3022                CONTINUE
 3021             CONTINUE
C ................................................................
C        OZONE 5
C ................................................................
               ELSEIF (PID .EQ. 5) THEN
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3051 ILAY1=1,NLAYER
                     DO 3052 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        OR(ILAY1)*DQDT(3,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 3051
 3052                CONTINUE
 3051             CONTINUE
C ................................................................
C        OZONE 6
C ................................................................
               ELSEIF (PID .EQ. 6) THEN
                  CALL I_DQDT(11,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3061 ILAY1=1,NLAYER
                     DO 3062 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*SECANT*
     &                        OR(ILAY1)*DQDT(11,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 3061
 3062                CONTINUE
 3061             CONTINUE
C ................................................................
C        OZONE 9
C ................................................................
               ELSEIF (PID .EQ. 9) THEN
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3091 ILAY1=1,NLAYER
                     DO 3092 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        DSQRT(SECANT*OR(ILAY1))*
     &                        DQDT(3,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 3091
 3092                CONTINUE
 3091             CONTINUE
C ................................................................
C        OZONE 11
C ................................................................
               ELSEIF (PID .EQ. 11) THEN
                  CALL I_DQDT(9,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3111 ILAY1=1,NLAYER
                     DO 3112 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        DSQRT(SECANT*OR(ILAY1))*
     &                        DQDT(9,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 3111
 3112                CONTINUE
 3111             CONTINUE
C ................................................................
C        OZONE 12
C ................................................................
               ELSEIF (PID .EQ. 12) THEN
                  CALL I_DQDT(9,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3121 ILAY1=1,NLAYER
                     DO 3122 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        SECANT*OR(ILAY1)*DQDT(9,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 3121
 3122                CONTINUE
 3121             CONTINUE
C ................................................................
C        OZONE 16
C ................................................................
               ELSEIF (PID .EQ. 16) THEN
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3161 ILAY1=1,NLAYER
                     DO 3162 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        SECANT*OR(ILAY1)*DQDT(3,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 3161
 3162                CONTINUE
 3161             CONTINUE
C ................................................................
C        OZONE 19
C ................................................................
               ELSEIF (PID .EQ. 19) THEN
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3191 ILAY1=1,NLAYER
                     DO 3192 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        ((SECANT*OR(ILAY1))**1.5)*
     &                        DQDT(3,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 3191
 3192                CONTINUE
 3191             CONTINUE
C ................................................................
C        OZONE 21
C ................................................................
               ELSEIF (PID .EQ. 21) THEN
                  CALL I_DQDT(3,DQDT_I,DQDT,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3211 ILAY1=1,NLAYER
                     DO 3212 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        ((SECANT*OR(ILAY1))**2)*
     &                        DQDT(3,ILAY1,ILAY2)
                        DKDT_M(ILAY1,ILAY2)=DKDT_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 3211
 3212                CONTINUE
 3211             CONTINUE
C ................................................................
C        OZONE - UNDEFINED
C ................................................................
               ELSEIF (PID .GT. 24) THEN
                  WRITE (*,*) 'FATAL ERROR IN C_DKDT:'
                  WRITE (*,*) '((IMOL .EQ. 3) .AND. (PID .GT. 24))'
                  WRITE (*,*) 'PROGRAM ABORTED'
                  STOP
               ENDIF
C ................................................................
            ENDIF
 600     CONTINUE

         ! NOTE: there is no tuning since it is an additive to
         !       the absorption coefficient which will disappear
         !       by differentiating

         ! now add DKDT_M to DKDT:
         ! add only on layers where K in SUBROUTINE C_K is
         ! not smaller than zero
         DO 801 ILAY1=1,NLAYER
            IF (K_NEG(IIRF,IMOL,ILAY1) .EQ. 0) THEN
               DO 802 ILAY2=ILAY1,NLAYER
                  DKDT(ILAY1,ILAY2)=DKDT(ILAY1,ILAY2) +
     &                              DKDT_M(ILAY1,ILAY2)
                  IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 801
 802           CONTINUE
            ENDIF
 801     CONTINUE

 700  CONTINUE

      RETURN
      END
C     END OF SUBROUTINE C_DKDT


C ====================================================================

      SUBROUTINE C_DKDW(IIRF,PREDID,NCOEFF,COEFF,
     &                    DQDW_I,DQDW,DKDW,SECANT,
     &                    ATEMP,AFIXED,AWATER,AOZONE,
     &                    RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEG,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     calculate derivatives dk/dW for IIRF
C     (k=absorption koefficients, T=temperature)
C
C     INPUT:
C       IIRF            IRF counter = actual IRF (channel, or channel ID)
C       PREDID          IDs of the predictors to be used (see COMMENT)
C       NCOEFF          number of coefficients to be used (see COMMENT)
C       COEFF           coefficients to be used (see COMMENT)
C       SECANT          observation secant
C       ATEMP,AFIXED,AWATER,AOZONE       atmosphere profile          
C       RTEMP,RFIXED,RWATER,ROZONE,PRES  reference atmosphere profile
C       K_NEG           indecees if K has been disregarded in SUBROUTINE C_K
C                       (=1) because it was negative, or used (=0)
C       VERBOS          screen output on/off (1=on,0=off)
C
C     INPUT AND OUTPUT:
C       DQDW_I          Q-profiles initialization flags
C       DQDW            derivatives of Q-profiles with respect to W
C
C     OUTPUT:
C       DKDW            profile of derivatives dk/dW
C
C     CALLING ROUNTINE
C       C_DTAU
C
C     COMMENTS:
C       1. PREDID,NCOEFF,COEFF contain the information for ALL channels,
C          no matter whether all channels are used or not.
C
C     HISTORY:
C       written 12/09/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE CALLING VARIABLES
      INTEGER IIRF,
     &        NCOEFF(MAXIRF,MAXMOL),
     &        DQDW_I(MAX_QP),
     &        PREDID(MAXIRF,MAXMOL,MAXPRD+1),
     &        K_NEG(MAXIRF,MAXMOL,NLAYER),
     &        VERBOS
      REAL*8  COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER),
     &        DKDW(NLAYER,NLAYER),
     &        DQDW(MAX_QP,NLAYER,NLAYER),
     &        SECANT,
     &        ATEMP(NLAYER),AFIXED(NLAYER),AWATER(NLAYER),
     &        AOZONE(NLAYER),RTEMP(NLAYER),RFIXED(NLAYER),
     &        RWATER(NLAYER),ROZONE(NLAYER),PRES(NLAYER)

C     LOCAL VARIABLES
      INTEGER IMOL,                  ! molecule counter
     &        ILAY1,ILAY2,           ! layer counters
     &        IPRED,                 ! predictor counter
     &        PID,                   ! predictor ID
     &        RCODE                  ! subroutine identifier, used in MISSIT
      REAL*8  DKADD,                 ! dk/dW-profile increment
     &        DKDW_M(NLAYER,NLAYER), ! dk/dW-profile of a certain molecule
     &        WR(NLAYER),            ! profile AWATER/RWATER
     &        WZ(NLAYER),            ! profile W_{Z} (see latex doc)
     &        TD(NLAYER),            ! profile ATEMP-RTEMP
     &        TR(NLAYER),            ! profile ATEMP/RTEMP
     &        TWZ2(NLAYER),          ! profile T_{Wz2} (see latex doc) 
     &        RDMMY1,RDMMY2,RDMMY3,RDMMY4,RDMMY5

      PARAMETER(RCODE=2) ! SUBROUTINE C_DKDT=1, C_DKDW=2, C_DKDO=3

C     ===========================================================
C     BEGIN
C     ===========================================================

C     ........ CALCULATE PREDICTORS: ...........
C     initialize dk/dW-profiles with zeros
      DO 101 ILAY1=1,NLAYER
         DO 102 ILAY2=1,NLAYER
            DKDW(ILAY1,ILAY2)=0.0
 102     CONTINUE
 101  CONTINUE

C     calculate profiles WR,WZ,TWZ2,TD,TR
      RDMMY1=0.0
      RDMMY2=0.0
      RDMMY4=0.0
      RDMMY5=0.0
      DO 103 ILAY1=NLAYER,1,-1
         TR(ILAY1)=ATEMP(ILAY1)/RTEMP(ILAY1)
         WR(ILAY1)=AWATER(ILAY1)/RWATER(ILAY1)
         TD(ILAY1)=ATEMP(ILAY1)-RTEMP(ILAY1)
         IF (ILAY1 .EQ. NLAYER) THEN
            TWZ2(ILAY1)=(AWATER(ILAY1)*ATEMP(ILAY1))/
     &           (RWATER(ILAY1)*RTEMP(ILAY1))
            WZ(ILAY1)=AWATER(ILAY1)/RWATER(ILAY1)
         ELSE
            RDMMY3=PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
            RDMMY1=RDMMY1+RDMMY3*AWATER(ILAY1)*ATEMP(ILAY1)
            RDMMY2=RDMMY2+RDMMY3*RWATER(ILAY1)*RTEMP(ILAY1)
            RDMMY4=RDMMY4+RDMMY3*AWATER(ILAY1)
            RDMMY5=RDMMY5+RDMMY3*RWATER(ILAY1)
            TWZ2(ILAY1)=RDMMY1/RDMMY2
            WZ(ILAY1)=RDMMY4/RDMMY5
         ENDIF
 103  CONTINUE

      DO 700 IMOL=2,MAXMOL ! there are no d(FIXED)/dW 
         ! reset DKDW_M
         DO 701 ILAY1=1,NLAYER
            DO 702 ILAY2=ILAY1,NLAYER
               DKDW_M(ILAY1,ILAY2)=0.0
 702        CONTINUE
 701     CONTINUE



C        loop all required predictors
         DO 600 IPRED=1,PREDID(IIRF,IMOL,1) ! PREDID(IIRF,IMOL,1) is the total
                                            ! number of predictors for (IIRF,IMOL)
            ! PID is introduced for programmer's convenience:
            PID = PREDID(IIRF,IMOL,1+IPRED) 

C           write (*,678) PID,IMOL 
C678        format(' C_DKDW: PID,IMOL=',I3,I3) 
            IF (IMOL .EQ. 2) THEN
C ................................................................
C        PREDICTORS FOR WATER
C ................................................................
C        WATER 1
C ................................................................
               IF (PID .EQ. 1) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2011 ILAY1=1,NLAYER
                     DO 2012 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*SECANT*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2011
 2012                CONTINUE
 2011             CONTINUE
C ................................................................
C        WATER 2
C ................................................................
               ELSEIF (PID .EQ. 2) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2021 ILAY1=1,NLAYER
                     DO 2022 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*0.5*
     &                        DSQRT(SECANT/WR(ILAY1))*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2021
 2022                CONTINUE
 2021             CONTINUE
C ................................................................
C        WATER 3
C ................................................................
               ELSEIF (PID .EQ. 3) THEN
                  CALL I_DQDW(5,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2031 ILAY1=1,NLAYER
                     DO 2032 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*0.5*
     &                        DSQRT(SECANT/WZ(ILAY1))*
     &                        DQDW(5,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2031
 2032                CONTINUE
 2031             CONTINUE
C ................................................................
C        WATER 4
C ................................................................
               ELSEIF (PID .EQ. 4) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2041 ILAY1=1,NLAYER
                     DO 2042 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*2.0*
     &                        (SECANT**2)*WR(ILAY1)*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2041
 2042                CONTINUE
 2041             CONTINUE
C ................................................................
C        WATER 5
C ................................................................
               ELSEIF (PID .EQ. 5) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  CALL I_DQDW(6,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2051 ILAY1=1,NLAYER
                     DO 2052 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*SECANT*(
     &                        WR(ILAY1)*DQDW(6,ILAY1,ILAY2) +
     &                        TWZ2(ILAY1)*DQDW(4,ILAY1,ILAY2))
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2051
 2052                CONTINUE
 2051             CONTINUE
C ................................................................
C        WATER 6
C ................................................................
               ELSEIF (PID .EQ. 6) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2061 ILAY1=1,NLAYER
                     DO 2062 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*3.0*
     &                        (SECANT**3)*(WR(ILAY1)**2)*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2061
 2062                CONTINUE
 2061             CONTINUE
C ................................................................
C        WATER 7
C ................................................................
               ELSEIF (PID .EQ. 7) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  CALL I_DQDW(5,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2071 ILAY1=1,NLAYER
                     DO 2072 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        SECANT*WR(ILAY1)*(
     &                        2.0*DQDW(4,ILAY1,ILAY2)/WZ(ILAY1) -
     &                        WR(ILAY1)*DQDW(5,ILAY1,ILAY2)/
     &                        (WZ(ILAY1)**2) )
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2071
 2072                CONTINUE
 2071             CONTINUE
C ................................................................
C        WATER 8
C ................................................................
               ELSEIF (PID .EQ. 8) THEN
                  CALL I_DQDW(5,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2081 ILAY1=1,NLAYER
                     DO 2082 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*2.0*
     &                        (SECANT**2)*WZ(ILAY1)*
     &                        DQDW(5,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2081
 2082                CONTINUE
 2081             CONTINUE
C ................................................................
C        WATER 9
C ................................................................
               ELSEIF (PID .EQ. 9) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2091 ILAY1=1,NLAYER
                     DO 2092 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        SECANT*TD(ILAY1)*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2091
 2092                CONTINUE
 2091             CONTINUE
C ................................................................
C        WATER 10
C ................................................................
               ELSEIF (PID .EQ. 10) THEN
!                 CALL MISSIT(RCODE,IMOL,PID)
!                 new: coded 9/24/98, T. Wehr
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2101 ILAY1=1,NLAYER
                     DO 2102 ILAY2=1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        0.25 * (SECANT**0.25) *
     &                        (WR(ILAY1)**(-0.75)) *
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2101
 2102                CONTINUE
 2101             CONTINUE
C ................................................................
C        WATER 11
C ................................................................
               ELSEIF (PID .EQ. 11) THEN
                  CALL I_DQDW(5,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2111 ILAY1=1,NLAYER
                     DO 2112 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        SECANT*DQDW(5,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2111
 2112                CONTINUE
 2111             CONTINUE
C ................................................................
C        WATER 12
C ................................................................
               ELSEIF (PID .EQ. 12) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  CALL I_DQDW(6,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2121 ILAY1=1,NLAYER
                     DO 2122 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*SECANT*(
     &                        (TWZ2(ILAY1)**2)*DQDW(4,ILAY1,ILAY2)+
     &                        2.0*WR(ILAY1)*TWZ2(ILAY1)*
     &                        DQDW(6,ILAY1,ILAY2) )
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2121
 2122                CONTINUE
 2121             CONTINUE
C ................................................................
C        WATER 13
C ................................................................
               ELSEIF (PID .EQ. 13) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2131 ILAY1=1,NLAYER
                     DO 2132 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        DQDW(4,ILAY1,ILAY2)/TR(ILAY1)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2131
 2132                CONTINUE
 2131             CONTINUE
C ................................................................
C        WATER 15
C ................................................................
               ELSEIF (PID .EQ. 15) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2151 ILAY1=1,NLAYER
                     DO 2152 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*0.5*
     &                        TD(ILAY1)*DSQRT(SECANT/WR(ILAY1))*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2151
 2152                CONTINUE
 2151             CONTINUE
C ................................................................
C        WATER 16
C ................................................................
               ELSEIF (PID .EQ. 16) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  CALL I_DQDW(5,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2161 ILAY1=1,NLAYER
                     DO 2162 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*SECANT*(
     &                        1.5*DSQRT(WR(ILAY1))*DQDW(4,ILAY1,ILAY2)/
     &                        WZ(ILAY1)-(WR(ILAY1)**1.5)*
     &                        DQDW(5,ILAY1,ILAY2)/(WZ(ILAY1)**2) )
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2161
 2162                CONTINUE
 2161             CONTINUE
C ................................................................
C        WATER 17
C ................................................................
               ELSEIF (PID .EQ. 17) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  CALL I_DQDW(5,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2171 ILAY1=1,NLAYER
                     DO 2172 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*(SECANT**2)*(
     &                        3.0*(WR(ILAY1)**2)*DQDW(4,ILAY1,ILAY2)/
     &                        WZ(ILAY1)-(WR(ILAY1)**3)*
     &                        DQDW(5,ILAY1,ILAY2)/(WZ(ILAY1)**2) )
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2171
 2172                CONTINUE
 2171             CONTINUE
C ................................................................
C        WATER 18
C ................................................................
               ELSEIF (PID .EQ. 18) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  CALL I_DQDW(5,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2181 ILAY1=1,NLAYER
                     DO 2182 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        SECANT*TD(ILAY1)* (
     &                        2.0*WR(ILAY1)*DQDW(4,ILAY1,ILAY2)/
     &                        WZ(ILAY1)-(WR(ILAY1)**2)*
     &                        DQDW(5,ILAY1,ILAY2)/(WZ(ILAY1)**2) )
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2181
 2182                CONTINUE
 2181             CONTINUE
C ................................................................
C        WATER 19
C ................................................................
               ELSEIF (PID .EQ. 19) THEN
                  CALL I_DQDW(5,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2191 ILAY1=1,NLAYER
                     DO 2192 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        DQDW(5,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2191
 2192                CONTINUE
 2191             CONTINUE
C ................................................................
C        WATER 23
C ................................................................
               ELSEIF (PID .EQ. 23) THEN
!                 CALL MISSIT(RCODE,IMOL,PID)
!                 new: coded 9/24/98, T. Wehr
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2231 ILAY1=1,NLAYER
                     DO 2232 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        SECANT*TD(ILAY1)*ABS(TD(ILAY1))*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2231
 2232                CONTINUE
 2231             CONTINUE

C ................................................................
C        WATER 24
C ................................................................
               ELSEIF (PID .EQ. 24) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2241 ILAY1=1,NLAYER
                     DO 2242 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*1.25*
     &                        (SECANT**1.25)*(WR(ILAY1)**0.25)*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2241
 2242                CONTINUE
 2241             CONTINUE
C ................................................................
C        WATER 25
C ................................................................
               ELSEIF (PID .EQ. 25) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2251 ILAY1=1,NLAYER
                     DO 2252 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*1.5*
     &                        (SECANT**1.5)*(WR(ILAY1)**0.5)*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2251
 2252                CONTINUE
 2251             CONTINUE
C ................................................................
C        WATER 26
C ................................................................
               ELSEIF (PID .EQ. 26) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2261 ILAY1=1,NLAYER
                     DO 2262 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*1.5*
     &                        (SECANT**1.5)*(WR(ILAY1)**0.5)*
     &                        TD(ILAY1)*DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2261
 2262                CONTINUE
 2261             CONTINUE
C ................................................................
C        WATER 27
C ................................................................
               ELSEIF (PID .EQ. 27) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2271 ILAY1=1,NLAYER
                     DO 2272 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*2.5*
     &                        (SECANT**2.5)*(WR(ILAY1)**1.5)*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2271
 2272                CONTINUE
 2271             CONTINUE
C ................................................................
C        WATER 28
C ................................................................
               ELSEIF (PID .EQ. 28) THEN
                  CALL MISSIT(RCODE,IMOL,PID)
C ................................................................
C        WATER 29
C ................................................................
               ELSEIF (PID .EQ. 29) THEN
                  CALL MISSIT(RCODE,IMOL,PID)
C ................................................................
C        WATER 30
C ................................................................
               ELSEIF (PID .EQ. 30) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  CALL I_DQDW(5,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2301 ILAY1=1,NLAYER
                     DO 2302 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        2*(SECANT**2)* (
     &                        WR(ILAY1)*DQDW(4,ILAY1,ILAY2)/
     &                        (WZ(ILAY1)**2) - (WR(ILAY1)**2)*
     &                        DQDW(5,ILAY1,ILAY2)/(WZ(ILAY1)**3) )
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2301
 2302                CONTINUE
 2301             CONTINUE
C ................................................................
C        WATER 31
C ................................................................
               ELSEIF (PID .EQ. 31) THEN
!                 CALL MISSIT(RCODE,IMOL,PID)
!                 new: coded 9/24/98, T. Wehr
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2311 ILAY1=1,NLAYER
                     DO 2312 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        DLOG(1.1*SECANT)*DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2311
 2312                CONTINUE
 2311             CONTINUE
C ................................................................
C        WATER 32
C ................................................................
               ELSEIF (PID .EQ. 32) THEN
                  CALL MISSIT(RCODE,IMOL,PID)
C ................................................................
C        WATER 33
C ................................................................
               ELSEIF (PID .EQ. 33) THEN
!                 CALL MISSIT(RCODE,IMOL,PID)
!                 new: coded 9/24/98, T. Wehr
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 2331 ILAY1=1,NLAYER
                     DO 2332 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        DEXP(SECANT/2.0)*DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2331
 2332                CONTINUE
 2331             CONTINUE
C ................................................................
C        WATER - UNDEFINED
C ................................................................
               ELSEIF (PID .GT. 32) THEN
                  WRITE (*,*) 'FATAL ERROR IN C_DKDW:'
                  WRITE (*,*) '((IMOL .EQ. 2) .AND. (PID .GT. 32))'
                  WRITE (*,*) 'PROGRAM ABORTED'
                  STOP
               ENDIF
C ................................................................
            ELSEIF (IMOL .EQ. 3) THEN
C ................................................................
C        PREDICTORS FOR OZONE
C ................................................................
C        OZONE 13
C ................................................................
               IF (PID .EQ. 13) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3131 ILAY1=1,NLAYER
                     DO 3132 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        DQDW(4,ILAY1,ILAY2)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 3131
 3132                CONTINUE
 3131             CONTINUE
C ................................................................
C        OZONE 15
C ................................................................
               ELSEIF (PID .EQ. 15) THEN
                  CALL I_DQDW(4,DQDW_I,DQDW,ATEMP,AFIXED,AWATER,
     &                 AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
                  DO 3151 ILAY1=1,NLAYER
                     DO 3152 ILAY2=ILAY1,NLAYER
                        DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                        (SECANT**2)*DQDW(4,ILAY1,ILAY2)*
     &                        AOZONE(ILAY1)/ROZONE(ILAY1)
                        DKDW_M(ILAY1,ILAY2)=DKDW_M(ILAY1,ILAY2)+DKADD
                        IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 3151
 3152                CONTINUE
 3151             CONTINUE
C ................................................................
C        OZONE - UNDEFINED
C ................................................................
               ELSEIF (PID .GT. 24) THEN
                  WRITE (*,*) 'FATAL ERROR IN C_DKDW:'
                  WRITE (*,*) '((IMOL .EQ. 3) .AND. (PID .GT. 24))'
                  WRITE (*,*) 'PROGRAM ABORTED'
                  STOP
               ENDIF
C ................................................................
            ENDIF
 600     CONTINUE

         ! NOTE: there is no tuning since it is an additive to
         !       the absorption coefficient which will disappear
         !       by differentiating

         ! now add DKDW_M to DKDW:
         ! add only on layers where K in SUBROUTINE C_K is
         ! not smaller than zero
         DO 801 ILAY1=1,NLAYER
            IF (K_NEG(IIRF,IMOL,ILAY1) .EQ. 0) THEN
               DO 802 ILAY2=ILAY1,NLAYER
                  DKDW(ILAY1,ILAY2)=DKDW(ILAY1,ILAY2) +
     &                              DKDW_M(ILAY1,ILAY2)
                  IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 801
 802           CONTINUE
            ENDIF
 801     CONTINUE

 700  CONTINUE

      RETURN
      END
C     END OF SUBROUTINE C_DKDW


C ====================================================================

      SUBROUTINE C_DKDO(IIRF,PREDID,NCOEFF,COEFF,
     &                    DQDO_I,DQDO,DKDO,SECANT,
     &                    ATEMP,AFIXED,AWATER,AOZONE,
     &                    RTEMP,RFIXED,RWATER,ROZONE,PRES,K_NEG,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     calculate derivatives dk/dO for IIRF
C     (k=absorption koefficients, T=temperature)
C
C     INPUT:
C       IIRF            IRF counter = actual IRF (channel, or channel ID)
C       PREDID          IDs of the predictors to be used (see COMMENT)
C       NCOEFF          number of coefficients to be used (see COMMENT)
C       COEFF           coefficients to be used (see COMMENT)
C       SECANT          observation secant
C       ATEMP,AFIXED,AWATER,AOZONE       atmosphere profile          
C       RTEMP,RFIXED,RWATER,ROZONE,PRES  reference atmosphere profile
C       K_NEG           indecees if K has been disregarded in SUBROUTINE C_K
C                       (=1) because it was negative, or used (=0)
C       VERBOS          screen output on/off (1=on,0=off)
C
C     INPUT AND OUTPUT:
C       DQDO_I          Q-profiles initialization flags
C       DQDO            derivatives of Q-profiles with respect to O
C
C     OUTPUT:
C       DKDO            profile of derivatives dk/dO
C
C     CALLING ROUNTINE
C       C_DTAU
C
C     COMMENTS:
C       1. PREDID,NCOEFF,COEFF contain the information for ALL channels,
C          no matter whether all channels are used or not.
C
C     HISTORY:
C       written 12/09/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE CALLING VARIABLES
      INTEGER IIRF,
     &        NCOEFF(MAXIRF,MAXMOL),
     &        DQDO_I(MAX_QP),
     &        PREDID(MAXIRF,MAXMOL,MAXPRD+1),
     &        K_NEG(MAXIRF,MAXMOL,NLAYER),
     &        VERBOS
      REAL*8  COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER),
     &        DKDO(NLAYER,NLAYER),
     &        DQDO(MAX_QP,NLAYER,NLAYER),
     &        SECANT,
     &        ATEMP(NLAYER),AFIXED(NLAYER),AWATER(NLAYER),
     &        AOZONE(NLAYER),RTEMP(NLAYER),RFIXED(NLAYER),
     &        RWATER(NLAYER),ROZONE(NLAYER),PRES(NLAYER)

C     LOCAL VARIABLES
      INTEGER IMOL,                  ! molecule counter
     &        ILAY1,ILAY2,           ! layer counters
     &        IPRED,                 ! predictor counter
     &        PID,                   ! predictor ID
     &        RCODE                  ! subroutine identifier, used in MISSIT
      REAL*8  DKADD,                 ! dk/dO-profile increment
     &        DKDO_M(NLAYER,NLAYER), ! dk/dO-profile of a certain molecule
     &        TR(NLAYER),            ! profile ... (see latex doc)
     &        TD(NLAYER),            ! profile ... (see latex doc)
     &        OR(NLAYER),            ! profile ... (see latex doc)
     &        OZ(NLAYER),            ! profile ... (see latex doc)
     &        TOZ(NLAYER),           ! profile ... (see latex doc)
     &        OZ2(NLAYER),           ! profile ... (see latex doc)
     &        TOZ2(NLAYER),          ! profile ... (see latex doc)
     &        RDMMY1,RDMMY2,RDMMY3,RDMMY4,RDMMY5,RDMMY6,SPRES

      PARAMETER(RCODE=3) ! SUBROUTINE C_DKDT=1, C_DKDW=2, C_DKDO=3

C     ===========================================================
C     BEGIN
C     ===========================================================

C     ........ CALCULATE PREDICTORS: ...........
C     initialize dk/dO-profiles with zeros
      DO 101 ILAY1=1,NLAYER
         DO 102 ILAY2=1,NLAYER
            DKDO(ILAY1,ILAY2)=0.0
 102     CONTINUE
 101  CONTINUE

C     calculate profiles OR,TR,OZ2,TD,OZ,TOZ
      RDMMY1=0.0
      RDMMY2=0.0
      RDMMY3=0.0
      RDMMY4=0.0
      RDMMY5=0.0
      RDMMY6=0.0
      DO 103 ILAY1=NLAYER,1,-1
         OR(ILAY1)=AOZONE(ILAY1)/ROZONE(ILAY1)
         TR(ILAY1)=ATEMP(ILAY1)/RTEMP(ILAY1)
         TD(ILAY1)=ATEMP(ILAY1)-RTEMP(ILAY1)
         IF (ILAY1 .EQ. NLAYER) THEN
            OZ(ILAY1)=0.0
            TOZ(ILAY1)=0.0
            OZ2(ILAY1)=AOZONE(ILAY1)/ROZONE(ILAY1)
            TOZ2(ILAY1)=(AOZONE(ILAY1)*ATEMP(ILAY1))/
     &                  (ROZONE(ILAY1)*RTEMP(ILAY1))
         ELSE
            SPRES=PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
            RDMMY1=RDMMY1+SPRES*OR(ILAY1+1)
            OZ(ILAY1)=RDMMY1
            RDMMY2=RDMMY2+SPRES*OR(ILAY1+1)*TR(ILAY1+1)
            TOZ(ILAY1)=RDMMY2
            RDMMY3=RDMMY3+SPRES*AOZONE(ILAY1)
            RDMMY4=RDMMY4+SPRES*ROZONE(ILAY1)
            OZ2(ILAY1)=RDMMY3/RDMMY4
            RDMMY5=RDMMY5+SPRES*AOZONE(ILAY1)*ATEMP(ILAY1)
            RDMMY6=RDMMY6+SPRES*ROZONE(ILAY1)*RTEMP(ILAY1)
            TOZ2(ILAY1)=RDMMY5/RDMMY6
         ENDIF
 103  CONTINUE

      IMOL=3 ! there is no d(FIXED)/dO and no d(WATER)/dO

      DO 701 ILAY1=1,NLAYER
         DO 702 ILAY2=ILAY1,NLAYER
            DKDO_M(ILAY1,ILAY2)=0.0
 702     CONTINUE
 701  CONTINUE

C     loop all required predictors
      DO 600 IPRED=1,PREDID(IIRF,IMOL,1) ! PREDID(IIRF,IMOL,1) is the total
                                         ! number of predictors for (IIRF,IMOL)

         ! PID is introduced for programmer's convenience:
         PID = PREDID(IIRF,IMOL,1+IPRED) 
C ................................................................
C        PREDICTORS FOR OZONE
C ................................................................
C        OZONE 1
C ................................................................
         IF (PID .EQ. 1) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3011 ILAY1=1,NLAYER
               DO 3012 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  0.5*DSQRT(SECANT/OR(ILAY1))*
     &                  DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3011
 3012          CONTINUE
 3011       CONTINUE
C ................................................................
C        OZONE 2
C ................................................................
         ELSEIF (PID .EQ. 2) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3021 ILAY1=1,NLAYER
               DO 3022 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  DQDO(7,ILAY1,ILAY2)/TR(ILAY1)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3021
 3022          CONTINUE
 3021       CONTINUE
C ................................................................
C        OZONE 3
C ................................................................
         ELSEIF (PID .EQ. 3) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3031 ILAY1=1,NLAYER
               DO 3032 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  SECANT*DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3031
 3032          CONTINUE
 3031       CONTINUE
C ................................................................
C        OZONE 4
C ................................................................
         ELSEIF (PID .EQ. 4) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            CALL I_DQDO(10,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3041 ILAY1=1,NLAYER
               DO 3042 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  DSQRT(SECANT)* (
     &                  1.5*DSQRT(OR(ILAY1))*
     &                  DQDO(7,ILAY1,ILAY2)/OZ2(ILAY1) -
     &                  (OR(ILAY1)**1.5)*DQDO(10,ILAY1,ILAY2)/
     &                  (OZ2(ILAY1)**2) )
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3041
 3042          CONTINUE
 3041       CONTINUE
C ................................................................
C        OZONE 5
C ................................................................
         ELSEIF (PID .EQ. 5) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3051 ILAY1=1,NLAYER
               DO 3052 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  TD(ILAY1)*DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3051
 3052          CONTINUE
 3051       CONTINUE
C ................................................................
C        OZONE 6
C ................................................................
         ELSEIF (PID .EQ. 6) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            CALL I_DQDO(11,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3061 ILAY1=1,NLAYER
               DO 3062 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  SECANT* (
     &                  TOZ2(ILAY1)*DQDO(7,ILAY1,ILAY2) +
     &                  OR(ILAY1)*DQDO(11,ILAY1,ILAY2) )
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3061
 3062          CONTINUE
 3061       CONTINUE
C ................................................................
C        OZONE 7
C ................................................................
         ELSEIF (PID .EQ. 7) THEN
!           new: coded 9/24/98, T. Wehr
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3071 ILAY1=1,NLAYER
               DO 3072 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  2.0*SECANT*SECANT*OR(ILAY1)*DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3071
 3072          CONTINUE
 3071       CONTINUE
C ................................................................
C        OZONE 8
C ................................................................
         ELSEIF (PID .EQ. 8) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            CALL I_DQDO(8,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3081 ILAY1=1,(NLAYER-1)
               DO 3082 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  (SECANT**1.5)* (
     &                  DSQRT(OZ(ILAY1))*DQDO(7,ILAY1,ILAY2) +
     &                  0.5*OR(ILAY1)*DQDO(8,ILAY1,ILAY2)/
     &                  DSQRT(OZ(ILAY1)) ) 
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3081
 3082          CONTINUE
 3081       CONTINUE
C ................................................................
C        OZONE 9
C ................................................................
         ELSEIF (PID .EQ. 9) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3091 ILAY1=1,NLAYER
               DO 3092 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  0.5*TD(ILAY1)*DSQRT(SECANT/OR(ILAY1))*
     &                  DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3091
 3092          CONTINUE
 3091       CONTINUE            
C ................................................................
C        OZONE 10
C ................................................................
         ELSEIF (PID .EQ. 10) THEN
            CALL I_DQDO(8,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3101 ILAY1=1,NLAYER
               DO 3102 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  SECANT*DQDO(8,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3101
 3102          CONTINUE
 3101       CONTINUE            
C ................................................................
C        OZONE 11
C ................................................................
         ELSEIF (PID .EQ. 11) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            CALL I_DQDO(9,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3111 ILAY1=1,NLAYER
               DO 3112 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  DSQRT(SECANT)* (
     &                  DSQRT(OR(ILAY1))*DQDO(9,ILAY1,ILAY2) +
     &                  0.5*TOZ(ILAY1)*DQDO(7,ILAY1,ILAY2)/
     &                  DSQRT(OR(ILAY1)) )
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3111
 3112          CONTINUE
 3111       CONTINUE            
C ................................................................
C        OZONE 12
C ................................................................
         ELSEIF (PID .EQ. 12) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            CALL I_DQDO(9,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3121 ILAY1=1,NLAYER
               DO 3122 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  SECANT* (
     &                  TOZ(ILAY1)*DQDO(7,ILAY1,ILAY2) +
     &                  OR(ILAY1)*DQDO(9,ILAY1,ILAY2) )
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3121
 3122          CONTINUE
 3121       CONTINUE            
C ................................................................
C        OZONE 14
C ................................................................
         ELSEIF (PID .EQ. 14) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            CALL I_DQDO(10,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3141 ILAY1=1,NLAYER
               DO 3142 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                 SECANT*(
     &                 2.0*OR(ILAY1)*DQDO(7,ILAY1,ILAY2)/OZ2(ILAY1)-
     &                 ((OR(ILAY1)/OZ2(ILAY1))**2)*DQDO(10,ILAY1,ILAY2))
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3141
 3142          CONTINUE
 3141       CONTINUE            
C ................................................................
C        OZONE 15
C ................................................................
         ELSEIF (PID .EQ. 15) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3151 ILAY1=1,NLAYER
               DO 3152 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  (SECANT**2)*DQDO(7,ILAY1,ILAY2)*
     &                  AWATER(ILAY1)/RWATER(ILAY1)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3151
 3152          CONTINUE
 3151       CONTINUE            
C ................................................................
C        OZONE 16
C ................................................................
         ELSEIF (PID .EQ. 16) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3161 ILAY1=1,NLAYER
               DO 3162 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  SECANT*TD(ILAY1)*DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3161
 3162          CONTINUE
 3161       CONTINUE            
C ................................................................
C        OZONE 17
C ................................................................
         ELSEIF (PID .EQ. 17) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3171 ILAY1=1,NLAYER
               DO 3172 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  1.25*(SECANT**1.25)*(OR(ILAY1)**0.25)*
     &                  DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3171
 3172          CONTINUE
 3171       CONTINUE            
C ................................................................
C        OZONE 18
C ................................................................
         ELSEIF (PID .EQ. 18) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3181 ILAY1=1,NLAYER
               DO 3182 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  1.5*(SECANT**1.5)*DSQRT(OR(ILAY1))*
     &                  DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3181
 3182          CONTINUE
 3181       CONTINUE            
C ................................................................
C        OZONE 19
C ................................................................
         ELSEIF (PID .EQ. 19) THEN
!           new: coded 9/24/98, T. Wehr
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3191 ILAY1=1,NLAYER
               DO 3192 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  1.5*(SECANT**1.5)*DSQRT(OR(ILAY1))*
     &                  TD(ILAY1)*DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3191
 3192          CONTINUE
 3191       CONTINUE            
C ................................................................
C        OZONE 20
C ................................................................
         ELSEIF (PID .EQ. 20) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3201 ILAY1=1,NLAYER
               DO 3202 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  2.5*(SECANT**2.5)*(OR(ILAY1)**1.5)*
     &                  DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3201
 3202          CONTINUE
 3201       CONTINUE            
C ................................................................
C        OZONE 21
C ................................................................
         ELSEIF (PID .EQ. 21) THEN
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3211 ILAY1=1,NLAYER
               DO 3212 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  2.0*(SECANT**2)*TD(ILAY1)*OR(ILAY1)*
     &                  DQDO(7,ILAY1,ILAY2)
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3211
 3212          CONTINUE
 3211       CONTINUE            
C ................................................................
C        OZONE 22
C ................................................................
         ELSEIF (PID .EQ. 22) THEN
            CALL MISSIT(RCODE,IMOL,PID)
C ................................................................
C        OZONE 23
C ................................................................
         ELSEIF (PID .EQ. 23) THEN
!           new: coded 9/24/98, T. Wehr
            CALL I_DQDO(7,DQDO_I,DQDO,ATEMP,AFIXED,AWATER,
     &           AOZONE,RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
            DO 3231 ILAY1=1,NLAYER
               DO 3232 ILAY2=ILAY1,NLAYER
                  DKADD=COEFF(IIRF,IMOL,IPRED,ILAY1)*
     &                  DQDO(7,ILAY1,ILAY2)/
     &                  (1.1*SECANT*OR(ILAY1))
                  DKDO_M(ILAY1,ILAY2)=DKDO_M(ILAY1,ILAY2)+DKADD
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 3231
 3232          CONTINUE
 3231       CONTINUE
C ................................................................
C        OZONE 24
C ................................................................
         ELSEIF (PID .EQ. 24) THEN
            CALL MISSIT(RCODE,IMOL,PID)
C ................................................................
C        OZONE - UNDEFINED
C ................................................................
         ELSEIF (PID .GT. 24) THEN
            WRITE (*,*) 'FATAL ERROR IN C_DKDO:'
            WRITE (*,*) '((IMOL .EQ. 3) .AND. (PID .GT. 24))'
            WRITE (*,*) 'PROGRAM ABORTED'
            STOP
         ENDIF

 600  CONTINUE

      ! NOTE: there is no tuning since it is an additive to
      !       the absorption coefficient which will disappear
      !       by differentiating

      ! now add DKDO_M to DKDO:
      ! add only on layers where K in SUBROUTINE C_K is
      ! not smaller than zero
      DO 801 ILAY1=1,NLAYER
         IF (K_NEG(IIRF,IMOL,ILAY1) .EQ. 0) THEN
            DO 802 ILAY2=ILAY1,NLAYER
               DKDO(ILAY1,ILAY2)=DKDO(ILAY1,ILAY2) +
     &                           DKDO_M(ILAY1,ILAY2)
               IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 801
 802        CONTINUE
         ENDIF
 801  CONTINUE

      RETURN
      END
C     END OF SUBROUTINE C_DKDO


      SUBROUTINE MISSIT(RCODE,IMOL,PID)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Notifies missing code for calculating a derivative and
C       stops program.
C
C     NOTE:
C       This routine reports incomplete code.
C
C     INPUT:
C       RCODE        code of the calling routine
C                    1:   C_DKDT
C                    2:   C_DKDW
C                    3:   C_DKDO
C       IMOL         molecule code
C                    1:   F
C                    2:   W
C                    3:   O
C       PID          predictor ID
C
C     OUTPUT:
C       no output
C
C     HISTORY:
C       written 12/03/97 by Tobias Wehr
C
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
C     SUBROUTINE CALLING PARAMETERS
      INTEGER RCODE,IMOL,PID
C     LOCAL PARAMETERS
      CHARACTER*17  DXNAME
      CHARACTER*5   MLNAME
      INTEGER       NOTED(MAXMOL,MAXPRD,3),   ! ONLY FOR CODE DEVELOPEMENT
     &              IMOL2,IPRD2,IRCODE,FIRSTC ! ONLY FOR CODE DEVELOPEMENT
      DATA          FIRSTC /1/                ! ONLY FOR CODE DEVELOPEMENT
      SAVE          NOTED,FIRSTC              ! ONLY FOR CODE DEVELOPEMENT

C     BEGIN

      ! THIS BLOCK ONLY FOR CODE DEVELOPMENT
      IF (FIRSTC .EQ. 1) THEN        
         FIRSTC=0                    
         DO 101 IMOL2=1,MAXMOL       
            DO 102 IPRD2=1,MAXPRD    
               DO 103 IRCODE=1,3
                  NOTED(IMOL2,IPRD2,IRCODE)=0
 103           CONTINUE
 102        CONTINUE                 
 101     CONTINUE                    
      ENDIF                          
            
      IF (RCODE .EQ. 1) THEN
         DXNAME='TEMP. DERIVATIVES'
      ELSEIF (RCODE .EQ. 2) THEN
         DXNAME='WATER DERIVATIVES'
      ELSEIF (RCODE .EQ. 3) THEN
         DXNAME='OZONE DERIVATIVES'
      ELSE
         WRITE (*,*) 'FATAL ERROR IN SUBROUTINE MISSIT:'
         WRITE (*,*) 'INVALID RCODE - PROGRAM ABORTED'
         STOP
      ENDIF
      
      IF (IMOL .EQ. 1) THEN
         MLNAME='FIXED'
      ELSEIF (IMOL .EQ. 2) THEN
         MLNAME='WATER'
      ELSEIF (IMOL .EQ. 3) THEN
         MLNAME='OZONE'
      ELSE
         WRITE (*,*) 'FATAL ERROR IN SUBROUTINE MISSIT'
         WRITE (*,*) 'INVALID IMOL - PROGRAM ABORTED'
         STOP
      ENDIF

      ! THIS BLOCK IS FOR THE FINAL CODE
      WRITE (*,*) 'FATAL ERROR: INCOMPLETE PROGRAM CODE'
      WRITE (*,10) DXNAME,MLNAME,PID
 10   FORMAT(' THE ',A17,' ARE NOT CODED FOR ',A5,' PID ',I2)
      WRITE (*,*) 'PROGRAM ABORTED'
      STOP

      ! THIS BLOCK ONLY FOR CODE DEVELOPMENT
c     IF (NOTED(IMOL,PID,RCODE) .EQ. 0) THEN
c        !IF (MLNAME .EQ. 'WATER') THEN
c        IF (RCODE .EQ. 3) THEN
c           WRITE (*,12) DXNAME,MLNAME,PID
c12         FORMAT(' THE ',A17,' ARE NOT CODED FOR ',A5,' PID ',I2)
c           NOTED(IMOL,PID,RCODE)=1
c        ENDIF
c     ENDIF
c     RETURN ! this line belongs to the "develpment"-code

      END
C     END OF SUBROUTINE MISSIT



C ====================================================================

      SUBROUTINE I_DQDT(QPID,DQDT_I,DQDT,
     &                  ATEMP,AFIXED,AWATER,AOZONE,
     &                  RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     initialize dQ/dT-profile number QPID
C
C     INPUT:
C       QPID            ID of Q-profile which shall be initialized
C       ATEMP,AFIXED,AWATER,AOZONE       atmosphere profile          
C       RTEMP,RFIXED,RWATER,ROZONE,PRES  reference atmosphere profile
C       VERBOS          screen output on/off (1=on,0=off)
C
C     INPUT AND OUTPUT:
C       DQDT_I          dQ/dT-profiles initialization flags
C       DQDT            dQ/dT-profiles 
C
C     CALLING ROUNTINE
C       C_DKDT
C
C     HISTORY:
C       written 12/5/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE CALLING VARIABLES
      INTEGER     QPID,DQDT_I(MAX_QP),VERBOS
      REAL*8      DQDT(MAX_QP,NLAYER,NLAYER), ! DQDT(?,i,j)=dQ(i)/dT(j)
     &            ATEMP(NLAYER),AFIXED(NLAYER),AWATER(NLAYER),
     &            AOZONE(NLAYER),RTEMP(NLAYER),RFIXED(NLAYER),
     &            RWATER(NLAYER),ROZONE(NLAYER),PRES(NLAYER)
C     LOCAL VARIABLES
      INTEGER     ILAY1,ILAY2
      REAL*8      RDMMY1,RDMMY2,RDMMY3

C     check if initialization flag is already set
      IF (DQDT_I(QPID) .EQ. 1) THEN
         RETURN
      ELSE
         IF (VERBOS .EQ. 1) THEN
            WRITE (*,20) QPID
 20         FORMAT(' initializing dQ/dT-profile ',I2)
         ENDIF
      ENDIF
      
C     calculate Q-profiles:
C     NOTE: ALL DQDT ARE ALREADY INITIALIZES AS ZEROS
C............................................................................
C     Q-PROFILE 1
C............................................................................
      IF (QPID .EQ. 1) THEN
         DO 1010 ILAY1=1,NLAYER
            DQDT(1,ILAY1,ILAY1)=1.0/RTEMP(ILAY1)
 1010    CONTINUE
         DQDT_I(1)=1
C............................................................................
C     Q-PROFILE 2
C............................................................................
      ELSEIF (QPID .EQ. 2) THEN
         DO 1021 ILAY1=1,(NLAYER-1)
            DO 1022 ILAY2=ILAY1,(NLAYER-1)
               DQDT(2,ILAY1,ILAY2)=PRES(ILAY2)*
     &              (PRES(ILAY2)-PRES(ILAY2+1))/RTEMP(ILAY2)
               IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1021
 1022       CONTINUE
 1021    CONTINUE
         DQDT_I(2)=1
C............................................................................
C     Q-PROFILE 3
C............................................................................
      ELSEIF (QPID .EQ. 3) THEN
         DO 1030 ILAY1=1,NLAYER
            DQDT(3,ILAY1,ILAY1)=1.0
 1030    CONTINUE
         DQDT_I(3)=1
C............................................................................
C     Q-PROFILE 6
C............................................................................
      ELSEIF (QPID .EQ. 6) THEN
         RDMMY1=0.0
         DO 1061 ILAY1=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
            RDMMY1=RDMMY1+RDMMY3*RWATER(ILAY1)*RTEMP(ILAY1)
            DO 1062 ILAY2=ILAY1,(NLAYER-1)
               RDMMY2=PRES(ILAY2)*(PRES(ILAY2)-PRES(ILAY2+1))
     &                *AWATER(ILAY2)
               DQDT(6,ILAY1,ILAY2)=RDMMY2/RDMMY1
               IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1061
 1062       CONTINUE
 1061    CONTINUE
         DQDT(6,NLAYER,NLAYER)=AWATER(NLAYER)/
     &                         (RWATER(NLAYER)*RTEMP(NLAYER))
         DQDT_I(6)=1
C............................................................................
C     Q-PROFILE 9
C............................................................................
      ELSEIF (QPID .EQ. 9) THEN
         DO 1091 ILAY1=1,NLAYER-1
            DO 1092 ILAY2=(ILAY1+1),NLAYER
               DQDT(9,ILAY1,ILAY2)=
     &              PRES(ILAY2-1)*(PRES(ILAY2-1)-PRES(ILAY2))*
     &              (AOZONE(ILAY2)/ROZONE(ILAY2))/RTEMP(ILAY2)
               IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1091
 1092       CONTINUE
 1091    CONTINUE
         DQDT_I(9)=1
C............................................................................
C     Q-PROFILE 11
C............................................................................
      ELSEIF (QPID .EQ. 11) THEN
         RDMMY1=0.0
         DO 1111 ILAY1=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
            RDMMY1=RDMMY1+RDMMY3*ROZONE(ILAY1)*RTEMP(ILAY1)
            DO 1112 ILAY2=ILAY1,(NLAYER-1)
               RDMMY2=PRES(ILAY2)*(PRES(ILAY2)-PRES(ILAY2+1))
     &                *AOZONE(ILAY2)
               DQDT(11,ILAY1,ILAY2)=RDMMY2/RDMMY1
               IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 1111
 1112       CONTINUE
 1111    CONTINUE
         DQDT(11,NLAYER,NLAYER)=AOZONE(NLAYER)/
     &                          (ROZONE(NLAYER)*RTEMP(NLAYER))
         DQDT_I(11)=1
C............................................................................
C     Q-PROFILES 4,5,7,8,10
C............................................................................
      ELSEIF ((QPID .EQ. 4) .OR. (QPID .EQ. 5) .OR. (QPID .EQ. 7)) THEN
         DQDT_I(QPID)=1
      ELSEIF ((QPID .EQ. 8) .OR. (QPID .EQ. 10)) THEN
         DQDT_I(QPID)=1
C............................................................................
      ELSE
         WRITE (*,*) 'FATAL ERROR IN I_DQDT: NO SUCH Q-PROFILE'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF

      RETURN
      END
C     END OF I_DQDT

C ====================================================================

      SUBROUTINE I_DQDW(QPID,DQDW_I,DQDW,
     &                  ATEMP,AFIXED,AWATER,AOZONE,
     &                  RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     initialize dQ/dW-profile number QPID
C
C     INPUT:
C       QPID            ID of Q-profile which shall be initialized
C       ATEMP,AFIXED,AWATER,AOZONE       atmosphere profile          
C       RTEMP,RFIXED,RWATER,ROZONE,PRES  reference atmosphere profile
C       VERBOS          screen output on/off (1=on,0=off)
C
C     INPUT AND OUTPUT:
C       DQDW_I          dQ/dW-profiles initialization flags
C       DQDW            dQ/dW-profiles
C
C     CALLING ROUNTINE
C       C_DKDW
C
C     HISTORY:
C       written 12/5/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE CALLING VARIABLES
      INTEGER     QPID,DQDW_I(MAX_QP),VERBOS
      REAL*8      DQDW(MAX_QP,NLAYER,NLAYER), ! DQDW(?,i,j)=dQ(i)/dW(j)
     &            ATEMP(NLAYER),AFIXED(NLAYER),AWATER(NLAYER),
     &            AOZONE(NLAYER),RTEMP(NLAYER),RFIXED(NLAYER),
     &            RWATER(NLAYER),ROZONE(NLAYER),PRES(NLAYER)
C     LOCAL VARIABLES
      INTEGER     ILAY1,ILAY2
      REAL*8      RDMMY1,RDMMY2,RDMMY3

C     check if initialization flag is already set
      IF (DQDW_I(QPID) .EQ. 1) THEN
         RETURN
      ELSE
         IF (VERBOS .EQ. 1) THEN
            WRITE (*,20) QPID
 20         FORMAT(' initializing dQ/dW-profile ',I2)
         ENDIF
      ENDIF
      
C     calculate Q-profiles:
C     NOTE: ALL DQDW ARE ALREADY INITIALIZES AS ZEROS
C............................................................................
C     Q-PROFILES 1 - 3
C............................................................................
      IF ((QPID .GE. 1) .AND. (QPID .LE. 3)) THEN
        DQDW_I(QPID)=1
C............................................................................
C     Q-PROFILE 4
C............................................................................
      ELSEIF (QPID .EQ. 4) THEN
         DO 1040 ILAY1=1,NLAYER
            DQDW(4,ILAY1,ILAY1)=1/RWATER(ILAY1)
 1040    CONTINUE
         DQDW_I(4)=1
C............................................................................
C     Q-PROFILE 5
C............................................................................
      ELSEIF (QPID .EQ. 5) THEN
         RDMMY1=0.0
         DO 1051 ILAY1=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
            RDMMY1=RDMMY1+RDMMY3*RWATER(ILAY1)
            DO 1052 ILAY2=ILAY1,(NLAYER-1)
               RDMMY2=PRES(ILAY2)*(PRES(ILAY2)-PRES(ILAY2+1))
               DQDW(5,ILAY1,ILAY2)=RDMMY2/RDMMY1
               IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 1051
 1052       CONTINUE
 1051    CONTINUE
         DQDW(5,NLAYER,NLAYER)=1/RWATER(NLAYER)
         DQDW_I(5)=1
C............................................................................
C     Q-PROFILE 6
C............................................................................
      ELSEIF (QPID .EQ. 6) THEN
         RDMMY1=0.0
         DO 1061 ILAY1=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
            RDMMY1=RDMMY1+RDMMY3*RWATER(ILAY1)*RTEMP(ILAY1)
            DO 1062 ILAY2=ILAY1,(NLAYER-1)
               RDMMY2=PRES(ILAY2)*(PRES(ILAY2)-PRES(ILAY2+1))
     &                *ATEMP(ILAY2)
               DQDW(6,ILAY1,ILAY2)=RDMMY2/RDMMY1
               IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 1061
 1062       CONTINUE
 1061    CONTINUE
         DQDW(6,NLAYER,NLAYER)=ATEMP(NLAYER)
     &                         /(RWATER(NLAYER)*RTEMP(NLAYER))
         DQDW_I(6)=1
C............................................................................
C     Q-PROFILES 7 - 11
C............................................................................
         ELSEIF ((QPID .GT. 6) .AND. (QPID .LE. 11)) THEN
            DQDW_I(QPID)=1
C............................................................................
      ELSE
         WRITE (*,*) 'FATAL ERROR IN I_DQDW: NO SUCH Q-PROFILE'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF

      RETURN
      END
C     END OF I_DQDW


C ====================================================================

      SUBROUTINE I_DQDO(QPID,DQDO_I,DQDO,
     &                  ATEMP,AFIXED,AWATER,AOZONE,
     &                  RTEMP,RFIXED,RWATER,ROZONE,PRES,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     initialize dQ/dO-profile number QPID
C
C     INPUT:
C       QPID            ID of Q-profile which shall be initialized
C       ATEMP,AFIXED,AWATER,AOZONE       atmosphere profile          
C       RTEMP,RFIXED,RWATER,ROZONE,PRES  reference atmosphere profile
C       VERBOS          screen output on/off (1=on,0=off)
C
C     INPUT AND OUTPUT:
C       DQDO_I          dQ/dO-profiles initialization flags
C       DQDO            dQ/dO-profiles
C
C     CALLING ROUNTINE
C       C_DKDO
C
C     HISTORY:
C       written 12/5/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE CALLING VARIABLES
      INTEGER     QPID,DQDO_I(MAX_QP),VERBOS
      REAL*8      DQDO(MAX_QP,NLAYER,NLAYER), ! DQDO(?,i,j)=dQ(i)/dO(j)
     &            ATEMP(NLAYER),AFIXED(NLAYER),AWATER(NLAYER),
     &            AOZONE(NLAYER),RTEMP(NLAYER),RFIXED(NLAYER),
     &            RWATER(NLAYER),ROZONE(NLAYER),PRES(NLAYER)
C     LOCAL VARIABLES
      INTEGER     ILAY1,ILAY2
      REAL*8      RDMMY1,RDMMY2,RDMMY3

C     check if initialization flag is already set
      IF (DQDO_I(QPID) .EQ. 1) THEN
         RETURN
      ELSE
         IF (VERBOS .EQ. 1) THEN
            WRITE (*,20) QPID
 20         FORMAT(' initializing dQ/dO-profile ',I2)
         ENDIF
      ENDIF
      

C     calculate Q-profiles:
C     NOTE: ALL DQDO ARE ALREADY INITIALIZES AS ZEROS
C............................................................................
C     Q-PROFILES 1 - 6
C............................................................................
      IF ((QPID .GE. 1) .AND. (QPID .LE. 6)) THEN
         DQDO_I(QPID)=1
C............................................................................
C     Q-PROFILE 7
C............................................................................
      ELSEIF (QPID .EQ. 7) THEN
         DO 1070 ILAY1=1,NLAYER
            DQDO(7,ILAY1,ILAY1)=1/ROZONE(ILAY1)
 1070    CONTINUE
         DQDO_I(7)=1
C............................................................................
C     Q-PROFILE 8
C............................................................................
      ELSEIF (QPID .EQ. 8) THEN
         DO 1081 ILAY1=1,(NLAYER-1)
            DO 1082 ILAY2=(ILAY1+1),NLAYER
               DQDO(8,ILAY1,ILAY2)=
     &              PRES(ILAY2-1)*(PRES(ILAY2-1)-PRES(ILAY2))/
     &              ROZONE(ILAY2)
               IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 1081
 1082       CONTINUE
 1081    CONTINUE
         DQDO_I(8)=1
C............................................................................
C     Q-PROFILE 9
C............................................................................
      ELSEIF (QPID .EQ. 9) THEN
         DO 1091 ILAY1=1,(NLAYER-1)
            DO 1092 ILAY2=(ILAY1+1),NLAYER
               DQDO(9,ILAY1,ILAY2)=
     &              PRES(ILAY2-1)*(PRES(ILAY2-1)-PRES(ILAY2))*
     &              ATEMP(ILAY2)/(RTEMP(ILAY2)*ROZONE(ILAY2))
               IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 1091
 1092       CONTINUE
 1091    CONTINUE
         DQDO_I(9)=1
C............................................................................
C     Q-PROFILE 10
C............................................................................
      ELSEIF (QPID .EQ. 10) THEN
         RDMMY1=0.0
         DO 1101 ILAY1=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
            RDMMY1=RDMMY1+RDMMY3*ROZONE(ILAY1)
            DO 1102 ILAY2=ILAY1,(NLAYER-1)
               RDMMY2=PRES(ILAY2)*(PRES(ILAY2)-PRES(ILAY2+1))
               DQDO(10,ILAY1,ILAY2)=RDMMY2/RDMMY1
               IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 1101
 1102       CONTINUE
 1101    CONTINUE
         DQDO(10,NLAYER,NLAYER)=1/ROZONE(NLAYER)
         DQDO_I(10)=1
C............................................................................
C     Q-PROFILE 11
C............................................................................
      ELSEIF (QPID .EQ. 11) THEN
         RDMMY1=0.0
         DO 1111 ILAY1=(NLAYER-1),1,-1
            RDMMY3=PRES(ILAY1)*(PRES(ILAY1)-PRES(ILAY1+1))
            RDMMY1=RDMMY1+RDMMY3*ROZONE(ILAY1)*RTEMP(ILAY1)
            DO 1112 ILAY2=ILAY1,(NLAYER-1)
               RDMMY2=PRES(ILAY2)*(PRES(ILAY2)-PRES(ILAY2+1))*
     &                ATEMP(ILAY2)
               DQDO(11,ILAY1,ILAY2)=RDMMY2/RDMMY1
               IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 1111
 1112       CONTINUE
 1111    CONTINUE
         DQDO(11,NLAYER,NLAYER)=ATEMP(NLAYER)/
     &                          (RTEMP(NLAYER)*ROZONE(NLAYER))
         DQDO_I(11)=1
C............................................................................
      ELSE
         WRITE (*,*) 'FATAL ERROR IN I_DQDO: NO SUCH Q-PROFILE'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF

      RETURN
      END
C     END OF I_DQDO

C ====================================================================

      SUBROUTINE C_DLST(CLCJAC,TOTIRF,ID_IRF,TAU,DTAUT,DTAUW,DTAUO)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       convert derivatives of layer transmittances to derivatives
C       of layer-to-space derivatives
C
C     INPUT:
C       TOTIRF      number of channels
C       ID_IRF      channel IDs
C       TAU         transmittances
C
C     INPUT AND OUTPUT:
C       DTAUT       d(tau)/dT (layer for INPUT, layer-to-space for OUTPUT)
C       DTAUW       d(tau)/dW (layer for INPUT, layer-to-space for OUTPUT)
C       DTAUO       d(tau)/dO (layer for INPUT, layer-to-space for OUTPUT)
C
C     CALLING ROUNTINE
C       KERNEL
C
C     HISTORY:
C       written 12/11/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE CALLING VARIABLES
      INTEGER        CLCJAC(3),          ! CALC JACOBIANS (T,W,O)
     &               TOTIRF,             ! TOTAL NUMBER OF IRFs
     &               ID_IRF(MAXIRF)      ! IDs OF IRFs 
      REAL*8         TAU(MAXIRF,NLAYER), ! LAYER TRANSMITTANCES
     &               DTAUT(MAXIRF,NLAYER,NLAYER),  ! dtau/dT, thermal upw.
     &               DTAUW(MAXIRF,NLAYER,NLAYER),  ! dtau/dW, thermal upw.
     &               DTAUO(MAXIRF,NLAYER,NLAYER)   ! dtau/dO, thermal upw.
C     LOCAL VARIABLES
      INTEGER        IIRF,IIRF_A,ILAY1,ILAY2
      REAL*8         TAULS(MAXIRF,NLAYER) ! layer-to-space transmittances 

C     BEGIN
      ! calculate the layer-to-space transmittances
      ! Note: the layer-to-space transmittance always refers to the
      !       transmittance above the bottom of the respective layer
      DO 100 IIRF=1,TOTIRF
         IIRF_A=ID_IRF(IIRF)
         ! for the top layer the layer-to-space transmittance is
         ! equal to the layer transmittance
         TAULS(IIRF_A,NLAYER)=TAU(IIRF_A,NLAYER)
         DO 101 ILAY1=(NLAYER-1),1,-1
            TAULS(IIRF_A,ILAY1)=TAULS(IIRF_A,ILAY1+1)*TAU(IIRF_A,ILAY1)
 101     CONTINUE
 100  CONTINUE

      DO 1000 IIRF=1,TOTIRF
         IIRF_A=ID_IRF(IIRF)
         ! all DTAUX(IIRF_A,NLAYER,ILAY2) remain unchanged
         IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
            ! TEMPERATURE
            DO 2011 ILAY1=(NLAYER-1),1,-1
               DO 2012 ILAY2=ILAY1,NLAYER
                  DTAUT(IIRF_A,ILAY1,ILAY2)=
     &                 DTAUT(IIRF_A,ILAY1,ILAY2)*TAULS(IIRF_A,ILAY1+1)+
     &                 DTAUT(IIRF_A,ILAY1+1,ILAY2)*TAU(IIRF_A,ILAY1)
                  IF (ILAY2 .EQ. (ILAY1+MILAYT)) GOTO 2011
 2012          CONTINUE
 2011       CONTINUE
         ENDIF
         IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
            ! WATER
            DO 2021 ILAY1=(NLAYER-1),1,-1
               DO 2022 ILAY2=ILAY1,NLAYER
                  DTAUW(IIRF_A,ILAY1,ILAY2)=
     &                 DTAUW(IIRF_A,ILAY1,ILAY2)*TAULS(IIRF_A,ILAY1+1)+
     &                 DTAUW(IIRF_A,ILAY1+1,ILAY2)*TAU(IIRF_A,ILAY1)
                  IF (ILAY2 .EQ. (ILAY1+MILAYW)) GOTO 2021
 2022          CONTINUE
 2021       CONTINUE
         ENDIF
         IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
            ! OZONE
            DO 2031 ILAY1=(NLAYER-1),1,-1
               DO 2032 ILAY2=ILAY1,NLAYER
                  DTAUO(IIRF_A,ILAY1,ILAY2)=
     &                 DTAUO(IIRF_A,ILAY1,ILAY2)*TAULS(IIRF_A,ILAY1+1)+
     &                 DTAUO(IIRF_A,ILAY1+1,ILAY2)*TAU(IIRF_A,ILAY1)
                  IF (ILAY2 .EQ. (ILAY1+MILAYO)) GOTO 2031
 2032          CONTINUE
 2031       CONTINUE
         ENDIF
 1000 CONTINUE
      RETURN
      END
C     END OF SUBROUTINE C_DLST

C ====================================================================

      SUBROUTINE ADDSEC(SEC1,SEC2,SECSUM)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       add two secant angles (convert secant to angles, add angles,
C       convert sum angle to secant, take absolute value)
C
C     INPUT:
C       SEC1        first secant     
C       SEC2        secont secant
C
C     OUTPUT:
C       SECSUM      abs(secant of sum of angles of 1st and 2nd secant)
C
C     HISTORY:
C       written 05/15/98 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE

C     SUBROUTINE CALLING VARIABLES
      REAL*8         SEC1,SEC2,SECSUM

C     BEGIN
      SECSUM = ABS(1.0/DCOS(DACOS(1.0/SEC1)+DACOS(1.0/SEC2)))

      RETURN
      END
C     END OF SUBROUTINE ADDSEC

C ====================================================================

