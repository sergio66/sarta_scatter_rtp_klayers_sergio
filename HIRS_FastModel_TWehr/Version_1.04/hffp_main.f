C FILE NAME: hffp_main.f
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
C HOW TO CALL THE PROGRAM
C
C     command line:
C     echo "controlfile.ctl" | hirs_main.x
C     or use the shell script hirs_main.sh
C
C HISTORY:
C     written by Tobias Wehr, 1997
C-----------------------------------------------------------------------

      PROGRAM HFFP
      IMPLICIT NONE

C     INCLUDE FILES
      include "hffp_glob_dec.f"
      include "hffp_init_dec.f"

C     LOCAL PARAMETERS
      INTEGER IJOB,           ! job (atmosphere) number counter
     &        CLCJAC(3),      ! CALCULATE JACOBIANS FOR T,W,O: 0=NO,1=YES(ANA-
                              ! LYTICAL JACOB.),2=YES(FINITE DIFFERENCES)
     &        JACTYP,         ! TYPE OF JACOBIANS: 1=DERIVATIVES OF RADIANCES
                              !                    2=DERIVATIVES OF T_BRIGHT
     &        CLCTBR,         ! CALCULATE T_BRIGHT (1=YES, 0=NO)
     &        OT_RTB          ! OFFSET TUNING OF RADIANCES (1) OR T_BRIGHT (2)
      REAL*8  CW_M(MAXIRF),   ! CHANNEL CENTER WAVE NUMBERS as used 
                              !   for calculation of TBMEAN
     &        RADTOT(MAXIRF), ! RADIANCES TOTAL (UPWELLING+REFLECED)
     &        TBMEAN(MAXIRF), ! BRIGTH. TEMPERATURE (calculated for 
                              !   respective mean center wavenumber)
     &        JACOBT(MAXIRF,NLAYER),   ! ANALYTIC JACOBIANS dY/dT (Y=RAD
     &        JACOBW(MAXIRF,NLAYER),   ! ANALYTIC JACOBIANS dY/dW    or
     &        JACOBO(MAXIRF,NLAYER)    ! ANALYTIC JACOBIANS dY/dO    TB)

      REAL*8  NEWTBR(MAXIRF), ! brightness temp. with disturbed atmosphere
                              ! for finite differences T_bright-Jacobians
     &        NEWRAD(MAXIRF), ! radiances with dist. atm. (for finite diff.
                              ! radiance-Jacobians
     &        FDJ_T(MAXIRF,NLAYER), ! f.d. Jacobians for temperature
     &        FDJ_W(MAXIRF,NLAYER), ! f.d. Jacobians for water
     &        FDJ_O(MAXIRF,NLAYER), ! f.d. Jacobians for ozone
     &        DISTPR(NLAYER), ! disturbed profile (either T,W or O)
     &        DELTAT,         ! temperature disturbance for f.d. Jac.
     &        DELTAW,         ! water profile disturb.  for f.d. Jac.
     &        DELTAO,         ! ozone profile disturb.  for f.d. Jac.
     &        DELTAX          ! dummy for abs. W or O disturbance for f.d.Jac.

      INTEGER VERBOS,          ! verbose parameter
     &        T_ADJ,           ! 1: the transmittance tuning parameters 
                               !    can be changed by SUBROUTINE C_K
                               !    (called by KERNEL), if found to be 
                               !    unphysical
                               ! 0: the transmittance tuning parameters 
                               !    will not be changed by SUBROUTINE C_K, 
                               !    (called by KERNEL), even if unphysical
     &        NOTECT           ! note if T_TAUX have changed in KERNEL
                               ! 1=yes, 0=no
      INTEGER ILAY1,ILAY2,IIRF,IIRF_A,IBENCH, ! misc. counter variables
     &        IDUMMY
      CHARACTER*80 DATMFN(MAXJOB)  ! dummy atmosphere filenames for f.d. Jac.
     &
      REAL*8  JTDMMY,JWDMMY,JODMMY ! dummies
      INTEGER DCLCJC(3)            ! dummy
      DATA DCLCJC/3*0/             ! DO NOT CHANGE THIS DATA STATEMENT

      DATA DELTAT/1.0/        ! set temp. disturbance to 1 Kelvin
      DATA DELTAW/1.0/        ! set water disturbance to 1 percent
      DATA DELTAO/1.0/        ! set ozone disturbance to 1 percent

C     ==== user can change this! ==========================
      DATA VERBOS/1/           ! SCREEN OUTPUT: 0=OFF, 1=ON
C     =====================================================

C     ====================================================================
C     BEGIN
C     ====================================================================

C     GREETING
      !IF (VERBOS .EQ. 1) THEN
      WRITE (*,*) '==================================='
      WRITE (*,*) 'HIRS FAST FORWARD PROGRAM (HFFP)'
      WRITE (*,*) HFFPVS
      WRITE (*,*) '==================================='
      !ENDIF

C     GET CONTROL FILE NAME
      WRITE (*,*) 'name of control file (including path): '
      READ (*,12) CTLDMY
 12   FORMAT(A80)
      CTLFIL=CTLDMY

C     GET INITIAL DATA
      CALL HFFP_INIT(CTLFIL,ATMFLN,SATNUM,NOJOBS,TOTIRF,ID_IRF,
     &               PREDID,NCOEFF,COEFF,T_TAUF,T_TAUW,T_TAUO,T_OFFS,
     &               T_SEC,
     &               SURF_T,SURF_P,SECANT,SUNSEC,SUREMI,SUNSOL,
     &               SOLRAD,CLCJAC,WRTOUT,VERBOS,JACTYP,CLCTBR,OT_RTB,
     &               T_ADJ)

      DO 2000 IBENCH=1,1 ! USE THIS LOOP FOR BENCHMARKING
         IF (IBENCH .GT. 1) WRITE (*,*) 'BENCHMARK LOOP ',IBENCH
         DO 1000 IJOB=1,NOJOBS
            IF (VERBOS .EQ. 1) THEN
               IF (NOJOBS .GT. 1) THEN
                  WRITE (*,*) 'CALCULATE ATMOSPHERE NO.',IJOB
               ENDIF
            ENDIF
            CALL GETATM(ATMFLN(IJOB),ATEMP,AFIXED,AWATER,AOZONE,VERBOS)
C           ======================
C           DO FORWARD CALCULATION
C           ======================
            CALL KERNEL(IJOB,SATNUM,TOTIRF,ID_IRF,
     &           PREDID,NCOEFF,COEFF,T_TAUF,T_TAUW,T_TAUO,T_OFFS,T_SEC,
     &           SURF_T(IJOB),SURF_P(IJOB),SECANT(IJOB),
     &           SUNSEC(IJOB),SUREMI(IJOB),SUNSOL(IJOB),
     &           SOLRAD,ATEMP,AFIXED,AWATER,AOZONE,CLCJAC,CW_M,RADTOT,
     &           TBMEAN,JACOBT,JACOBW,JACOBO,VERBOS,NOTECT,JACTYP,
     &           CLCTBR,OT_RTB,T_ADJ)
C           WRITE RADIATION AND BRIGHTNESS TEMP. RESULTS TO OUTPUT
            IF (WRTOUT .EQ. 1) THEN
               CALL ROUT(CTLFIL,TOTIRF,ID_IRF,RADTOT,TBMEAN,CW_M,IJOB,
     &                   VERBOS,CLCTBR)
            ENDIF
C           ============================
C           FINITE DIFFERENCES JACOBIANS
C           ============================
            ! temperature Jacobians
            IF ((CLCJAC(1) .EQ. 2) .OR. (CLCJAC(1) .EQ. 3)) THEN
               IF (VERBOS .EQ. 1) THEN
                  WRITE (*,*) 'CALC FIN.DIF. JACOBIANS dY/dT'
               ENDIF
               DO 1011 ILAY1=1,NLAYER
                  DO 1012 ILAY2=1,NLAYER
                     DISTPR(ILAY2)=ATEMP(ILAY2)
 1012             CONTINUE
                  DISTPR(ILAY1)=DISTPR(ILAY1)+DELTAT
                  CALL KERNEL(IJOB,SATNUM,TOTIRF,ID_IRF,
     &                 PREDID,NCOEFF,COEFF,T_TAUF,T_TAUW,T_TAUO,T_OFFS,
     &                 T_SEC,
     &                 SURF_T(IJOB),SURF_P(IJOB),SECANT(IJOB),
     &                 SUNSEC(IJOB),SUREMI(IJOB),SUNSOL(IJOB),
     &                 SOLRAD,DISTPR,AFIXED,AWATER,AOZONE,DCLCJC,
     &                 CW_M,NEWRAD,NEWTBR,JTDMMY,JWDMMY,JODMMY,0,
     &                 IDUMMY,JACTYP,CLCTBR,OT_RTB,T_ADJ)
                  DO 1013 IIRF=1,TOTIRF
                     IIRF_A=ID_IRF(IIRF)
                     IF (JACTYP .EQ. 1) THEN
                        FDJ_T(IIRF_A,ILAY1)=
     &                       (NEWRAD(IIRF_A)-RADTOT(IIRF_A))/DELTAT
                     ELSE
                        FDJ_T(IIRF_A,ILAY1)=
     &                       (NEWTBR(IIRF_A)-TBMEAN(IIRF_A))/DELTAT
                     ENDIF
 1013             CONTINUE
 1011          CONTINUE
            ENDIF
            ! water Jacobians
            IF ((CLCJAC(2) .EQ. 2) .OR. (CLCJAC(2) .EQ. 3)) THEN
               IF (VERBOS .EQ. 1) THEN
                  WRITE (*,*) 'CALC FIN.DIF. JACOBIANS dY/dW'
               ENDIF
               DO 1021 ILAY1=1,NLAYER
                  DO 1022 ILAY2=1,NLAYER
                     DISTPR(ILAY2)=AWATER(ILAY2)
 1022             CONTINUE
                  DELTAX=DISTPR(ILAY1)*DELTAW/100.0
                  DISTPR(ILAY1)=DISTPR(ILAY1)+DELTAX
                  CALL KERNEL(IJOB,SATNUM,TOTIRF,ID_IRF,
     &                 PREDID,NCOEFF,COEFF,T_TAUF,T_TAUW,T_TAUO,T_OFFS,
     &                 T_SEC,
     &                 SURF_T(IJOB),SURF_P(IJOB),SECANT(IJOB),
     &                 SUNSEC(IJOB),SUREMI(IJOB),SUNSOL(IJOB),
     &                 SOLRAD,ATEMP,AFIXED,DISTPR,AOZONE,DCLCJC,
     &                 CW_M,NEWRAD,NEWTBR,JTDMMY,JWDMMY,JODMMY,0,
     &                 IDUMMY,JACTYP,CLCTBR,OT_RTB,T_ADJ)
                  DO 1023 IIRF=1,TOTIRF
                     IIRF_A=ID_IRF(IIRF)
                     IF (JACTYP .EQ. 1) THEN
                        FDJ_W(IIRF_A,ILAY1)=
     &                       (NEWRAD(IIRF_A)-RADTOT(IIRF_A))/DELTAX
                     ELSE
                        FDJ_W(IIRF_A,ILAY1)=
     &                       (NEWTBR(IIRF_A)-TBMEAN(IIRF_A))/DELTAX
                     ENDIF
 1023             CONTINUE
 1021          CONTINUE
            ENDIF
            ! ozone Jacobians
            IF ((CLCJAC(3) .EQ. 2) .OR. (CLCJAC(3) .EQ. 3)) THEN
               IF (VERBOS .EQ. 1) THEN
                  WRITE (*,*) 'CALC FIN.DIF. JACOBIANS dY/dO'
               ENDIF
               DO 1031 ILAY1=1,NLAYER
                  DO 1032 ILAY2=1,NLAYER
                     DISTPR(ILAY2)=AOZONE(ILAY2)
 1032             CONTINUE
                  DELTAX=DISTPR(ILAY1)*DELTAO/100.0
                  DISTPR(ILAY1)=DISTPR(ILAY1)+DELTAX
                  CALL KERNEL(IJOB,SATNUM,TOTIRF,ID_IRF,
     &                 PREDID,NCOEFF,COEFF,T_TAUF,T_TAUW,T_TAUO,T_OFFS,
     &                 T_SEC,
     &                 SURF_T(IJOB),SURF_P(IJOB),SECANT(IJOB),
     &                 SUNSEC(IJOB),SUREMI(IJOB),SUNSOL(IJOB),
     &                 SOLRAD,ATEMP,AFIXED,AWATER,DISTPR,DCLCJC,
     &                 CW_M,NEWRAD,NEWTBR,JTDMMY,JWDMMY,JODMMY,0,
     &                 IDUMMY,JACTYP,CLCTBR,OT_RTB,T_ADJ)
                  DO 1033 IIRF=1,TOTIRF
                     IIRF_A=ID_IRF(IIRF)
                     IF (JACTYP .EQ. 1) THEN
                        FDJ_O(IIRF_A,ILAY1)=
     &                       (NEWRAD(IIRF_A)-RADTOT(IIRF_A))/DELTAX
                     ELSE
                        FDJ_O(IIRF_A,ILAY1)=
     &                       (NEWTBR(IIRF_A)-TBMEAN(IIRF_A))/DELTAX
                     ENDIF
 1033             CONTINUE
 1031          CONTINUE
            ENDIF
C           === end of finite differences Jacobians ===
C           write Jacobians to files
            IF (WRTOUT .EQ. 1)  THEN
               IF ((CLCJAC(1) .NE. 0) .OR. (CLCJAC(2) .NE. 0) .OR.
     &              (CLCJAC(3) .NE. 0)) THEN
                  CALL WR_JAC(CTLFIL,TOTIRF,ID_IRF,CLCJAC,IJOB,
     &                 JACOBT,JACOBW,JACOBO,FDJ_T,FDJ_W,FDJ_O,VERBOS)
               ENDIF
            ENDIF

 1000    CONTINUE
 2000 CONTINUE 

      IF ((VERBOS .EQ. 1) .AND. (NOTECT .EQ. 1)) THEN
         WRITE (*,*) 'WARNING: TUNING PARAMETERS HAVE BEEN CHANGED.'
         WRITE (*,*) '         USED SET OF TUNING PARAMETERS:'
         WRITE (*,*) '         CHANNEL   FIXED       WATER       OZONE'
         DO 3001 IIRF=1,TOTIRF
            IIRF_A=ID_IRF(IIRF)
            WRITE (*,3002) IIRF_A,
     &                     T_TAUF(IIRF_A),T_TAUW(IIRF_A),T_TAUO(IIRF_A)
 3001    CONTINUE
 3002    FORMAT('         ',I5,'  ',F12.6,F12.6,F12.6)
      ENDIF

      !IF (VERBOS .EQ. 1) THEN
      WRITE (*,*) 'PROGRAM SUCCESSFULLY COMPLETED'
      !ENDIF
      END
C     END OF MAIN PROGRAM


C ====================================================================

      SUBROUTINE HFFP_INIT(CTLFIL,ATMFLN,SATNUM,NOJOBS,TOTIRF,ID_IRF,
     &               PREDID,NCOEFF,COEFF,T_TAUF,T_TAUW,T_TAUO,T_OFFS,
     &               T_SEC,
     &               SURF_T,SURF_P,SECANT,SUNSEC,SUREMI,SUNSOL,
     &               SOLRAD,CLCJAC,WRTOUT,VERBOS,JACTYP,CLCTBR,OT_RTB,
     &               T_ADJ)

C --------------------------------------------------------------------
C
C     PURPOSE:
C       - Initialize all parameters. 
C       - This subroutine is only used by the hffp_main program, 
C         not by the NASA-DAO system.
C
C     INPUT
C       CTLFIL   name of control file (see http_init_dec.f)
C
C     OUTPUT 
C       everything else, as described in hffp_init_dec.f
C 
C     NOTE: 
C       ATEMP,AFIXED,AWATER,AOZONE can also be defined (read from
C       file) within the main subroutine (KERNEL)
C
C     HISTORY:
C       written 11/19/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      include "hffp_glob_dec.f"
C     CALLING PARAMETERS
      include "hffp_init_dec.f"
      INTEGER        VERBOS
      INTEGER        CLCJAC(3)     ! CALCULATE JACOBIANS FOR T,W,O
                     ! 0=NO,1=YES(ANALYTICAL J.),2=YES(FINITE DIFFERENCES)
      INTEGER        JACTYP,CLCTBR,OT_RTB,T_ADJ
      
C     BEGIN
      IF (VERBOS .EQ. 1) THEN
         WRITE (*,*) 'INIT ', HFFPVS
      ENDIF
      IF ((VERBOS .EQ. 1) .OR. (MATLAB .EQ. 1)) THEN
         WRITE(*,13) CTLFIL
 13      FORMAT(' control file: ',A80)
      ENDIF
C     READ CONTROL FILE
      CALL RDCTL(CTLFIL,SATNUM,NOJOBS,TOTIRF,ID_IRF,SURF_T,SURF_P,
     &           T_TAUF,T_TAUW,T_TAUO,T_OFFS,T_SEC,ATMFLN,SECANT,SUNSEC,
     &           SUREMI,SUNSOL,CLCJAC,WRTOUT,JACTYP,CLCTBR,OT_RTB,
     &           T_ADJ)
C     READ PREDICTOR IDs AND COEFFICIENTS
C     read predictor IDs from file
      CALL RDPRID(SATNUM,PREDID,VERBOS)
C     read coefficients from file
      CALL RDCOEF(SATNUM,NCOEFF,COEFF,VERBOS)
C     read in solar radiances
      CALL RDSOLR(SATNUM,SOLRAD,VERBOS)
C     =================================================================

      RETURN
      END
C     END OF SUBROUTINE HFFP_INIT

C ====================================================================

      SUBROUTINE RDCTL(CTLFIL,SATNUM,NOJOBS,TOTIRF,ID_IRF,
     &                 SURF_T,SURF_P,
     &                 T_TAUF,T_TAUW,
     &                 T_TAUO,T_OFFS,T_SEC,
     &                 ATMFLN,SECANT,SUNSEC,SUREMI,SUNSOL,
     &                 CLCJAC,WRTOUT,JACTYP,CLCTBR,OT_RTB,T_ADJ)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     read control file, all returned values are read from control file
C
C     INPUT:
C       CTLFIL           name of control file
C
C     OUTPUT:
C       all other parameter, explanation see main program
C
C     HISTORY:
C       written 8/26/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE ARGUMENT PARAMETER
      CHARACTER*80   CTLFIL,ATMFLN(MAXJOB)
      INTEGER        SATNUM,NOJOBS,TOTIRF,ID_IRF(MAXIRF),
     &               CLCJAC(3),WRTOUT,JACTYP,CLCTBR,OT_RTB,T_ADJ
      REAL*8         SURF_T(MAXJOB),SURF_P(MAXJOB),
     &               T_TAUF(MAXIRF),T_TAUW(MAXIRF),T_TAUO(MAXIRF),
     &               T_OFFS(MAXIRF),T_SEC(MAXIRF)
      REAL*8         SECANT(MAXJOB),SUNSEC(MAXJOB),SUREMI(MAXJOB),
     &               SUNSOL(MAXJOB)
C     LOCAL PARAMETER
      CHARACTER*40   CTLLNE,REQLNE,REQLNA
      INTEGER        FP,IERR,IIRF,IATM
      REAL*8         RANGLE,PI
      PARAMETER(PI=3.141592654)

 1000 FORMAT(A40)
 1010 FORMAT('ERROR: keyword not found: ', A40, /,
     &     '       I read:            ', A40 )
      FP=10
      OPEN (UNIT=FP,FILE=CTLFIL,FORM='FORMATTED',
     &     STATUS='OLD',IOSTAT=IERR)
C     READ SATELLITE NUMBER
      REQLNE='SATELLITE NUMBER'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         READ (FP,*) SATNUM
      ENDIF
C     READ TOTAL NUMBER OF IRF
      REQLNE='TOTAL NUMBER OF CHANNELS'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         READ (FP,*) TOTIRF
      ENDIF
C     READ ID NUMBERS OF IRF
      REQLNE='ID NUMBERS OF CHANNELS'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1020 IIRF=1,TOTIRF
            READ (FP,*) ID_IRF(IIRF)
 1020    CONTINUE
      ENDIF
C     READ TUNING PARAMETERS FOR TAU(FIXED), TAU(H2O), TAU(O3)
      REQLNE='TAU(FIXED) TUNING PARAMETERS'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1030 IIRF=1,TOTIRF
            READ (FP,*) T_TAUF(ID_IRF(IIRF))
 1030    CONTINUE
      ENDIF
      REQLNE='TAU(H2O) TUNING PARAMETERS'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1040 IIRF=1,TOTIRF
            READ (FP,*) T_TAUW(ID_IRF(IIRF))
 1040    CONTINUE
      ENDIF
      REQLNE='TAU(O3) TUNING PARAMETERS'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1050 IIRF=1,TOTIRF
            READ (FP,*) T_TAUO(ID_IRF(IIRF))
 1050    CONTINUE
      ENDIF
C     READ OFFSET TUNING PARAMERTERS
      REQLNE='OFFSET TUNING PARAMETER'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1060 IIRF=1,TOTIRF
            READ (FP,*) T_OFFS(ID_IRF(IIRF))
 1060    CONTINUE
      ENDIF
C     READ CHANNEL ZENITH ANGLE OFFSET AND CONVERT TO SECANT OFFSET
      REQLNE='CHANNEL ZENITH ANGLE OFFSET'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1061 IIRF=1,TOTIRF
            READ (FP,*) RANGLE
            T_SEC(ID_IRF(IIRF))=1/DCOS(2*PI*RANGLE/360.0) ! convert to secant
            !T_SEC(ID_IRF(IIRF))=RANGLE ! wrong!!! convert to secant!!!
 1061    CONTINUE
      ENDIF
C     READ WHETHER OFFSET TUNES RADIANCES OR BRIGHTNESS TEMPERATURES
      REQLNE='OFFSET TUNES RADIANCES (1) OR T_BRIGHT (2)'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         READ (FP,*) OT_RTB
         IF ((OT_RTB .NE. 1) .AND. (OT_RTB .NE. 2)) THEN
            WRITE (*,*) 'ERROR: OFFSET TUNES TYPE MUST BE 1 OR 2'
            STOP
         ENDIF
      ENDIF
C     READ T_ADJ: TAU TUNING PARAMETERS SHOULD BE ADJUSTED OR NOT
      REQLNE='ADJUST TAU TUNING PARAMETERS (0=no, 1=yes)'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         READ (FP,*) T_ADJ
      ENDIF
C     READ IF BRIGHTNESS TEMPERATURES ARE TO BE CALCULATED
      REQLNE='CALCULATE T_BRIGHT (0=NO, 1=YES)'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         READ (FP,*) CLCTBR
         IF ((CLCTBR .NE. 0) .AND. (CLCTBR .NE. 1)) THEN
            WRITE (*,*) 'ERROR: CALC. T_BRIGHT PARAMETER MUST BE 0 OR 1'
            STOP
         ENDIF
      ENDIF
C     READ FLAG FOR CALCULATING JACOBIANS
      REQLNE='CALCULATE JACOBIANS FOR TWO (0=NO,1=AJ,2=FD,3=AJ&FD)'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         READ (FP,*) CLCJAC(1),CLCJAC(2),CLCJAC(3)
         IF ((CLCJAC(1) .GT. 3) .OR. (CLCJAC(2) .GT. 3) .OR.
     &        (CLCJAC(3) .GT. 3)) THEN
            WRITE (*,*) 'ERROR: ILLEGAL PARAMETER IN CONTROL FILE'
            WRITE (*,*) 'FOR JACOBIANS: ',CLCJAC(1),CLCJAC(2),CLCJAC(3)
            STOP
         ENDIF
      ENDIF
C     READ IF JACOBIANS BE DERIVATIVES OF RADIANCES OR BRIGHTN.TEMP.
      REQLNE='JACOBIANS BE DERIVATIVES OF RAD (1) OR T_BRIGHT (2)'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         READ (FP,*) JACTYP
         IF ((JACTYP .NE. 1) .AND. (JACTYP .NE. 2)) THEN
            WRITE (*,*) 'ERROR: JACOBIANS TYPE MUST BE 1 OR 2'
            STOP
         ENDIF
      ENDIF
C     READ IF OUTPUT SHALL BE WRITTEN TO FILE
      REQLNE='WRITE OUTPUT TO FILE (1=YES, 0=NO)'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         READ (FP,*) WRTOUT
      ENDIF
C     READ NUMBER OF ATMOSPHERES
      REQLNE='NUMBER OF ATMOSPHERES'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         READ (FP,*) NOJOBS
         IF (NOJOBS .GT. MAXJOB) THEN
            WRITE (*,*) 'PROGRAM ERROR: NUMBER OF ATMOSPHERES TOO LARGE'
            WRITE (*,*) 'recompile program with larger value of MAXJOB'
            WRITE (*,*) '   in file hffp_glob_dec.f'
            WRITE (*,*) 'current value for MAXJOB (max. number of'
            WRITE (*,*) '   atmospheres) is ',MAXJOB
            WRITE (*,*) 'PROGRAM ABORTED'
            STOP
         ENDIF
      ENDIF
C     READ FILE NAMES OF ATMOSPHERIC PROFILE
      REQLNE='ATMOSPHERE PROFILE FILE NAMES'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1070 IATM=1,NOJOBS
            READ (FP,*) ATMFLN(IATM)
 1070    CONTINUE
      ENDIF
C     READ SURFACE PRESSURES
      REQLNE='SURFACE PRESSURES (unit mbar)'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1080 IATM=1,NOJOBS
            READ (FP,*) SURF_P(IATM)
 1080    CONTINUE
      ENDIF
C     READ SURFACE TEMPERATURES
      REQLNE='SURFACE TEMPERATURES (unit K)'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1090 IATM=1,NOJOBS
            READ (FP,*) SURF_T(IATM)
 1090    CONTINUE
      ENDIF
C     READ OBSERVATION ZENITH ANGLE AND CONVERT TO SECANT
      REQLNE='OBSERVATION ZENITH ANGLES'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         REQLNA=REQLNE
         REQLNE='OBSERVATION ZENITH SECANTS'
         IF (CTLLNE .NE. REQLNE) THEN
            WRITE (*,1010) REQLNE, CTLLNE
            WRITE (*,*) 'possible alternative for required keyword'
            WRITE (*,*) 'would be ',REQLNA
            STOP
         ELSE
            DO 1100 IATM=1,NOJOBS
               READ (FP,*) SECANT(IATM)
 1100       CONTINUE
         ENDIF
      ELSE
         DO 1110 IATM=1,NOJOBS
            READ (FP,*) RANGLE
            SECANT(IATM) = 1 / DCOS(2*PI*RANGLE/360.0)
 1110    CONTINUE
      ENDIF
C     READ SUN ZENITH ANGLE
      REQLNE='SUN ZENITH ANGLES'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         REQLNA=REQLNE
         REQLNE='SUN ZENITH SECANTS'
         IF (CTLLNE .NE. REQLNE) THEN
            WRITE (*,1010) REQLNE, CTLLNE
            WRITE (*,*) 'possible alternative for required keyword'
            WRITE (*,*) 'would be ',REQLNA
            STOP
         ELSE
            DO 1120 IATM=1,NOJOBS
               READ(FP,*) SUNSEC(IATM)
 1120       CONTINUE
         ENDIF
      ELSE
         DO 1130 IATM=1,NOJOBS
            READ (FP,*) RANGLE
            IF (RANGLE .GE. 90.0) THEN 
                                ! angle is too large (>=90 deg)
               SUNSEC(IATM)=MAXSNS
            ELSE
               SUNSEC(IATM) = 1 / DCOS(2*PI*RANGLE/360.0)
               IF (SUNSEC(IATM) .GE. MAXSNS) THEN 
                                ! too large angle (< 90deg .AND. close to 90deg)
                  SUNSEC(IATM)=MAXSNS
               ENDIF
            ENDIF
 1130    CONTINUE
      ENDIF
C     READ SOLID ANGLE OF THE SUN
                                ! as an average value is can be used:  PI*2.15976E-5
      REQLNE='SUN SOLID ANGLE'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1140 IATM=1,NOJOBS
            READ (FP,*) SUNSOL(IATM)
 1140    CONTINUE
      ENDIF
C     READ SURFACE EMISSIVITY
      REQLNE='SURFACE EMISSIVITIES'
      CALL GETCLN(FP,CTLLNE)
      IF (CTLLNE .NE. REQLNE) THEN
         WRITE (*,1010) REQLNE, CTLLNE
         STOP
      ELSE
         DO 1150 IATM=1,NOJOBS
            READ (FP,*) SUREMI(IATM)
 1150    CONTINUE
      ENDIF

      CLOSE (UNIT=FP)

C     SOME CONSISTENCY CHECKS
      IF (CLCTBR .EQ. 0) THEN
         IF (OT_RTB .EQ. 2) THEN
            WRITE (*,*) 'WARNING: SWITCHING BRIGHTNESS TEMPERATURE'
            WRITE (*,*) '         CALCULATION "ON" because the offset-'
            WRITE (*,*) '         tuning refers to brightness temp.'
            CLCTBR=1
         ENDIF
         IF (JACTYP .EQ. 2) THEN
            IF ((CLCJAC(1) .GT. 0) .OR. (CLCJAC(2) .GT. 0) .OR. 
     &          (CLCJAC(3) .GT. 0)) THEN
               WRITE (*,*) 'WARNING: SWITCHING BRIGHTNESS TEMPERATURE'
               WRITE (*,*) '         CALCULATION "ON" because the Jaco-'
               WRITE (*,*) '         bians are derivat. of brigh.temp.'
               CLCTBR=1
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END
C     END OF SUBROUTINE RDCTL

C ====================================================================

      SUBROUTINE GETATM(ATMFLN,TEMP,FIXED,WATER,OZONE,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     read both the atmosphere for which the tau shall be calculated
C     and the reference atmosphere. The pressure is read from the ref. prof.!
C
C     INPUT:
C       ATMFLN          atmosphere file name
C                       if ATMFLN .EQ. 'REFERENCE ATMOSPHERE' then 
C                         read reference atmosphere
C
C     OUTPUT:
C       all other parameter, explanation see main program
C
C     COMMENT:
C       this subroutine reads the atmosphere file format as used
C       for the creation of the predictors/coefficients
C
C     HISTORY:
C       written 10/14/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     SUBROUTINE ARGUMENT PARAMETER
      CHARACTER*80   ATMFLN

      REAL*8         TEMP(NLAYER),FIXED(NLAYER),
     &               WATER(NLAYER),OZONE(NLAYER)
      INTEGER        VERBOS
C     LOCAL VARIABLES
      INTEGER        FP,IERR,ILAY,IDUMMY
      REAL*8         RDMMY1,RDMMY2,RDMMY3
      CHARACTER*130  RDSTR
      
C     BEGIN
      ILAY=0
      FP=10
      OPEN(UNIT=FP,FILE=ATMFLN,STATUS='OLD',FORM='FORMATTED',
     &     IOSTAT=IERR)
      ! BEGIN OF READ-IN LOOP
 10   CONTINUE 
      READ(FP,15,END=20) RDSTR
 15   FORMAT(A130)
      IF ((RDSTR(1:1) .NE. '!') .AND. (RDSTR(1:1) .NE. 'r')) THEN
         ILAY=ILAY+1
         ! The user can easily change the requirements of the
         ! atmosphere input files by changing this subroutine,
         ! for example deleting the first 4 dummy parameters
         ! in the next command line
         READ(RDSTR,*) IDUMMY,RDMMY1,RDMMY2,RDMMY3,TEMP(ILAY),
     &                 FIXED(ILAY),WATER(ILAY),OZONE(ILAY)
      ENDIF
      GO TO 10
 20   CONTINUE
      ! END OF READ-IN LOOP
      IF (VERBOS .EQ. 1) THEN
         WRITE (*,21) ILAY
 21      FORMAT(' read in ',I4,' atm layer from atmosphere file')
      ENDIF
      CLOSE(FP)

      IF (ILAY .NE. NLAYER) THEN
         WRITE (*,*) 'FATAL ERROR:'
         WRITE (*,*) 'DID NOT READ ENOUGH DATA FROM ATMOSPHERE FILE'
         WRITE (*,22) ATMFLN
 22      FORMAT(' CHECK YOUR ATMOSPHERE INPUT FILE ',A80)
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF
      RETURN
      END
C     END OF SUBROUTINE GETATM

C ====================================================================

      SUBROUTINE ROUT(CTLFIL,TOTIRF,ID_IRF,RAD,TBR,CW_M,IJOB,VERBOS,
     &                CLCTBR)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Write calculated radiances to output file
C       The output file name will be constructed from the 
C       controlfile name
C
C     INPUT:
C       CTLFIL      controlfile name
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       RAD         radiances
C       TBR         brightness temperatures
C       CW_M        mean center wave numbers, as used for the
C                   calculation of TBR
C       IJOB        KERNEL loop number (used for output
C                   file name construction)
C       VERBOS      verbose (screen output yes/no)
C       CLCTBR      brightness temp. calculated or not (1=yes, 0=no)
C
C     OUTPUT
C       output is a file with the following content:
C       1st line: 
C         TOTIRF
C       2-(TOTIRF+1)-th line: 
C         ID_IRF(IIRF) RAD(ID_IRF(IIRF)) TBR(ID_IRF(IIRF)) CW_M(ID_IRF(IIRF))
C
C     CALLING ROUNTINE
C       KERNEL
C
C     HISTORY:
C       written 10/15/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
      INCLUDE "hffp_aux_dec.f"

C     CALLING VARIABLES
      CHARACTER*80   CTLFIL              ! name of controlfile
      INTEGER        TOTIRF,             ! total number of channels
     &               IJOB,               ! forward calculations loop counter
     &               ID_IRF(MAXIRF),     ! channel IDs
     &               VERBOS,             ! verbose (1=yes,0=no)
     &               CLCTBR              ! T_bright calc (1=yes,0=no)
      REAL*8         RAD(MAXIRF),        ! radiances
     &               TBR(MAXIRF),        ! brightness temperature
     &               CW_M(MAXIRF)        ! mean center wave numbers, as used
                                         ! for calculating TBR
C     LOCAL VARIABLES
      CHARACTER*80   OUTFLN,CDUMMY,EXTNUM
      INTEGER        FP,IIRF

C     BEGIN

      WRITE (EXTNUM,1) IJOB   ! this will fail if IJOB is larger than
 1    FORMAT(I10)             ! possible in the FORMAT statement
      ! remove leading blanks in EXTNUM
 2    IF (EXTNUM(1:1) .EQ. ' ') THEN 
         EXTNUM(1:79)=EXTNUM(2:80)
         GOTO 2
      ENDIF
      CALL CATSTR('.rad. ',EXTNUM,OUTFLN) !NOTE: place blank after .rad
      CALL CATSTR(CTLFIL,OUTFLN,CDUMMY) 
      CALL CATSTR(OUTPTH,CDUMMY,OUTFLN)
   
      IF (VERBOS .EQ. 1) THEN
         WRITE (*,*) 'creating output file for radiances: '
         WRITE (*,*) OUTFLN
      ENDIF
      FP=10
      ! use STATUS='NEW' if you want to avoid over-writing old files
      OPEN(UNIT=FP,FILE=OUTFLN,FORM='FORMATTED')
      ! OPEN(UNIT=FP,FILE=OUTFLN,FORM='FORMATTED',STATUS='NEW')
      WRITE (FP,201) TOTIRF
      DO 100 IIRF=1,TOTIRF
         ! write to file
         IF (CLCTBR .EQ. 1) THEN
            WRITE (FP,200) ID_IRF(IIRF),
     &           RAD(ID_IRF(IIRF)),TBR(ID_IRF(IIRF)),
     &           CW_M(ID_IRF(IIRF)) 
         ELSE
            WRITE (FP,202) ID_IRF(IIRF),RAD(ID_IRF(IIRF)) 
         ENDIF
 100  CONTINUE

      IF (MATLAB .EQ. 1) THEN
         WRITE (*,*) 'RAD_fast= ...'
         WRITE (*,20) RAD(1),RAD(2),RAD(3),RAD(4),RAD(5),
     &        RAD(6),RAD(7),RAD(8),RAD(9),RAD(10),
     &        RAD(11),RAD(12),RAD(13),RAD(14),RAD(15),
     &        RAD(16),RAD(17),RAD(18),RAD(19)
         WRITE (*,*) 'BRIGHTNESS TEMPERATURES:'
         WRITE (*,20) TBR(1),TBR(2),TBR(3),TBR(4),TBR(5),
     &        TBR(6),TBR(7),TBR(8),TBR(9),TBR(10),
     &        TBR(11),TBR(12),TBR(13),TBR(14),TBR(15),
     &        TBR(16),TBR(17),TBR(18),TBR(19)
 20      FORMAT('[',F20.15,',',F20.15,',',F20.15,',',F20.15,',',
     &        F20.15,',',F20.15,',',F20.15,',',F20.15,',',F20.15,',',
     &        F20.15,',',F20.15,',',F20.15,',',F20.15,',',F20.15,',',
     &        F20.15,',',F20.15,',',F20.15,',',F20.15,',',F20.15,'];')      
      ENDIF
 200  FORMAT(I5,' ',E20.12,' ',E20.12,' ',F20.12)
 201  FORMAT(I5)
 202  FORMAT(I5,' ',E20.12,' -1.0 -1.0')
      
      CLOSE(UNIT=FP)
      END 
C     END OF SUBROUTINE ROUT


C ====================================================================

      SUBROUTINE WR_JAC(CTLFIL,TOTIRF,ID_IRF,CLCJAC,IJOB,
     &                  JACOBT,JACOBW,JACOBO,FDJ_T,FDJ_W,FDJ_O,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Write calculated Jacobians to output file
C       The output file name will be constructed from the 
C       controlfile name
C
C     INPUT:
C       CTLFIL      controlfile name
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       CLCJAC      switches for calculating Jacobians
C       IJOB        KERNEL loop number (used for output
C                   file name construction)
C       JACOBT      ANALYTIC TEMP. JACOBIANS
C       JACOBW      ANALYTIC WATER JACOBIANS
C       JACOBO      ANALYTIC OZONE JACOBIANS
C       FDJ_T       FINITE DIFFERENCES TEMP. JACOBIANS
C       FDJ_W       FINITE DIFFERENCES WATER JACOBIANS
C       FDJ_O       FINITE DIFFERENCES OZONE JACOBIANS
C       VERBOS      VERBOSE PARAMETER (SCREEN OUTPUT ON/OFF)
C
C     OUTPUT
C       output are files with the following names
C       1st line: 
C         TOTIRF
C       2-(TOTIRF+1)-th line: 
C         ID_IRF(IIRF) RAD(ID_IRF(IIRF)) TBR(ID_IRF(IIRF)) CW_M(ID_IRF(IIRF))
C
C     CALLING ROUNTINE
C       main program
C
C     HISTORY:
C       written 10/15/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
      INCLUDE "hffp_aux_dec.f"

C     CALLING VARIABLES
      CHARACTER*80  CTLFIL
      INTEGER       TOTIRF,ID_IRF(MAXIRF),CLCJAC(3),IJOB,VERBOS
      REAL*8        JACOBT(MAXIRF,NLAYER),JACOBW(MAXIRF,NLAYER),
     &              JACOBO(MAXIRF,NLAYER),FDJ_T(MAXIRF,NLAYER),
     &              FDJ_W(MAXIRF,NLAYER),FDJ_O(MAXIRF,NLAYER)

C     LOCAL VARIABLES
      CHARACTER*80  OUTFLN,EXTNUM,CDUMMY

C     BEGIN
      WRITE (EXTNUM,1) IJOB  ! this will fail if IJOB is larger than
 1    FORMAT(I10)            ! possible in the FORMAT statement
      ! remove leading blanks in EXTNUM
 2    IF (EXTNUM(1:1) .EQ. ' ') THEN 
         EXTNUM(1:79)=EXTNUM(2:80)
         GOTO 2
      ENDIF

C     ANALYTIC JACOBIANS
      ! temperature Jacobians
      IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
         CALL CATSTR('.ajt. ',EXTNUM,OUTFLN) !NOTE: place blank after .ajt
         CALL CATSTR(CTLFIL,OUTFLN,CDUMMY) 
         CALL CATSTR(OUTPTH,CDUMMY,OUTFLN)
         CALL WRJACF(OUTFLN,TOTIRF,ID_IRF,JACOBT,VERBOS)
      endif
      ! water Jacobians
      IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
         CALL CATSTR('.ajw. ',EXTNUM,OUTFLN) !NOTE: place blank after .ajw
         CALL CATSTR(CTLFIL,OUTFLN,CDUMMY) 
         CALL CATSTR(OUTPTH,CDUMMY,OUTFLN)
         CALL WRJACF(OUTFLN,TOTIRF,ID_IRF,JACOBW,VERBOS)
      endif
       ! ozone Jacobians
      IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
         CALL CATSTR('.ajo. ',EXTNUM,OUTFLN) !NOTE: place blank after .ajo
         CALL CATSTR(CTLFIL,OUTFLN,CDUMMY) 
         CALL CATSTR(OUTPTH,CDUMMY,OUTFLN)
         CALL WRJACF(OUTFLN,TOTIRF,ID_IRF,JACOBO,VERBOS)
      endif
C     FINITE DIFFERENCES JACOBIANS
      ! temperature Jacobians
      IF ((CLCJAC(1) .EQ. 2) .OR. (CLCJAC(1) .EQ. 3)) THEN
         CALL CATSTR('.fjt. ',EXTNUM,OUTFLN) !NOTE: place blank after .fjt
         CALL CATSTR(CTLFIL,OUTFLN,CDUMMY) 
         CALL CATSTR(OUTPTH,CDUMMY,OUTFLN)
         CALL WRJACF(OUTFLN,TOTIRF,ID_IRF,FDJ_T,VERBOS)
      endif
      ! water Jacobians
      IF ((CLCJAC(2) .EQ. 2) .OR. (CLCJAC(2) .EQ. 3)) THEN
         CALL CATSTR('.fjw. ',EXTNUM,OUTFLN) !NOTE: place blank after .fjw
         CALL CATSTR(CTLFIL,OUTFLN,CDUMMY) 
         CALL CATSTR(OUTPTH,CDUMMY,OUTFLN)
         CALL WRJACF(OUTFLN,TOTIRF,ID_IRF,FDJ_W,VERBOS)
      endif
       ! ozone Jacobians
      IF ((CLCJAC(3) .EQ. 2) .OR. (CLCJAC(3) .EQ. 3)) THEN
         CALL CATSTR('.fjo. ',EXTNUM,OUTFLN) !NOTE: place blank after .fjo
         CALL CATSTR(CTLFIL,OUTFLN,CDUMMY) 
         CALL CATSTR(OUTPTH,CDUMMY,OUTFLN)
         CALL WRJACF(OUTFLN,TOTIRF,ID_IRF,FDJ_O,VERBOS)
      endif
        

      RETURN
      END
C     END OF SUBROUTINE WR_JAC


C ====================================================================

      SUBROUTINE WRJACF(OUTFLN,TOTIRF,ID_IRF,JAC,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       writes Jacobian matrix to file
C
C     INPUT
C       OUTFLN     file name
C       TOTIRF     total number of channels
C       ID_IRF     channel IDs
C       JAC        Jacobians
C
C     OUTPUT
C       output is a file of the struckture:
C
C     HISTORY
C       written 12/21/1997 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
      
C     CALLING VARIABLES
      CHARACTER*80   OUTFLN
      INTEGER        TOTIRF,ID_IRF(MAXIRF),VERBOS
      REAL*8         JAC(MAXIRF,NLAYER)

C     LOCAL VARIABLES
      INTEGER        FP,IIRF,ILAY,IIRF_A

C     BEGIN

      IF (VERBOS .EQ. 1) THEN
         WRITE (*,99) OUTFLN
 99      FORMAT(' creating Jacobian output file ',A80)
      ENDIF

c     ! only for debugging
c     write (*,*) 'Jacobians: instead ',OUTFLN
c     write (*,*) 'use filename:'
c     read (*,12) OUTFLN
c12   format(a80)

      FP=10
      OPEN(UNIT=FP,FILE=OUTFLN,FORM='FORMATTED')
      WRITE (FP,1) TOTIRF
 1    FORMAT('totirf=',I2,';')
      WRITE (FP,2) TOTIRF
 2    FORMAT('id_irf=zeros(1,',I2,');')
      WRITE (FP,3) TOTIRF,NLAYER
 3    FORMAT('jac=zeros(',I2,',',I3,');')
      DO 10 IIRF=1,TOTIRF
         IIRF_A=ID_IRF(IIRF)
         DO 20 ILAY=1,NLAYER
            WRITE (FP,4) IIRF,ILAY,JAC(IIRF_A,ILAY)
 4          FORMAT('jac(',I2,',',I3,')=',E30.23,';')
 20      CONTINUE
 10   CONTINUE
      CLOSE(UNIT=FP)
      RETURN
      END
C     END OF SUBROUTINE WRJACF
