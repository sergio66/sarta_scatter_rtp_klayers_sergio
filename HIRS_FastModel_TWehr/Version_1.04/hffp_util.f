C FILE NAME: hffp_util.f
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
C     utilities for HFFP
C ====================================================================

C ====================================================================

      SUBROUTINE KERNEL(IJOB,SATNUM,TOTIRF,ID_IRF,
     &               PREDID,NCOEFF,COEFF,T_TAUF,T_TAUW,T_TAUO,T_OFFS,
     &               T_SEC,
     &               SURF_T,SURF_P,SECANT,SUNSEC,SUREMI,SUNSOL,
     &               SOLRAD,ATEMP,AFIXED,AWATER,AOZONE,CLCJAC,
     &               CW_M,RADTOT,TBMEAN,JACOBT,JACOBW,JACOBO,VERBOS,
     &               NOTECT,JACTYP,CLCTBR,OT_RTB,T_ADJ)
      IMPLICIT NONE
      include "hffp_glob_dec.f"
      include "hffp_ref_dec.f"
      include "hffp_aux_dec.f"
      include "hffp_kernel_dec.f"
C --------------------------------------------------------------------
C
C     PURPOSE:
C       this is the main "forward" loop for HFFP
C
C     INPUT/OUTPUT: see include file hffp_kernel_dec.f
C
C     HISTORY:
C       written   01/14/98 by Tobias Wehr
C --------------------------------------------------------------------

C     LOCAL PARAMETER (PART A)
      REAL*8   CW_IRF(MAXIRF,NLAYER+1), ! FILTER CHANNEL CENTER WAVENO.
     &         TAU(MAXIRF,NLAYER), ! LAYER TRANSMITTANCES
     &         TAU_RS(MAXIRF,NLAYER), ! LAYER TRANSMITTANCES REFL. SOLAR
     &         PLANCK(MAXIRF,NLAYER+1), ! PLANCK FCT. FOR EACH CHANNEL AND
                                        ! LAYER. (NLAYER+1) REFERS TO THE
                                        ! SURFACE VALUE, NOT TO AN ATMOSPHERE 
                                        ! LAYER!!!
                                        ! for calculation of TBMEAN
                                        ! (NLAYER+1) REFERS TO THE SURFACE
     &         RAD_UP(MAXIRF),          ! RADIANCES UPWELLING THERMAL
     &         RAD_RS(MAXIRF),          ! RADIANCES REFLECTED SOLAR
     &         RAD_RT(MAXIRF),          ! RADIANCES REFLECTED THERMAL
     &         ALPHA              ! MULTIPLIER OF K OF THE SURFACE LAYER
      INTEGER  ISURF,             ! INDEX OF THE SURFACE LAYER
     &         K_NEG(MAXIRF,MAXMOL,NLAYER),  ! indicator if K is negative
                                             !   (see SUBROUTINE C_K)
     &         K_NEGS(MAXIRF,MAXMOL,NLAYER)  ! (same for solar)
       ! NOTE: the derivatives of the transmittances (variables DTAUT,
       !       DTAUW,DTAUO,DTAUTS,DTAUWS,DTAUOS) will be calculated as 
       !       derivatives of layer transmittances first and changed 
       !       afterward into derivatives of layer-to-space transmittances
      REAL*8   DTAUT(MAXIRF,NLAYER,NLAYER),  ! dtau/dT, thermal upw.
     &         DTAUW(MAXIRF,NLAYER,NLAYER),  ! dtau/dW, thermal upw.
     &         DTAUO(MAXIRF,NLAYER,NLAYER),  ! dtau/dO, thermal upw.
     &         DTAUTS(MAXIRF,NLAYER,NLAYER), ! dtau/dT, refl. solar
     &         DTAUWS(MAXIRF,NLAYER,NLAYER), ! dtau/dW, refl. solar
     &         DTAUOS(MAXIRF,NLAYER,NLAYER), ! dtau/dO, refl. solar
     &         DPLNCK(MAXIRF,NLAYER+1)       ! DERIVATIVE OF PLANCK FCT.
C     ! THE FOLLOWING 6 PARAMETERS ARE ONLY USED FOR DEBUGGING (COMPARISON
C     ! OF ANALYTICALLY DIFFERENCIATED TAU AND FINITE DIFFERENCES
!     REAL*8   FTAUT(MAXIRF,NLAYER,NLAYER),  ! dtau/dT, thermal upw.
!    &         FTAUW(MAXIRF,NLAYER,NLAYER),  ! dtau/dW, thermal upw.
!    &         FTAUO(MAXIRF,NLAYER,NLAYER),  ! dtau/dO, thermal upw.
!    &         FTAUTS(MAXIRF,NLAYER,NLAYER), ! dtau/dT, refl. solar
!    &         FTAUWS(MAXIRF,NLAYER,NLAYER), ! dtau/dW, refl. solar
!    &         FTAUOS(MAXIRF,NLAYER,NLAYER)  ! dtau/dO, refl. solar

C     LOCAL PARAMETER (PART B)
      REAL*8   JACUPT(MAXIRF,NLAYER), ! d(rad)/dT for thermal upwelling
     &         JACUPW(MAXIRF,NLAYER), ! d(rad)/dW for thermal upwelling
     &         JACUPO(MAXIRF,NLAYER), ! d(rad)/dO for thermal upwelling
     &         JACRST(MAXIRF,NLAYER), ! d(rad)/dT for refl. solar
     &         JACRSW(MAXIRF,NLAYER), ! d(rad)/dW for refl. solar
     &         JACRSO(MAXIRF,NLAYER)  ! d(rad)/dO for refl. solar
      INTEGER  IIRFH,      ! channel number counter ("IIRF HELP")
     &         ILAYH       ! layer number counter ("ILAY HELP")

C     BEGIN
      IF (VERBOS .EQ. 1) THEN
         WRITE (*,*) 'RUN  ',HFFPVS
      ENDIF

C     CALCULATE TRANSMITTANCES
      CALL C_TAU(IJOB,SATNUM,TOTIRF,ID_IRF,PREDID,
     &     NCOEFF,COEFF,SECANT,SUNSEC,
     &     T_TAUF,T_TAUW,T_TAUO,T_SEC,TAU,TAU_RS,
     &     ATEMP,AFIXED,AWATER,AOZONE,RTEMP,RFIXED,RWATER,ROZONE,
     &     PRES,K_NEG,K_NEGS,VERBOS,T_ADJ,NOTECT)
      ! NOTE: channel N refers to TAU(N,*)
      ! for the not-desired channels M, TAU(M,*) is not defined

C     CALULATE PARAMETERS OF LOWEST LAYER
      ! find out about the surface layer
      CALL LOWLAY(SURF_P,ISURF,ALPHA)
         
C     SURFACE TEMPERATURE: if the surface temperature is <0.0 in
      ! the control file, it will be set to the surface layer
      ! temperature
      IF (SURF_T .LT. 0.0) THEN
         SURF_T = ATEMP(ISURF)
         IF (VERBOS .EQ. 1) THEN
            WRITE (*,*) 'set T(surface) to ',SURF_T
         ENDIF
      ENDIF

C     CALCULATE CENTER WAVENUMBERS OF THE CHANNELS
      ! NOTE: CW_IRF will have only defined values for the
      !       desired channels, all other values are undefined.
      !       Indecees: channel N refers to CW_IRF(N,*)
      !       CW_IRF(*,NLAYER+1) refers to the surface temperature
      ! Now calculate center wavenumbers for the desired channels:
      ! use GETCWN (use GETCWS only for test runs)
      CALL GETCWN(SATNUM,CFFPTH,CFFILE,TOTIRF,
     &     ID_IRF,CW_IRF,ATEMP,SURF_T,VERBOS)

C     CHANGE SURFACE LAYER TAU(*,ISURF)
      ! NOTE: ALPHA is a scaling factor for the absorption coefficent k,
      !       i.e., since TAU=exp(-k), ALPHA is a power of TAU !
      DO 100 IIRFH=1,TOTIRF
         TAU(ID_IRF(IIRFH),ISURF)=TAU(ID_IRF(IIRFH),ISURF)**ALPHA
         TAU_RS(ID_IRF(IIRFH),ISURF)=TAU_RS(ID_IRF(IIRFH),ISURF)
     &        **ALPHA
 100  CONTINUE
      ! NOTE: all TAU below the surface layer are left unchanged
      !       in order to save computation time. So, be careful
      !       with the radiative transfer code, that you
      !       do not use the layer below the surface layer !!!

C     CALCULATE PLANCK FUNCTIONS FOR EACH LAYER AND CHANNEL
      ! NOTE: PLANCK(*,NLAYER+1) is the surface value
      CALL CPLNCK(TOTIRF,ID_IRF,ATEMP,CW_IRF,ISURF,
     &     SURF_T,PLANCK)

C     =================================================================
C       DO NOT DELETE THIS PART OF THE CODE SINCE IT IS USED
C       IN THE PROCESS OF FITTING THE REFLECTED THERMAL APPROXIMATION
C       FOR NEW INSTRUMENTS (E.G., NOAA-15++)
C       (comment it out if not needed)
C     =================================================================
c     ! output for matlab:
c     ! "QUICK HACK" USED FOR FITTING REFLECTED THERMAL
c     ! NOTE: run it for a control file with all 48 profiles, no Jacobians
c     IF (IJOB .EQ. 1) WRITE (*,*) 'planck=zeros(19,48,100);'
c     DO 101 IIRFH=1,19
c        DO 102 ILAYH=1,100
c           WRITE (*,103) IIRFH,IJOB,ILAYH,PLANCK(IIRFH,ILAYH)
c102     CONTINUE
c101  CONTINUE
c103  FORMAT('planck(',I2,',',I2,',',I3,')=',E25.18,';')


C     CALCULATE UPWELLING THERMAL RADIANCES
      ! NOTE: PLANCK(*,NLAYER+1) is the surface value
      CALL RADUPW(TOTIRF,ID_IRF,ISURF,TAU,SUREMI,
     &     PLANCK,RAD_UP)

C     CALCULATE REFLECTED SOLAR RADIANCES
      CALL SOLRFL(TOTIRF,ID_IRF,ISURF,TAU_RS,SUREMI,
     &     SUNSOL,SUNSEC,SOLRAD,RAD_RS)

C     CALCULATE REFLECTED THERMAL
      ! approximated for emissivity=0.975 and surface pressure approx. 1 atm
      CALL THRFL(SATNUM,TOTIRF,ID_IRF,ISURF,PLANCK,TAU,
     &           SECANT,T_SEC,SUREMI,RAD_RT)
C     ! make a test printout in matlab format:
C     WRITE (*,*) 'radrt=zeros(1,19);'
C     DO 301 IIRFH=1,TOTIRF
C        WRITE (*,302) IIRFH,RAD_RT(IIRFH)
C301  CONTINUE
C302  FORMAT('radrt(',I2,')=',E25.18,';')

C     ADD UP RADIANCES
      DO 300 IIRFH=1,TOTIRF
         RADTOT(ID_IRF(IIRFH))=RAD_UP(ID_IRF(IIRFH))+
     &        RAD_RS(ID_IRF(IIRFH)) + RAD_RT(ID_IRF(IIRFH))
 300  CONTINUE

      IF (OT_RTB .EQ. 1) THEN
C        TUNE RADIANCES (see also "IF (OT_RTB .EQ. 2)" below)
         CALL TUNETB(TOTIRF,ID_IRF,T_OFFS,RADTOT)
      ENDIF

      IF (CLCTBR .EQ. 1) THEN
C        CALCULATE BRIGHTNESS TEMPERATURES FOR THE
C        MEAN CENTER WAVENUMBER FOR EACH CHANNEL
         CALL C_TBR(TOTIRF,ID_IRF,CW_IRF,RADTOT,TBMEAN,CW_M)
      ENDIF

      IF (OT_RTB .EQ. 2) THEN
C        TUNE BRIGHTNESS TEMP. (see also "IF (OT_RTB .EQ. 1)" above)
         CALL TUNETB(TOTIRF,ID_IRF,T_OFFS,TBMEAN)
      ENDIF

C     -------------------------------------------------------
C     ANALYTIC JACOBIANS:
C     -------------------------------------------------------
      IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3) .OR.
     &    (CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3) .OR.
     &    (CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
         IF (VERBOS .EQ. 1) THEN
            IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
               WRITE (*,*) 'CALC ANALYTIC JACOBIANS dY/dT'
            ENDIF
            IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
               WRITE (*,*) 'CALC ANALYTIC JACOBIANS dY/dW'
            ENDIF
            IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
               WRITE (*,*) 'CALC ANALYTIC JACOBIANS dY/dO'
            ENDIF
         ENDIF
         IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
            ! calculate deriva tives of Planck function
            ! NOTE: DPLNCK(*,NLAYER+1) is the surface value
            IF (VERBOS .EQ. 1) THEN
               WRITE (*,*) 'CALCULATE d(Planck)/dT'
            ENDIF

            ! calculate d(Planck)/dT
            ! either: use CDPLCK for analytical derivatives
            ! or: use CDPLCK_FD for finite differences (slowly!)
            ! both have the same I/O data
            ! (Rem: the finite diff. subroutine was developed
            !  to check the analytical subroutine)
            CALL CDPLCK(TOTIRF,ID_IRF,ATEMP,CW_IRF,ISURF,
     &           SURF_T,DPLNCK)
c           CALL CDPLCK_FD(TOTIRF,ID_IRF,ATEMP,CW_IRF,ISURF,
c    &           SURF_T,DPLNCK)
            ! for plotting use output routine:
c           CALL WRITE_DPLANCK(TOTIRF,ID_IRF,DPLNCK)

         ENDIF
         ! calculate derivatives of layer transmittances
         ! for thermal upwelling and reflected solar
         CALL C_DTAU(CLCJAC,SATNUM,TOTIRF,ID_IRF,PREDID,
     &        NCOEFF,COEFF,SECANT,T_SEC,SUNSEC,
     &        TAU,TAU_RS,
     &        ATEMP,AFIXED,AWATER,AOZONE,
     &        RTEMP,RFIXED,RWATER,ROZONE,PRES,
     &        DTAUT,DTAUW,DTAUO,DTAUTS,DTAUWS,DTAUOS,
     &        K_NEG,K_NEGS,VERBOS)
C     ALTERNATIVELY, C_DTAU_FD CAN BE USED TO CALCULATE D(TAU)/DX
C     USING FINITE DIFFERNCES INSTEAD OF ANALYTICAL DERIVATIVES:
!        CALL C_DTAU_FD(CLCJAC,SATNUM,TOTIRF,ID_IRF,PREDID,
!    &        NCOEFF,COEFF,SECANT,T_SEC,SUNSEC,
!    &        ATEMP,AFIXED,AWATER,AOZONE,
!    &        RTEMP,RFIXED,RWATER,ROZONE,PRES,
!    &        DTAUT,DTAUW,DTAUO,DTAUTS,DTAUWS,DTAUOS,
!    &        VERBOS,IJOB,T_TAUF,T_TAUW,T_TAUO)
C     USE THE FOLLOWING FOR COMPARISON (WRITE TO FTAUX INSTEAD DTAUX)
!        CALL C_DTAU_FD(CLCJAC,SATNUM,TOTIRF,ID_IRF,PREDID,
!    &        NCOEFF,COEFF,SECANT,T_SEC,SUNSEC,
!    &        ATEMP,AFIXED,AWATER,AOZONE,
!    &        RTEMP,RFIXED,RWATER,ROZONE,PRES,
!    &        FTAUT,FTAUW,FTAUO,FTAUTS,FTAUWS,FTAUOS,
!    &        VERBOS,IJOB,T_TAUF,T_TAUW,T_TAUO)
!        ! write DTAUT,DTAUW,DTAUO to file (output for debugging)
!        CALL WRITE_DTAUX(TOTIRF,ID_IRF,
!    &        DTAUT,DTAUW,DTAUO,DTAUTS,DTAUWS,DTAUOS,
!    &        FTAUT,FTAUW,FTAUO,FTAUTS,FTAUWS,FTAUOS)


         ! CHANGE SURFACE LAYER DTAUs 
         ! NOTE: the d(tau)/dX are calculated for all 100 layer, not
         !       considering, that the surface layer index might be > 1.
         !       Therefore, all d(tau)/dX below the surface are garbage!
         DO 400 IIRFH=1,TOTIRF
            DO 401 ILAYH=NLAYER,ISURF,-1
               DTAUT(ID_IRF(IIRFH),ISURF,ILAYH)=ALPHA*
     &              DTAUT(ID_IRF(IIRFH),ISURF,ILAYH)
               DTAUW(ID_IRF(IIRFH),ISURF,ILAYH)=ALPHA*
     &              DTAUW(ID_IRF(IIRFH),ISURF,ILAYH)
               DTAUO(ID_IRF(IIRFH),ISURF,ILAYH)=ALPHA*
     &              DTAUO(ID_IRF(IIRFH),ISURF,ILAYH)
               DTAUTS(ID_IRF(IIRFH),ISURF,ILAYH)=ALPHA*
     &              DTAUTS(ID_IRF(IIRFH),ISURF,ILAYH)
               DTAUWS(ID_IRF(IIRFH),ISURF,ILAYH)=ALPHA*
     &              DTAUWS(ID_IRF(IIRFH),ISURF,ILAYH)
               DTAUOS(ID_IRF(IIRFH),ISURF,ILAYH)=ALPHA*
     &              DTAUOS(ID_IRF(IIRFH),ISURF,ILAYH)
 401        CONTINUE
 400     CONTINUE

         ! calculate derivatives of layer-to-space transmittances
         ! for the thermal upwelling
         ! NOTE: DTAUT,DTAUW,DTAUO will be changed from derivatives
         !       of layer transmittances to derivatives of layer-
         !       to-space transmittances!
         CALL C_DLST(CLCJAC,TOTIRF,ID_IRF,TAU,DTAUT,DTAUW,DTAUO)
         ! calculate derivatives of layer-to-space transmittances
         ! for the reflected solar
         ! NOTE: DTAUTS,DTAUWS,DTAUOS will be changed from derivatives
         !       of layer transmittances to derivatives of layer-
         !       to-space transmittances!
         CALL C_DLST(CLCJAC,TOTIRF,ID_IRF,TAU_RS,DTAUTS,DTAUWS,
     &        DTAUOS)
         ! NOTE: DTAUTS,DTAUWS,DTAUOS are now deriv. of layer-to-sp. tr.
         ! calculate the Jacobians for upwelling thermal radiation
         IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
            CALL C_J_UP(TOTIRF,ID_IRF,TAU,DTAUT,PLANCK,DPLNCK,
     &           JACUPT,1,ISURF,SUREMI)
         ENDIF
         IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
            CALL C_J_UP(TOTIRF,ID_IRF,TAU,DTAUW,PLANCK,DPLNCK,
     &           JACUPW,2,ISURF,SUREMI)
         ENDIF
         IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
            CALL C_J_UP(TOTIRF,ID_IRF,TAU,DTAUO,PLANCK,DPLNCK,
     &           JACUPO,2,ISURF,SUREMI)
         ENDIF
         ! calculate Jacobians for reflected solar
         IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
            CALL C_J_RS(TOTIRF,ID_IRF,ISURF,DTAUTS,SUREMI,
     &                  SUNSOL,SUNSEC,SOLRAD,JACRST)
         ENDIF
         IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
            CALL C_J_RS(TOTIRF,ID_IRF,ISURF,DTAUWS,SUREMI,
     &                  SUNSOL,SUNSEC,SOLRAD,JACRSW)
         ENDIF
         IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
            CALL C_J_RS(TOTIRF,ID_IRF,ISURF,DTAUOS,SUREMI,
     &                  SUNSOL,SUNSEC,SOLRAD,JACRSO)
         ENDIF
         ! (not yet coded: calculate Jacobians for reflected thermal)
         ! add up Jacobians 
         IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
            DO 501 IIRFH=1,TOTIRF
               DO 502 ILAYH=1,NLAYER
                  JACOBT(ID_IRF(IIRFH),ILAYH)=
     &                 JACUPT(ID_IRF(IIRFH),ILAYH) +
     &                 JACRST(ID_IRF(IIRFH),ILAYH)
 502           CONTINUE
 501        CONTINUE
         ENDIF
         IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
            DO 503 IIRFH=1,TOTIRF
               DO 504 ILAYH=1,NLAYER
                  JACOBW(ID_IRF(IIRFH),ILAYH)=
     &                 JACUPW(ID_IRF(IIRFH),ILAYH) +
     &                 JACRSW(ID_IRF(IIRFH),ILAYH)
 504           CONTINUE
 503        CONTINUE
         ENDIF
         IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
            DO 505 IIRFH=1,TOTIRF
               DO 506 ILAYH=1,NLAYER
                  JACOBO(ID_IRF(IIRFH),ILAYH)=
     &                 JACUPO(ID_IRF(IIRFH),ILAYH) +
     &                 JACRSO(ID_IRF(IIRFH),ILAYH)
 506           CONTINUE
 505        CONTINUE
         ENDIF

C        Until here the Jacobians are radiance Jacobians.
C        If desired, convert now to bright.temp. Jacobians.
C        NOTE: JACOBT, JACOBW, JACOBO will be changed!
         IF (JACTYP .EQ. 2) THEN
            IF ((CLCJAC(1) .EQ. 1) .OR. (CLCJAC(1) .EQ. 3)) THEN
               CALL C_JTBR(TOTIRF,ID_IRF,CW_M,RADTOT,JACOBT)
            ENDIF
            IF ((CLCJAC(2) .EQ. 1) .OR. (CLCJAC(2) .EQ. 3)) THEN
               CALL C_JTBR(TOTIRF,ID_IRF,CW_M,RADTOT,JACOBW)
            ENDIF
            IF ((CLCJAC(3) .EQ. 1) .OR. (CLCJAC(3) .EQ. 3)) THEN
               CALL C_JTBR(TOTIRF,ID_IRF,CW_M,RADTOT,JACOBO)
            ENDIF
         ENDIF

      ENDIF ! end of analytic Jacobians calculations

      RETURN
      END 
C     END OF SUBROUTINE KERNEL

C=====================================================================

      SUBROUTINE FNDBLK(CASTR,IFIND)
C     find first nonblank character in a string (written by Sergio)
      IMPLICIT NONE

      CHARACTER*80 CASTR
      INTEGER IFIND

      IFIND=1
 10   CONTINUE
      IF ((CASTR(IFIND:IFIND) .NE. ' ') .AND. (IFIND .LT. 80)) THEN
        IFIND=IFIND+1
        GOTO 10
        END IF

      RETURN
      END

C ====================================================================

      SUBROUTINE CATSTR(STR1,STR2,SUMSTR)
C     COMBINE STRINGS STR1 AND STR2 TO SUMSTR, E.G.:
C     STR1='file.' STR2='dat' =>  SUMSTR='file.dat'
C     TRAILING BLANKS IN STR1 ARE REMOVED
      IMPLICIT NONE

      CHARACTER*80 STR1,STR2,SUMSTR      ! ALL LIMITED IN LENGTH!!!
      INTEGER      LNGTH1,LNGTH2,ISTRLG,I

C     BEGIN
C     FIND 1ST BLANKS     
      CALL FNDBLK(STR1,LNGTH1)
      LNGTH1=LNGTH1-1
      CALL FNDBLK(STR2,LNGTH2)
      LNGTH2=LNGTH2-1

      ISTRLG=LNGTH1+LNGTH2
      IF (ISTRLG .GT. 80) THEN
         WRITE (*,*) 'ERROR IN SUBROUTINE CATSTR'
         WRITE (*,*) '  STR1="',STR1,'"'
         WRITE (*,*) '  STR2="',STR2,'"'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ELSE
         SUMSTR(1:LNGTH1)             =STR1(1:LNGTH1)
         SUMSTR((LNGTH1+1):(ISTRLG+1))=STR2(1:LNGTH2)
         IF (ISTRLG .LT. 80) THEN
            DO 10 I=(ISTRLG+1),80
               SUMSTR(I:I)=' '
 10         CONTINUE
         ENDIF
      ENDIF

      RETURN
      END

C ====================================================================

      SUBROUTINE SATTOI(SATNUM,ISAT)
C --------------------------------------------------------------------
C
C     PURPOSE
C       get the index of NOAANR (from hffp_aux_dec.f) or return -1 if 
C       failed
C    
C     INPUT
C       SATNUM   satellite ID (e.g. 7 for NOAA7)
C
C     OUTPUT
C       ISAT     satellite index (index of NOAANR)
C
C     NOTE
C       this is VERY specific on NOAANR in hffp_aux_dec.f !
C
C     HISTORY
C       written 11/13/97 by T.Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_aux_dec.f"

      INTEGER  SATNUM,ISAT,I

      DO 100 I=1,7
         IF (SATNUM .EQ. NOAANR(I)) THEN
            ISAT=I
            GOTO 101
         ELSE
            ISAT=-1
         ENDIF
 100  CONTINUE

 101  IF (ISAT .EQ. -1) THEN
         WRITE (*,*) 'FATAL ERROR IN SUBROUTINE SATTOI'
         WRITE (*,*) '  ISAT NOT FOUND IN NOAANR'
         WRITE (*,*) '  PRGRAMMING ERROR'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF

      RETURN
      END ! OF SUBROUTINE SATTOI


C ====================================================================

      SUBROUTINE GETCLN(FP,CTLLNE)
C --------------------------------------------------------------------
C     
C     PURPOSE:
C       repeat reading one line from the file until the first
C       read charactar is not a comment sign "!"
C
C     INPUT:
C       FP      file pointer of an open file
C
C     OUTPUT:
C       CTLLNE  string variable
C
C     HISTORY:
C       written 11/5/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
C     CALLING PARAMETERS
      CHARACTER*40   CTLLNE
      INTEGER        FP
      
 100  FORMAT(A40)
      CTLLNE='! THIS IS A COMMENT LINE'
 200  IF (CTLLNE(1:1) .EQ. '!') THEN
         READ(FP,100) CTLLNE
         GOTO 200
      ENDIF
      
      RETURN
      END ! OF SUBROUTINE GETCLN


C ====================================================================

      SUBROUTINE RDPRID(SATNUM,PREDID,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     read the identifiers (IDs) of the predictors for all channels (IRFs)
C
C     INPUT:
C       SATNUM          satellite number
C       VERBOS          screen output on/off (1=on,0=off)
C
C     OUTPUT:
C       PREDID          IDs of the predictors (3-dim matrix)
C
C     COMMENT:
C       this subroutine reads all IDs whether the channel (IRF)
C       shall be calculated or not
C
C     HISTORY:
C       written 10/17/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
      INCLUDE "hffp_aux_dec.f"

C     SUBROUTINE ARGUMENT PARAMETER
      INTEGER        SATNUM,VERBOS
      INTEGER        PREDID(MAXIRF,MAXMOL,MAXPRD+1) ! ORDER OF MOL: F,W,O
C     LOCAL PARAMETER
      CHARACTER*80   FILENM           ! FILE NAME
      INTEGER        FP,IERR          ! FILE POINTER
      CHARACTER*80   CDMMY1
      INTEGER        IIRF,IDUMMY,IMOL,IPRD,N,I1,I2,I3,ISAT

C     initialize PREDID with -1 numbers
      DO 1 I1=1,MAXIRF
         DO 2 I2=1,MAXMOL
            DO 3 I3=1,(1+MAXPRD)
               PREDID(I1,I2,I3)=-1
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE

      ! get file name of predictor ID-file
      CALL SATTOI(SATNUM,ISAT)
      CDMMY1=PIDFLS(ISAT)
      CALL CATSTR(PCPATH,CDMMY1,FILENM)

      IF ((VERBOS .EQ. 1) .OR. (MATLAB .EQ. 1)) THEN
         WRITE (*,30) FILENM
 30      FORMAT(' read predictor ID-file ',A80)
      ENDIF
         
      FP=10
      OPEN(UNIT=FP,FILE=FILENM,STATUS='OLD',FORM='FORMATTED',
     &     IOSTAT=IERR)
      READ(FP,*) IDUMMY
      !write (*,*) IDUMMY
      IF (IDUMMY .NE. MAXIRF) THEN
         WRITE (*,*) 'FATAL ERROR IN RDPRID: IDUMMY .NE. MAXIRF'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF
      DO 300 IIRF=1,MAXIRF
         READ (FP,*) IDUMMY
         IF (IDUMMY .NE. IIRF) THEN
            WRITE (*,*) 'FATAL ERROR IN RDPRID: IDUMMY .NE. IIRF'
            WRITE (*,*) 'PROGRAM ABORTED'
            STOP
         ENDIF
         DO 200 IMOL=1,MAXMOL
            READ (FP,*) PREDID(IIRF,IMOL,1)
            BACKSPACE(UNIT=FP)  
            N=PREDID(IIRF,IMOL,1)
            READ (FP,*) IDUMMY, (PREDID(IIRF,IMOL,IPRD),IPRD=2,(N+1))
            !WRITE (*,*) (PREDID(IIRF,IMOL,IPRD),IPRD=1,N+1)
 200     CONTINUE
         READ (FP,*) IDUMMY
         IF (IDUMMY .NE. -1) THEN
            WRITE (*,*) 'FATAL ERROR IN RDPRID: IDUMMY .NE. -1'
            WRITE (*,*) 'PROGRAM ABORTED'
            STOP
         ENDIF
 300  CONTINUE
      CLOSE(FP)

      RETURN
      END


C ====================================================================

      SUBROUTINE RDCOEF(SATNUM,NCOEFF,COEFF,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     read the coefficients for all channels (IRFs)
C
C     INPUT:
C       SATNUM          satellite number
C       VERBOS          screen output on/off (1=on,0=off)
C
C     OUTPUT:
C       NCOEFF          number of coefficients
C       COEFF           coefficients (4-dim matrix)
C
C     COMMENT:
C       this subroutine reads all coefficients, whether the channel (IRF)
C       shall be calculated or not
C
C     HISTORY:
C       written 9/1/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
      INCLUDE "hffp_aux_dec.f"

C     SUBROUTINE ARGUMENT PARAMETER
      INTEGER        SATNUM,VERBOS
      INTEGER        NCOEFF(MAXIRF,MAXMOL) ! number of coefficients
      REAL*8         COEFF(MAXIRF,MAXMOL,MAXPRD,NLAYER) ! ORDER OF MOL: F,W,O
C     LOCAL PARAMETER
      CHARACTER*80   FILENM           ! FILE NAME
      INTEGER        FP,IERR          ! FILE POINTER
      CHARACTER*80   CDMMY1
      INTEGER        IIRF,IDUMMY,IMOL,ICOEF,ILAY,ISAT

      ! get file name for coefficents-file
      CALL SATTOI(SATNUM,ISAT)
      CDMMY1=COFFLS(ISAT)
      CALL CATSTR(PCPATH,CDMMY1,FILENM)

      IF ((VERBOS .EQ. 1) .OR. (MATLAB .EQ. 1)) THEN
         WRITE (*,30) FILENM
 30      FORMAT(' read coefficients-file ',A80)
      ENDIF

      FP=10
      OPEN(UNIT=FP,FILE=FILENM,STATUS='OLD',FORM='FORMATTED',
     &     IOSTAT=IERR)
      READ(FP,*) IDUMMY
      IF (IDUMMY .NE. MAXIRF) THEN
         WRITE (*,*) 'FATAL ERROR IN RDCOEF: IDUMMY .NE. MAXIRF'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF

      DO 800 IIRF=1,MAXIRF
         ! READ ON IRF-BLOCK
         READ (FP,*) IDUMMY
         IF (IDUMMY .NE. IIRF) THEN
            WRITE (*,*) 'FATAL ERROR IN RDCOEF: IDUMMY .NE. IIRF'
            WRITE (*,*) 'PROGRAM ABORTED'
            STOP
         ENDIF
         DO 700 IMOL=1,MAXMOL
            ! read number of coefficients for this IRF and MOL
            READ(FP,*) NCOEFF(IIRF,IMOL)
            IF (NCOEFF(IIRF,IMOL) .GT. 0) THEN
               ! read coefficients
               DO 600 ILAY=1,NLAYER
                  READ(FP,*) (COEFF(IIRF,IMOL,ICOEF,ILAY),
     &                        ICOEF=1,NCOEFF(IIRF,IMOL))
 600           CONTINUE
            ENDIF
            READ(FP,*) IDUMMY
            IF (IDUMMY .NE. -1) THEN
               WRITE (*,*) 'FATAL ERROR IN RDCOEF: IDUMMY .NE. -1'
               WRITE (*,*) 'PROGRAM ABORTED'
               STOP
            ENDIF
 700     CONTINUE
 800  CONTINUE
      
      CLOSE(FP)
      RETURN
      END

C ====================================================================

      SUBROUTINE RDSOLR(SATNUM,SOLRAD,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE
C       read the convolved solar radiances from file
C
C     INPUT
C       SATNUM     satellite number
C       VERBOS     screen output on/off (1=on,0=off)
C
C     OUTPUT
C       SOLRAD     vector of convolved solar radiances for 
C                  all channels
C
C     NOTES
C       1. this will read all convolved solar radiances, whether all
C          channels will be used later or not
C       2. def "convolved solar radiances": incoming solar radiance
C          convolved with the instrument response functions 
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
      INCLUDE "hffp_aux_dec.f"
C     SUBROUTINE CALLING VARIABLES
      INTEGER        SATNUM,VERBOS
      REAL*8         SOLRAD(MAXIRF)
C     LOCAL VARIABLES
      INTEGER        FP,ISAT,IIRF,IDUMMY
      CHARACTER*80   FILENM,CDMMY

C     BEGIN
      CALL SATTOI(SATNUM,ISAT)
      FP=10
      CDMMY=SOLRDF(ISAT)
      CALL CATSTR(CFFPTH,CDMMY,FILENM)
      IF (VERBOS .EQ. 1) THEN
         WRITE (*,100) FILENM
 100     FORMAT(' read convolved solar radiances from file ',A80)
      ENDIF

      OPEN(UNIT=FP,FILE=FILENM,FORM='FORMATTED',STATUS='OLD')
      READ(FP,*) IDUMMY
      IF (IDUMMY .NE. MAXIRF) THEN
         WRITE (*,*) 'FATAL ERROR IN RDSOLR: IDUMMY .NE. MAXIRF'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF
      DO 200 IIRF=1,MAXIRF
         READ(FP,*) SOLRAD(IIRF)
 200  CONTINUE
      CLOSE (UNIT=FP)
      !write (*,*) 'subroutine rdsolr: solrad=',solrad
      RETURN
      END ! OF SUBROUTINE RDSOLR

C ====================================================================

      SUBROUTINE GETCWN(SATNUM,CFFPTH,CFFILE,TOTIRF,ID_IRF,
     &                  CW_IRF,ATEMP,SURF_T,VERBOS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     interpolate actual filter center wavenumbers from tabulated
C     center wavenumbers for the actual temperature profile
C
C     INPUT:
C       SATNUM          satellite number (e.g. 7 for NOAA7)
C       CFFPTH          path for center wavenumbers files
C       CFFILE          center wavenumbers files
C       TOTIRF          total number of channels
C       ID_IRF          channel identifiers (channel numbers)
C       ATEMP           the atmosphere temperatures
C       SURF_T          the surface temperature
C       VERBOS          screen output on/off (1=on,0=off)
C
C     OUTPUT
C       CW_IRF          center wavenumbers as function of channel ID and layer
C                       the values for the surface temperature will be written
C                       to CW_IRF(*,NLAYER+1)
C
C     CALLING ROUNTINE
C       KERNEL
C
C     HISTORY:
C       written 10/6/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
C     SUBROUTINE CALLING VARIABLES
      CHARACTER*80 CFFPTH,CFFILE(7)   ! path and c.f. files names
      INTEGER SATNUM,                 ! satellite number (7,9,10,11,12,14)
     &        TOTIRF,ID_IRF(MAXIRF),  ! total no channels, channel IDs
     &        VERBOS                  ! verbose
      REAL*8  CW_IRF(MAXIRF,NLAYER+1),! center wavenumbers
     &        ATEMP(NLAYER),          ! temperature profile
     &        SURF_T                  ! surface temperature

C     LOCAL VARIABLES
      CHARACTER*80 CFFN               ! actually used file name
      INTEGER FP,IERR,                ! file pointer
     &        NOTMP,ITMP,IIRF,IXTMP   ! no of temperatures, index variables
      INTEGER MAXTTB                  ! maximum number of temperatures
      PARAMETER(MAXTTB=21)
      REAL*8  TTAB(MAXTTB),           ! tabel of T from centerfrq. file
     &        CFTAB(MAXIRF,MAXTTB)    ! tabel of centerwaven. for [channel,T]
      INTEGER ILAY                    ! layer index counter
      REAL*8  TDIFF1,TDIFF2,
     &        TLOCAL                  ! "local" temperature is either one of
                                      ! the layer temp.s or the surface temp.
      INTEGER OK
      
C     BEGIN
      IF (SATNUM .EQ. 7) THEN
         CALL CATSTR(CFFPTH,CFFILE(1),CFFN)
      ELSEIF (SATNUM .EQ. 9) THEN
         CALL CATSTR(CFFPTH,CFFILE(2),CFFN)
      ELSEIF (SATNUM .EQ. 10) THEN
         CALL CATSTR(CFFPTH,CFFILE(3),CFFN)
      ELSEIF (SATNUM .EQ. 11) THEN
         CALL CATSTR(CFFPTH,CFFILE(4),CFFN)
      ELSEIF (SATNUM .EQ. 12) THEN
         CALL CATSTR(CFFPTH,CFFILE(5),CFFN)
      ELSEIF (SATNUM .EQ. 14) THEN
         CALL CATSTR(CFFPTH,CFFILE(6),CFFN)
      ELSEIF (SATNUM .EQ. 15) THEN
         CALL CATSTR(CFFPTH,CFFILE(7),CFFN)
      ELSE
         WRITE (*,*) 'FATAL ERROR IN GETCFQ: INVALID SATNUM'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF

C     read table of center wavenumbers and temperatures
      FP=10
      OPEN(UNIT=FP,FILE=CFFN,FORM='FORMATTED',STATUS='OLD',
     &     IOSTAT=IERR)
C     read number of filter channels
      READ (FP,*) NOTMP ! NOTMP is misused here as the "number of channels"
      IF (NOTMP .NE. MAXIRF) THEN
         WRITE (*,*) 'FATAL ERROR IN SUBROUTINE GETCFQ:'
         WRITE (*,12) MAXIRF
 12      FORMAT(' NO OF FILTER CHANNELS IN C.FRQ. FILE NOT EQUAL ',I3)
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF
C     read number of temperatures
      READ (FP,*) NOTMP
      IF (NOTMP .GT. MAXTTB) THEN
         WRITE (*,*) 'FATAL ERROR IN SUBROUTINE GETCFQ:'
         WRITE (*,*) 'NUMBER OF TEMPERATURES LARGER THAN MAXTTB'
         WRITE (*,11) NOTMP
 11      FORMAT(' ADJUST PARAMETER IN SUBROUTINE GETCFQ TO MAXTTB=',I5)
         WRITE (*,*) 'RECOMPILE PROGRAM AND START AGAIN'
         STOP
      ENDIF

      DO 1010 ITMP=1,NOTMP
         READ (FP,*) TTAB(ITMP)
         !WRITE (*,*) TTAB(ITMP)
 1010 CONTINUE

      DO 1020 IIRF=1,MAXIRF
         DO 1030 ITMP=1,NOTMP
            READ (FP,*) CFTAB(IIRF,ITMP)
            !WRITE (*,*) CFTAB(IIRF,ITMP)
 1030    CONTINUE
 1020 CONTINUE
      CLOSE (UNIT=FP)       

C     INTERPOLATE THE CENTER WAVENUMBER

      ! NOTE: in order to save a little computer time, we could 
      ! introduce here a mean temperature instead of using a different
      ! temperature for each layer. This could be done for channel 1-18
      ! but not for channel 19

      ! LOOP ALL LAYERS
      IXTMP=1
      DO 1040 ILAY=1,NLAYER+1
         OK=0
         IF (ILAY .EQ. NLAYER+1) THEN
            TLOCAL=SURF_T
         ELSE
            TLOCAL=ATEMP(ILAY)
         ENDIF
         IF (TLOCAL .LT. TTAB(1)) THEN
            ! EXTRAPOLATE
            WRITE (*,*) 'WARNING: EXTRAPOLATING T FOR CENTER WAVENUMBER'
            ! WRITE (*,*) 'TLOCAL,TTAB(1)=',TLOCAL,TTAB(1)
            DO 1050 IIRF=1,TOTIRF
               CALL INTPCW(TTAB(1),TTAB(2),TLOCAL,
     &                     CFTAB(ID_IRF(IIRF),1),
     &                     CFTAB(ID_IRF(IIRF),2),
     &                     CW_IRF(ID_IRF(IIRF),ILAY))
 1050       CONTINUE
            OK=1
         ELSE
            ! TRY TO INTERPOLATE
 1060       CONTINUE ! this is a GOTO destination 
            IF ( (TLOCAL .GE. TTAB(IXTMP)) .AND. 
     &           (TLOCAL .LE. TTAB(IXTMP+1)) ) THEN
               DO 1070 IIRF=1,TOTIRF
                  CALL INTPCW(TTAB(IXTMP),TTAB(IXTMP+1),TLOCAL,
     &                        CFTAB(ID_IRF(IIRF),IXTMP),
     &                        CFTAB(ID_IRF(IIRF),IXTMP+1),
     &                        CW_IRF(ID_IRF(IIRF),ILAY))
 1070          CONTINUE
               OK=1
            ELSE
               IF (IXTMP .EQ. 1) THEN
                  IXTMP=IXTMP+1
                  GOTO 1060
               ELSE
                  ! calculate gradient for next search
                  TDIFF1=ABS(TLOCAL-TTAB(IXTMP+1))
                  TDIFF2=ABS(TLOCAL-TTAB(IXTMP-1))
                  IF (TDIFF1 .LE. TDIFF2) THEN
                     IXTMP=IXTMP+1
                  ELSE
                     IXTMP=IXTMP-1
                  ENDIF
                  IF (IXTMP .LT. NOTMP) THEN
                     GOTO 1060
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         ! if attempt of extrapolation to lower T and interpolation has failed:
         IF (OK .NE. 1) THEN
            ! TLOCAL IS LARGER THAN TABEL VALUES:
            IF (TLOCAL .GT. TTAB(NOTMP)) THEN
               WRITE (*,*) 'WARNING: EXTRAPOLATING T FOR CENTER WAVENO.'
               DO 1080 IIRF=1,TOTIRF
                  CALL INTPCW(TTAB(NOTMP-1),TTAB(NOTMP),TLOCAL,
     &                 CFTAB(ID_IRF(IIRF),NOTMP-1),
     &                 CFTAB(ID_IRF(IIRF),NOTMP),
     &                 CW_IRF(ID_IRF(IIRF),ILAY))
 1080          CONTINUE
               OK=1
            ELSE
               WRITE (*,*) 'FATAL SOFTWARE ERROR IN GETCFQ:'
               WRITE (*,*) 'PROGRAM ABORTED'
               STOP
            ENDIF
         ENDIF

 1040 CONTINUE
         
      END ! of SUBROUTINE GETCFQ


C ====================================================================

      SUBROUTINE INTPCW(T1,T2,TN,F1,F2,FN)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     interpolates filter channel center wavenumber FN for the 
C     temperature TN from the temperatures T1,T2 and wavenumbers F1,F2
C     Interpolation: approximately it is: (a=constant)
C                    center wavenumber = a * exp(-1/temperature)
C                    because of the dependency of the center wavenumber on 
C                    the Planck function
C     Therefore the wavenumber is linear interpolated with respect to
C     THETA=exp(-1/T)
C
C     INPUT:
C       T1,T2        temperatures
C       F1,F2        center waveno. for the temperatures T1,T2
C       TN           temp. for which the center waveno. shall be interpolated
C
C     OUTPUT
C       FN           interpolated center waveno. at temperature TN
C
C     CALLING ROUNTINE
C       GETCWN
C
C     HISTORY:
C       written 10/9/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8  T1,T2,TN,F1,F2,FN,TH1,TH2,THN

      ! calculate THETAs
      TH1=DEXP(-1/T1)
      TH2=DEXP(-1/T2)
      THN=DEXP(-1/TN)
      ! interpolate
      FN=F1+(F2-F1)*((THN-TH1)/(TH2-TH1))

      END ! of SUBROUTINE INTPCW

C ====================================================================

      SUBROUTINE LOWLAY(SURF_P,ISURF,ALPHA)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     Def.: we call the layer, in which the surface is located the
C           "surface layer"
C     1. find the index of the surface layer (ISURF)
C        so that the first tau to be considered in the radiative transfer
C        will be TAU(*,ISURF)
C     2. find the absorption coefficient (k) scaling factor of the 
C        surface layer DSURF, with which the k of the surface layer has 
C        to be multiplied.
C        Usually, it is 0 < ALPHA <= 1 ,but if the surface pressure 
C        is larger than the largest pressure in the reference profile. 
C        Then, ALPHA will be extrapolated and is > 1.
C
C        The layer transmittance (TAU) of the surface layer is then
C        TAU(*,ISURF)**ALPHA
C
C        ASSUMPTION:
C        We assume the that the absorptino coefficient k is a linear 
C        function of the pressure p:
C          x p = k = -log (tau)
C
C
C     INPUT:
C       SURF_P       surface pressure in mbar (same unit as PRESBD)
C
C     OUTPUT
C       ISURF        surface layer index
C       ALPHA        surface layer transmittance scaling factor
C
C     CALLING ROUNTINE
C       GETCWN
C
C     HISTORY:
C       written 10/14/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
      INCLUDE "hffp_ref_dec.f"

C     SUBROUTINE CALLING PARAMETERS
      REAL*8      SURF_P,               ! surface pressure
     &            ALPHA                 ! surface layer abs.coeff. scaling factor
      INTEGER     ISURF                 ! surface layer index
C     LOCAL VARIABLES
      REAL*8      ALPHA1,ALPHA2
      PARAMETER(ALPHA1=1.0)
      PARAMETER(ALPHA2=0.0)
      
      ISURF=1
      IF (SURF_P .GT. PRESBD(1)) THEN
         WRITE (*,*) 'WARNING FROM SUBROUTINE LOWLAY:'
         WRITE (*,*) 'SURFACE PRESSURE LARGER THAN MODEL PROFILE'
         WRITE (*,*) '  SURF_P=',SURF_P,' PRESBD(1)=',PRESBD(1)
         WRITE (*,*) '  EXTRAPOLATING SURFACE LAYER ABS.COEF.'
         ! extrapolate DSURF
         ISURF=1
         CALL INTLIN(PRESBD(1),PRESBD(2),SURF_P,ALPHA1,ALPHA2,ALPHA)
      ELSE
 100     CONTINUE ! this is a GOTO destination
         IF ( (SURF_P .LE. PRESBD(ISURF)) .AND. 
     &        (SURF_P .GT. PRESBD(ISURF+1)) ) THEN
            CALL INTLIN(PRESBD(ISURF),PRESBD(ISURF+1),
     &                  SURF_P,ALPHA1,ALPHA2,ALPHA)
         ELSE
            ISURF=ISURF+1
            IF (ISURF .GT. NLAYER) THEN
               WRITE(*,*)'FATAL ERROR IN SUBROUTIN LOWLAY:'
               WRITE(*,*)'SOMETHING IS WRONG WITH YOUR SURFACE PRESSURE'
               WRITE(*,*)'IT IS SMALLER THAN THE LOWEST REFERENCE MODEL'
               WRITE(*,*)'PRESSURE! YOUR SURFACE PRESSURE IS ',SURF_P
               WRITE(*,*)'PROGRAM ABORTED'
               STOP
            ENDIF
            GOTO 100
         ENDIF
      ENDIF

      END ! of SUBROUTINE LOWLAY


C ====================================================================

      SUBROUTINE INTLIN(X1,X2,XN,Y1,Y2,YN)
C --------------------------------------------------------------------
C
C     PURPOSE:
C     linear interpolation
C
C     INPUT:
C       X1,X2        tabulated independend values
C       XN           independend value for which YN shall be interpolated
C       Y1,Y2        tabulated dependend values
C
C     OUTPUT
C       YN           interpolated value YN(XN)
C
C     HISTORY:
C       written 10/14/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8  X1,X2,XN,Y1,Y2,YN

      ! interpolate
      YN=Y1+(Y2-Y1)*((X1-XN)/(X1-X2))

      END ! of SUBROUTINE INTPCW

C ====================================================================

      SUBROUTINE CPLNCK(TOTIRF,ID_IRF,ATEMP,CW_IRF,ISURF,SURF_T,PLANCK)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Calculate the Planck function for every channel and every
C       layer starting at layer ISURF
C
C     INPUT:
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       ATEMP       layer temperatures
C       CW_IRF      center wavenumbers of the filter channels for each layer,
C                   CW_IRF(*,NLAYER+1) after the surface temperature 
C       ISURF       surface layer index
C       SURF_T      surface temperature
C
C     OUTPUT
C       PLANCK      value of Planck function for each channel and layer
C                   NOTE: PLANCK(*,NLAYER+1) contains the surface value
C
C     CALLING ROUNTINE
C       KERNEL
C
C     NOTE:
C       1. CW_IRF is a matrix of MAXIRF by NLAYER, i.e., every channel
C          has a different center wavenumber in every layer
C          CW_IRF(*,NLAYER+1) is the cw value for the surface temperature
C          ISURF is used in order to calculate Planck only for the surface
C          layer and the layers above. (Saves calculation time.)
C       2. the Planck function value for the surface temperature SURF_T
C          will be written to PLANCK(*,NLAYER+1)
C
C     HISTORY:
C       written 10/14/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     CALLING VARIABLES
      INTEGER     TOTIRF,ID_IRF(MAXIRF),  ! # of channels, channel IDs
     &            ISURF                   ! surface layer index
      REAL*8      ATEMP(NLAYER),          ! temperature profile
     &            SURF_T,                 ! surface temperature
     &            CW_IRF(MAXIRF,NLAYER+1),! channel center wavenumbers
     &            PLANCK(MAXIRF,NLAYER+1) ! Planck function values
C     LOCAL VARIABLES
      INTEGER     IIRF,ILAY
C     NATURAL CONSTANTS
      REAL*8 C1,C2 ! Planck function constants in wavenumber-units
      DATA C1/1.1911E-8/  ! C1 constant of Planck function (wavenumber units)
      DATA C2/1.4387863/  ! C2 constant of Planck function (wavenumber units)
        ! Planck function:
        ! B = C1*WAVENUMBER**3 / (EXP(C2*WAVENUMBER/TEMPERATURE) -1)
        ! units:
        ! [C1]= W m-2 sr-1 (cm-1)-4  
        ! [C2]= K (cm-1)-1
        ! [B] = W m-2 sr-1 (cm-1)-1

C     BEGIN
      DO 100 ILAY=ISURF,NLAYER+1 ! NLAYER+1 refers to the surface value
         DO 200 IIRF=1,TOTIRF
            ! calculate Planck function value
            IF (ILAY .EQ. NLAYER+1) THEN
               ! for the surface temperature
               PLANCK(ID_IRF(IIRF),ILAY)=C1*CW_IRF(ID_IRF(IIRF),ILAY)**3 
     &               / ( DEXP(C2*CW_IRF(ID_IRF(IIRF),ILAY)/SURF_T) - 1 )
            ELSE
               ! for the atmospheric temperatures
               PLANCK(ID_IRF(IIRF),ILAY)=C1*CW_IRF(ID_IRF(IIRF),ILAY)**3 
     &            / (DEXP(C2*CW_IRF(ID_IRF(IIRF),ILAY)/ATEMP(ILAY)) - 1)
            ENDIF
 200     CONTINUE       
 100  CONTINUE

      END ! of SUBROUTINE CPLNCK

C ====================================================================

      SUBROUTINE RADUPW(TOTIRF,ID_IRF,ISURF,TAU,SUREMI,PLANCK,RAD_UP)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Radiative transfer calculation for the upwelling radiation
C
C     INPUT:
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       ISURF       surface layer index
C       TAU         layer transmittances for each channel and layer
C       SUREMI      surface emissivity of this atmosphere (= for this job)
C       PLANCK      planck values for each channel and layer
C
C     OUTPUT
C       RAD_UP      upwelling radiance for each channel and layer
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

C     CALLING VARIABLES
      INTEGER    TOTIRF,ID_IRF(MAXIRF),  ! # of channels, channel IDs
     &           ISURF                   ! surface layer index
      REAL*8     TAU(MAXIRF,NLAYER),     ! layer transmittances
     &           SUREMI,                 ! surface emissivity
     &           PLANCK(MAXIRF,NLAYER+1),! Planck function values
     &           RAD_UP(MAXIRF)          ! upwelling radiances
C     LOCAL VARIABLES
      INTEGER    IIRF,ILAY
      REAL*8     TAU_LS,                 ! layer-to-space transmittance
     &           RAD                     ! radiance

C     BEGIN
      DO 100 IIRF=1,TOTIRF
         TAU_LS=1.
         RAD=0.
         DO 200 ILAY=NLAYER,ISURF,-1
            ! CALCULATE ATMOSPHERE EMISSION AND TRANSFER
            ! Note: TAU_LS is the layer-to-space transmittance of ILAY+1
            RAD=RAD+(1-TAU(ID_IRF(IIRF),ILAY))
     &              *PLANCK(ID_IRF(IIRF),ILAY)*TAU_LS
            TAU_LS=TAU_LS*TAU(ID_IRF(IIRF),ILAY)
 200     CONTINUE
         ! ADD SURFACE EMISSION AND TRANSFER
         RAD_UP(ID_IRF(IIRF))=RAD + SUREMI*
     &                        PLANCK(ID_IRF(IIRF),NLAYER+1)*TAU_LS
 100  CONTINUE
      END ! OF SUBROUTINE RADUPW

C ====================================================================

      SUBROUTINE SOLRFL(TOTIRF,ID_IRF,ISURF,TAU_RS,SUREMI,
     &                  SUNSOL,SUNSEC,SOLRAD,RAD_RS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Radiative transfer calculation for the reflected solar radiation,
C       which is an additional term for the total radiation
C
C     INPUT:
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       ISURF       surface layer index
C       SUREMI      surface emissivity
C       SUNSOL      sun solid angle
C       SUNSEC      sun secant
C       SOLRAD      convolved solar radiance (this is the incoming
C                   radiation into the atmosphere, convolved with
C                   the instrument response functions)
C       TAU_RS      refleced solar transmittances
C
C     OUTPUT
C       RAD_RS      reflected solar radiances
C
C     CALLING ROUNTINE
C       KERNEL
C
C     HISTORY:
C       written 11/14/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     CALLING VARIABLES
      INTEGER    TOTIRF,ID_IRF(MAXIRF),  ! # of channels, channel IDs
     &           ISURF                   ! surface layer index
      REAL*8     RAD_RS(MAXIRF),         ! reflected solar radiance
     &           TAU_RS(MAXIRF,NLAYER),  ! reflected solar transmittance
     &           SUREMI,                 ! surface emissivity
     &           SUNSOL,                 ! solid angle of the sun
     &           SUNSEC,                 ! sun secant
     &           SOLRAD(MAXIRF)          ! convolved incoming solar radiance
C     LOCAL VARIABLES
      INTEGER    IIRF,ILAY
      REAL*8     TAU_LS,       ! layer-to-space transmittance (surface layer)
     &           PI
      PARAMETER(PI=3.1415926535897931)
C     BEGIN
      DO 100 IIRF=1,TOTIRF
         IF (ID_IRF(IIRF) .GE. SNIIRF) THEN
            ! calculate layer-to-space transmittance for the surface layer
            TAU_LS=1.0
            DO 200 ILAY=ISURF,NLAYER
               TAU_LS=TAU_LS*TAU_RS(ID_IRF(IIRF),ILAY)
 200        CONTINUE
            ! calculate solar reflected radiance
            RAD_RS(ID_IRF(IIRF))=(1-SUREMI)*SUNSOL*SOLRAD(ID_IRF(IIRF))
     &           *TAU_LS/(PI*SUNSEC)
         ELSE
            RAD_RS(ID_IRF(IIRF))=0.0
         ENDIF
 100  CONTINUE

      RETURN
      END 
C     END OF SUBROUTINE SOLRFL

C ====================================================================

      SUBROUTINE THRFL(SATNUM,TOTIRF,ID_IRF,ISURF,PLANCK,TAU,
     &                 SECANT,T_SEC,SUREMI,RAD_RT)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       calculation of the reflected thermal radiation
C       which is an additional term to the total radiation
C       Note: the refl. thermal radiation is calculated after
C             model calculations for 48 sample profiles with
C             a surface reflectivity of 0.975 and a surface
C             pressure of about standard surface pressure
C
C     INPUT:
C       SATNUM      satellite number (one of: 7,9,10,11,...)
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       ISURF       surface layer index
C       PLANCK      Planck function
C       TAU         layer transmittances for each channel and layer
C       SECANT      secant
C       T_SEC       channel dependend secant offset (tuning parameter)
C
C     OUTPUT
C       RAD_RT      reflected solar radiances
C
C     CALLING ROUNTINE
C       KERNEL
C
C     HISTORY:
C       written 01/25/98 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
      INCLUDE "hffp_rt_dec.f"

C     CALLING VARIABLES
      INTEGER  SATNUM,                 ! satellite number (7,9,10,11,...)
     &         TOTIRF,ID_IRF(MAXIRF),  ! # of channels, channel IDs
     &         ISURF                   ! surface layer index
      REAL*8   PLANCK(MAXIRF,NLAYER+1),! Planck function
     &         TAU(MAXIRF,NLAYER),     ! layer transmittances
     &         SECANT,                 ! secant
     &         T_SEC(MAXIRF),          ! channel depend. secant offset
     &         SUREMI,                 ! surface emissivity
     &         RAD_RT(MAXIRF)          ! reflected thermal radiance

C     VARIABLES DECLARED IN THE INCLUDED FILE hffp_rt_dec.f:
C     (see included file to make it sure!)
C     ! integer of the layer used for the Planck(of layer X)-predictors:
C     INTEGER IPLNCK(19,N)  ! 19 channels, N satellites (noaa 7,9,10,11,...)
C     REAL*8  AC(5,19,N)    ! 5 predictors, 19 channels, N satellites
C     INTEGER IAC           ! help variable for assigning AC data

C     LOCAL VARIABELS
      INTEGER   IIRF,IIRF_A,ISAT,ILAY
      REAL*8    FI,         ! F_i (see latex doc)
     &          TAUSS,      ! surface-to-space transmittance
     &          C_SEC       ! = SECANT + T_SEC, "tuned" secant

C     BEGIN
      CALL SATTOI(SATNUM,ISAT)
      DO 100 IIRF=1,TOTIRF
         IIRF_A=ID_IRF(IIRF)
         C_SEC=SECANT+T_SEC(IIRF_A)
         IF ((IIRF_A .EQ. 1) .OR. (IIRF_A .EQ. 2)) THEN
            RAD_RT(IIRF_A)=0.0
         ELSE
            ! create FI
            FI = AC(1,IIRF_A,ISAT) +
     &           AC(2,IIRF_A,ISAT) / C_SEC +
     &           AC(3,IIRF_A,ISAT) * 
     &              PLANCK(IIRF_A,IPLNCK(IIRF_A,ISAT)) +
     &           AC(4,IIRF_A,ISAT) * 
     &              PLANCK(IIRF_A,IPLNCK(IIRF_A,ISAT)) / C_SEC +
     &           AC(5,IIRF_A,ISAT) * PLANCK(IIRF_A,ISURF) / 
     &              PLANCK(IIRF_A,IPLNCK(IIRF_A,ISAT))
            ! calculate TAUSS (surface-to-space transmittance)
            TAUSS=1.0
            DO 200 ILAY=NLAYER,ISURF,-1
               TAUSS=TAUSS*TAU(IIRF_A,ILAY)
 200        CONTINUE
            ! create RD
            RAD_RT(IIRF_A)=
     &           (1-SUREMI)*PLANCK(IIRF_A,ISURF)*TAUSS*(1-TAUSS)*FI
         ENDIF
 100  CONTINUE
      END
C     END OF SUBROUTINE THRFL

C ====================================================================

      SUBROUTINE C_J_RS(TOTIRF,ID_IRF,ISURF,DTAU,SUREMI,
     &                  SUNSOL,SUNSEC,SOLRAD,JACRS)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Calculate analytic jacobians (radiative transfer part) for 
C       reflected solar radiatin
C
C     INPUT:
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       ISURF       surface layer index
C       SUREMI      surface emissivity
C       SUNSOL      sun solid angle
C       SUNSEC      sun secant
C       SOLRAD      convolved solar radiance (this is the incoming
C                   radiation into the atmosphere, convolved with
C                   the instrument response functions)
C       DTAU        d(tau)/dX for refl. solar; X= T or W or O
C
C     OUTPUT
C       RAD_RS      reflected solar radiances
C
C     CALLING ROUNTINE
C       KERNEL
C
C     HISTORY:
C       written 11/14/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
C     CALLING VARIABLES
      INTEGER    TOTIRF,ID_IRF(MAXIRF),  ! # of channels, channel IDs
     &           ISURF                   ! surface layer index
      REAL*8     JACRS(MAXIRF,NLAYER),   ! d(reflected solar radiance)/dX
     &           DTAU(MAXIRF,NLAYER,NLAYER),! d(tau)/dX (X=T or W or O)
     &           SUREMI,                 ! surface emissivity
     &           SUNSOL,                 ! solid angle of the sun
     &           SUNSEC,                 ! sun secant
     &           SOLRAD(MAXIRF)          ! convolved incoming solar radiance
C     LOCAL VARIABLES
      INTEGER    IIRF,ILAY2
      REAL*8     PI
      PARAMETER(PI=3.1415926535897931)
C     BEGIN
 
      DO 100 IIRF=1,TOTIRF
         IF (ID_IRF(IIRF) .GE. SNIIRF) THEN
            DO 201 ILAY2=ISURF,NLAYER
               JACRS(ID_IRF(IIRF),ILAY2)=DTAU(ID_IRF(IIRF),ISURF,ILAY2)*
     &              (1-SUREMI)*SUNSOL*SOLRAD(ID_IRF(IIRF))/(PI*SUNSEC)
 201        CONTINUE
            IF (ISURF .GT. 1) THEN
               DO 202 ILAY2=1,(ISURF-1)
                  JACRS(ID_IRF(IIRF),ILAY2)=0.0
 202           CONTINUE
            ENDIF
         ELSE
            DO 300 ILAY2=1,NLAYER
               JACRS(ID_IRF(IIRF),ILAY2)=0.0
 300        CONTINUE
         ENDIF
 100  CONTINUE

      RETURN
      END
C     END OF C_J_RS

C ====================================================================

      SUBROUTINE C_TBR(TOTIRF,ID_IRF,CW_IRF,RAD,TBRIGT,CW_M)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Calculate the brightness temperatures from radiances, assuming
C       a mean (over all 100 layer) channel center wavenumber.
C
C     INPUT:
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       CW_IRF      filter channel center wavenumbers
C       RAD         radiances
C
C     OUTPUT
C       TBRIGT      brightness temperature for mean atmosph. temp and 
C                   mean filter channel center wavenumber
C       CW_M        mean channel center frequencies, as used for the 
C                   calculation of TBRIGT
C
C     CALLING ROUNTINE
C       KERNEL
C
C     HISTORY:
C       written 10/20/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     CALLING VARIABLES
      INTEGER  TOTIRF,                 ! total number of channels
     &         ID_IRF(MAXIRF)          ! channel IDs
      REAL*8   CW_IRF(MAXIRF,NLAYER+1),! center wavenumbers
     &         RAD(MAXIRF),            ! radiances
     &         TBRIGT(MAXIRF)          ! brightness temperatures
      REAL*8   CW_M(MAXIRF) ! mean center wavenumber [cm-1] for each channel
C     LOCAL VARIABLES
      INTEGER  ILAY,IIRF,
     &         CW_SET       ! 0 if CW_M not yet calculated, 1 if already calc.
C     NATURAL CONSTANTS
      REAL*8      C1,C2     ! constants for Planck function (wavenumber units)
      DATA C1/1.1911E-8/  ! C1 constant of Planck function (wavenumber units)
      DATA C2/1.4387863/  ! C2 constant of Planck function (wavenumber units)
        ! Planck function:
        ! B = C1*WAVENUMBER**3 / (EXP(C2*WAVENUMBER/TEMPERATURE) -1)
        ! units:
        ! [C1]= W m-2 sr-1 (cm-1)-4  
        ! [C2]= K (cm-1)-1
        ! [B] = W m-2 sr-1 (cm-1)-1
C     INITIALIZE
      DATA CW_SET/0/
C     SAVE STATEMENT
      SAVE CW_SET

C     BEGIN
      IF (CW_SET .EQ. 0) THEN
         DO 100 IIRF=1,TOTIRF
            CW_M(ID_IRF(IIRF))=0.
 100     CONTINUE
         ! use the mean center wave number of all layers for each channel
         DO 200 ILAY=1,NLAYER
            DO 300 IIRF=1,TOTIRF
               CW_M(ID_IRF(IIRF))=CW_M(ID_IRF(IIRF))
     &              +CW_IRF(ID_IRF(IIRF),ILAY)/NLAYER
 300        CONTINUE
 200     CONTINUE
         CW_SET=1
      ENDIF
  
      IF (MATLAB .EQ. 1) THEN
         WRITE (*,*) 'frequencies for brightness temperatures:' 
         WRITE (*,*) 'frq= ...'
         WRITE (*,301) CW_M(1),CW_M(2),CW_M(3),CW_M(4),CW_M(5),
     &        CW_M(6),CW_M(7),CW_M(8),CW_M(9),CW_M(10),
     &        CW_M(11),CW_M(12),CW_M(13),CW_M(14),CW_M(15),
     &        CW_M(16),CW_M(17),CW_M(18),CW_M(19)
 301     FORMAT('[',F15.10,',',F15.10,',',F15.10,',',F15.10,',',
     &        F15.10,',',F15.10,',',F15.10,',',F15.10,',',F15.10,',',
     &        F15.10,',',F15.10,',',F15.10,',',F15.10,',',F15.10,',',
     &        F15.10,',',F15.10,',',F15.10,',',F15.10,',',F15.10,'];')
      ENDIF

      DO 400 IIRF=1,TOTIRF
         TBRIGT(ID_IRF(IIRF)) =
     &      (C2*CW_M(ID_IRF(IIRF))) /
     &      DLOG( (C1*CW_M(ID_IRF(IIRF))**3)/RAD(ID_IRF(IIRF)) + 1)        
 400  CONTINUE

      END ! OF SUBROUTINE C_TBR

C ====================================================================

      SUBROUTINE C_JTBR(TOTIRF,ID_IRF,CW_M,RAD,JACOB)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       convert radiance-Jacobians to brightness temperature-Jacobians,
C       assuming a mean (over all 100 layer) channel center wavenumber.
C
C     INPUT:
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       CW_M        mean channel center frequencies
C       RAD         radiances
C
C     INPUT/OUTPUT
C       JACOB       as input:  radiance Jacobians 
C                   as output: brightness temperature Jacobians
C
C     CALLING ROUNTINE
C       KERNEL
C
C     HISTORY:
C       written 12/20/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     CALLING VARIABLES
      INTEGER  TOTIRF,        ! total number of channels
     &         ID_IRF(MAXIRF) ! channel IDs
      REAL*8   CW_M(MAXIRF),  ! mean center wavenumber [cm-1] for each channel
     &         RAD(MAXIRF),   ! radiances
     &         JACOB(MAXIRF,NLAYER) ! Jacobians
C     LOCAL VARIABLES
      INTEGER  ILAY,IIRF
      REAL*8   NU,DMMY1,DMMY2
C     NATURAL CONSTANTS
      REAL*8      C1,C2     ! constants for Planck function (wavenumber units)
      DATA C1/1.1911E-8/  ! C1 constant of Planck function (wavenumber units)
      DATA C2/1.4387863/  ! C2 constant of Planck function (wavenumber units)
        ! Planck function:
        ! B = C1*WAVENUMBER**3 / (EXP(C2*WAVENUMBER/TEMPERATURE) -1)
        ! units:
        ! [C1]= W m-2 sr-1 (cm-1)-4  
        ! [C2]= K (cm-1)-1
        ! [B] = W m-2 sr-1 (cm-1)-1

C     BEGIN
      DO 100 IIRF=1,TOTIRF
         NU=CW_M(ID_IRF(IIRF)) ! for programmer's convenience
         DMMY1=C1*(NU**3)/RAD(ID_IRF(IIRF)) + 1.0
         DMMY2=((DLOG(DMMY1))**2) * DMMY1 * RAD(ID_IRF(IIRF))**2
         DMMY1=C1*C2*(NU**4)/DMMY2
         DO 200 ILAY=1,NLAYER
            JACOB(ID_IRF(IIRF),ILAY)=DMMY1*JACOB(ID_IRF(IIRF),ILAY)
 200     CONTINUE
 100  CONTINUE

      END ! OF SUBROUTINE C_JTBR



C ====================================================================

      SUBROUTINE GETCWS(SATNUM,CFFPTH,CFFILE,TOTIRF,ID_IRF,
     &                  CW_IRF,ATEMP,SURF_T,VERBOS)
C --------------------------------------------------------------------
C
C SAME AS GETCWN, BUT USE CENTER WAVE NUMBERS AS USED IN THE
C MATLAB PROGRAMS
C
C This routine is only used for debugging and consistency checks
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
C     SUBROUTINE CALLING VARIABLES
      CHARACTER*80 CFFPTH,CFFILE(7)   ! path and c.f. files names
      INTEGER SATNUM,                 ! satellite number (7,9,10,11,...)
     &        TOTIRF,ID_IRF(MAXIRF),  ! total no channels, channel IDs
     &        VERBOS                  ! verbose parameter 
      REAL*8  CW_IRF(MAXIRF,NLAYER+1),! center wavenumbers
     &        ATEMP(NLAYER),          ! temperature profile
     &        SURF_T                  ! surface temperature
C     LOCAL VARIABLES
      INTEGER I

      IF ((VERBOS .EQ. 1) .OR. (MATLAB .EQ. 1)) THEN
         WRITE (*,1) SATNUM
 1       FORMAT(' satellite number: ',I3)
      ENDIF
      IF (SATNUM .EQ. 7) THEN
         DO 100 I=1,NLAYER+1
            CW_IRF( 1,I) =  667.965015229724710
            CW_IRF( 2,I) =  678.775099981398850
            CW_IRF( 3,I) =  691.704633440290190
            CW_IRF( 4,I) =  704.529456728980680
            CW_IRF( 5,I) =  717.135056268601150
            CW_IRF( 6,I) =  733.657589867001430
            CW_IRF( 7,I) =  749.520100911458370
            CW_IRF( 8,I) =  897.008265904017890
            CW_IRF( 9,I) = 1025.681349981398900
            CW_IRF(10,I)= 1220.023646763392900
            CW_IRF(11,I)= 1362.740001860119000
            CW_IRF(12,I)= 1481.082141694568500
            CW_IRF(13,I)= 2181.978887648809600
            CW_IRF(14,I)= 2206.521740141369300
            CW_IRF(15,I)= 2239.589529854910800
            CW_IRF(16,I)= 2270.503034319196600
            CW_IRF(17,I)= 2356.928350539434600
            CW_IRF(18,I)= 2512.998523530505700
            CW_IRF(19,I)= 2645.690627325148900
 100     CONTINUE
      ELSEIF (SATNUM .EQ. 9) THEN
         DO 200 I=1,NLAYER+1
            CW_IRF( 1,I)=   667.571777343750000
            CW_IRF( 2,I)=   679.742803664434520
            CW_IRF( 3,I)=   691.052952357700860
            CW_IRF( 4,I)=   703.401602608816920
            CW_IRF( 5,I)=   717.137570335751430
            CW_IRF( 6,I)=   732.783054896763360
            CW_IRF( 7,I)=   749.512343633742600
            CW_IRF( 8,I)=   897.807245163690480
            CW_IRF( 9,I)=  1030.093447730654800
            CW_IRF(10,I)=  1218.675833565848300
            CW_IRF(11,I)=  1365.180832635788600
            CW_IRF(12,I)=  1479.900652204241200
            CW_IRF(13,I)=  2190.341924758184600
            CW_IRF(14,I)=  2209.358444940476100
            CW_IRF(15,I)=  2242.571207682291500
            CW_IRF(16,I)=  2273.603073846726100
            CW_IRF(17,I)=  2359.490815662202300
            CW_IRF(18,I)=  2516.095505487351100
            CW_IRF(19,I)=  2656.655680338541500
 200     CONTINUE
      ELSEIF (SATNUM .EQ. 10) THEN
         DO 300 I=1,NLAYER+1
            CW_IRF( 1,I)=   667.708054315476150
            CW_IRF( 2,I)=   680.140183221726150
            CW_IRF( 3,I)=   691.221729096912210
            CW_IRF( 4,I)=   704.382771809895870
            CW_IRF( 5,I)=   716.652413504464330
            CW_IRF( 6,I)=   732.825395856584810
            CW_IRF( 7,I)=   750.381958007812500
            CW_IRF( 8,I)=   898.240129743303560
            CW_IRF( 9,I)=  1027.364513578869000
            CW_IRF(10,I)=  1222.061715262276700
            CW_IRF(11,I)=  1363.226225353422700
            CW_IRF(12,I)=  1479.989972795758800
            CW_IRF(13,I)=  2189.303001767113100
            CW_IRF(14,I)=  2206.324404761904600
            CW_IRF(15,I)=  2239.139846075148900
            CW_IRF(16,I)=  2267.840413411458500
            CW_IRF(17,I)=  2358.712913876488100
            CW_IRF(18,I)=  2511.958821614583500
            CW_IRF(19,I)=  2653.427304222470400
 300     CONTINUE
      ELSEIF (SATNUM .EQ. 11) THEN
         DO 400 I=1,NLAYER+1
            CW_IRF( 1,I)=   668.958403087797590
            CW_IRF( 2,I)=   678.859151204427120
            CW_IRF( 3,I)=   689.697329566592320
            CW_IRF( 4,I)=   703.317574637276830
            CW_IRF( 5,I)=   716.793718610491060
            CW_IRF( 6,I)=   732.045494442894320
            CW_IRF( 7,I)=   749.344953264508940
            CW_IRF( 8,I)=   900.091849190848170
            CW_IRF( 9,I)=  1030.897553943452300
            CW_IRF(10,I)=   795.727047874813930
            CW_IRF(11,I)=  1360.570707775297700
            CW_IRF(12,I)=  1477.118710472470200
            CW_IRF(13,I)=  2189.655215308779600
            CW_IRF(14,I)=  2209.417433965773900
            CW_IRF(15,I)=  2238.981224423363100
            CW_IRF(16,I)=  2267.639195033482300
            CW_IRF(17,I)=  2416.028134300595400
            CW_IRF(18,I)=  2511.164515904017700
            CW_IRF(19,I)=  2657.778483072916500
 400     CONTINUE
      ELSEIF (SATNUM .EQ. 12) THEN
         DO 500 I=1,NLAYER+1
            CW_IRF( 1,I)=   667.572236560639910
            CW_IRF( 2,I)=   680.161417643229130
            CW_IRF( 3,I)=   689.959394182477690
            CW_IRF( 4,I)=   704.124407087053560
            CW_IRF( 5,I)=   716.244457426525290
            CW_IRF( 6,I)=   732.040254138764910
            CW_IRF( 7,I)=   751.663219633556540
            CW_IRF( 8,I)=   898.355927966889910
            CW_IRF( 9,I)=  1025.184256417410800
            CW_IRF(10,I)=  1219.729125976562500
            CW_IRF(11,I)=  1366.707327706473300
            CW_IRF(12,I)=  1477.133155459449400
            CW_IRF(13,I)=  2188.081752232142700
            CW_IRF(14,I)=  2210.308593750000000
            CW_IRF(15,I)=  2238.129987444196600
            CW_IRF(16,I)=  2267.342680431547700
            CW_IRF(17,I)=  2361.262230282738100
            CW_IRF(18,I)=  2513.613862537202300
            CW_IRF(19,I)=  2646.569533575148900
 500     CONTINUE
      ELSEIF (SATNUM .EQ. 14) THEN
         DO 600 I=1,NLAYER+1
            CW_IRF( 1,I)=   668.896164667038650
            CW_IRF( 2,I)=   679.374488467261930
            CW_IRF( 3,I)=   689.593218122209810
            CW_IRF( 4,I)=   703.592032296316920
            CW_IRF( 5,I)=   714.603288922991060
            CW_IRF( 6,I)=   732.176751999627980
            CW_IRF( 7,I)=   749.569205147879420
            CW_IRF( 8,I)=   898.301720028831820
            CW_IRF( 9,I)=  1028.048659551711400
            CW_IRF(10,I)=   795.956202915736640
            CW_IRF(11,I)=  1360.438563755580400
            CW_IRF(12,I)=  1478.257585797991200
            CW_IRF(13,I)=  2190.976946149553400
            CW_IRF(14,I)=  2207.041503906250000
            CW_IRF(15,I)=  2235.962646484375000
            CW_IRF(16,I)=  2267.851469494047700
            CW_IRF(17,I)=  2419.864641462053400
            CW_IRF(18,I)=  2511.431663876488100
            CW_IRF(19,I)=  2641.804815383184600
 600     CONTINUE
      ELSE
         WRITE (*,*) 'FATAL ERROR IN SUBROUTINE GETCWS'
         WRITE (*,*) 'SATNUM NOT VALID'
         WRITE (*,*) 'PROGRAM ABORTED'
         STOP
      ENDIF
      END ! of SUBROUTINE GETCFQ

C ====================================================================

      SUBROUTINE TUNETB(TOTIRF,ID_IRF,T_OFFS,RAD_TB)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Tunes calculated radiances or brightness temperatures with 
C       the offset-tuning parameters
C
C     INPUT:
C       TOTIRF    total number of channels
C       ID_IRF    channel IDs
C       T_OFFS    offset tuning parameters (either for rad or T_bright)
C
C     INPUT/OUTPUT
C       RAD_TB       radiance or brightness temperatures
C
C     CALLING ROUNTINE
C       KERNEL
C
C     HISTORY:
C       written 11/4/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     CALLING VARIABLES
      INTEGER   TOTIRF,         ! TOTAL NUMBER OF IRFs (IRF=channel)
     &          ID_IRF(MAXIRF)  ! IDs OF IRFs
      REAL*8    T_OFFS(MAXIRF), ! RAD. OR T_BRIGHT TUNING OFFSET
     &          RAD_TB(MAXIRF)  ! RADIANCES OR BRIGHTNESS TEMPERATURES
C     LOCAL VARIABLES
      INTEGER   IIRF

C     BEGIN
      DO 100 IIRF=1,TOTIRF
         RAD_TB(ID_IRF(IIRF))=RAD_TB(ID_IRF(IIRF))+T_OFFS(ID_IRF(IIRF))
 100  CONTINUE
      
      END ! OF SUBROUTINE TUNETB

C ====================================================================

      SUBROUTINE CDPLCK(TOTIRF,ID_IRF,ATEMP,CW_IRF,ISURF,SURF_T,DPLNCK)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Calculate the temperature derivative of the Planck function 
C       for every channel and every layer starting at layer ISURF
C       Methode: analytical
C
C     INPUT:
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       ATEMP       layer temperatures
C       CW_IRF      center wavenumbers of the filter channels for each layer,
C                   CW_IRF(*,NLAYER+1) after the surface temperature 
C       ISURF       surface layer index
C       SURF_T      surface temperature
C
C     OUTPUT
C       DPLNCK      value of Planck function for each channel and layer
C                   NOTE: PLANCK(*,NLAYER+1) contains the surface value
C
C     CALLING ROUNTINE
C       KERNEL
C
C     NOTE:
C       1. CW_IRF is a matrix of MAXIRF by NLAYER, i.e., every channel
C          has a different center wavenumber in every layer
C          CW_IRF(*,NLAYER+1) is the cw value for the surface temperature
C          ISURF is used in order to calculate Planck only for the surface
C          layer and the layers above. (Saves calculation time.)
C       2. the Planck function value for the surface temperature SURF_T
C          will be written to DPLNCK(*,NLAYER+1)
C
C     HISTORY:
C       written 12/02/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     CALLING VARIABLES
      INTEGER     TOTIRF,ID_IRF(MAXIRF),  ! # of channels, channel IDs
     &            ISURF                   ! surface layer index
      REAL*8      ATEMP(NLAYER),          ! temperature profile
     &            SURF_T,                 ! surface temperature
     &            CW_IRF(MAXIRF,NLAYER+1),! channel center wavenumbers
     &            DPLNCK(MAXIRF,NLAYER+1) ! Planck function values
C     LOCAL VARIABLES
      INTEGER     IIRF,ILAY
C     NATURAL CONSTANTS
      REAL*8 C1,C2 ! Planck function constants in wavenumber-units
      DATA C1/1.1911E-8/  ! C1 constant of Planck function (wavenumber units)
      DATA C2/1.4387863/  ! C2 constant of Planck function (wavenumber units)
        ! Planck function:
        ! B = C1*WAVENUMBER**3 / (EXP(C2*WAVENUMBER/TEMPERATURE) -1)
        ! temperature derivative of Planck function:
        ! dB/dT=C1*C2*WAVENUMBER**4 * EXP(C2*WAVENUMBER/TEMPERATURE)
        !       / ( TEMPERATURE*(EXP(C2*WAVENUMBER/TEMPERATURE)-1) )**2
        ! units:
        ! [C1]= W m-2 sr-1 (cm-1)-4  
        ! [C2]= K (cm-1)-1
        ! [B] = W m-2 sr-1 (cm-1)-1
C     BEGIN
      DO 100 ILAY=ISURF,NLAYER+1 ! NLAYER+1 refers to the surface value
         DO 200 IIRF=1,TOTIRF
            ! calculate Planck function value
            IF (ILAY .EQ. NLAYER+1) THEN
               ! for the surface temperature
               DPLNCK(ID_IRF(IIRF),ILAY)=
     &              C1*C2*CW_IRF(ID_IRF(IIRF),ILAY)**4
     &              *DEXP(C2*CW_IRF(ID_IRF(IIRF),ILAY)/SURF_T) /
     &              (SURF_T*
     &               (DEXP(C2*CW_IRF(ID_IRF(IIRF),ILAY)/SURF_T)-1)
     &              )**2
            ELSE
               ! for the atmospheric temperatures
               DPLNCK(ID_IRF(IIRF),ILAY)=
     &              C1*C2*CW_IRF(ID_IRF(IIRF),ILAY)**4
     &              *DEXP(C2*CW_IRF(ID_IRF(IIRF),ILAY)/ATEMP(ILAY)) /
     &              (ATEMP(ILAY)*
     &               (DEXP(C2*CW_IRF(ID_IRF(IIRF),ILAY)/ATEMP(ILAY))-1)
     &              )**2
            ENDIF
 200     CONTINUE       
 100  CONTINUE

      RETURN
      END 
C     END OF SUBROUTINE CPLNCK

C ====================================================================

      SUBROUTINE CDPLCK_FD(TOTIRF,ID_IRF,ATEMP,CW_IRF,ISURF,SURF_T,
     &                     DPLNCK)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Calculate the temperature derivative of the Planck function 
C       for every channel and every layer starting at layer ISURF
C       Methode: finite differences
C
C     INPUT:
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       ATEMP       layer temperatures
C       CW_IRF      center wavenumbers of the filter channels for each layer,
C                   CW_IRF(*,NLAYER+1) after the surface temperature 
C       ISURF       surface layer index
C       SURF_T      surface temperature
C
C     OUTPUT
C       DPLNCK      value of Planck function for each channel and layer
C                   NOTE: PLANCK(*,NLAYER+1) contains the surface value
C
C     CALLING ROUNTINE
C       KERNEL
C
C     NOTE:
C       1. CW_IRF is a matrix of MAXIRF by NLAYER, i.e., every channel
C          has a different center wavenumber in every layer
C          CW_IRF(*,NLAYER+1) is the cw value for the surface temperature
C          ISURF is used in order to calculate Planck only for the surface
C          layer and the layers above. (Saves calculation time.)
C       2. the Planck function value for the surface temperature SURF_T
C          will be written to DPLNCK(*,NLAYER+1)
C
C     HISTORY:
C       written 12/31/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     CALLING VARIABLES
      INTEGER     TOTIRF,ID_IRF(MAXIRF),  ! # of channels, channel IDs
     &            ISURF                   ! surface layer index
      REAL*8      ATEMP(NLAYER),          ! temperature profile
     &            SURF_T,                 ! surface temperature
     &            CW_IRF(MAXIRF,NLAYER+1),! channel center wavenumbers
     &            DPLNCK(MAXIRF,NLAYER+1) ! Planck function values
C     LOCAL VARIABLES
      INTEGER     IIRF,IIRF_A,ILAY1,ILAY2
      REAL*8      PLANCK(MAXIRF,NLAYER+1),PLNCKD(MAXIRF,NLAYER+1),
     &            DTEMP(NLAYER),DSURFT,DELTAT

      PARAMETER(DELTAT=1.0) ! temperature disturbance (delta T)
      
C     BEGIN
      CALL CPLNCK(TOTIRF,ID_IRF,ATEMP,CW_IRF,ISURF,SURF_T,PLANCK)
      DO 1000 ILAY2=1,NLAYER+1
         IF (ILAY2 .EQ. (NLAYER+1)) THEN
            DSURFT=SURF_T+DELTAT
            CALL CPLNCK(TOTIRF,ID_IRF,ATEMP,CW_IRF,ISURF,DSURFT,PLNCKD)
         ELSE
            DO 2001 ILAY1=1,NLAYER
               DTEMP(ILAY1)=ATEMP(ILAY1)
 2001       CONTINUE
            DTEMP(ILAY2)=ATEMP(ILAY2)+DELTAT
            CALL CPLNCK(TOTIRF,ID_IRF,DTEMP,CW_IRF,ISURF,SURF_T,PLNCKD)
         ENDIF
         DO 3000 IIRF=1,TOTIRF
            IIRF_A=ID_IRF(IIRF)
            DPLNCK(IIRF_A,ILAY2)=
     &           (PLNCKD(IIRF_A,ILAY2)-PLANCK(IIRF_A,ILAY2))/DELTAT
 3000    CONTINUE
 1000 CONTINUE

      RETURN
      END
C     END OF SUBROUTINE CPLNCK_FD

C ====================================================================

      SUBROUTINE C_J_UP(TOTIRF,ID_IRF,TAU,DTAUX,PLANCK,DPLNCK,JACUPX,
     &                  DMODE,ISURF,SUREMI)
C --------------------------------------------------------------------
C
C     PURPOSE:
C       Calculate the Jacobians for the upwelling thermal radiation
C
C     INPUT:
C       TOTIRF      total number of channels
C       ID_IRF      channel IDs
C       TAU         layer transmittances
C       DTAUX       d(TAU_LS)/dX  (X is either T,W,O), TAU_LS is
C                                   layer-to-space transmittance
C       PLANCK      Planck function values
C       DPLNCK      d(Planck)/dT
C       DMODE       1: if calculate d/dT, 2: if calculate d/dW or d/dO
C                   (if d/dT then d(Planck)/dT-term must be added)
C       ISURF       surface layer index
C       SUREMI      surface emissivity
C
C     OUTPUT
C       JACUPX      Jacobians for upwelling thermal: d(radiation)/dX
C
C     CALLING ROUNTINE
C       KERNEL
C
C     NOTE:
C       1. CW_IRF is a matrix of MAXIRF by NLAYER, i.e., every channel
C          has a different center wavenumber in every layer
C          CW_IRF(*,NLAYER+1) is the cw value for the surface temperature
C          ISURF is used in order to calculate Planck only for the surface
C          layer and the layers above. (Saves calculation time.)
C       2. the Planck function value for the surface temperature SURF_T
C          will be written to DPLNCK(*,NLAYER+1)
C
C     HISTORY:
C       written 12/02/97 by Tobias Wehr
C
C --------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"

C     CALLING VARIABLES
      INTEGER    TOTIRF,                      ! TOTAL NUMBER OF IRFs
     &           ID_IRF(MAXIRF),              ! IDs OF IRFs 
     &           DMODE,                       ! MODE IF d/dT
     &           ISURF                        ! SURFACE LAYER INDEX
      REAL*8     TAU(MAXIRF,NLAYER),          ! LAYER TRANSMITTANCES
     &           DTAUX(MAXIRF,NLAYER,NLAYER), ! dTAU/dX
     &           PLANCK(MAXIRF,NLAYER+1),     ! PLANCK FCT.
     &           DPLNCK(MAXIRF,NLAYER+1),     ! T-DERIVATIVE OF PLANCK FCT.
      ! NOTE: for PLANCK and DPLNCK the element NLAYER+1 refers to the
      !       surface, which is NOT needed in this subroutine!
     &           JACUPX(MAXIRF,NLAYER),       ! JACOBIANS d(rad)/dX
     &           SUREMI
C     LOCAL VARIABLES
      INTEGER    IIRF,IIRF_A,ILAY1,ILAY2
      REAL*8     DRDX,TAU_LS,TAUOLS

C     BEGIN
      DO 1000 IIRF=1,TOTIRF

         IIRF_A=ID_IRF(IIRF)
         DO 2002 ILAY2=NLAYER,ISURF,-1
            ! start (for ILAY1) in the highest layer NLAYER
            DRDX=PLANCK(IIRF_A,NLAYER)*(-1.0)*DTAUX(IIRF_A,NLAYER,ILAY2)
            IF ((DMODE .EQ. 1) .AND. (ILAY2 .EQ. NLAYER)) THEN
               TAU_LS=TAU(IIRF_A,NLAYER) ! layer-to-space tau for ILAY1
               DRDX=DRDX+(1.0-TAU_LS)*DPLNCK(IIRF_A,ILAY2)
               !remember: DPLNCK are the diagonal elements of the matrix
            ENDIF
            ! loop through all other layers down to ILAY2
            DO 2001 ILAY1=(NLAYER-1),ISURF,-1 !ILAY2,-1
               DRDX=DRDX+
     &              PLANCK(IIRF_A,ILAY1)*
     &              ( DTAUX(IIRF_A,ILAY1+1,ILAY2)-
     &                DTAUX(IIRF_A,ILAY1,ILAY2) )
               IF ((DMODE .EQ. 1) .AND. (ILAY1 .EQ. ILAY2)) THEN
                  TAUOLS=TAU_LS                   ! layer-to-space tau of i+1
                  TAU_LS=TAU_LS*TAU(IIRF_A,ILAY1) ! layer-to-space tau of i
                  DRDX=DRDX+DPLNCK(IIRF_A,ILAY2)*(TAUOLS-TAU_LS)
              ENDIF
 2001       CONTINUE
            ! surface layer
            DRDX=DRDX + SUREMI*
     &           PLANCK(IIRF_A,NLAYER+1)*DTAUX(IIRF_A,ISURF,ILAY2)
            JACUPX(IIRF_A,ILAY2)=DRDX
 2002    CONTINUE
         IF (ISURF .GT. 1) THEN
            DO 2003 ILAY2=1,(ISURF-1)
               JACUPX(IIRF_A,ILAY2)=0.0
 2003       CONTINUE
         ENDIF

 1000 CONTINUE

      RETURN
      END
C     END OF SUBROUTINE C_J_UP

C ==================================================================
      SUBROUTINE WRITE_DPLANCK(TOTIRF,ID_IRF,DPLNCK)
      ! THIS SUBROUTINE IS FOR DEBUGGING ONLY
      IMPLICIT NONE
      INCLUDE "hffp_glob_dec.f"
C     VARIABLES
      INTEGER     TOTIRF,ID_IRF(MAXIRF)
      REAL*8      DPLNCK(MAXIRF,NLAYER+1)
C     LOCAL VARIABLES
      INTEGER     IIRF,IIRF_A,ILAY
C     BEGIN
      
      DO 100 IIRF=1,TOTIRF
         IIRF_A=ID_IRF(IIRF)
         DO 200 ILAY=1,NLAYER+1
            WRITE (*,300) IIRF,ILAY,DPLNCK(IIRF_A,ILAY)
 200     CONTINUE
 100  CONTINUE
 300  FORMAT('dp(',I2,',',I3,')=',E25.15,';')
      RETURN
      END
C     END OF SUBROUTINE WRITE_DPLANCK

C==================================================

      SUBROUTINE WRITE_DTAUX(TOTIRF,ID_IRF,
     &     DTAUT,DTAUW,DTAUO,DTAUTS,DTAUWS,DTAUOS,
     &     FTAUT,FTAUW,FTAUO,FTAUTS,FTAUWS,FTAUOS)

C     THIS SUBROUTINE IS FOR WRITING DTAUX & FTAUX TO OUTPUT FILE
C     AND ONLY USED FOR DEBUGGING

      IMPLICIT NONE
      include "hffp_glob_dec.f"
C     CALLING PARAMETERS
      INTEGER  TOTIRF,
     &         ID_IRF(MAXIRF)
      REAL*8   DTAUT(MAXIRF,NLAYER,NLAYER),  ! dtau/dT, thermal upw.
     &         DTAUW(MAXIRF,NLAYER,NLAYER),  ! dtau/dW, thermal upw.
     &         DTAUO(MAXIRF,NLAYER,NLAYER),  ! dtau/dO, thermal upw.
     &         DTAUTS(MAXIRF,NLAYER,NLAYER), ! dtau/dT, refl. solar
     &         DTAUWS(MAXIRF,NLAYER,NLAYER), ! dtau/dW, refl. solar
     &         DTAUOS(MAXIRF,NLAYER,NLAYER)  ! dtau/dO, refl. solar
      REAL*8   FTAUT(MAXIRF,NLAYER,NLAYER),  ! dtau/dT, thermal upw.
     &         FTAUW(MAXIRF,NLAYER,NLAYER),  ! dtau/dW, thermal upw.
     &         FTAUO(MAXIRF,NLAYER,NLAYER),  ! dtau/dO, thermal upw.
     &         FTAUTS(MAXIRF,NLAYER,NLAYER), ! dtau/dT, refl. solar
     &         FTAUWS(MAXIRF,NLAYER,NLAYER), ! dtau/dW, refl. solar
     &         FTAUOS(MAXIRF,NLAYER,NLAYER)  ! dtau/dO, refl. solar
C     LOCAL PARAMETERS
      INTEGER       FP,IIRF,ILAY1,ILAY2,IIRF_A
      CHARACTER*80  FILENM,CDMMY1,CDMMY2,PATH,CIIRF

C     BEGIN
      FP=10

      !PATH='/scratch/HIRS/OUTPUT/ '  ! on squirrel
      !PATH='/salsify/data/Wehr/hirs/OUTPUT/ ' ! on salsify
      PATH='OUTPUT/' 
      DO 5000 IIRF=1,TOTIRF
         IIRF_A=ID_IRF(IIRF)
         IF (IIRF_A .LT. 10) THEN
            WRITE (CDMMY1,10) IIRF_A
         ELSE
            WRITE (CDMMY1,11) IIRF_A
         ENDIF
 10      FORMAT(I1)
 11      FORMAT(I2)
         CALL CATSTR(CDMMY1,' ',CIIRF)
         
         CALL CATSTR(PATH,'dtaudt_ ',CDMMY1)
         CALL CATSTR(CDMMY1,CIIRF,CDMMY2)
         CALL CATSTR(CDMMY2,'.m ',FILENM)
         WRITE (*,*) 'CREATING FILE ',FILENM
         OPEN(UNIT=FP,FILE=FILENM,FORM='FORMATTED')
         WRITE (FP,*) 'ajac=zeros(100,100);'
         WRITE (FP,*) 'fjac=zeros(100,100);'
         ! first 2 indecees is layer, 3rd is IIRF, 4th is anal. or fin.diff.
         DO 1002 ILAY1=1,NLAYER
            DO 1003 ILAY2=1,NLAYER
               WRITE (FP,201) ILAY1,ILAY2,DTAUT(IIRF_A,ILAY1,ILAY2)
               WRITE (FP,202) ILAY1,ILAY2,FTAUT(IIRF_A,ILAY1,ILAY2)
 1003       CONTINUE
 1002    CONTINUE
         CLOSE(FP)
         CALL CATSTR(PATH,'dtaudtS_ ',CDMMY1)
         CALL CATSTR(CDMMY1,CIIRF,CDMMY2)
         CALL CATSTR(CDMMY2,'.m ',FILENM)
         WRITE (*,*) 'CREATING FILE ',FILENM
         OPEN(UNIT=FP,FILE=FILENM,FORM='FORMATTED')
         WRITE (FP,*) 'ajac=zeros(100,100);'
         WRITE (FP,*) 'fjac=zeros(100,100);'
         ! first 2 indecees is layer, 3rd is IIRF, 4th is anal. or fin.diff.
         DO 1012 ILAY1=1,NLAYER
            DO 1013 ILAY2=1,NLAYER
               WRITE (FP,201) ILAY1,ILAY2,DTAUTS(IIRF_A,ILAY1,ILAY2)
               WRITE (FP,202) ILAY1,ILAY2,FTAUTS(IIRF_A,ILAY1,ILAY2)
 1013       CONTINUE
 1012    CONTINUE
         CLOSE(FP)

         CALL CATSTR(PATH,'dtaudw_ ',CDMMY1)
         CALL CATSTR(CDMMY1,CIIRF,CDMMY2)
         CALL CATSTR(CDMMY2,'.m ',FILENM)
         WRITE (*,*) 'CREATING FILE ',FILENM
         OPEN(UNIT=FP,FILE=FILENM,FORM='FORMATTED')
         WRITE (FP,*) 'ajac=zeros(100,100);'
         WRITE (FP,*) 'fjac=zeros(100,100);'
         ! first 2 indecees is layer, 3rd is IIRF, 4th is anal. or fin.diff.
         DO 2002 ILAY1=1,NLAYER
            DO 2003 ILAY2=1,NLAYER
               WRITE (FP,201) ILAY1,ILAY2,DTAUW(IIRF_A,ILAY1,ILAY2)
               WRITE (FP,202) ILAY1,ILAY2,FTAUW(IIRF_A,ILAY1,ILAY2)
 2003       CONTINUE
 2002    CONTINUE
         CLOSE(FP)
         CALL CATSTR(PATH,'dtaudwS_ ',CDMMY1)
         CALL CATSTR(CDMMY1,CIIRF,CDMMY2)
         CALL CATSTR(CDMMY2,'.m ',FILENM)
         WRITE (*,*) 'CREATING FILE ',FILENM
         OPEN(UNIT=FP,FILE=FILENM,FORM='FORMATTED')
         WRITE (FP,*) 'ajac=zeros(100,100);'
         WRITE (FP,*) 'fjac=zeros(100,100);'
         ! first 2 indecees is layer, 3rd is IIRF, 4th is anal. or fin.diff.
         DO 2012 ILAY1=1,NLAYER
            DO 2013 ILAY2=1,NLAYER
               WRITE (FP,201) ILAY1,ILAY2,DTAUWS(IIRF_A,ILAY1,ILAY2)
               WRITE (FP,202) ILAY1,ILAY2,FTAUWS(IIRF_A,ILAY1,ILAY2)
 2013       CONTINUE
 2012    CONTINUE
         CLOSE(FP)
         
         CALL CATSTR(PATH,'dtaudo_ ',CDMMY1)
         CALL CATSTR(CDMMY1,CIIRF,CDMMY2)
         CALL CATSTR(CDMMY2,'.m ',FILENM)
         WRITE (*,*) 'CREATING FILE ',FILENM
         OPEN(UNIT=FP,FILE=FILENM,FORM='FORMATTED')
         WRITE (FP,*) 'ajac=zeros(100,100);'
         WRITE (FP,*) 'fjac=zeros(100,100);'
         ! first 2 indecees is layer, 3rd is IIRF, 4th is anal. or fin.diff.
         DO 3002 ILAY1=1,NLAYER
            DO 3003 ILAY2=1,NLAYER
               WRITE (FP,201) ILAY1,ILAY2,DTAUO(IIRF_A,ILAY1,ILAY2)
               WRITE (FP,202) ILAY1,ILAY2,FTAUO(IIRF_A,ILAY1,ILAY2)
 3003       CONTINUE
 3002    CONTINUE
         CLOSE(FP)
         CALL CATSTR(PATH,'dtaudoS_ ',CDMMY1)
         CALL CATSTR(CDMMY1,CIIRF,CDMMY2)
         CALL CATSTR(CDMMY2,'.m ',FILENM)
         WRITE (*,*) 'CREATING FILE ',FILENM
         OPEN(UNIT=FP,FILE=FILENM,FORM='FORMATTED')
         WRITE (FP,*) 'ajac=zeros(100,100);'
         WRITE (FP,*) 'fjac=zeros(100,100);'
         ! first 2 indecees is layer, 3rd is IIRF, 4th is anal. or fin.diff.
         DO 3012 ILAY1=1,NLAYER
            DO 3013 ILAY2=1,NLAYER
               WRITE (FP,201) ILAY1,ILAY2,DTAUOS(IIRF_A,ILAY1,ILAY2)
               WRITE (FP,202) ILAY1,ILAY2,FTAUOS(IIRF_A,ILAY1,ILAY2)
 3013       CONTINUE
 3012    CONTINUE
         CLOSE(FP)

 5000 CONTINUE

 201  FORMAT('ajac(',I3,',',I3,')=',E25.15,';')
 202  FORMAT('fjac(',I3,',',I3,')=',E25.15,';')

      RETURN
      END
C     END OF SUBROUTINE WRITE_DTAUX
