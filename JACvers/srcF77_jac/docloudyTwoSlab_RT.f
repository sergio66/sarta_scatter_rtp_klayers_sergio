      SUBROUTINE docloudyTwoSlab_RT(I, IPROF, FREQ, LBOT, NWANTC, INDCHN, 
     $  TEMP,TSURF,TAU,TAUZ, TAUZSN, 
     $  BLMULT, EMIS, FCLEAR, COSDAZ, SECANG, SECSUN, DOSUN, SUNFAC, HSUN, RHOSUN, 
     $  RHOTHR, LABOVE, COEFF, LCTOP1, LCBOT1, LBLAC1, LCTOP2, LCBOT2, LBLAC2,       ! DEPEND ON CPRTOP1,CPRBOT1
     $  TCBOT1, TCTOP1, TCBOT2, TCTOP2, CFRAC1, CFRAC2, CLRB1, CLRT1, CLRB2, CLRT2,  ! DEPEND ON CPRTOP2,CPRBOT2
     $  CFRA12, CFRA1X, CFRA2X, 
     $  CEMIS1, CRHOS1, CRHOT1, CEMIS2, CRHOS2, CRHOT2, TEMPC1, TEMPC2,
     $  MASEC1, MASUN1, CFRCL1, G_ASY1, NEXTO1, NSCAO1, 
     $  MASEC2, MASUN2, CFRCL2, G_ASY2, NEXTO2, NSCAO2,
     $  POLYNOM_BACKSCAT1, POLYNOM_BACKSCAT2,ISCALING1,ISCALING2,
     $  QUICKINDNTE, NCHNTE, CLISTN, COEFN, SUNCOS, SCOS1, CO2TOP,
     $  RAD, DOJAC, TAU4, PLANCK_RAD4, PURE_RAD4, RTHERM4_SOLAR4, DBTDT,
     $               NLTEJACPRED5T,NLTEJACPRED7Q,
     $               RAACLDOD4,CLDEFFOD1,CLDEFFOD2)

      IMPLICIT NONE
      include "incFTC.f"

c input
       REAL CO2TOP             ! top layers CO2 mixing ratio
       INTEGER I,IPROF         ! channel list and profile number
       INTEGER INDCHN(MXCHAN)  ! array indices for all channels
       INTEGER NWANTC          ! number of wanted channels (-1=all)
       REAL   FREQ(MXCHAN)     ! chan center frequency
       INTEGER   LBOT          ! bottom layer index number
       REAL   TEMP(MAXLAY)     ! prof layer average temperature
       REAL  TSURF             ! surface temperature
       REAL BLMULT                ! bottom layer fractional multiplier
       REAL    TAU(MAXLAY,MXCHAN) ! chan layer effective optical depth
       REAL   TAUZ(MAXLAY,MXCHAN) ! chan surface-to-space trans
       REAL TAUZSN(MAXLAY,MXCHAN) ! sun space-to-surface-to-space OD
       REAL   EMIS(MXCHAN) ! chan surface emissivity

       REAL SUNCOS         ! cosine of sun zenith angle
       REAL SCOS1          ! cosine of sun zenith angle at layer1
       REAL FCLEAR            ! clear (no cloud) fraction of FOV
       LOGICAL DOSUN       ! do sun calc?
       REAL SECANG(MAXLAY) ! secant of effective     local path angle
       REAL SECSUN(MAXLAY) ! secant of effective sun local path angle
       REAL SUNFAC         ! sun solid angles times cosine at surface
       REAL   HSUN(MXCHAN) ! sun radiance (direct from sun)
       REAL RHOSUN(MXCHAN) ! chan reflectivity for sun
       REAL RHOTHR(MXCHAN) ! chan reflectivity for downwelling thermal
       INTEGER LABOVE(MXCHAN) ! chan downwelling thermal layer above
       REAL  COEFF(NFCOEF,MXCHAN)        ! coefs for chan "F" factor
       REAL COSDAZ         ! cosine(solazi - satazi) {COS Delta AZimuth}

       REAL CFRA12            ! cloud1+2(both) fraction of FOV
       REAL CFRA1X            ! cloud1(exclusively) fraction of FOV
       REAL CFRA2X            ! cloud2(exclusively) fraction of FOV

       INTEGER LCBOT1         ! layer containing cloud bottom
       INTEGER LCTOP1         ! layer containing cloud top
       INTEGER LCBOT2         ! layer containing cloud bottom
       INTEGER LCTOP2         ! layer containing cloud top
       LOGICAL LBLAC1  ! black cloud1? {Mie cloud if false}
       LOGICAL LBLAC2  ! black cloud2? {Mie cloud if false}
       REAL TCBOT1            ! temperature at cloud bottom
       REAL TCTOP1            ! temperature at cloud top, DEPENDS ON CPRTOP but only for black cloud
       REAL TCBOT2            ! temperature at cloud bottom
       REAL TCTOP2            ! temperature at cloud top, DEPENDS ON CPRTOP but only for black cloud
       REAL CFRAC1            ! cloud1(total) fraction of FOV
       REAL CFRAC2            ! cloud2(total) fraction of FOV
       REAL  CLRB1            ! frac of layer at bottom of cloud clear
       REAL  CLRT1            ! frac of layer at top of cloud clear, DEPENDS ON CPRTOP but only for black cloud
       REAL  CLRB2            ! frac of layer at bottom of cloud clear
       REAL  CLRT2            ! frac of layer at top of cloud clear, DEPENDS ON CPRTOP but only for black cloud
       REAL TEMPC1            ! cloud1 frac layer (above cloud) mean temp
       REAL TEMPC2            ! cloud2 frac layer (above cloud) mean temp

       REAL CEMIS1(MXCHAN) ! chan surface emissivity cloud1
       REAL CRHOS1(MXCHAN) ! chan solar reflectivity cloud1
       REAL CRHOT1(MXCHAN) ! chan thermal reflectivity cloud1
       REAL CEMIS2(MXCHAN) ! chan surface emissivity cloud2
       REAL CRHOS2(MXCHAN) ! chan solar reflectivity cloud2
       REAL CRHOT2(MXCHAN) ! chan thermal reflectivity cloud2

       REAL MASEC1            ! mean cloud view angle secant
       REAL MASUN1            ! mean cloud sun-only angle secant
       REAL CFRCL1(MAXLAY)    ! fraction of cloud in layer
       REAL G_ASY1(MXCHAN)    ! "g" asymmetry
       REAL NEXTO1(MXCHAN)    ! nadir extinction optical depth
       REAL NSCAO1(MXCHAN)    ! nadir scattering optical depth

       REAL MASEC2            ! mean cloud view angle secant
       REAL MASUN2            ! mean cloud sun-only angle secant
       REAL CFRCL2(MAXLAY)    ! fraction of cloud in layer
       REAL G_ASY2(MXCHAN)    ! "g" asymmetry
       REAL NEXTO2(MXCHAN)    ! nadir extinction optical depth
       REAL NSCAO2(MXCHAN)    ! nadir scattering optical depth

       REAL POLYNOM_BACKSCAT1(4), POLYNOM_BACKSCAT2(4) !!! backscatter coeffs
       INTEGER ISCALING1,ISCALING2                     !!! backscatter scaling

       INTEGER QUICKINDNTE(MXCHAN)       ! list of non-LTE channels
       INTEGER NCHNTE                    ! number of non-LTE channels
       INTEGER CLISTN(MXCNTE)            ! non-LTE channel list
       REAL  COEFN(NNCOEF,MXCNTE)        ! non-LTE coefficients

c output
       REAL  NLTEJACPRED5T(5,MXCHAN),NLTEJACPRED7Q(MXCHAN)
       REAL    RAD(MXCHAN) ! chan radiance
       LOGICAL DOJAC       ! are we planning on jacs???
       REAL TAU4(4,MAXLAY,MXCHAN) ! chan layer effective optical depth for CLR,CLD1,CLD2,CLD12       
       REAL PLANCK_RAD4(4,MAXLAY,MXCHAN) ! -radiance(L) + planck(TL)          for CLR,CLD1,CLD2,CLD12
       REAL PURE_RAD4(4,MAXLAY,MXCHAN)   ! +radiance(L)                       for CLR,CLD1,CLD2,CLD12
       REAL RTHERM4_SOLAR4(4,MAXLAY,MXCHAN)    ! downwell solar/ thermal background term at surface (about same for all 4 calcs but whatever)
       REAL DBTDT(MAXLAY,MXCHAN)           ! dBT(T,L)/dT
       REAL RAACLDOD4(4,MAXLAY,MXCHAN)     ! cloud ods for the 4 streams
       REAL CLDEFFOD1, CLDEFFOD2  ! calrad1 computes K1=NEXTO1(I) - NSCAO1(I)*(1.0+G_ASY1(I))/2.0 = total OD overr CLDL layers

c local
       REAL RTHERM         ! downwell thermal background term at surface
       REAL C1C2V4         ! rad constant c1c2 times freq^4
       REAL C1V3           ! rad constant c1 times freq^3
       REAL C2V            ! rad constant c2 times freq
       INTEGER L, III,IJUNK
       REAL RADNTE                       ! temporary NLTE rad                        

       REAL RPLNCK(MAXLAY) ! layer Planck
       REAL RSURFE         ! surface emission
       REAL RSURFC         ! black cloud surface emission
       REAL  TRANL(MAXLAY) ! clear air layer transmittance
       REAL  TRANZ(MXCHAN) ! clear air layer-to-space transmittance
       REAL  TRANS(MXCHAN) ! clear air total reflected solar trans

       REAL   RAD0         ! radiance no clouds
       REAL  RADC1         ! radiance cloud1
       REAL  RADC2         ! radiance cloud2
       REAL RADC12         ! radiance cloud1+cloud2

       REAL RADLAY(MAXLAY)        ! chan layer radiance                for CLDONLY
       REAL CLDTAU(MAXLAY)        ! chan layer effective optical depth for CLDONLY
       REAL VSTORE(6)             ! temporary storage for various variables
       REAL QIKEXP, RJUNK1, RJUNK2, RAD2BT

c************************************************************************
C     Radiation constants for current channel
      C1V3=C1*(FREQ(I)**3)
      C2V=C2*FREQ(I)

      IF (DOJAC) THEN 
        C1C2V4=C1*C2*(FREQ(I)**4)
      END IF

c     IF (NWANTC .GT. 0) THEN
c         write(*,'(A,I5,7(F12.4))') 'start rads/tau',I,FREQ(I),TSURF,TEMP(LBOT),
c    $             TAU(LBOT,I),TAUZ(LBOT-1,I),TAUZSN(LBOT,I),TAUZSN(LBOT-1,I)
c     END IF

C     Calculate Planck & clear airs trans for full layers
      DO L=1,LBOT-1
         RPLNCK(L)=C1V3/( EXP( C2V/TEMP(L) ) - 1.0 )
         TRANL(L)=QIKEXP( -TAU(L,I) )
      ENDDO

C     Note: TEMP(LBOT) already adjusted for bottom fractional layer
      RPLNCK(LBOT)=C1V3/( EXP( C2V/TEMP(LBOT) ) - 1.0 )

      IF (DOJAC) THEN
C       calculate planck derivatives        
        DO L=1,LBOT
           DBTDT(L,I)=C1C2V4 * EXP( C2V/TEMP(L) ) / ((TEMP(L)*( EXP( C2V/TEMP(L) ) - 1.0 ))**2)
        END DO
        L = LBOT+1
        DBTDT(L,I)=C1C2V4 * EXP( C2V/TSURF ) / ((TSURF*( EXP( C2V/TSURF ) - 1.0 ))**2)
      END IF

C     Calculate clear airs trans for bottom fractional layer
      RJUNK1=-TAU(LBOT,I)*BLMULT
      TRANL(LBOT)=QIKEXP( RJUNK1 )
      TRANL(LBOT)=QIKEXP( RJUNK1 )
      TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LBOT-1,I) )
      TRANS(I)=QIKEXP( BLMULT*(TAUZSN(LBOT-1,I)-TAUZSN(LBOT,I)) -
     $         TAUZSN(LBOT-1,I) )

C     Planck for surface
      RSURFE=EMIS(I)*C1V3/( EXP( C2V/TSURF ) - 1.0 )

C     Calculate clear radiance
      IF (FCLEAR .GT. 0.0) THEN
         CALL CALRAD0( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $       TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $       RHOTHR, LABOVE, COEFF, RAD0, DOJAC, CLDTAU, RADLAY, RTHERM )
      ELSE
        RAD0=0.0
      ENDIF
      IF (DOJAC) THEN
        TAU4(1,:,I) = TAU(:,I)   !!! CLDTAU is a dummy
        RTHERM4_SOLAR4(1,:,I) = RTHERM*SECANG     !!! I think this is a little odd to include since writeoutjac does 
                                                  !!! - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * dTAU_DG1(1:NLAY,1:NCHAN)
                                                  !!! and that dTAU_DG1 etc derivative already should have factor of 1/mu
        PURE_RAD4(1,:,I)   = (+RADLAY)
        PLANCK_RAD4(1,:,I) = (-RADLAY + RPLNCK)
        !PLANCK_RAD4(1,:,I) = (-RADLAY + RPLNCK)*SECANG
        RAACLDOD4(1,:,I) = TAU(:,I)  !!! so this is clear OD, which we need to keep substracting from (clr_cld) OD ... need to reset this to zero at the end
      END IF

C     Store original values
      VSTORE(1)=TRANL(LCTOP2)
      VSTORE(2)=TRANZ(I)
      VSTORE(3)=TRANS(I)
      VSTORE(4)=RHOTHR(I)
      VSTORE(5)=RHOSUN(I)
      VSTORE(6)=RPLNCK(LCTOP2)
C     Updates for new surface if bottom cloud2 is black
      IF (CFRAC2 .GT. 0.0 .AND. LBLAC2) THEN
         RJUNK1=-TAU(LCTOP2,I)*CLRT2
         TRANL(LCTOP2)=QIKEXP( RJUNK1 )
         TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LCTOP2-1,I) )
         TRANS(I)=QIKEXP( CLRT2*(TAUZSN(LCTOP2-1,I)-
     $          TAUZSN(LCTOP2,I)) - TAUZSN(LCTOP2-1,I) )
         RSURFC=CEMIS2(I)*C1V3/( EXP( C2V/TCTOP2 ) - 1.0 )
         RHOTHR(I)=CRHOT2(I)
         RHOSUN(I)=CRHOS2(I)
         RPLNCK(LCTOP2)=C1V3/( EXP( C2V/TEMPC2 ) - 1.0 )
c             RSURFC=C1V3/( EXP( C2V/TEMPC2 ) - 1.0 )
          ENDIF

C    Calculate bottom cloud2 radiance
      IF (CFRA2X .GT. 0.0) THEN
         IF (LBLAC2) THEN
            CALL CALRAD0( DOSUN, I, LCTOP2, RPLNCK, RSURFC, SECANG,
     $          TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, RADC2, DOJAC, CLDTAU, RADLAY, RTHERM )
         ELSE
            CALL CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $          TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, CFRCL2, MASEC2, MASUN2, COSDAZ,
     $          NEXTO2, NSCAO2, G_ASY2, LCTOP2, LCBOT2, 
     $          POLYNOM_BACKSCAT2, ISCALING2, 
     $          RADC2, DOJAC, CLDTAU, RADLAY, RTHERM, CLDEFFOD2 )
         ENDIF
      ELSE
         RADC2=0.0
      ENDIF
      RAACLDOD4(3,:,I) = 0
      IF (DOJAC) THEN
        TAU4(3,:,I) = CLDTAU
        RTHERM4_SOLAR4(3,:,I) = RTHERM*SECANG     !!! I think this is a little odd to include since writeoutjac does 
                                                  !!! - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * dTAU_DG1(1:NLAY,1:NCHAN)
                                                  !!! and that dTAU_DG1 etc derivative already should have factor of 1/mu
        PURE_RAD4(3,:,I)   = (+RADLAY)
        PLANCK_RAD4(3,:,I) = (-RADLAY + RPLNCK)
        !PLANCK_RAD4(3,:,I) = (-RADLAY + RPLNCK)*SECANG
        IF (CFRA2X .GT. 0.0) RAACLDOD4(3,:,I) = CLDTAU - RAACLDOD4(1,:,I)   !!! so this is cld1 OD
      END IF

C      Calculate combined cloud1+cloud2 radiance
      IF (CFRA12 .GT. 0.0) THEN
         IF (LBLAC2) THEN
            CALL CALRAD1( DOSUN, I, LCTOP2, RPLNCK, RSURFC, SECANG,
     $          TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, COSDAZ,
     $          NEXTO1, NSCAO1, G_ASY1, LCTOP1, LCBOT1, 
     $          POLYNOM_BACKSCAT1, ISCALING1, 
     $          RADC12, DOJAC, CLDTAU, RADLAY, RTHERM, CLDEFFOD1 )
         ELSE
            CALL CALRAD2( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $          TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, COSDAZ,
     $          CFRCL1, MASEC1, MASUN1, NEXTO1, NSCAO1, G_ASY1, LCTOP1, LCBOT1, 
     $          CFRCL2, MASEC2, MASUN2, NEXTO2, NSCAO2, G_ASY2, LCTOP2, LCBOT2, 
     $          POLYNOM_BACKSCAT1, POLYNOM_BACKSCAT2, ISCALING1, ISCALING2,
     $          RADC12, DOJAC, CLDTAU, RADLAY, RTHERM )
         ENDIF
      ELSE
         RADC12=0.0
      ENDIF
      RAACLDOD4(4,:,I) = 0
      IF (DOJAC) THEN
        TAU4(4,:,I) = CLDTAU
        RTHERM4_SOLAR4(4,:,I) = RTHERM*SECANG     !!! I think this is a little odd to include since writeoutjac does 
                                                  !!! - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * dTAU_DG1(1:NLAY,1:NCHAN)
                                                  !!! and that dTAU_DG1 etc derivative already should have factor of 1/mu
        PURE_RAD4(4,:,I)   = (+RADLAY)
        PLANCK_RAD4(4,:,I) = (-RADLAY + RPLNCK)
        !PLANCK_RAD4(4,:,I) = (-RADLAY + RPLNCK)*SECANG
        IF (CFRA12 .GT. 0) RAACLDOD4(4,:,I) = CLDTAU - RAACLDOD4(1,:,I)   !!! so this is cld1+cld2 OD
      END IF

C     Restore original values
      TRANL(LCTOP2)=VSTORE(1)
      TRANZ(I)=VSTORE(2)
      TRANS(I)=VSTORE(3)
      RHOTHR(I)=VSTORE(4)
      RHOSUN(I)=VSTORE(5)
      RPLNCK(LCTOP2)=VSTORE(6)
C     Updates for new surface if top cloud1 is black
      IF (CFRAC1 .GT. 0.0 .AND. LBLAC1) THEN
        RJUNK1=-TAU(LCTOP1,I)*CLRT1
        TRANL(LCTOP1)=QIKEXP( RJUNK1 )
        TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LCTOP1-1,I) )
        TRANS(I)=QIKEXP( CLRT1*(TAUZSN(LCTOP1-1,I)-
     $          TAUZSN(LCTOP1,I)) - TAUZSN(LCTOP1-1,I) )
        RSURFC=CEMIS1(I)*C1V3/( EXP( C2V/TCTOP1 ) - 1.0 )
        RHOTHR(I)=CRHOT1(I)
        RHOSUN(I)=CRHOS1(I)
        RPLNCK(LCTOP1)=C1V3/( EXP( C2V/TEMPC1 ) - 1.0 )
c       RSURFC=C1V3/( EXP( C2V/TEMPC1 ) - 1.0 )
      ENDIF

C     Calculate top cloud1 radiance
      IF (CFRA1X .GT. 0.0) THEN
         IF (LBLAC1) THEN
            CALL CALRAD0( DOSUN, I, LCTOP1, RPLNCK, RSURFC, SECANG,
     $          TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, RADC1, DOJAC, CLDTAU, RADLAY, RTHERM )
         ELSE
            CALL CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $          TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $          RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, COSDAZ,
     $          NEXTO1, NSCAO1, G_ASY1, LCTOP1, LCBOT1, 
     $          POLYNOM_BACKSCAT1, ISCALING1, 
     $          RADC1, DOJAC, CLDTAU, RADLAY, RTHERM, CLDEFFOD1 )
         ENDIF
      ELSE
         RADC1=0.0
      ENDIF
      RAACLDOD4(2,:,I) = 0
      IF (DOJAC) THEN
        TAU4(2,:,I) = CLDTAU
        RTHERM4_SOLAR4(2,:,I) = RTHERM*SECANG     !!! I think this is a little odd to include since writeoutjac does 
                                                  !!! - RTHERM4_SOLAR4(1,1:NLAY,1:NCHAN) * dTAU_DG1(1:NLAY,1:NCHAN)
                                                  !!! and that dTAU_DG1 etc derivative already should have factor of 1/mu
        PURE_RAD4(2,:,I)   = (+RADLAY)
        PLANCK_RAD4(2,:,I) = (-RADLAY + RPLNCK)
        !PLANCK_RAD4(2,:,I) = (-RADLAY + RPLNCK)*SECANG
        IF (CFRA1X .GT. 0) RAACLDOD4(2,:,I) = CLDTAU - RAACLDOD4(1,:,I)   !!! so this is cld2 OD
      END IF
      RAACLDOD4(1,:,I) = 0                           !!! so now say this is zero as there is no cloud in the first stream

c see test_cld_jacs_rads.m       
c      DO IJUNK = 1,LBOT
c        write(*,'(A,I4,13(ES16.7))') 'check taucldOD and radpure',IJUNK,TAU4(1,IJUNK,I),
c     $                          RAACLDOD4(1,IJUNK,I),RAACLDOD4(2,IJUNK,I),RAACLDOD4(3,IJUNK,I),RAACLDOD4(4,IJUNK,I),
c     $                          PURE_RAD4(1,IJUNK,I),PURE_RAD4(2,IJUNK,I),PURE_RAD4(3,IJUNK,I),PURE_RAD4(4,IJUNK,I),
c     $                          PLANCK_RAD4(1,IJUNK,I),PLANCK_RAD4(2,IJUNK,I),PLANCK_RAD4(3,IJUNK,I),PLANCK_RAD4(4,IJUNK,I)
c      END DO

c      Total the clear & various cloudy radiances
      RAD(I)=RAD0*FCLEAR + RADC1*CFRA1X + RADC2*CFRA2X + RADC12*CFRA12

c rdinfo.f : these are 8 nice chans to look at
c       sarta  fin=input.rtp  fout=output.rtp  listp=1,2,3 listc=445,449,1092,1291,1614,2070,2333,2353 listj=-1
c          would give you radiances at                           790,791,1042,1231,1419,2350,2616,2637 cm-1
c                                                                 CO2Q    O3   WIN  WV  NLTE  WIN  HDO
      IF (  ((NWANTC .GT. 0) .AND. DEBUG) .OR. ((NWANTC .GT. 0) .AND. (NWANTC .LE. 8)) ) THEN
        write(*,'(A,I5,3(F10.4),A,4(F10.4),A,4(ES10.3),A,ES10.3,F12.4)') ' [Index,Freq,Emis,Tsuf ]',I,FREQ(I),EMIS(I),TSURF,
     $           ' [ Cldfrac fc,c1x,c2x,c12 ]',FCLEAR,CFRA1X,CFRA2X,CFRA12,
     $           ' [ CldRad  rc,r1x,r2x,r12 ]',RAD0,RADC1,RADC2,RADC12,' [ final rad/BT ]',RAD(I),rad2bt(FREQ(I),RAD(I))
      END IF

ccc this block for testing
C       IF (I .EQ. 1291) THEN
c         print *,'chan1291 : iPROF,rad0,radc1,radc2,radc12,FINAL=',
c     $      IPROF,RAD0,RADC1,RADC2,RADC12,RAD(I)
C         print *,'chan1291 : IPROF,rad0,FCLEAR,CFRA1X,CFRA2X,CFRA12=',
C     $      IPROF,RAD0,FCLEAR,CFRA1X,CFRA2X,CFRA12
c         PRINT *,'CLOUD1 emis,temp = ',CEMIS1(I),TCTOP1
c         PRINT *,'CLOUD2 emis,temp = ',CEMIS2(I),TCTOP2
c       endif

c Tang
c       IF (I .EQ. 1291) THEN
c         write(*,'(I8,I8,13(F12.4))') 
c     c    I,IPROF,FCLEAR,CFRA1X,CFRA2X,CFRA12,TCTOP1,TCTOP2,CLDEFFOD1,CLDEFFOD2,RAD0,RADC1,RADC2,RADC12,RAD(I)
c       endif

ccc

C         -----------------
C         Calculate non-LTE
C         -----------------
C         comment: the nonLTE calculation does not consider cloud effects,
C         but clouds are generally below the altitude where nonLTE occurs.
      IF (DOSUN) THEN
CC FORGET THIS  III = INTERSECT(I,CLISTN(1:NCHNTE),NCHNTE)
        III = QUICKINDNTE(I)
        IF (III .GT. 0) THEN
          RADNTE = RAD(I)
          CALL YCALNTE ( INDCHN, TEMP, SUNCOS, SCOS1, SECANG(1),
     $                        NCHNTE, CLISTN, COEFN, CO2TOP, RADNTE, III,
     $                        DOJAC,NLTEJACPRED5T,NLTEJACPRED7Q,I)
c          IF (NWANTC .GT. 0) print *,'NLTE',I,III,RAD(I),RADNTE
          RAD(I) = RADNTE           
        END IF
      ENDIF

      !!!! so we already do 
      !!!!      PLANCK_RAD4(2,:,I) = (-RADLAY + RPLNCK)*SECANG  etc == (-RADLAY + RPLNCK) /mu
      !!!!      PURERAD4(2,:,I)    = (+RADLAY)         *SECANG  etc == (RADLAY) /mu

! orig code Apr 2023
!      IF (DOJAC) THEN 
!        PLANCK_RAD4(1,1:LBOT,I) = PLANCK_RAD4(1,1:LBOT,I)/SECANG(1:LBOT)
!        PLANCK_RAD4(2,1:LBOT,I) = PLANCK_RAD4(2,1:LBOT,I)/SECANG(1:LBOT)
!        PLANCK_RAD4(3,1:LBOT,I) = PLANCK_RAD4(3,1:LBOT,I)/SECANG(1:LBOT)
!        PLANCK_RAD4(4,1:LBOT,I) = PLANCK_RAD4(4,1:LBOT,I)/SECANG(1:LBOT)
!        PURE_RAD4(1,1:LBOT,I) = PURE_RAD4(1,1:LBOT,I)/SECANG(1:LBOT)
!        PURE_RAD4(2,1:LBOT,I) = PURE_RAD4(2,1:LBOT,I)/SECANG(1:LBOT)
!        PURE_RAD4(3,1:LBOT,I) = PURE_RAD4(3,1:LBOT,I)/SECANG(1:LBOT)
!        PURE_RAD4(4,1:LBOT,I) = PURE_RAD4(4,1:LBOT,I)/SECANG(1:LBOT)
!      END IF

      RETURN
      END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL FUNCTION rad2bt(f,r)

      IMPLICIT NONE
      include "incFTC.f"

      REAL f,r

      REAL x1,x2

      x1 = c2*f
      x2 = log(1 + c1*f**3/r)
      rad2bt = x1/x2
      RETURN
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
