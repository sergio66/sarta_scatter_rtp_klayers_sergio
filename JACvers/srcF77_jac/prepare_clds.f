      SUBROUTINE PREPARE_CLDS(
     $    LBLAC1, CTYPE1, CFRAC1, CPSIZ1, CPRTO1, CPRBO1, CNGWA1,  ! from GETCLD
     $    XCEMI1, XCRHO1, CSTMP1,                                  ! from GETCLD
     $    LBLAC2, CTYPE2, CFRAC2, CPSIZ2, CPRTO2, CPRBO2, CNGWA2,  ! from GETCLD
     $    XCEMI2, XCRHO2, CSTMP2, CFRA12, FCLEAR, CFRA1X, CFRA2X,  ! from GETCLD
     $    NCHAN, IPROF, LBOT, PSURF, PLEV, PLAY, TEMP, SECANG, SECSUN,
     $    MIETYP, MIENPS, MIEPS, MIEABS, MIEEXT, MIEASY,
     $    LCBOT1, LCTOP1, CLRB1, CLRT1, TCBOT1, TCTOP1, MASEC1, MASUN1,
     $    CFRCL1, G_ASY1, NEXTO1, NSCAO1, TEMPC1, 
     $    LCBOT2, LCTOP2, CLRB2, CLRT2, TCBOT2, TCTOP2, MASEC2, MASUN2,
     $    CFRCL2, G_ASY2, NEXTO2, NSCAO2, TEMPC2,
     $    DOJAC, 
     $    JACA_G_ASY1, JACA_NEXTO1, JACA_NSCAO1, JACA_FINAL_1, 
     $    JACS_G_ASY1, JACS_NEXTO1, JACS_NSCAO1, JACS_FINAL_1,
     $    JACTOP_CFRCL1,JACBOT_CFRCL1,JACTOP_CFRCL1_v,JACBOT_CFRCL1_v,
     $    JACA_G_ASY2, JACA_NEXTO2, JACA_NSCAO2, JACA_FINAL_2, 
     $    JACS_G_ASY2, JACS_NEXTO2, JACS_NSCAO2, JACS_FINAL_2,
     $    JACTOP_CFRCL2,JACBOT_CFRCL2,JACTOP_CFRCL2_v,JACBOT_CFRCL2_v)

      IMPLICIT NONE

      include "incFTC.f"

c input from GETCLD
       LOGICAL LBLAC1  ! black cloud1? {Mie cloud if false}
       LOGICAL LBLAC2  ! black cloud2? {Mie cloud if false}
       INTEGER CTYPE1         ! cloud1 type code number
       INTEGER CTYPE2         ! cloud2 type code number
       REAL CFRAC1            ! cloud1(total) fraction of FOV
       REAL CFRAC2            ! cloud2(total) fraction of FOV
       REAL CFRA12            ! cloud1+2(both) fraction of FOV
       REAL CNGWA1            ! cloud1 non-gases water
       REAL CNGWA2            ! cloud1 non-gases water
       REAL CPRBO1            ! cloud1 bottom pressure
       REAL CPRBO2            ! cloud2 bottom pressure
       REAL CPRTO1            ! cloud1 top pressure
       REAL CPRTO2            ! cloud2 top pressure
       REAL CPSIZ1            ! cloud1 particle size
       REAL CPSIZ2            ! cloud2 particle size
       REAL XCEMI1(MXEMIS)    ! cloud1 emissivity
       REAL XCEMI2(MXEMIS)    ! cloud2 emissivity
       REAL XCRHO1(MXEMIS)    ! cloud1 reflectivity
       REAL XCRHO2(MXEMIS)    ! cloud2 reflectivity
       REAL CSTMP1            ! cloud1 top/surf temperature
       REAL CSTMP2            ! cloud2 top/surf temperature
       REAL CFRA1X            ! cloud1(exclusively) fraction of FOV
       REAL CFRA2X            ! cloud2(exclusively) fraction of FOV
       REAL FCLEAR            ! clear (no cloud) fraction of FOV

       REAL SECANG(MAXLAY)        ! local path angle secant
       REAL SECSUN(MAXLAY) ! secant of effective sun local path angle
       INTEGER  NCHAN          ! # of selected channels
       INTEGER  IPROF      ! profile loop counter
       INTEGER   LBOT             ! bottom layer index number
       REAL  PSURF                ! surface pressure
       REAL PLEV(MAXLAY+1)
       REAL PLAY(MAXLAY)   ! layer mean pressure
       REAL   TEMP(MAXLAY) ! prof layer average temperature
       
       INTEGER MIETYP(NMIETY)      ! mie type
       REAL  MIEPS(MXMIEA,NMIETY)        ! Mie particle size for table
       INTEGER MIENPS(NMIETY)            ! number of particle sizes
       REAL MIEABS(MXCHAN,MXMIEA,NMIETY) ! Mie absorption table
       REAL MIEEXT(MXCHAN,MXMIEA,NMIETY) ! Mie extinction table
       REAL MIEASY(MXCHAN,MXMIEA,NMIETY) ! Mie asymmetry table
       LOGICAL DOJAC

c output
C      for CCPREP cloud1
       INTEGER LCBOT1         ! layer containing cloud bottom
       INTEGER LCTOP1         ! layer containing cloud top
       REAL  CLRB1            ! frac of layer at bottom of cloud clear
       REAL  CLRT1            ! frac of layer at top of cloud clear
       REAL TCBOT1            ! temperature at cloud bottom
       REAL TCTOP1            ! temperature at cloud top
       REAL MASEC1            ! mean cloud view angle secant
       REAL MASUN1            ! mean cloud sun-only angle secant
       REAL CFRCL1(MAXLAY)    ! fraction of cloud in layer
       REAL G_ASY1(MXCHAN)    ! "g" asymmetry
       REAL NEXTO1(MXCHAN)    ! nadir extinction optical depth
       REAL NSCAO1(MXCHAN)    ! nadir scattering optical depth
C
C      for CCPREP cloud2
       INTEGER LCBOT2         ! layer containing cloud bottom
       INTEGER LCTOP2         ! layer containing cloud top
       REAL  CLRB2            ! frac of layer at bottom of cloud clear
       REAL  CLRT2            ! frac of layer at top of cloud clear
       REAL TCBOT2            ! temperature at cloud bottom
       REAL TCTOP2            ! temperature at cloud top
       REAL MASEC2            ! mean cloud view angle secant
       REAL MASUN2            ! mean cloud sun-only angle secant
       REAL CFRCL2(MAXLAY)    ! fraction of cloud in layer
       REAL G_ASY2(MXCHAN)    ! "g" asymmetry
       REAL NEXTO2(MXCHAN)    ! nadir extinction optical depth
       REAL NSCAO2(MXCHAN)    ! nadir scattering optical depth
       REAL TEMPC1            ! cloud1 frac layer (above cloud) mean temp
       REAL TEMPC2            ! cloud2 frac layer (above cloud) mean temp

       REAL JACA_G_ASY1(MXCHAN),JACA_NEXTO1(MXCHAN),JACA_NSCAO1(MXCHAN),JACA_FINAL_1(MXCHAN)  !! amount jacs
       REAL JACS_G_ASY1(MXCHAN),JACS_NEXTO1(MXCHAN),JACS_NSCAO1(MXCHAN),JACS_FINAL_1(MXCHAN)  !! sze jacs
       REAL JACA_G_ASY2(MXCHAN),JACA_NEXTO2(MXCHAN),JACA_NSCAO2(MXCHAN),JACA_FINAL_2(MXCHAN)  !! amount jacs
       REAL JACS_G_ASY2(MXCHAN),JACS_NEXTO2(MXCHAN),JACS_NSCAO2(MXCHAN),JACS_FINAL_2(MXCHAN)  !! sze jacs
       REAL JACTOP_CFRCL1(MAXLAY),JACBOT_CFRCL1(MAXLAY)    ! derivatives of fraction of cloud in layer
       REAL JACTOP_CFRCL2(MAXLAY),JACBOT_CFRCL2(MAXLAY)    ! derivatives of fraction of cloud in layer
       REAL JACTOP_CFRCL1_v(MAXLAY,MXCHAN),JACBOT_CFRCL1_v(MAXLAY,MXCHAN)    ! spectral derivatives of fraction of cloud in laye
       REAL JACTOP_CFRCL2_v(MAXLAY,MXCHAN),JACBOT_CFRCL2_v(MAXLAY,MXCHAN)    ! spectral derivatives of fraction of cloud in layer
c local
       INTEGER IERR1, IERR2
       INTEGER INDMI1  ! index in MIETYP for CTYPE1
       INTEGER INDMI2  ! index in MIETYP for CTYPE2

c************************************************************************

C        Check and prepare (top) cloud1
       IF (CFRAC1 .GT. 0.0) THEN
         IF (LBLAC1) THEN
           CALL BKPREP(IPROF, 1, CTYPE1, CFRAC1, CPRTO1,
     $          LBOT, PSURF, PLEV, PLAY, TEMP, LCTOP1, TCTOP1,
     $          TEMPC1, CLRT1)
           IF (CSTMP1 .GT. 0.0) TCTOP1=CSTMP1
         ELSE
C           Determine which lookup table to use
            CALL GETMIE(CTYPE1,MIETYP,INDMI1,IERR1)
C           Prepare selected lookup table for given cpsize
            CALL CCPREP( NCHAN, LBOT, INDMI1, MIENPS,
     $          CNGWA1, CPSIZ1, CPRTO1, CPRBO1, PLEV, TEMP, SECANG,
     $          SECSUN, MIEPS, MIEABS, MIEEXT, MIEASY, LCBOT1, LCTOP1,
     $          CLRB1, CLRT1, TCBOT1, TCTOP1, MASEC1, MASUN1,
     $          CFRCL1, G_ASY1, NEXTO1, NSCAO1, 
     $          DOJAC, JACA_G_ASY1, JACA_NEXTO1, JACA_NSCAO1, JACA_FINAL_1, 
     $                 JACS_G_ASY1, JACS_NEXTO1, JACS_NSCAO1, JACS_FINAL_1,
     $                 JACTOP_CFRCL1,JACBOT_CFRCL1,JACTOP_CFRCL1_v,JACBOT_CFRCL1_v)
          ENDIF
       ENDIF

C      Check and prepare (bottom) cloud2
       IF (CFRAC2 .GT. 0.0) THEN
          IF (LBLAC2) THEN
             CALL BKPREP(IPROF, 2, CTYPE2, CFRAC2, CPRTO2,
     $          LBOT, PSURF, PLEV, PLAY, TEMP, LCTOP2, TCTOP2,
     $          TEMPC2, CLRT2)
           IF (CSTMP2 .GT. 0.0) TCTOP2=CSTMP2
          ELSE
C          Determine which lookup table to use
           CALL GETMIE(CTYPE2,MIETYP,INDMI2,IERR2)
C          Prepare lookup data for cloud2
           CALL CCPREP( NCHAN, LBOT, INDMI2, MIENPS,
     $          CNGWA2, CPSIZ2, CPRTO2, CPRBO2, PLEV, TEMP, SECANG,
     $          SECSUN, MIEPS, MIEABS, MIEEXT, MIEASY, LCBOT2, LCTOP2,
     $          CLRB2, CLRT2, TCBOT2, TCTOP2, MASEC2, MASUN2,
     $          CFRCL2, G_ASY2, NEXTO2, NSCAO2, 
     $          DOJAC, JACA_G_ASY2, JACA_NEXTO2, JACA_NSCAO2, JACA_FINAL_2,
     $                 JACS_G_ASY2, JACS_NEXTO2, JACS_NSCAO2, JACS_FINAL_2,
     $                 JACTOP_CFRCL2,JACBOT_CFRCL2,JACTOP_CFRCL2_v,JACBOT_CFRCL2_v)
         ENDIF
       ELSE
C         Safe default for non-existant cloud2
          LCTOP2=1
       ENDIF

       RETURN 
       END
