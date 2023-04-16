C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore County (UMBC)
C
C    AIRS
C
C    SCATUMBC
C
!F77====================================================================


!ROUTINE NAME:
C    SCATUMBC


!ABSTRACT:
C    Subroutine to model upto 2 complex clouds using the PCLSAM method.


!CALL PROTOCOL
C    none (main program)


!INPUT PARAMETERS:
C    none


!OUTPUT PARAMETERS:
C    none


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S)
C    sarta_pclsam


!ROUTINES CALLED:
C    CALRAD : calc radiance


!FILES ACCESSED:
C    NONE

!COMMON BLOCKS
C    COMLEV : layer boundary pressure levels
c    py_fittedcoe : Ping Yang scattering coeffs


!DESCRIPTION:
C    Dec 2005 version of the SARTA (Stand-Alone Rapid
C    Transmittance Algorith with RTP I/O) by
C    L.L.Strow, S.Hannon, and H.Mottler
C
C    Computes radiances for the layers profiles contained in the
C    input RTP file.  This is the main program, and consists
C    primarily of calls to the external routines to do most of the
C    computing.


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    This program is only intended as a demo of the fast model.


!ROUTINE HISTORY:
C Date         Programmer    Comments
C ----------- -------------- -------------------------------------------
C 25 Jan 2009 Sergio Machado Created

!END====================================================================

       SUBROUTINE SCATUMBC (
     $              NCHAN,FREQ,TEMP,TAU,TAUZ,TAUZSN,EMIS,TSURF,
     $              SVA,SECANG,SUNFAC,HSUN,RHOSUN,RHOTHR,LABOVE,COEFF,
     $              DOSUN,LBOT,NLAY,FCLEAR,CFRA12,
     $              BLMULT, SUNCOS,COSDAZ,SOLARTOA,PLAY,
     $              CTYPE1,CNGWA1,CPSIZ1,CPRBO1,CPRTO1,
     $              CFRAC1,LBLAC1,CFRA1X,TCTOP1,TEMPC1,
     $              CFRCL1, MASEC1, MASUN1,
     $              NEXTO1, NSCAO1, G_ASY1, LCTOP1, LCBOT1, 
     $              CTYPE2,CNGWA2,CPSIZ2,CPRBO2,CPRTO2,
     $              CFRAC2,LBLAC2,CFRA2X,TCTOP2,TEMPC2,
     $              CFRCL2, MASEC2, MASUN2, 
     $              NEXTO2, NSCAO2, G_ASY2, LCTOP2, LCBOT2, 
     $              IPYCODE,RAD)

C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE

C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
       include 'incFTC.f'
       include 'rtpdefs.f'
       include 'pyang2cld.param'

C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
c      many

C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------

       REAL   TEMP(MAXLAY) ! prof layer average temperature
       INTEGER CTYPE1         ! cloud1 type code number
       INTEGER CTYPE2         ! cloud2 type code number
       INTEGER NLAY        ! number of layers in profile
       REAL CNGWA1            ! cloud1 non-gases water or <tau550nm>
       REAL CNGWA2            ! cloud1 non-gases water or <tau550nm>
       REAL CPSIZ1            ! cloud1 particle size
       REAL CPSIZ2            ! cloud2 particle size

C      FOR PING YANGs CODE
       INTEGER IPYCODE     !-1 for ONE dust cloud and/or ONE black cloud
                           !-1 for TWO black clouds
                           !-1 for TWO complex clouds (dust/ice/water)
                           !+1 for ONE ice/water cloud and/or ONE black cloud

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
C
       INTEGER  NCHAN         ! # of selected channels
       REAL   FREQ(MXCHAN)    ! chan center frequency
       REAL BLMULT                ! bottom layer fractional multiplier
       INTEGER LABOVE(MXCHAN) ! chan downwelling thermal layer above
       REAL  COEFF(NFCOEF,MXCHAN)        ! coefs for chan "F" factor

C      for surface
       INTEGER   LBOT             ! bottom layer index number
       REAL SECANG(MAXLAY)        ! local path angle secant
       REAL    SVA         ! satellite viewing angle (degrees)
C
C      for CALT
       REAL    TAU(MAXLAY,MXCHAN) ! chan layer effective optical depth
       REAL   TAUZ(MAXLAY,MXCHAN) ! chan surface-to-space trans
C
C      for CALRAD
       REAL SUNFAC         ! sun solid angles times cosine at surface
       REAL RPLNCK(MAXLAY) ! layer Planck
       REAL RSURFE         ! surface emission
       REAL RSURFC         ! black cloud surface emission
       REAL  TRANL(MAXLAY) ! clear air layer transmittance
       REAL  TRANZ(MXCHAN) ! clear air layer-to-space transmittance
       REAL  TRANS(MXCHAN) ! clear air total reflected solar trans
       REAL  TSURF         ! surface temperature
       REAL   EMIS(MXCHAN) ! chan surface emissivity
       REAL RHOSUN(MXCHAN) ! chan reflectivity for sun
       REAL RHOTHR(MXCHAN) ! chan reflectivity for downwelling thermal
       REAL    RAD(MXCHAN) ! chan radiance
C      For clear/cloudy radiances
       REAL   RAD0         ! radiance no clouds
       REAL  RADC1         ! radiance cloud1
       REAL  RADC2         ! radiance cloud2
       REAL RADC12         ! radiance cloud1+cloud2
C
C      for RDSUN
       REAL   HSUN(MXCHAN) ! sun radiance (direct from sun)
C
C      Other variables for the sun
       REAL SUNCOS         ! cosine of sun zenith angle
       REAL COSDAZ         ! cosine(solazi - satazi) {COS Delta AZimuth}
       REAL TAUZSN(MAXLAY,MXCHAN) ! sun space-to-surface-to-space OD
       LOGICAL DOSUN       ! do sun calc?
C
C      Basic cloud info
       REAL CEMIS1            ! cloud1 emissivity
       REAL CEMIS2            ! cloud2 emissivity
       REAL CFRAC1            ! cloud1(total) fraction of FOV
       REAL CFRAC2            ! cloud2(total) fraction of FOV
       REAL CFRA1X            ! cloud1(exclusively) fraction of FOV
       REAL CFRA2X            ! cloud2(exclusively) fraction of FOV
       REAL CFRA12            ! cloud1+2(both) fraction of FOV
       REAL CPRBO1            ! cloud1 bottom pressure
       REAL CPRBO2            ! cloud2 bottom pressure
       REAL CPRTO1            ! cloud1 top pressure
       REAL CPRTO2            ! cloud2 top pressure
       REAL  CRHO1            ! cloud1 reflectivity
       REAL  CRHO2            ! cloud2 reflectivity
       REAL FCLEAR            ! clear (no cloud) fraction of FOV
       REAL TEMPC1            ! cloud1 frac layer (above cloud) mean temp
       REAL TEMPC2            ! cloud2 frac layer (above cloud) mean temp
C
C      for GETMIE
       LOGICAL LBLAC1  ! black cloud1? {Mie cloud if false}
       LOGICAL LBLAC2  ! black cloud2? {Mie cloud if false}
C
C      for CCPREP cloud1
       INTEGER LCBOT1         ! layer containing cloud bottom
       INTEGER LCTOP1         ! layer containing cloud top
       REAL  CLRT1            ! frac of layer at top of cloud clear
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
       REAL  CLRT2            ! frac of layer at top of cloud clear
       REAL TCTOP2            ! temperature at cloud top
       REAL MASEC2            ! mean cloud view angle secant
       REAL MASUN2            ! mean cloud sun-only angle secant
       REAL CFRCL2(MAXLAY)    ! fraction of cloud in layer
       REAL G_ASY2(MXCHAN)    ! "g" asymmetry
       REAL NEXTO2(MXCHAN)    ! nadir extinction optical depth
       REAL NSCAO2(MXCHAN)    ! nadir scattering optical depth
C
C      used locally only
       INTEGER      I      ! loop counter
       INTEGER      L      ! loop counter
       REAL RJUNK1         ! junk/work
       REAL C1V3           ! rad constant c1 times freq^3
       REAL C2V            ! rad constant c2 times freq
       REAL VSTORE(6)      ! temporary storage for various variables
C
C      Profile data structure
C
C      Boundary pressure levels
       COMMON /COMLEV/ PLEV
       REAL PLEV(MAXLAY+1)
C
C      for function QIKEXP
       REAL QIKEXP

C      FOR PING YANGs CODE
       REAL PLAY(MAXLAY)   ! layer mean pressure
       REAL SOLARTOA

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none
C
C***********************************************************************
C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
C
       DO I=1,NCHAN

C        Radiation constants for current channel
         C1V3=C1*(FREQ(I)**3)
         C2V=C2*FREQ(I)

C        Calculate Planck & clear airs trans for full layers
         DO L=1,LBOT-1
            RPLNCK(L)=C1V3/( EXP( C2V/TEMP(L) ) - 1.0 )
            TRANL(L)=QIKEXP( -TAU(L,I) )
         ENDDO
C        Note: TEMP(LBOT) already adjusted for bottom fractional layer
         RPLNCK(LBOT)=C1V3/( EXP( C2V/TEMP(LBOT) ) - 1.0 )

C        Calculate clear airs trans for bottom fractional layer
         RJUNK1=-TAU(LBOT,I)*BLMULT
         TRANL(LBOT)=QIKEXP( RJUNK1 )
         TRANL(LBOT)=QIKEXP( RJUNK1 )
         TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LBOT-1,I) )
         TRANS(I)=QIKEXP( BLMULT*(TAUZSN(LBOT-1,I)-TAUZSN(LBOT,I)) -
     $      TAUZSN(LBOT-1,I) )

C        Planck for surface
         RSURFE=EMIS(I)*C1V3/( EXP( C2V/TSURF ) - 1.0 )

C        Calculate clear radiance
         IF (FCLEAR .GT. 0.0) THEN
            CALL CALRAD0( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $      TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $       RHOTHR, LABOVE, COEFF, RAD0 )
          ELSE
            RAD0=0.0
         ENDIF

C        Store original values
         VSTORE(1)=TRANL(LCTOP2)
         VSTORE(2)=TRANZ(I)
         VSTORE(3)=TRANS(I)
         VSTORE(4)=RHOTHR(I)
         VSTORE(5)=RHOSUN(I)
         VSTORE(6)=RPLNCK(LCTOP2)
C        Updates for new surface if bottom cloud2 is black
         IF (CFRAC2 .GT. 0.0 .AND. LBLAC2) THEN
            RJUNK1=-TAU(LCTOP2,I)*CLRT2
            TRANL(LCTOP2)=QIKEXP( RJUNK1 )
            TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LCTOP2-1,I) )
            TRANS(I)=QIKEXP( CLRT2*(TAUZSN(LCTOP2-1,I)-
     $            TAUZSN(LCTOP2,I)) - TAUZSN(LCTOP2-1,I) )
            RSURFC=CEMIS2*C1V3/( EXP( C2V/TCTOP2 ) - 1.0 )
            RHOTHR(I)=CRHO2
            RHOSUN(I)=CRHO2
               RPLNCK(LCTOP2)=C1V3/( EXP( C2V/TEMPC2 ) - 1.0 )
            ENDIF

C        Calculate bottom cloud2 radiance
         IF (CFRA2X .GT. 0.0) THEN
            IF (LBLAC2) THEN
               CALL CALRAD0( DOSUN, I, LCTOP2, RPLNCK, RSURFC, 
     $         SECANG, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $         RHOTHR, LABOVE, COEFF, RADC2 )
            ELSE
               CALL CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $         TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $         RHOTHR, LABOVE, COEFF, CFRCL2, MASEC2, MASUN2, COSDAZ,
     $         NEXTO2, NSCAO2, G_ASY2, LCTOP2, LCBOT2, RADC2 )
            ENDIF
         ELSE
            RADC2=0.0
         ENDIF

C        Calculate combined cloud1+cloud2 radiance
         IF (CFRA12 .GT. 0.0) THEN
            IF (LBLAC2) THEN
               CALL CALRAD1( 
     $            DOSUN, I, LCTOP2, RPLNCK, RSURFC, SECANG,
     $            TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $            RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, 
     $            COSDAZ, NEXTO1, NSCAO1, G_ASY1, LCTOP1, LCBOT1, 
     $            RADC12 )
            ELSE
               CALL CALRAD2( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $         TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $         RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, NEXTO1,
     $         NSCAO1, G_ASY1, LCTOP1, LCBOT1, CFRCL2, MASEC2, 
     $         MASUN2, COSDAZ, NEXTO2, NSCAO2, G_ASY2, 
     $         LCTOP2, LCBOT2, RADC12 )
            ENDIF
         ELSE
            RADC12=0.0
         ENDIF

C        Restore original values
         TRANL(LCTOP2)=VSTORE(1)
         TRANZ(I)=VSTORE(2)
         TRANS(I)=VSTORE(3)
         RHOTHR(I)=VSTORE(4)
         RHOSUN(I)=VSTORE(5)
         RPLNCK(LCTOP2)=VSTORE(6)
C        Updates for new surface if top cloud1 is black
         IF (CFRAC1 .GT. 0.0 .AND. LBLAC1) THEN
            RJUNK1=-TAU(LCTOP1,I)*CLRT1
            TRANL(LCTOP1)=QIKEXP( RJUNK1 )
            TRANZ(I)=QIKEXP( RJUNK1 - TAUZ(LCTOP1-1,I) )
            TRANS(I)=QIKEXP( CLRT1*(TAUZSN(LCTOP1-1,I)-
     $            TAUZSN(LCTOP1,I)) - TAUZSN(LCTOP1-1,I) )
            RSURFC=CEMIS1*C1V3/( EXP( C2V/TCTOP1 ) - 1.0 )
            RHOTHR(I)=CRHO1
            RHOSUN(I)=CRHO1
            RPLNCK(LCTOP1)=C1V3/( EXP( C2V/TEMPC1 ) - 1.0 )
         ENDIF

C        Calculate top cloud1 radiance
         IF (CFRA1X .GT. 0.0) THEN
            IF (LBLAC1) THEN
               CALL CALRAD0( DOSUN, I, LCTOP1, RPLNCK, RSURFC, 
     $         SECANG, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $         RHOTHR, LABOVE, COEFF, RADC1 )
            ELSE
               CALL CALRAD1( DOSUN, I, LBOT, RPLNCK, RSURFE, SECANG,
     $         TAU, TRANL, TRANZ, SUNFAC, HSUN, TRANS, RHOSUN,
     $         RHOTHR, LABOVE, COEFF, CFRCL1, MASEC1, MASUN1, COSDAZ,
     $         NEXTO1, NSCAO1, G_ASY1, LCTOP1, LCBOT1, RADC1 )
            ENDIF
         ELSE
            RADC1=0.0
         ENDIF

ccc this block for testing
c       if (I .EQ. 1) then
c       print *,'chan1 rad0,radc1,radc2,radc12=',
c     $ RAD0,RADC1,RADC2,RADC12
c       endif
ccc

C        Total the clear & various cloudy radiances
         RAD(I)=RAD0*FCLEAR + RADC1*CFRA1X + RADC2*CFRA2X +
     $      RADC12*CFRA12

c         print *,RAD0,RADC1,RADC2,RADC12,FCLEAR,CFRA1X,CFRA2X,CFRA12
c         stop
         ENDDO ! channels

       RETURN
       END
