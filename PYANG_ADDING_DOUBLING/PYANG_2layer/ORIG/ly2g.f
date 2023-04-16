c------------------------------------------------------------------------------------------
c NAME:  two_layer_model
c
c PURPOSE:
c        Computes 587-2350 cm-1 radiance and brightness temperature at top of 
c        atmosphere for a scattering atmosphere with at most two cloud layers
c
c LANGUAGE: Fortran 95
c
c CREATION HISTORY:
c        Feb 2009    - turned into F77 by Sergio DeSouza-Machado, UMBC
c
c        16 May 2006 - Corrected error when testing whether to compute downwelling
c                      radiance at surface using 2 quadrature points
c
c         7 Apr 2006 - Added Leslie Moy's improvement in computing the total 
c                      downwelling flux at the surface in clear sky conditions 
c                      only (referred to in earlier code as DFLUX)
c
c         6 Apr 2006 - Incorporated new ice R/T databases based on Baum's bulk ice 
c                      scattering property models
c
c         3 Apr 2006 - Added variable file format for opening coefficient files so 
c                      compatible with pgf90, g95 and Intel f90 compilers. Also added
c                      error handling for load_rt_tables.
c
c         6 Mar 2006 - Previous modification was wrong. Changed keyword for FORM in
c                      OPEN statements back to 'UNFORMATTED'.
c
c         2 Mar 2006 - Changed FORM in OPEN statements in load_rt_tables from
c                      'UNFORMATTED' to 'BINARY' in order to be compatible with different
c                      compilers
c
c        28 Feb 2006 - Significantly speeded up R/T table interpolation; modified 
c                      routines cldice, cldwater, intpth, and intpmu
c
c        27 Feb 2006 - Removed routines get_cloud_prop, get_CloudsMM5 and 
c                      parseClouds, which now are in a separate module (cloud_properties);
c                      also removed routines planck and ebbt, which are in module planckf
c
c        13 Dec 2005 - Loads R/T tables into memory. Argument itab was
c                      removed from routine ly2g
c
c        28 Nov 2005 - Modified by TG to use same source function integration 
c                      scheme as the SOI model 
c
c        Written by Tom Greenwald  UW/CIMSS  19 Oct - 28 Oct 2005
c------------------------------------------------------------------------------------------

      MODULE two_layer_model

c---------------
c Module usage
c---------------
       USE type_kinds, ONLY : fp_kind
       USE SOI_thermal_source
       USE planckf
       USE open_file_format

        IMPLICIT NONE
  
c        CHARACTER*11, PUBLIC, PARAMETER :: FILE_FORMAT = 'BINARY     ' c 'UNFORMATTED' for pgf90 & g95 
                                                                       c 'BINARY     ' for Intel f90

        INTEGER, PUBLIC, PARAMETER                ::  THETMAX = 9      c Number of zenith angles
        INTEGER, PUBLIC, PARAMETER                ::  NUM_LIQ = 2938       c Number of wavenumbers
        INTEGER, PUBLIC, PARAMETER                ::  NUM_ICE = 3151       c Number of wavenumbers
        INTEGER, PUBLIC, PARAMETER                ::  TAUMAX_ICE = 25  c No. of ice optical depths
        INTEGER, PUBLIC, PARAMETER                ::  DMAX_ICE = 18    c No. of ice effective diameters
        INTEGER, PUBLIC, PARAMETER                ::  TAUMAX_LIQ = 22  c No. of liquid optical depths
        INTEGER, PUBLIC, PARAMETER                ::  DMAX_LIQ = 13    c No. of liquid effective diameters
        INTEGER, PUBLIC, PARAMETER                ::  NLEVEL_MAX = 101
        INTEGER, PUBLIC, PARAMETER                ::  NSPEC_MAX = 5    c Max number of hydrometeor species
        INTEGER, PUBLIC, PARAMETER                ::  NANG_MAX = 3
        REAL( fp_kind ), PUBLIC, PARAMETER        ::  ZERO = 0.0_fp_kind
        REAL( fp_kind ), PUBLIC, PARAMETER        ::  ONE = 1.0_fp_kind
        REAL( fp_kind ), PUBLIC, PARAMETER        ::  TWO = 2.0_fp_kind
        REAL( fp_kind ), PUBLIC, PARAMETER        ::  PI = 3.14159265_fp_kind

        INTEGER, PUBLIC :: TAU_idx
        INTEGER, PUBLIC :: TAU_idx_
        INTEGER, PUBLIC :: D_idx
        INTEGER, PUBLIC :: D_idx_
        INTEGER, PUBLIC :: THET_idx
        INTEGER, PUBLIC :: THET_idx_
        REAL( fp_kind ), PUBLIC, DIMENSION( THETMAX ) :: THET
        REAL, PUBLIC, DIMENSION( 2 * THETMAX, NUM_ICE, TAUMAX_ICE * DMAX_ICE ) :: RT_ICE
        REAL, PUBLIC, DIMENSION( 2 * THETMAX, NUM_LIQ, TAUMAX_LIQ * DMAX_LIQ ) :: RT_LIQ

c----------------
c Visibilities
c----------------

        PUBLIC  :: ly2g, &
                   LayerTemp, &
                   load_rt_tables
                   
        PRIVATE :: cldice, &
                   cldwater, &
                   prod, &
                   trp, &
                   intpth, &
                   intpmu

      CONTAINS
c-----------------------------------------------------------------------
c
c  This is a significant rewrite of the 2-layer multiple scattering 
c  model (Niu et al. 2005) as a callable subroutine so it can more 
c  easily be "plugged into" other code and accept arbitrary pressure 
c  level profiles and gas optical depth profiles.
c  Model is limited to between 587.3 and 2349.5 cm-1 due to limits of 
c  R/T function tables.
c
c  INPUT:
c
c     nlev      :   Number of atmospheric levels
c     temp      :   Profile of temperatures at levels (Units of Kelvin)
c     gtran     :   Profile of layer gas transmittances
c     gtau      :   Profile of layer gas optical depths (slant path)
c     tauvis4   :   Visible optical depth for upper cloud (layer 4)
c     de4       :   Effective particle diameter for upper cloud (layer 4)
c     ncld4     :   Level number for upper cloud top
c     tauvis2   :   Visible optical depth for lower cloud (layer 2)
c     de2       :   Effective particle diameter for lower cloud (layer 2)
c     ncld2     :   Level number for lower cloud top
c     cld_indt  :   Cloud type identifier
c                   (0: Clear; 1: single-layer ice; 2: single-layer liquid;
c                    3: ice + liquid; 4: ice + ice)
c     waven     :   Wavenumber of incident radiation (inverse cm)
c     tsfc      :   Surface temperature (Units of Kelvin)
c     emiss     :   Surface emissivity  (Unitless)
c     sfc_type  :   Surface type ( 1: Water, 2: Land ) 
c     zen       :   Zenith angle for outgoing radiance (Units of 0-90 deg.)
c     
c  OUTPUT:
c
c     rad       :   Radiance at top of atmosphere (Units of Watts / 
c                                                  meter**2 steradian inverse cm) 
c     tb        :   Brightness temperature at top of atmosphere (Units of Kelvin)
c
c-----------------------------------------------------------------------
c  References:
c
c     Niu, J., P. Yang, H.-L. Huang, J. E. Davies, J. Li, B. Baum, and 
c          Y. X. Hu, 2005: A fast infrared radiative transfer model for 
c          overlapping clouds. To appear in Journal of Quantitative 
c          Spectroscopy & Radiative Transfer.
c
c
c  Written by Tom Greenwald  UW/CIMSS   10/19/2005 - 10/28/2005
c-----------------------------------------------------------------------
        
      SUBROUTINE ly2g( nlev,      !  c Input
                       temp,      !  c Input
                       gtau,      !  c Input
                       nang,      !  c Input
                       qwgt,      !  c Input
                       tauvis4,   !  c Input
                       de4,       !  c Input
                       ncld4,     !  c Input/output
                       tauvis2,   !  c Input
                       de2,       !  c Input
                       ncld2,     !  c Input/output
                       cld_indt,  !  c Input
                       waven,     !  c Input
                       tsfc,      !  c Input
                       emiss,     !  c Input
                       sfc_type,  !  c Input
                       zen,       !  c Input
                       rad,       !  c Output
                       tb     )    ! c Output

c------------
c Arguments
c------------

        INTEGER nlev
        REAL temp(*), gtau(*), qwgt(*)
        INTEGER nang,sfc_type
        REAL tauvis2,de2,tauvis4,de4,
        INTEGER nld2,ncld4,cld_indt,waven,tsdc,emiss,zen,rad,tb

c-------------------
c Local variables
c-------------------
        REAL T1( nang_max )
        REAL T2( nang_max )
        REAL T3( nang_max )
        REAL T4( nang_max )
        REAL T5
        REAL T15
        REAL T25
        REAL T35
        REAL T45
        REAL TCL
        REAL TCH
        REAL RC2
        REAL RC4
        REAL Itop
        REAL ISC
        REAL TBItop
        REAL I5inf
        REAL Iinf5(NANG_MAX )
        REAL I34
        REAL I43(NANG_MAX)
        REAL I23
        REAL I32(NANG_MAX)
        REAL I12
        REAL I21(NANG_MAX)
        REAL I01
        REAL I10(NANG_MAX)
        REAL BS
        REAL BDair(NLEVEL_MAX - 1)
        REAL reflect
        REAL TRANSMI_C( NLEVEL_MAX - 1, NANG_MAX )
        REAL slopd, B_TOP, B_BOT, ssa, B_AVG, B0, B1, mu
        REAL B_UP(NLEVEL_MAX - 1, NANG_MAX), B_DN(NLEVEL_MAX - 1, NANG_MAX)
        REAL tauir2, tauir4
        REAL zang
        REAL downwelling_tot
        REAL gtran( NLEVEL_MAX, NANG_MAX )
        DOUBLE PRECISION ff              
        INTEGER i, j, k, l
        INTEGER nlayer
        CHARACTER*256 FNCIRUS, FNWATER
        LOGICAL first
        DATA first / .true. /

        REAL( fp_kind ), PARAMETER :: ONE_EIGHTY = 180.0_fp_kind
        REAL( fp_kind ), PARAMETER :: ZEROPT5 = 0.5_fp_kind
        
c--------------
c  Range check
c--------------
        IF ( zen > 80.0 ) THEN
c         print *, ' WARNING: Input zenith angle is too large in ly2g, setting to 80 degrees'
          zang = 80.0 ! c Maximum zenith angle in RT tables
        ELSE
          zang = zen 
        END IF
c--------------------------------------------------
c Check for consistency between nang and cld_indt
c If problem, reset nang
c--------------------------------------------------
c        IF ( nang /= 3 .AND. cld_indt == 0 ) THEN
c           print *, ' ERROR - inconsistency between nang and cld_indt' 
c           nang = 3
c           STOP
c        ENDIF

c-------------------------
c Compute transmittances
c-------------------------
        DO i = 1, nlev
          DO j = 1, nang
            gtran( i, j ) = exp( -gtau( i, j ) )
          END DO
        END DO

c ----------------------------------------------------------------
c      INITIALIZATION
c ----------------------------------------------------------------
        nlayer = nlev - 1

c---------------------------------
c Let's do the radiative transfer
c---------------------------------
c Initialize layer reflectances and transmittances
        RC2 = ZERO   ! c Lower cloud layer reflectance
        RC4 = ZERO  ! c Upper cloud layer reflectance
        TBItop = ZERO! c Downwelling radiance at top of atmosphere
        TCL = ONE   ! c Lower cloud layer transmittance
        TCH = ONE   ! c Upper cloud layer transmittance
        T1 =  ONE
        T2  = ONE
        T3  = ONE
        T4  = ONE
        T5  = ONE
        T15 = ONE
        DO j = 1, nlayer
          DO i = 1, nang
            TRANSMI_C( j, i ) = gtran( j, i )
          END DO
        END DO
c-------------------------------------------------------------------
c Note:  Even when there are no cloud layers, or where there
c        is only 1 cloud layer, NCLD2 and NCLD4 must be set to
c        values which won't screw up calculations with the logic
c        of the loop ranges in subsequent code.  JED 20-Jan-2005.
c-------------------------------------------------------------------
        IF ( cld_indt == 4 ) THEN

          CALL cldice( waven, TAUvis4, De4, zang, RC4, TCH )
          CALL cldice( waven, TAUvis2, De2, zang, RC2, TCL )
          TAUIR2 = -LOG( TCL )
          TAUIR4 = -LOG( TCH )

        ELSE IF( cld_indt == 3 ) THEN

          CALL cldice( waven, TAUvis4, De4, zang, RC4, TCH )
          CALL cldwater( waven, TAUvis2, De2, zang, RC2, TCL)

          TAUIR2 = -LOG( TCL )
          TAUIR4 = -LOG( TCH )

        ELSE IF( cld_indt == 2 ) THEN

          CALL cldwater( waven, TAUvis2, De2, zang, RC2, TCL )  
          NCLD4 = NCLD2 - 1 ! c Dummy value
          TAUIR2 = -LOG( TCL )
          TAUIR4 = ZERO

        ELSE IF( cld_indt == 1 ) THEN

          CALL cldice( waven, TAUvis4, De4, zang, RC4, TCH )
          NCLD2 = NCLD4 + 1 ! c Dummy value
          TAUIR4 = -LOG( TCH )
          TAUIR2 = ZERO

        ELSE IF( cld_indt == 0 ) THEN

          NCLD4 = 2 ! c Dummy value
          NCLD2 = 4 ! c Dummy value
          TAUIR2 = ZERO
          TAUIR4 = ZERO

        ENDIF

c------------------------
c   Inserting cloud
c------------------------
40      DO i = 1, nang
          TRANSMI_C( NCLD4, i ) = TCH * TRANSMI_C( NCLD4, i )
          TRANSMI_C( NCLD2, i ) = TCL * TRANSMI_C( NCLD2, i )
        ENDDO

c----------------------------------------------------------
c  Compute atmospheric transmittance for each of 5 layers
c----------------------------------------------------------
        T5 = PROD( TRANSMI_C( :, 1), nlayer, 1 , NCLD4 - 1 )
        DO i = 1, nang
          T4( i ) = TRANSMI_C( NCLD4, i )
          T3( i ) = PROD( TRANSMI_C( :, i), nlayer, NCLD4 + 1, NCLD2 - 1 )
          T2( i ) = TRANSMI_C( NCLD2, i )
          T1( i ) = PROD( TRANSMI_C( :, i), nlayer, NCLD2 + 1, NLAYER )
        END DO

        T45 = T4( 1 ) * T5
        T35 = T3( 1 ) * T45
        T25 = T2( 1 ) * T35
        T15 = T1( 1 ) * T25

c--------------------------------------
c   Surface emission contribution
c--------------------------------------
        ff = waven
        CALL planck( ff, tsfc, bs )

c---------------------------------------
c  Compute atmospheric thermal sources
c---------------------------------------
        ssa = ZERO
        mu = cos( PI * zang / ONE_EIGHTY )
        DO i = 1, NLAYER
          CALL PLANCK( FF, TEMP( i ), B_TOP )
          CALL PLANCK( FF, TEMP( i + 1 ), B_BOT )
          DO j = 1, nang
            IF ( i == NCLD2 ) THEN
              slopd = gtau( i, j ) + TAUIR2 / mu
            ELSE IF ( i == NCLD4 ) THEN
              slopd = gtau( i, j ) + TAUIR4 / mu
            ELSE
              slopd = gtau( i, j )
            END IF
            call thermal_source( slopd, ssa, TRANSMI_C( i, j ), B_TOP, B_BOT, B_UP( i, j ), B_DN( i, j ) )
          END DO
        END DO

c-----------------------
c Surface emission term 
c-----------------------
        ISC = emiss * BS * PROD( TRANSMI_C( :, 1), nlayer, 1, NLAYER )

c---------------------------------------------------
c Compute upwelling radiances at layer interfaces
c---------------------------------------------------
        I01 = ZERO
        DO j = NCLD2 + 1, NLAYER          
           I01 = I01 + B_UP( j, 1 ) * PROD( TRANSMI_C( :, 1 ), nlayer, NCLD2 + 1, j - 1 )
        END DO
        CALL PLANCK( FF, TEMP( NCLD2  ), B0 ) 
        CALL PLANCK( FF, TEMP( NCLD2 + 1 ), B1 )
        B_AVG = ZEROPT5 * ( B0 + B1 )                
        I12 = B_UP( NCLD2, 1 ) - B_AVG * RC2
        DO i = 1, nang
          I21( i ) = B_DN( NCLD2, i ) - B_AVG * RC2
        END DO
        I23 = ZERO
        DO j = NCLD4 + 1, NCLD2 - 1
          I23 = I23 + B_UP( j, 1 ) * PROD( TRANSMI_C( :, 1 ), nlayer, NCLD4 + 1, j - 1 )
        END DO
        CALL PLANCK( FF, TEMP( NCLD4  ), B0 ) 
        CALL PLANCK( FF, TEMP( NCLD4 + 1 ), B1 )
        B_AVG = ZEROPT5 * ( B0 + B1 )                
        I34 = B_UP( NCLD4, 1 ) - B_AVG * RC4
        DO i = 1, nang
          I43( i ) = B_DN( NCLD4, i ) - B_AVG * RC4
        END DO
        I5inf = ZERO
        DO j = 1, NCLD4 - 1
          I5inf = I5inf + B_UP( j, 1 ) * PROD( TRANSMI_C( :, 1 ), nlayer, 1, j - 1 )
        END DO

c------------------------------
c   Compute I10, I32, I54
c------------------------------
        DO i = 1, nang
          Iinf5( i ) = ZERO
          DO j = 1, NCLD4 - 1
             Iinf5( i ) = Iinf5( i ) + B_DN( j, i ) * PROD( TRANSMI_C( :, i ), nlayer, j + 1, NCLD4 - 1 )
          END DO
        END DO
        DO i = 1, nang
          I32( i ) = ZERO
          DO j = NCLD4 + 1, NCLD2 - 1
            I32( i ) = I32( i ) + B_DN( j, i ) * PROD( TRANSMI_C( :, i ), nlayer, j + 1, NCLD2 - 1 ) 
          END DO
        END DO
        DO i = 1, nang
          I10( i ) = ZERO
          DO j = NCLD2 + 1, nlayer
            I10( i ) = I10( i ) + B_DN( j, i ) * PROD( TRANSMI_C( :, i ), nlayer, j + 1, NLAYER )
          END DO
        END DO

        IF ( sfc_type == 1 ) THEN
          reflect = ONE - emiss c specular (water) reflectance
        ELSE
          reflect = ( ONE - emiss ) / PI c rough surface (land) reflectance
        END IF

c---------------------------------------------------------------------------------------------
c If clear-sky profile, then compute downwelling flux at surface using two-point quadrature 
c Adopted from Leslie Moy's approach in gifts_dflux routine
c---------------------------------------------------------------------------------------------
        IF ( cld_indt == 0 .and. nang == 3 ) THEN
          downwelling_tot = ZERO
          DO i = 2, nang
            downwelling_tot =  downwelling_tot + ( I10( i ) + I21( i ) * T1( i ) + I32( i ) * T1( i ) &
                               * T2( i ) + I43( i ) * T1( i ) * T2( i ) * T3( i ) + &
                               Iinf5( i ) * T1( i ) * T2( i ) * T3( i ) * T4( i ) ) * qwgt( i - 1 )
          END DO
          downwelling_tot = TWO * pi * downwelling_tot
        ELSE 
          downwelling_tot =  I10( 1 ) + I21( 1 ) * T1( 1 ) + I32( 1 ) * T1( 1 ) * T2( 1 ) + I43( 1 ) &
                             * T1( 1 ) * T2( 1 ) * T3( 1 ) + Iinf5( 1 ) * T1( 1 ) * T2( 1 ) * T3( 1 ) * T4( 1 ) 
        END IF       
c---------------------------------------------------------
c  Add up all terms to get radiance at top of atmosphere
c---------------------------------------------------------
        Itop = ISC + &
               I01 * T25 + &
               I12 * T35 + &
               I23 * T45 + &
               I34 * T5 + &
               I5inf + &
               RC4 * Iinf5( 1 ) * T5 + &
               RC2 * T35 * ( I32( 1 ) + I43( 1 ) * T3( 1 ) + Iinf5( 1 ) * T3( 1 ) * T4( 1 ) ) + &
               reflect * T15 * downwelling_tot

        rad = Itop c Units of Watts / meter**2 steradian inverse cm) 

c----------------------------------------------------------------
c Convert to equivalent Black body temperature (Units of Kelvin)
c----------------------------------------------------------------
        call ebbt( ff, Itop, tb )


   END SUBROUTINE ly2g

