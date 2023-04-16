!------------------------------------------------------------------------------------------
! NAME:  two_layer_model
!
! PURPOSE:
!        Computes 587-2350 cm-1 radiance and brightness temperature at top of 
!        atmosphere for a scattering atmosphere with at most two cloud layers
!
! LANGUAGE: Fortran 95
!
! CREATION HISTORY:
!
!
!        16 May 2006 - Corrected error when testing whether to compute downwelling
!                      radiance at surface using 2 quadrature points
!
!         7 Apr 2006 - Added Leslie Moy's improvement in computing the total 
!                      downwelling flux at the surface in clear sky conditions 
!                      only (referred to in earlier code as DFLUX)
!
!         6 Apr 2006 - Incorporated new ice R/T databases based on Baum's bulk ice 
!                      scattering property models
!
!         3 Apr 2006 - Added variable file format for opening coefficient files so 
!                      compatible with pgf90, g95 and Intel f90 compilers. Also added
!                      error handling for load_rt_tables.
!
!         6 Mar 2006 - Previous modification was wrong. Changed keyword for FORM in
!                      OPEN statements back to 'UNFORMATTED'.
!
!         2 Mar 2006 - Changed FORM in OPEN statements in load_rt_tables from
!                      'UNFORMATTED' to 'BINARY' in order to be compatible with different
!                      compilers
!
!        28 Feb 2006 - Significantly speeded up R/T table interpolation; modified 
!                      routines cldice, cldwater, intpth, and intpmu
!
!        27 Feb 2006 - Removed routines get_cloud_prop, get_CloudsMM5 and 
!                      parseClouds, which now are in a separate module (cloud_properties);
!                      also removed routines planck and ebbt, which are in module planckf
!
!        13 Dec 2005 - Loads R/T tables into memory. Argument itab was
!                      removed from routine ly2g
!
!        28 Nov 2005 - Modified by TG to use same source function integration 
!                      scheme as the SOI model 
!
!        Written by Tom Greenwald  UW/CIMSS  19 Oct - 28 Oct 2005
!------------------------------------------------------------------------------------------

      MODULE two_layer_model

!---------------
! Module usage
!---------------
       USE type_kinds, ONLY : fp_kind
       USE SOI_thermal_source
       USE planckf
       USE open_file_format

        IMPLICIT NONE
  
!        CHARACTER*11, PUBLIC, PARAMETER :: FILE_FORMAT = 'BINARY     ' ! 'UNFORMATTED' for pgf90 & g95 
                                                                       ! 'BINARY     ' for Intel f90

        INTEGER, PUBLIC, PARAMETER                ::  THETMAX = 9      ! Number of zenith angles
        INTEGER, PUBLIC, PARAMETER                ::  NUM_LIQ = 2938       ! Number of wavenumbers
        INTEGER, PUBLIC, PARAMETER                ::  NUM_ICE = 3151       ! Number of wavenumbers
        INTEGER, PUBLIC, PARAMETER                ::  TAUMAX_ICE = 25  ! No. of ice optical depths
        INTEGER, PUBLIC, PARAMETER                ::  DMAX_ICE = 18    ! No. of ice effective diameters
        INTEGER, PUBLIC, PARAMETER                ::  TAUMAX_LIQ = 22  ! No. of liquid optical depths
        INTEGER, PUBLIC, PARAMETER                ::  DMAX_LIQ = 13    ! No. of liquid effective diameters
        INTEGER, PUBLIC, PARAMETER                ::  NLEVEL_MAX = 101
        INTEGER, PUBLIC, PARAMETER                ::  NSPEC_MAX = 5    ! Max number of hydrometeor species
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

!----------------
! Visibilities
!----------------

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
!-----------------------------------------------------------------------
!
!  This is a significant rewrite of the 2-layer multiple scattering 
!  model (Niu et al. 2005) as a callable subroutine so it can more 
!  easily be "plugged into" other code and accept arbitrary pressure 
!  level profiles and gas optical depth profiles.
!  Model is limited to between 587.3 and 2349.5 cm-1 due to limits of 
!  R/T function tables.
!
!  INPUT:
!
!     nlev      :   Number of atmospheric levels
!     temp      :   Profile of temperatures at levels (Units of Kelvin)
!     gtran     :   Profile of layer gas transmittances
!     gtau      :   Profile of layer gas optical depths (slant path)
!     tauvis4   :   Visible optical depth for upper cloud (layer 4)
!     de4       :   Effective particle diameter for upper cloud (layer 4)
!     ncld4     :   Level number for upper cloud top
!     tauvis2   :   Visible optical depth for lower cloud (layer 2)
!     de2       :   Effective particle diameter for lower cloud (layer 2)
!     ncld2     :   Level number for lower cloud top
!     cld_indt  :   Cloud type identifier
!                   (0: Clear; 1: single-layer ice; 2: single-layer liquid;
!                    3: ice + liquid; 4: ice + ice)
!     waven     :   Wavenumber of incident radiation (inverse cm)
!     tsfc      :   Surface temperature (Units of Kelvin)
!     emiss     :   Surface emissivity  (Unitless)
!     sfc_type  :   Surface type ( 1: Water, 2: Land ) 
!     zen       :   Zenith angle for outgoing radiance (Units of 0-90 deg.)
!     
!  OUTPUT:
!
!     rad       :   Radiance at top of atmosphere (Units of Watts / 
!                                                  meter**2 steradian inverse cm) 
!     tb        :   Brightness temperature at top of atmosphere (Units of Kelvin)
!
!-----------------------------------------------------------------------
!  References:
!
!     Niu, J., P. Yang, H.-L. Huang, J. E. Davies, J. Li, B. Baum, and 
!          Y. X. Hu, 2005: A fast infrared radiative transfer model for 
!          overlapping clouds. To appear in Journal of Quantitative 
!          Spectroscopy & Radiative Transfer.
!
!
!  Written by Tom Greenwald  UW/CIMSS   10/19/2005 - 10/28/2005
!-----------------------------------------------------------------------
        
      SUBROUTINE ly2g( nlev,      &  ! Input
                       temp,      &  ! Input
!                       gtran,     &  ! Input
                       gtau,      &  ! Input
                       nang,      &  ! Input
                       qwgt,      &  ! Input
                       tauvis4,   &  ! Input
                       de4,       &  ! Input
                       ncld4,     &  ! Input/output
                       tauvis2,   &  ! Input
                       de2,       &  ! Input
                       ncld2,     &  ! Input/output
                       cld_indt,  &  ! Input
                       waven,     &  ! Input
                       tsfc,      &  ! Input
                       emiss,     &  ! Input
                       sfc_type,  &  ! Input
                       zen,       &  ! Input
                       rad,       &  ! Output
                       tb     )      ! Output

!------------
! Arguments
!------------

        INTEGER,                            INTENT( IN )    :: nlev
        REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )    :: temp
!        REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )    :: gtran
        REAL( fp_kind ), DIMENSION( :, : ),    INTENT( IN )    :: gtau
        INTEGER,                            INTENT( INOUT )    :: nang
        REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )    :: qwgt        
        REAL( fp_kind ),                    INTENT( IN )    :: tauvis2
        REAL( fp_kind ),                    INTENT( IN )    :: de2
        INTEGER,                            INTENT( INOUT ) :: ncld2                 
        REAL( fp_kind ),                    INTENT( IN )    :: tauvis4
        REAL( fp_kind ),                    INTENT( IN )    :: de4
        INTEGER,                            INTENT( INOUT ) :: ncld4
        INTEGER,                            INTENT( IN )    :: cld_indt
        REAL( fp_kind ),                    INTENT( IN )    :: waven
        REAL( fp_kind ),                    INTENT( IN )    :: tsfc
        REAL( fp_kind ),                    INTENT( IN )    :: emiss
        INTEGER,                            INTENT( IN )    :: sfc_type
        REAL( fp_kind ),                    INTENT( IN )    :: zen
        REAL( fp_kind ),                    INTENT( OUT )   :: rad
        REAL( fp_kind ),                    INTENT( OUT )   :: tb

!-------------------
! Local variables
!-------------------
        REAL( fp_kind ) :: T1( nang_max )
        REAL( fp_kind ) :: T2( nang_max )
        REAL( fp_kind ) :: T3( nang_max )
        REAL( fp_kind ) :: T4( nang_max )
        REAL( fp_kind ) :: T5
        REAL( fp_kind ) :: T15
        REAL( fp_kind ) :: T25
        REAL( fp_kind ) :: T35
        REAL( fp_kind ) :: T45
        REAL( fp_kind ) :: TCL
        REAL( fp_kind ) :: TCH
        REAL( fp_kind ) :: RC2
        REAL( fp_kind ) :: RC4
        REAL( fp_kind ) :: Itop
        REAL( fp_kind ) :: ISC
        REAL( fp_kind ) :: TBItop
        REAL( fp_kind ) :: I5inf
        REAL( fp_kind ), DIMENSION( NANG_MAX ) :: Iinf5
        REAL( fp_kind ) :: I34
        REAL( fp_kind ), DIMENSION( NANG_MAX ) :: I43
        REAL( fp_kind ) :: I23
        REAL( fp_kind ), DIMENSION( NANG_MAX ) :: I32
        REAL( fp_kind ) :: I12
        REAL( fp_kind ), DIMENSION( NANG_MAX ) :: I21
        REAL( fp_kind ) :: I01
        REAL( fp_kind ), DIMENSION( NANG_MAX ) :: I10
        REAL( fp_kind ) :: BS
        REAL( fp_kind ), DIMENSION( NLEVEL_MAX - 1 ) :: BDair
        REAL( fp_kind ) :: reflect
        REAL( fp_kind ), DIMENSION( NLEVEL_MAX - 1, NANG_MAX ) :: TRANSMI_C
        REAL( fp_kind ) :: slopd, B_TOP, B_BOT, ssa, B_AVG, B0, B1, mu
        REAL( fp_kind ), DIMENSION( NLEVEL_MAX - 1, NANG_MAX ) :: B_UP, B_DN
        REAL( fp_kind ) :: tauir2, tauir4
        REAL( fp_kind ) :: zang
        REAL( fp_kind ) :: downwelling_tot
        REAL( fp_kind ), DIMENSION( NLEVEL_MAX, NANG_MAX ) :: gtran
        REAL*8 :: ff              
        INTEGER :: i, j, k, l
        INTEGER :: nlayer
        CHARACTER( 256 ) :: FNCIRUS, FNWATER
        LOGICAL :: first
        DATA first / .true. /

        REAL( fp_kind ), PARAMETER :: ONE_EIGHTY = 180.0_fp_kind
        REAL( fp_kind ), PARAMETER :: ZEROPT5 = 0.5_fp_kind
        
!--------------
!  Range check
!--------------
        IF ( zen > 80.0 ) THEN
!         print *, ' WARNING: Input zenith angle is too large in ly2g, setting to 80 degrees'
          zang = 80.0   ! Maximum zenith angle in RT tables
        ELSE
          zang = zen 
        END IF
!--------------------------------------------------
! Check for consistency between nang and cld_indt
! If problem, reset nang
!--------------------------------------------------
!        IF ( nang /= 3 .AND. cld_indt == 0 ) THEN
!           print *, ' ERROR - inconsistency between nang and cld_indt' 
!           nang = 3
!           STOP
!        ENDIF

!-------------------------
! Compute transmittances
!-------------------------
        DO i = 1, nlev
          DO j = 1, nang
            gtran( i, j ) = exp( -gtau( i, j ) )
          END DO
        END DO

! ----------------------------------------------------------------
!      INITIALIZATION
! ----------------------------------------------------------------
        nlayer = nlev - 1

!---------------------------------
! Let's do the radiative transfer
!---------------------------------
! Initialize layer reflectances and transmittances
        RC2 = ZERO     ! Lower cloud layer reflectance
        RC4 = ZERO    ! Upper cloud layer reflectance
        TBItop = ZERO  ! Downwelling radiance at top of atmosphere
        TCL = ONE     ! Lower cloud layer transmittance
        TCH = ONE     ! Upper cloud layer transmittance
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
!-------------------------------------------------------------------
! Note:  Even when there are no cloud layers, or where there
!        is only 1 cloud layer, NCLD2 and NCLD4 must be set to
!        values which won't screw up calculations with the logic
!        of the loop ranges in subsequent code.  JED 20-Jan-2005.
!-------------------------------------------------------------------
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
          NCLD4 = NCLD2 - 1   ! Dummy value
          TAUIR2 = -LOG( TCL )
          TAUIR4 = ZERO

        ELSE IF( cld_indt == 1 ) THEN

          CALL cldice( waven, TAUvis4, De4, zang, RC4, TCH )
          NCLD2 = NCLD4 + 1   ! Dummy value
          TAUIR4 = -LOG( TCH )
          TAUIR2 = ZERO

        ELSE IF( cld_indt == 0 ) THEN

          NCLD4 = 2   ! Dummy value
          NCLD2 = 4   ! Dummy value
          TAUIR2 = ZERO
          TAUIR4 = ZERO

        ENDIF

!------------------------
!   Inserting cloud
!------------------------
40      DO i = 1, nang
          TRANSMI_C( NCLD4, i ) = TCH * TRANSMI_C( NCLD4, i )
          TRANSMI_C( NCLD2, i ) = TCL * TRANSMI_C( NCLD2, i )
        ENDDO

!----------------------------------------------------------
!  Compute atmospheric transmittance for each of 5 layers
!----------------------------------------------------------
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

!--------------------------------------
!   Surface emission contribution
!--------------------------------------
        ff = waven
        CALL planck( ff, tsfc, bs )

!---------------------------------------
!  Compute atmospheric thermal sources
!---------------------------------------
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

!-----------------------
! Surface emission term 
!-----------------------
        ISC = emiss * BS * PROD( TRANSMI_C( :, 1), nlayer, 1, NLAYER )

!---------------------------------------------------
! Compute upwelling radiances at layer interfaces
!---------------------------------------------------
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

!------------------------------
!   Compute I10, I32, I54
!------------------------------
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
          reflect = ONE - emiss ! specular (water) reflectance
        ELSE
          reflect = ( ONE - emiss ) / PI ! rough surface (land) reflectance
        END IF

!---------------------------------------------------------------------------------------------
! If clear-sky profile, then compute downwelling flux at surface using two-point quadrature 
! Adopted from Leslie Moy's approach in gifts_dflux routine
!---------------------------------------------------------------------------------------------
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
!---------------------------------------------------------
!  Add up all terms to get radiance at top of atmosphere
!---------------------------------------------------------
        Itop = ISC + &
               I01 * T25 + &
               I12 * T35 + &
               I23 * T45 + &
               I34 * T5 + &
               I5inf + &
               RC4 * Iinf5( 1 ) * T5 + &
               RC2 * T35 * ( I32( 1 ) + I43( 1 ) * T3( 1 ) + Iinf5( 1 ) * T3( 1 ) * T4( 1 ) ) + &
               reflect * T15 * downwelling_tot

        rad = Itop ! Units of Watts / meter**2 steradian inverse cm) 

!----------------------------------------------------------------
! Convert to equivalent Black body temperature (Units of Kelvin)
!----------------------------------------------------------------
        call ebbt( ff, Itop, tb )


   END SUBROUTINE ly2g

   SUBROUTINE load_rt_tables( input_filename_ice, & ! Input
                              input_filename_liq  ) ! Input

      IMPLICIT NONE
      
!------------
! Arguments
!------------
      CHARACTER( * ), INTENT( IN )  :: input_filename_ice
      CHARACTER( * ), INTENT( IN )  :: input_filename_liq

!------------------
! Local variables
!------------------
      INTEGER :: i, IOSTAT

      IOSTAT = 0

      OPEN( 16, FILE = input_filename_ice, RECL = NUM_ICE * THETMAX * 2 * 4, ACCESS = 'DIRECT', &
            FORM = TRIM( FILE_FORMAT ), STATUS = 'OLD', IOSTAT = IOSTAT )
      IF ( IOSTAT > 0 ) THEN
        print *, ' ** OPEN ERROR in load_rt_tables for file ', trim(input_filename_ice), ' *** '
        STOP
      END IF 

      DO i = 1, TAUMAX_ICE * DMAX_ICE
        READ( 16, REC = i, ERR = 99 ) RT_ICE( :, :, i )
      END DO      
      CLOSE( 16 )

      OPEN( 16, FILE = input_filename_liq, RECL = NUM_LIQ * THETMAX * 2 * 4, ACCESS = 'DIRECT', &
            FORM = TRIM( FILE_FORMAT ), STATUS = 'OLD', IOSTAT = IOSTAT )
      IF ( IOSTAT > 0 ) THEN
        print *, ' ** OPEN ERROR in load_rt_tables for file ', trim(input_filename_liq), ' *** '
        STOP
      END IF 

      DO i = 1, TAUMAX_LIQ * DMAX_LIQ
        READ( 16, REC = i, ERR = 99 ) RT_LIQ( :, :, i )
      END DO      
      CLOSE( 16 )

      RETURN
      
 99   PRINT *, ' ** READ ERROR in load_rt_tables ** '
      STOP

   END SUBROUTINE load_rt_tables

   REAL FUNCTION prod( ARR,       & ! Input
                       NCLMN,     & ! Input
                       STARTCLMN, & ! Input
                       ENDCLMN )    ! Input

!-----------
! Arguments
!-----------
       REAL( fp_kind ), DIMENSION( NCLMN ), INTENT( IN ) :: ARR
       INTEGER,                             INTENT( IN ) :: NCLMN
       INTEGER,                             INTENT( IN ) :: STARTCLMN
       INTEGER,                             INTENT( IN ) :: ENDCLMN

!------------------
! Local variables 
!------------------
       REAL( fp_kind ) :: MULTI
       INTEGER         :: i

       IF ( STARTCLMN > ENDCLMN ) THEN
         PROD = 1
         RETURN
       ELSE
         MULTI = ONE
         DO i = STARTCLMN, ENDCLMN
           MULTI = MULTI * ARR( i )
         END DO
         PROD = MULTI
       ENDIF

    END FUNCTION prod

      SUBROUTINE LayerTemp( nlevel,     &  ! Input
                            P_AIRS,     &  ! Input
                            T_AIRS,     &  ! Input
                            PSURF,      &  ! Input
                            HPA4,       &  ! Input
                            HPA2,       &  ! Input
                            temp_layer, &  ! Output
                            ncld4,      &  ! Output
                            ncld2 )        ! Output
!--------------------------------------------------------------------
!
!     Subroutine to provide layer effective temperatures from 
!     temperatures at specified pressure levels.  HPA{4,2} are 
!     used to compute the cloud layers ncld{4,2}; if HPA4 or HPA2 
!     are -ve, we assume ncld{4,2} = IABS(HPA{4,2}), otherwise 
!     we take them as pressure levels and retun the layer numbers
!     within which they occur.
!
!     INPUTS:  (UPPERCASE variables)
!     P_AIRS      AIRS pressure levels (hPa)
!     T_AIRS      temperature at AIRS pressure levels (K)
!     (the above are both 1D arrays of length NLEVEL)
!     PSURF       surface pressure (hPa)
!     HPA4        layer 4 cloud position in hPa or -ve layer number
!     HPA2        layer 2 cloud position in hPa or -ve layer number
!     
!     outputs: (lowercase variables)
!     temp_layer  layer effective temperatures
!     (the above is 1D array of length NLEVEL)
!     nlayer      last layer number (i.e. surface layer)
!     ncld4       cloud4 layer number 
!     ncld2       cloud2 layer number
!
!     USES:  
!     TRP  linear interpolation function from NSWC (attached)
!
!     NOTE:  To really do this properly we need to follow the example
!            of Paul van Delst, see:
!
!        http://www.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/
!------------------------------------------------------------------------

        IMPLICIT NONE

!---------------
!  Arguments
!---------------
        INTEGER,                         INTENT( IN )  :: nlevel
        REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: P_AIRS
        REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: T_AIRS
        REAL( fp_kind ),                 INTENT( IN )  :: PSURF
        REAL( fp_kind ),                 INTENT( IN )  :: HPA4
        REAL( fp_kind ),                 INTENT( IN )  :: HPA2
        REAL( fp_kind ), DIMENSION( : ), INTENT( OUT ) :: temp_Layer
        INTEGER,                         INTENT( OUT ) :: ncld4
        INTEGER,                         INTENT( OUT ) :: ncld2

!------------------  
! Local variables 
!------------------
        REAL( fp_kind ), DIMENSION( NLEVEL_MAX )  :: rlayer
        REAL( fp_kind ) :: pbar
        LOGICAL :: first
        INTEGER :: i
        INTEGER :: NLAYER

!--------
! Data
!--------
        DATA first / .TRUE. /

!-------
! SAVE
!-------
        SAVE first, rlayer

        nlayer = nlevel - 1

        IF ( first ) THEN
          DO i = 1, NLEVEL
            rlayer( i ) = float( i )
          ENDDO
          first = .FALSE.
        END IF

! Fixed a problem here - TG 11/14/2005        
        DO i = 1, nlayer - 1
          pbar = ( p_airs( i + 1 ) + p_airs( i ) ) / TWO
          temp_layer( i ) = TRP( PBAR, NLEVEL, P_AIRS, T_AIRS )
        END DO

        pbar = ( psurf + p_airs( nlayer ) ) / TWO
        temp_layer( nlayer ) = TRP( PBAR, NLEVEL, P_AIRS, T_AIRS )

! - assign the cloud layer numbers
        IF ( HPA4 < 0 ) THEN
          ncld4 = INT( ABS( HPA4 ) )
        ELSE
          ncld4 = INT( TRP( HPA4, NLEVEL, P_AIRS, RLAYER ) )
        END IF
        IF ( HPA2 < 0 ) THEN
          ncld2 = INT( ABS( HPA2 ) )
        ELSE
          ncld2 = INT( TRP( HPA2, NLEVEL, P_AIRS, RLAYER ) )
        END IF

      END SUBROUTINE LayerTemp


      FUNCTION trp( A, &  ! Input
                    N, &  ! Input
                    X, &  ! Input
                    Y )   ! Input

        IMPLICIT NONE
       
!--------------
!  Arguments
!--------------

        REAL( fp_kind ),                 INTENT( IN ) :: A
        INTEGER,                         INTENT( IN ) :: N
        REAL( fp_kind ), DIMENSION( N ), INTENT( IN ) :: X
        REAL( fp_kind ), DIMENSION( N ), INTENT( IN ) :: Y

!------------------
! Local variables
!------------------
        INTEGER            :: NM1
        INTEGER            :: IL
        INTEGER            :: IR
        INTEGER            :: I
        REAL( fp_kind )    :: R

!------------
! Functions
!------------
        REAL    :: trp

        NM1 = N-1
        IF ( A < X( 2 ) ) GO TO 50
        IF ( A >= X( NM1 ) ) GO TO 40
        IL = 2
        IR = NM1
!
!     BISECTION SEARCH
!
   10   I = ( IL + IR ) / 2
        IF ( I == IL ) GO TO 60
        IF ( A - X( I ) ) 20, 60, 30
   20   IR = I
        GO TO 10
   30   IL = I
        GO TO 10
!
!     A.LT.X(2) .OR. A.GE.X(N-1)
!
   40   I = NM1
        GO TO 60
   50   I = 1
!
!     EVALUATION
!
   60   R = ( A - X( I ) ) / ( X( I + 1 ) - X( I ) )
        TRP = Y( I ) + R * ( Y( I + 1 ) - Y( I ) )

      END FUNCTION trp
 

       !--------------------------------------------------------
       ! PROGRAM      CLDICE SUBROUTINE
       !--------------------------------------------------------
       ! PURPOSE      Reading ice cloud reflectance, transmittance
       !              and interpolating them if necessary
       !--------------------------------------------------------
       ! VERSION:     1.0  JIANGUO NIU  2004/07/12
       !--------------------------------------------------------
       ! DISCRIPTION  CLD LOOKUP TABLE IS COMPUTED IN THE 
       !              WAVENUMBER REGION FROM 587.30 cm-1 
       !              TO 2349.50 cm-1 WITH THE INTERVAL OF
       !              0.6 cm-1, 2938 data points.
       !--------------------------------------------------------
       ! SUBROUTINES  INTPTH INTPMU
       !--------------------------------------------------------
       ! INPUT PARAMETERS:
       ! RQWV         REQUIRED WAVENUMBER TO INTERPOLATE 
       ! TAU_vis      CLD VISIBLE OPTICAL THICKNESS
       ! De           CLD PARTICLE EFFECTIVE SIZE
       ! THETA        SATELITTE ZENITH ANGLE
       
       ! OUTPUTS:
       ! R            CLD REFLECTANCE
       ! T            CLD TRANSMITTANCE
       
       ! NUM          TOTAL WAVENUMBERS IN R&T LOOKUP TABLE
       ! TAUMAX       TOTAL DIFFERENT TAUvis
       ! DMAX         TOTAL DIFFERENT De
       ! THETMAX      TOTAL DIFFERENT THETA ANGLES
       !--------------------------------------------------------

       SUBROUTINE CLDICE( RQWV,    &  ! Input
                          TAU_vis, &  ! Input
                          De,      &  ! Input
                          THETA,   &  ! Input
                          R,       &  ! Output
                          T )         ! Output

          IMPLICIT NONE

!-------------
!  Arguments
!-------------
          REAL( fp_kind ), INTENT( IN )  :: RQWV
          REAL( fp_kind ), INTENT( IN )  :: TAU_vis
          REAL( fp_kind ), INTENT( IN )  :: De
          REAL( fp_kind ), INTENT( IN )  :: THETA
          REAL( fp_kind ), INTENT( OUT ) :: R 
          REAL( fp_kind ), INTENT( OUT ) :: T 

!-----------------
! Local variables
!-----------------
          INTEGER :: indexT
          INTEGER :: indexT_
          INTEGER :: indexD
          INTEGER :: indexD_
          INTEGER :: i
          INTEGER :: index
          REAL( fp_kind ), DIMENSION( TAUMAX_ICE )       :: TAU_v
          REAL( fp_kind ), DIMENSION( DMAX_ICE )         :: D
          REAL( fp_kind ), DIMENSION( NUM_ICE )          :: nu
          LOGICAL                                        :: first
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RTT1
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RTT2
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RTD1
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RTD2
          REAL( fp_kind )                                :: RD1
          REAL( fp_kind )                                :: TD1
          REAL( fp_kind )                                :: RD2
          REAL( fp_kind )                                :: TD2
          REAL( fp_kind )                                :: RT1
          REAL( fp_kind )                                :: TT1
          REAL( fp_kind )                                :: RT2
          REAL( fp_kind )                                :: TT2
          REAL( fp_kind )                                :: Rtemp1
          REAL( fp_kind )                                :: Ttemp1
          REAL( fp_kind )                                :: Rtemp2
          REAL( fp_kind )                                :: Ttemp2
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RT
          REAL( fp_kind )                                :: R1
          REAL( fp_kind )                                :: T1
          REAL( fp_kind )                                :: R2
          REAL( fp_kind )                                :: T2
          INTEGER                                        :: index_wn
          INTEGER                                        :: index_wn1
          INTEGER                                        :: index_wn2

!----------
!  Data
!----------
          DATA D / 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., &
                   160., 170., 180. /
          DATA TAU_v / 0.01, 0.05, 0.1, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, &
                       5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0 / 
          DATA first / .TRUE. /
          SAVE

! Compute wavenumber grid used in tables first time through only
          IF ( first ) THEN
            DO i = 1, NUM_ICE
              nu( i ) = 100.0_fp_kind + 1.0_fp_kind * ( i - 1 )
            END DO
            DO i = 1, THETMAX
              THET( i ) =  10.0_fp_kind * ( i - 1 ) 
            END DO
            first = .FALSE.
          END IF
  
          IF ( TAU_vis < TAU_v( 1 ) ) THEN
            PRINT*, 'INPUT TAU SMALLER THAN LIMIT (0.01)'
            GOTO 1000
          END IF
           
          IF ( TAU_vis > TAU_v( TAUMAX_ICE ) ) THEN
            TAU_idx = TAUMAX_ICE
            TAU_idx_= TAUMAX_ICE
            GOTO 100
          END IF
          DO i = 1, TAUMAX_ICE
            IF (TAU_v( i ) == TAU_vis ) THEN
              TAU_idx = i
              TAU_idx_= i
              GOTO 100    
            ELSE IF ( TAU_v( i ) > TAU_vis ) THEN
              TAU_idx = i
              TAU_idx_= i-1
              GOTO 100
            END IF
          END DO
         
100       IF ( De > D( DMAX_ICE ) ) THEN     ! Check added by TG 3/1/06
            D_idx = DMAX_ICE                 !
            D_idx_= DMAX_ICE                 !
            GOTO 110                         !
          END IF                             !
          DO i = 1, DMAX_ICE
            IF ( D( i ) == De ) THEN
              D_idx = i
              D_idx_= i
              GOTO 110
            ELSE IF ( D( i ) > De ) THEN
              D_idx = i
              D_idx_= i-1
              GOTO 110
            END IF
          END DO
       
 110      CONTINUE
!-------------------------
! Get wavenumber indices
!-------------------------
!          index_wn = int( ( RQWV - 587.3_fp_kind ) / 0.6_fp_kind + 1.5_fp_kind )
          index_wn = int( ( RQWV - 100.0_fp_kind ) + 1.5_fp_kind )
          IF ( RQWV > nu( NUM_ICE ) ) index_wn = NUM_ICE
          index_wn1 = index_wn
          IF ( index_wn == NUM_ICE ) THEN
            index_wn2 = index_wn1
          ELSE
            index_wn2 = index_wn1 + 1
          END IF

       !--------------------------------------------------------
       !  TAUvis and De are on the grid of cld table
       !--------------------------------------------------------
          IF ( ( TAU_idx == TAU_idx_ ) .AND. ( D_idx == D_idx_ ) ) THEN
            index = DMAX_ICE * ( TAU_idx - 1 ) + D_idx
            RT = RT_ICE( :, index_wn1, index )
            CALL INTPTH( THETA, RT, R1, T1 )
            RT = RT_ICE( :, index_wn2, index )
            CALL INTPTH( THETA, RT, R2, T2 )
            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF
       !--------------------------------------------------------
       ! TAUvis on the grid of table, De not on
       !--------------------------------------------------------
          IF ( (TAU_idx == TAU_idx_ ) .AND. ( D_idx /= D_idx_ ) ) THEN
            indexD_ = DMAX_ICE * ( TAU_idx - 1 ) + D_idx_
            indexD  = DMAX_ICE * ( TAU_idx - 1 ) + D_idx
            RTD1 = RT_ICE( :, index_wn1, indexD_ )
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            RTD2 = RT_ICE( :, index_wn1, indexD )
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            R1 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            T1 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1

            RTD1 = RT_ICE( :, index_wn2, indexD_ )
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            RTD2 = RT_ICE( :, index_wn2, indexD )
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            R2 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            T2 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1

            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF
       !--------------------------------------------------------
       ! BOTH TAUvis and De not on the grid of cld table
       !--------------------------------------------------------
          IF ( ( TAU_idx /= TAU_idx_ ) .AND. ( D_idx /= D_idx_ ) ) THEN
            indexT_= DMAX_ICE * ( TAU_idx_ - 1 ) + D_idx_
            indexT = DMAX_ICE * ( TAU_idx_ - 1 ) + D_idx
            indexD_= DMAX_ICE * ( TAU_idx - 1 ) + D_idx_
            indexD = DMAX_ICE * ( TAU_idx - 1 ) + D_idx
            RTT1 = RT_ICE( :, index_wn1, indexT_ )
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            RTT2 = RT_ICE( :, index_wn1, indexT )
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            Rtemp1 = ( De - D( D_idx_ ) ) * ( RT2 - RT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RT1
            Ttemp1 = ( De - D( D_idx_ ) ) * ( TT2 - TT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TT1
            RTD1 = RT_ICE( :, index_wn1, indexD_ )
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            RTD2 = RT_ICE( :, index_wn1, indexD )
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            Rtemp2 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            Ttemp2 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1
            R1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Rtemp2 - Rtemp1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Rtemp1
            T1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Ttemp2 - Ttemp1 ) /  &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Ttemp1
            RTT1 = RT_ICE( :, index_wn2, indexT_ )
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            RTT2 = RT_ICE( :, index_wn2, indexT )
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            Rtemp1 = ( De - D( D_idx_ ) ) * ( RT2 - RT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RT1
            Ttemp1 = ( De - D( D_idx_ ) ) * ( TT2 - TT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TT1
            RTD1 = RT_ICE( :, index_wn2, indexD_ )
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            RTD2 = RT_ICE( :, index_wn2, indexD )
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            Rtemp2 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            Ttemp2 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1
            R2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Rtemp2 - Rtemp1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Rtemp1
            T2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Ttemp2 - Ttemp1 ) /  &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Ttemp1
            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R , T )

          END IF
       !--------------------------------------------------------
       ! TAUvis not on the grid, De on the grid of cld table
       !--------------------------------------------------------
          IF ( ( TAU_idx /= TAU_idx_ ) .AND. ( D_idx == D_idx_ ) ) THEN
            indexT_= DMAX_ICE * ( TAU_idx_ - 1 ) + D_idx
            indexT = DMAX_ICE * ( TAU_idx  - 1 ) + D_idx
            RTT1 = RT_ICE( :, index_wn1, indexT_ )
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            RTT2 = RT_ICE( :, index_wn1, indexT )
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            R1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * (RT2 - RT1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + RT1
            T1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( TT2 - TT1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + TT1

            RTT1 = RT_ICE( :, index_wn2, indexT_ )
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            RTT2 = RT_ICE( :, index_wn2, indexT )
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            R2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * (RT2 - RT1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + RT1
            T2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( TT2 - TT1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + TT1

            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF

1000      RETURN

       END SUBROUTINE cldice

!------------------------------------------------------------
!  PROGRAM    INTPMU SUBROUTINE
!------------------------------------------------------------
!  PURPOSE    SUBROUTINE TO MAKE R & T INTERPOLATION
!             FOR CORRESPONDING REQUIRED WAVENUMBER
!------------------------------------------------------------
!  VERSION    1.0      JIANGUO NIU      2004/08/04
!------------------------------------------------------------
!  DESCRIPTION 
!------------------------------------------------------------
!  VARIABLE   DESCRIPTION
!  nu         WAVENUMBER GRID IN CLD DATABASE
!  inu1       Start wavenumber index
!  inu2       End wavenumber index
!  NUM        TOTAL WAVENUMBERS IN R&T DATABASE
!  RQWV       REQUIRED WAVENUMBER GRID FOR INTERPOLATION
!  R1         ORIGINAL R DATA FOR WAVENUMBER INTERPOLATION
!  T1         ORIGINAL T DATA FOR WAVENUMBER INTERPOLATION
!  R2         ORIGINAL R DATA FOR WAVENUMBER INTERPOLATION
!  T2         ORIGINAL T DATA FOR WAVENUMBER INTERPOLATION
!  R          OUTPUT R ON THE REQUIRED WAVENUMBER GRID
!  T          OUTPUT T ON THE REQUIRED WAVENUMBER GRID
!------------------------------------------------------------
       SUBROUTINE INTPMU( nu,    & ! Input
                          inu1,  & ! Input
                          inu2,  & ! Input 
                          R1,    & ! Input
                          T1,    & ! Input
                          R2,    & ! Input
                          T2,    & ! Input
                          RQWV,  & ! Input
                          R,     & ! Output
                          T )      ! Output

          IMPLICIT NONE

!--------------
! Arguments
!--------------

          REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: nu
          INTEGER,                         INTENT( IN )  :: inu1
          INTEGER,                         INTENT( IN )  :: inu2
          REAL( fp_kind ),                 INTENT( IN )  :: R1
          REAL( fp_kind ),                 INTENT( IN )  :: T1
          REAL( fp_kind ),                 INTENT( IN )  :: R2
          REAL( fp_kind ),                 INTENT( IN )  :: T2
          REAL( fp_kind ),                 INTENT( IN )  :: RQWV
          REAL( fp_kind ),                 INTENT( OUT ) :: R
          REAL( fp_kind ),                 INTENT( OUT ) :: T

!------------------
! Local variables
!------------------
          
          INTEGER :: i, jT    

!
! Don't bother to interpolate if input wavennumber occurs at either endpoint of lookup table
!
          IF ( nu( inu1 ) == RQWV ) THEN
             R = R1
             T = T1
             RETURN
          END IF
          IF ( nu( inu2 ) == RQWV ) THEN
             R = R2
             T = T2
             RETURN
          END IF

          IF ( inu2 /= inu1 ) THEN
            R = R1 + ( R2 - R1 ) * ( RQWV - nu( inu1 ) ) / ( nu( inu2 ) - nu( inu1 ) )
            T = T1 + ( T2 - T1 ) * ( RQWV - nu( inu1 ) ) / ( nu( inu2 ) - nu( inu1 ) )
          ELSE
            R = R2
            T = T2
          END IF

       END SUBROUTINE intpmu

!----------------------------------------------------------------
!  PROGRAM    INTPTH  SUBROUTINE
!----------------------------------------------------------------
!  PURPOSE    VIEWING ANGLE INTERPOLATION FOR R&T OF
!             ICE & WATER CLOUD IN DATABASE
!----------------------------------------------------------------
!  VERSION    1.0 JIANGUO NIU  2004/AUGUST/02
!----------------------------------------------------------------
!  VARIABLE   DESCRIPTION
!----------------------------------------------------------------
!  THETA      INPUT, SATELITTE VIEWING ANGLE
!  RT         INPUT, R&T ARRAY
!  R          OUTPUT, REFLECTANCES AT REQUIRED THETA
!  T          OUTPUT, TREANMITTANCES AT REQUIRED THETA
!----------------------------------------------------------------
       SUBROUTINE INTPTH( THETA,   & ! Input
                          RT,      & ! Input
                          R,       & ! Output
                          T )        ! Ouput

          IMPLICIT NONE
          
!-----------
! Arguments 
!-----------           
          REAL( fp_kind ),                           INTENT( IN )  :: THETA
          REAL( fp_kind ), DIMENSION( 2 * THETMAX ), INTENT( IN )  :: RT
          REAL( fp_kind ),                           INTENT( OUT ) :: R
          REAL( fp_kind ),                           INTENT( OUT ) :: T
          
!------------------
! Local variables
!------------------
          INTEGER :: i, k
          REAL( fp_kind ) :: SLP

          DO i = 1, THETMAX
            IF ( THET( i ) == THETA ) THEN
              THET_idx  = i
              THET_idx_ = i
              GOTO 210
            ELSE IF ( THET( i ) > THETA ) THEN
              THET_idx  = i
              THET_idx_ = i - 1
              GOTO 210
            END IF
          END DO

210       IF ( THET_idx == THET_idx_ ) THEN
            R = RT( THET_idx )
            T = RT( THETMAX + THET_idx )
          ELSE

! This code does not interpolate properly - TG 5/23/2006
!            R = ( THETA - THET( THET_idx_ ) ) * ( RT( THET_idx ) &
!                      - RT( THET_idx_ ) ) / ( THET( THET_idx ) - THETA &
!                      - THET( THET_idx_ ) ) + RT( THET_idx_ )
! This code does not interpolate properly - TG 5/23/2006
!            T = ( THETA - THET( THET_idx_ ) ) * ( RT( THETMAX + THET_idx ) &
!                      - RT( THETMAX + THET_idx_ ) ) / ( THET( THET_idx ) - THETA &
!                      - THET( THET_idx_ ) ) + RT( THETMAX + THET_idx_ )
            
            SLP = ( THETA - THET( THET_idx_ ) ) / ( THET( THET_idx ) - THET( THET_idx_ ) )
            R = ( 1.0_fp_kind - SLP ) * RT( THET_idx_ ) + SLP * RT( THET_idx )
            T = ( 1.0_fp_kind - SLP ) * RT( THETMAX + THET_idx_) + SLP * RT( THETMAX + THET_idx )

!    print *, ' R = ', r
!    print *, ' rt = ', rt(thet_idx), rt(thet_idx_)
!    print *, ' T = ', t
!    print *, ' THETA = ', theta, thet(thet_idx), thet(thet_idx_)
!    print *, ' rt = ', rt(thetmax+thet_idx), rt(thetmax+thet_idx_)

          END IF

       END SUBROUTINE intpth


       !--------------------------------------------------------
       ! PROGRAM    CLDWATER SUBROUTINE
       !--------------------------------------------------------
       ! PURPOSE    READ WATER CLD REFLECTANCE, TRANSMITTANCE
       !            SND INTERPOLATING THEM IF NECESSARY
       !--------------------------------------------------------
       ! VERSION:   1.0  JIANGUO NIU  2004/07/28
       !--------------------------------------------------------
       ! DISCRIPTION  CLD LOOKUP TABLE IS COMPUTED IN THE
       !              WAVENUMBER REGION FROM 587.30 cm-1
       !              TO 2349.50 cm-1 WITH THE INTERVAL OF
       !              0.6 cm-1, 2938 DATA POINTS.
       !--------------------------------------------------------
       ! SUBROUTINES  INTPTH INTPMU
       !--------------------------------------------------------
       ! indexT : INDEX FOR READING  TAU_vis
       ! indexT_: INDEX FOR READING  TAU_vis at 'indexT-1'
       ! INPUT PARAMETERS:
       ! FNAME        R & T LOOKUP TABLE NAME
       ! RQNUM        REQUIRED WAVENUMBERS IN MAIN ROUTINE
       ! TAU_vis      CLD VISIBLE OPTICAL THICKNESS
       ! De           CLD PARTICAL EFFECTIVE SIZE
       ! THETA        SATELITTE ZENITH ANGLE
       
       ! OUTPUTS PARAMETERS:
       ! R            CLD REFLECTANCE
       ! T            CLD TRANSMITTANCE
       
       ! NUM          TOTAL WAVENUMBERS IN R&T LOOKUP TABLE
       ! TAUMAX       TOTAL DIFFERENT TAUvis
       ! DMAX         TOTAL DIFFERENT De
       ! THETMAX      TOTAL DIFFERENT THETA ANGLES
       !--------------------------------------------------------
       SUBROUTINE CLDWATER( RQWV,    & ! Input
                            TAU_vis, & ! Input
                            De,      & ! Input
                            THETA,   & ! Input
                            R,       & ! Output
                            T )        ! Output
       

          IMPLICIT NONE
          
!------------
! Arguments
!------------
          REAL( fp_kind ),    INTENT( IN ) :: RQWV
          REAL( fp_kind ),    INTENT( IN ) :: TAU_vis
          REAL( fp_kind ),    INTENT( IN ) :: De
          REAL( fp_kind ),    INTENT( IN ) :: THETA
          REAL( fp_kind ),    INTENT( OUT ) :: R
          REAL( fp_kind ),    INTENT( OUT ) :: T

!-------------------
! Local variables
!-------------------
          INTEGER :: indexT
          INTEGER :: indexT_
          INTEGER :: indexD
          INTEGER :: indexD_
          INTEGER :: i
          INTEGER :: index
          REAL( fp_kind ), DIMENSION( TAUMAX_LIQ )  :: TAU_v
          REAL( fp_kind ), DIMENSION( DMAX_LIQ )    :: D
          REAL( fp_kind ), DIMENSION( NUM_LIQ )     :: nu
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RTT1
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RTT2
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RTD1
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RTD2
          REAL( fp_kind )                                :: RD1
          REAL( fp_kind )                                :: TD1
          REAL( fp_kind )                                :: RD2
          REAL( fp_kind )                                :: TD2
          REAL( fp_kind )                                :: RT1
          REAL( fp_kind )                                :: TT1
          REAL( fp_kind )                                :: RT2
          REAL( fp_kind )                                :: TT2
          REAL( fp_kind )                                :: Rtemp1
          REAL( fp_kind )                                :: Ttemp1
          REAL( fp_kind )                                :: Rtemp2
          REAL( fp_kind )                                :: Ttemp2
          REAL( fp_kind ), DIMENSION( 2 * THETMAX )      :: RT
          REAL( fp_kind )                                :: R1
          REAL( fp_kind )                                :: T1
          REAL( fp_kind )                                :: R2
          REAL( fp_kind )                                :: T2
          INTEGER                                        :: index_wn
          INTEGER                                        :: index_wn1
          INTEGER                                        :: index_wn2
          LOGICAL first

!---------------------
! Data initialization
!---------------------
          DATA  TAU_v / 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, &
                        14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 110.0 /
          DATA  D     / 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 25.0, 37.0 , 50.0, 75.0, 100.0 /
          DATA  first / .TRUE. /

          SAVE  first, nu

          IF ( first ) THEN
            DO i = 1, NUM_LIQ
              nu( i ) = 587.3_fp_kind + 0.6_fp_kind * ( i - 1 )
            END DO 
            DO i = 1, THETMAX
              THET( i ) = 10.0_fp_kind * ( i - 1 )
            END DO 
            first = .FALSE.
          END IF

          IF ( TAU_vis < TAU_v( 1 ) ) THEN
            PRINT*,' INPUT TAU SMALLER THAN LOWER LIMIT (0.5)'
            RETURN
          END IF           
          IF ( TAU_vis > TAU_v( TAUMAX_LIQ ) ) THEN
            TAU_idx  = TAUMAX_LIQ
            TAU_idx_ = TAUMAX_LIQ
            GOTO 100
          END IF
          DO i = 1, TAUMAX_LIQ
            IF ( TAU_v( i ) == TAU_vis ) THEN
              TAU_idx  = i
              TAU_idx_ = i
              GOTO 100    
            ELSE IF ( TAU_v( i ) > TAU_vis ) THEN
              TAU_idx  = i
              TAU_idx_ = i - 1
              GOTO 100
            END IF
          END DO
       
100       IF ( De > D( DMAX_LIQ ) ) THEN     ! Check added by TG 3/1/06
            D_idx = DMAX_LIQ                 !
            D_idx_= DMAX_LIQ                 !
            GOTO 110                         !
          END IF                             !
          DO i = 1, DMAX_LIQ
            IF ( D( i ) == De ) THEN
              D_idx  = i
              D_idx_ = i
              GOTO 110
            ELSE IF ( D( i ) > De ) THEN
              D_idx  = i
              D_idx_ = i - 1
              GOTO 110
            END IF
          END DO

110       CONTINUE

!-------------------------
! Get wavenumber indices
!-------------------------
          index_wn = int( ( RQWV - 587.3_fp_kind ) / 0.6_fp_kind + 1.5_fp_kind )
          IF ( RQWV > nu( NUM_LIQ ) ) index_wn = NUM_LIQ
          index_wn1 = index_wn
          IF ( index_wn == NUM_LIQ ) THEN
            index_wn2 = index_wn1
          ELSE
            index_wn2 = index_wn1 + 1
          END IF

       !--------------------------------------------------------
       ! WATER TAU & De ARE ON THE GRID OF ITS TABLE
       !--------------------------------------------------------
          IF ( ( TAU_idx == TAU_idx_ ) .AND. ( D_idx == D_idx_ ) ) THEN
            index = DMAX_LIQ * ( TAU_idx - 1 ) + D_idx
            RT = RT_LIQ( :, index_wn1, index )
            CALL INTPTH( THETA, RT, R1, T1 )
            RT = RT_LIQ( :, index_wn2, index )
            CALL INTPTH( THETA, RT, R2, T2 )
            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF
       !--------------------------------------------------------
       ! WATER TAU ON THE GRID OF TABLE, De NOT
       !--------------------------------------------------------
          IF ( ( TAU_idx == TAU_idx_ ) .AND. ( D_idx /= D_idx_ ) ) THEN
            indexD_ = DMAX_LIQ * ( TAU_idx - 1 ) + D_idx_
            indexD  = DMAX_LIQ * ( TAU_idx - 1 ) + D_idx
            RTD1 = RT_LIQ( :, index_wn1, indexD_ )
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            RTD2 = RT_LIQ( :, index_wn1, indexD )
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            R1 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            T1 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1

            RTD1 = RT_LIQ( :, index_wn2, indexD_ )
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            RTD2 = RT_LIQ( :, index_wn2, indexD )
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            R2 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            T2 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1

            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF
       !--------------------------------------------------------
       ! BOTH TAU & De ARE NOT ON THE GRID OF TABLE
       !--------------------------------------------------------
          IF ( ( TAU_idx /= TAU_idx_ ) .AND. ( D_idx /= D_idx_ ) ) THEN
            indexT_ = DMAX_LIQ * ( TAU_idx_ - 1 ) + D_idx_
            indexT  = DMAX_LIQ * ( TAU_idx_ - 1 ) + D_idx
            indexD_ = DMAX_LIQ * ( TAU_idx - 1 ) + D_idx_
            indexD  = DMAX_LIQ * ( TAU_idx - 1 ) + D_idx
            RTT1 = RT_LIQ( :, index_wn1, indexT_ )
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            RTT2 = RT_LIQ( :, index_wn1, indexT )
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            Rtemp1 = ( De - D( D_idx_ ) ) * ( RT2 - RT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RT1
            Ttemp1 = ( De - D( D_idx_ ) ) * ( TT2 - TT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TT1
            RTD1 = RT_LIQ( :, index_wn1, indexD_ )
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            RTD2 = RT_LIQ( :, index_wn1, indexD)
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            Rtemp2 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            Ttemp2 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1
            R1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Rtemp2 - Rtemp1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Rtemp1
            T1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Ttemp2 - Ttemp1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Ttemp1

            RTT1 = RT_LIQ( :, index_wn2, indexT_ )
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            RTT2 = RT_LIQ( :, index_wn2, indexT )
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            Rtemp1 = ( De - D( D_idx_ ) ) * ( RT2 - RT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RT1
            Ttemp1 = ( De - D( D_idx_ ) ) * ( TT2 - TT1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TT1
            RTD1 = RT_LIQ( :, index_wn2, indexD_ )
            CALL INTPTH( THETA, RTD1, RD1, TD1 )
            RTD2 = RT_LIQ( :, index_wn2, indexD)
            CALL INTPTH( THETA, RTD2, RD2, TD2 )
            Rtemp2 = ( De - D( D_idx_ ) ) * ( RD2 - RD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + RD1
            Ttemp2 = ( De - D( D_idx_ ) ) * ( TD2 - TD1 ) / ( D( D_idx ) - D( D_idx_ ) ) + TD1
            R2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Rtemp2 - Rtemp1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Rtemp1
            T2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( Ttemp2 - Ttemp1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + Ttemp1

            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF
       !--------------------------------------------------------
       ! WATER De ON THE GRID, TAU NOT ON IT
       !--------------------------------------------------------
          IF ( ( TAU_idx /= TAU_idx_ ) .AND. ( D_idx == D_idx_ ) ) THEN
            indexT_= DMAX_LIQ * ( TAU_idx_ - 1 ) + D_idx
            indexT = DMAX_LIQ * ( TAU_idx - 1 ) + D_idx
            RTT1 = RT_LIQ( :, index_wn1, indexT_ )
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            RTT2 = RT_LIQ( :, index_wn1, indexT )
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            R1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( RT2 - RT1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + RT1
            T1 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( TT2 - TT1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + TT1

            RTT1 = RT_LIQ( :, index_wn2, indexT_ )
            CALL INTPTH( THETA, RTT1, RT1, TT1 )
            RTT2 = RT_LIQ( :, index_wn2, indexT )
            CALL INTPTH( THETA, RTT2, RT2, TT2 )
            R2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( RT2 - RT1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + RT1
            T2 = ( TAU_vis - TAU_v( TAU_idx_ ) ) * ( TT2 - TT1 ) / &
                        ( TAU_v( TAU_idx ) - TAU_v( TAU_idx_ ) ) + TT1

            CALL INTPMU( nu, index_wn1, index_wn2, R1, T1, R2, T2, RQWV, R, T )
          END IF

       END SUBROUTINE cldwater

  END MODULE two_layer_model

