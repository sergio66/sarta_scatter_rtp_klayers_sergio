USE INCFTC

TYPE CLOUD
	! clear flag/code
    INTEGER :: clrflag              ! clear flag/code
                
    ! cloud1 data
    INTEGER :: ctype                ! cloud type code
    REAL :: cfrac                ! cloud fraction 
    REAL :: cemis(MXEMIS)       ! cloud top emissivity
    REAL :: crho(MXEMIS)        ! cloud top reflectivity
    REAL :: cprtop               ! cloud top pressure
    REAL :: cprbot               ! cloud bottom pressure
    REAL :: cngwat               ! cloud non-gas water
    REAL :: cpsize               ! cloud particle size
    REAL :: cstemp               ! cloud surface temperature

    ! cloud2 data
    INTEGER :: ctype2               ! cloud2 type code
    REAL :: cfrac2               ! cloud2 fraction 
    REAL :: cemis2(MXEMIS)      ! cloud2 top emissivity
    REAL :: crho2(MXEMIS)       ! cloud2 top reflectivity
    REAL :: cprtop2              ! cloud2 top pressure
    REAL :: cprbot2              ! cloud2 bottom pressure
    REAL :: cngwat2              ! cloud2 non-gas water
    REAL :: cpsize2              ! cloud2 particle size
    REAL :: cstemp2              ! cloud2 surface temperature
    REAL :: cfrac12              ! cloud1+2 fraction

END TYPE CLOUD

END