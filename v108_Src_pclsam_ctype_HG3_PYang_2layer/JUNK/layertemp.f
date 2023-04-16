
      SUBROUTINE LayerTemp( nlevel, P_AIRS, T_AIRS, PSURF, HPA4, HPA2, 
     $                      temp_layer, ncld4, ncld2)
c--------------------------------------------------------------------
c
c     Subroutine to provide layer effective temperatures from 
c     temperatures at specified pressure levels.  HPA{4,2} are 
c     used to compute the cloud layers ncld{4,2}; if HPA4 or HPA2 
c     are -ve, we assume ncld{4,2} = IABS(HPA{4,2}), otherwise 
c     we take them as pressure levels and retun the layer numbers
c     within which they occur.
c
c     INPUTS:  (UPPERCASE variables)
c     P_AIRS      AIRS pressure levels (hPa)
c     T_AIRS      temperature at AIRS pressure levels (K)
c     (the above are both 1D arrays of length NLEVEL)
c     PSURF       surface pressure (hPa)
c     HPA4        layer 4 cloud position in hPa or -ve layer number
c     HPA2        layer 2 cloud position in hPa or -ve layer number
c     
c     outputs: (lowercase variables)
c     temp_layer  layer effective temperatures
c     (the above is 1D array of length NLEVEL)
c     nlayer      last layer number (i.e. surface layer)
c     ncld4       cloud4 layer number 
c     ncld2       cloud2 layer number
c
c     USES:  
c     TRP  linear interpolation function from NSWC (attached)
c
c     NOTE:  To really do this properly we need to follow the example
c            of Paul van Delst, see:
c
c        http://www.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/
c------------------------------------------------------------------------

        IMPLICIT NONE
        INCLUDE 'pyang2cld.param'

c---------------
c  Arguments
c---------------
        INTEGER nlevel
        REAL p_airs(*)
        REAL t_airs(*)
        REAL psurf
        REAL HPA4
        REAL HPA2
        REAL temp_Layer(*)
        INTEGER ncld4, ncld2

c------------------  
c Local variables 
c------------------
        REAL rlayer(NLEVEL_MAX)
        REAL pbar
        LOGICAL first
        INTEGER i
        INTEGER NLAYER
        REAL trp   !!!function
c--------
c Data
c--------
        DATA first / .TRUE. /

c-------
c SAVE
c-------
        SAVE first, rlayer

        nlayer = nlevel - 1

        IF ( first ) THEN
          DO i = 1, NLEVEL
            rlayer( i ) = float( i )
          ENDDO
          first = .FALSE.
        END IF

c Fixed a problem here - TG 11/14/2005        
        DO i = 1, nlayer - 1
          pbar = ( p_airs( i + 1 ) + p_airs( i ) ) / TWO
          temp_layer( i ) = TRP( PBAR, NLEVEL, P_AIRS, T_AIRS )
        END DO

        pbar = ( psurf + p_airs( nlayer ) ) / TWO
        temp_layer( nlayer ) = TRP( PBAR, NLEVEL, P_AIRS, T_AIRS )

c - assign the cloud layer numbers
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

      RETURN
      END

