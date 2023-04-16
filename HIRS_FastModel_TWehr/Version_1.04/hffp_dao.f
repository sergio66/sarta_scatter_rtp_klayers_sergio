C FILE NAME: hffp_dao.f
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
C     HFFP-DAO-interface stuff
C ====================================================================

      SUBROUTINE HFFP_K(SATNUM,TOTIRF,ID_IRF,
     &               PREDID,NCOEFF,COEFF,T_TAUF,T_TAUW,T_TAUO,T_OFFS,
     &               T_SEC,
     &               SURF_T,SURF_P,SECANT,SUNSEC,SUREMI,
     &               SOLRAD,ATEMP,AWATER,AOZONE,CLCJAC,
     &               CW_M,RADTOT,TBMEAN,JACOBT,JACOBW,JACOBO,
     &               NOTECT)
      IMPLICIT NONE
      include "hffp_glob_dec.f"
      include "hffp_ref_dec.f"
      include "hffp_aux_dec.f"
      include "hffp_kernel_dec.f"
C --------------------------------------------------------------------
C
C     PURPOSE:
C       this is the HFFP interface called by DAO
C
C     INPUT/OUTPUT: see include file hffp_kernel_dec.f
C
C     HISTORY:
C       written   04/27/98 by Tobias Wehr
C --------------------------------------------------------------------

C     LOCAL PARAMETER 
      INTEGER  ILAYH       ! layer number counter

C     BEGIN     
      IJOB=1     ! IJOB is always 1     
      VERBOS=0   ! verbose parameter (1=verbose, 0=silent)
      JACTYP=2   ! Jacobians Type; 2 = T_bright Jac. (not radiance Jac.)
      CLCTBR=1   ! calculate brightness temperatures on
      OT_RTB=2   ! offset tuning factors tune brightn.temp. (not radiances)
      T_ADJ=0    ! do not adjust tuning parameters
      SUNSOL=6.7851E-05        ! sun solid angle
      DO 100 ILAYH=1,NLAYER    ! copy reference to fixed gas profile
         AFIXED(ILAYH)=RFIXED(ILAYH)
         ! avoid possible division by zero in predictor calculations
         IF (AWATER(ILAYH) .le. 1.0E-16) THEN
            AWATER(ILAYH) = 1.0E-16
         ENDIF
         IF (AOZONE(ILAYH) .le. 1.0E-16) THEN
            AOZONE(ILAYH) = 1.0E-16
         ENDIF
 100  CONTINUE

      CALL KERNEL(IJOB,SATNUM,TOTIRF,ID_IRF,
     &            PREDID,NCOEFF,COEFF,T_TAUF,T_TAUW,T_TAUO,T_OFFS,
     &            T_SEC,
     &            SURF_T,SURF_P,SECANT,SUNSEC,SUREMI,SUNSOL,
     &            SOLRAD,ATEMP,AFIXED,AWATER,AOZONE,CLCJAC,
     &            CW_M,RADTOT,TBMEAN,JACOBT,JACOBW,JACOBO,VERBOS,
     &            NOTECT,JACTYP,CLCTBR,OT_RTB,T_ADJ)

      RETURN
      END
C     END OF SUBROUTINE HFFP_K

C ====================================================================

      subroutine HFFP_Ipres(presrt)
C     This subroutine copies PRES (from hffp_ref_dec.f) to
C     presrt
      implicit none
      include "hffp_glob_dec.f"
      include "hffp_ref_dec.f"
      real,dimension(:),intent(out)::presrt
      integer                      ::ilay

C     BEGIN
      if (size(presrt) .ne. nlayer) then
         write (*,*) 'FATAL ERROR in Subroutine HFFP_Ipres'
         write (*,*) '  size(presrt) .ne. nlayer'
         write (*,*) '  program aborted'
         stop
      endif
      do ilay=1,NLAYER
         ! copy in inverse order and convert atm to hPa
         presrt(ilay)=PRES(NLAYER-ilay+1)*1013.25
      enddo
      return
      end
C     END OF SUBROUTINE HFFP_IPRES

C ====================================================================

      subroutine HFFP_int_h2o(h2o,h2o_p,t_air,h2o_nr)
C --------------------------------------------------------------------
C
C     PURPOSE
C       convert a specific humidity point profile
C       to a water amount layer profile
C
C     PHYSICS
C       Calculating the amount:
C       amount = integral (vmr p)/R T  dz
C         substituting dz = dp/(-g rho)        (ideal gas law)
C       amount = - integral (vmr p)/(R T g rho) dp
C         assuming vmr,R,T,g,rho are constant within one layer:
C       amount = (p_1^2-p_2^2) vmr / 2 g rho R T
C         where  vmr = volume mixing ratio [dimensionless]
C                      vmr = partial pressure (H2O) / total pressure
C                      note that the total pressure is the pressure
C                      of dry air plus water (=pressure of wet air)
C                  g = gravity acceleration at layer [m/s^2]
C                rho = air density [kg/m^3]
C                  R = gas constant [J/(mol*K)]
C                  T = layer temperature [K]
C             amount = layer h2o amount [moles/m^2] which will
C                      be converted to [1000moles/cm^2] afterwards
C
C       The vmr can be calculated from the specific humidity H by
C            vmr = 1/(1+X),  X= ((1-H)/H) * (m_w/m_a)
C       where m_w is the molar mass of water, m_a is the molar mass of dry air
C
C       Note: unit conversions: 
C       Pa = N/m^2, J = N*m, N = kg*m/s^2
C
C     INPUT
C       h2o      specific humidity
C                units:  mass_water [kg] / mass_air [kg]
C                        mass_air = mass_dryair + mass_water
C       h2o_p    pressure values for h2o
C       h2o_nr   dimension of h2o and h2o_p
C
C     OUTPUT
C       h2o      water amounts in units of 1000*moles/cm^2
C
C     HISTORY
C       written 05/13/98 by T.Wehr
C
C --------------------------------------------------------------------
      implicit none
      include "hffp_glob_dec.f"
      include "hffp_ref_dec.f"
C     INPUT VARIABLES
      real,dimension(:),intent(in)     :: h2o_p  ! pressures
      real,dimension(:),intent(in)     :: t_air  ! temperatures
      integer,          intent(in)     :: h2o_nr ! number of layers
C     OUTPUT VARIABLES
      real,dimension(:),intent(inout)  :: h2o    ! h2o 
C     LOCAL VARIABLES
      integer                          :: ilay
      real,dimension(NLAYER)           :: glayer ! av. earth acc. at layer
      real                             :: vmr,rho,rhop_0,amount
      real                             :: m_air,m_h2o
      real                             :: R_gas,N_A

      data rhop_0    /2.6867E+19/   ! [particles/cm^3] density air
      data m_air     /28.968E-3/    ! [kg/mol]         molar mass air
      data m_h2o     /18.0153E-3/   ! [kg/mol]         molar mass water
      data R_gas     /8.31441/      ! [J/(mol*K)]      gas constant
      data N_A       /6.022045E23/  ! [1/mol]          Avogadro number

      ! These glayer data are calculated for the HFFP reference
      ! profile, assuming  glayer = g0 * (r^2/R^2)
      !    where r = av. Earth radius = 6370 km
      !          R = r + layer altitude
      !          g = 9.81 m/s^2
      !
      ! glayer: index 1 = highest pressure, index 100 = lowest pressure
      ! (that means inverse order compared to h2o indecees !)
      data glayer
     & / 9.811768 , 9.811058 , 9.810342 , 9.809620 , 9.808893 , 
     &   9.808160 , 9.807421 , 9.806676 , 9.805925 , 9.805168 , 
     &   9.804406 , 9.803638 , 9.802863 , 9.802081 , 9.801294 , 
     &   9.800499 , 9.799698 , 9.798892 , 9.798079 , 9.797258 , 
     &   9.796431 , 9.795598 , 9.794758 , 9.793910 , 9.793055 , 
     &   9.792194 , 9.791325 , 9.790448 , 9.789563 , 9.788672 , 
     &   9.787772 , 9.786867 , 9.785952 , 9.785028 , 9.784099 , 
     &   9.783160 , 9.782212 , 9.781257 , 9.780293 , 9.779319 , 
     &   9.778337 , 9.777345 , 9.776345 , 9.775331 , 9.774297 , 
     &   9.773243 , 9.772170 , 9.771073 , 9.769954 , 9.768813 , 
     &   9.767645 , 9.766453 , 9.765236 , 9.763989 , 9.762714 , 
     &   9.761408 , 9.760073 , 9.758705 , 9.757304 , 9.755867 , 
     &   9.754393 , 9.752880 , 9.751326 , 9.749730 , 9.748087 , 
     &   9.746395 , 9.744649 , 9.742846 , 9.740982 , 9.739051 , 
     &   9.737052 , 9.734979 , 9.732828 , 9.730592 , 9.728264 , 
     &   9.725841 , 9.723310 , 9.720665 , 9.717898 , 9.714992 , 
     &   9.711935 , 9.708706 , 9.705272 , 9.701607 , 9.697680 , 
     &   9.693454 , 9.688890 , 9.683937 , 9.678532 , 9.672597 , 
     &   9.666048 , 9.658838 , 9.650979 , 9.642440 , 9.633138 , 
     &   9.622940 , 9.611612 , 9.598797 , 9.583865 , 9.565517 /

C     BEGIN

      ! subroutine requires h2o_nr .eq. NLAYER
      if (h2o_nr .ne. NLAYER) then
         write (*,*) 'FATAL ERROR in HFFP_int_h2o'
         write (*,*) 'h2o_nr .ne. NLAYER'
         write (*,*) 'program aborted'
         stop
      endif
      if (size(t_air) .ne. NLAYER) then
         write (*,*) 'FATAL ERROR in HFFP_int_h2o'
         write (*,*) 'size(t_air) .ne. NLAYER'
         write (*,*) 'program aborted'
         stop
      endif

      do ilay=1,nlayer
         ! calculate volume mixing ratio
         vmr=1.0/(1.0+ ((1.0-h2o(ilay))/h2o(ilay))*m_h2o/m_air)

         !vmr= h2o(ilay)*m_air/m_h2o ! Wrong
         ! calculate air density in [kg/m^3] as function of p and T
         ! and convert from [particles/cm^3] to [kg/m^3]
         rho=rhop_0 * (m_air/N_A) * 1.0E6
     &              * (273.15/t_air(ilay)) * (h2o_p(ilay)/1013.25)
         ! integrate (incl units conversion to [1000 moles / cm^2] )
         ! conversion factors:
         !
         !   PARAMETER    UNIT CONVERSION     CONVERSION FACTOR
         !   vmr          remains 1                   1
         !   presbd       hPa -> Pa                 100**2
         !   glayer       remains m/s^2               1
         !   rho          remains kg/m^3              1
         !   R_gas        remains J/(mol*K)           1
         !   t_air        remains K                   1
         !   
         !   after considering these conversions the unit of 
         !   the result would be moles/m^2
         !
         !   converting to 1000 moles/m^2           1/1000
         !   converting to 1000 moles/cm^2          1/100**2
         !   -------------------------------------------------
	 !   total                                  1/1000
         !
         amount= vmr
     &     *( presbd(nlayer-ilay+1)*presbd(nlayer-ilay+1) 
     &       -presbd(nlayer-ilay+2)*presbd(nlayer-ilay+2) )
     &     /(2000.0*glayer(nlayer-ilay+1)*rho*R_gas*t_air(ilay))
         ! copy amount to output variable
         h2o(ilay)=amount
         
         !write (*,*) amount,h2o(ilay),h2o_p(ilay)
      enddo

      RETURN
      END
