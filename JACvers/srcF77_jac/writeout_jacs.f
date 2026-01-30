c testing
c       FCLEAR = 0.0; CFRA1X = 1.0; CFRA2X = 0.0; CFRA12 = 0.0
c       FCLEAR = 0.0; CFRA1X = 0.0; CFRA2X = 1.0; CFRA12 = 0.0
c       FCLEAR = 0.0; CFRA1X = 0.0; CFRA2X = 0.0; CFRA12 = 1.0
c       FCLEAR = 1.0; CFRA1X = 0.0; CFRA2X = 0.0; CFRA12 = 0.0

!!!!!! V1 orig code, all of April 2023
!        L2S4above = 0  !!! L2S4above(4,MAXLAY,MXCHAN)
!        L2S4above(:,1:NLAY,1:NCHAN) =  L2S4(:,1:NLAY,1:NCHAN)

!!!!!! V2 new code, May 3, 2023, T is very good
!        L2S4above = 0  !!! L2S4above(4,MAXLAY,MXCHAN)
!        L2S4above(:,2:NLAY,1:NCHAN) =  L2S4(:,1:NLAY-1,1:NCHAN)
!        L2S4above(:,NLAY+1,1:NCHAN) =  L2S4(:,NLAY,    1:NCHAN)

!!!!!! V3 orig code, May2023
        L2S4above = 0  !!! L2S4above(4,MAXLAY,MXCHAN)
        L2S4above(:,1:NLAY+1,1:NCHAN) =  L2S4(:,1:NLAY+1,1:NCHAN)

!!! all of SARTA ODs are in terms of angles ie TAU_SARTA = TAU_NADIR/cos(theta) = mu TAU_NADIR
!!!                                         so d(TAU_SARTA)/dX = mu d(TAU_NADIR)/dX
!!!                                         so d(TAU_NADIR)/dX = 1/mu  d(TAU_SARTA)/dX = sec(theta) d(TAU_SARTA)/dX
!      ONESSS(1:NCHAN) = 1.0
!      ONESSS = 1
!      RAMU(1:NLAY,1:NCHAN) = matmul( RESHAPE(SECANG(1:NLAY), (/NLAY,1/)), ONESSS) 
!      RAMU(1:NLAY,1:NCHAN) = 1    !! we uncompressed tau(x) = 1/cos(x) tau(0) = 1/mu tau(0)
!                                  !! so dtau(x)/dX = 1/mu dtau(0)/dX which is what we want

!!!!      RTHERM4_SOLAR4 = 0    !!! DEBUG
c%%%%%%%%%%%%%%%%%%%%%%%%%
c note : we know "blmult" in sarta_pclsam.f and will use it to adjust jacobians for lowest layer
      include "writeout_gas1_12jacs.f"
      include "writeout_TZjacs_WGTFCN.f"
      include "writeout_cldjacs.f"
c%%%%%%%%%%%%%%%%%%%%%%%%%
