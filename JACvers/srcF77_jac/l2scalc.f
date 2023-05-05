      SUBROUTINE L2Scalc(DOJAC,DOSUN,SECANG,SECSUN,NLAY,NCHAN,TAU4,L2S4,WGT4)

      IMPLICIT NONE

      include "incFTC.f"

c input
      LOGICAL DOSUN              ! do sun calc?
      LOGICAL DOJAC              ! do jac calc?
      REAL SECANG(MAXLAY)        ! local path angle secant
      REAL SECSUN(MAXLAY)        ! secant of effective sun local path angle
      INTEGER NLAY               ! number of layers in profile
      INTEGER NCHAN              ! number of layers in profile
      REAL TAU4(4,MAXLAY,MXCHAN) ! chan layer effective optical depth for CLR,CLD1,CLD2,CLD12 at angle theta

c output
      REAL L2S4(4,MAXLAY,MXCHAN)  ! L2S tramnsmittance
      REAL WGT4(4,MAXLAY,MXCHAN)  ! WGT functiom

c local
      INTEGER iL
      INTEGER I         ! which channel

!! V1  orig code, all of April 2023, ST was bad  
!      L2S4 = 0          ! initialize
!      DO iL = 2, NLAY
!        L2S4(:,iL,1:NCHAN) = L2S4(:,iL-1,1:NCHAN) + TAU4(:,IL-1,1:NCHAN)
!      END DO

!! V2 new code, May 3, 2023
!! recall 1 = topmost layer, LBOT = GND, ST is good, T(z) is very good
      L2S4 = 0          ! initialize
      iL = 1
      L2S4(:,iL,1:NCHAN) = TAU4(:,IL,1:NCHAN)
      DO iL = 2, NLAY
        L2S4(:,iL,1:NCHAN) = L2S4(:,iL-1,1:NCHAN) + TAU4(:,IL,1:NCHAN)
      END DO

!! V3 new code, May 4, 2023
      L2S4 = 0          ! initialize
      iL = 1
      L2S4(:,iL,1:NCHAN) = 0
      DO iL = 2, NLAY
        L2S4(:,iL,1:NCHAN) = L2S4(:,iL-1,1:NCHAN) + TAU4(:,IL-1,1:NCHAN)
      END DO
      IL = NLAY + 1
      L2S4(:,iL,1:NCHAN) = L2S4(:,iL-1,1:NCHAN) + TAU4(:,IL-1,1:NCHAN)

c************************************************************************
      L2S4(:,1:NLAY+1,1:NCHAN) = exp(-L2S4(:,1:NLAY+1,1:NCHAN))

      WGT4 = 0          ! initialize
      WGT4(:,1:NLAY,1:NCHAN) = 1-exp(-TAU4(:,1:NLAY,1:NCHAN))

      RETURN
      END
