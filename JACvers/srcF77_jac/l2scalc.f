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
      REAL TAU4(4,MAXLAY,MXCHAN) ! chan layer effective optical depth for CLR,CLD1,CLD2,CLD12

c output
      REAL L2S4(4,MAXLAY,MXCHAN),WGT4(4,MAXLAY,MXCHAN)

c local
      INTEGER iL
      INTEGER I               ! which channel
      
      L2S4 = 0          ! initialize
      WGT4 = 0          ! initialize

      DO iL = 2, NLAY
        L2S4(:,iL,1:NCHAN) = L2S4(:,iL-1,1:NCHAN) + TAU4(:,IL-1,1:NCHAN)
      END DO
      L2S4(:,1:NLAY,1:NCHAN) = exp(-L2S4(:,1:NLAY,1:NCHAN))

      WGT4(:,1:NLAY,1:NCHAN) = 1-exp(-TAU4(:,1:NLAY,1:NCHAN))

      RETURN
      END
