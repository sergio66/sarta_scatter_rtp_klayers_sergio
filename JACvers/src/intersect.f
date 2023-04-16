! this subroutine looks for intersections of two arrays
! assumes they are sorted
! copied from /home/sergio/KCARTA/SRCv1.22_f90/n_rtp.f90

      FUNCTION intersect(iX,iaX,iN)
! this sees if integer iX lies in array iaX of length iN

      IMPLICIT NONE

      integer iaX(*)
      integer iX,iN,intersect

      integer iJ,iK

      iJ = -1
      iK = 1

c      print *, iaX
c      print *,'check for iX = ',iX,' in array of length iN = ',iN

      do while ( (iK <= iN) .and. (iJ < 0))
c        print *,iK,iJ,iX,iN
        if (iX .EQ. iaX(iK)) THEN
           iJ = iK
        else
           iK = iK + 1
        end if
      end do

      intersect = iJ

      END FUNCTION intersect

