      subroutine readtram(iounit, n, freq, ntrlev, ncoef, iret)
      implicit none

#include "amsurta.com"

C***********************************************************************
C*
C*  Program Name:      readtr
C*
C*  Description:       read rapid algorithm transmittance coefficients
C*
C*  Input Variables:
C*    Name       Type      Units       Description
C*    ----       ----      -----       -----------
C*   iounit     integer             I/O unit for the coefficient data set
C*  Output Variables:
C*    Name       Type      Units       Description
C*    ----       ----      -----       -----------
C*    freq       real    wave number frequency of channel
C*                       gHz for AMSU
C*
C*  Routines Called:
C*
C*  Files Used:
C*
C*  Accuracy:
C*
C*  Error Handling:
C*     The return code iret is non-zero if there is an error in input.
C*      0 for normal return
C*      1 for bad data
C*     -1 for end of file
C*     -2 for read error
C*     -3 if file elements not the same size as opentr/header of file
C* 
C*  References:
C*
C*  History:
C*      Version   Date          Short Descriptions
C*      VER01.00 03/10/91    NEW
C*
C*  Comments:
C*    1.   Data set format changes depending on the flag unformatted.
C*          If unformatted=.false., the subroutine reads a formatted data set.
C*          If unformatted=.true., the subroutine reads an unformatted data
C*           set.   This makes I/O faser and saves disk space.   But the
C*           unformatted data set becomes machine dependent.
C*
C***********************************************************************
C Input variables
c ---------------
      integer*4    iounit
      integer*4    n
      integer*4    ntrlev   ! # of levels
      integer*4    ncoef    ! # of coefficients
C***********************************************************************
C Output variables
c ----------------
      integer*4    iret
      real*4       freq
C
C Local variables
c ---------------

      integer*4  i, L, ichan, nlev, ncof

C***********************************************************************
C***********************************************************************

      maxchan = n
      if(maxchan.gt.MAXRTACHAN) then
        print 3000, n, MAXRTACHAN
        goto 720
      endif

C     Read data
C     ---------

      read ( iounit, 900, end=700, err=710 )
     $     freq, ichan, nlev, ncof
c     print 300, n, freq, ichan, nlev, ncof

      if(nlev.ne.ntrlev) then
        print 3010, n, freq, nlev, ntrlev
        goto 720
      endif
      if(ncof.ne.ncoef) then
        print 3020, n, freq, ncof, ncoef
        goto 720
      endif

      do L = 1, nlev
        read(iounit,910, end=700, err=710 )
     $       ( trcoef(i,L,n), i = 1, ncof )
c       print 310, ( trcoef(i,L,n), i = 1, ncof )
      enddo
      iret  = 0
      return

C     end of input reached
C     --------------------
 700  iret  = -1
      return
C
C     I/O error
C     ---------
 710  iret  = -2
      return

C     Weird stuff in file sizes
c     -------------------------

 720  iret  = -3
      l = 7
      call softexit('READTRAM',l)

      return
 900  format( f10.3, 3i5)
 910  format( 7e11.4 )
 3000 format('readtram: *** TOO MANY CHANNELS in',
     1 ' AMSU or MHS RTA file ',i7,' > ',i5,' ***')
 3010 format('readtram: i=',i3,' f=',f8.3,'  # levels is wrong in'
     1 ' AMSU or MHS RTA file ',i7,' <> ',i5,' ***')
 3020 format('readtram: i=',i3,' f=',f8.3,'  # coeficients is wrong in'
     1 ' AMSU or MHS RTA file ',i7,' <> ',i5,' ***')
 300  format('n=',i3,' f=',f10.5,3i4)
 310  format(7e10.3)
      end


