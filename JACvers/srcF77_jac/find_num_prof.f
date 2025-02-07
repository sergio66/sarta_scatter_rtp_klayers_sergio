      SUBROUTINE find_num_prof(FIN,NUMCHAN,NUMPROF)

       IMPLICIT NONE

       include 'incFTC.f'
       include 'rtpdefs.f'

       CHARACTER*90 FIN
       INTEGER NUMPROF,NUMCHAN

       INTEGER rtpopen       ! function rtpopen
       INTEGER STATUS        ! status of RTP file open
       CHARACTER*1 MODE      ! mode for rtpopen: "c"=create, "r"=read
C      Structures (see "rtpdefs.f")
       RECORD /RTPHEAD/ HEAD            ! header data
       RECORD /RTPATTR/ HATT(MAXNATTR)  ! header attributes
       RECORD /RTPPROF/ PROF
       RECORD /RTPATTR/ PATT(MAXNATTR)  ! profile attributes

       INTEGER  IOPCI  ! I/O unit ("profile channel") for input file
       INTEGER ISTAT
       INTEGER IPROF
       INTEGER rtpread  ! for calling read rtp interface routine
       INTEGER rtpclose ! for calling close rtp interface routine

       MODE='r'
       STATUS=rtpopen(FIN, MODE, HEAD, HATT, PATT, IOPCI)

       IPROF = 0  ! initialize profile counter

 8888  IPROF = IPROF + 1
       ISTAT=rtpread(IOPCI, PROF)
       IF (ISTAT .EQ. -1) THEN
         GOTO 9999  ! reached End Of File
       ELSE
         GOTO 8888
       END IF
       
 9999  ISTAT=rtpclose(IOPCI)       

       NUMPROF = IPROF - 1
       NUMCHAN = HEAD%NCHAN

       RETURN
       END
