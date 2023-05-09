      SUBROUTINE write_header_jac_files(LISTJ,NWANTJ,NUMPROF,NUMCHAN,FREQ,
     $      IOUNTZ,IOUNG1,IOUNG2,IOUNG3,IOUNG4,IOUNG5,IOUNG6,IOUNG9,
     $      IOUNG11,IOUNG12,IOUNG103,IOUNWGT,IOUNCLD,
     $      caJacTZ,caJACWGT,caJACG1,caJACG2,caJACG3,caJACG4,caJACG5,
     $      caJACG6,caJACG9,caJACG11,caJACG12,caJACG103,caJACCLD)

       IMPLICIT NONE
       include "incFTC.f"

c input
       INTEGER NWANTJ          ! number of wanted jacs (default 0=none)
       INTEGER  LISTJ(MAXPRO)  ! list of wanted channels
       INTEGER NUMPROF,NUMCHAN ! number of channels, profiles in rtp file
       REAL   FREQ(MXCHAN)    ! chan center frequency
       INTEGER IOUNTZ,IOUNG1,IOUNG2,IOUNG3,IOUNG4,IOUNG5,IOUNG6,
     $         IOUNG9,IOUNG11,IOUNG12,IOUNG103,IOUNWGT,IOUNCLD,iFileErr
       CHARACTER*180 caJacTZ,caJACWGT,caJACG1,caJACG2,caJACG3,caJACG4,caJACG5,
     $               caJACG6,caJACG9,caJACG11,caJACG12,caJACG103,caJACCLD

c local
       INTEGER iC, INTERSECT

       write(*,'(A)') 'opening (multiple) jacobian files ... they may end up being huge!!!!'

       IF (INTERSECT(1,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG1 = ',caJacG1
         OPEN(UNIT=IOUNG1,FILE=caJacG1,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG1) NUMPROF
         WRITE(IOUNG1) NUMCHAN
         WRITE(IOUNG1) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(2,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG2 = ',caJacG2
         OPEN(UNIT=IOUNG2,FILE=caJacG2,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG2) NUMPROF
         WRITE(IOUNG2) NUMCHAN
         WRITE(IOUNG2) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(3,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG3 = ',caJacG3
         OPEN(UNIT=IOUNG3,FILE=caJacG3,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG3) NUMPROF
         WRITE(IOUNG3) NUMCHAN
         WRITE(IOUNG3) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(4,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG4 = ',caJacG4
         OPEN(UNIT=IOUNG4,FILE=caJacG4,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG4) NUMPROF
         WRITE(IOUNG4) NUMCHAN
         WRITE(IOUNG4) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(5,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG5 = ',caJacG5
         OPEN(UNIT=IOUNG5,FILE=caJacG5,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG5) NUMPROF
         WRITE(IOUNG5) NUMCHAN
         WRITE(IOUNG5) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(6,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG6 = ',caJacG6
         OPEN(UNIT=IOUNG6,FILE=caJacG6,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG6) NUMPROF
         WRITE(IOUNG6) NUMCHAN
         WRITE(IOUNG6) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(9,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG9 = ',caJacG9
         OPEN(UNIT=IOUNG9,FILE=caJacG9,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG9) NUMPROF
         WRITE(IOUNG9) NUMCHAN
         WRITE(IOUNG9) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(11,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG11 = ',caJacG11
         OPEN(UNIT=IOUNG11,FILE=caJacG11,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG11) NUMPROF
         WRITE(IOUNG11) NUMCHAN
         WRITE(IOUNG11) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(12,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG12 = ',caJacG12
         OPEN(UNIT=IOUNG12,FILE=caJacG12,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG12) NUMPROF
         WRITE(IOUNG12) NUMCHAN
         WRITE(IOUNG12) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(103,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacG103 = ',caJacG103
         OPEN(UNIT=IOUNG103,FILE=caJacG103,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNG103) NUMPROF
         WRITE(IOUNG103) NUMCHAN
         WRITE(IOUNG103) (FREQ(iC),iC=1,NUMCHAN)
       END IF

       IF (INTERSECT(100,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacTZ = ',caJacTZ
         OPEN(UNIT=IOUNTZ,FILE=caJacTZ,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNTZ) NUMPROF
         WRITE(IOUNTZ) NUMCHAN
         WRITE(IOUNTZ) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(200,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacWGT = ',caJacWGT
         OPEN(UNIT=IOUNWGT,FILE=caJacWGT,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNWGT) NUMPROF
         WRITE(IOUNWGT) NUMCHAN
         WRITE(IOUNWGT) (FREQ(iC),iC=1,NUMCHAN)
       END IF
       IF (INTERSECT(300,LISTJ(1:NWANTJ),NWANTJ) .GT. 0) THEN
         write(*,'(A,A)') 'opening jac file caJacCLD = ',caJacCLD
         OPEN(UNIT=IOUNCLD,FILE=caJacCLD,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
         WRITE(IOUNCLD) NUMPROF
         WRITE(IOUNCLD) NUMCHAN
         WRITE(IOUNCLD) (FREQ(iC),iC=1,NUMCHAN)
       END IF

       RETURN
       END
