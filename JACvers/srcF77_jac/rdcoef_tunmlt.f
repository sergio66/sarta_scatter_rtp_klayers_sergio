      SUBROUTINE RDCOEF_TUNMLT(IOUN, NCHAN, INDCHN, SETCHN, 
     $  NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7,
     $ CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
     $  COEF1,  COEF2,  COEF3,  COEF4,  COEF5,  COEF6,  COEF7,
     $   FREQ, LABOVE,  COEFF, INDCO2, COFCO2, INDSO2, COFSO2,
     $ INDHNO, COFHNO, INDN2O, COFN2O, INDNH3, COFNH3, 
     $ INDHDO, COFHDO,
     $ INDH2O,  WAZOP, WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN,
     $ QUICKCLIST1, QUICKCLIST2, QUICKCLIST3, QUICKCLIST4, QUICKCLIST5, 
     $ QUICKCLIST6, QUICKCLIST7,
     $ NWANTC, LISTC, LSTCHN, RINDCHN)

       IMPLICIT NONE

       include 'incFTC.f'

       INTEGER NWANTC         ! number of wanted channels (-1=all)
       INTEGER  LISTC(MAXPRO) ! list of wanted channels
       INTEGER LSTCHN(MXCHAN)  ! list of selected channels

       INTEGER   IOUN         ! I/O unit number
       INTEGER  NCHAN          ! # of selected channels
       REAL     FCHAN(MXCHAN)  ! chan center frequency
       INTEGER INDCHN(MXCHAN)  ! array indices for all channels
       INTEGER RINDCHN(MXCHAN) ! list of locations of chans eg chID1291 = 1231 cm-1 but this is location 1520 in L1C
       REAL   FREQ(MXCHAN)    ! chan center frequency

       INTEGER LABOVE(MXCHAN) ! chan downwelling thermal layer above

C      for RDCOEF             ! Info for selected channels only
       INTEGER SETCHN(MXCHAN) ! set # for each channel
       INTEGER  NCHN1         ! # of set1 channels
       INTEGER  NCHN2         ! # of set2 channels
       INTEGER  NCHN3         ! # of set3 channels
       INTEGER  NCHN4         ! # of set4 channels
       INTEGER  NCHN5         ! # of set5 channels
       INTEGER  NCHN6         ! # of set6 channels
       INTEGER  NCHN7         ! # of set7 channels

       INTEGER CLIST1(MXCHN1) ! list of set1 channels
       INTEGER CLIST2(MXCHN2) ! list of set2 channels
       INTEGER CLIST3(MXCHN3) ! list of set3 channels
       INTEGER CLIST4(MXCHN4) ! list of set4 channels
       INTEGER CLIST5(MXCHN5) ! list of set5 channels
       INTEGER CLIST6(MXCHN6) ! list of set6 channels
       INTEGER CLIST7(MXCHN7) ! list of set7 channels

       INTEGER QUICKCLIST1(MXCHAN) ! list of set1 channels
       INTEGER QUICKCLIST2(MXCHAN) ! list of set2 channels
       INTEGER QUICKCLIST3(MXCHAN) ! list of set3 channels
       INTEGER QUICKCLIST4(MXCHAN) ! list of set4 channels
       INTEGER QUICKCLIST5(MXCHAN) ! list of set5 channels
       INTEGER QUICKCLIST6(MXCHAN) ! list of set6 channels
       INTEGER QUICKCLIST7(MXCHAN) ! list of set7 channels

       REAL  COEF1(N1COEF,MAXLAY,MXCHN1) ! coefs for set1 chans
       REAL  COEF2(N2COEF,MAXLAY,MXCHN2) ! coefs for set2 chans
       REAL  COEF3(N3COEF,MAXLAY,MXCHN3) ! coefs for set3 chans
       REAL  COEF4(N4COEF,MAXLAY,MXCHN4) ! coefs for set4 chans
       REAL  COEF5(N5COEF,MAXLAY,MXCHN5) ! coefs for set5 chans
       REAL  COEF6(N6COEF,MAXLAY,MXCHN6) ! coefs for set6 chans
       REAL  COEF7(N7COEF,MAXLAY,MXCHN7) ! coefs for set7 chans
       REAL  COEFF(NFCOEF,MXCHAN)        ! coefs for chan "F" factor
       INTEGER INDCO2(MXCHAN)            ! chan indices for CO2 pert
       REAL COFCO2(  NCO2,MAXLAY,MXCHNC) ! coefs for CO2 pert
       INTEGER INDSO2(MXCHAN)            ! chan indices for SO2 pert
       REAL COFSO2(  NSO2,MAXLAY,MXCHNS) ! coefs for SO2 pert
       INTEGER INDHDO(MXCHAN)            ! chan indices for HDO pert
       REAL COFHDO(  NHDO,MAXLAY,MXCHND) ! coefs for HDO pert
       INTEGER INDHNO(MXCHAN)            ! chan indices for HNO3 pert
       REAL COFHNO( NHNO3,MAXLAY,MXCHNH) ! coefs for HNO3 pert
       INTEGER INDN2O(MXCHAN)            ! chan indices for N2O pert
       REAL COFN2O(  NN2O,MAXLAY,MXCHNN) ! coefs for N2O pert
       INTEGER INDNH3(MXCHAN)            ! chan indices for NH3 pert
       REAL COFNH3(  NNH3,MAXLAY,MXCHNA) ! coefs for NH3 pert
       INTEGER INDH2O(MXCHAN)            ! chan indices for OPTRAN H2O
       REAL   WAZOP(MXOWLY)              ! OPTRAN water l-to-s amounts
       REAL  WAVGOP(NOWAVG,MXOWLY)       ! OPTRAN raw predictor averages
       REAL COFH2O(  NH2O,MXOWLY,MXCHNW) ! coefs for OPTRAN H2O
       REAL     FX(MAXLAY)               ! fixed gases adjustment

       INTEGER NCHNTE                    ! number of non-LTE channels
       INTEGER CLISTN(MXCNTE)            ! non-LTE channel list
       REAL  COEFN(NNCOEF,MXCNTE)        ! non-LTE coefficients

c      local
       INTEGER I,III
       INTEGER intersect

C      ------------------------
C      Read the coef data files
C      ------------------------
       CALL RDCOEF( IOUN, NCHAN, INDCHN, SETCHN,
     $  NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7,
     $ CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
     $  COEF1,  COEF2,  COEF3,  COEF4,  COEF5,  COEF6,  COEF7,
     $   FREQ, LABOVE,  COEFF, INDCO2, COFCO2, INDSO2, COFSO2,
     $ INDHNO, COFHNO, INDN2O, COFN2O, INDNH3, COFNH3, 
     $ INDHDO, COFHDO,
     $ INDH2O,  WAZOP, WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN)
C
C      Get and apply multipler tuning to coefficients {note: ignores HNO3}
       CALL TUNMLT( IOUN, NCHAN, INDCHN, SETCHN,
     $  NCHN1,  NCHN2,  NCHN3,  NCHN4,  NCHN5,  NCHN6,  NCHN7,
     $ CLIST1, CLIST2, CLIST3, CLIST4, CLIST5, CLIST6, CLIST7,
     $  COEF1,  COEF2,  COEF3,  COEF4,  COEF5,  COEF6,  COEF7,
     $   FREQ, LABOVE,  COEFF, INDCO2, COFCO2, INDSO2, COFSO2,
     $ INDHNO, COFHNO, INDN2O, COFN2O,
     $ INDH2O,  WAZOP, WAVGOP, COFH2O, FX, NCHNTE, CLISTN, COEFN )

       IF (NWANTC .GT. 0) THEN
         write(*,'(A)') '     ---------------------------------------------------------------------------------------------'
         write(*,'(A,I5,A)') 'after opnrtp, NCHAN = ',NCHAN,' .... here are the channels after opnrtp .... '
         write(*,'(A)') '              I          LSTCHN(I)      RINDCHN(I)    INDCHN(LSTCHN(I))  BREAKOUT     FCHAN(I) '
         write(*,'(A)') '                                                                     SETCHN(LSTCHN(I)          '
         write(*,'(A)') '     ---------------------------------------------------------------------------------------------'
         DO I = 1,NCHAN
           write(*,'(5(I15),F20.7)') I,LSTCHN(I),RINDCHN(I),INDCHN(LSTCHN(I)),SETCHN(LSTCHN(I)),FCHAN(I)
         END DO
         write(*,'(A)') '     ---------------------------------------------------------------------------------------------'
       END IF

      !!!! suppose this is AIRS and there are 2834 chans
      !!! see incFTC.f : so if you want all 2834 chans computed, then NCHN1 == MXCHN1 = 1461 for FWO, NCHN2 = MXCHN2 = 325 for FOW etc etc etc
      !!!      ie 1461   325     396     85    210    217    140     which sums to 2834
      !!!
      !!! but eg in my retrieval code CRODGERS_FAST_CLOUD, I have chose LW chnas only about 420 between 15 um, window, O3, WV
      !!! ie h.nchan = 426, h.ichan = >> h.ichan(1:15)' = 25 52 62 69 70 71 73 75 77 78 79 80 82 83 84 ...
      !!! then NCHNX will say of the ODs in setX, how many of these IDs should be computed!!!!!
      !!! but in the example from retrievals, this ends up being 277 56 69 67 7 0 0   which sums to 476     since no SW channels (set 6,7 = none chose) 
!      print *,NCHN1,NCHN2,NCHN3,NCHN4,NCHN5,NCHN6,NCHN7
      !!! and NCHAN witll b the sum of above == 476 = h.ichan
!      print *,NCHAN

      !!!! so there will be ZEROS where h.ichan does not exist, and the index where it does so for example this will look like
      !!!!                       0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           1           0           0           0           0     index 25
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           0           0
      !!!!           0           0           0           0           2           0     index 52
      !!!!           0           0           0           0           0           0
      !!!!           0           0           3           0           0           0     index 62
      !!!!           0           0           0           4           5           6     index 69 70 71
      !!!!           0           7           0           8           0           9     index 73 75 77
      !!!!          10          11          12           0          13          14     index 78 79 80 82 83 ... 

      QUICKCLIST1(INDCHN(CLIST1(1:NCHN1))) = (/(I,I=1,NCHN1)/)
      QUICKCLIST2(INDCHN(CLIST2(1:NCHN2))) = (/(I,I=1,NCHN2)/)
      QUICKCLIST3(INDCHN(CLIST3(1:NCHN3))) = (/(I,I=1,NCHN3)/)
      QUICKCLIST4(INDCHN(CLIST4(1:NCHN4))) = (/(I,I=1,NCHN4)/)
      QUICKCLIST5(INDCHN(CLIST5(1:NCHN5))) = (/(I,I=1,NCHN5)/)
      QUICKCLIST6(INDCHN(CLIST6(1:NCHN6))) = (/(I,I=1,NCHN6)/)
      QUICKCLIST7(INDCHN(CLIST7(1:NCHN7))) = (/(I,I=1,NCHN7)/)

       if (DEBUG) then
         print *,'NCHN4 = ',NCHN4
         print *,'SETCHN(CLIST4(1:NCHN4)) = ',SETCHN(CLIST4(1:NCHN4))
         print *,'CLIST4(1:NCHN4) = ',CLIST4(1:NCHN4)
         print *,'INDCHN(CLIST4(1:NCHN4)) = ',INDCHN(CLIST4(1:NCHN4))
         print *,'QLIST4(INDCHN(CLIST4(1:NCHN4))) = ',QUICKCLIST4(INDCHN(CLIST4(1:NCHN4)))
         do I = 1,NCHN4
           print *,I,CLIST4(I),INDCHN(CLIST4(I)),QUICKCLIST4(INDCHN(CLIST4(I)))
         end do
       end if

      if (DEBUG) then
        print *,NCHN1,NCHN2,NCHN3,NCHN4,NCHN5,NCHN6,NCHN7
        print *,NCHAN
        print *,'INDCHN(1:NCHAN) = ',INDCHN(1:NCHAN)
        DO I = 1,NCHAN
           !! print i,h.ichan(i),h.vchan(i)
           print *,I,LSTCHN(I),SETCHN(I),FREQ(I)
         END DO
         print *,'CLIST1 = ',CLIST1
         print *,'CLIST2 = ',CLIST2
!        print *,'CLIST3 = ',CLIST3
!        print *,'CLIST4 = ',CLIST4
!        print *,'CLIST5 = ',CLIST5
!        print *,'CLIST6 = ',CLIST6
!        print *,'CLIST7 = ',CLIST7
       END IF

       if (DEBUG) then
         print *,'CLIST2 = ',CLIST2(1:NCHN2)
         print *,'INDCHN(CLIST2) = ',INDCHN(CLIST2(1:NCHN2))
         DO I = 1,NCHAN
           III = intersect(I,INDCHN(CLIST2(1:NCHN2)), NCHN2)
           IF (III .GT. 0) print *,I,III,INDCHN(CLIST2(III))
         END DO
         STOP
       end if

      !!! of the above INDCHN, which match up with set OD1
      !!! CLIST1 = 25 52 62 69 70 71 73 75 77 78 79 80 82 83 84 .... 1785
      !!! CLIST2 = 149 150 196 203 206 207 217 220 221 223 224 ... 1221 1236

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    SO print *,INDCHN(CLIST1) should be INDCHN([25 52 62 69 70 71 73 75 77 ...]) = 1 2 3 4 5 6 ....
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       print *,INDCHN(CLIST2(1:NCHN2))
!         63  64  110 117 120 121
!         131 134 135 137 138 139
!         142 143 144 145 146 147
!         148 149 150 164 166 167
!         168 169 170 172 179 180
!         242 243 244 245 246 247
!         248 249 250 251 252 253
!         254 255 256 257 259 260
!         261 262 263 264 265 266
!         267 268            

       RETURN
       END
