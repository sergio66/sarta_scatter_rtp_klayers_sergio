       SUBROUTINE  RDSUN_SET_NFAKE(IOUN, INDCHN, HSUN, 
     $         NCHN1,NCHN2,NCHN3,CLIST1,CLIST2,CLIST3,
     $         NFAKE,INDFAK,QUICKINDFAK)

       IMPLICIT NONE
       include 'incFTC.f'

cinout
       INTEGER   IOUN         ! I/O unit number
       INTEGER INDCHN(MXCHAN)  ! array indices for all channels

       INTEGER CLIST1(MXCHN1) ! list of set1 channels
       INTEGER CLIST2(MXCHN2) ! list of set2 channels
       INTEGER CLIST3(MXCHN3) ! list of set3 channels

       INTEGER  NCHN1         ! # of set1 channels
       INTEGER  NCHN2         ! # of set2 channels
       INTEGER  NCHN3         ! # of set3 channels
       INTEGER  intersect

C      output for FAKETZ and HSUN
       INTEGER  NFAKE              ! # of channels to "fake"
       INTEGER INDFAK(MXCHAN)      ! indices of channels to fake
       INTEGER QUICKINDFAK(MXCHAN) ! list of set1 channels
       REAL   HSUN(MXCHAN) ! sun radiance (direct from sun)

c      local
       INTEGER I,III

C      --------------------------
C      Read in the solar radiance
C      --------------------------

       CALL RDSUN(IOUN, INDCHN, HSUN)
C
C      -----------------------------------------------
C      All channels from sets 1, 2, and 3 are to use a
C      fake effective sun angle layer-to-space trans
c      at the end, NFAKE = NCHN1 + NCHN2 + NCHN3
C      -----------------------------------------------
       NFAKE=0
C 
       DO I=1,NCHN1
          NFAKE=NFAKE + 1
          INDFAK(NFAKE)=INDCHN( CLIST1(I) )
       ENDDO
C
       DO I=1,NCHN2
          NFAKE=NFAKE + 1
          INDFAK(NFAKE)=INDCHN( CLIST2(I) )
       ENDDO
C
       DO I=1,NCHN3
          NFAKE=NFAKE + 1
          INDFAK(NFAKE)=INDCHN( CLIST3(I) )
       ENDDO

      DO I = 1,NFAKE
        III = intersect(I,INDFAK(1:NFAKE), NFAKE)  !! so I = INDFAK(III)
        QUICKINDFAK(I) = III
      END DO

      RETURN
      END
