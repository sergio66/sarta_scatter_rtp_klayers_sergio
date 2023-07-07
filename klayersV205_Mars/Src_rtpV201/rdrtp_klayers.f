C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              RDRTP
C
!F77====================================================================


!ROUTINE NAME: RDRTP


!ABSTRACT:
C    Read a profile from a previously openned RTP file


!CALL PROTOCOL:
C    RDRTP(ISTAT, IOPCI, NGASF, GLISTF, INDEXF, LAT,  LON,
C       NIN, PIN, LNPIN, TIN, MRIN, ZIN, LZ, PSURF, ZSURF, PROF)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   IOPCI   input RTP file I/O number   none
C    INTEGER   NGASF   # of found wanted gases     none
C    INT arr   GLISTF  list of found wanted gases  none (HITRAN IDs)
C    INT arr   INDEXF  indices of found gases      none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   ISTAT   I/O error status            none
C    REAL      LAT     latitude                    degrees
C    REAL      LON     longitude                   degrees
C    INTEGER   NIN     number of profile levels    none
C    REAL arr  PIN     air pressure                atm
C    REAL arr  LNPIN   log of air pressure         log of atm
C    REAL arr  TIN     air temperature             Kelvin
C    REAL arr  MRIN    gas mixing ratios           ppmv
C    REAL arr  ZIN     level altitudes             meters
C    LOGICAL   LZ      true/false ZIN supplied     none
C    REAL      PSURF   surface pressure            mb
C    REAL      ZSURF   surface altitude            meters
C    STRUCT    PROF    profile data structure      (see attributes)


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): KLAYERS


!ROUTINES CALLED: none


!FILES ACCESSED:
C    Input RTP file with I/O number IOPCI


!COMMON BLOCKS: none


!DESCRIPTION:
C    Reads a single profile from a previously openned RTP file.
C    The routine expects to find the data specified in the header
C    of the input RTP file.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 31 Jan 2001 Scott Hannon      created
C  9 Feb 2001 Scott Hannon      ensure returned profile is bottom-up
C 23 Mar 2001 Scott Hannon      add output vars PSURF & ZSURF
C 30 Aug 2001 Scott Hannon      Bug fix: top-down profiles now work
C  1 May 2002 Scott Hannon      Add error trap for NIN > MAXLEV
C 06 Oct 2004 Scott Hannon      Add error trop for |LAT| > 90
C 30 Aug 2011 Scott Hannon      Add error trap for NIN < 10

!END====================================================================

C      =================================================================
       SUBROUTINE RDRTP(ISTAT, IOPCI, NGASF, GLISTF, INDEXF, LAT, LON,
     $    NIN, PIN, LNPIN, TIN, MRIN, ZIN, LZ, PSURF, ZSURF, PROF)
C      =================================================================


C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
      include 'incLAY.f'
      include 'rtpdefs.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input parameters:
       INTEGER IOPCI
       INTEGER  NGASF
       INTEGER  GLISTF( MXGAS)
       INTEGER  INDEXF( MXGAS)
C
C      Output parameters:
       INTEGER ISTAT
       REAL    LAT
       REAL    LON
       INTEGER    NIN
       REAL    PIN(  MXIN)
       REAL  LNPIN(  MXIN)
       REAL    TIN(  MXIN)
       REAL   MRIN(  MXIN, MXGAS)
       REAL    ZIN(  MXIN)
       LOGICAL LZ
       REAL  PSURF
       REAL  ZSURF
C      Profile data structure
       RECORD /RTPPROF/ PROF


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER IG,J
       INTEGER IL
       INTEGER ILR
       INTEGER rtpread
       REAL P1
       REAL PNIN
       REAL ZJUNK


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************

C      ------------------------
C      Read the current profile
C      ------------------------
       ISTAT=rtpread(IOPCI, PROF)
C
       IF (ISTAT .EQ. -1) GOTO 9999  ! reached end of file
C
C      -----------------------------------
C      Pull out the data needed by KLAYERS
C      -----------------------------------
       LAT=PROF.plat
       LON=PROF.plon
       NIN=PROF.nlevs
       PSURF=PROF.spres
       ZSURF=PROF.salti
C
C      Can not correctly compute gravity if LAT is bad
       IF (ABS(LAT) .GT. 90.01) THEN
          WRITE(IOERR,1005) LAT
 1005     FORMAT('ERROR! Bad profile latitude LAT=',1PE11.4)
          STOP
       ENDIF
C      Note: LON is not currently used and thus need not be checked
C
C      Can not correctly read RTP if NIN=PROF.nlevs > MAXLEV
       IF (NIN .GT. MAXLEV) THEN
          WRITE(IOERR,1010) NIN, MAXLEV
 1010     FORMAT('ERROR! Profile NLEVS=',I4,' > RTP MAXLEV=',I4)
          STOP
       ENDIF
C
C      Assume a problem if NIN < 10
       IF (NIN .LT. 10) THEN
          WRITE(IOERR,1015) NIN
 1015     FORMAT('ERROR! Profile NLEVS=',I4,' < 10')
          STOP
       ENDIF
C
C      Loop over the levels
C      Note: returned profile data must be sorted bottom-up
       ZJUNK=-9999.0
       P1=PROF.plevs(1)
       PNIN=PROF.plevs(NIN)
C
       IF (PNIN .GT. P1) THEN
C         Input profile sorted top-down (Pmin to Pmax)
          DO IL=1,NIN
             ILR=NIN + 1 - IL  ! reverse order of loop index
             PIN(ILR)=PROF.plevs(IL)
             LNPIN(ILR)=LOG( PIN(ILR) )
             TIN(ILR)=PROF.ptemp(IL)
             ZIN(ILR)=PROF.palts(IL)
             ZJUNK=MAX( ZJUNK, ZIN(ILR) )
C            Loop over the wanted gases
             DO IG=1,NGASF
                MRIN(ILR,IG)=PROF.gamnt(IL,INDEXF(IG))
             ENDDO
          ENDDO
C
       ELSE
C         Input profile sorted bottom-up (Pmax to Pmin)
          DO IL=1,NIN
             PIN(IL)=PROF.plevs(IL)
             LNPIN(IL)=LOG( PIN(IL) )
             TIN(IL)=PROF.ptemp(IL)
             ZIN(IL)=PROF.palts(IL)
             ZJUNK=MAX( ZJUNK, ZIN(IL) )
C            Loop over the wanted gases
             DO IG=1,NGASF
                MRIN(IL,IG)=PROF.gamnt(IL,INDEXF(IG))
             ENDDO
          ENDDO
       ENDIF
C
C      Check altitude
       IF (ZJUNK .GT. 0.0) THEN
          LZ=.TRUE.
       ELSE
          LZ=.FALSE.
       ENDIF
       IF (PROF.palts(1) .LT. -998) LZ=.FALSE.
C

ccc for testing ccc
c       print *, 'LZ=',LZ
c       DO IL=1,NIN
c          WRITE(6,1234) IL, PIN(IL), LNPIN(IL), TIN(IL), ZIN(IL),
c     $       (MRIN(IL,IG),IG=1,NGASF)
c 1234     FORMAT(I4,2(1PE12.5),0PF8.3,6(1PE12.5))
c       ENDDO
ccccccccc

C
 9999  RETURN
       END
