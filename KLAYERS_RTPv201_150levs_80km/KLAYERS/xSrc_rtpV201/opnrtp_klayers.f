C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              OPNRTP
C
!F77====================================================================


!ROUTINE NAME: OPNRTP


!ABSTRACT:
C    Open RTP input & output files.


!CALL PROTOCOL:
C    OPNRTP(FIN, FOUT, COMMNT, NWANT, LISTG, MASSW, POMIN, POMAX,
C       NGASI, GLISTI, NGASF, GLISTF, GUNITSF, MASSF, INDEXF,
C       IOPCI, IOPCO)



!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    CHAR*80   FIN     input RTP file name         none
C    CHAR*80   FOUT    output RTP file name        none
C    CHAR*256  COMMNT  klayers comment string      none
C    INTEGER   NWANT   # of wanted gases           none
C    INT arr   LISTG   list of wanted gases        none (HITRAN ID #s)
C    REAL arr  MASSW   molec mass of wanted gases  AMU
C    REAL      POMIN   min pressure in output      mb
C    REAL      POMAX   max pressure in output      mb


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   NGASI   # of gases in input file    none
C    INT arr   GLISTI  list of gases in input      none (HITRAN ID #s)
C    INTEGER   NGASF   # of found wanted gases     none
C    INT arr   GLISTF  list of found wanted gases  none (HITRAN ID #s)
C    INT arr   GUNITF  units code # of found gases none
C    REAL arr  MASSF   molec. mass of found gases  AMU
C    INT arr   INDEXF  indices of found gases      none
C    INTEGER   IOPCI   input RTP file I/O unit     none
C    INTEGER   IOPCO   output RTP file I/O unit    none


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): klayers_rtp


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    IOPCI : input RTP file I/O unit ("profile channel")
C    IOPCO : output RTP file I/O unit ("profile channel")


!COMMON BLOCKS:
C    COMUNI : List of allowed input gas amount units


!DESCRIPTION:
C    Opens the input RTP file & reads the header info.  opens the
C    output RTP file and writes the header.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 31 Jan 2001 Scott Hannon      Created
C 19 Feb 2001 Scott Hannon      Add GUNITF
C 27 Feb 2001 Scott Hannon      Added MASSW and MASSF
C 27 Mar 2001 Scott Hannon      Add BITS2N & N2BITS for processing new
C                                  rtp var HEAD.pfields
C  3 Jun 2001 Scott Hannon      Add "BAD" checks for HEAD.mlevs, etc
C 05 Aug 2003 Scott Hannon      Change FIN & FOUT to CHAR*80 (was 70)
C 16 Feb 2007 Scott Hannon      Update HEAD.gunit for GUCCLD
C 07 Nov 2008 Scott Hannon      Minor update for rtpV201

!END====================================================================


C      =================================================================
       SUBROUTINE OPNRTP(FIN, FOUT, COMMNT, NWANT, LISTG, MASSW,
     $    POMIN, POMAX, NGASI, GLISTI, NGASF, GLISTF, GUNITF, MASSF,
     $    INDEXF, IOPCI, IOPCO)
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
C      Input
       CHARACTER*80 FIN        ! input RTP filename
       CHARACTER*80 FOUT       ! output RTP filename
       CHARACTER*256 COMMNT    ! klayers info comment string
       INTEGER  NWANT          ! number of wanted gases
       INTEGER LISTG( MXGAS)   ! list of all wanted gases
       REAL MASSW( MXGAS)      ! molecular mass of all wanted gases
       REAL POMIN              ! min pressure (mb) for output file
       REAL POMAX              ! max pressure (mb) for output file
C
C      Output
       INTEGER  NGASI          ! number of gases in input file
       INTEGER GLISTI( MXGAS)  ! list of gas IDs in input file
       INTEGER  NGASF          ! number of wanted gases found in input
       INTEGER GLISTF( MXGAS)  ! list of wanted gases found in input
       INTEGER GUNITF( MXGAS)  ! gas amount units code number
       REAL     MASSF( MXGAS)  ! molecular mass
       INTEGER INDEXF( MXGAS)  ! indices of wanted gases found in input
       INTEGER  IOPCI  ! I/O unit ("profile channel") for input file
       INTEGER  IOPCO  ! I/O unit ("profile channel") for output file


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER      I
       INTEGER     IC
       INTEGER     IG
       INTEGER    IND
       INTEGER IPATTW(MAXGAS)
       INTEGER IPATTJ(MAXGAS)
       INTEGER      J
       INTEGER      K
       INTEGER  LENNB
       INTEGER  NHATT        ! counter for # of header attributes
       INTEGER  NJUNK        ! counter of junk/unwanted gas patt
       INTEGER  NMOVE        ! number of attributes to move
       INTEGER  NPATT        ! counter for # of profile attributes
       INTEGER STATUS        ! status of RTP file open
       CHARACTER*1 MODE      ! mode for rtpopen: "c"=create, "r"=read
c       CHARACTER*2 IDSTR     ! string for ID numbers
c       CHARACTER*4 C4STR
       CHARACTER*30 CGASES   ! string to be added to end of COMMNT
       LOGICAL LJUNK
C
C      for N2BITS
       INTEGER*4 NUMBER
       LOGICAL LFLAGS(32)
C
C      Structures (see "rtpdefs.f")
       INTEGER rtpopen
       RECORD /RTPHEAD/ HEAD            ! header data
       RECORD /RTPATTR/ HATT(MAXNATTR)  ! header attributes
       RECORD /RTPATTR/ PATT(MAXNATTR)  ! profile attributes


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************

C      -------------------
C      Open RTP input file
C      -------------------
       MODE='r'
       STATUS=rtpopen(FIN, MODE, HEAD, HATT, PATT, IOPCI)
C
ccc
c       print *, 'read open status = ', STATUS
ccc
       IF (HEAD.mlevs   .EQ. BAD) HEAD.mlevs=0
       IF (HEAD.memis   .EQ. BAD) HEAD.memis=0
ccc Removed for rtpV201
c       IF (HEAD.mrho    .EQ. BAD) HEAD.mrho=0
c       IF (HEAD.mwmemis .EQ. BAD) HEAD.mwmemis=0
c       IF (HEAD.mwmstb  .EQ. BAD) HEAD.mwmstb=0
ccc

C      ---------------------
C      Check pfields & ptype
C      ---------------------
       IF (HEAD.ptype .EQ. LAYPRO) THEN
          WRITE(IOERR,2010)
 2010     FORMAT('ERROR: KLAYERS can only process RTP ptype=LEVPRO')
          STOP
       ENDIF
C
       NUMBER=HEAD.pfields
       CALL N2BITS(NUMBER,LFLAGS)
ccc
c      print *, NUMBER
c      print *, LFLAGS
ccc
       IF (.NOT. LFLAGS(1)) THEN
          WRITE(IOERR,2012)
 2012     FORMAT('ERROR: input RTP file has no profile data!')
          STOP
       ENDIF
C

C      ----------------------------
C      Pull out input file gas info
C      ----------------------------
       NGASI=HEAD.ngas
       DO I=1,NGASI
          GLISTI(I)=HEAD.glist(I)
       ENDDO
C

C      ---------------------------------------------
C      Detemine which wanted gases are in input file
C      ---------------------------------------------
C      Loop over wanted gases
       NGASF=0
       DO I=1,NWANT
C         Loop over gases in input file
          DO J=1,NGASI
C            Check if current gas in input file is wanted
             IF (GLISTI(J) .EQ. LISTG(I)) THEN
C               Increment count of found gases
                NGASF=NGASF + 1
C               ID of found gas
                GLISTF(NGASF)=HEAD.glist(J)
C               Gas amount units code number
                GUNITF(NGASF)=HEAD.gunit(J)
C               Molecular mass
                MASSF(NGASF)=MASSW(I)
C               Index of found gas
                INDEXF(NGASF)=J
C
C               Check gunits to see if input gas amount units are OK
                CALL CHKGUC( GUNITF(NGASF), IND )
                IF (IND .LT. 1) THEN
                   WRITE(IOERR,1010) J, GLISTF(NGASF), GUNITF(NGASF)
 1010              FORMAT('ERROR! index=',I2,', gas ID=',I3,
     $             ', KLAYERS can not handle units code=',I3)
                   STOP
                ENDIF
C
             ENDIF
          ENDDO
       ENDDO
C

C      -----------------------------
C      Add gas info to end of COMMNT
C      -----------------------------
C      Assign CGASES string (30 char max)
       CGASES='; Gases from RTP input='
C
C      Add gas info to end of COMMNT
       I=LENNB(CGASES)  ! length excluding trailing blanks
       J=LENNB(COMMNT)
       K=I + 3*NGASF ! number of chars to add to COMMNT
       IF (K .LE. 256-J) THEN
          COMMNT(J+1:J+I)=CGASES(1:I)
          J=J + I  ! new index of last non-black char
C         Loop over found gases
          DO K=1,NGASF
C           Write gas ID to end of COMMNT
            IF (GLISTF(K) .LT. 10) THEN
               WRITE(COMMNT(J+1:J+2),1021) GLISTF(K)
 1021          FORMAT(I2)
               J=J + 2
            ELSE
               IF (GLISTF(K) .LT. 100) THEN
                  WRITE(COMMNT(J+1:J+3),1022) GLISTF(K)
 1022             FORMAT(I3)
                  J=J + 3
               ELSE
                  WRITE(COMMNT(J+1:J+4),1023) GLISTF(K)
 1023             FORMAT(I4)
                  J=J + 4
               ENDIF
            ENDIF
          ENDDO
       ENDIF
C
C      -----------------------------------
C      Update header attributes for output
C      -----------------------------------
C      Add klayers comment to header attributes
C      Count the number of header attributes
       I=1
       DO WHILE (ICHAR(HATT(I).fname) .NE. 0 .AND. I .LE. MAXNATTR)
          I=I + 1
       ENDDO
       NHATT=I
       J=LENNB(COMMNT)
       HATT(NHATT).fname='header'  // CHAR(0)
       HATT(NHATT).aname='klayers' // CHAR(0)
       HATT(NHATT).atext=COMMNT(1:J) // CHAR(0)
C      Add a char(0) to end of attributes if less than maxnattr
       IF (NHATT .LT. MAXNATTR) THEN
          HATT(NHATT + 1).fname=CHAR(0)
       ENDIF
C
C      Update pmin & pmax
       HEAD.pmin=POMIN
       HEAD.pmax=POMAX
C
C      Update ngas & glist & gunit
       HEAD.ngas=NWANT
       DO I=1,NWANT
          HEAD.glist(I)=LISTG(I)
          HEAD.gunit(I)=GUCOUT
          IF (LISTG(I) .GE. IDCLD1) HEAD.gunit(I)=GUCCLD
       ENDDO
C
C      Update mlevs
       HEAD.mlevs=NBPLEV
C
C      Update ptype for output
       HEAD.ptype=LAYPRO  ! layer profile
C

C      ------------------------------------
C      Update profile attributes for output
C      ------------------------------------
C      Update attribute entry (if any) for wanted gases.
C      Remove attribute entry for not-wanted gases.
c
C      Initialize array of attribute indices of wanted gases
       DO IG=1,NWANT
          IPATTW(IG)=-1
       ENDDO
C      Initialize array of attribute indices of junk/unwanted gases
       NJUNK=0
       DO IG=1,MAXGAS
          IPATTJ(IG)=-1
       ENDDO
C
cccccccc this block for debugging ccccccccccccccccccccc
c       I=1
c       WHILE (ICHAR(PATT(I).fname) .NE. 0 .AND. I .LE. MAXNATTR)
c          print *, '-----  input patt number ', I, '-------'
c          print *, PATT(I).fname
c          print *, PATT(I).aname
c          print *, PATT(I).atext
c          I=I+1
c       ENDDO
ccccccccccccccccccccccccccccc

C      --------------------------------------
C      Check input rtp gas_i units attributes
C      --------------------------------------
       I=1
       DO WHILE (ICHAR(PATT(I).fname) .NE. 0 .AND. I .LE. MAXNATTR)
C
          IF ( PATT(I).fname(1:4) .EQ. 'gas_' .AND.
     $         PATT(I).aname(1:5) .EQ. 'units') THEN
C
ccc Do not bother pulling out gas ID
cC            Pull out gas ID (integer) from fname (string)
c             IC=ICHAR( PATT(I).fname(6:6) )
c             IF (IC .GE. ICHAR('0') .AND. IC .LE. ICHAR('9')) THEN
c                READ( PATT(I).fname(5:6), *) IG
c             ELSE
c                READ( PATT(I).fname(5:5), *) IG
c             ENDIF
ccc
C
C            ---------------------------
C            Check & update wanted gases
C            ---------------------------
             LJUNK=.TRUE.
ccc Get rid of all gas_# units attributes
c             DO J=1,NWANT  ! loop over wanted gases
c                IF (IG .EQ. LISTG(J)) THEN
cC
cC                  Check that units are ppmv
c                   C4STR=PATT(I).atext(1:4)
c                   CALL UPCASE(C4STR)
c                   IF (C4STR .NE. 'PPMV') THEN
c                      WRITE(IOERR,1030) IG, PATT(I).atext(1:20)
c 1030                 FORMAT('ERROR! units for gas ID ',I2,
c     $                ' are ',A20,' instead of PPMV.')
c                      STOP
c                   ENDIF
cC
cC                  Update units
c                   IPATTW(J)=I  ! index of attribute
c                   LJUNK=.FALSE.
cC                  Update 'units'
c                   PATT(I).atext='kilomoles/cm^2' // CHAR(0)
cc                   PATT(I).atext='molecules/cm^2' // CHAR(0)
cC
c                ENDIF
c             ENDDO
ccc
C
C            ------------------------------------------
C            Keep track of junk/unwanted gas attributes
C            ------------------------------------------
             IF (LJUNK) THEN
                NJUNK=NJUNK + 1  ! count of junk/unwanted attributes
                IPATTJ(NJUNK)=I  ! index of this unwanted attribute
             ENDIF
C
          ENDIF
          I=I + 1
       ENDDO
       NPATT=I - 1  ! count of profile attributes in input file
C
C      ---------------------------------------------------------------
C      Add gas_i units attributes for wanted gases not in original rtp
C      ---------------------------------------------------------------
ccc Do not add any gas_# units attriubtes
c       DO J=1,NWANT
c          IF (IPATTW(J) .LT. 0) THEN
c
c       print *,'doing new entry for wanted gas J=',J
c
cC            Create ID number string
c             IF (LISTG(J) .GT. 9) THEN
c                WRITE(IDSTR,1040) LISTG(J)
c 1040           FORMAT(I2)
c             ELSE
c                WRITE(IDSTR,1041) LISTG(J), CHAR(0)
c 1041           FORMAT(I1,A1)
c             ENDIF
c             IF (NJUNK .GT. 0) THEN
cC               Re-cycle an old junk units attribute
c                IC=IPATTJ(NJUNK)  ! index of a junk attribute
c                PATT(IC).fname='gas_' // IDSTR // CHAR(0)
c                PATT(IC).atext='kilomoles/cm^2' // CHAR(0)
cc                PATT(IC).atext='molecules/cm^2' // CHAR(0)
c                NJUNK=NJUNK - 1
c             ELSE
cC               Create a new units attribute
c                IC=NPATT + 1
c                IF (IC .GT. MAXNATTR) THEN
c                   WRITE(IOERR,1050)
c 1050              FORMAT('ERROR! output profile attributes '
c     $             'exceeds MAXNATTR')
c                   STOP
c                ELSE
c                   PATT(IC).fname='gas_' // IDSTR // CHAR(0)
c                   PATT(IC).aname='units' // CHAR(0)
c                   PATT(IC).atext='kilomoles/cm^2' // CHAR(0)
cc                   PATT(IC).atext='molecules/cm^2' // CHAR(0)
c                ENDIF
c                NPATT=NPATT + 1
c             ENDIF
c          ENDIF
c       ENDDO
ccc
C
C      ----------------------------------------------------
C      Remove any remaining unwanted gas_i units attributes
C      ----------------------------------------------------
       IF (NJUNK .GT. 0) THEN
C         New number of profile attributes for output
          NPATT=NPATT - NJUNK
C
C         Count number of junk attributes between index 1 & NPATT
          NMOVE=0
          DO I=1,NJUNK
             IF (IPATTJ(I) .LE. NPATT) NMOVE=NMOVE + 1
          ENDDO
C
C         Loop over good attributes to move
          IC=NPATT + 1  ! initialize good index to be moved
          DO I=1,NMOVE
             J=IPATTJ(I)  ! index of junk attr to overwrite
C
C            Find index of non-junk attribute to be moved
             DO K=NMOVE+I,NJUNK
                IF (IC .EQ. IPATTJ(K)) IC=IC + 1
              ENDDO
C
C            Move good attribute from index IC to index J
             PATT(J).fname=PATT(IC).fname
             PATT(J).aname=PATT(IC).aname
             PATT(J).atext=PATT(IC).atext
             IC=IC + 1  ! increment
          ENDDO
       ENDIF
       IF (NPATT .LT. MAXNATTR) PATT(NPATT + 1).fname = CHAR(0)
C
cccccc this block for debugging ccccccccccccccccccccccc
c       I=1
c       WHILE (ICHAR(PATT(I).fname) .NE. 0 .AND. I .LE. MAXNATTR)
c          print *, '-----  output patt number ', I, '-------'
c          print *, PATT(I).fname
c          print *, PATT(I).aname
c          print *, PATT(I).atext
c          I=I+1
c       ENDDO
ccccccccccccccccccccccccccccc

C      --------------------
C      Open RTP output file
C      --------------------
       MODE='c'
       STATUS=rtpopen(FOUT, MODE, HEAD, HATT, PATT, IOPCO)

ccc
c       print *, 'read open status = ', STATUS
ccc

C
       RETURN
       END
