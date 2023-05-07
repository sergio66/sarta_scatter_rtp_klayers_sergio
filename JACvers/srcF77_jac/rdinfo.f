c rdinfo processes command line arguments
c
c       sarta  fin=input.rtp  fout=output.rtp  listp=1,2,3 listc=445,449,1092,1291,1614,2070,2333,2353 listj=-1
c          would give you radiances at                           790,791,1042,1231,1419,2350,2616,2637 cm-1
c                                                                 CO2Q    O3   WIN  WV  NLTE  WIN  HDO
c
c
c to compile
c   Absoft/Linux: f77 -N109 -o klayers $(SRC) -lU77
c   SGI Irix: no special compiler options are needed
C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              RDINFO_sarta
C
!F77====================================================================


!ROUTINE NAME: RDINFO


!ABSTRACT:
C    Get info about the sarta run: the names of input & output
C    files, the channel list, and list of profile numbers.


!CALL PROTOCOL:
C       SUBROUTINE RDINFO(FIN, FOUT, LRHOT, NWANTP, LISTP, NWANTC, LISTC,
C     $             NWANTJ, LISTJ, NUMCHAN, NUMPROF, caJacTZ, caJacG1, caJacG3, caJacWgt)

!INPUT PARAMETERS:
C    none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    CHAR*80   FIN     input filename              none
C    CHAR*80   FOUT    output filename             none
C    LOGICAL   LRHOT   force RHO for refl thermal? none
C    INTEGER   NWANTP  Number of desired profiles  none
C    INT arr   LISTP   List of desired prof nums   none
C    INTEGER   NWANTC  Number of desired channels  none
C    INT arr   LISTC   List of desired channels nums   none
C    INTEGER   NWANTJ  Number of desired jacs      none
C    INT arr   LISTJ   List of desired jac nums    none

!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): SARTA_rtp


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    none


!COMMON BLOCKS:
C      none


!DESCRIPTION:
C    Gets various info about sarta run.
C
C    Each command line argument is of the form <variable>=<value>
C    Each <variable>=<value> string must be 80 char or less.  It is
C    necessary to enclose values in quotes unless they contain
C    blanks.  The recognized command-line variables are:
C
C    fin : name of input file
C
C    fout : name of output file
C
C    lrhot : force reflected thermal rho?; logical. true/false of T/F
C       If true, the refl therm will use rho=(1-emis)/pi rather than
C       the rho (if any) from the input file.
C
C    listp : list of desired profile numbers (all other profiles will
C       be ignored).  If "listp" is not specified, SARTA will process
C       all profiles.
C
C       The listp profile numbers may be specified either as a
C       sequence of integers separated by a comma, or alternately as
C       a quoted string containing integers separated by a blank space.
C       Examples:
C          listp=1,2,3,4,5
C          listp='1 2 3 4 5'
C       Due to the 80 char limit, the maximum number of entries
C       in listp is limited.  (Eg 15 four digit numbers, or
C       25 two digit numbers.  MAXPRO is the hardcoded limit.)
C
C    listc : list of desired channel numbers (all other channels will
C       be ignored).  If "listc" is not specified, SARTA will process
C       all channels in the rtp file (eg if h.nchan = 400, all 400 chans will be used).
C
C       The listc profile numbers may be specified either as a
C       sequence of integers separated by a comma, or alternately as
C       a quoted string containing integers separated by a blank space.
C       Examples:
C          listc=1,2,3,4,5
C          listc='1 2 3 4 5'
C       Due to the 80 char limit, the maximum number of entries
C       in listc is limited.  (Eg 15 four digit numbers, or
C       25 two digit numbers.  MAXPRO is the hardcoded limit.)
C
C    listj : list of desired jacs (default 0 = none)
C
C       The listj profile numbers may be specified either as a
C       sequence of integers separated by a comma, or alternately as
C       a quoted string containing integers separated by a blank space.
C       Examples:
C          listj= -1,1,2,3,4,5
C          listj='-1 1 2 3 4 5'
C       Due to the 80 char limit, the maximum number of entries
C       in listc is limited.  (Eg 15 four digit numbers, or
C       25 two digit numbers.  MAXPRO is the hardcoded limit.)

!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 13 Feb 2001 H.Motteler/S.Hannon Re-write of KLAYERS version
C 28 Nov 2001 Scott Hannon      Remove command-line argument "nwantp"
C  5 Dec 2001 Scott Hannon      Remove unused local var LENNB
C 05 Aug 2003 Scott Hannon      Correct FIN & FOUT to CHAR*80 (not 70)
C 06 Feb 2004 Scott Hannon      Add LRHOT argument and associated code


!END====================================================================


C      =================================================================
       SUBROUTINE RDINFO(FIN, FOUT, LRHOT, NWANTP, LISTP, NWANTC, LISTC,
     $             NWANTJ, LISTJ, NUMCHAN, NUMPROF, 
     $     caJacTZ,caJACWGT,caJACG1,caJACG2,caJACG3,caJACG4,caJACG5,
     $     caJACG6,caJACG9,caJacG11,caJACG12,caJacG103)
C      =================================================================

C      use unix_library
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
       include 'incFTC.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      From "util.f"
C      subroutine UPCASE = converts a string to upper case
C      function STR2BO = converts true/false string to boolean (LOGICAL)


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input:
C      none
C
C      Output:
       CHARACTER*80 FIN
       CHARACTER*80 FOUT
       LOGICAL  LRHOT
       INTEGER NWANTP
       INTEGER  LISTP(MAXPRO)
       INTEGER NWANTC
       INTEGER  LISTC(MAXPRO)
       INTEGER NWANTJ,NUMPROF,NUMCHAN,XNUMPROF,XNUMCHAN
       INTEGER  LISTJ(MAXPRO)
       CHARACTER*180 caJacTZ,caJACWGT,caJACG1,caJACG2,caJACG3,caJACG4,caJACG5,
     $               caJACG6,caJACG9,caJACG11,caJACG12,caJACG103

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER I
       INTEGER IARGC
       INTEGER IP
       INTEGER IPJUNK(MXGAS+1)  ! junk gas id work array
       INTEGER J
       INTEGER K
       INTEGER NARGS   ! number of arguments
       INTEGER SORTED  ! flag for sorting

       CHARACTER*80 BUF
       CHARACTER*80 VAL
       CHARACTER*80 VAR

       LOGICAL LLISTP
       LOGICAL LLISTC
       LOGICAL LLISTJ

       INTEGER NWANTJX
       INTEGER  LISTJX(MAXPRO)
       
       LOGICAL STR2BO

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************

C      ------------
C      Set defaults
C      ------------
       FIN='sarta_in.rtp'                 ! input filename
       FOUT='sarta_out.rtp'               ! output filename
       NWANTP = -1   ! do sarta for all profiles found in input file
       NWANTC = -1   ! do sarta for all channels found in input file
       NWANTJ = 0    ! do sarta rads only (no jacs)

       NUMPROF = NWANTP
       NUMCHAN = NWANTC
       caJacTZ  = ' '
       caJacG1  = ' '
       caJacG3  = ' '
       caJacWGT = ' '
       LRHOT=.FALSE. ! use input rho for reflected thermal
C
C      -----------------------------------------------------------------
C      Loop on program parameters
C      --------------------------
C      Determine the number of command-line arguments
       NARGS=IARGC()
C
C      Loop over the command-line arguments
       LLISTP=.FALSE.
       DO I = 1, NARGS
C
C         Pull out the ith argument
          CALL GETARG(I, BUF)
C
C         Find the "=" character in the command-line argument string
          J=INDEX(BUF, '=')
C
          IF (J .NE. 0) THEN
C
C            Name of variable
             VAR = BUF(1:J-1)
             CALL UPCASE(VAR)
C
C            Specified value
             VAL = BUF(J+1:LEN(BUF))
C
C            Big "IF" to set parameters
C            ----------------------------
             IF (VAR(1:3) .EQ. 'FIN') THEN
                FIN=VAL

             ELSEIF (VAR(1:4) .EQ. 'FOUT') THEN
                FOUT=VAL

             ELSEIF (VAR(1:5) .EQ. 'LRHOT') THEN
                LRHOT=STR2BO(VAL)

             ELSEIF (VAR(1:5) .EQ. 'LISTP') THEN
                LLISTP=.TRUE.

C               Read the indices of the desired profiles
                K=1
 10             IF (K .GT. MAXPRO) THEN
                   WRITE(6,1017)
 1017              FORMAT('ERROR! bad LISTP, ',
     $             'either an unrecognized value or too many entries')
                   STOP
                ENDIF
                READ(VAL,*,END=19) (IPJUNK(IP),IP=1,K)
                K=K + 1  ! increment count of profiles
                GOTO 10  ! loop to next entry
 19             CONTINUE
                K=K - 1  ! number of profiles
                DO IP=1,K
                   LISTP(IP)=IPJUNK(IP)
                ENDDO
                NWANTP=K
C

             ELSEIF (VAR(1:5) .EQ. 'LISTC') THEN
                LLISTC=.TRUE.

C               Read the indices of the desired channels
                K=1
 20             IF (K .GT. MAXPRO) THEN
                   WRITE(6,2017)
 2017              FORMAT('ERROR! bad LISTC, ',
     $             'either an unrecognized value or too many entries')
                   STOP
                ENDIF
                READ(VAL,*,END=29) (IPJUNK(IP),IP=1,K)
                K=K + 1  ! increment count of profiles
                GOTO 20  ! loop to next entry
 29             CONTINUE
                K=K - 1  ! number of profiles
                DO IP=1,K
                   LISTC(IP)=IPJUNK(IP)
                ENDDO
                NWANTC=K

             ELSEIF (VAR(1:5) .EQ. 'LISTJ') THEN
                LLISTJ=.TRUE.

C               Read the indices of the desired channels
                K=1
 30             IF (K .GT. MAXPRO) THEN
                   WRITE(6,3017)
 3017              FORMAT('ERROR! bad LISTJ, ',
     $             'either an unrecognized value or too many entries')
                   STOP
                ENDIF
                READ(VAL,*,END=39) (IPJUNK(IP),IP=1,K)
                K=K + 1  ! increment count of profiles
                GOTO 30  ! loop to next entry
 39             CONTINUE
                K=K - 1  ! number of profiles
                DO IP=1,K
                   LISTJ(IP)=IPJUNK(IP)
                ENDDO
                NWANTJ=K

             ELSE
                WRITE(6,1020) VAR(1:6)
 1020           FORMAT('Unknown command-line argument: ',A6)
                STOP

             ENDIF

          ENDIF
       ENDDO  ! end of loop over command-line arguments
C      -----------------------------------------------------------------


C      -------------------------------------
C      Sort prof numbers & check for repeats
C      -------------------------------------
       IF ((NWANTP .EQ. 1) .AND. (LISTP(1) .EQ. -1)) NWANTP = -1

       IF (NWANTP .GT. 0) THEN
C
C         Sort in ascending order
          SORTED=1  ! initialize flag for first pass
 120       IF (SORTED .EQ. 1) THEN
             SORTED=0  ! initialize flag for this loop
             DO K=1,NWANTP-1
                IF (LISTP(K) .GT. LISTP(K+1)) THEN
                   IP=LISTP(K)
                   LISTP(K)=LISTP(K+1)
                   LISTP(K+1)=IP
                   SORTED=1  ! set flag to indicate ordering was altered
                ENDIF
             ENDDO
             GOTO 120
          ENDIF
C
C         Check for repeats
          DO K=1,NWANTP-1
             IF (LISTP(K) .EQ. LISTP(K+1)) THEN
                WRITE(6,1045) LISTP(K)
 1045           FORMAT('ERROR! profile ',I2,
     $          ' appears more than once in LISTP')
                STOP
             ENDIF
          ENDDO

       ENDIF
C

C      -------------------------------------
C      Sort chan numbers & check for repeats
C      -------------------------------------
       IF ((NWANTC .EQ. 1) .AND. (LISTC(1) .EQ. -1)) NWANTC = -1

       IF ((NWANTC .EQ. 1) .AND. (LISTC(1) .EQ. -100)) THEN
         write (*,'(A)') 'LISCT = -100 so doing AIRS Channel ID 445,449,1092,1291,1614,2070,2333,2353'
         NWANTC = 8
         LISTC = 0
         LISTC(1:NWANTC) = (/445,449,1092,1291,1614,2070,2333,2353/)
       END IF

       IF (NWANTC .GT. 0) THEN
C
C         Sort in ascending order
          SORTED=1  ! initialize flag for first pass
 60       IF (SORTED .EQ. 1) THEN
             SORTED=0  ! initialize flag for this loop
             DO K=1,NWANTC-1
                IF (LISTC(K) .GT. LISTC(K+1)) THEN
                   IP=LISTC(K)
                   LISTC(K)=LISTC(K+1)
                   LISTC(K+1)=IP
                   SORTED=1  ! set flag to indicate ordering was altered
                ENDIF
             ENDDO
             GOTO 60
          ENDIF
C
C         Check for repeats
          DO K=1,NWANTC-1
             IF (LISTC(K) .EQ. LISTC(K+1)) THEN
                WRITE(6,2045) LISTC(K)
 2045           FORMAT('ERROR! channel ',I2,
     $          ' appears more than once in LISTC')
                STOP
             ENDIF
          ENDDO

       ENDIF

C      -------------------------------------
C      Sort jac numbers & check for repeats
C      -------------------------------------
       IF ((NWANTJ .EQ. 1) .AND. (LISTJ(1) .EQ. 0)) NWANTJ = 0
       IF (NWANTJ .GT. 0) THEN
C
C         Sort in ascending order
          SORTED=1  ! initialize flag for first pass
 90       IF (SORTED .EQ. 1) THEN
             SORTED=0  ! initialize flag for this loop
             DO K=1,NWANTJ-1
                IF (LISTJ(K) .GT. LISTJ(K+1)) THEN
                   IP=LISTJ(K)
                   LISTJ(K)=LISTJ(K+1)
                   LISTJ(K+1)=IP
                   SORTED=1  ! set flag to indicate ordering was altered
                ENDIF
             ENDDO
             GOTO 90
          ENDIF
C
C         Check for repeats
          DO K=1,NWANTJ-1
             IF (LISTJ(K) .EQ. LISTJ(K+1)) THEN
                WRITE(6,3045) LISTJ(K)
 3045           FORMAT('ERROR! channel ',I2,
     $          ' appears more than once in LISTJ')
                STOP
             ENDIF
          ENDDO

          NUMPROF = NWANTP
          NUMCHAN = NWANTC
          IF ((NWANTP .EQ. -1) .OR. (NWANTC .EQ. -1)) THEN
            CALL find_num_prof(FIN,XNUMCHAN,XNUMPROF)
            IF ((NWANTP .EQ. -1) .AND. (NWANTC .EQ. -1)) THEN
               NUMPROF = XNUMPROF
               NUMCHAN = XNUMCHAN
            ELSEIF ((NWANTP .EQ. -1) .AND. (NWANTC .GT. 0)) THEN
              NUMPROF = XNUMPROF            
            ELSEIF ((NWANTP .GT. 0) .AND. (NWANTC .EQ. -1)) THEN
              NUMCHAN = XNUMCHAN
            END IF
          ENDIF

C remember : 0 = NO jacs, -1 = all jacs (T,WV,O3), 1,3 = only G1 G2 G3 G4 G5 G6 G9 G12, 100 = only T+ST, 200 = WGT
          NWANTJX  = NWANTJ
          LISTJX   = LISTJ
          caJacTZ  = trim(trim(FOUT) // '_jacTZ')
          caJACG1  = trim(trim(FOUT) // '_jacG1')
          caJACG2  = trim(trim(FOUT) // '_jacG2')
          caJACG3  = trim(trim(FOUT) // '_jacG3')
          caJACG4  = trim(trim(FOUT) // '_jacG4')
          caJACG5  = trim(trim(FOUT) // '_jacG5')
          caJACG6  = trim(trim(FOUT) // '_jacG6')
          caJACG9  = trim(trim(FOUT) // '_jacG9')
          caJACG11 = trim(trim(FOUT) // '_jacG11')
          caJACG12 = trim(trim(FOUT) // '_jacG12')
          caJACG103= trim(trim(FOUT) // '_jacG103')
          caJACWGT = trim(trim(FOUT) // '_WGTFCN')

          IF ((NWANTJX .EQ. 1) .AND. (LISTJX(1) .EQ. -1)) THEN
            LISTJ = 0
            NWANTJ = 4
            LISTJ(1) = 1
            LISTJ(2) = 3
            LISTJ(3) = 100
            LISTJ(4) = 200

            LISTJ = 0
            NWANTJ = 8
            LISTJ(1) = 1
            LISTJ(2) = 2
            LISTJ(3) = 3
            LISTJ(4) = 4
            LISTJ(5) = 5
            LISTJ(6) = 6
            LISTJ(7) = 100
            LISTJ(8) = 200

            LISTJ = 0
            NWANTJ = 10
            LISTJ(1) = 1
            LISTJ(2) = 2
            LISTJ(3) = 3
            LISTJ(4) = 4
            LISTJ(5) = 5
            LISTJ(6) = 6
            LISTJ(7) = 9
            LISTJ(8) = 12
            LISTJ(9)  = 100
            LISTJ(10) = 200
          END IF            
       ENDIF


       RETURN
       END
