c rdinfo processes command line arguments
c
c       klayers fin=test.rtp fout=out.rtp nwant=3 listg=1,2,3
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
C              RDINFO_cl
C
!F77====================================================================


!ROUTINE NAME: RDINFO


!ABSTRACT:
C    Get info about the user supplied profile.  Gets info from
C    command line arguments else uses the hardcoded defaults.
C    This version for command-line arguments is a total re-write,
C    and even a few of the passed parameters are different.


!CALL PROTOCOL:
C    RDINFO(FIN, FOUT, COMMNT, FNAFGL, MNAFGL,
C       NWANT, LISTG, MASSW, TOFF, SCALEH, LDRY, LSVP, LSPLIN, LPMAX)


!INPUT PARAMETERS:
C    none


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    CHAR*80   FIN     user profile file name      none
C    CHAR*80   FOUT    user profile file name      none
C    CHAR*256  COMMNT  KLAYERS info comment        none
C    CHAR*80   FNAFGL  AFGL profile filename       none
C    INTEGER   MNAFGL  AFGL profile model number   none (1-6)
C    INTEGER   NWANT   number of wanted gases      none
C    INT arr   LISTG   List of wanted gases        none (HITRAN ID #s)
C    REAL arr  MASSW   Masses of wanted gases      AMU
C    REAL      TOFF    offset to apply to temp.    Kelvin
C    REAL      SCALEH  Pres vs Alt scale height    kilometers
C    LOGICAL   LDRY    T/F mix ratios for dry air  none
C    LOGICAL   LSVP    T/F check water saturation  none
C    LOGICAL   LSPLIN  T/F do spline interpolation none
C    LOGICAL   LPMAX   T/F extend profile to Pmax  none


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): KLAYERS


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    Input unit IOUNIT (used by RDUPRO)
C    Standard input (keyboard) unit 5


!COMMON BLOCKS:
C      Allowed gas IDs from "cbgids.f"
C      COMMON /COMGID/ GIDS, GMASS
C      INTEGER GIDS(NGIDS)
C      REAL GMASS(NGIDS)


!DESCRIPTION:
C    Gets various info about the user supplied profile to use
C    with kLAYERS.  The info comes from the command line arguments.
C    Default values are used for all non-specified arguments.
C
C    Each command line argument is of the form <variable>=<value>
C    Each <variable>=<value> must be 80 char or less.  It is not
C    necessary to enclose string values in quotes. The recognized
C    command-line variables are:
C
C    fin : name of input file; string
C
C    fout : name of output file; string
C
C    fnafgl : name of AFGL file; string
C
C    mnafgl : AFGL model number; integer.  The six models are:
C       1 = Sergio (testing)
C       2 = Midlatitude Summer
C       3 = Midlatitude Winter
C       4 = Subarctic Summer
C       5 = Subarctic Winter
C       6 = U.S. Standard 
C
C    nwant : number of desired gases; integer 1 to MXGAS, or -1=all.
C       If nwant=-1, all gases will be used and listg must not
C       be specified.  Otherwise listg must be specified and it must
C       contain nwant number of entries.
C
C    listg : list of integer gas IDs separated by a comma.
C       Alternately the IDs may be specified using a quoted string
C       containing integers separated by a blank space.  Examples:
C          listg=1,2,3,4,5
C          listg='1 2 3 4 5'
C       Note: due to the 80 char limit, the max number of entries
C       in listg is limited to between 25 to 37 gas IDs.  If you
C       want more gases than 25 and less than all of them, you'll
C       need to edit this routine to change the default listg.
C
C    toff : temperature offset to apply to entire profile; real
C
C    scaleh : scale height for Pres=10*e^(Alt/scaleh).  This
C        is only used if no surface altitude is supplied in the
C        input profile.
C
C    lpmax : extend profile to PMAX?; logical. true/false or T/F
C        If lpmax=true, the profile will be extended down to
C        PMAX (ie PLEV(1), the max pressure of the layering grid)
C        instead of stopping at whichever layer contains PSURF.
C
C    ldry : mixing ratios are for dry air?; logical.  true/false or T/F
C
C    lsvp : check water sat. vapor pres.?; logical.  true/false or T/F
C
C    lsplin : do spline interpolations?; logical.  true/false or T/F
C
C    Note: the meaning of "dry" and "wet" mixing ratios:
C       "wet"=the mixing ratios are exactly as given. The effects of
C           water on the make-up of the air has already been factored
C           into the supplied values.
C       "dry"=the mixing ratios are with respect to "dry air"; that is,
C           air without any water vapor present. The supplied values
C           should be adjusted for the presence of water vapor in the
C           actual profile.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C  2 Feb 2001 H.Motteler/S.Hannon Created; re-write of Howard Motteler's
C                                  kargs.f demo program.
C 21 Feb 2001 Scott Hannon      Changed error/warning messages from
C                                  I/O unit 6 to units IOERR & IOINFO.
C 27 Feb 2001 Scott Hannon      Added MASS and MASSW; correct bug in
C                                  LISTG assign when NWANT=-1
C 23 Mar 2001 Scott Hannon      Replace output var and command-line
C                                  argument PZMAX with SCALEH.
C  5 Sep 2001 Scott Hannon      Fix bug so MASS is real not integer
C 12 Sep 2001 Scott Hannon      Replace hardcoded default FNAFGL
C                                  with new parameter DFAFGL
C 21 Nov 2001 Scott Hannon      Add COMLEV for use in COMMNT; redo
C                                  COMMNT using VKLAYE & NAMGRD;
C                                  increase CJUNK & CJUNK from 30 to 40
C 27 Nov 2001 Scott Hannon      Add "lpmax" command line argument
C 05 Aug 2003 Scott Hannon      Change FIN & FOUT to CHAR*80 (was70)
C 23 Jun 2005 Scott Hannon      "trace" version for CO2,SO2,HNO3,N2O
C 16 Feb 2007 Scott Hannon      Updated check gas IDs for v2.05 clouds

!END====================================================================


C      =================================================================
       SUBROUTINE RDINFO(FIN, FOUT, COMMNT, FNAFGL, MNAFGL,
     $    NWANT, LISTG, MASSW, TOFF, SCALEH, LDRY, LSVP, LSPLIN, LPMAX)
C      =================================================================


C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE


C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
       include 'incLAY.f'


C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      From "util.f"
C      function LENNB = length of string excluding trailing blanks
C      function STR2BO = converts true/false string to boolean (LOGICAL)
C      subroutine UPCASE = converts a string to upper case


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input:
C      none
C
C      Output:
       CHARACTER*80 FIN
       CHARACTER*80 FOUT
       CHARACTER*256 COMMNT
       CHARACTER*70 FNAFGL
       INTEGER MNAFGL
       INTEGER  NWANT
       INTEGER  LISTG( MXGAS)
       REAL MASSW( MXGAS)
       REAL   TOFF
       REAL   SCALEH
       LOGICAL   LDRY
       LOGICAL   LSVP
       LOGICAL LSPLIN
       LOGICAL  LPMAX
C
C      Allowed gas IDs from "cbgids.f"
       COMMON /COMGID/ GIDS, GMASS
       INTEGER GIDS(NGIDS)
       REAL GMASS(NGIDS)

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER I
       INTEGER ID
       INTEGER IDCMAX
       INTEGER IG
       INTEGER IGJUNK(NGIDS+MXCLD)  ! junk gas id work array ..... THIS IS AN UPDATE on ../../klayersV205/Src_rtpV201/
       INTEGER IARGC
       INTEGER IWAT
       INTEGER J
       INTEGER K
       INTEGER LENNB
       INTEGER NARGS

       REAL MASS

       CHARACTER*80 BUF
       CHARACTER*80 VAL
       CHARACTER*80 VAR

       CHARACTER*1 CLDRY
       CHARACTER*1 CLSVP
       CHARACTER*1 CLSPLI
       CHARACTER*1 CLPMAX

       CHARACTER*40 CJUNK  ! junk work string
       CHARACTER*40 CJUNK2 ! 2nd junk work string

       LOGICAL STR2BO
       LOGICAL LNWANT
       LOGICAL LLISTG

C      Boundary pressure levels
       COMMON /COMLEV/ NAMGRD, PLEV
       CHARACTER*40 NAMGRD
       REAL PLEV(NBPLEV)

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
       FIN='klayers_in.rtp'         ! input filename
       FOUT='klayers_out.rtp'       ! output filename
       FNAFGL=DFAFGL                ! default AFGL file (see incLAY.f)
       MNAFGL=1                     ! AFGL model 1 = SERGIO TEST
       TOFF=0                       ! Temperature offset (usually 0)
       SCALEH=11.1                  ! scale height (km), see https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
       LDRY=.TRUE.                  ! Mixing ratios are for dry air?
       LSVP=.TRUE.                  ! Check water SVP?
       LSPLIN=.FALSE.               ! Do spline interpolation?
       LPMAX=.FALSE.                ! Extend profile to Pmax?
C
C      Default for desired output gases
C      Note: Set NWANT=-1 to use all gases
c      for earth
c       NWANT=8                      ! AIRS "trace" RTA needs 8 gases:
c       LISTG(1)=1                   !  1 = H2O
c       LISTG(2)=2                   !  2 = CO2
c       LISTG(3)=3                   !  3 = O3
c       LISTG(4)=4                   !  4 = N2O
c       LISTG(5)=5                   !  5 = CO
c       LISTG(6)=6                   !  6 = CH4
c       LISTG(7)=9                   !  9 = SO2
c       LISTG(8)=12                  ! 12 = HNO3

c      for Mars
       NWANT=6                      ! AIRS "trace" RTA needs 8 gases:
       LISTG(1)=1                   !  1 = H2O
       LISTG(2)=2                   !  2 = CO2
       LISTG(3)=3                   !  3 = N2
       LISTG(4)=4                   !  4 = O2
       LISTG(5)=5                   !  5 = CO
       LISTG(6)=6                   !  6 = NO
C
C      Comment string
C      Note: additional info will be tacked on to the end
       CJUNK=VKLAYE
       CJUNK2=NAMGRD
       I=LENNB(CJUNK)
       J=LENNB(CJUNK2)
       COMMNT='KLAYERS src=' // CJUNK(1:I) // '; grid=' // CJUNK2(1:J)

C      -----------------------------------------------------------------
C      Loop on program parameters
C      --------------------------
C      Determine the number of command-line arguments
       NARGS=IARGC()
C
C      Loop over the command-line arguments
       LNWANT=.FALSE.
       LLISTG=.FALSE.
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

             ELSEIF (VAR(1:6) .EQ. 'FNAFGL') THEN
                FNAFGL=VAL

             ELSEIF (VAR(1:6) .EQ. 'MNAFGL') THEN
                READ(VAL,*) MNAFGL

             ELSEIF (VAR(1:5) .EQ. 'NWANT') THEN
                LNWANT=.TRUE.
C               Make sure VAL is a 1 or 2 digit number
                K=LENNB(VAL)
                IF (VAL(1:1) .NE. '-' .AND.
     $              (ICHAR(VAL(1:1)) .LT. ICHAR('1') .OR.
     $               ICHAR(VAL(1:1)) .GT. ICHAR('9'))) K=3
                IF (K .EQ. 2 .AND. 
     $              (ICHAR(VAL(2:2)) .LT. ICHAR('0') .OR.
     $               ICHAR(VAL(2:2)) .GT. ICHAR('9'))) K=3
                IF (K .GT. 2) THEN
                   WRITE(IOERR,1010)
 1010              FORMAT('ERROR! bad NWANT, ',
     $             'either an unrecognized value or string too long')
                   STOP
                ENDIF
                K=NWANT
                READ(VAL,*) NWANT
                IF (NWANT .GT. MXGAS .OR.
     $          (NWANT .NE. -1 .AND. NWANT .LT. 1)) THEN
                   WRITE(IOERR,1011) NWANT, MXGAS
 1011              FORMAT('ERROR! NWANT=',I6,
     $             ', must have 1<=NWANT<=',I2,
     $             ', or NWANT=-1 (all)')
                   STOP
                ENDIF
                IF (LLISTG) THEN
                   IF (NWANT .EQ. -1) THEN
                      WRITE(IOERR,1012)
 1012                 FORMAT('ERROR! LISTG not allowed if NGAS=-1')
                      STOP
                   ELSEIF (K .NE. NWANT) THEN
                      WRITE(IOERR,1015) NWANT, K
 1015                 FORMAT('ERROR! NWANT=',I2,
     $                ' but LISTG specified ',I2,' gases.')
                      STOP
                   ENDIF
                ENDIF

             ELSEIF (VAR(1:5) .EQ. 'LISTG') THEN
                LLISTG=.TRUE.
C
C               Read all the gas IDs
                K=1
 10             READ(VAL,*,END=19) (IGJUNK(IG),IG=1,K)
                IF (K .GT. 25) THEN
                   WRITE(IOERR,1017)
 1017              FORMAT('ERROR! bad LISTG, ',
     $             'either an unrecognized value or more than 25 IDs')
                   STOP
                ELSEIF (K .GT. MXGAS) THEN
                   WRITE(IOERR,1018) MXGAS
 1018              FORMAT('ERROR! number of entries in LISTG ',
     $             'exceeds MXGAS=',I2)
                   STOP
                ENDIF
                K=K + 1  ! increment count of gases
                GOTO 10  ! loop to next entry
 19             CONTINUE
                K=K - 1  ! number of gases
                DO IG=1,K
                   LISTG(IG)=IGJUNK(IG)
                ENDDO
C
C               Check NWANT (if already set)
                IF (LNWANT) THEN
                   IF (NWANT .EQ. -1) THEN
                      WRITE(IOERR,1012)
                      STOP
                   ELSEIF (K .NE. NWANT) THEN
                      WRITE(IOERR,1015) NWANT, K
                      STOP
                   ENDIF
                ELSE
                   NWANT=K
                ENDIF

             ELSEIF (VAR(1:4) .EQ. 'TOFF') THEN
                READ(VAL,*) TOFF

             ELSEIF (VAR(1:6) .EQ. 'SCALEH') THEN
                READ(VAL,*) SCALEH

             ELSEIF (VAR(1:4) .EQ. 'LDRY') THEN
                LDRY=STR2BO(VAL)

             ELSEIF (VAR(1:4) .EQ. 'LSVP') THEN
                LSVP=STR2BO(VAL)

             ELSEIF (VAR(1:6) .EQ. 'LSPLIN') THEN
                LSPLIN=STR2BO(VAL)

             ELSEIF (VAR(1:5) .EQ. 'LPMAX') THEN
                LPMAX=STR2BO(VAL)

             ELSE
                WRITE(IOERR,1020) VAR
 1020           FORMAT('Unknown command-line argument: ',A6)
                STOP

             ENDIF

          ENDIF
       ENDDO  ! end of loop over command-line arguements
C      -----------------------------------------------------------------
C
       IF (NWANT .NE. -1 .AND. LNWANT .AND. .NOT. LLISTG) THEN
          WRITE(IOERR,1030) NWANT
 1030     FORMAT('ERROR! NWANT=',I2,' but LISTG not specifed')
          STOP
       ENDIF
C

C      -------------
C      Check gas IDs
C      -------------
       IF (NWANT .EQ. -1) THEN
C         Want all gases
          NWANT=MXGAS
          DO IG=1,NGIDS
             LISTG(IG)=GIDS(IG)
             MASSW(IG)=GMASS(IG)
          ENDDO
C         Also want all clouds
          DO IG=1,MXCLD
             LISTG(NGIDS+IG)=IDCLD1 + IG - 1
             MASSW(NGIDS+IG)=1.0
          ENDDO
C
       ELSE
C         Initialize IGJUNK count to zero
          DO IG=1,MXGAS
             IGJUNK(IG)=0
          ENDDO
C
C         Max cloud ID
          IDCMAX=IDCLD1 + MXCLD-1

C         Loop over LISTG
          DO IG=1,NWANT
C
C            current ID
             ID = LISTG(IG)

             IF (ID .GE. IDCLD1) THEN
                IF (ID .LE. IDCMAX) THEN
                   MASSW(IG)=1.0
                   I=NGIDS + ID - IDCLD1 + 1
                   IGJUNK(I)=IGJUNK(I) + 1
                ELSE
                   WRITE(IOERR,1042) IG, ID
                ENDIF
             ELSE
C
                CALL CHKGAS(ID,I,MASS)
C
                IF (I .GT. 0) THEN
                MASSW(IG)=MASS
                IGJUNK(I)=IGJUNK(I) + 1
                ELSE
                   WRITE(IOERR,1042) IG, ID
 1042              FORMAT('ERROR! LISTG(',I2,')=',I6,
     $             ' is not an allowed gas ID')
                   STOP
                ENDIF
             ENDIF
C
C            Check for repeats
             IF (IGJUNK(I) .GT. 1) THEN
                WRITE(IOERR,1045) ID
 1045           FORMAT('ERROR! gas ',I2,
     $          ' appears more than once in LISTG')
                STOP
             ENDIF
          ENDDO

       ENDIF
C
C      ---------------
C      Check for water
C      ---------------
       IF (LDRY) THEN
          IWAT=0
          DO I=1,NWANT
             IF (LISTG(I) .EQ. 1) IWAT=I
          ENDDO
          IF (IWAT .EQ. 0) THEN
             WRITE(IOERR,1050)
 1050        FORMAT('ERROR!, LDRY=TRUE but LISTG omits water (ID=1)!')
             STOP
          ENDIF
       ENDIF
C      

C      --------------------------
C      Add info to end of comment
C      --------------------------
C      Add TOFF to end of comment if not zero
       IF (ABS(TOFF) .GE. 0.01) THEN
          CJUNK='; TOFF='
          WRITE(CJUNK2,1060) TOFF
 1060     FORMAT(F6.2)
          I=LENNB(COMMNT)
          J=LENNB(CJUNK)
          K=LENNB(CJUNK2)
          IF (I .LE. 256-J-K) THEN
             COMMNT=COMMNT(1:I) // CJUNK(1:J) // CJUNK2(1:K)
          ENDIF
       ELSEIF (TOFF .NE. 0.0) THEN
          WRITE(IOINFO,1063) TOFF
 1063     FORMAT('WARNING! TOFF=',1PE10.3,' too small so reset to 0')
          TOFF=0.0
       ENDIF
C
C      Add AFGL model number to end of comment string
       CJUNK='; MNAFGL='
       WRITE(CJUNK2,1066) MNAFGL
 1066  FORMAT(I1)
       I=LENNB(COMMNT)
       J=LENNB(CJUNK)
       K=LENNB(CJUNK2)
       IF (I .LE. 256-J-K) THEN
          COMMNT=COMMNT(1:I) // CJUNK(1:J) // CJUNK2(1:K)
       ENDIF
C
C      Add state of logical vars to end of comment string
       CLDRY='F'
       CLSVP='F'
       CLSPLI='F'
       CLPMAX='F'
       IF (LDRY) CLDRY='T'
       IF (LSVP) CLSVP='T'
       IF (LSPLIN) CLSPLI='T'
       IF (LPMAX) CLPMAX='T'
       CJUNK='; LDRY=' // CLDRY // ', LSVP=' // CLSVP //
     $    ', LSPLIN=' // CLSPLI // ', LPMAX=' // CLPMAX
       I=LENNB(CJUNK)  ! length excluding trailing blanks
       J=LENNB(COMMNT)
       IF (J .LE. 256-I) THEN
          COMMNT=COMMNT(1:J) // CJUNK(1:I)
       ENDIF

       RETURN
       END
