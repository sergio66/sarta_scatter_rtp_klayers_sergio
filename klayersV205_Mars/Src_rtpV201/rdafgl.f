C***********************************************************************
C
C  SUBROUTINE     RDAFGL        
C
C  PURPOSE        TO READ AFGL MODEL PROFILES
C
C  VERSION        1.0 April 1995
C                 based on version 3.y of LAYERS routine PROFIL
C                 by D.P. EDWARDS   10/05/93
C 
C  I/O UNITS      IND    I/P       AFGL MODEL PROFILE DATA SET
C
C Modified: 17 July 2000 by Scot Hannon; adjust mixing ratios
C of CFC-11,12,13 and CCl4.
C
C Revised: 9 October 2000 by Scott to fix a bug introduced by the
C last modification which mult'ed mixing ratio of last desired gas
C by 1.9*2.3*17*0.74=54.9746 if not all gases were used.
C
C Modified: 24 July 2001 by Scott Hannon.  Moved the adjustments
C from the previous two mods to new routine ADAFGL.
C
C Update: 16 Feb 2007 Scott Hannon - add dummy cloud profile data
C***********************************************************************
C
       SUBROUTINE RDAFGL(FNAFGL, IND, NMODEL, NWANT, IDWANT, NLEVAF,
     $    NAFGL, IDAFGL, LATAF, HAF, PAF, TAF, DAF, MIXAF)
C
C-----------------------------------------------------------------------
C
      include 'incLAY.f'
C
C      Input/Output Parameters
       INTEGER    IND               ! AFGL file I/O unit number
       INTEGER IDWANT( MXGAS)       ! gas IDs of wanted gases
       INTEGER  NWANT               ! number of gases wanted
       INTEGER NLEVAF               ! AFGL number of levels
       INTEGER  NAFGL               ! number of AFGL gas profiles
       INTEGER IDAFGL( MXGAS)       ! AFGL gas IDs
       INTEGER NMODEL               ! AFGL model number
       REAL  LATAF                  ! AFGL profile latitude
       REAL  MIXAF(  MXIN, MXGAS)   ! AFGL profile mixing ratio
       REAL    HAF(  MXIN)          ! AFGL profile height
       REAL    PAF(  MXIN)          ! AFGL profile pressure
       REAL    TAF(  MXIN)          ! AFGL profile temperature
       REAL    DAF(  MXIN)          ! AFGL profile density
       CHARACTER*70 FNAFGL          ! AFGL filename
C
C      Local variables
       INTEGER      I
       INTEGER   IGAS
       INTEGER      L
       REAL    XAF(  MXIN)
       CHARACTER*1 COMMNT
       CHARACTER*1 WORDIN
       CHARACTER*7    KEY
       CHARACTER*7  MODEL
       LOGICAL   SKIP

C-----------------------------------------------------------------------
C
       COMMNT='!'
C
C      Open AFGL file
       OPEN(IND,FILE=FNAFGL,FORM='FORMATTED',STATUS='OLD')
C
C      Read the number of levels in the model profiles
   20  READ(IND,1000) WORDIN
       IF (WORDIN .EQ. COMMNT) GOTO 20
       IF (WORDIN .NE. COMMNT) BACKSPACE(IND)
       READ(IND,*) NLEVAF
C
C      Load "MODEL" variable with name of string to look for
       WRITE(MODEL,1010) NMODEL
 1010  FORMAT('MODEL',I2)
C
C      Read model number; skip until reach desired model
   30  READ(IND,1000) WORDIN
       IF (WORDIN .EQ. COMMNT) GOTO 30
       IF (WORDIN .NE. COMMNT) BACKSPACE(IND)
   40  READ(IND,1020) KEY
       IF (KEY .NE. MODEL) GOTO 40

C      READ IN CONSTITUENT GAS DENSITIES FROM AFGL DATA 
C      Height[km], Pres[mb], Temp[K], AirDen[cm-3],
C      GasDen[ppmv], in order of ID:
C      1(H2O), 2(CO2), 3(O3), 4(N2O), 5(CO), 6(CH4), 7(O2) etc
C
C      Read profile latitude
   50  READ(IND,1000) WORDIN
       IF (WORDIN .EQ. COMMNT) GOTO 50
       IF (WORDIN .NE. COMMNT) BACKSPACE(IND)
       READ(IND,*) LATAF
C
C      Read level altitudes
   60  READ(IND,1000) WORDIN
       IF (WORDIN .EQ. COMMNT) GOTO 60
       IF (WORDIN .NE. COMMNT) BACKSPACE(IND)
       READ(IND,*) (HAF(I),I=1,NLEVAF)
C
C      Read pressure levels
   70  READ(IND,1000) WORDIN
       IF (WORDIN .EQ. COMMNT) GOTO 70
       IF (WORDIN .NE. COMMNT) BACKSPACE(IND)
       READ(IND,*) (PAF(I),I=1,NLEVAF)
C
C      Read level temperatures
   80  READ(IND,1000) WORDIN
       IF (WORDIN .EQ. COMMNT) GOTO 80
       IF (WORDIN .NE. COMMNT) BACKSPACE(IND)
       READ(IND,*) (TAF(I),I=1,NLEVAF)
C
C      Read level air densities
   90  READ(IND,1000) WORDIN
       IF (WORDIN .EQ. COMMNT) GOTO 90
       IF (WORDIN .NE. COMMNT) BACKSPACE(IND)
       READ(IND,*) (DAF(I),I=1,NLEVAF)
C
C      ---------------------
C      Read model gas mixing ratios
       NAFGL=0
  100  READ(IND,1000) WORDIN
       IF (WORDIN(1:1) .EQ. COMMNT) GOTO 100
       IF (WORDIN(1:1) .NE. COMMNT) BACKSPACE(IND)
       READ(IND,1020) KEY
       IF (KEY .EQ. 'MODEND ') GOTO 110
C
       BACKSPACE(IND)
       READ(IND,*) IGAS
C
       SKIP=.TRUE.
       DO I=1,NWANT
          IF (IGAS .EQ. IDWANT(I)) THEN
             SKIP=.FALSE.
             NAFGL=NAFGL+1
             IDAFGL(NAFGL)=IGAS
          ENDIF
       ENDDO
       IF (SKIP) THEN
          READ(IND,*) (XAF(I),I=1,NLEVAF)
       ELSE
          READ(IND,*) (MIXAF(I,NAFGL),I=1,NLEVAF)
C
C         Adjust AFGL mixing ratios if necessary
          CALL ADAFGL(NAFGL,IDAFGL,NLEVAF,MIXAF)
C
       ENDIF
C
       GOTO 100
C      --------
C
 110   CONTINUE

CCCC UNCOMMENT C      READ MINOR GAS CONCENTRATIONS 
CCCC UNCOMMENT C      Gas Den[ppmv] for minor gases in order of ID
CCCC UNCOMMENT C      from 8(NO), 9(SO2), etc through to 28(PH3) etc
CCCC UNCOMMENT C
CCCC UNCOMMENT   110  READ(IND,1020) KEY
CCCC UNCOMMENT        IF (KEY .NE. 'MINGAS ') GOTO 110 
CCCC UNCOMMENT C
CCCC UNCOMMENT C      -----------------------------
CCCC UNCOMMENT C      Read in minor gas mixng ratio     
CCCC UNCOMMENT   120  READ(IND,1000) WORDIN
CCCC UNCOMMENT        IF (WORDIN .EQ. COMMNT) GOTO 120
CCCC UNCOMMENT        IF (WORDIN .NE. COMMNT) BACKSPACE(IND)
CCCC UNCOMMENT        READ(IND,1020) KEY
CCCC UNCOMMENT        IF (KEY .EQ. 'DATEND ') GOTO 130 
CCCC UNCOMMENT C
CCCC UNCOMMENT        BACKSPACE(IND)     
CCCC UNCOMMENT        READ(IND,*) IGAS
CCCC UNCOMMENT C
CCCC UNCOMMENT        SKIP=.TRUE.
CCCC UNCOMMENT        DO I=1,NWANT
CCCC UNCOMMENT           IF (IGAS .EQ. IDWANT(I)) THEN
CCCC UNCOMMENT              SKIP=.FALSE.
CCCC UNCOMMENT              NAFGL=NAFGL+1
CCCC UNCOMMENT              IDAFGL(NAFGL)=IGAS
CCCC UNCOMMENT           ENDIF
CCCC UNCOMMENT        ENDDO
CCCC UNCOMMENT        IF (SKIP) THEN
CCCC UNCOMMENT           READ(IND,*) (XAF(I),I=1,NLEVAF)
CCCC UNCOMMENT        ELSE
CCCC UNCOMMENT           READ(IND,*) (MIXAF(I,NAFGL),I=1,NLEVAF)
CCCC UNCOMMENT C
CCCC UNCOMMENT C         Adjust AFGL mixing ratios if necessary
CCCC UNCOMMENT           CALL ADAFGL(NAFGL,IDAFGL,NLEVAF,MIXAF)
CCCC UNCOMMENT C
CCCC UNCOMMENT        ENDIF
CCCC UNCOMMENT C
CCCC UNCOMMENT        GOTO 120
CCCC UNCOMMENT C      --------
CCCC UNCOMMENT C

C      Close the AFGL file
  130  CLOSE(IND)
C
 1020  FORMAT(A7) 
 1000  FORMAT(A1)
C
C
C      Load dummy cloud data
       DO I=1,NWANT
          IF (IDWANT(I) .GE. IDCLD1) THEN
             NAFGL=NAFGL+1
             IDAFGL(NAFGL)=IDWANT(I)
             DO L=1,NLEVAF
                MIXAF(L,NAFGL)=0.0
             ENDDO
          ENDIF
       ENDDO
C
       RETURN
       END
