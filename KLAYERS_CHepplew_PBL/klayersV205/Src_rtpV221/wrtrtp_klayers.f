C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              WRTRTP
C
!F77====================================================================


!ROUTINE NAME: WRTRTP


!ABSTRACT:
C    Write a profile to a previously openned RTP file


!CALL PROTOCOL:
C    WRTRTP(IP, ISTAT, IOPCO, LAYBOT, NGASES, LISTG, GORDER, HBOUND,
C       PLEV, PLAY, TLAY, ALAY, TOFF, PROF)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   IP      profile count so far        none
C    INTEGER   IOPCO   input RTP file I/O number   none
C    INTEGER   LAYBOT  bottom layer number         none
C    INTEGER   NGASES  # of gases                  none
C    INT arr   LISTGF  list of gases               none (HITRAN IDs)
C    INT arr   GORDER  order of gases in ALAY      none
C    REAL arr  HBOUND  layer boundary altitudes    meters
C    REAL arr  PLEV    layer boundary pressures    millibars
C    REAL arr  PLAY    layer average pressure      atmospheres
C    REAL arr  TLAY    layer average air temp.     Kelvin
C    REAL arr  ALAY    layer column densities      kilomoles/cm^2
ccc
cC    REAL arr  ALAY    layer column densities      molecules/cm^2
ccc
C    STRUCT    PROF    RTP profile structure       (see attributes)


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   ISTAT   I/O error status            none


!INPUT/OUTPUT PARAMETERS: none


!RETURN VALUES: none


!PARENT(S): KLAYERS


!ROUTINES CALLED: none


!FILES ACCESSED:
C    Output RTP file with I/O number IOPCO
C    unit 6: standard output (error message)


!COMMON BLOCKS: none


!DESCRIPTION:
C    Writes a single profile to a previously openned RTP file.


!ALGORITHM REFERENCES: see DESCRIPTION


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C 31 Jan 2001 Scott Hannon      created
C  9 Feb 2001 Scott Hannon      output data top-down
C  9 Mar 2001 Scott Hannon      add if/elseif for GUCOUT 1 & 2
C 23 Mar 2001 Scott Hannon      add LAYBOT to input vars; add local vars
C                                  NLAY & NLEV; change loops 1 to
C                                  NPBLEV into 1 to NLEV; change loops
C                                  1 to MYNLAY into 1 to NLAY.
C 16 Feb 2007 Scott Hannon      Add AMULT and update for v2.05 clouds


!END====================================================================

C      =================================================================
       SUBROUTINE WRTRTP(IP, ISTAT, IOPCO, LAYBOT, NGASES, LISTG,
     $    GORDER, HBOUND, PLEV, PLAY, TLAY, ALAY, TOFF, PROF)
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
       INTEGER     IP
       INTEGER  IOPCO
       INTEGER LAYBOT
       INTEGER NGASES
       INTEGER  LISTG( MXGAS)
       INTEGER GORDER( MXGAS)
       REAL HBOUND(MYNLAY+1)
       REAL   PLEV(NBPLEV)
       REAL   PLAY(MYNLAY)
       REAL   TLAY(MYNLAY)
       REAL   ALAY(MYNLAY, MXGAS)
       REAL   TOFF
C      Profile data structure
       RECORD /RTPPROF/ PROF
C
C      Output parameters:
       INTEGER ISTAT


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER IG
       INTEGER IL
       INTEGER ILR
       INTEGER NLAY
       INTEGER NLEV
       INTEGER rtpwrite
       REAL AMULT


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************

C      Number of layers and levels used by this profile
       NLAY=MYNLAY - LAYBOT + 1
       NLEV=NLAY + 1
C      Note: data in all "REAL" input vars is in sorted bottom
C      of atmos to top, with bottom of atmos at index LAYBOT.

C      -------------------------
C      Load up the new PROF data
C      -------------------------
       PROF.nlevs=NLEV
C
C      Loop over the boundary levels
       DO IL=1,NLEV
          ILR=NBPLEV + 1 - IL  ! reverse order of loop index
          PROF.plevs(IL)=PLEV(ILR)
          PROF.palts(IL)=HBOUND(ILR)
       ENDDO
C
C      Loop over the layers
       DO IL=1,NLAY
          ILR=MYNLAY + 1 - IL  ! reverse order of loop index
ccc Removed for rtpV201
c          PROF.plays(IL)=PLAY(ILR)*1013.25  ! convert atm to mb
ccc
          PROF.ptemp(IL)=TLAY(ILR) + TOFF
       ENDDO
C
C      Load up layer integrated gas amounts
       IF (GUCCLD .NE. 5) THEN
C         This routine only works with GUCCLD=5 (grams per meter^2)
          WRITE(IOERR,1000) GUCCLD
 1000     FORMAT('Error! GUCCLD=',I4,' but wrtrtp hardcoded for 5')
       ENDIF
       IF (GUCOUT .EQ. 1) THEN
C         Want molecules/cm^2; convert from kilomoles/cm^2
C         Note: this is the unit used by AIRS level2 processing
          DO IG=1,NGASES
             IF (LISTG(IG) .GE. IDCLD1) THEN
C               Convert kilomoles/cm^2 to g/m^2 for 1.0 g/mole
                AMULT=1.0E+7
             ELSE
                AMULT=6.02214199E+26
             ENDIF
             DO IL=1,NLAY
                ILR=MYNLAY + 1 - IL  ! reverse order of loop index
                PROF.gamnt(IL,IG)=ALAY(ILR,GORDER(IG))*AMULT
             ENDDO
          ENDDO

       ELSEIF (GUCOUT .EQ. 2) THEN
C         Want kilomoles/cm^2
C         Note: this is the unit used internally by KLAYERS
          DO IG=1,NGASES
             IF (LISTG(IG) .GE. IDCLD1) THEN
                AMULT=1E+7
             ELSE
                AMULT=1.0
             ENDIF
             DO IL=1,NLAY
                ILR=MYNLAY + 1 - IL  ! reverse order of loop index
                PROF.gamnt(IL,IG)=ALAY(ILR,GORDER(IG))*AMULT
             ENDDO
          ENDDO
       ELSE
          WRITE(IOERR,1005) GUCOUT
 1005     FORMAT('ERROR! wrtrtp can not process GUCOUT=',I4)
          STOP
       ENDIF
C

C      -------------------------
C      Write the current profile
C      -------------------------
       ISTAT=rtpwrite(IOPCO, PROF)
C
       IF (ISTAT .EQ. -1) THEN
          WRITE(IOERR,1010) IP
 1010     FORMAT('ERROR! unable to write PROF data for prof ',I5)
ccc
c Do not stop, let klayers_rtp handle ISTAT=-1 error
c          STOP
ccc
       ENDIF
C
cccc for debugging ccccccccccc
c       ILR=MIN(NGASES,6)
c       DO IL=1,NLAY
c          WRITE(6,2010) IL, PROF.palts(IL), PROF.plevs(IL),
cccc
c     $       PROF.ptemp(IL),(PROF.gamnt(IL,IG),IG=1,ILR)
c 2010     FORMAT(I4,9(1X,1PE11.4))
cccc Replaced with above for rtpV201
cc     $       PROF.plays(IL), PROF.ptemp(IL),
cc     $       (PROF.gamnt(IL,IG),IG=1,ILR)
cc 2010     FORMAT(I4,10(1X,1PE11.4))
cccc
c       ENDDO
ccccccccccccccc
C
       RETURN
       END
