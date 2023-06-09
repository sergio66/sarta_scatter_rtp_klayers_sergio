C=======================================================================
C=======================================================================
C
C    University of Maryland Baltimore Country (UMBC)
C
C    AIRS
C
C    SETEMS
C
C copied from /home/sergio/SARTA_CLOUDY/v107_Src_pclsam_ctype_100layer_rtpV201/SAFECOPY/setems.f
C removed NRHO,FRHO from argument list, and made them variables (as now NRHO = NEMIS, FRHO = FEMIS)
!F77====================================================================


!ROUTINE NAME:
C    SETEMS

!ABSTRACT:
C    Assign the emissivity and reflectivity data for every channel
C    being used.


!CALL PROTOCOL
C    RDEMIS ( NCHAN, NEMIS, FREQ, XEMIS, FEMIS, XRHO, 
C       LRHOT, EMIS, RHOSUN, RHOTHR )


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    INTEGER   NCHAN   number of channels          none
C    INTEGER   NEMIS   number of emis pts          none
C    REAL arr  FREQ    channel frequencies         cm^-1
C    REAL arr  XEMIS   raw emis points             none (0 to 1)
C    REAL arr  FEMIS   raw emis freq points        cm^-1
C    REAL arr  XRHO    raw reflec points           1/steradian
C    LOGICAL   LRHOT   force refl therm rho?       none

!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  EMIS    surface emissivity          none
C    REAL arr  RHOSUN  reflectivity for solar      1/steradian
C    REAL arr  RHOTHR  reflectivity for thermal    1/steradian


!INPUT/OUTPUT PARAMETERS:
C    none


!RETURN VALUES:
C    none


!PARENT(S):
C    USEFAST


!ROUTINES CALLED:
C    none


!FILES ACCESSED:
C    incFTC.f : include file of parameter statements accessed during
C       compilation only.


!COMMON BLOCKS
C    none


!DESCRIPTION:
C    August 2000 version of the 100 layer AIRS Fast Transmittance
C    Code by L.Strow/S.Hannon.
C
C    Loops over the passed in EMIS & RHO points and interpolates
C    or extrapolates them onto the channel freqs being used.
C    Interpolations are linear in wavelength (not frequency).


!ALGORITHM REFERENCES:
C    none


!KNOWN BUGS AND LIMITATIONS:
C    none


!ROUTINE HISTORY:
C    Date        Programmer     Comments
C    ----------- -------------- ----------------------------------------
C    14 Feb 2001 Scott Hannon   Created based on rdemis.f
C    06 Feb 2004 Scott Hannon   Add LRHOT argument and associated code


!END====================================================================

C      =================================================================
       SUBROUTINE SETEMS ( NCHAN, NEMIS, FREQ, XEMIS, FEMIS,
     $    XRHO, LRHOT, EMIS, RHOSUN, RHOTHR )
C      =================================================================


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
C      none


C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input
       INTEGER  NCHAN
       INTEGER  NEMIS
       REAL   FREQ(MXCHAN)
       REAL  XEMIS(MXEMIS)
       REAL  FEMIS(MXEMIS)
       REAL   XRHO(MXEMIS)
       LOGICAL  LRHOT
C
C      Output
       REAL   EMIS(MXCHAN)
       REAL RHOSUN(MXCHAN)
       REAL RHOTHR(MXCHAN)


C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       INTEGER   NRHO
       REAL   FRHO(MXEMIS)
       INTEGER      I  ! generic looping variable
       INTEGER      J  ! generic looping variable
       INTEGER SORTED  ! conrol flag for sorts
       REAL DELELM(MXEMIS)  ! delta for emis wavelength
       REAL DELEMS(MXEMIS)  ! delta for emis
       REAL DELRHO(MXEMIS)  ! delta for rho
       REAL DELRLM(MXEMIS)  ! delta for rho wavelength
       REAL     DX          ! delta X
       REAL   ELAM(MXEMIS)  ! emis wavelength
       REAL   RLAM(MXEMIS)  ! rho wavelength
       REAL RJUNK1  ! generic junk/work variable
       REAL RJUNK2  ! generic junk/work variable


C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none


C***********************************************************************
C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C***********************************************************************
       NRHO = NEMIS
       DO I = 1,NEMIS
         FRHO(I) = FEMIS(I)
       END DO

C
C      ---------------------------
C      Sort the emis & rho by freq
C      ---------------------------
C      Emissivity
       SORTED=1
 10    IF (SORTED .EQ. 1) THEN
          SORTED=0
          DO I=1,NEMIS-1
             IF (FEMIS(I) .GT. FEMIS(I+1)) THEN
                RJUNK1=FEMIS(I)
                RJUNK2=XEMIS(I)
                FEMIS(I)=FEMIS(I+1)
                XEMIS(I)=XEMIS(I+1)
                FEMIS(I+1)=RJUNK1
                XEMIS(I+1)=RJUNK2
                SORTED=1
             ENDIF
          ENDDO
          GOTO 10
       ENDIF
C
C      Reflectivity
       SORTED=1
 11    IF (SORTED .EQ. 1) THEN
          SORTED=0
          DO I=1,NRHO-1
             IF (FRHO(I) .GT. FRHO(I+1)) THEN
                RJUNK1=FRHO(I)
                RJUNK2=XRHO(I)
                FRHO(I)=FRHO(I+1)
                XRHO(I)=XRHO(I+1)
                FRHO(I+1)=RJUNK1
                XRHO(I+1)=RJUNK2
                SORTED=1
             ENDIF
          ENDDO
          GOTO 11
       ENDIF
C

C      ----------------------------------------
C      Calc wavelength and interpolation deltas
C      ----------------------------------------
C      Raw emissivity points
       ELAM(1)=1.0E+4/FEMIS(1) ! wavelength in microns
       DO J=2,NEMIS
          ELAM(J)=1.0E+4/FEMIS(J)
          DELELM(J-1)=ELAM(J) - ELAM(J-1)
          DELEMS(J-1)=XEMIS(J) - XEMIS(J-1)
       ENDDO
C
C      Raw reflectivity points
       RLAM(1)=1.0E+4/FRHO(1) ! wavelength in microns
       DO J=2,NRHO
          RLAM(J)=1.0E+4/FRHO(J)
          DELRLM(J-1)=RLAM(J) - RLAM(J-1)
          DELRHO(J-1)=XRHO(J) - XRHO(J-1)
       ENDDO
C

C      ------------------------------------
C      Interpolate emis & rho onto channels
C      ------------------------------------
C      Loop over channels (note: chan freqs may be in any order)
       DO I=1,NCHAN
C
C         Emissivity
          IF ( FREQ(I) .LE. FEMIS(1) ) THEN
             EMIS(I)=XEMIS(1)
          ELSEIF ( FREQ(I) .GE. FEMIS(NEMIS) ) THEN
             EMIS(I)=XEMIS(NEMIS)
          ELSE
C            Determine the index of the upper bounding FEMIS
             J=2
 20          IF ( FEMIS(J) .LT. FREQ(I) ) THEN
                J=J + 1
                GOTO 20
             ENDIF
C            Convert to lower boundary index
             J=J - 1
C
             DX=(( 1.0E+4/FREQ(I) - ELAM(J) )/DELELM(J))
             EMIS(I) = XEMIS(J) + DX*DELEMS(J)
          ENDIF
C
C         Reflectivity
          IF ( FREQ(I) .LE. FRHO(1) ) THEN
             RHOSUN(I)=XRHO(1)
          ELSEIF ( FREQ(I) .GE. FRHO(NRHO) ) THEN
             RHOSUN(I)=XRHO(NRHO)
          ELSE
C            Determine the index of the upper bounding FRHO
             J=2
 21          IF ( FRHO(J) .LT. FREQ(I) ) THEN
                J=J + 1
                GOTO 21
             ENDIF
C            Convert to lower boundary index
             J=J - 1
C
             DX=(( 1.0E+4/FREQ(I) - RLAM(J) )/DELRLM(J))
             RHOSUN(I) = XRHO(J) + DX*DELRHO(J)
          ENDIF
C
          IF (LRHOT) THEN
C            Force reflected thermal rho = (1-emis)/pi
             RHOTHR(I)=(1.0 - EMIS(I))/PI
          ELSE
C            Reflected thermal uses the same rho as solar
             RHOTHR(I)=RHOSUN(I)
          ENDIF
C
       ENDDO

       RETURN
       END
