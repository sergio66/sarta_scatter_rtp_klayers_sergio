C=======================================================================
C=======================================================================
C
C              University of Maryland Baltimore County [UMBC]
C
C              AIRS
C
C              INTEG
C
!F77====================================================================

!ROUTINE NAME: INTEG

!ABSTRACT:
C    Integrates AIRS sub-layers to determine the 100 AIRS layers
C    profile values for use with the Fast Transmittance Code.

!CALL PROTOCOL:
C    INTEG(LTAIR, LAYBOT, NSUB, NFINE, NGASES,
C          PFINE, PSUB, TSUB, TGSUB, MRSUB,
C          DZSUB, PLAY,  TLAY,  TGLAY, ALAY, DZLAY)

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  DZSUB   thickness of each sublayer  m
C    INTEGER   LAYBOT  bottom layer number
C    LOGICAL   LTAIR   Use TSUB for ALAY calc?     none
C    INTEGER   NFINE   number of fine sublevels    none
C    INTEGER   NGASES  number of gases             none
C    INT arr   NSUB    # of sub-levels each layer  none
C    REAL arr  MRSUB   avg mix ratio each layer    ppmv
C    REAL arr  PFINE   fine sub-level pressures    mb
C    REAL arr  PSUB    avg pres of each sub-layer  mb
C    REAL arr  TGSUB   avg gas temp of sub-layers  K
C    REAL arr  TSUB    avg air temp of sub-layers  K

!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL arr  ALAY    avg gas amount, AIRS layer  k.mol/cm^2
C    READ arr  DZLAY   layer thickness             m
C    REAL arr  PLAY    avg air pressure, AIRS lay  atm
C    REAL arr  TGLAY   avg gas temp, AIRS layers   K
C    REAL arr  TLAY    avg air temp, AIRS layers   K

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S): KLAYERS

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C    Integrates (sums) the sub-layers for each layer using averaged
C    sub-layer values for temperature, mixing ratio, and pressure.
C
C    The ideal gas law (PV=nRT) is used to compute "air" density (a
C    quite reasonable approximation for atmospheric pressures) and
C    then multiplied by the pathlength and mixing ratio to determine
C    absorber "amount". 
C
C    An average layer temperature is computed using the "air amount"
C    in the sub-layers as a weight.
C
C    An average layer pressure is computed using only the profile
C    independent 100 AIRS layer pressure boundaries. The equation 
C    used, delta_P/delta_ln(P), comes from assuming P = B*exp(A*z)
C    and letting Pavg = integral{P*dz}/integral{dz}.

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS:
C    No error checking of the input profile 

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C Mar  8 1995 Scott Hannon/UMBC created
C Mar 28 1995 Scott Hannon/UMBC added layer thickness DZLAY
C Jun 23 1995 Scott Hannon      Correct some comments
C Oct  9 2000 Scott Hannon      Change loops from MXGAS to NGASES
C 23 Mar 2001 Scott Hannon      Add LAYBOT to input vars
C  9 May 2001 Scott Hannon      Change loop 1 to MYNLAY into LAYBOT
C                                  to MYNLAY; I forgot to implement
C                                  this change 23 Mar as intended; this
C                                  bug caused problems only when
C                                  looping over multiple profiles)

!END====================================================================


C      =================================================================
       SUBROUTINE INTEG(LTAIR, LAYBOT, NSUB, NFINE, NGASES, 
     $    PFINE, PSUB, TSUB,
     $    TGSUB, MRSUB, DZSUB, PLAY, TLAY, TGLAY, ALAY, DZLAY)
C      =================================================================
C
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
C      none

C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------
C      Input parameters:
       INTEGER NSUB(MYNLAY), NFINE, NGASES
       INTEGER LAYBOT
       LOGICAL LTAIR
       REAL PFINE(NSUBLV), PSUB(NSUBLV-1), TSUB(NSUBLV-1),
     $    MRSUB(NSUBLV-1,MXGAS), DZSUB(NSUBLV-1), TGSUB(NSUBLV-1,MXGAS)
C
C      Output parameters:
       REAL PLAY(MYNLAY), TLAY(MYNLAY), TGLAY(MYNLAY,MXGAS),
     $    ALAY(MYNLAY,MXGAS), DZLAY(MYNLAY)

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
       REAL ASUM(MXGAS),TSUM,DENAV,RJUNK,RJUNK2,AJUNK,TGSUM(MXGAS),
     $    DZSUM
       INTEGER I,J,K,L,M

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************
C
C      DENAV=density_ref particles/(atm.cm^3) * 10E-6/ppmv * 100 cm/m *
C         1/(N_avag*1000) k.mol/particles * atm/(1013.25 mb) * Tref (K)
C      with:
C         density_ref=2.6867E+19 = N_avag * ( n/V = P/(R*T) ) at STP
C         N_avag=6.022045E+23
C         Tref=273.15
C
       DENAV=1.2027E-12
C
C      -----------------------------
C      Loop over the 100 AIRS layers
C      -----------------------------
       J=0
       DO L=LAYBOT,MYNLAY
C
          DO M=1,NGASES
             ASUM(M)=0.0
             TGSUM(M)=0.0
          ENDDO
          TSUM=0.0
          DZSUM=0.0
          AJUNK=0.0
C
C         ------------------------
C         Loop over the sub-layers
C         ------------------------
          DO I=1,NSUB(L)
             K=I+J
C
             DZSUM=DZSUM+DZSUB(K)
C
C            Weight the temperature sums by the amount
             RJUNK=DENAV*DZSUB(K)*PSUB(K)
             AJUNK=AJUNK + (RJUNK/TSUB(K))
             TSUM=TSUM+RJUNK
             IF (LTAIR) THEN
c                DO M=1,MXGAS
                DO M=1,NGASES
                   RJUNK2=RJUNK*MRSUB(K,M)
                   ASUM(M)=ASUM(M) + (RJUNK2/TSUB(K))
                   TGSUM(M)=TGSUM(M)+RJUNK2
                ENDDO
             ELSE
c                DO M=1,MXGAS
                DO M=1,NGASES
                   RJUNK2=RJUNK*MRSUB(K,M)
                   ASUM(M)=ASUM(M) + (RJUNK2/TGSUB(K,M))
                   TGSUM(M)=TGSUM(M)+RJUNK2
                ENDDO
             ENDIF
C
          ENDDO
C
          TLAY(L)=TSUM/AJUNK
          DZLAY(L)=DZSUM
          DO M=1,NGASES
             TGLAY(L,M)=TGSUM(M)/ASUM(M)
             ALAY(L,M)=ASUM(M)
          ENDDO
C
C         Because the FTC code requires constant layer pressure, the
C         average layer pressure must be profile independent (ie it
C         only depends upon the profile layer boundary pressures).
          PLAY(L)=( PFINE(J+1+NSUB(L))-PFINE(J+1) )/
     +         (1013.25*LOG( PFINE(J+1+NSUB(L))/PFINE(J+1) ))
C
          J=J+NSUB(L)
       ENDDO
C
       RETURN
       END
