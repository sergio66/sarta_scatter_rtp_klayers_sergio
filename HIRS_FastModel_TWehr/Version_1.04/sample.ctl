!=======================================================================
!
!    Joint Center for Earth Systems Technology (JCET)
!    University of Maryland Baltimore County   (UMBC)
!
!    HIRS Fast Radiative Transfer Program
!
!    Version number see HFFPVS in file hffp_glob_dec.f
!
!                                   Tobias Wehr       wehr@umbc.edu
!                                   L.Larrabee Strow  strow@umbc.edu
!
!=======================================================================
!
! sample.ctl is a template control file for running the program
! 
! The program will read input informations from this control file.
! It looks for certain keywords, followed by parameters. The 
! parameters can be real numbers, integer numbers or strings.
! Strings are always encapsulated with ''. 
!
! The order of the keywords must not be changed.
!
! Comment lines in a control file are allowed only before a
! keyword line. Comment lines must have a ! symbol in the first
! column. Blank lines are not allowed. At the end of a control
! file must be a carriage return.
!
! ...............................................................
!
! valid satellite numbers are 7,9,10,11,12,14,15
SATELLITE NUMBER
15
!
TOTAL NUMBER OF CHANNELS
19
! 
! For each channel the channel ID is required, in ascending
! order. If the total number of channels is N, then N lines
! must be given here with N (different) channel IDs
ID NUMBERS OF CHANNELS
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
!
! the fixed-transmittance tuning parameters for each channel
! set to 1.0 for no tuning
TAU(FIXED) TUNING PARAMETERS
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
!
! the water-transmittance tuning parameter for each channel
! set to 1.0 for no tuning
TAU(H2O) TUNING PARAMETERS
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
!
! the ozone-transmittance tuning parameters for each channel
! set to 1.0 for no tuning
TAU(O3) TUNING PARAMETERS
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
1.0
!
! the radiance tuning parameter for each channel. This is 
! an offset! Set to 0.0 for no tuning
! NOTE: we are tuning the radiances, not the brightness
!       temperatures!
OFFSET TUNING PARAMETER
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
!
! channel dependend zenith angle offset (tuning parameter)
CHANNEL ZENITH ANGLE OFFSET
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
0.0
!
! the offset tuning parameters tunes either the radiances
! or the brightness temperatures
OFFSET TUNES RADIANCES (1) OR T_BRIGHT (2)
1
!
! set this parameter to zero!
ADJUST TAU TUNING PARAMETERS (0=no, 1=yes)
0
!
! calculate brightness temperatures (must be 1 if 
! offset tunes T_bright)
CALCULATE T_BRIGHT (0=NO, 1=YES)
0
!
! calculate Jacobians for temperature (T), water (W)
! and/or ozone (O): 0=don't calculate Jacobians, 
! 1=calculate analytic Jacobians, 2=calculate finite diff.J.
CALCULATE JACOBIANS FOR TWO (0=NO,1=AJ,2=FD,3=AJ&FD)
3 3 3
!
! Jacobians might be calculated as derivatives of either 
! radiances or brightness temperatures. If no Jacobians 
! are calculated this is a dummy.
JACOBIANS BE DERIVATIVES OF RAD (1) OR T_BRIGHT (2)
1
!
! output will be written to file only if requested!
WRITE OUTPUT TO FILE (1=YES, 0=NO)
1
!
! ============================================================
! The next part is the definition of the atmospheres for which
! the radiances shall be calculated
! ============================================================
!
! number of atmospheres you want to calculate
NUMBER OF ATMOSPHERES
2
!
! Below, for every atmosphere one parameter line has to be
! given
!
ATMOSPHERE PROFILE FILE NAMES
'./INPUT/prof1'
'./INPUT/prof2'
SURFACE PRESSURES (unit mbar)
1100.0
1100.0
!
! If the surface temperature is set to <0.0 the program
! will use the temperature of the surface layer as 
! the temperature of the surface
SURFACE TEMPERATURES (unit K)
-1.0
-1.0
!
! The next gives the observation angle. We have two options
! here: either to give the zenith angle in degree or the
! secant of the angle. Depending on the, the keyword has
! to be either
! OBSERVATION ZENITH ANGLES
! or
! OBSERVATION ZENITH SECANTS
! followed by the parameter lines.
! In this sample we use the secant
OBSERVATION ZENITH SECANTS
1.0
1.0
!
! The next are the sun zenith angles. Also two options
! here: either angles in degree or in secants. The respective
! keywords are
! SUN ZENITH SECANTS
! or
! SUN ZENITH ANGLES
SUN ZENITH ANGLES
0.0
0.0
SUN SOLID ANGLE
6.7851E-05
6.7851E-05
SURFACE EMISSIVITIES
0.975
0.975
!
! The control file must end with a comment line or carr.return
