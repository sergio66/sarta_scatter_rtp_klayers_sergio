! Microwave ocean surface emissivity model constants
! *COMMON/Constants/*
!      common/emismw/FreqGHz,Pangdeg,Pcc,Pc2,Ps2,emc
      real*4 emc(59)
      data emc/       
     & 17.535    ,    -0.61767 ,   0.008948   , 3.1842,
     & 0.019189  ,    -0.010873,   0.00025818 , 68.396,
     & -0.40643  ,    0.022832 ,   -0.00053061, 4.7629, 
     & 0.1541    ,    -0.033717,   0.00084428 , 78.287,
     & -0.0043463,    5.3125   ,   -0.011477  , 3.141592654,
     & -1.0      ,    1.95E-5  ,    2.55,
     & -6.37182  ,  0.0253918  ,3.57569e-05,
     & 9.42928   ,-0.0332839 ,-6.47724e-05,
     & -3.29282  , 0.00965450 , 2.81588e-05,
     & 0.252676  , 0.00343867 ,-1.56362e-05,
     & -0.000156669,  1.39485e-05 ,-4.07633e-08,
     & -0.141316 ,-0.00356556  ,1.42869e-05,
     & -2.40701  , -0.0563888  ,0.000325227,
     & 2.96005   , 0.0704675 ,-0.000426440,
     & -0.751252 ,  -0.0191934 , 0.000125937,
     & -0.288253 , -0.00102655 , 2.26701e-06,
     & -0.00119072, -2.63165e-05 , 1.14597e-07,
     & 0.406300  , 0.00200031 ,-7.81635e-06/

!Coefficients post October 29 change
!     & -6.37182  ,  0.0253918  ,3.57569e-05,
!     & 9.42928   ,-0.0332839 ,-6.47724e-05,
!     & -3.29282  , 0.00965450 , 2.81588e-05,
!     & 0.252676  , 0.00343867 ,-1.56362e-05,
!     & -0.000156669,  1.39485e-05 ,-4.07633e-08,
!     & -0.141316 ,-0.00356556  ,1.42869e-05,
!     & -2.40701  , -0.0563888  ,0.000325227,
!     & 2.96005   , 0.0704675 ,-0.000426440,
!     & -0.751252 ,  -0.0191934 , 0.000125937,
!     & -0.288253 , -0.00102655 , 2.26701e-06,
!     & -0.00119072, -2.63165e-05 , 1.14597e-07,
!     & 0.406300  , 0.00200031 ,-7.81635e-06/

!Coefficients post April 24 change
!     & -2.35188  ,    -0.132072,    0.000658883,
!     & 3.44125   ,  0.187539   , -0.000931400,
!     & -1.14758  , -0.0626640  ,  0.000308404,
!     & -0.232072 ,  0.0119520  , -4.01439e-05,
!     & 0.00549094,  1.51716e-06, -9.59848e-09,
!     & 0.238199  , -0.00881164 ,  2.71126e-05,
!     & 11.6962   , -0.346877   ,  0.00121357,
!     & -15.9742  ,  0.443442   , -0.00156718,
!     & 5.11736   , -0.127518   ,  0.000457518,
!     & -0.773612 ,  0.000775719,  1.80312e-06,
!     & 0.0103726 , -0.000133524,  3.65630e-07,
!     & 0.625144  ,  0.00625598 , -2.56224e-05/

! Coefficients pre-April 24 change
!                                                   -.271499E+01,
!     & -.285110E-01,   .178220E-03,   .394354E+01,  .433341E-01,
!     & -.262158E-03,  -.126083E+01,  -.160313E-01,   .910478E-04,
!     & -.609792E-01,   .368813E-02,  -.128188E-04,   .545718E-02,
!     & -.138377E-04,  -.192805E-09,   .980869E-01,  -.325372E-02,
!     & .103577E-04,   .445209E+01,  -.804326E-03,   .493464E-04,
!     & -.669273E+01,  -.798760E-02,  -.327148E-04,   .253287E+01,
!     & .700837E-02,  -.837905E-05,  -.519445E+00,  -.304075E-02,
!     & .140881E-04,   .690632E-02,  -.116684E-04,  -.123170E-07,
!     & .504731E+00,   .337427E-02,  -.160390E-04/
!
! Explanation:
! FreqGHz: Observation frequency in GHz
! Angdeg: local zenith angle
! Ci: cosine of local zenith angle
! CiCi: cosine squared of local zenith angle
! SiSi: sine squared of local zenith angle
! emc(38): Emissivity model data
! Permittivity model data (Lamkaouchi model)
!   [1-3]: Temperature polynomial coefficients for Tau1 - Lamkaouchi (1996)
!   [4-7]: Temperature polynomial coefficients for Tau2 - Lamkaouchi (1996)
!  [8-11]: Temperature polynomial coefficients for Del1 - Lamkaouchi (1996)
! [12-15]: Temperature polynomial coefficients for Del2 - Lamkaouchi (1996)
! [16-17]: Temperature polynomial coefficients for static permittivity - Lamkaouchi (1996)
! [18-19]: Temperature polynomial coefficients for infinite freq. permittivity - Lamkaouchi (1996)
! Pi is stored for good measure
!    [20]: Stored value of Pi  - temporary, use RTTOV pi when available.
! Large scale correction model version 1: This does *NOT* correct for
! hemispherical scattering and is *NO LONGER USED*
!    [21]: Angle coefficient for large scale correction - see English (1997)
!    [22]: Windspeed coefficient for large scale correction - see English (1997)
!    [23]: Constant for large scale correction - see English (1997)
!    [24]: Reference frequency for large scale correction - see English (1997)
!    [25]: Normalisation frequency for large scale correction - see English (1997)
!    [26]: Scaling factor for large scale correction - see English (1997)
! Bragg scattering correction coefficients
!    [27]: Scaling factor for small scale correction - see English (1997)
! Foam model coefficients for Monahan model
!    [28]: First coefficient in Monahan foam model (neutral stability)  - see English (1997)
!    [29]: Second coefficient in Monahan foam model (neutral stability) - see English (1997)
! Alternative permittivity model (Liebe)
!    [30]: a1 in Liebe's dielectric model - see Liebe (1989)
!    [31]: b1 in Liebe's dielectric model - see Liebe (1989)
!    [32]: c1 in Liebe's dielectric model - see Liebe (1989)
!    [33]: c2 in Liebe's dielectric model - see Liebe (1989)
!    [34]: d1 in Liebe's dielectric model - see Liebe (1989)
!    [35]: d2 in Liebe's dielectric model - see Liebe (1989)
!    [36]: d3 in Liebe's dielectric model - see Liebe (1989)
!    [37]: e1 in Liebe's dielectric model - see Liebe (1989)
!    [38]: e2 in Liebe's dielectric model - see Liebe (1989)
! Version 2 of large scale correction which *DOES»* take account of
! hemispherical scattering.
! 1.) Mixed polarisation mode (nominal V at nadir)
!    [39]: Term a00 in mixed pol of large scale correction model
!    [40]: Term a01 in mixed pol mode of large scale correction model
!    [41]: Term a02 in mixed pol mode of large scale correction model
!    [42]: Term a10 in mixed pol mode of large scale correction model
!    [43]: Term a11 in mixed pol mode of large scale correction model
!    [44]: Term a12 in mixed pol mode of large scale correction model
!    [45]: Term a20 in mixed pol mode of large scale correction model
!    [46]: Term a21 in mixed pol mode of large scale correction model
!    [47]: Term a22 in mixed pol mode of large scale correction model
!    [48]: Term a30 in mixed pol mode of large scale correction model
!    [49]: Term a31 in mixed pol mode of large scale correction model
!    [50]: Term a32 in mixed pol mode of large scale correction model
!    [51]: Term a40 in mixed pol mode of large scale correction model
!    [52]: Term a41 in mixed pol mode of large scale correction model
!    [53]: Term a42 in mixed pol mode of large scale correction model
! 2.) Vertical polarisation mode
!    [54]: Term a00 in vertical pol mode of large scale correction model
!    [55]: Term a01 in vertical pol mode of large scale correction model
!    [56]: Term a02 in vertical pol mode of large scale correction model
!    [57]: Term a10 in vertical pol mode of large scale correction model
!    [58]: Term a11 in vertical pol mode of large scale correction model
!    [59]: Term a12 in vertical pol mode of large scale correction model
!    [60]: Term a20 in vertical pol mode of large scale correction model
!    [61]: Term a21 in vertical pol mode of large scale correction model
!    [62]: Term a22 in vertical pol mode of large scale correction model
!    [63]: Term a30 in vertical pol mode of large scale correction model
!    [64]: Term a31 in vertical pol mode of large scale correction model
!    [65]: Term a32 in vertical pol mode of large scale correction model
!    [66]: Term a40 in vertical pol mode of large scale correction model
!    [67]: Term a41 in vertical pol mode of large scale correction model
!    [68]: Term a42 in vertical pol mode of large scale correction model 
! 3. ) Horizontal polarisation mode
!    [69]: Term a00 in horizontal pol mode of large scale correction model
!    [70]: Term a01 in horizontal pol mode of large scale correction model
!    [71]: Term a02 in horizontal pol mode of large scale correction model
!    [72]: Term a10 in horizontal pol mode of large scale correction model
!    [73]: Term a11 in horizontal pol mode of large scale correction model
!    [74]: Term a12 in horizontal pol mode of large scale correction model
!    [75]: Term a20 in horizontal pol mode of large scale correction model
!    [76]: Term a21 in horizontal pol mode of large scale correction model
!    [77]: Term a22 in horizontal pol mode of large scale correction model
!    [78]: Term a30 in horizontal pol mode of large scale correction model
!    [79]: Term a31 in horizontal pol mode of large scale correction model
!    [80]: Term a32 in horizontal pol mode of large scale correction model
!    [81]: Term a40 in horizontal pol mode of large scale correction model
!    [82]: Term a41 in horizontal pol mode of large scale correction model
!    [83]: Term a42 in horizontal pol mode of large scale correction model
!    [84]: Windspeed coefficient in mixed polarisation high U, theta correction
!    [85]: View angle coefficient in mixed polarisation high U, theta correction
!    [86]: Constant coefficient in mixed polarisation high U, theta correction
!    [87]: Windspeed coefficient in vertical polarisation high U, theta correction
!    [88]: View angle coefficient in vertical polarisation high U, theta correction
!    [89]: Constant coefficient in vertical polarisation high U, theta correction
!    [90]: Windspeed coefficient in horizontal polarisation high U, theta correction
!    [91]: View angle coefficient in horizontal polarisation high U, theta correction
!    [92]: Constant coefficient in horizontal polarisation high U, theta correction

