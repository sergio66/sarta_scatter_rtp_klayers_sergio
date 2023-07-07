addpath /asl/matlib/h4tools
addpath /home/sergio/KCARTA/MATLAB

[hm1,~,pm1,~] = rtpread('../Data/mars_sergio_rtp.ip.rtp');
fid = fopen('/asl/ftp/pub/sergio/Mars/testMarsSim_ip.txt','w');
nlevs = pm1.nlevs;
fprintf(fid,'%% surf temp = %8.6f K  \n',pm1.stemp);
fprintf(fid,'%% surf pres = %8.6e mb \n',pm1.spres);
fprintf(fid,'%% \n');
fprintf(fid,'%% ppmv : N2 = 0.0259, O2 = 0.00160, CO = 0.00060, NO = 0.00010 \n');
fprintf(fid,'%% \n');
fprintf(fid,'%%  plevs(z)       T(z)       WV(z)      CO2(z)    \n');
fprintf(fid,'%%    mb            K         ppm        ppm       \n');
fprintf(fid,'%% ----------------------------------------------- \n');
data = [pm1.plevs(1:nlevs)  pm1.ptemp(1:nlevs)  pm1.gas_1(1:nlevs)  pm1.gas_2(1:nlevs)];
  whos data
fprintf(fid,'%8.6e %8.6f %8.6e %8.6e \n',data');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
units of .op.rtp are molecules/cm2

Example : for lowest layer in tropics, CO2 = 400 ppm  L = 250 m
pV = n R T ==? q = (n/V)L = (p/RT)L = 1013 mb * (100 mb--> N/m2) / 8.31 / 300 K * 250 m * 400 * 1e-6 (ppm --> ratio) = 101300/8.31/300*250*400*1e-6 = 4.06 mol/m2
 4.06*6.023e23*1e-4 (convert to molecules/m2, and then to molecules/cm2) = 2.477e20 molecules/cm2

format short e; [pRTP.ptemp(lowest) pRTP.gas_2(lowest)]
  299.01  2.3329e+20

so pav (mb), L (m),  MR (ppmv), T (K)  -->> (pav * L /T) * ppm * 100/8.31*1e-6*6.023e23*1e-4   molecules/cm2 == (pav * L/T) * ppm * 7.2479e+14
%}

%{
lowest layer for what I gave Giuliano   ppm = MR * 1e6
   pav          L          T            WV         CO2         N2            O2           CO           NO
   1.289381e+01 315.218750 248.789474 2.484778e+12 1.124066e+16 3.064560e+14 1.893164e+13 7.099367e+12 1.183228e+12
                                                          ^^^^^
(12 * 315/249) * (0.95 * 1e6) * 7.2479e+14 = 1.0453e+22      |
   OOOPOPS 1e6 too small oooeeerrrr   -----------------------|
%}

if ~exist('pm')
  [hm,~,pm,~] = rtpread('../Data/mars_sergio_rtp.op.rtp');
end
fid = fopen('/asl/ftp/pub/sergio/Mars/testMarsSim_op.txt','w');
nlevs = pm.nlevs;
nlays = pm.nlevs-1;
pN = pm.plevs(1:100)-pm.plevs(2:101);
pD = log(pm.plevs(1:100)./pm.plevs(2:101));
pm.plays = pN./pD;
fprintf(fid,'%% surf temp = %8.6f K  \n',pm.stemp);
fprintf(fid,'%% surf pres = %8.6e mb \n',pm.spres);
fprintf(fid,'%% \n');

fprintf(fid,'%%   plays       palts(z)      T(z)         WV(z)       CO2(z)        N2(z)        O2(z)       CO(z)        NO(z)     \n');
fprintf(fid,'%%    mb           m             K       <--------------------- molecules/cm2 --------------------------------------> \n');
fprintf(fid,'%% -----------------------------------------------------------------------------------------------------------------  \n');
data = [pm.plays(1:nlays) pm.palts(1:nlays) pm.ptemp(1:nlays)  pm.gas_1(1:nlays)  pm.gas_2(1:nlays) pm.gas_3(1:nlays) pm.gas_4(1:nlays) pm.gas_5(1:nlays) pm.gas_6(1:nlays)];
  whos data
fprintf(fid,'%8.6e %8.6f %8.6f %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e \n',data');
fclose(fid);

if ~exist('dEarth')
  [dEarth,w] = readkcstd('earth_kcarta_122.dat');
end
if ~exist('dMars')
  [dMars,w] = readkcstd('mars_kcarta.dat');
end
fid = fopen('/asl/ftp/pub/sergio/Mars/testMarsSim_rad_0p0025.txt','w');
data = [w; dMars']; whos data
fprintf(fid,'%15.6f %16.6e \n',data);
fclose(fid);

if ~exist('qc')
  addpath /home/sergio/MATLABCODE
  [fc,qc] = quickconvolve(w,[dEarth dMars],0.5,0.5);
end
fid = fopen('/asl/ftp/pub/sergio/Mars/testMarsSim_rad_0p5.txt','w');
data = [fc; qc(:,2)']; whos data
fprintf(fid,'%15.6f %16.6e \n',data);
fclose(fid);

figure(1); plot(fc,rad2bt(fc,qc)); hl = legend('Earth Tropical','Mars Std','location','best','fontsize',10);
  xlabel('Wavenumber cm-1'); ylabel('BT(K)');

junk1 = load('/asl/ftp/pub/sergio/Mars/testMarsSim_rad_0p0025.txt'); figure(3); plot(junk1(:,1),junk1(:,2))
junk2 = load('/asl/ftp/pub/sergio/Mars/testMarsSim_rad_0p5.txt'); figure(3); plot(junk1(:,1),junk1(:,2),'b',junk2(:,1),junk2(:,2),'r')
