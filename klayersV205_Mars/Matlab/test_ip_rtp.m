addpath /asl/matlib/h4tools

%% see https://www.grc.nasa.gov/www/k-12/airplane/atmosmrm.html

h = 0:250:120000;
p = zeros(size(h));
T = zeros(size(h));

%% T in C, h in meters, p in kPa
oo = find(h <= 7000);
 T(oo) = -31 - 0.000998 * h (oo);       %% -dT/dh = 0.998 K/km
 p(oo) = .699 * exp(-0.00009 * h(oo));

oo = find(h > 7000);
  T(oo) = -23.4 - 0.00222 * h(oo) + 1e-8*h(oo).^2;      %% -dT/dh = 2.2 K/km, I did the quadratic modification 
  p(oo) = .699 * exp(-0.00009 * h(oo));

T = T + 273;
n = p ./ (0.1921 * T);   %% density using eqn of state

p = p*1000/100;  %% KPa to Pa = N/m2 to mb
T = T;           %% in K
n = n/1e6;       %% per cm3

figure(1); plot(T,h/1000); ylabel('h(km)'); xlabel('T(K)');
figure(2); plot(p,h/1000); ylabel('h(km)'); xlabel('p(mb)');
figure(3); plot(n,h/1000); ylabel('h(km)'); xlabel('n(molecules/cm3)');

%% An Overview of the AIRS Radiative Transfer Model, 
%% IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. 41, NO. 2, FEBRUARY 2003

%% All component gas transmittances, with the exception of
%% water vapor for 600 channels, are parameterized on this grid
%% which spans the range 1100–0.005 hPa. The 101 atmospheric
%% pressure levels which divide the atmosphere into 100 layers are
%% defined as
%% (19) Plev(i) = (Ai^2 + Bi + C).^7/2
%% where is level number, and A,B,C are constants. By
%% fixing P(1) = 1100, p(38) = 300, p(101) = 0.005, , and
%% hPa, we can then solve for these three constants. This relation gives us smoothly varying layers and is fine enough to
%% not limit the accuracy of the radiative transfer equation

%{
for Mars we want P(1) = 13 mb; P(40) = 0.3; P(101) = 0.001 mb
A + B + C = 13^2/7
1600A + 40B + C = 0.3^2/7
10201A + 101B + C = 0.01^2/7

(1       1 1)(A)   (13   )
(1600   40 1)(B) = (0.3  ) ^ (2/7)
(10201 101 1)(C)   (0.001)
%}
hp = 38; matr = [1 1 1; hp*hp hp 1; 10201 101 1]; yvec = [1100 300 0.005].^(2/7);  ABC = matr\yvec'  %% for AIRS
figure(3); A = ABC(1); B = ABC(2); C = ABC(3); ii = 1:101; P = (A*ii.^2 + B*ii + C).^(7/2); semilogy(ii,P,'o-')

hp = 80; matr = [1 1 1; hp*hp hp 1; 10201 101 1]; yvec = [13 0.3 0.001].^(2/7);  ABC = matr\yvec'     %% for Mars
figure(4); A = ABC(1); B = ABC(2); C = ABC(3); ii = 1:101; P = (A*ii.^2 + B*ii + C).^(7/2); semilogy(ii,P,'o-')

TP = interp1(log(p),T,log(P),[],'extrap'); figure(3); semilogy(T,p,TP,P); set(gca,'ydir','reverse');
HP = interp1(log(p),h,log(P),[],'extrap'); figure(4); semilogy(h,p,HP,P); set(gca,'ydir','reverse');
NP = interp1(log(p),n,log(P),[],'extrap'); figure(2); semilogy(n,p,NP,P); set(gca,'ydir','reverse');

%% https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
%% Surface density: ~0.020 kg/m3
%% Atmospheric composition (by volume) : PPMV = 10^6 VMR
%%     Major      : Carbon Dioxide (CO2) - 95.1% ; Nitrogen (N2) - 2.59%
%%                  Argon (Ar) - 1.94%; Oxygen (O2) - 0.16%; Carbon Monoxide (CO) - 0.06% 
%%     Minor (ppm): Water (H2O) - 210; Nitrogen Oxide (NO) - 100; Neon (Ne) - 2.5;
%%                  Hydrogen-Deuterium-Oxygen (HDO) - 0.85; Krypton (Kr) - 0.3; 
%% 		 			    Xenon (Xe) - 0.08
%% 
%% when we write ../Data/glmars.dat we see surface denity 3.31970e+17 molecules/cm3
%%  so surface density = (3.32e17*1e6) * (44/1000) / 6.023e23   where I convert to molecules/m3 and CO2 molar mass = 44g/mol
%%  0.0243 kg/m3 YAY
%%
%% /asl/rta/kcarta_sergio/KCDATA/RefProf_Mars/refgas2 says lowemost layer has 1.8660463e-05 kmol/cm2
%% so lowest layer = 300 m, 0.95 VMR ==> q = (3.32e17*1e6) * 0.75 * 300 molecules/m2 = (3.32e17*1e6) * 0.75 * 300 / 1e4 molecules/cm2 
%%                                         = (3.32e17*1e6) * 0.75 * 300 / 1e4 /6e26 kmol/cm2 = 1.24e-5  YAY

hx.ptype = 0;
hx.pfields = 1;
hx.pmin = 0;
hx.pmax = 10;
hx.ngas = 2;
hx.glist = [1 2]';
hx.gunit = [10 10]';
hx.vcmin = 605;
hx.vcmax = 2805;

px.rlat = [0   -20]; px.plat = px.rlat;
px.rlon = [-90 +90]; px.plon = px.rlon;
px.nlevs = [50 50];
px.plevs = zeros(px.nlevs(1),2);
px.ptemp = zeros(px.nlevs(1),2);
px.gas_1 = zeros(px.nlevs(1),2);
px.gas_2 = zeros(px.nlevs(1),2);
px.satzen = [10 20];
px.solzen = [10 120];
px.spres = [8 12];
px.salti = [1000 100];

px.plevs(:,1) = P(1:2:100);
px.plevs(:,2) = P(2:2:101);

  sinsin = 1 * cos((1:50)'/50*2*pi);
px.ptemp(:,1) = TP(1:2:100) + sinsin';
  sinsin = 2 * sin((1:50)'/25*2*pi);
px.ptemp(:,2) = TP(1:2:100) - sinsin';
px.stemp(1) = max(px.ptemp(:,1)) + 1;
px.stemp(2) = max(px.ptemp(:,2)) - 1;

% RECALL PPMV = 10^6 VMR
% May 26, 2021. I messsed up MR and ppm!!! these are MR  
% px.gas_1(:,1) = 100*1e-6; px.gas_1(:,2) = 50*1e-6;
% px.gas_2(:,1) = 0.95;     px.gas_2(:,2) = 0.96;

% June 7, 2021
px.gas_1(:,1) = 200;          px.gas_1(:,2) = 50;
px.gas_2(:,1) = 0.97*1e6;     px.gas_2(:,2) = 0.96*1e6;

px.rlat  = single(px.rlat);
px.plat  = single(px.plat);
px.rlon  = single(px.rlon);
px.plon  = single(px.plon);
px.gas_1 = single(px.gas_1);
px.gas_2 = single(px.gas_2);
px.ptemp = single(px.ptemp);
px.plevs = single(px.plevs);
px.stemp = single(px.stemp);
px.spres = single(px.spres);
px.satzen  = single(px.satzen);
px.solzen  = single(px.solzen);

rtpwrite('test_ip.rtp',hx,[],px,[]);

semilogy(px.ptemp,px.plevs); set(gca,'ydir','reverse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

klayerser = ['!../BinV201/klayers_mars_wetwater fin=test_ip.rtp fout=test_op.rtp'];
eval(klayerser)

[hx2,~,px2,~] = rtpread('test_op.rtp');
pN = px2.plevs(1:100,:)-px2.plevs(2:101,:);
pD = log(px2.plevs(1:100,:)./px2.plevs(2:101,:));
px2.plays = pN./pD;
px2.plays(101,:) = 0;

figure(1)
semilogy(px.ptemp(:,1),px.plevs(:,1),'co-',px2.ptemp(1:px2.nlevs(1),1),px2.plays(1:px2.nlevs(1),1),'b',px.ptemp(:,2),px.plevs(:,2),'mo-',px2.ptemp(1:px2.nlevs(2),2),px2.plays(1:px2.nlevs(2),2),'r')
  set(gca,'ydir','reverse');

figure(2)
loglog(px2.gas_1(1:px2.nlevs(1),1),px2.plays(1:px2.nlevs(1),1),px2.gas_1(1:px2.nlevs(2),2),px2.plays(1:px2.nlevs(2),2)); set(gca,'ydir','reverse');

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
RH = layeramt2RH(hx2,px2);
figure(3); semilogy(RH,px2.plays(1:100,:))

[ppmv1] = layers2ppmv(hx2,px2,1:2,1); figure(4); semilogy(ppmv1,px2.plays(1:96,:))
loglog(px.gas_1(:,1),px.plevs(:,1),'co-',ppmv1(:,1),px2.plays(1:96,1),'b',px.gas_1(:,2),px.plevs(:,2),'mo-',ppmv1(:,2),px2.plays(1:96,2),'r')

[ppmv2] = layers2ppmv(hx2,px2,1:2,2); figure(4); semilogy(ppmv2,px2.plays(1:96,:))
loglog(px.gas_2(:,1),px.plevs(:,1),'co-',ppmv2(:,1),px2.plays(1:96,1),'b',px.gas_2(:,2),px.plevs(:,2),'mo-',ppmv2(:,2),px2.plays(1:96,2),'r')

semilogx(P,HP/1000,'co-',px2.plays(1:px2.nlevs(1)-1,1),px2.palts(1:px2.nlevs(1)-1,1)/1000,'b.-',...
         P,HP/1000,'mo-',px2.plays(1:px2.nlevs(1)-1,1),px2.palts(1:px2.nlevs(1)-1,1)/1000,'r'); xlabel('p(mb)'); ylabel('H(km)')
