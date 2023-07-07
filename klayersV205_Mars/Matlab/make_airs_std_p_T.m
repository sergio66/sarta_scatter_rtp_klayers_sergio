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

%% original in May 21, 2021
T = T + 273;
n = p ./ (0.1921 * T);   %% density using eqn of state

%% fixed?? June 7, 2021
%% p V = n R T ==> n/V = p/RT  moles/m3 = p*1000*6.023e23/RT molecules/m3         p in kPa so x1000 to put in Pa
n = p*1000*6.023e23 ./ (8.31 *T);

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

iWrite = +1;
if iWrite > 0
  disp('<<< writing to Grid/cbplev_mars0.f ... after this you need to edit it and put the header info from Grid/cbplev_mars_header.f at the INSERT HERE >>>')
  disp('<<< writing to Grid/cbplev_mars0.f ... after this you need to edit it and put the header info from Grid/cbplev_mars_header.f at the INSERT HERE >>>')
  disp('<<< writing to Grid/cbplev_mars0.f ... after this you need to edit it and put the header info from Grid/cbplev_mars_header.f at the INSERT HERE >>>')
  disp(' ')
  disp('so eg combine Grid/cbplev_mars_header.f with Grid/cbplev_marsX.f')
  disp(' ')
  fid = fopen('../Grid/cbplev_marsX.f','w');
  fprintf(fid,'      DATA (PLEV(I), I = 1,101) \n');
  fprintf(fid,'     $ / \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'     $ %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',P(ind));
  end
  ind = 101;
  fprintf(fid,'     $ %11.5f \n',P(ind));
  fprintf(fid,'     $ / \n');
%  fprintf(fid,'      END \n');  
  fclose(fid);
end

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

iWrite = +1;
if iWrite > 0
  disp('<<< writing to Data/glmars.dat >>>')
  disp(' ')
  fid = fopen('../Data/glmars.dat','w');
  fprintf(fid,'! No. of levels in these profiles: \n')
  fprintf(fid,'101 \n');
  fprintf(fid,'!\n');
  fprintf(fid,'MODEL 1 : SERGIO \n');
  fprintf(fid,'!\n')
  fprintf(fid,'! Latitude (deg) \n');
  fprintf(fid,'0.0 \n');
  fprintf(fid,'!\n');

  fprintf(fid,'! Altitude (km) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',HP(ind)/1000);
  end
  ind = 101;
  fprintf(fid,'  %11.5f \n',HP(ind)/1000);
  fprintf(fid,'!\n');

  fprintf(fid,'! Pressure (mb) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',P(ind));
  end
  ind = 101;
  fprintf(fid,'  %11.5f \n',P(ind));
  fprintf(fid,'!\n');

  fprintf(fid,'! Temperature (K) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',TP(ind));
  end
  ind = 101;
  fprintf(fid,'  %11.5f \n',TP(ind));
  fprintf(fid,'!\n');

  fprintf(fid,'! Density (cm-3) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5e, %11.5e, %11.5e, %11.5e, %11.5e,  \n',NP(ind));
  end
  ind = 101;
  fprintf(fid,'  %11.5e \n',NP(ind));
  fprintf(fid,'!\n');

  gas = ones(size(TP)) * 210;
  fprintf(fid,'1 : WV (ppmv) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',gas(ind));
  end
  ind = 101;
  fprintf(fid,'  %11.5f \n',gas(ind));
  fprintf(fid,'!\n');

  gas = ones(size(TP)) * 95/100 * 1e6;
  fprintf(fid,'2 : CO2 (ppmv) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',gas(ind));
  end
  ind = 101;
  fprintf(fid,'  %11.5f \n',gas(ind));
  fprintf(fid,'!\n');

  gas = ones(size(TP)) * 2.59/100 * 1e6;
  fprintf(fid,'3 : N2 (ppmv) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',gas(ind));
  end
  ind = 101;
  fprintf(fid,'  %11.5f \n',gas(ind));
  fprintf(fid,'!\n');

  gas = ones(size(TP)) * 0.16/100 * 1e6;
  fprintf(fid,'4 : O2 (ppmv) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',gas(ind));
  end
  ind = 101;
  fprintf(fid,'  %11.5f \n',gas(ind));
  fprintf(fid,'!\n');

  gas = ones(size(TP)) * 0.06/100 * 1e6;
  fprintf(fid,'5 : CO (ppmv) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',gas(ind));
  end
  ind = 101;
  fprintf(fid,'  %11.5f \n',gas(ind));
  fprintf(fid,'!\n');

  gas = ones(size(TP)) * 100;
  fprintf(fid,'6 : NO (ppmv) \n');
  for ii = 1 : 20
    ind = (1:5) + (ii-1)*5;
    fprintf(fid,'  %11.5f, %11.5f, %11.5f, %11.5f, %11.5f,  \n',gas(ind));
  end
  ind = 101;
  fprintf(fid,'  %11.5f \n',gas(ind));
  fprintf(fid,'!\n');

  fprintf(fid,'MODEND\n');
  fprintf(fid,'!\n');
  fprintf(fid,'! CONSTITUENT PROFILES FOR THE MINOR ABSORBING ATMOSPHERIC GASES\n');
  fprintf(fid,'!\n');
  fprintf(fid,'MINGAS\n');
  fprintf(fid,'!\n');
  fprintf(fid,'DATAEND\n');

  fclose(fid);
end
