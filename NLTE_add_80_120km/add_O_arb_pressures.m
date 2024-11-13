function [saber_O_pall,saber_O,saber_p] = add_O_arb_pressures(hall,pall)

addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD

%% this takes 
%%   hall,pall = rtp levels tructures
%% and outputs
%%   saber_O = saber atomic oxygen profile      at AFGL 50 levels from glatm_16Aug2010.dat
%%   saber_p = the 50 pressure levels           at AFGL 50 levels from glatm_16Aug2010.dat
%%   saber_O_pall = saber atomic oxygen profile at pall pressure levels

afgl34 = quick_read_afgl(34,1);

[yy,mm,dd,hh] = tai2utcSergio(pall.rtime);
yyx = yy + (mm-1)/12 + (dd-1)/30/12;

atomicO = read_netcdf_lls('/home/sergio/KCARTA/NONLTE2/sergio/VT_48PROFILES_FOR_MANUEL/AtomicO_SABER_Mlynczak/atox_athy_night_YY2014_V1.02.nc');
atomicO.ppmv =  toppmv(atomicO.pressure*ones(1,length(atomicO.lat)),atomicO.ktemp,atomicO.qatox,16,12);
atomicO.lon = wrapTo180(atomicO.lon);
daysINyear  = [0 31 28 31 30 31 30 31 31 30 31 30 31];
cdaysINyear = cumsum(daysINyear);
for ii = 1 : 12
  boo = find(atomicO.day > cdaysINyear(ii) & atomicO.day <= cdaysINyear(ii+1));
  atomicO.month(boo) = ii;
  %figure(1); scatter_coast(atomicO.lon(boo),atomicO.lat(boo),10,atomicO.day(boo)); title(num2str(ii)); pause
end

%% match pall(rlat/rlon/month) to atomicO(lat/lon/month)

for ii = 1 : 12
  fprintf(1,'month = %2i \n',ii)
  booO = find(atomicO.month == ii);
  booP = find(mm == ii);
  [r,closestind_2to1,closestind_1to2,ind1_maxdist,ind2_maxdist] = haversine3_matrix(pall.rlat(booP),pall.rlon(booP),atomicO.lat(booO),atomicO.lon(booO),50000000);
  figure(1); scatter_coast(pall.rlon(booP),pall.rlat(booP),10,ones(size(booP))); title('ECMWF SAF'); ylim([-90 +90])
  figure(2); scatter_coast(atomicO.lon(booO(closestind_2to1)),atomicO.lat(booO(closestind_2to1)),10,ones(size(booP))); title('SABER'); ylim([-90 +90])
  figure(3); plot(pall.rlon(booP),pall.rlat(booP),'b.',atomicO.lon(booO(closestind_2to1)),atomicO.lat(booO(closestind_2to1)),'r.')
  for jj = 1 : length(booP)
    [saber_O(:,booP(jj)),saber_p] = add_afgl_g34(atomicO.ppmv(:,booO(closestind_2to1(jj))),atomicO.pressure,afgl34);
  end
  pause(0.1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : length(pall.stemp)
  pjunk = pall.plevs(:,ii);
  saber_O_pall(:,ii) = interp1(log(saber_p),saber_O(:,ii),log(pjunk));
end
