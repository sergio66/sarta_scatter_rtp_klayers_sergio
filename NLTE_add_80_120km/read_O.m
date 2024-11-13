addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/PLOTTER

afgl34 = quick_read_afgl(34,1);

atomicO = read_netcdf_lls('/home/sergio/KCARTA/NONLTE2/sergio/VT_48PROFILES_FOR_MANUEL/AtomicO_SABER_Mlynczak/atox_athy_night_YY2013_V1.02.nc');
atomicO.ppmv =  toppmv(atomicO.pressure*ones(1,length(atomicO.lat)),atomicO.ktemp,atomicO.qatox,16,12);

plot(nanmean(atomicO.alt,2),atomicO.pressure); set(gca,'ydir','reverse'); grid; xlabel('Hgt [km]'); ylabel('P [mb]')


daysINyear  = [0 31 28 31 30 31 30 31 31 30 31 30 31];
cdaysINyear = cumsum(daysINyear);
for ii = 1 : 12
  boo = find(atomicO.day > cdaysINyear(ii) & atomicO.day <= cdaysINyear(ii+1));
  atomicO.month(boo) = ii;
  figure(1); scatter_coast(atomicO.lon(boo),atomicO.lat(boo),10,atomicO.day(boo)); title(num2str(ii)); pause
end
figure(1); clf; simplemap(atomicO.lat,atomicO.lon,atomicO.month',5); colormap jet
  title('SABER months, 5 deg grid avg')

figure(3); clf; simplemap(atomicO.lat,atomicO.lon,atomicO.qatox(1,:)',5); colormap jet
figure(3); clf; loglog(nanmean(atomicO.qatox,2),atomicO.pressure); set(gca,'ydir','reverse')
  xlabel('Atomic O [VMR]'); ylabel('P [mb]'); grid

figure(3); clf; simplemap(atomicO.lat,atomicO.lon,atomicO.ppmv(1,:)',5); colormap jet
figure(3); clf; loglog(nanmean(atomicO.ppmv,2),atomicO.pressure,'b',afgl34.qstd,afgl34.pstd,'r'); set(gca,'ydir','reverse')
  xlabel('Atomic O [PPMV]'); ylabel('P [mb]'); hl = legend('SABER','AFGL34','location','best','fontsize',10); grid
figure(3); clf; loglog(nanmean(atomicO.ppmv,2),atomicO.pressure,'bx-',afgl34.qstd,afgl34.pstd,'r',...
                       nanmin(atomicO.ppmv,[],2),atomicO.pressure,'g',nanmax(atomicO.ppmv,[],2),atomicO.pressure,'k'); set(gca,'ydir','reverse')
  xlabel('Atomic O [PPMV]'); ylabel('P [mb]'); hl = legend('mean SABER','AFGL34','min SABR','max SABER','location','best','fontsize',10); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[meanO_X,pjunk] = add_afgl_g34(nanmean(atomicO.ppmv,2),atomicO.pressure,afgl34);
[minO_X,pjunk]  = add_afgl_g34(nanmin(atomicO.ppmv,[],2),atomicO.pressure,afgl34);
[maxO_X,pjunk]  = add_afgl_g34(nanmax(atomicO.ppmv,[],2),atomicO.pressure,afgl34);
figure(3); clf; 
  loglog(nanmean(atomicO.ppmv,2),atomicO.pressure,'bx-',afgl34.qstd,afgl34.pstd,'r',...
         nanmin(atomicO.ppmv,[],2),atomicO.pressure,'gx-',nanmax(atomicO.ppmv,[],2),atomicO.pressure,'kx-',...
         meanO_X,pjunk,'b',minO_X,pjunk,'g',maxO_X,pjunk,'k')
           set(gca,'ydir','reverse')
  xlabel('Atomic O [PPMV]'); ylabel('P [mb]'); hl = legend('mean SABER','AFGL34','min SABR','max SABER','location','best','fontsize',10); grid
xticks(10.^[-7:2:+6]); grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rlat = -90 : 10 : +90;
dayofmonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
cdayofmonth = cumsum(dayofmonth);
for ii = 1 : 12
  boo = find(atomicO.day > cdayofmonth(ii) & atomicO.day <= cdayofmonth(ii+1));
  atomicO.monthly(:,ii) = nanmean(atomicO.qatox(:,boo),2);
  for jj = 1 : length(rlat)-1
    boo = find(atomicO.day > cdayofmonth(ii) & atomicO.day <= cdayofmonth(ii+1) & ...
           atomicO.lat >= rlat(jj) & atomicO.lat < rlat(jj+1));
    for pp = 1 : 16
      atomicO.monthly_O(pp,ii,jj) = nanmean(atomicO.qatox(pp,boo));  %% 96 km
      atomicO.monthly_Z(pp,ii,jj) = nanmean(atomicO.alt(pp,boo));  %% 96 km
    end
  end
end

[mm,nn,oo] = size(atomicO.monthly_Z);
for ii = 1 : nn
  for jj = 1 : oo
    z = squeeze(atomicO.monthly_Z(:,ii,jj));
    q = squeeze(atomicO.monthly_O(:,ii,jj));
    boo = find(isfinite(z) & isfinite(q));
    if length(boo) > 2
      atomicO.q96(ii,jj)     = interp1(z(boo),q(boo),96); 
    end
    boo = find(z >= 95 & z <= 97);
    if length(boo) > 0
      atomicO.q96mean(ii,jj) = nanmean(q(boo));
    end
  end
end

figure(4); clf; pcolor(atomicO.monthly);     xlabel('Month'); ylabel('P [mb]'); colorbar; shading interp
figure(5); clf; pp = 3; pcolor(1:12,-85:10:+85,squeeze(atomicO.monthly_O(pp,:,:))'); xlabel('Month'); ylabel('Latitude'); colorbar; shading interp
figure(5); clf; pp = 3; pcolor(1:12,-85:10:+85,squeeze(atomicO.q96(:,:))'); xlabel('Month'); ylabel('Latitude'); colorbar; shading interp
figure(5); clf; pp = 3; pcolor(1:12,-85:10:+85,squeeze(atomicO.q96mean(:,:))'); xlabel('Month'); ylabel('Latitude'); colorbar; shading interp
  ylim([-1 +1]*40); caxis([0.009 0.026]); colormap jet



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
