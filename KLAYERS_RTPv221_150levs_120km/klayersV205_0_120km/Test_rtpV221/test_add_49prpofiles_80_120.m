addpath /home/sergio/KCARTA/NONLTE2/sergio/VT_48PROFILES_FOR_MANUEL/AtomicO_SABER_Mlynczak/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/EOF_PCA

[h0,ha,p0,pa] = rtpread('/asl/packages/klayersV205/Data/adafgl_16Aug2010_ip.rtp');
[h,ha,p,pa] = rtpread('/asl/s1/sergio/RTP_pin_feb2002/pin_feb2002_sea_airsnadir_ip.so2.latlon.const_emiss.rtp');
for ii = 1 : 49
  fprintf(1,'ii p.nlevs ptop = %3i %3i %12.8e \n',ii,p.nlevs(ii),p.plevs(p.nlevs(ii),ii))
end

[h,p] = subset_rtp(h,p,[],[],6:48);
for ii = 1 : 43
  fprintf(1,'ii p.nlevs ptop = %3i %3i %12.8e \n',ii,p.nlevs(ii),p.plevs(p.nlevs(ii),ii))
end

%% [h,p] = subset_rtp(h,p,[],[],[1 2]);

p.gas_1 = toppmv(p.plevs,p.ptemp,p.gas_1,18,21);
p.gas_2 = toppmv(p.plevs,p.ptemp,p.gas_2,44,21);
p.gas_3 = toppmv(p.plevs,p.ptemp,p.gas_3,48,21);
p.gas_5 = toppmv(p.plevs,p.ptemp,p.gas_5,28,21);
p.gas_6 = toppmv(p.plevs,p.ptemp,p.gas_6,16,21);
p.gas_9 = toppmv(p.plevs,p.ptemp,p.gas_9,66,21);
h.gunit = 10 * ones(size(h.gunit));

p.rtime = utc2taiSergio(2002,01,01,12) * ones(size(p.stemp));
more_p = linspace(log10(2.5e-5),log10(1e-3),10); more_p = exp10(fliplr(more_p)); more_p = more_p';
[profnew,proftop] = add_othergases_arb_pressures(h,p,more_p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a fake dataset with 4 spatial points. Note that we are arranging
% our array with time as the first dimension, increasing with row
% number, so each column is a time series. We could use any number of
% spatial locations (number of columns), and they don't have to be
% arranged any particular way. They could be all points in a rectangular
% lon-lat grid, for example, strung out out into a single list of points
% by flattening the 2-D array. To begin with, it is easiest to see how
% the calculation works by keeping the number of spatial points very
% small, hence our choice of 4 points; but in general there can even be
% more spatial points than temporal samples, and in fact there is
% nothing special about "space" and "time", or the choice of which one
% is the first dimension and which is the second.

% plev series x spatial locations

wah = detrend(proftop.ptemp);
figure(7); 
  subplot(121); semilogy(proftop.ptemp,proftop.plevs); title('proftop(p)'); 
    set(gca,'ydir','reverse')
  subplot(122); semilogy(wah,proftop.plevs); title('detrend(p)')
    set(gca,'ydir','reverse')

[L, EOFs, EC, error] = driver_eof1(proftop.ptemp);
[L, EOFs, EC, error] = driver_eof1(detrend(proftop.ptemp));
figure(8); semilogy(L,'o-'); title('Sergio 49 regr eigenvalues')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9)

ch = load('saf_sub300_reg49_plv147_400ppm_v3');
newT = ch.pcat.ptemp(:,7:349);
newP = ch.pcat.plevs(:,7:349);

%% this shows above 100 mb, the plevs are about the same
semilogx(nanstd(newP,[],2),nanmean(newP,2)); set(gca,'ydir','reverse')
plot(nanstd(newP,[],2),nanmean(newP,2)); set(gca,'ydir','reverse'); 

ua = find(nanmean(newP,2) < 50);
ua = find(nanmean(newP,2) < 1);
ua = find(nanmean(newP,2) < 0.01);  %% see psaf below, min(psaf.plevs) == 0.01)
loglog(nanstd(newP(ua,:),[],2),nanmean(newP(ua,:),2)); set(gca,'ydir','reverse')
  xlabel('std dev P');  ylabel('P [mb]')

wah = detrend(newT(ua,:));
figure(9); 
  subplot(121); semilogy(newT(ua,:),newP(ua,:)); title('proftop(p)'); 
    set(gca,'ydir','reverse')
  subplot(122); semilogy(wah,newP(ua,:)); title('detrend(p)')
    set(gca,'ydir','reverse')

pah = newT(ua,:);
[L, EOFs, EC, error] = driver_eof1(pah);
[L, EOFs, EC, error] = driver_eof1(detrend(pah));
figure(10); semilogy(L,'o-'); title('Chris 300 eigenvalues')
figure(10); semilogy(L/max(L),'o-'); title('Chris 300 eigenvalues')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hsaf,~,psaf,~] = rtpread('/home/sergio/MATLABCODE_Git/REGR_PROFILES_SARTA/ECMWF_SAF_137Profiles/save_SAF_25000_profiles_25-May-2016_xmb_400ppmv.ip.rtp');
plot(psaf.plevs(1,:))
[min(psaf.plevs(1,:)) max(psaf.plevs(1,:)) mean(psaf.plevs(1,:)) std(psaf.plevs(1,:))]

wv = load('/home/sergio/MATLABCODE_Git/REGR_PROFILES_SARTA/ECMWF_SAF_137Profiles/ECMWF_RTTOV_91_25000_and_704_Profiles/Vers137_2024/rtp_q.mat');
