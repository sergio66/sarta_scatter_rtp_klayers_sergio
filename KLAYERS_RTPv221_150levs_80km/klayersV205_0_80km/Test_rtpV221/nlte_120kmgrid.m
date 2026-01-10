if ~exist('pcat')
  test_300_49
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first find z(p) ie the heights at the plevs in Chris original mat file
pcatx.palts = zeros(size(pcatx.plevs));
for ii = 1 : iNumProf
  klayers_nlevs = p120.nlevs(ii);
  klayers_plevs = p120.plevs(1:klayers_nlevs,ii);
  klayers_palts = p120.palts(1:klayers_nlevs,ii);
  max_klayers_alts(ii) = max(klayers_palts);
  min_klayers_plevs(ii) = min(klayers_plevs);

  nlevs = pcatx.nlevs(ii);
  plevs = pcatx.plevs(1:nlevs,ii);
  palts = interp1(log(klayers_plevs),klayers_palts,log(plevs),[],'extrap');
  pcatx.palts(1:nlevs,ii) = palts;
end
figure(1); clf; semilogy(max_klayers_alts/1000,min_klayers_plevs)
  xlabel('Max klayers0\_120 Z [km]'); ylabel('Min klayers0\_120 p [mb]')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now interpolate everything to a uniform H grid

Hmax = floor(nanmin(p120.palts(1,:))/1000);
fprintf(1,'Hmax = %8.2f km \n',Hmax)

Hmax = floor(nanmin(pcatx.palts(1,:))/1000);
fprintf(1,'Hmax = %8.2f km \n',Hmax)

zgrid = 0 : Hmax;
zgrid = zgrid * 1000;

hnlte = hcatx;
pnlte = pcatx;
pnlte = rmfield(pnlte,'plevs');
pnlte = rmfield(pnlte,'ptemp');
pnlte = rmfield(pnlte,'palts');
pnlte = rmfield(pnlte,'gas_1');
pnlte = rmfield(pnlte,'gas_2');
pnlte = rmfield(pnlte,'gas_3');
pnlte = rmfield(pnlte,'gas_5');
pnlte = rmfield(pnlte,'gas_6');
pnlte = rmfield(pnlte,'gas_9');
pnlte = rmfield(pnlte,'gas_1_press');
pnlte = rmfield(pnlte,'altitude');

for ii = 1 : iNumProf
  nlevs = pcatx.nlevs(ii);
  palts = pcatx.palts(1:nlevs,ii);
  plevs = pcatx.plevs(1:nlevs,ii);
  ptemp = pcatx.ptemp(1:nlevs,ii);
  gas_1 = pcatx.gas_1(1:nlevs,ii);
  gas_2 = pcatx.gas_2(1:nlevs,ii);
  gas_3 = pcatx.gas_3(1:nlevs,ii);
  gas_5 = pcatx.gas_5(1:nlevs,ii);
  gas_6 = pcatx.gas_6(1:nlevs,ii);
  gas_9 = pcatx.gas_9(1:nlevs,ii);
  
  boo = find(isfinite(palts) & isfinite(plevs) & gas_1 >= 0 & gas_2 >= 0 & gas_3 >= 0 & gas_6 >= 0);
  palts = palts(boo);
  plevs = plevs(boo);
  ptemp = ptemp(boo);
  gas_1 = gas_1(boo);
  gas_2 = gas_2(boo);
  gas_3 = gas_3(boo);
  gas_5 = gas_5(boo);
  gas_6 = gas_6(boo);
  gas_9 = gas_9(boo);

  plevsX = interp1(palts,log(plevs),zgrid,[],'extrap'); plevsX = exp(plevsX);
  ptempX = interp1(palts,ptemp,zgrid,[],'extrap'); ptempX = ptempX;
  gas_1X = interp1(palts,log(gas_1),zgrid,[],'extrap'); gas_1X = exp(gas_1X);
  gas_2X = interp1(palts,log(gas_2),zgrid,[],'extrap'); gas_2X = exp(gas_2X);
  gas_3X = interp1(palts,log(gas_3),zgrid,[],'extrap'); gas_3X = exp(gas_3X);
  gas_5X = interp1(palts,log(gas_5),zgrid,[],'extrap'); gas_5X = exp(gas_5X);
  gas_6X = interp1(palts,log(gas_6),zgrid,[],'extrap'); gas_6X = exp(gas_6X);
  gas_9X = interp1(palts,log(gas_9),zgrid,[],'extrap'); gas_9X = exp(gas_9X);
  
  pnlte.plevs(:,ii) = plevsX;
  pnlte.ptemp(:,ii) = ptempX;
bad = find(ptempX < 125);
if length(bad) > 0
  plot(ptemp,palts,zgrid,ptempX)
  plot(ptemp,palts/1000,'b.',ptempX,zgrid/1000,'ro'); legend('Chris','Sergio');
  error('bad T')
end

  pnlte.gas_1(:,ii) = max(0,gas_1X);
  pnlte.gas_2(:,ii) = max(0,gas_2X);
  pnlte.gas_3(:,ii) = max(0,gas_3X);
  pnlte.gas_5(:,ii) = max(0,gas_5X);
  pnlte.gas_6(:,ii) = max(0,gas_6X);
  pnlte.gas_9(:,ii) = max(0,gas_9X);
  pnlte.palts(:,ii) = zgrid;
  pnlteX.nlevs(ii) = length(zgrid);
end

figure(1); clf; pcolor(1:iNumProf,pcatx.plevs,pcatx.ptemp); colormap jet; shading interp; colorbar; title('Torig'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); cx = caxis;
figure(2); clf; pcolor(1:iNumProf,pnlte.plevs,pnlte.ptemp); colormap jet; shading interp; colorbar; title('Tnew'); set(gca,'ydir','reverse'); set(gca,'yscale','log');  caxis(cx);

figure(3); clf; pcolor(1:iNumProf,pcatx.plevs,(pcatx.gas_1)); colormap jet; shading interp; colorbar; title('WVorig'); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); cx = caxis;
figure(4); clf; pcolor(1:iNumProf,pnlte.plevs,(pnlte.gas_1)); colormap jet; shading interp; colorbar; title('WVnew'); set(gca,'ydir','reverse'); set(gca,'yscale','linear');  caxis(cx);

raLat = -90:10:+90;
for ii = 1 : length(raLat)-1
  boo = find(pcat.rlat >= raLat(ii) & pcat.rlat < raLat(ii+1));

  numfound(ii) = length(boo);
  meanT(:,ii)  = nanmean(p120.ptemp(:,boo),2); 
  stdT(:,ii)   = nanstd(p120.ptemp(:,boo),[],2); 
  meanWV(:,ii) = nanmean(p120.gas_1(:,boo),2); 
  stdWV(:,ii)  = nanstd(p120.gas_1(:,boo),[],2); 

  daT  = p120.ptemp(:,boo);
  daWV = p120.gas_1(:,boo);
  for jj = 1 : length(boo)
    nn = pcat.nlevs(boo(jj)) - 1;
    daT(nn:101,jj) = nan;
    daWV(nn:101,jj) = nan;
  end
  daWV = daWV./(nanmean(daWV,2)*ones(1,length(boo)));
  numfound(ii) = length(boo);
  meanT(:,ii)  = nanmean(daT,2); 
  stdT(:,ii)   = nanstd(daT,[],2); 
  meanWV(:,ii) = nanmean(daWV,2); 
  stdWV(:,ii)  = nanstd(daWV,[],2); 
end

figure(5); pcolor(meanvaluebin(raLat),nanmean(p120.plevs,2),meanT); shading interp; colorbar; xlabel('Latitude'); ylabel('P [mb]'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('Mean T'); 
figure(6); pcolor(meanvaluebin(raLat),nanmean(p120.plevs,2),stdT);  shading interp; colorbar; xlabel('Latitude'); ylabel('P [mb]'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('Std T'); 
figure(6); colormap jet; set(gca,'fontsize',10)
figure(5); colormap jet; set(gca,'fontsize',10)

figure(7); pcolor(meanvaluebin(raLat),nanmean(p120.plevs,2),meanWV); shading interp; colorbar; xlabel('Latitude'); ylabel('P [mb]'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('Mean WV'); 
figure(8); pcolor(meanvaluebin(raLat),nanmean(p120.plevs,2),stdWV);  shading interp; colorbar; xlabel('Latitude'); ylabel('P [mb]'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('Std WV'); 
figure(7); colormap jet; set(gca,'fontsize',10)
figure(8); colormap jet; set(gca,'fontsize',10)

figure(9); errorbar(1:110,mean(pnlte.plevs,2),std(pnlte.plevs,[],2)); set(gca,'yscale','log')
  ylabel('Mean/std plevs'); xlabel('Z [km]')

addpath /home/sergio/MATLABCODE/EOF_PCA
figure(10);
boo = find(zgrid > 50*1000);
pah = pnlte.ptemp(boo,:);
L = driver_eof1(pah);
LL = L/sum(L); plot(1:length(LL),cumsum(LL),'.-')
  axis([0 10 0 1])

figure(11); subplot(121); 
  semilogy(pnlte.ptemp(boo,:),pnlte.plevs(boo,:)); set(gca,'ydir','reverse'); set(gca,'fontsize',10)
  title('Sergio')
  ax = axis;
figure(11); subplot(122); 
  boo = 1:floor(mean(i005));
  boo = 1:floor(mean(i1));
  semilogy(pcat.ptemp(boo,:),pcat.plevs(boo,:)); set(gca,'ydir','reverse'); set(gca,'fontsize',10)
  axis(ax);
  title('Chris')

%{
comment = 'see /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_120km/klayersV205_0_120km/Test_rtpV221/test_300_49.m';
comment = 'see /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_120km/klayersV205_0_120km/Test_rtpV221/driver_test_300_49_make_nlte_Z.m';
save pZ_0_120km.mat hnlte pnlte comment
%}
