addpath /home/sergio/MATLABCODE_Git
addpath /home/sergio/MATLABCODE_Git/COLORMAP

addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools

frtp = 'cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.op.rtp'; iProf = 45; %% satzen =  0  deg
frtp = 'cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.op.rtp'; iProf = 1;  %% satzen = -50  deg
frtp = 'newdayx_1_100_12150.op.rtp'; iProf = 59;    %% DCC
frtp = 'cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.op.rtp'; iProf = 12150/2+45;  %% satzen = +50  deg
frtp = 'newdayx_1_100_12150.op.rtp'; iProf = 113;   %% almost clear
frtp = 'newdayx_1_100_12150.op.rtp'; iProf = 1;     %% same as profile 1 from cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.op.rtp

if ~exist('porig')

  iExtraGas = [];
  iExtraGas = [2 4 5 6 9 12];
  iExtraGas = [9];
  iExtraGas = input('Enter which gas ID to check [2 4 5 6 9 12] : ');
  if length(iExtraGas) == 0
    iExtraGas = 2;
  end

  sartaer = ['!date; time ../bin/airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)];
  eval(sartaer);
  [h,ha,porig,pa] = rtpread('newdayx.rtp');
  fprintf(1,'nlevs = %3i scanang = %8.6f \n',porig.nlevs,porig.satzen)
  %jacx = quicksartajac(h,porig,1,1,-1);
  t1 = datetime("now");
  jacx = quicksartajac(h,porig,1,1,-1,iExtraGas);
  t2 = datetime("now");
end
nlays = porig.nlevs-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time with and without jacs
sartaer = ['!ls -lt ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug; date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)];
eval(sartaer);
sartaer = ['!ls -lt ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug; date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)  ' listj=100'];
sartaer = ['!ls -lt ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug; date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)  ' listj=-1'];
t3 = datetime("now");
eval(sartaer);
t4 = datetime("now");
fprintf(1,'time to run FINITE DIFF vs ANALYTIC = %8.6f %8.6f \n',etime(datevec(t2),datevec(t1)),etime(datevec(t4),datevec(t3)))
[h,ha,pnew,pa] = rtpread('newdayx.rtp');

figure(1); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc)-rad2bt(h.vchan,pnew.rcalc)); xlim([640 1640]); title('BT DIFF orig-new')

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_WGTFCN',200);  dd = squeeze(d(1,:,:)); % whos dd d
figure(2); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('WGT ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacTZ',100);  dd = squeeze(d(1,:,:)); % whos dd d
figure(1); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc)-rad2bt(h.vchan,pnew.rcalc),'kx-',h.vchan,jacx.jac(:,7),'b',h.vchan,dd(:,nlays+1),'r'); xlim([640 1640]); 
  hl = legend('BT DIFF orig-new','finite diff','analytic','location','best','fontsize',10);

figure(3); clf; pcolor(h.vchan,1:nlays,jacx.tjac(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ FINITE DIFF'); colormap jet; caxis([0 1]/10); xlim([640 1640])
figure(4); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])
figure(5); clf; pcolor(h.vchan,1:nlays,jacx.tjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ FINITE/ANALYTIC'); colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])
ratio = jacx.tjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)'); bad = find(abs(jacx.tjac(:,1:nlays)') < 1e-6); ratio(bad) = nan;
figure(5); clf; pcolor(h.vchan,1:nlays,ratio-1); colorbar; shading interp; set(gca,'ydir','reverse'); title('FINITE/ANALYTIC TZ-1'); colormap(usa2); caxis([-1 +1]); xlim([640 1640])

figure(15); clf; 
  figure(15); hold on; plot(h.vchan,10*sum(jacx.tjac(:,1:nlays)'),'r.-',h.vchan,10*sum(dd(:,1:nlays)'),'m'); hold on
figure(16); clf; 
  figure(16); hold on; plot(h.vchan,(eps+sum(jacx.tjac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'r'); hold on; axis([640 1640 0 +2]); 

%%%%%%%%%%%%%%%%%%%%%%%%%

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacG1',1);  dd = squeeze(d(1,:,:)); % whos dd d
figure(6); clf; pcolor(h.vchan,1:nlays,jacx.wvjac(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 FINITE DIFF'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(7); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(8); clf; pcolor(h.vchan,1:nlays,jacx.wvjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 FINITE/ANALYTIC'); colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])
ratio = jacx.wvjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)'); bad = find(abs(jacx.wvjac(:,1:nlays)') < 1e-6); ratio(bad) = nan;
figure(8); clf; pcolor(h.vchan,1:nlays,ratio-1); colorbar; shading interp; set(gca,'ydir','reverse'); title('FINITE/ANALYTIC G1-1'); colormap(usa2); caxis([-1 +1]); xlim([640 1640])

figure(15); hold on; plot(h.vchan,sum(jacx.wvjac(:,1:nlays)'),'b.-',h.vchan,sum(dd(:,1:nlays)'),'c'); hold on; 
  axis([640 1640 -10 +10])
  hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','location','best','fontsize',10);
figure(16); hold on; plot(h.vchan,(eps+sum(jacx.wvjac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'b'); hold on; axis([640 1640 0 +2]); 

%%%%%%%%%%%%%%%%%%%%%%%%%
[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacG3',3);  dd = squeeze(d(1,:,:)); % whos dd d
figure(9); clf; pcolor(h.vchan,1:nlays,jacx.o3jac(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 FINITE DIFF'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(10); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(11); clf; pcolor(h.vchan,1:nlays,jacx.o3jac(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 FINITE/ANALYTIC'); colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])
ratio = jacx.o3jac(:,1:nlays)' ./ (eps + dd(:,1:nlays)'); bad = find(abs(jacx.o3jac(:,1:nlays)') < 1e-6); ratio(bad) = nan;
figure(11); clf; pcolor(h.vchan,1:nlays,ratio-1); colorbar; shading interp; set(gca,'ydir','reverse'); title('FINITE/ANALYTIC G3-1'); colormap(usa2); caxis([-1 +1]); xlim([640 1640])

figure(15); hold on; plot(h.vchan,sum(jacx.o3jac(:,1:nlays)'),'k.-',h.vchan,sum(dd(:,1:nlays)'),'g'); hold off; 
  axis([640 1640 -20 +20])
  hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','sum(O3jac),finite diff','sum(O3jac),analytic','location','best','fontsize',8);
figure(16); hold on; plot(h.vchan,(eps+sum(jacx.o3jac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'k'); hold on; axis([640 1640 0 +3]); title('ratio FINITE/ANALYTIC'); hl = legend('TZ','WV','O3','location','best','fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%
if length(iExtraGas) == 1
  [w,d,iaProf,iaNumLay] = readsarta_jac(['newdayx.rtp_jacG' num2str(iExtraGas)],iExtraGas);  dd = squeeze(d(1,:,:)); % whos dd d
  str = ['gNjac = jacx.G' num2str(iExtraGas) 'jac;']; 
  eval(str);
  figure(12); clf; pcolor(h.vchan,1:nlays,gNjac(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title(['G' num2str(iExtraGas) ' FINITE DIFF']); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
  figure(13); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title(['G' num2str(iExtraGas) ' ANALYTIC']); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
  figure(14); clf; pcolor(h.vchan,1:nlays,gNjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title(['G' num2str(iExtraGas) ' FINITE/ANALYTIC']); 
    colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])
  ratio = gNjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)'); bad = find(abs(gNjac(:,1:nlays)') < 1e-6); ratio(bad) = nan;
  figure(14); clf; pcolor(h.vchan,1:nlays,ratio-1); colorbar; shading interp; set(gca,'ydir','reverse'); 
  title(['FINITE/ANALYTIC G' num2str(iExtraGas) '-1']); colormap(usa2); caxis([-1 +1]); xlim([640 1640])

  if iExtraGas == 5
    figure(12); xlim([2000 2700])
    figure(13); xlim([2000 2700])
    figure(14); xlim([2000 2700])
  end
  
  strG = ['G' num2str(iExtraGas)];
  figure(15); hold on; plot(h.vchan,sum(gNjac(:,1:nlays)'),'k.-',h.vchan,sum(dd(:,1:nlays)'),'g'); hold off; 
    axis([640 1640 -20 +20])
    hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','sum(O3jac),finite diff','sum(O3jac),analytic',...
                ['sum(G' num2str(iExtraGas) 'jac),finitediff'],['sum(G' num2str(iExtraGas) 'jac),analytic'],'location','best','fontsize',8);  
  figure(16); hold on; plot(h.vchan,(eps+sum(gNjac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'g'); hold on; axis([640 1640 0 +3]); title('ratio FINITE/ANALYTIC \newline fg(cyan=column/sum(lay) WV)'); 
     hl = legend('TZ','WV','O3',strG,'location','best','fontsize',10);
end
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(15); hold off
figure(16); hold on; plot(h.vchan,jacx.jac(:,1)'./nansum(jacx.wvjac(:,1:nlays)'),'c'); 
figure(16); hold off

iAX = input('Change x-axis to [640 2780] ? (-1 default /+1) : ');
if length(iAX) == 0
  iAX = -1;
end
if iAX > 0
  for ii = 1 : 16
    figure(ii); xlim([640 2800]);
  end
  figure(01); figure(02);             pause(1)
  figure(03); figure(04); figure(05); pause(1)
  figure(06); figure(07); figure(08); pause(1)
  figure(09); figure(10); figure(11); pause(1)
  figure(12); figure(13); figure(14); pause(1)
end
