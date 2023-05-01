addpath /home/sergio/MATLABCODE_Git
addpath /home/sergio/MATLABCODE_Git/COLORMAP

addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools

frtp = 'cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.op.rtp';
iProf = 45; %% satzen =  0  deg
iProf = 1;  %% satzen = -50 deg

if ~exist('porig')
  sartaer = ['!date; time ../bin/airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)];
  eval(sartaer);
  [h,ha,porig,pa] = rtpread('newdayx.rtp');
  fprintf(1,'nlevs = %3i scanang = %8.6f \n',porig.nlevs,porig.satzen)
  jacx = quicksartajac(h,porig,1,1,-1);
end
nlays = porig.nlevs-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time with and without jacs
sartaer = ['!ls -lt ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug; date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)];
eval(sartaer);
sartaer = ['!ls -lt ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug; date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)  ' listj=100'];
sartaer = ['!ls -lt ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug; date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)  ' listj=-1'];
eval(sartaer);
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

figure(15); clf; 
  figure(15); hold on; plot(h.vchan,10*sum(jacx.tjac(:,1:nlays)'),'r.-',h.vchan,10*sum(dd(:,1:nlays)'),'m'); hold on
figure(16); clf; 
  figure(16); hold on; plot(h.vchan,(eps+sum(jacx.tjac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'b'); hold on; axis([640 1640 0 +2]); 

%%%%%%%%%%%%%%%%%%%%%%%%%

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacG1',1);  dd = squeeze(d(1,:,:)); % whos dd d
figure(6); clf; pcolor(h.vchan,1:nlays,jacx.wvjac(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 FINITE DIFF'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(7); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(8); clf; pcolor(h.vchan,1:nlays,jacx.wvjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 FINITE/ANALYTIC'); colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])

figure(15); hold on; plot(h.vchan,sum(jacx.wvjac(:,1:nlays)'),'b.-',h.vchan,sum(dd(:,1:nlays)'),'c'); hold on; 
  axis([640 1640 -10 +10])
  hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','location','best','fontsize',10);
figure(16); hold on; plot(h.vchan,(eps+sum(jacx.wvjac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'r'); hold on; axis([640 1640 0 +2]); 

%%%%%%%%%%%%%%%%%%%%%%%%%
[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacG3',3);  dd = squeeze(d(1,:,:)); % whos dd d
figure(9); clf; pcolor(h.vchan,1:nlays,jacx.o3jac(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 FINITE DIFF'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(10); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(11); clf; pcolor(h.vchan,1:nlays,jacx.o3jac(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 FINITE/ANALYTIC'); colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])

figure(15); hold on; plot(h.vchan,sum(jacx.o3jac(:,1:nlays)'),'k.-',h.vchan,sum(dd(:,1:nlays)'),'g'); hold off; 
  axis([640 1640 -20 +20])
  hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','sum(O3jac),finite diff','sum(O3jac),analytic','location','best','fontsize',8);
figure(16); hold on; plot(h.vchan,(eps+sum(jacx.o3jac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'k'); hold on; axis([640 1640 0 +3]); title('ratio FINITE/ANALYTIC'); hl = legend('TZ','WV','O3','location','best','fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(15); hold off
figure(16); hold off
