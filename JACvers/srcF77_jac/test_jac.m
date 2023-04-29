addpath /home/sergio/MATLABCODE_Git
addpath /home/sergio/MATLABCODE_Git/COLORMAP

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time with and without jacs
sartaer = ['!ls -lt ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug; date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)];
eval(sartaer);
sartaer = ['!ls -lt ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug; date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=' frtp ' fout=newdayx.rtp listp=' num2str(iProf)  ' listj=-1'];
eval(sartaer);
[h,ha,pnew,pa] = rtpread('newdayx.rtp');

figure(1); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc)-rad2bt(h.vchan,pnew.rcalc)); xlim([640 1640]); title('BT DIFF orig-new')

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_WGTFCN',200);  dd = squeeze(d(1,:,:)); % whos dd d
figure(2); clf; pcolor(h.vchan,1:97,dd(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('WGT ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacTZ',100);  dd = squeeze(d(1,:,:)); % whos dd d
figure(1); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc)-rad2bt(h.vchan,pnew.rcalc),'kx-',h.vchan,jacx.jac(:,7),'b',h.vchan,dd(:,98),'r'); xlim([640 1640]); 
  hl = legend('BT DIFF orig-new','finite diff','analytic','location','best','fontsize',10);
figure(3); clf; pcolor(h.vchan,1:97,jacx.tjac(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ FINITE DIFF'); colormap jet; caxis([0 1]/10); xlim([640 1640])
figure(4); clf; pcolor(h.vchan,1:97,dd(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])
figure(9); clf; 
  figure(9); hold on; plot(h.vchan,10*sum(jacx.tjac(:,1:97)'),'r.-',h.vchan,10*nansum(dd(:,1:97)'),'m'); hold on
figure(10); clf; 
  figure(10); hold on; plot(h.vchan,sum(jacx.tjac(:,1:97)')./nansum(dd(:,1:97)'),'b'); hold on; axis([640 1640 0 +2]); 

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacG1',1);  dd = squeeze(d(1,:,:)); % whos dd d
figure(5); clf; pcolor(h.vchan,1:97,jacx.wvjac(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 FINITE DIFF'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(6); clf; pcolor(h.vchan,1:97,dd(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(9); hold on; plot(h.vchan,sum(jacx.wvjac(:,1:97)'),'b.-',h.vchan,sum(dd(:,1:97)'),'c'); hold on; 
  axis([640 1640 -10 +10])
  hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','location','best','fontsize',10);
figure(10); hold on; plot(h.vchan,sum(jacx.wvjac(:,1:97)')./nansum(dd(:,1:97)'),'r'); hold on; axis([640 1640 0 +2]); 

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacG3',3);  dd = squeeze(d(1,:,:)); % whos dd d
figure(7); clf; pcolor(h.vchan,1:97,jacx.o3jac(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 FINITE DIFF'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(8); clf; pcolor(h.vchan,1:97,dd(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(9); hold on; plot(h.vchan,sum(jacx.o3jac(:,1:97)'),'k.-',h.vchan,sum(dd(:,1:97)'),'g'); hold off; 
  axis([640 1640 -20 +20])
  hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','sum(O3jac),finite diff','sum(O3jac),analytic','location','best','fontsize',10);
figure(10); hold on; plot(h.vchan,sum(jacx.o3jac(:,1:97)')./nansum(dd(:,1:97)'),'k'); hold on; axis([640 1640 0 +3]); title('ratio FINITE/ANALYTIC'); hl = legend('TZ','WV','O3','location','best','fontsize',10);

figure(9); hold off
figure(10); hold off
