addpath /home/sergio/MATLABCODE/COLORMAP

if ~exist('porig')
  sartaer = ['!date; time ../bin/airs_l1c_2834_cloudy_may19_prod_debug fin=cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.op.rtp fout=newdayx.rtp listp=1'];
  eval(sartaer);
  [h,ha,porig,pa] = rtpread('newdayx.rtp');
  fprintf(1,'nlevs = %3i \n',porig.nlevs)
  jacx = quicksartajac(h,porig,1,1,-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sartaer = ['!date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.op.rtp fout=newdayx.rtp listp=1 listj=-1'];
eval(sartaer);
[h,ha,pnew,pa] = rtpread('newdayx.rtp');

figure(1); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc)-rad2bt(h.vchan,pnew.rcalc)); xlim([640 1640]); title('BT DIFF orig-new')

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_WGTFCN',200);  dd = squeeze(d(1,:,:)); whos dd d
figure(2); clf; pcolor(h.vchan,1:97,dd(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('WGT ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacTZ',100);  dd = squeeze(d(1,:,:)); whos dd d
figure(1); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc)-rad2bt(h.vchan,pnew.rcalc),'kx-',h.vchan,jacx.jac(:,7),'b',h.vchan,dd(:,98),'r'); xlim([640 1640]); 
  hl = legend('BT DIFF orig-new','finite diff','analytic','location','best','fontsize',10);
figure(3); clf; pcolor(h.vchan,1:97,jacx.tjac(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ FINITE DIFF'); colormap jet; caxis([0 1]/10); xlim([640 1640])
figure(4); clf; pcolor(h.vchan,1:97,dd(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])

[w,d,iaProf,iaNumLay] = readsarta_jac('newdayx.rtp_jacG1',1);  dd = squeeze(d(1,:,:)); whos dd d
figure(5); clf; pcolor(h.vchan,1:97,jacx.wvjac(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 FINITE DIFF'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(6); clf; pcolor(h.vchan,1:97,dd(:,1:97)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
