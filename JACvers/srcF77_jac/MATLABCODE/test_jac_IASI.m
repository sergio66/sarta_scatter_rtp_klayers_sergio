addpath /home/sergio/MATLABCODE_Git
addpath /home/sergio/MATLABCODE_Git/COLORMAP

addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools

frtp = 'TEST_RTP/mktemp_N2dzVSCD_airs_001_kl.rtp_1'; iProf = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('porig')
  rmer = ['!/bin/rm jactest.rtp jactest.rtp_jac* jactest.rtp_WGTFCN ']; eval(rmer);

  iExtraGas = [];
  iExtraGas = [2 4 5 6 9 12];
  iExtraGas = [9];
  iExtraGas = input('Enter which gas ID to check [2 4 5 6 9 12] : ');
  if length(iExtraGas) == 0
    iExtraGas = 2;
  end

  disp('MAKE SURE THE quicksartajac FINITE JACS use < dQ = 0.01;  dT = 0.1 > ... if you make them smaller, perversely the accuracy gets worse!!!!')
  disp('MAKE SURE THE quicksartajac FINITE JACS use < dQ = 0.01;  dT = 0.1 > ... if you make them smaller, perversely the accuracy gets worse!!!!')
  disp('MAKE SURE THE quicksartajac FINITE JACS use < dQ = 0.01;  dT = 0.1 > ... if you make them smaller, perversely the accuracy gets worse!!!!')

  sartaer = ['!date; time /asl/packages/sartaV108/BinV201/sarta_iasi_may09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4  fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)];
  eval(sartaer);
  [h,ha,porig,pa] = rtpread('jactest.rtp');
  fprintf(1,'nlevs = %3i satzen = %8.6f solzen = %8.6f\n',porig.nlevs,porig.satzen,porig.solzen)
  %jacx = quicksartajac(h,porig,1,3,-1);
  t1 = datetime("now");
  jacx = quicksartajac(h,porig,1,3,-1,iExtraGas);
  t2 = datetime("now");
end
nlays = porig.nlevs-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time with and without jacs

sarta_exec = '../bin/jac_sarta_iasi_may09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4';

fprintf(1,'sarta_exec = %s \n',sarta_exec);

sartaer = ['!ls -lt ' sarta_exec ';  date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)];
eval(sartaer);
sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=100'];
sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=-1'];
t3 = datetime("now");
eval(sartaer);
t4 = datetime("now");
fprintf(1,'time to run FINITE DIFF vs ANALYTIC = %8.6f %8.6f \n',etime(datevec(t2),datevec(t1)),etime(datevec(t4),datevec(t3)))
[h,ha,pnew,pa] = rtpread('jactest.rtp');

figure(1); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc)-rad2bt(h.vchan,pnew.rcalc)); xlim([640 1640]); title('BT DIFF orig-new')

[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_WGTFCN',200);  dd = squeeze(d(1,:,:)); % whos dd d
figure(2); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('WGT ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])
wgtfcn_sarta_fast = dd(:,1:nlays);
%%%%%%%%%%%%%%%%%%%%%%%%%

[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacTZ',100);  dd = squeeze(d(1,:,:)); % whos dd d
figure(1); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc)-rad2bt(h.vchan,pnew.rcalc),'kx-',h.vchan,jacx.jac(:,7),'b',h.vchan,dd(:,nlays+1),'r'); xlim([640 1640]); 
  hl = legend('BT DIFF orig-new','STEMP JAC finite diff','STEMP JACanalytic','location','best','fontsize',10);
stempjac_sarta_fast = dd(:,nlays+1);
ptempjac_sarta_fast = dd(:,1:nlays);

figure(3); clf; pcolor(h.vchan,1:nlays,jacx.tjac(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ FINITE DIFF'); colormap jet; caxis([0 1]/10); xlim([640 1640])
figure(4); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])
figure(5); clf; pcolor(h.vchan,1:nlays,jacx.tjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ FINITE/ANALYTIC'); colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])
ratio = jacx.tjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)'); bad = find(abs(jacx.tjac(:,1:nlays)') < 1e-6); ratio(bad) = nan;
figure(5); clf; pcolor(h.vchan,1:nlays,ratio-1); colorbar; shading interp; set(gca,'ydir','reverse'); title('FINITE/ANALYTIC TZ - 1'); colormap(usa2); caxis([-1 +1]); xlim([640 1640])

figure(15); clf; 
  figure(15); hold on; plot(h.vchan,10*sum(jacx.tjac(:,1:nlays)'),'r.-',h.vchan,10*sum(dd(:,1:nlays)'),'m'); hold on
figure(16); clf; 
  figure(16); hold on; plot(h.vchan,(eps+sum(jacx.tjac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'r'); hold on; axis([640 1640 0 +2]); 

%%%%%%%%%%%%%%%%%%%%%%%%%

[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacG1',1);  dd = squeeze(d(1,:,:)); % whos dd d
g1jac_sarta_fast = dd(:,1:nlays);
figure(6); clf; pcolor(h.vchan,1:nlays,jacx.wvjac(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 FINITE DIFF'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(7); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(8); clf; pcolor(h.vchan,1:nlays,jacx.wvjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 FINITE/ANALYTIC'); colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])
ratio = jacx.wvjac(:,1:nlays)' ./ (eps + dd(:,1:nlays)'); bad = find(abs(jacx.wvjac(:,1:nlays)') < 1e-6); ratio(bad) = nan;
figure(8); clf; pcolor(h.vchan,1:nlays,ratio-1); colorbar; shading interp; set(gca,'ydir','reverse'); title('FINITE/ANALYTIC G1 - 1'); colormap(usa2); caxis([-1 +1]); xlim([640 1640])

figure(15); hold on; plot(h.vchan,sum(jacx.wvjac(:,1:nlays)'),'b.-',h.vchan,sum(dd(:,1:nlays)'),'c'); hold on; 
  axis([640 1640 -10 +10])
  hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','location','best','fontsize',10);
figure(16); hold on; plot(h.vchan,(eps+sum(jacx.wvjac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'b'); hold on; axis([640 1640 0 +2]); 

%%%%%%%%%%%%%%%%%%%%%%%%%
[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacG3',3);  dd = squeeze(d(1,:,:)); % whos dd d
g3jac_sarta_fast = dd(:,1:nlays);
figure(9); clf; pcolor(h.vchan,1:nlays,jacx.o3jac(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 FINITE DIFF'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(10); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
figure(11); clf; pcolor(h.vchan,1:nlays,jacx.o3jac(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title('G3 FINITE/ANALYTIC'); colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])
ratio = jacx.o3jac(:,1:nlays)' ./ (eps + dd(:,1:nlays)'); bad = find(abs(jacx.o3jac(:,1:nlays)') < 1e-6); ratio(bad) = nan;
figure(11); clf; pcolor(h.vchan,1:nlays,ratio-1); colorbar; shading interp; set(gca,'ydir','reverse'); title('FINITE/ANALYTIC G3 - 1'); colormap(usa2); caxis([-1 +1]); xlim([640 1640])

figure(15); hold on; plot(h.vchan,sum(jacx.o3jac(:,1:nlays)'),'k.-',h.vchan,sum(dd(:,1:nlays)'),'g'); hold off; 
  axis([640 1640 -20 +20])
  hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','sum(O3jac),finite diff','sum(O3jac),analytic','location','best','fontsize',8);
figure(16); hold on; plot(h.vchan,(eps+sum(jacx.o3jac(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'k'); hold on; axis([640 1640 0 +3]); title('ratio FINITE/ANALYTIC'); 
hl = legend('TZ','WV','O3','location','best','fontsize',10);

figure(17); clf; plot(h.vchan,jacx.jac(:,1),'b.-',h.vchan,sum(g1jac_sarta_fast,2),'r')
   title('Comparing G1 (WV) jacs'); hl = legend('SARTA Finite','SARTA analytic','location','best','fontsize',10); xlim([1000 1500])
   xlabel('Wavenumber cm-1'); ylabel('sum( dBT/dX X)')
figure(18); clf; plot(h.vchan,jacx.jac(:,3),'b.-',h.vchan,sum(g3jac_sarta_fast,2),'r')
   title('Comparing G3 (O3) jacs'); hl = legend('SARTA Finite','SARTA analytic','location','best','fontsize',10); xlim([900 1200])
   xlabel('Wavenumber cm-1'); ylabel('sum( dBT/dX X)')
figure(19); clf; plot(h.vchan,jacx.jac(:,8),'b.-',h.vchan,sum(ptempjac_sarta_fast,2),'r')
   title('Comparing TZ jacs'); hl = legend('SARTA Finite','SARTA analytic','location','best','fontsize',10); xlim([650 1000])
   xlabel('Wavenumber cm-1'); ylabel('sum( dBT/dX X)')

%%%%%%%%%%%%%%%%%%%%%%%%%
if length(iExtraGas) == 1
  [w,d,iaProf,iaNumLay] = readsarta_jac(['jactest.rtp_jacG' num2str(iExtraGas)],iExtraGas);  dd = squeeze(d(1,:,:)); % whos dd d
  str = ['gNjac_sarta_fast = jacx.G' num2str(iExtraGas) 'jac;']; 
  eval(str);
  figure(12); clf; pcolor(h.vchan,1:nlays,gNjac_sarta_fast(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title(['G' num2str(iExtraGas) ' FINITE DIFF']); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
  figure(13); clf; pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title(['G' num2str(iExtraGas) ' ANALYTIC']); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])
  figure(14); clf; pcolor(h.vchan,1:nlays,gNjac_sarta_fast(:,1:nlays)' ./ (eps + dd(:,1:nlays)')); colorbar; shading interp; set(gca,'ydir','reverse'); title(['G' num2str(iExtraGas) ' FINITE/ANALYTIC']); 
    colormap(usa2); caxis([-1 +1]*2); xlim([640 1640])
  ratio = gNjac_sarta_fast(:,1:nlays)' ./ (eps + dd(:,1:nlays)'); bad = find(abs(gNjac_sarta_fast(:,1:nlays)') < 1e-6); ratio(bad) = nan;
  figure(14); clf; pcolor(h.vchan,1:nlays,ratio-1); colorbar; shading interp; set(gca,'ydir','reverse'); 
  title(['FINITE/ANALYTIC G' num2str(iExtraGas) ' - 1']); colormap(usa2); caxis([-1 +1]); xlim([640 1640])

  if iExtraGas == 5
    figure(12); xlim([2000 2700])
    figure(13); xlim([2000 2700])
    figure(14); xlim([2000 2700])
  end
  
  strG = ['G' num2str(iExtraGas)];
  figure(15); hold on; plot(h.vchan,sum(gNjac_sarta_fast(:,1:nlays)'),'k.-',h.vchan,sum(dd(:,1:nlays)'),'g'); hold off; 
    axis([640 1640 -20 +20])
    hl = legend('10*sum(TZjac),finite diff','10*sum(TZjac),analytic','sum(WVjac),finite diff','sum(WVjac),analytic','sum(O3jac),finite diff','sum(O3jac),analytic',...
                ['sum(G' num2str(iExtraGas) 'jac),finitediff'],['sum(G' num2str(iExtraGas) 'jac),analytic'],'location','best','fontsize',8);  
  figure(16); hold on; plot(h.vchan,(eps+sum(gNjac_sarta_fast(:,1:nlays)'))./(eps+sum(dd(:,1:nlays)')),'g'); hold on; axis([640 1640 0 +3]); 
  figure(16); hold on; plot(h.vchan,jacx.jac(:,1)'./nansum(jacx.wvjac(:,1:nlays)'),'c'); 
     title('ratio FINITE/ANALYTIC \newline (cyan=column/sum(lay) WV)'); 
     hl = legend('TZ','WV','O3',strG,'finite WV jac column/sumlay','location','best','fontsize',10);
  
  if iExtraGas == 2
    figure(20); clf; plot(h.vchan,jacx.jac(:,2),'b.-',h.vchan,sum(gNjac_sarta_fast,2),'r')
     title('Comparing G2 (CO2) jacs'); hl = legend('SARTA Finite','SARTA analytic','location','best','fontsize',10); xlim([640 1000])
     xlabel('Wavenumber cm-1'); ylabel('sum( dBT/dX X)')
  elseif iExtraGas == 4
    figure(20); clf; plot(h.vchan,jacx.jac(:,4),'b.-',h.vchan,sum(gNjac_sarta_fast,2),'r')
     title('Comparing G4 (N2O) jacs'); hl = legend('SARTA Finite','SARTA analytic','location','best','fontsize',10); xlim([1000 1500])
     xlabel('Wavenumber cm-1'); ylabel('sum( dBT/dX X)')
  elseif iExtraGas == 5
    figure(20); clf; plot(h.vchan,jacx.jac(:,5),'b.-',h.vchan,sum(gNjac_sarta_fast,2),'r')
     title('Comparing G5 (CO) jacs'); hl = legend('SARTA Finite','SARTA analytic','location','best','fontsize',10); xlim([2150 2250])
     xlabel('Wavenumber cm-1'); ylabel('sum( dBT/dX X)')
  elseif iExtraGas == 6
    figure(20); clf; plot(h.vchan,jacx.jac(:,6),'b.-',h.vchan,sum(gNjac_sarta_fast,2),'r')
     title('Comparing G6 (CH4) jacs'); hl = legend('SARTA Finite','SARTA analytic','location','best','fontsize',10); xlim([1000 1500])
     xlabel('Wavenumber cm-1'); ylabel('sum( dBT/dX X)')
  end

end
%%%%%%%%%%%%%%%%%%%%%%%%%
[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacCLD',300);  dd = squeeze(d(1,:,:)); % whos dd d
cldjac_sarta_fast = dd(:,1:12);
figure(21); 
for icld = 1 : 12
  switch icld
    case 1; cldstr = 'cfrac1';
    case 2; cldstr = 'cngwat1';
    case 3; cldstr = 'cpsize1';
    case 4; cldstr = 'cprtop1';
    case 5; cldstr = 'cprbot1';
    case 6; cldstr = 'cfrac2';
    case 7; cldstr = 'cngwat2';
    case 8; cldstr = 'cpsize2';
    case 9; cldstr = 'cprtop2';
    case 10; cldstr = 'cprbot2';
    case 11; cldstr = 'cfrac12';
    case 12; cldstr = 'stemp';
  end  
  if icld <= 11
    plot(h.vchan,jacx.cldjac(:,icld),'b.-',h.vchan,cldjac_sarta_fast(:,icld),'r'); title(cldstr);
  elseif icld == 12
    plot(h.vchan,jacx.jac(:,7),'b.-',h.vchan,cldjac_sarta_fast(:,icld),'r'); title(cldstr);
  end
  xlim([640 1640])
  pause(1.0);
  %pause
end
print_cloud_params(h,porig,1);

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(15); hold off
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
  figure(15); figure(16); figure(17); figure(18); figure(19); figure(20); pause(1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
handO3jac = load('o3_jac_column_sum.mat');
i1041 = find(h.vchan >= 1041,1);
figure(20); clf;
plot(jacx.o3jac(i1041,1:97),1:97,'b',g3jac_sarta_fast(i1041,1:97),1:97,'r',handO3jac.jacBT,1:97,'k'); set(gca,'ydir','reverse')
hl = legend('SARTA finite diff','SARTA analystic','hand jac','location','best');
[sum(jacx.o3jac(i1041,1:97)) sum(g3jac_sarta_fast(i1041,1:97)) sum(handO3jac.jacBT)]
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
