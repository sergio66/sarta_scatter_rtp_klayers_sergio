function [stempjac_sarta_fast,ptempjac_sarta_fast,g1jac_sarta_fast,cldjac_sarta_fast] = quick_sarta_jac(h,ha,porig,pa,iProf,sarta_exec,iCloudorClear)

addpath /home/sergio/MATLABCODE_Git
addpath /home/sergio/MATLABCODE_Git/COLORMAP

addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time with and without jacs

if nargin == 5
  sarta_exec = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod_debug_save3';
  sarta_exec = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod_debug';
  sarta_exec = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
  iCloudorClear = +1;
elseif nargin == 6
  iCloudorClear = +1;
end

fprintf(1,'sarta_exec = %s \n',sarta_exec);

nlays = porig.nlevs(iProf)-1;

frtp = 'jactest.op.rtp';

%iCloudorClear = +1;
%iCloudorClear = -1;
if iCloudorClear < 0
  disp('quick_sarta_jac.m : resetting clouds to zero!!')
  disp('quick_sarta_jac.m : resetting clouds to zero!!')
  porig.ctype  = -9999 * ones(size(porig.stemp));
  porig.ctype2 = -9999 * ones(size(porig.stemp));
  porig.cngwat  = -9999 * ones(size(porig.stemp));
  porig.cngwat2 = -9999 * ones(size(porig.stemp));
  porig.cfrac  = -9999 * ones(size(porig.stemp));
  porig.cfrac2 = -9999 * ones(size(porig.stemp));
  porig.cfrac12 = -9999 * ones(size(porig.stemp));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do finite diff jacs 

dT = 0.1; dQ = 0.001;
dT = 0.1; dQ = 0.100;
dT = 0.1; dQ = 0.010;

rtpwrite(frtp,h,ha,porig,pa);
sartaer = ['!ls -lt ' sarta_exec ';  date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)];
eval(sartaer);
[h,ha,pnew0,pa] = rtpread('jactest.rtp');

pjunk = porig;
pjunk.stemp = pjunk.stemp + dT;
rtpwrite(frtp,h,ha,pjunk,pa);
sartaer = ['!ls -lt ' sarta_exec ';  date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)];
eval(sartaer);
[h,ha,pST1,pa] = rtpread('jactest.rtp');

pjunk = porig;
pjunk.ptemp = pjunk.ptemp + dT;
rtpwrite(frtp,h,ha,pjunk,pa);
sartaer = ['!ls -lt ' sarta_exec ';  date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)];
eval(sartaer);
[h,ha,pT1,pa] = rtpread('jactest.rtp');

pjunk = porig;
pjunk.gas_1 = pjunk.gas_1 * (1 + dQ);
rtpwrite(frtp,h,ha,pjunk,pa);
sartaer = ['!ls -lt ' sarta_exec ';  date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)];
eval(sartaer);
[h,ha,pG1,pa] = rtpread('jactest.rtp');

%    case 1; cldstr = 'cfrac1';
%    case 2; cldstr = 'cngwat1';
%    case 3; cldstr = 'cpsize1';
%    case 4; cldstr = 'cprtop1';
%    case 5; cldstr = 'cprbot1';

[hnew,pjunk] = replicate_rtp_headprof(h,porig,iProf,5);
pjunk.cfrac(1) = min(pjunk.cfrac(1)*(1 + dQ),1);
pjunk.cngwat(2) = pjunk.cngwat(2)*(1 + dQ);
pjunk.cpsize(3) = pjunk.cpsize(3)*(1 + dQ);
pjunk.cprtop(4) = pjunk.cprtop(4)*(1 + dQ);
pjunk.cprbot(5) = pjunk.cprbot(5)*(1 + dQ);
rtpwrite(frtp,h,ha,pjunk,pa);
sartaer = ['!ls -lt ' sarta_exec ';  date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp >& ugh'];
eval(sartaer);
[h,ha,pCLD1,pa] = rtpread('jactest.rtp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% anallytic jacs

rtpwrite(frtp,h,ha,porig,pa);

sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=-1'];             %% ALL
sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=100'];            %% T
sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=100,200,1'];      %% T,WGT,WV
sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=100,200,300,1'];  %% T,WGT,CLD,WV

t3 = datetime("now");
eval(sartaer);
t4 = datetime("now");
fprintf(1,'time to run ANALYTIC JAC = %8.6f \n',etime(datevec(t4),datevec(t3)))
[h,ha,pnew,pa] = rtpread('jactest.rtp');

figure(1); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc(:,iProf))-rad2bt(h.vchan,pnew.rcalc)); xlim([640 1640]); title('BT DIFF orig-new')

[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_WGTFCN',200);  dd = squeeze(d(1,:,:)); % whos dd d
figure(2); clf; 
pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); 
title('WGT ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])
wgtfcn_sarta_fast = dd(:,1:nlays);
%%%%%%%%%%%%%%%%%%%%%%%%%

[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacTZ',100);  dd = squeeze(d(1,:,:)); % whos dd d
figure(1); clf; plot(h.vchan,dd(:,nlays+1),'r'); title('STEMP jac')
xlim([640 1640]);  

stempjac_sarta_fast = dd(:,nlays+1);
ptempjac_sarta_fast = dd(:,1:nlays);

figure(3); clf; 
pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); 
title('TZ ANALYTIC'); colormap jet; caxis([0 1]/10); xlim([640 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%

[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacG1',1);  dd = squeeze(d(1,:,:)); % whos dd d
g1jac_sarta_fast = dd(:,1:nlays);

figure(4); clf; 
pcolor(h.vchan,1:nlays,dd(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); 
title('G1 ANALYTIC'); colormap(usa2); caxis([-1 1]/2); xlim([640 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%
[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacCLD',300);  dd = squeeze(d(1,:,:)); % whos dd d
cldjac_sarta_fast = dd(:,1:12);

%% figure(21); 
%% for icld = 1 : 12
%%   switch icld
%%     case 1; cldstr = 'cfrac1';
%%     case 2; cldstr = 'cngwat1';
%%     case 3; cldstr = 'cpsize1';
%%     case 4; cldstr = 'cprtop1';
%%     case 5; cldstr = 'cprbot1';
%%     case 6; cldstr = 'cfrac2';
%%     case 7; cldstr = 'cngwat2';
%%     case 8; cldstr = 'cpsize2';
%%     case 9; cldstr = 'cprtop2';
%%     case 10; cldstr = 'cprbot2';
%%     case 11; cldstr = 'cfrac12';
%%     case 12; cldstr = 'stemp';
%%   end  
%%   if icld <= 11
%%     plot(h.vchan,jacx.cldjac(:,icld),'b.-',h.vchan,cldjac_sarta_fast(:,icld),'r'); title(cldstr);
%%   elseif icld == 12
%%     plot(h.vchan,jacx.jac(:,7),'b.-',h.vchan,cldjac_sarta_fast(:,icld),'r'); title(cldstr);
%%   end
%%   xlim([640 1640])
%%   pause(1.00);
%%   %pause
%% end
%% print_cloud_params(h,porig,1);

%%%%%%%%%%%%%%%%%%%%%%%%%

pause(0.5)
figure(1); clf
stdiff = (rad2bt(h.vchan,pST1.rcalc)-rad2bt(h.vchan,pnew0.rcalc))/dT;
Tdiff  = (rad2bt(h.vchan,pT1.rcalc) -rad2bt(h.vchan,pnew0.rcalc))/dT;
WVdiff = (rad2bt(h.vchan,pG1.rcalc)-rad2bt(h.vchan,pnew0.rcalc))/log(1 + dQ);

plot(h.vchan,stempjac_sarta_fast,'rx-',h.vchan,stdiff,'m',...
     h.vchan,nansum(ptempjac_sarta_fast,2),'gx-',h.vchan,Tdiff,'k',...
     h.vchan,nansum(g1jac_sarta_fast,2)/10,'bx-',h.vchan,WVdiff/10,'c','linewidth',2)
hl = legend('ST Analytic','ST FiniteDiff','T Analytic','T FiniteDiff','WV/10 Analytic','WV/10 FiniteDiff','location','best','fontsize',8);
plotaxis2;
if iCloudorClear > 0
  title(['Cloudy sky jacs for Profile ' num2str(iProf)]);
else
  title(['Clear sky jacs for Profile ' num2str(iProf)]);
end

figure(2); clf
Cdiff(:,1) = (rad2bt(h.vchan,pCLD1.rcalc(:,1))-rad2bt(h.vchan,pnew0.rcalc))/log(1 + dQ);
Cdiff(:,2) = (rad2bt(h.vchan,pCLD1.rcalc(:,2))-rad2bt(h.vchan,pnew0.rcalc))/log(1 + dQ);
Cdiff(:,3) = (rad2bt(h.vchan,pCLD1.rcalc(:,3))-rad2bt(h.vchan,pnew0.rcalc))/log(1 + dQ);
Cdiff(:,4) = (rad2bt(h.vchan,pCLD1.rcalc(:,4))-rad2bt(h.vchan,pnew0.rcalc))/log(1 + dQ);
Cdiff(:,5) = (rad2bt(h.vchan,pCLD1.rcalc(:,5))-rad2bt(h.vchan,pnew0.rcalc))/log(1 + dQ);
plot(h.vchan,cldjac_sarta_fast(:,1),'b',h.vchan,cldjac_sarta_fast(:,2),'g',h.vchan,cldjac_sarta_fast(:,3),'c',h.vchan,cldjac_sarta_fast(:,4),'k',h.vchan,cldjac_sarta_fast(:,5),'r',...
     h.vchan,Cdiff(:,1),'b--',h.vchan,Cdiff(:,2),'g--',h.vchan,Cdiff(:,3),'c--',h.vchan,Cdiff(:,4),'k--',h.vchan,Cdiff(:,5),'r--'); plotaxis2; 
  hl = legend('cfrac','cngwat','cpsize','cprtop','cprbot','location','best','fontsize',10);
  title('CLD jacs solid=analytic, dashed=finite diff')
rmer = ['!/bin/rm jactest.rtp* jactest.op.rtp']; eval(rmer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
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
%}

%%%%%%%%%%%%%%%%%%%%%%%%%
%{
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
%}
%%%%%%%%%%%%%%%%%%%%%%%%%
