addpath /home/sergio/MATLABCODE_Git
addpath /home/sergio/MATLABCODE_Git/COLORMAP

addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools

frtp = '/asl/s1/sergio/RTP_pin_feb2002/pin_feb2002_sea_airsnadir_g80_op.so2.rtp';
frtp = '/asl/s1/sergio/RTP_pin_feb2002/pin_feb2002_sea_airsnadir_op.so2.latlon.rtp';
[h,ha,p0,pa] = rtpread(frtp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time with and without jacs

sarta_exec = '../../bin/jac_airs_l1c_2834_cloudy_may19_prod';

fprintf(1,'sarta_exec = %s \n',sarta_exec);

iProf = input('Enter Profile number : ');
if length(iProf) == 0
 iProf = 1;
end

rmer = ['!/bin/rm jactest.rtp_*']; eval(rmer);

sartaer = ['!ls -lt ' sarta_exec ';  date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)];
eval(sartaer);
sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=100'];
sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=-1'];
eval(sartaer)

[h,ha,pnew,pa] = rtpread('jactest.rtp');

figure(1); clf; plot(h.vchan,rad2bt(h.vchan,pnew.rcalc)); xlim([640 1640]); title('BT')

nlays = p0.nlevs(iProf)-1;

[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_WGTFCN',200);  dWGT = squeeze(d(1,:,:)); % whos dWGT d
figure(2); clf; pcolor(h.vchan,1:nlays,dWGT(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('WGT ANALYTIC'); colormap jet; xlim([640 1640])
wgtfcn_sarta_fast = dWGT(:,1:nlays);
caxis([0 1]/10); 
%%%%%%%%%%%%%%%%%%%%%%%%%

[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacTZ',100);  dTz = squeeze(d(1,:,:)); % whos dTz d
stempjac_sarta_fast = dTz(:,nlays+1);
ptempjac_sarta_fast = dTz(:,1:nlays);
figure(3); clf; pcolor(h.vchan,1:nlays,dTz(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('TZ ANALYTIC'); colormap jet; xlim([640 1640])
caxis([0 1]/10); 
%%%%%%%%%%%%%%%%%%%%%%%%%

[w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacG1',1);  dG1 = squeeze(d(1,:,:)); % whos dG1 d
g1jac_sarta_fast = dG1(:,1:nlays);
figure(4); clf; pcolor(h.vchan,1:nlays,dG1(:,1:nlays)'); colorbar; shading interp; set(gca,'ydir','reverse'); title('G1 ANALYTIC'); colormap(usa2); xlim([640 1640])
caxis([-1 1]/2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0.plays = plevs2plays(p0.plevs);
for ii = 1 : h.nchan
  moo = dWGT(ii,:);
  miaow = find(moo == max(moo),1);
  wgtpeak(ii) = miaow;
  Ppeak(ii) = p0.plays(miaow,iProf);
  Tpeak(ii) = p0.ptemp(miaow,iProf);
end
figure(5); plot(h.vchan,rad2bt(h.vchan,pnew.rcalc),'b',h.vchan,Tpeak,'r')
  hl = legend('rad2bt(v,rad)','T(z) at wgt peak','location','best','fontsize',10);
  title(['Profile ' num2str(iProf) ' of 49'])
figure(6); plot(h.vchan,Ppeak); set(gca,'ydir','reverse'); title('Wgt Function Peak');
