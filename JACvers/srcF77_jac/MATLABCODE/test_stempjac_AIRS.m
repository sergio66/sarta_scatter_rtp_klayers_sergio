addpath /home/sergio/MATLABCODE_Git
addpath /home/sergio/MATLABCODE_Git/COLORMAP
addpath /home/sergio/MATLABCODE_Git/CONVERT_GAS_UNITS/

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
 iProf = -1;
end

rmer = ['!/bin/rm jactest.rtp_*']; eval(rmer);

sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=100'];
%sartaer = ['!ls -lt ' sarta_exec '; date; time ' sarta_exec ' fin=' frtp ' fout=jactest.rtp listp=' num2str(iProf)  ' listj=-1'];
eval(sartaer)

[h,ha,pnew,pa] = rtpread('jactest.rtp');

figure(1); clf; plot(h.vchan,rad2bt(h.vchan,pnew.rcalc)); xlim([640 1640]); title('BT')

if iProf > 0
  nlays = p0.nlevs(iProf)-1;
else
  nlays = p0.nlevs - 1;
end

if iProf > 0
  [w,d,iaProf,iaNumLay] = readsarta_jac('jactest.rtp_jacTZ',100);  
  for ii = 1 : length(p0.stemp)
    dTz = squeeze(d(ii,:,:)); % whos dTz d
    stempjac_sarta_fast(:,ii)   = dTz(:,nlays(ii)+1);
    ptempjac_sarta_fast(ii,:,:) = dTz(:,1:nlays(ii));
  end
  mmw = mmwater_rtp(h,p0);
  figure(2); [Y,I] = sort(p0.stemp,'desc'); plot(p0.stemp(I),stempjac_sarta_fast(1291,I)); xlabel('stemp'); ylabel('ST jac');
  figure(3); [Y,I] = sort(mmw,'desc');      plot(mmw(I),stempjac_sarta_fast(1291,I));      xlabel('mmw');   ylabel('ST jac');
    P = polyfit(mmw(I),stempjac_sarta_fast(1291,I),1); Y = polyval(P,mmw(I));
    plot(mmw(I),stempjac_sarta_fast(1291,I),mmw(I),Y);      xlabel('mmw');   ylabel('ST jac'); title(['dBT/dSKT = ' num2str(P(1)) ' mmw + ' num2str(P(2)) ])
end

%%%%%%%%%%%%%%%%%%%%%%%%%

