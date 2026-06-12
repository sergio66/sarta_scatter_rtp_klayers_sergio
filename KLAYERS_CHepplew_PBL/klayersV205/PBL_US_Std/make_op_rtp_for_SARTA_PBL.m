addpath /home/sergio/git/matlabcode
addpath /home/sergio/git/matlabcode/CONVERT_GAS_UNITS
addpath /home/sergio/git/matlabcode/matlibSergio/matlib2025/h4tools
addpath /home/sergio/git/matlabcode/matlibSergio/matlib2025/rtptools

% this is to compare CO2 ppmv and make potential many 370,385,400 ppmv for SARTA breakouts

klayers = '../BinV221/klayers_pbl_wetwater_test';
co2 = 370;
co2 = 385;
co2 = 400;

fip0 = ['/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/regr49_1100_with_co2_' num2str(co2) 'ppm_9gases_unitemiss.ip.rtp'];
fipX = ['us_std_for_pbl_breakouts' num2str(co2) '.ip.rtp'];
fopX = ['us_std_for_pbl_breakouts' num2str(co2) '.op.rtp'];

[h,ha,p,pa] = rtpread(fip0);
%[h,p] = subset_rtp(h,p,[],[],49);

rtpwrite(fipX,h,ha,p,pa);

klayerser = ['!' klayers ' fin=' fipX ' fout=' fopX];
eval(klayerser);

[hsergio,ha,psergio,pa] = rtpread(fopX);
[hchris,ha,pchris,pa] = rtpread('/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/regr49_pbl.op.rtp'); %% by chris  LOOKS LIKE 370 ppm

compare_two_structures(psergio,pchris);

[sppmvLAY,sppmvAVG,sppmvMAX,spavgLAY,tavgLAY,sppmv500,sppmv75,sppmvSURF] = layers2ppmv(hsergio,psergio,1:49,2);
[cppmvLAY,cppmvAVG,cppmvMAX,cpavgLAY,tavgLAY,cppmv500,cppmv75,cppmvSURF] = layers2ppmv(hchris,pchris,1:49,2);
semilogy(nanmean(sppmvLAY,2),nanmean(spavgLAY,2),nanmean(cppmvLAY,2),nanmean(cpavgLAY,2)); set(gca,'ydir','reverse'); ylim([10 1100]); legend('sergio','chris');
