function [Z,P,JAC,p,fout,rout,dz,dp] = comparejacs_AIRSL2_vs_PBL(iAIRS_or_CrIS,i100_L2_or_PBL,iRTP);

%% see ~/git/MATLABCODE_Git/AveragingKernel_DOF/driver_vertical_resolution.m

if iAIRS_or_CrIS == 1
  str = 'AIRS';
elseif iAIRS_or_CrIS == -1
  str = 'CrIS';
elseif iAIRS_or_CrIS == 0
  str = 'IASI';
end

if i100_L2_or_PBL == 1
  str = [str ' L2-100'];
elseif i100_L2_or_PBL == -1
  str = [str ' PBL-100'];
end

fprintf(1,'doing %s \n',str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if i100_L2_or_PBL > 0
  if iAIRS_or_CrIS  > 0  
    sarta = ['/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 '];
  elseif iAIRS_or_CrIS  < 0 
    sarta = ['/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_jan25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_newFeb2025 '];
  elseif iAIRS_or_CrIS  == 0
    sarta = ['/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_sarta_iasi_jan25_H2020_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4 '];
  end
  klayers = ['/asl/packages/klayersV205/BinV201_chip/klayers_airs'];
elseif i100_L2_or_PBL < 0
  if iAIRS_or_CrIS  > 0
    sarta = ['/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_feb25_H2020_PBL '];
  elseif iAIRS_or_CrIS  < 0
    sarta = ['/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_feb25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new_PBL '];
  elseif iAIRS_or_CrIS  == 0
    sarta = ['/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_sarta_iasi_jan25_H2020_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4_PBL '];
  end
  klayers = [' /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_PBL/Feb05_2025/klayersV205/BinV201_chip/klayers_airs'];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/set_rtp.m
%cp /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/Aux_jacs_AIRS_August2018_2002_2017/LATS40_avg/Desc/latbin1_40.rp.rtp RTP/.
%use_this_rtp = 'RTP/latbin1_40.rp.rtp';

cd JUNK
cd /home/sergio/git/MATLABCODE_Git/AveragingKernel_DOF/JUNK/
if iRTP < 0
  error('need valid profile');
elseif iRTP > 0
  %% legal profile number, set correctly
  fip = '/asl/s1/sergio/alldata/RTP_pin_feb2002/regr49_maybe.ip.rtp';
  fop = '../LATBIN_1_40/junk.op.rtp';
  klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ughklayers'];

  fprintf(1,'runner_klayers = %s \n',klayerser)
  disp('running klayers')
  eval(klayerser);
  disp('ran klayers')

  runner_sarta = [' fin=' fop ' fout=blah.rad.rtp listp=' num2str(iRTP) ' listj=-1'];
end

runner_sarta = ['!time /bin/rm blah.rad.rtp*; time ' sarta ' ' runner_sarta ' >& ugh'];

fprintf(1,'runner_sarta = %s \n',runner_sarta)
disp('running sarta')
eval(runner_sarta);
disp('ran sarta')

%disp(' run this on chip!!!! run this on chip!!!! run this on chip!!!!')
%disp(' in /home/sergio/git/MATLABCODE_Git/AveragingKernel_DOF/JUNK/')

pause(0.1);
%disp('ret to continue'); pause
cd ../

%%%%%%%%%%%%%%%%%%%%%%%%% do this on chip %%%%%%%%%%%%%%%%%%%%%%%%%
if  iAIRS_or_CrIS > 0
  load /home/sergio/asl/s1/sergio/alldata/AIRSPRODUCTS_JACOBIANS/STD/temp_jac.mat
end

[h,ha,p,pa] = rtpread('JUNK/blah.rad.rtp');
p.plays = plevs2plays(p.plevs);
[w,jacWV, iaProf,iaNumLev] = readsarta_jac('JUNK/blah.rad.rtp_jacG1',1);     jacWV = squeeze(jacWV);
[w,jacO3, iaProf,iaNumLev] = readsarta_jac('JUNK/blah.rad.rtp_jacG3',3);     jacO3 = squeeze(jacO3);
[w,jacWGT,iaProf,iaNumLev] = readsarta_jac('JUNK/blah.rad.rtp_WGTFCN',200);  jacWGT = squeeze(jacWGT);
[w,jacTZ, iaProf,iaNumLev] = readsarta_jac('JUNK/blah.rad.rtp_jacTZ',100);   jacTZ = squeeze(jacTZ);
iaNumLay = iaNumLev - 1;

jacST = jacTZ(:,iaNumLev) * ones(1,4);

jacWV  = jacWV(:,1:iaNumLay);
jacO3  = jacO3(:,1:iaNumLay);
jacWGT = jacWGT(:,1:iaNumLay);
jacTZ  = jacTZ(:,1:iaNumLay);

jacWV  = fliplr(jacWV);
jacO3  = fliplr(jacO3);
jacWGT = fliplr(jacWGT);
jacTZ  = fliplr(jacTZ);

if iAIRS_or_CrIS > 0
  f2378 = instr_chans('airs',1);
  NeDT2378 = instr_chans('airs',2);
  NeDT = interp1(f2378,NeDT2378,h.vchan,[],'extrap');  
  g = dogoodchan;
elseif iAIRS_or_CrIS < 0
  f2211 = instr_chans('cris',1);
  NeDT2211 = instr_chans('cris',2);
  NeDT = interp1(f2211,NeDT2211,h.vchan,[],'extrap');
  g = 1 : 2223;
elseif iAIRS_or_CrIS == 0
  f8462 = instr_chans('iasi',1);
  %NeDT8462 = instr_chans('iasi',2);
  NeDT8462 = 0.5 * ones(length(f8462));
  NeDT = interp1(f8462,NeDT8462,h.vchan,[],'extrap');
  g = 1 : length(f8462);
  g = g(1:h.nchan);
end

JAC = [jacWV jacO3 jacTZ jacWGT jacST];  %% WV O3 T WGT SURF
foutkc = fout;
fout = w;
rout = p.rcalc; %% did jacobians for one profile so no need to subset

figure(7); clf; dz = abs(diff(p.palts))/1000;   semilogy(dz(1:p.nlevs-1),p.plays(1:p.nlevs-1)); set(gca,'ydir','reverse'); xlabel('dz'); ylabel('Plays (mb)'); ylim([10 1000])
figure(8); clf; Z = meanvaluebin(p.palts)/1000; plot(dz(1:p.nlevs-1),Z(1:p.nlevs-1)); xlabel('dz'); ylabel('Hlays (km)'); ylim([0 50])
P = p.plays(1:p.nlevs-1);
dp = abs(diff(p.plevs));
dp = dp(1:p.nlevs-1);

Z = flipud(Z);
P = flipud(P);
dz = flipud(dz);
dp = flipud(dp);
