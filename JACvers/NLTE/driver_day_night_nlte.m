addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER

%{
see ../src/calnte.f
       THIGH = (TEMP(1) + TEMP(2) + TEMP(3) + TEMP(4) + TEMP(5))/5.0
       PRED1 = 1.0
       PRED2 = SCOS1
       PRED3 = SCOS1*SCOS1
       PRED4 = SCOS1*VSEC1
       PRED5 = SCOS1*THIGH
       PRED6 = SUNCOS

C    REAL      SUNCOS  solzen cosine at surface    none
C    REAL      SCOS1   solzen cosine at layer1     none
C    REAL      VSEC1   satzen secant at layer1     none

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if I load in rtp file, subset first 48 profiles and do 
%%    T_high = sum(p.ptemp(01:05,:))/5
%% then profile 16 has a very low T_high!!!  (180 K vs 220 K)
%% Chris H throws that out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts.csens    = 'AIRS_L1C';
opts.prod     = '2022';
opts.run_date = 'apr2021';
opts.range    = 'day';
opts.range    = 'night';
opts.regset   = 'r49';
opts.vers     ='hmmmm_';

%% /home/chepplew/data/sarta/prod_2022/airs_l1c/apr2021/nonLTE/xnonlte_meanratiominusone_airs_l1c.txt exists
iDoSARTA = -1;
if iDoSARTA > 0
  [freq,coefall,brmsall] = fit_extended_nonLTE_all_generic_x(opts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('p')
  rtpfile = '/asl/s1/sergio/RTP_pin_feb2002/pin_feb2002_sea_airsnadir_op.rtp';
  %% see /home/chepplew/gitLib/ftc_dev/fit_nonLTE/merge_RnonLTE_generic.m
  rtpfile=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/'...
               'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1100_400ppm_unitemiss.op.rtp'];

  [h,ha,p,pa] = rtpread(rtpfile);
  profind = 1:48;    %% fpr some reason, never did profile 49
  %profind = setdiff(profind,16);  %% Chris H throws out this "cold at TOA" profile
  [h,p] = subset_rtp(h,p,[],[],profind);
  fprintf(1,'length(p.stemp) == %2i  length(profind) == %2i \n',length(p.stemp),length(profind))
  T_high0 = sum(p.ptemp(1:5,:))/5; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% let solzen  --> sangles = solar zenith = solar angle at GND
%% let scanang --> vangles = scanang at satellite
%%
%% let satzen --> vang1 = vang at TOA
%% let solang --> sang1 = satellite angle at TOA

%% April 27, 2021, private slack channel with me and Chris
%% nlte_B and lte_Bfunny4um
iVersKC = 0; %% SARTA NLTE 2019, n all solang <= 90
iVersKC = 1; %% SARTA NLTE 2021, n all solang <= 120
if ~exist('nlte')
  read_kcarta_nlte_data   %% does airs, though can easily switch to others, uses palt(1) here also
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define day vs night
day = find(sangles < 90);
day = find(sangles <= 90);
night = find(sangles > 90);

ind_day_for_sang2_vang2 = [];
ind_night_for_sang2_vang2 = [];
for ii = 1 : oo
  ind_day_for_sang2_vang2   = [ind_day_for_sang2_vang2   ((ii-1)*length(vangles) + day)];
  ind_night_for_sang2_vang2 = [ind_night_for_sang2_vang2 ((ii-1)*length(vangles) + night)];
end

iMax = 48; %% see all profiles
iMax = 06; %% boring, only see some
[Y,I] = sort(sangles);
for ii = 1 : iMax
  figure(1); clf
  plot(fairs,nanmean(squeeze(tnlte(:,:,ii)),2),'r',fairs,nanmean(squeeze(tlte(:,:,ii)),2),'b')
  plot(fairs,nanmean(squeeze(tnlte(:,:,ii)),2)-nanmean(squeeze(tlte(:,:,ii)),2),'r')
  plot(fairs,nanmean(squeeze(tnlte(:,:,ii)),2)-nanmean(squeeze(tlte(:,:,ii)),2),'rx-',...
       fairs,nanmean(squeeze(tnlte(:,day,ii)),2)-nanmean(squeeze(tlte(:,day,ii)),2),'b',...
       fairs,nanmean(squeeze(tnlte(:,night,ii)),2)-nanmean(squeeze(tlte(:,night,ii)),2),'k')
  plotaxis2;
  title(num2str(ii))
  hl = legend('all','day','night','location','best'); 

  figure(2); clf
  ix = find(fairs >= 2361,1);
  blah = squeeze(tnlte(ix,:,ii)-tlte(ix,:,ii));
  plot(sangles(I),blah(I)); xlabel('Solzen angles'); ylabel('NLTE-LTE'); 
  title(num2str(ii)); plotaxis2;
  pause(0.1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fit_dayNnight
