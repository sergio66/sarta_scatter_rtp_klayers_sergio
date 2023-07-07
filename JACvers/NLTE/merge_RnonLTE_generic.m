function [iok] = merge_RnonLTE_generic(csens, prod, run_date, regset, extn)

% copied from /home/chepplew/gitLib/ftc_dev/fit_nonLTE/merge_RnonLTE_generic.m


%
% INPUTS:
%   csens:    the sensor to work with: {'airs_l1b','airs_l1c','cris_hr','iasi'}
%   prod:     production year e.g. '2018' etc
%   run_date: date of ODs or sim. e.g. 'dec2018' etc
%   regset:   the regression set to work with: {'r49','saf704'}
%   extn:     Extended (0-120.deg solzen) 1=yes, 0=No. (0-90.deg solzen) 
%
% This is a matlab program to merge the individual profile LTE and nonLTE
% radiance data files into one matlab output file containing *all* the
% data needed to do a fast model for deltaRnonLTE = R_nonLTE - R_LTE.
%
% Besides the matlab files with the rad and angle data, we might
% also want to use the profile temperatures as a predictor, as well
% as reference profile.  This is pulled in from separate a RTP file.
%
% NB the VT*.mat convolved files give 352 channels for AIRS, but for CrIS NSR,FSR
%    MedR, IASI, the full channel grid is supplied. So trim channels except for AIRS.
%
% Created: Scott Hannon, 9 December 2003 (based on "mergetherm.m")
% Update: 25 February 2005, S.Hannon - minor changes for Jan/Feb 2005 data
% Update: 14 March 2005, S.Hannon - minor changes for March 2005 X data
% Update: 06 Sep 2005, S.Hannon - minor changes for Sept 2005 data (no "X")
% Update: 10 Oct 2005, S.Hannon - minor changes for Oct 2005 data
% Update: 15 Mar 2007, S.Hannon - minor changes for IASI Mar 2007 data
% Update: 03 Oct 2008, S.Hannon - minor changes for IASI Sep 2008 data
% Update: 02 Aug 2016, C Hepplewhite - Modified to work with Sergio's kcarta runs
% 20 Mar 2018. CLH: updated for prod_2018 and HITRAN2016 kcarta.
% 1 May 2018.  CLH. Tweaks for v120 nonLTE kcarta. 
% Dec 2018.    CLH. adapt for use with specified sesnor, production run and kCARTA.
%                   prod_run, kcarta_run, all_sensors, and path defs.
% Mar 2021     CLH  updated paths to source dpath.
% Apr 2021     CLH: added IASI-NG. Extended solar zenith angles to past-eclipse.
%                   19 solar x 6 view angles. Updated paths *again*.
% May 2021     CLH  Added CRIS_LR
% Jul 2022     CLH  added option 'extn' whether to use std fitting (0-90 deg solzen)
%                   or extended (0-120 deg solzen).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/chepplew/gitLib/ftc_dev/run

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the variables csens, regset, run_date, prod_run as needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( nargin ~= 5) 
  disp('Wrong number inputs: using defaults'); 
  csens    = 'CRIS_HR'; % 'CRIS_HR'; % 'AIRS_L1C'; % 'CRIS_LR';
  regset   = 'r49';
  prod     = '2022';
  run_date = 'jul2022';
  extn     = 0;
end

% Check kCARTA calcs date
allruns = {'mar2018','sep2018','dec2018','mar2020','mar2021','apr2021', ...
           'may2021','jul2022'};
if(~ismember(run_date,allruns)) error('Wrong run date');return; end

% Check Production run.
allprods = {'2018','2019','2020','2021','2022'};
if(~ismember(prod,allprods)); error('Wrong production'); return; end
prod_run = ['prod_' prod];

% check sensor
csens = upper(csens);
allsens = {'CRIS_LR','CRIS_HR','AIRS_L1B','AIRS_L1C','IASI','IASI_NG','CHIRP'};
if( ~ismember(csens,allsens) )
   error('Wrong instrument.')
   return
end

% check choice of regression profile set
regset = upper(regset);
if( ~ismember(regset,{'R49','SAF704'}) )
  error('Incorrection option for regression profiles. Options are R49, SAF704');
  return
end

% Check extn
if(0 > extn | extn > 1)
  error('invalid value for extn [0,1]'); return;
end
if(extn == 0)
  fno_sffx = '';
elseif(extn == 1)
  fno_sffx = 'x';
end

% display configuration
fprintf(1,'Sensor: %s. prod run: %s. build: %s. extended? %i\n', csens,prod_run,run_date,extn)

nlay=97;
psurf=1013.25;

% comment: string to be included in output file
%%%comment='02-Aug-2016 CrIS HiRes 4g coefficients for deltaRnonLTE CO2 400ppm';

% rtpfile : Name of the matlab file with the reg & ref profiles
% The file must be a "layers" contain the variables:
%    head.ptype = 1 ("layers" profile)
%    prof.nlevs (1 x nprof+1)
%    prof.spres (1 x nprof+1)
%    prof.ptemp (nlay x nprof+1)
% where the reference profile is number nprof+1

% set the source directories and files:
if(strcmp(regset,'R49'))
  nprof = 48;
  switch run_date
    case 'mar2018'
      %dpath   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/';
      dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2016_Mar2018/'];
      rtpfile=[dpath 'regr49_1100_400ppm_unitemiss.op.rtp'];

    case 'sep2018'
       dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
             'REGR49_400ppm_H2016_Sept2018_AIRS2645/'];
       rtpfile=[dpath 'regr49_1100_400ppm_unitemiss.op.rtp'];

    case 'dec2018'
      dpath=['/home/sergio/KCARTA/NONLTE_PRODUCTION/VT_48PROFILES_120_400ppmv_v120_' ...
             'H16_NLTEH16_Dec2018/Results/CONV_Results/'];
      rtpfile=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/'...
               'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1100_400ppm_unitemiss.op.rtp'];
      % updated for AIRS L1C
      dpath=['/home/sergio/KCARTA/NONLTE/VT_48PROFILES_120_400ppmv_v121_H16_NLTEH16_Apr2019/'...
             'Results/CONV_Results_Apr23_2019/'];

    case 'mar2020'
       dpath = ['/asl/s1/sergio/NLTE_CALCS/VT_48PROFILES_120_SRCv1.21_400ppmv_H16_Mar2020/' ...
               'CONV_Results/'];
       rtpfile=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/'...
            'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1100_400ppm_unitemiss.op.rtp'];

    case 'mar2021'
       dpath = ['/home/sergio/KCARTA/NONLTE_PRODUCTION/' ...
                'VT_48PROFILES_120_400ppmv_v121_H16_NLTEH16_Apr2021_NewNLTEProfiles_SAVETHISGOOD/' ...
		'Results/CONV_Results/'];
       rtpfile=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/'...
            'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1100_400ppm_unitemiss.op.rtp'];

    case 'apr2021'
      % extended modeling with 19 solzen to 120.deg
      %dpath = ['/asl/s1/sergio/NLTE_CALCS/' ...
      dpath = ['/home/chepplew/data/KCARTA/NLTE_CALCS/' ...
               'VT_48PROFILES_120_SRCv1.21_400ppmv_H16_Apr2021_NewNLTEProfiles/CONV_Results/'];
      rtpfile=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/'...
            'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1100_400ppm_unitemiss.op.rtp'];
     
    case 'may2021'
      dpath = ['/home/sergio/KCARTA/NONLTE_PRODUCTION/' ...
               'VT_48PROFILES_120_400ppmv_v121_H16_NLTEH16_Jan2020/Results/CONV_Results/']; 
      rtpfile = ['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/'...
            'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1100_400ppm_unitemiss.op.rtp'];

    case 'jul2022'
      % USES: apr2021 extended modeling with 19 solzen to 120.deg
      dpath = ['/asl/s1/sergio/NLTE_CALCS/' ...
               'VT_48PROFILES_120_SRCv1.21_400ppmv_H16_Apr2021_NewNLTEProfiles/CONV_Results/'];
      rtpfile=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/'...
            'REGR49_400ppm_H2016_Dec2018_AIRS2834/regr49_1100_400ppm_unitemiss.op.rtp'];

  end
  pnums   = [1:48]';
  comment = [csens ' r49 400ppm CO2 H2016 ' fno_sffx '_' run_date];
  outdr   = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) '/' run_date '/'];
end
if(strcmp(regset,'SAF704'))
  %dpath   = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/SAF704/';
  dpath=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/' ...
         'SAF704_400ppm_H2016_Mar2018/'];
  rtpfile = [dpath 'save_SAF_704_profiles_29-Apr-2016_1100mb.op.rtp'];
  pnums = [1:703]';
  comment = [csens ' SAF704 400ppm CO2 HIT2016'];
end

%rtpfile='pin_feb2002_reg_op.rtp';
%rtpfile='/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49/regr49_1013_385ppm.op.rtp';
%rtpfile=['/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/'...
%         'regr49_1100_400ppm.op.rtp'];
%rtpfile=['/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/'...
%         'REGR49_400ppm_H2016_Mar2018/regr49_1100_400ppm.op.rtp'];

% raddir: dir with radiance data files
%raddir='/carrot/s1/sergio/VT_48PROFILES_120_SRCv1.12_370ppmv/CONV_Results';
%raddir='/home/sergio/KCARTA/NONLTE/VT_48PROFILES_120_385ppmv_v118/Results/CONV_Results/';
%raddir='/home/sergio/KCARTA/NONLTE/VT_48PROFILES_120_400ppmv_v118/Results/CONV_Results/';
%raddir='/home/sergio/KCARTA/NONLTE/VT_48PROFILES_120_400ppmv_v118/Results/CONV_Results/';
% Ths update 1-May-2018 Sergio update to nonLTE alogirthm
%raddir=['/home/sergio/KCARTA/NONLTE/VT_48PROFILES_120_400ppmv_v120_H16_NLTEH16/'...
%         'Results/CONV_Results/'];

% outfile: Name prefix of output matlab file to create
%outfile='/home/chepplew/data/sarta/prod_2018/p400/nlte_v120_rad_cris_hrg4_p400';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The main executable code begins below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the profile temperature data
[head, hattr, prof, pattr] = rtpread(rtpfile);
%
% Check that the RTP data looks plausible
if (head.ptype ~= 1)
   error('RTP file must contain "layers" type profiles')
end
ii=nprof + 1;
if (length(prof.nlevs) < ii)
   error('RTP file contains too few profiles')
end
if (length(prof.nlevs) > ii)
   disp('WARNING: RTP file contains more profiles than expected')
end
if ( min(prof.nlevs(1:ii))-1 < nlay)
   error('Some RTP profiles have fewer layers than expected')
end
if ( min(prof.spres(1:ii)) < psurf)
   error('Some RTP profiles have lower surface pressure than expected')
end
% Note: prof.spres < psurf is bad, but prof.spres > psurf is OK
% Pull out the temperature data
temp = prof.ptemp(1:nlay,1:nprof);
tref = prof.ptemp(1:nlay,nprof+1);

%plays=prof.plays(1:nlay,1);       %%%% prof.plays does not exist anymore! so do the following
pN    = prof.plevs(1:nlay,:)-prof.plevs(2:nlay+1,:);
pD    = log(prof.plevs(1:nlay,:) ./ prof.plevs(2:nlay+1,:));
plays = zeros(nlay,size(prof.plevs,2));
plays(1:nlay,:) = pN ./ pD;

plevs = prof.plevs(1:(nlay+1),1);
palts = 0.5*(prof.palts(1:nlay,1:nprof) + prof.palts(2:(nlay+1),1:nprof));
whos temp tref plays plevs palts;
clear head prof hattr pattr

% Check all profiles are available
junk = dir([dpath 'vt*.mat']);
if(length(junk) ~= nprof)
  error('The number of vailable profiles does not match expectations');
  return;
end

% Read in the first nLTE data file
% Note: "xrad" files used xfconvkc convolver
iprof = 1;
xx    = load([dpath,'vt' int2str(iprof) '.mat']);

% Sergio has computed all channels for CrIS, IASI but only some SW channels for AIRS
switch csens
  case 'AIRS_L1C'
    hinfo   = hdfinfo('/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables_lls_new.hdf');
    fchn    = hdfread(hinfo.SDS(2));
    idchn   = hdfread(hinfo.SDS(1));
    %yy = load('~chepplew/projects/sarta/prod_2019/airs_l1c/airs_l1c_chans_for_nte.mat');
    %freq   = yy.RnonLTE.vchan;
    %idchan = yy.RnonLTE.ichan;
    %ichan  = [1:length(idchan)];
    freq   = xx.fairs;
    [~, idchan] = intersect(fchn, freq);
    ichans = [1:length(freq)];
    
  case 'CRIS_LR'
    %band3  = textread('/home/chepplew/gitLib/ftc_dev/chanLists/cris_lr_g4_list_all');
    %freq   = intersect(xx.fcris,band3(:,2))'; %'
    %ichan  = band3(:,1)';
    freq   = xx.lo_fcris;
    idchan = [1:length(freq)];    % NB. does not record guard channels

  case 'CRIS_HR'
    %band3  = textread('/home/chepplew/gitLib/ftc_dev/chanLists/cris_hr_g4_B3_list');
    %freq   = intersect(xx.fcris,band3(:,2))'; %'
    %ichan  = band3(:,1)';
    freq   = xx.hi_fcris;
    idchan = [1:length(freq)];    % NB. does not record guard channels
  case 'IASI'
    freq   = xx.fiasi;
    idchan = [1:8461];
  case 'IASI_NG'
    freq   = xx.ng_fiasi;
    idchan = [1:length(freq)]';
  case 'CHIRP'
    freq   = xx.med_fcris;
    idchan = [1:length(freq)];
end
% ---------------------------------------------------------------------------
% fitting to daytime only using extended modelling requires restricting range
%   to solzen <=90 otherwise <=120 for extended.
% ---------------------------------------------------------------------------
usolz  = unique(xx.solar);
switch extn
  case 0
    iwsz   = find(xx.solar <= 90);
  case 1
    iwsz   = find(xx.solar <= 120);
end
nchan  = length(freq); 
nang   = length(xx.solar(iwsz));           % number of angle combos exp: 36 or 72
nsang  = length(unique(xx.solar(iwsz)));   % number of solar angles exp: 6 or 12
nvang  = length(unique(xx.viewer));        % number of view  angles exp: 6.

% Declare output arrays
tprof     = zeros(nlay, nang*nprof);
zprof     = zeros(nlay, nang*nprof);
radLTE    = zeros(nchan, nang*nprof);
radnonLTE = zeros(nchan, nang*nprof);
sangles   = zeros(1,nang*nprof);
vangles   = zeros(1,nang*nprof);
whos tprof zprof radLTE radnonLTE sangles vangles;

% Start to fill arrays for iprof=1
% NB some (cris -> hi_ med_) field names changed after adding CHIRP.
ind = 1:nang;

switch csens
  case 'AIRS_L1C'
    radLTE(:,ind)    = xx.lte_airs_Bfunny4um(:,iwsz);
    radnonLTE(:,ind) = xx.nlte_airsB(:,iwsz);
  case 'CRIS_LR'
    radLTE(:,ind)    = xx.lo_lte_cris_Bfunny4um(idchan,iwsz);
    radnonLTE(:,ind) = xx.lo_nlte_crisB(idchan,iwsz);
  case 'CRIS_HR'
    radLTE(:,ind)    = xx.hi_lte_cris_Bfunny4um(idchan,iwsz);
    radnonLTE(:,ind) = xx.hi_nlte_crisB(idchan,iwsz);
    %radLTE(:,ind)    = xx.hi_lte_cris_Bfunny4um(idchan,iwsz);
    %radnonLTE(:,ind) = xx.hi_nlte_crisB(idchan,iwsz);
  case 'IASI'
    radLTE(:,ind)    = xx.lte_iasi_Bfunny4um(idchan,iwsz);
    radnonLTE(:,ind) = xx.nlte_iasiB(idchan,iwsz);
  case 'IASI_NG'
    radLTE(:,ind)    = xx.ng_lte_iasi_Bfunny4um(idchan,iwsz);
    radnonLTE(:,ind) = xx.ng_nlte_iasiB(idchan,iwsz);
  case 'CHIRP'
    radLTE(:,ind)    = xx.med_lte_cris_Bfunny4um(idchan,iwsz);
    radnonLTE(:,ind) = xx.med_nlte_crisB(idchan,iwsz);
end
%%%
tprof(:,ind) = temp(:,iprof)*ones(1,nang);
zprof(:,ind) = palts(:,iprof)*ones(1,nang);
sangles(ind) = xx.solar(iwsz);
vangles(ind) = xx.viewer(iwsz);
clear xx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over the other profiles
xprf = 0;
%%for iprof=2:nprof
%%for iprof=[2:16 20:48]  % #17,18,19 corrupt kC data
for iprof=[2:15 17:48]   % omit #16 too cold
  xprf = xprf+1;
  ind = ind + nang;
  xx  = load([dpath,'vt' int2str(iprof) '.mat']);
  tprof(:,ind) = temp(:,xprf)*ones(1,nang);
  zprof(:,ind) = palts(:,xprf)*ones(1,nang);
%
  switch csens 
    case 'AIRS_L1C'
      radLTE(:,ind)    = xx.lte_airs_Bfunny4um(:,iwsz);
      radnonLTE(:,ind) = xx.nlte_airsB(:,iwsz);
    case 'CRIS_LR'
      radLTE(:,ind)    = xx.lo_lte_cris_Bfunny4um(idchan,iwsz);
      radnonLTE(:,ind) = xx.lo_nlte_crisB(idchan,iwsz);
    case 'CRIS_HR'
      radLTE(:,ind)    = xx.hi_lte_cris_Bfunny4um(idchan,iwsz);
      radnonLTE(:,ind) = xx.hi_nlte_crisB(idchan,iwsz);
      %radLTE(:,ind)    = xx.hi_lte_cris_Bfunny4um(idchan,iwsz);
      %radnonLTE(:,ind) = xx.hi_nlte_crisB(idchan,iwsz);
    case 'IASI'
      radLTE(:,ind)    = xx.lte_iasi_Bfunny4um(idchan,iwsz);
      radnonLTE(:,ind) = xx.nlte_iasiB(idchan,iwsz);
    case 'IASI_NG'
      radLTE(:,ind)    = xx.ng_lte_iasi_Bfunny4um(idchan,iwsz);
      radnonLTE(:,ind) = xx.ng_nlte_iasiB(idchan,iwsz);
    case 'CHIRP'
      radLTE(:,ind)    = xx.med_lte_cris_Bfunny4um(idchan,iwsz);
      radnonLTE(:,ind) = xx.med_nlte_crisB(idchan,iwsz);    
  end
%
   sangles(ind) = xx.solar(iwsz);
   vangles(ind) = xx.viewer(iwsz);
   fprintf(1,'.');
end
fprintf('\n');

% Clean up unwanted variables
clear xx iprof temp

% with corrupt VT files, arrays are zero-filled toward the end
%  use max value of index counter ind to trim.
ind_end = max(ind);
radLTE(:,ind_end+1:end)    = [];
radnonLTE(:,ind_end+1:end) = [];
sangles(ind_end+1:end)     = [];
vangles(ind_end+1:end)     = [];
tprof(:,ind_end+1:end)     = [];
zprof(:,ind_end+1:end)     = [];
nprof                      = xprf+1;
palts(:,nprof+1:end)       = [];
plays(:,nprof+1:end)       = [];

% equalize values for the unwanted channels
inuke = find(freq < 2205.999 | freq > 2393.001);
radnonLTE(inuke,:) = radLTE(inuke,:);

% Save output
if(~exist([outdr '/nonLTE'])) 
  disp(['Creating output directory: ' outdr '/nonLTE/']);
  mkdir([outdr '/nonLTE/']);
end
outfn   = [ fno_sffx 'nlte_merged.mat'];
savVars = {'freq', 'idchan', 'tprof', 'tref', 'sangles', 'vangles', 'radLTE',...
           'radnonLTE', 'nprof' 'nlay' 'nang' 'nsang' 'nvang' 'comment' 'plays',...
           'plevs', 'palts', 'zprof'};
disp(['Saving data to file: ' [outdr '/nonLTE/' outfn]]);

save([outdr '/nonLTE/' outfn], savVars{:});

%clear
disp('Merge_RnonLTE_generic : Completed')

iok = nchan;

%%% end of program %%%
