function [freq,coefall,brmsall] = fit_extended_nonLTE_all_generic(opts)

%% cp -a   /home/chepplew/gitLib/ftc_dev/fit_nonLTE/src/fit_extended_nonLTE_all_generic_x.m .

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do a fit of deltaRnonLTE for a single channel
%
% INPUTS:
%  opts : struct with fields:
%    csens: {'CRIS_LR','CRIS_MR','CrIS_HR','AIRS_L1B','AIRS_L1C','IASI','IASI_NG','CHIRP'}
%    prod:  {'2018','2019','2020','2021'};
%    run_date: {'mar2018','sep2018','dec2018','mar2020','mar2021','apr2021'};
%    range:    {'day','night','c'} for 90 > solzen > 90 deg OR cosine(solzen).
%    regset    {'r49','saf704'};
%    vers      string. Additional subscript for output file name
% 
% first: merge_RnonLTE_generic.m (with extn = 1)
% then:  nonlte_change_with_co2_generic.m (with extn = 1)
% remember for CrIS to reorder_guard_chans.m
%
% Apr 2021 C. Hepplewhite: New version to stretch solar angles to 120-deg
%                          for nonLTE after sunset. ref <TBD>
%                          new variable xang{1,s,i} from sang{1,s,i}
% Jun 2022 C Hepplewhite:  version fits daylight and twilight separately
%                          using variable: range={'day','night'}
%                          first run: prod_2021/iasi/apr2021/
% Aug 2022 C Hepplewhite:  added 'x' in front of data file names for compliance
%                          with updated merge_RnonLTE and meanratio.... scripts.
% Mar 2023 CLH:  Changed input parameter list to structure and mods for cosine(solzen)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% cd /home/chepplew/gitLib/ftc_dev/run  % fit_nonLTE

addpath /asl/matlib/science        % vaconv.m saconv.m (was /asl/matlab2012/science)
addpath /asl/matlib/aslutil        % rad2bt.m
% addpath /home/chepplew/gitLib/ftc_dev/fit_nonLTE/src

disp('Initializing...')

% =================================
% Check input parameter and fields
% =================================
if( nargin ~= 1)
  error('One input plz: structure');
  return
end
if(~isstruct(opts))
  error('Input parameter needs to be a structure')
  return
end
%
if(~isfield(opts,'csens'))
  error('Sensor not defined')
  return;
else
  csens = upper(opts.csens);
end
%
if(~isfield(opts,'prod'))
  error('prod year not defined')
  return
else
  prod = opts.prod;
end
%
if(~isfield(opts,'run_date'))
  error('run_date not defined')
  return
else
  run_date = opts.run_date;
end
%
if(~isfield(opts,'regset'))
  error('regset not defined')
  return
else
  regset = upper(opts.regset);
end
%
if(~isfield(opts,'range'))
  error('range not set')
  return
else
  range = lower(opts.range);
end
%
if(~isfield(opts,'vers'))
  error('vers not set')
  return
else
  vers = lower(opts.vers);
end

% %band   = 'all';   % no longer used here.

% Which sensors
allsens = {'CRIS_LR','CRIS_MR','CrIS_HR','AIRS_L1B','AIRS_L1C','IASI','IASI_NG','CHIRP'};
if( ~ismember(csens, allsens) )
   error('Wrong instrument. Options are: cris_lr, cris_mr, cris_hr, airs_l1b, airs_l1c, iasi')
   return
end

% check choice of regression profile set
regset = upper(regset);
if( ~ismember(regset,{'R49','SAF704'}) )
  error('Incorrection option for regression profiles. Options are R49, SAF704');
  return
end

% Check kCARTA calcs date
allruns = {'mar2018','sep2018','dec2018','mar2020','mar2021','apr2021'};
if(~ismember(run_date,allruns)) error('Wrong run date');return; end

% Check Production run.
allprods = {'2018','2019','2020','2021','2022'};
if(~ismember(prod,allprods)); error('Wrong production'); return; end
prod_run = ['prod_' prod];

% solzen angle range
allranges={'day','night','c'};
if(~ismember(range,allranges)); error('Wrong range selection'); return;end

% display settings
fprintf(1, ...
  'Processing: %s, prod: %s, calcs: %s, range: %s, regset: %s, vers: %s\n', ...
   csens,prod,run_date,range,regset,vers);

% source and destination directories:
srcdr = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) '/' ...
         run_date '/nonLTE/'];
outdr = ['/home/chepplew/data/sarta/' prod_run '/' lower(csens) '/' ...
         run_date '/fitc/' lower(regset) '/'];
outdr = ['/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/NLTE/' prod_run '/' lower(csens) '/' ...
         run_date '/fitc/' lower(regset) '/'];
	 
if(strcmp(regset,'R49'))
  switch csens
    case 'CRIS_HR'
      sat_alt   = 834000;      % km
      datafile  = [srcdr 'xnlte_merged.mat'];
      ratiofile = [srcdr 'xnonlte_meanratiominusone_cris_hr.txt'];
      csname    = ['CRIS_HR_' regset ];
      fnout     = [outdr csname '_cutcoef_xnte_' vers];
    case 'AIRS_L1C'
      sat_alt   = 710000;      % km
      datafile  = [srcdr 'xnlte_merged.mat'];
      ratiofile = [srcdr 'xnonlte_meanratiominusone_airs_l1c.txt'];
      csname    = ['AIRS_L1C_' regset ];
      fnout     = [outdr csname '_cutcoef_xnte_' vers];      % was ... band]
    case 'IASI'
      sat_alt   = 817000;      % km
      datafile  = [srcdr 'xnlte_merged.mat'];
      ratiofile = [srcdr 'xnte_meanratiominusone_iasi.txt'];
      csname    = ['IASI_' regset ];
      fnout     = [outdr csname '_rawcoef_xnte_' vers];
      % /asl/s2/hannon/SpinachHome/Fit_deltaR_nonLTE/nonlte_meanratiominusone_iasi.txt
    case 'CHIRP'
      sat_alt   = 834000;      % km
      datafile  = [srcdr 'xnlte_merged.mat'];
      ratiofile = [srcdr 'xnonlte_meanratiominusone_chirp.txt'];
      csname    = ['CHIRP_' regset ];
      fnout     = [outdr csname '_rawcoef_xnte_' vers];
  end
  % prod_2018
  %datafile='/home/chepplew/data/sarta/prod_2018/p400/SO2/cris_hrg4_so2_data_long.mat';
  % prof 2021
  %datafle = '/home/chepplew/data/sarta/prod_2022/airs_l1c/apr2021/nonLTE/xnlte_merged.mat';
end

switch range
  case 'day'
    fnout = [fnout '_le90'];
  case 'night'
    fnout = [fnout '_gt90'];
  case 'c'
    fnout = [fnout '_cos'];
end

%
ratiodata = load(ratiofile);

%ratiodeltaco2 = 15.0; % 385-370 ppmv
ratiodeltaco2 = 1.0;  % 400 based on 400 ppmv

% Lowest layer to consider
lmax = 45;

% Number of nonLTE coeffcients
ncoef = 6;

% Min number of profiles
nprofmin = 6;

% Assign weight min & max
% Warning: do not weight too heavily; a small deltaRnonLTE
% radiance may still be a significant fraction of the total
% radiance if the atmos is cold.
wmax = 2;
wmin = 1;

% Load in the convolved data
%eval(['load ../Data/raddata_iasi_sep2008']);
%load('/asl/s1/chepplew/projects/sarta/cris_hr/nlte_raddata_cris_hrg4_jun2016');
%load('/asl/s1/chepplew/projects/sarta/cris_hr/nlte_raddata_cris_hrg4_400ppm');
%load('/home/chepplew/data/sarta/prod_2018/p400/nlte_rad_cris_hrg4_p400');
% prod_2018 H2016 nonLTE v120
load(datafile);

deltaR = radnonLTE - radLTE;
nchan  = length(freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ix = find(freq  >= 2330,1);  %% SERGIO DEBUG as in fit_dayNnight.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% transform solar zenith angle to new scale
%%%xangles = sangles*90/120;

% Find indices of profiles with solar angle <= 90 and > 90.
switch range
  case 'day'
    iwnt = find(sangles <  90.00);   % daylight sergio
    iwnt = find(sangles <= 90.00);   % daylight
  case 'night'
    iwnt = find(sangles > 90.00);  % extended range (twilight)
end

% Calc sun & view angles 
% Use top layer of atmosphere to model region of nonLTE
junk = zprof(1,:);
nobs = length(junk);
[vang1] = vaconv( vangles, sat_alt*ones(1,nobs), junk ); 
[sang1] = saconv( sangles, junk );

% Declare output variables
ncoefall = zeros(nchan,1);
coefall  = zeros(nchan,ncoef);
rmsall   = 1.0E+16*ones(nchan,1);
%
bcoefall = zeros(nchan,ncoef);
brmsall  = 1.0E+16*ones(nchan,1);

% Convert radiances to brightness temperature
btR  = real(rad2bt(freq,radLTE));
btnR = real(rad2bt(freq,radnonLTE));
dbt  = btnR - btR;

% Calculate upper atmosphere mean temperatures
%% tprof comes from /home/chepplew/gitLib/ftc_dev/fit_nonLTE/merge_RnonLTE_generic.m
%% tprof comes from /home/chepplew/gitLib/ftc_dev/fit_nonLTE/merge_RnonLTE_generic.m
%% and claims for iprof=[2:15 17:48]
[nlay,nprofang] = size(tprof);
layhb = 1:5;
thb   = mean(tprof(layhb,:));

disp('Looping over channels ...');
% Loop over the channels
for ichan = 1:nchan
   %clf
   % Pull out the data for this channel
   if strfind(opts.range,'day')
     ind  = find( abs(dbt(ichan,:)) > 0.02 ); 
   else
     ind  = find( abs(dbt(ichan,:)) > 0.02 ); 
     ind = 1 : length(dbt(ichan,:));
   end

if ichan == 184
  fprintf(1,'length(ind) for chan %4i freq %8.6f \n',[ichan freq(ichan)])
  whos dbt ind
end

   % remove solar angle 120 radiances
   %%%%%%%ind  = setdiff(ind,i120);
   % include samples with dbt>0.02 and sangle in range.
   ind  = intersect(ind,iwnt);
   npts = length(ind);
   %
   if (npts > nprofmin)
      dbti  = dbt(ichan,ind);
      deltaRi=deltaR(ichan,ind);
      Ri    = radLTE(ichan,ind);
      nRi   = radnonLTE(ichan,ind);
      freqi = freq(ichan)*ones(1,npts);
      sangi = sangles(ind);
      vangi = vangles(ind);
      scosi = cos(sangi.*pi/180);
      vseci = 1./cos(vangi.*pi/180);
      thbi  = thb(ind);
      scos1 = cos(sang1(ind).*pi/180);
      vsec1 = 1./cos(vang1(ind).*pi/180);
      %
      freq(ichan);
      x=1:npts;
      %
      btr=zeros(lmax,npts);
      bref=zeros(lmax);
      for i=1:lmax
         btr(i,:) = bt2rad(freqi',tprof(i,ind)')'; %'
         bref(i)  = bt2rad(freq(ichan),tref(i));
      end
      %
      %%%%%
      % Define weight
      rmax=max(abs(deltaRi));
      rmin=min(abs(deltaRi));
      m=(wmax - wmin)/(rmax - rmin);
      b=wmax - m*rmax;
      weight=m*deltaRi + b;
      %
%%%
%      % Loop over the layers
%      for laybtx=1:lmax
%%
%         % Calc some layer dependent predictors
%         jj=laybtx:(laybtx+5);
%         tx=mean( tprof(jj,ind) );
%%%
         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         term1=weight;
         term2=weight .* scos1;
         term3=weight .* (scos1.^2);
         term4=weight .* scos1 .* vsec1;
         term5=weight .* scos1 .* thbi;
         term6=weight .* scosi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         term1b=weight;
         term2b=weight .* scos1;
         term3b=weight .* (scos1.^2);
         term4b=weight .* scos1 .* vseci;
         term5b=weight .* scos1 .* thbi;
         term6b=weight .* scosi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %
         nterm=1;
         A=[term1'];  %'
         B=[term1b']; %'
         %
         if (npts > 12)
           A=[A, term2'];  %'
           B=[B, term2b']; %'
           nterm=2;
           %
           if (npts > 25)
             A=[A, term3'];  %'
             B=[B, term3b']; %'
             nterm=3;
             %
             if (npts > 50)
               A=[A, term4'];  %'
               B=[B, term4b']; %'
               nterm=4;
               %
               if (npts > 100)
                 A=[A, term5'];  %'
                 B=[B, term5b']; %'
                 nterm=5;
                 %
                 if (npts > 200)
                   A=[A, term6'];  %'
                   B=[B, term6b']; %'
                   nterm=6;
                 end
               end
             end
           end
         end
         %
         f=weight.*deltaRi;
         coef=A\f'; %'
         calc=A*coef;
         %
         bcoef=B\f'; %'
         bcalc=B*bcoef;
         %
         % Calculated total nonLTE radiance
         rjunk  = Ri + calc'./weight;   %'  calc nRi
         rjunkb = Ri + bcalc'./weight;  %' bcalc nRi
         %
         % Calc dBT for nonLTE
         bjunk  = real(rad2bt(freqi,rjunk))  - btR(ichan,ind); %'
         bjunkb = real(rad2bt(freqi,rjunkb)) - btR(ichan,ind); %'
         bjunkt = dbt(ichan,ind)'; %'
         %
         wmean=mean(weight);
         dif=weight'.*(calc./f' - 1)./wmean;
         rms=sqrt(mean(dif.^2));
         dif=weight'.*(bcalc./f' - 1)./wmean;
         brms=sqrt(mean(dif.^2));
%{
         subplot(311)
           plot(x,100*(calc-f')./f','ro'),grid
           title(['deltaRnonLTE % error chan ' int2str(ichan)]);
         subplot(312)
           plot(x,bjunk,'ro',x,bjunkt,'k.'),grid
         subplot(313)
         % Show as error in deltaBTnonLTE
           plot(x,bjunk-bjunkt','ro', x,bjunkb-bjunkt','cx'),grid
         %
         pause
%}
         ncoefall(ichan) = nterm;
         coefall(ichan,1:nterm) = coef;
         rmsall(ichan) = rms;
         %
         bcoefall(ichan,1:nterm) = bcoef;
         brmsall(ichan) = brms;

%%%
%      end % End loop over layers
   end  % end if npts > nprofmin

   if ichan == ix
     coef
     switch range
       case 'day'
         saver = ['save debug_fit_extended_nonLTE_all_generic_x_ichan_' num2str(ichan) '_filtered_with_tprof.mat A f coef weight term1 term2 term3 term4 term5 term6 tprof vangles sangles vang1 sang1 zprof thbi thb'];
       case 'night'
         saver = ['save debug_night_fit_extended_nonLTE_all_generic_x_ichan_' num2str(ichan) '_filtered_with_tprof.mat A f coef weight term1 term2 term3 term4 term5 term6 tprof vangles sangles vang1 sang1 zprof thbi thb'];
     end
     whos A f coef weight term1 term2 term3 term4 term5 term6 tprof vangles sangles vang1 sang1 zprof thbi thb
     fprintf(1,'saver = %s \n',saver);
     eval(saver);
   end

end % End of loop over channels
whos *all dbt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine which channels to write out
amaxc    = max( abs(coefall),[],2 );
dbtmax   = max(dbt, [], 2);
dbtmean  = mean(dbt, 2);
iwant    = find( amaxc > 1E-7 & dbtmax > 0.03 & dbtmean > 0.003 & ncoefall >= 2);

%{
clf;
plot(freq,dbtmax,'b',freq,dbtmean,'c',freq,amaxc,'r',freq,ncoefall,'m');grid on;
hold on;axis([2200 2450 -0.5 35]);legend('dbtmax','dbtmean','amaxc','ncoef')
plot(freq(iwant),dbtmax(iwant),'bo')
disp('review selection')
pause
%}

disp(['iwant : ' num2str(length(iwant)) ' idchan : ' num2str(length(idchan))])
idchan  = idchan(iwant);
freq    = freq(iwant);
coefall = coefall(iwant,:);
dbtmean = dbtmean(iwant);
brmsall = brmsall(iwant);

% Create list file
junk = [idchan(:), freq(:), dbtmean(:)];
%junk = [idchan(iwant)', freq(iwant), dbtmean(iwant)];

if(exist([fnout '_list']))
  disp('list file already exists, i wont overwrite it')
else
  disp(['Writing list file: ' fnout '_list']);
  %%%%save([fnout '_list'],'junk', '-ascii')
end

% Write out coefficients
ierr  = 0;
fname = [fnout '_6term.dat'];
if(exist(fname))
  disp('6-term file already exists, i wont overwrite it')
else
  disp(['Saving 6-term nonLTE coefficients to: ' fname]); 
  %%%ierr = wrtcoef_nte(fname, idchan(iwant), freq(iwant), coefall(iwant,:));
  ierr = wrtcoef_nte(fname, idchan, freq, coefall);
  if (ierr ~= 0)
     disp('ERROR writing 6 term coef file')
  end
end
%
fname = [fnout '_7term.dat'];

% Option to write 7th-term as zeros: set LTERM7=0.
LTERM7 = logical(1);

% for AIRS need to sort frequency before comparing
[~, isf] = sort(freq);
[~, isr] = sort(ratiodata(:,2));

if(exist(fname))
  disp('7-term file already exists, i wont overwrite it')
else

  if(LTERM7)
    addpath /asl/packages/airs_decon/source
    [ix iy] = seq_match(ratiodata(isr,2), freq(isf));
    %%%%[fx iy] = intersect(ratiodata(:,2), freq);
    %whos freq ratiodata fx ij ix iy
    junk  = max( abs(ratiodata(isr(ix),2) - freq(isf(iy))) );
    if (junk > 0.01)
      disp('ERROR ratio data freq does not match coef data freq')
    else
      %%%coefall7 = [coefall ratiodata(idchan,3)/ratiodeltaco2];
      coefall7 = [coefall ratiodata(isr(ix),3)/ratiodeltaco2];
      disp(['Saving 7-term nonLTE coefficients to: ' fname]);
      ierr = wrtcoef_nte(fname, idchan, freq, coefall7);
      if (ierr ~= 0)
         disp('ERROR writing 7 term coef file')
      end
    end

  else
    ns    = size(coefall,1);
    vect0 = zeros(ns,1);
    coefall7 = [coefall vect0];
    disp(['Saving 7-term nonLTE coefficients to: ' fname]);
    ierr = wrtcoef_nte(fname, idchan, freq, coefall7);
    if (ierr ~= 0)
       disp('ERROR writing 7 term coef file')
    end
  end
end


%%% end of program %%%
