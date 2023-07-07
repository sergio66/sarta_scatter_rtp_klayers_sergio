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

%% deltarad(Nx1) = Pred(Nx6) Coeff(6x1)        
%%   where N = number sangles x num vangles x number profiles = 12 x 6 x 48 = 3456
%% 12 = length(unique(sangles(day))), 6 = length(unique(sangles(vangles)))), 48 = numprofs
%% dr = M c ==> M' dr = M'M c ==> c = inv(M'M) M'dr

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
duO3  = dobson_rtp(h,p);
duO3  = dobson_gas_rtp(h,p,3);
duN2O = dobson_gas_rtp(h,p,4);
duCO  = dobson_gas_rtp(h,p,5);

% l2s = load('/home/sergio/KCARTA/L2SComparisons/l2s_kc122_H16_605_2830.mat');
% [fc,qc] = quickconvolve(l2s.w,l2s.d(:,[3 4 5]),1,1);   
% plot(fc,exp(-qc(:,1)),'b',fc,exp(-qc(:,2)),'r',fc,exp(-qc(:,3)),'k','linewidth',2); xlim([2200 2300])
% hl = legend('3','4','5','location','best');  %% so mostly N2O, very little CO, no )3 in this region! amazing!

clear coefN coefD

iVersFit = 0; %% Scott/Sergio/Strow/Manuel 2005 GRL
iVersFit = 1; %% Marco/Manuel 2019
iVersFit = 2; %% Sergio 2023

iVersFit

if iVersFit == 0
  sarta = read_sarta_nlte_coef();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DAY
%T_highD = ones(length(day),1) * sum(p.ptemp(01:05,:))/5; T_highD = T_highD'; T_highD = T_highD(:);
T_highD = ones(length(day),1) * sum(p.ptemp(01:05,:))/5; T_highD = T_highD(:);
T_medD = ones(length(day),1) * sum(p.ptemp(06:15,:))/10; T_medD  = T_medD(:);
T_lowD = ones(length(day),1) * sum(p.ptemp(16:25,:))/10; T_lowD  = T_lowD(:);
o3D    = ones(length(day),1) * duO3;                       o3D     = o3D(:);
if strfind(opts.range,'day')
  disp('using Chris ThighD'); T_highD = chris2D.thbi';
end
  disp('using Chris ThighD'); T_highD = chris2D.thbi';
suncosD = cos(sangles(day)*pi/180)' * ones(1,length(profind));        suncosD = suncosD(:);
scos1D  = cos(sang1(day)*pi/180)' * ones(1,length(profind));          scos1D  = scos1D(:);
vsec1D  = sec(vang1(day)*pi/180)' * ones(1,length(profind));          vsec1D  = vsec1D(:);
scos2D  = cos(sang2(ind_day_for_sang2_vang2)*pi/180)';
vsec2D  = sec(vang2(ind_day_for_sang2_vang2)*pi/180)';
scos1D = scos2D;
vsec1D = vsec2D;

unityD = ones(size(T_highD));

if iVersFit == 0
  %% Scott/Sergio/Strow/Manuel 2005 GRL
  pred1D = unityD;
  pred2D = scos1D; 
  pred3D = pred2D.*pred2D; 
  pred4D = pred2D.*vsec1D; 
  pred5D = pred2D.*T_highD; 
  pred6D = suncosD;
  predsD = [pred1D pred2D pred3D pred4D pred5D pred6D];
elseif iVersFit == 1
  %% Marco/Manuel 2019
  pred1D = unityD;
  pred2D = scos1D; 
  pred3D = sqrt(scos1D); 
  pred4D = scos1D.*vsec1D; 
  pred5D = (scos1D.*vsec1D).^2; 
  pred6D = scos1D.*T_highD; 
  pred7D = scos1D.*T_medD; 
  pred8D = vsec1D.*T_highD; 
  pred9D = vsec1D.*T_medD; 
  predsD = [pred1D pred2D pred3D pred4D pred5D pred6D pred7D pred8D pred9D];
elseif iVersFit == 2
  %% Sergio 2023
  pred1D = unityD;
  pred2D = scos1D; 
  pred3D = pred2D.*pred2D; 
  pred4D = pred2D.*vsec1D; 
  pred5D = pred2D.*T_highD/100; 
  pred6D = suncosD;
  predsD = [pred1D pred2D pred3D pred4D pred5D pred6D];
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NIGHT
if iVersKC > 0
  T_highN = ones(length(night),1) * sum(p.ptemp(01:05,:))/5; T_highN = T_highN(:);
  T_medN = ones(length(night),1) * sum(p.ptemp(06:15,:))/10; T_medN = T_medN(:);
  T_lowN = ones(length(night),1) * sum(p.ptemp(16:25,:))/10; T_lowN = T_lowN(:);
  o3N    = ones(length(night),1) * duO3;                     o3N     = o3N(:); o3N = o3N.^(2.1125);

%  if strfind(opts.range,'night')
%  disp('using Chris ThighN'); T_highN = chris2N.thbi';
%  end
  disp('using Chris ThighN'); T_highN = chris2N.thbi';
  suncosN = cos(sangles(night)*pi/180)' * ones(1,length(profind));       suncosN = suncosN(:);  %suncosN = abs(suncosN);
  scos1N  = cos(sang1(night)*pi/180)' * ones(1,length(profind));         scos1N  = scos1N(:);   %scos1N = abs(scos1N);
  vsec1N  = sec(vang1(night)*pi/180)' * ones(1,length(profind));         vsec1N  = vsec1N(:);
  scos2N  = cos(sang2(ind_night_for_sang2_vang2)*pi/180)';
  vsec2N  = sec(vang2(ind_night_for_sang2_vang2)*pi/180)';
  scos1N = scos2N;
  vsec1N = vsec2N;
  o3N    = ones(length(night),1) * duO3/100;                        o3N  = o3N(:); o3N = o3N.^(2.07625).*vsec1N;
  n2oN   = ones(length(night),1) * duN2O/100;                       n2oN = n2oN(:); n2oN = n2oN.^(0.25).*vsec1N.*(0.5);
  coN    = ones(length(night),1) * duCO/100;                        coN  = coN(:); coN = coN.^(6.7).*vsec1N;
  
  unityN = ones(size(T_highN));
  
  if iVersFit == 0
    %% Scott/Sergio/Strow/Manuel 2005 GRL
    pred1N = unityN;
    pred2N = scos1N; 
    pred3N = pred2N.*pred2N; 
    pred4N = pred2N.*vsec1N; 
    pred5N = pred2N.*T_highN; 
    pred6N = suncosN;
    predsN = [pred1N pred2N pred3N pred4N pred5N pred6N];
  elseif iVersFit == 1
    %% Marco/Manuel 2019
    pred1N = unityN;
    pred2N = scos1N; 
    pred3N = sqrt(scos1N); 
    pred4N = scos1N.*vsec1N; 
    pred5N = (scos1N.*vsec1N).^2; 
    pred6N = scos1N.*T_highN; 
    pred7N = scos1N.*T_medN; 
    pred8N = vsec1N.*T_highN; 
    pred9N = vsec1N.*T_medN; 
    predsN = [pred1N pred2N pred3N pred4N pred5N pred6N pred7N pred8N pred9N];
  elseif iVersFit == 2  
    %% Sergio 2023
    pow = 2;
    pow = 1;
    pow10 = 0.5;
  
    pred1N = unityN;
    pred2N = T_highN/250;         pred2N = pred2N.^pow; 
    pred3N = (T_highN-T_medN)/25; pred3N = pred3N.^pow; 
    pred4N = (T_highN-T_lowN)/25; pred4N = pred4N.^pow; 
    pred5N = suncosN; 
    pred6N = vsec1N;
    pred7N = pred5N.*pred2N.*vsec1N.^2;
    pred8N = pred5N.*pred3N.*vsec1N.^2; 
    pred9N = pred5N.*pred4N.*vsec1N.^2; 
    %pred7N = pred5N.*pred2N./vsec1N;
    %pred8N = pred5N.*pred3N./vsec1N;
    %pred9N = pred5N.*pred4N./vsec1N;
    pred10N = (abs(suncosN.*vsec1N)).^pow10;
    pred11N = o3N;
    pred12N = n2oN;
    pred13N = coN;
    predsN = [pred1N pred2N pred3N pred4N pred5N pred6N pred7N pred8N pred9N pred10N pred11N pred12N pred13N];
  end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wmax = 2;
wmin = 1;

iDoALLorONE = -1; %% do one chan
iDoALLorONE = +1; %% do all chan

if iDoALLorONE > 0
  %%%%%%%%%%%% this does all chans %%%%%%%%%%%%
  ixall = 1:mm;
  iPlot = -1;
elseif iDoALLorONE < 0
  %%%%%%%%%%%% this does ONE chan %%%%%%%%%%%%
  ix = find(fairs >= 2361,1);
  ix = find(fairs >= 2330,1); 
  if strfind(opts.range,'day')
    sarta_debug = load('debug_fit_extended_nonLTE_all_generic_x_ichan_184_filtered_with_tprof.mat'); %% comes from running "fit_extended_nonLTE_all_generic_x.m"  with ix > 0
  elseif strfind(opts.range,'night')
    sarta_debug = load('debug_night_fit_extended_nonLTE_all_generic_x_ichan_184_filtered_with_tprof.mat'); %% comes from running "fit_extended_nonLTE_all_generic_x.m"  with ix > 0
  end
  
  plot(1:length(profind),sarta_debug.tprof(1,1:114:length(profind)*114),'bx-',1:length(profind),p.ptemp(1,:),'r.-')   %% 114 = length(unique(vangles)) * length(unique(sangles)) = 6 x 19; length(profind)*114 = 47*(6*19)
  plot(1:length(profind),sarta_debug.tprof(1,1:114:length(profind)*114)-p.ptemp(1,:),'r.-')   %% 114 = length(unique(vangles)) * length(unique(sangles)) = 6 x 19; length(profind)*114 = 47*(6*19)
  
  plot(1:length(profind),sort(sarta_debug.tprof(1,1:114:length(profind)*114)),'bx-',1:length(profind),sort(p.ptemp(1,:)),'r.-')   %% 114 = length(unique(vangles)) * length(unique(sangles)) = 6 x 19; length(profind)*114 = 47*(6*19)
  plot(1:length(profind),sort(sarta_debug.tprof(1,1:114:length(profind)*114)) - sort(p.ptemp(1,:)),'r.-')   %% 114 = length(unique(vangles)) * length(unique(sangles)) = 6 x 19; length(profind)*114 = 47*(6*19)
  iPlot = +1;
  ixall = ix;
end

for iii = 1 : length(ixall)
  ii = ixall(iii);

  %{
  %% Chris has this, but I think he has already subset for day, so ignore this piece of code!
  deltaradAll = squeeze(nlte(ii,:,:) - lte(ii,:,:));
  deltaradAll = deltaradAll(:);
  deltaRi = deltaradAll;
  weightDAll = ones(size(deltaRi));
  rmax  = max(abs(deltaRi));
  rmin  = min(abs(deltaRi));
  m     = (wmax - wmin)/(rmax - rmin);
  b     = wmax - m*rmax;
  weightDAll = m*deltaRi + b;
  weightD = weightDAll(ind_day_for_sang2_vang2);
  %}

  %% I have this
  deltaBT_D = squeeze(tnlte(ii,day,:) - tlte(ii,day,:));
  deltaBT_D = deltaBT_D(:);
  %% now Chris actually goes to look for deltaBT_D > 0.2 before proceeding, but I'm not doing that (for these tests)
  deltaradD = squeeze(nlte(ii,day,:) - lte(ii,day,:));
  deltaradD = deltaradD(:);
  deltaRi = deltaradD;
  weightD = ones(size(deltaRi));
  rmax  = max(abs(deltaRi));
  rmin  = min(abs(deltaRi));
  m     = (wmax - wmin)/(rmax - rmin);
  b     = wmax - m*rmax;
  weightD = m*deltaRi + b;


%  weightD = ones(size(deltaradD));

  %whos T_high suncos scos1 vsec1 deltarad unity preds
  %% deltarad = preds x coef ==>
  [~,numD] = size(predsD);
  coefD(:,ii) = ((weightD*ones(1,numD)).*predsD) \ (weightD.*deltaradD);

  tpred = lte(ii,day,:);      tpred = tpred(:); tpred = tpred+predsD*coefD(:,ii); tpred = rad2bt(fairs(ii),tpred);
  tactual = tnlte(ii,day,:);  tactual = tactual(:);
  tltex   = tlte(ii,day,:);   tltex   = tltex(:);

  moo = unique(suncosD);
  for mmm = 1 : length(moo)
    aooo = find(suncosD == moo(mmm));
    biasdiffD(mmm) = mean(tactual(aooo)-tpred(aooo));
    stddiffD(mmm) = std(tactual(aooo)-tpred(aooo));
    meanD(mmm) = mean(tactual(aooo)-tltex(aooo));
    stdD(mmm)  = std(tactual(aooo)-tltex(aooo));
  end
  meanallD(ii) = mean(tactual-tltex);
  stdallD(ii)  = std(tactual-tltex);
  biasallD(ii) = mean(tactual-tpred);
  stddallD(ii) = std(tactual-tpred);

  if iPlot > 0
    figure(1); clf;  plot(1:length(deltaradD),deltaradD,'b',1:length(deltaradD),predsD*coefD(:,ii),'r'); 
      title(['DAY ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
    figure(2); clf;  plot(1:length(deltaradD),deltaradD - predsD*coefD(:,ii),'b')
      title(['DAY ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
    
    blahD = deltaradD;
    figure(3); clf;  plot(acos(suncosD)*180/pi,blahD,'bx',acos(suncosD)*180/pi,predsD*coefD(:,ii),'r.'); xlabel('Solzen angles'); ylabel('NLTE-LTE : kcarta vs pred');
      title(['DAY ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
    figure(4); clf;  plot(acos(suncosD)*180/pi,blahD-predsD*coefD(:,ii),'b.'); xlabel('Solzen angles'); ylabel('NLTE-LTE : kcarta-pred');
      title(['DAY ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
    figure(5); clf;  plot(acos(suncosD)*180/pi,100*(blahD-predsD*coefD(:,ii))./(blahD),'b.'); xlabel('Solzen angles'); ylabel('100*(kcarta-pred)/kcarta');
      title(['DAY ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
  
    figure(6); clf;  plot(acos(suncosD)*180/pi,tactual-tltex,'b.',acos(suncosD)*180/pi,tactual-tpred,'r.'); xlabel('Solzen angles'); ylabel('NLTE-LTE : BT kcarta-pred');
      title(['DAY ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
  
    figure(7); clf;  plot(acos(suncosD)*180/pi,tactual-tpred,'r.'); xlabel('Solzen angles'); ylabel('NLTE-LTE : BT kcarta-pred');
      title(['DAY ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
  
    figure(8); clf;  
      subplot(211); errorbar(acos(moo)*180/pi,meanD,stdD,'r.');                           ylabel('d(BT) kcarta signal'); plotaxis2;
      title(['DAY NLTE-LTE dBT(K) : ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
      subplot(212); errorbar(acos(moo)*180/pi,biasdiffD,stddiffD,'r.'); xlabel('Solzen angles'); ylabel('d(BT) actual-pred'); plotaxis2;
  end

  if iPlot == 1 & strfind(opts.range,'day')
    figure_9_debug_day
    error(';lkaf;lakf')
  end

  %pause(0.1)
  %%%%%%%%%%%%%%%%%%%%%%%%%    
  if iVersKC > 0
    deltaradN = squeeze(nlte(ii,night,:) - lte(ii,night,:));
    deltaradN = deltaradN(:);
    weightN = ones(size(deltaradD));

    deltaRi = deltaradN;
    rmax  = max(abs(deltaRi));
    rmin  = min(abs(deltaRi));
    m     = (wmax - wmin)/(rmax - rmin);
    b     = wmax - m*rmax;
    weightN = m*deltaRi + b;

    %whos T_high suncos scos1 vsec1 deltarad unity preds
    %% deltarad = preds x coef ==>
    [~,numN] = size(predsN);
    coefN(:,ii) = ((weightN*ones(1,numN)).*predsN) \ (weightN.*deltaradN);
  
    tpred = lte(ii,night,:);      tpred = tpred(:); tpred = tpred+predsN*coefN(:,ii); tpred = rad2bt(fairs(ii),tpred);
    tactual = tnlte(ii,night,:);  tactual = tactual(:);
    tltex   = tlte(ii,night,:);   tltex   = tltex(:);
  
    moo = unique(suncosN);
    for mmm = 1 : length(moo)
      aooo = find(suncosN == moo(mmm));
      biasdiffN(mmm) = mean(tactual(aooo)-tpred(aooo));
      stddiffN(mmm) = std(tactual(aooo)-tpred(aooo));
      meanN(mmm) = mean(tactual(aooo)-tltex(aooo));
      stdN(mmm)  = std(tactual(aooo)-tltex(aooo));
    end
    meanallN(ii) = mean(tactual-tltex);
    stdallN(ii)  = std(tactual-tltex);
    biasallN(ii) = mean(tactual-tpred);
    stddallN(ii) = std(tactual-tpred);
  
    if iPlot > 0
      figure(1); hold;  plot(1:length(deltaradN),deltaradN,'g',1:length(deltaradN),predsN*coefN(:,ii),'y'); 
        title(['DAY/NIGHT ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
      figure(2); hold;  plot(1:length(deltaradN),deltaradN - predsN*coefN(:,ii),'g')
        title(['DAY/NIGHT ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
    
      blahN = deltaradN;
      figure(3); hold;  plot(acos(suncosN)*180/pi,blahN,'bx',acos(suncosN)*180/pi,predsN*coefN(:,ii),'r.'); xlabel('Solzen angles'); ylabel('NLTE-LTE : kcarta vs pred');
        title(['DAY/NIGHT ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ]); plotaxis2;     
      figure(4); hold;  plot(acos(suncosN)*180/pi,blahN-predsN*coefN(:,ii),'b.'); xlabel('Solzen angles'); ylabel('NLTE-LTE : kcarta-pred');
        title(['DAY/NIGHT ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ]); plotaxis2;     
      figure(5); hold;  plot(acos(suncosN)*180/pi,100*(blahN-predsN*coefN(:,ii))./(blahN),'b.'); xlabel('Solzen angles'); ylabel('100*(kcarta-pred)/kcarta');
        title(['DAY/NIGHT ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ]); plotaxis2;     
      ylim([-1 +1]*100);
    
      figure(6); hold;  plot(acos(suncosN)*180/pi,tactual-tltex,'b.',acos(suncosN)*180/pi,tactual-tpred,'r.'); xlabel('Solzen angles'); ylabel('NLTE-LTE : BT kcarta-pred');
        title(['DAY/NIGHT ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ]); plotaxis2;
        hl = legend('True NLTE-LTE','True NLTE-Pred NLTE','location','best'); 
    
      figure(7); hold;  plot(acos(suncosN)*180/pi,tactual-tpred,'r.'); xlabel('Solzen angles'); ylabel('NLTE-LTE : BT kcarta-pred');
        title(['DAY/NIGHT ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ]); plotaxis2; ylim([-1 +1]*2)
    
      figure(8); 
        subplot(211); hold; errorbar(acos(moo)*180/pi,meanN,stdN,'r.');                           ylabel('d(BT) kcarta signal'); plotaxis2;
        title(['DAY/NIGHT NLTE-LTE dBT(K) : ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
        subplot(212); hold; errorbar(acos(moo)*180/pi,biasdiffN,stddiffN,'r.'); xlabel('Solzen angles'); ylabel('d(BT) actual-pred'); plotaxis2;
    
        clf;
        subplot(211); errorbar(acos(moo)*180/pi,meanN,stdN,'r.');                           ylabel('d(BT) kcarta signal'); plotaxis2;
        title(['NIGHT NLTE-LTE dBT(K) : ch ' num2str(ii) ' fairs = ' num2str(fairs(ii)) ])
        subplot(212); errorbar(acos(moo)*180/pi,biasdiffN,stddiffN,'r.'); xlabel('Solzen angles'); ylabel('d(BT) : actual-pred'); plotaxis2;
    end       %% if iPlot > 0

    if iPlot == 1 & strfind(opts.range,'night')
      figure_9_debug_night
      error(';lkaf;lakf')
    end
 
  end         %% if iVersKC > 0  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iPlot < 0
  if iVersKC > 0

    figure(9); clf
    errorbar(fairs,meanallD,stdallD,'b'); hold on
    errorbar(fairs,biasallD,stddallD,'c'); hold on
    errorbar(fairs,meanallN,stdallN,'r'); hold on
    errorbar(fairs,biasallN,stddallN,'m'); hold off
    hl = legend('Signal(D)=NLTE(D)-LTE(D) kc','fit(D) bias/std','Signal(N)=NLTE(N)-LTE(N) kc','fit(N) bias/std','location','best');
    xlim([2200 2400]);

    figure(10); clf
    subplot(211); plot(fairs,meanallD,'b',fairs,biasallD,'c',fairs,meanallN,'r',fairs,biasallN,'m');
    subplot(212); plot(fairs,stdallD,'b', fairs,stddallD,'c',fairs,stdallN,'r',fairs,stddallN,'m');
    xlim([2200 2400]);  
  
    figure(10); clf
    subplot(211); plot(fairs,meanallD,'b',fairs,meanallN*10,'r'); plotaxis2; ylabel('mean'); title('KCARTA NLTE-LTE'); hl = legend('D','N*10'); xlim([2200 2400]); 
       ylim([-5 10])
    subplot(212); plot(fairs,stdallD,'b',fairs,stdallN*10,'r');   plotaxis2; ylabel('std');                            hl = legend('D','N*10'); xlim([2200 2400]);
       ylim([0 5])

    figure(11); clf
    subplot(211); plot(fairs,biasallD,'c',fairs,biasallN,'m'); plotaxis2; ylabel('mean'); title('fit-KCARTA (NLTE-LTE)'); hl = legend('D','N'); xlim([2200 2400]);
      ylim([-0.4 +0.2])
    subplot(212); plot(fairs,stddallD,'c',fairs,stddallN,'m'); plotaxis2; ylabel('std');                                  hl = legend('D','N'); xlim([2200 2400]);
      ylim([0 0.75])
  
    figure(12); plot(fairs,coefD,'linewidth',2); xlim([2200 2400]); plotaxis2; title('Day coeffs');
    figure(13); plot(fairs,coefN,'linewidth',2); xlim([2200 2400]); plotaxis2; title('Night Coeffs')

  else
    figure(9); clf
    errorbar(fairs,meanallD,stdallD,'b'); xlim([2200 2400]); hold on; 
    errorbar(fairs,biasallD,stddallD,'c'); xlim([2200 2400]); hold off
  
    figure(10); clf
    subplot(211); plot(fairs,meanallD,'b',fairs,biasallD,'c'); xlim([2200 2400]);
    subplot(212); plot(fairs,stdallD,'b', fairs,stddallD,'c'); xlim([2200 2400]);
  
    figure(10); clf
    subplot(211); plot(fairs,meanallD,'b'); xlim([2200 2400]);
    subplot(212); plot(fairs,stdallD,'b'); xlim([2200 2400]);
    figure(11); clf
    subplot(211); plot(fairs,biasallD,'c'); plotaxis2; xlim([2200 2400]);
    subplot(212); plot(fairs,stddallD,'c'); plotaxis2; xlim([2200 2400]);
  
    figure(12); plot(fairs,coefD,'linewidth',2); xlim([2200 2400]); plotaxis2; title('Day coeffs');
  end
end

figure(14); plot(suncosD,scos1D,'.',suncosN,scos1N,'.'); plotaxis2; 
xlabel('suncos = solar angle at GND'); ylabel('scos = solar ang at TOA')

figure(15); plot(sarta.vchan,sarta.cofn,'linewidth',2); plotaxis2; title('Official SARTA coeffs')
