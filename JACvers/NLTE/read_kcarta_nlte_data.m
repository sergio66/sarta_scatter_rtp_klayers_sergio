if iVersKC == 0
  boo0 = ['/home/sergio/KCARTA/NONLTE_PRODUCTION/VT_48PROFILES_120_400ppmv_v121_H16_NLTEH16_Mar2020_NewConvolvers/Results/CONV_Results/'];  
  boo0 = ['/home/sergio/KCARTA/NONLTE_PRODUCTION/VT_48PROFILES_120_400ppmv_v120_H16_NLTEH16_Dec2018/Results/CONV_Results/'];                %% since I compare to SARTA May19
  %                  solar: [0 40 60 80 85 90 0 40 60 80 85 90 0 40 60 80 85 90 0 40 60 80 85 90 0 40 60 80 85 90 0 40 60 80 85 90]
  %                  viewer: [0 0 0 0 0 0 10 10 10 10 10 10 20 20 20 20 20 20 30 30 30 30 30 30 40 40 40 40 40 40 50 50 50 50 50 50]
elseif iVersKC == 1
  boo0 = ['/asl/s1/sergio/NLTE_CALCS/VT_48PROFILES_120_SRCv1.21_400ppmv_H16_Apr2021_NewNLTEProfiles/CONV_Results/'];
 %                    solar: [0 10 20 30 40 50 60 70 80 85 87 90 92 94 96 98 100 105 120 0 10 20 30 40 50 60 70 80 85 87 90 92 94 96 98 100 105 120 0 10 20 30 40 ... ]
 %                   viewer: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 20 20 20 20 20 20 20 20 20 20 20 20 ... ]
end

for iii = 1 : length(profind)
  ii = profind(iii);
  if mod(iii,10) == 0
    fprintf(1,'x');
  else
    fprintf(1,'.');
  end
  boo = [boo0 '/vt' num2str(ii) '.mat'];

  x = load(boo);
  nlte(iii,:,:) = x.nlte_airs;
  lte(iii,:,:)  = x.lte_airs_Bfunny4um;

end
fprintf(1,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% let solzen  --> sangles = solar zenith = solar angle at GND
%% let scanang --> vangles = scanang at satellite
%%
%% let satzen --> vang1 = vang at TOA
%% let solang --> sang1 = satellite angle at TOA

sathgt = 710000;  %% AIRS should be 705000 but Chris uses 710000????
fairs   = x.fairs;
sangles = x.solar;
vangles = x.viewer;
junk    = mean(0.5*(p.palts(1,:)+p.palts(2,:))) * ones(size(sangles));   %% recall palts is at pressure levels, while we want pressure layers


chris1 =  load('/home/chepplew/data/sarta/prod_2022/airs_l1c/apr2021//nonLTE/xnlte_merged.mat','zprof');
chris2D = load('debug_fit_extended_nonLTE_all_generic_x_ichan_184_filtered_with_tprof.mat','sangles','vangles','vang1','sang1','thbi');
chris2N = load('debug_night_fit_extended_nonLTE_all_generic_x_ichan_184_filtered_with_tprof.mat','sangles','vangles','vang1','sang1','thbi');
if strfind(opts.range,'day')
  chris2 = chris2D;
else
  chris2 = chris2N;
end
plot(chris1.zprof(1,1:114:114*48))
plot(1:48,chris1.zprof(1,1:114:114*48),1:48,p.palts(1,:))
plot(1:48,chris1.zprof(1,1:114:114*48),'b.-',1:48,0.5*(p.palts(1,:)+p.palts(2,:)))

junk_palt1 = chris1.zprof(1,1:114:114*48);
junk    = mean(junk_palt1) * ones(size(sangles));

clear x boo ii

nlte = permute(nlte,[2 3 1]);
lte  = permute(lte,[2 3 1]);

tnlte = rad2bt(fairs,nlte);
tlte  = rad2bt(fairs,lte);

[mm,nn,oo] = size(tnlte);
fprintf(1,'size(tlnte) = numchan %4i x vangXsang %4i x numprof %4i \n',[mm,nn,oo])
fprintf(1,'[nn length(unique(vangles))*length(unique(sangles))] = %4i %4i %4i \n',[nn length(unique(vangles))*length(unique(sangles))])

vang1 = vaconv(vangles,sathgt*ones(size(vangles)),junk);   %% satellite zenith angle at ground
sang1 = saconv(sangles,junk);                              %% solar ang at TOA

vang2 = [];
sang2 = [];
for ii = 1 : oo
  junk2 = 0.5*(p.palts(1,ii)+p.palts(2,ii)) * ones(size(sangles));
  junk2 = junk_palt1(ii) * ones(size(sangles));  %% use Chris H mystery heights
  mooV = vaconv(vangles,sathgt*ones(size(vangles)),junk2); 
  mooS = saconv(sangles,junk2);
  vang2 = [vang2 mooV];
  sang2 = [sang2 mooS];
end
clear moo*

plot(1:114,chris2.vang1(1:114),'b.-',1:114,vang1)
plot(1:114,chris2.vang1(1:114)-vang1,'b',1:114,chris2.vang1(1:114)-vang2(1:114),'r')
plot(1:114,chris2.sang1(1:114),'b.-',1:114,sang1)
plot(1:114,chris2.sang1(1:114)-sang1,'b',1:114,chris2.sang1(1:114)-sang2(1:114),'r')
