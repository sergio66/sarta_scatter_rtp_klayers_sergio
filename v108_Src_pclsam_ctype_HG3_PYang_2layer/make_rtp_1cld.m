weff = [20]; %% typical water particles
deff = [30];
cfrac = 1.0; cfrac = 0.75; cfrac = 1.0;
cngwat = 500;
cprtop = 133;
cprbot = 380;
iMax = 3;

% filename of klayers executable 
klayers = '/asl/packages/klayers/Bin/klayers_airs'; 

%% should have the same database; one for 2378 chans, other for 2834 chans
sarta0= ... 
  '/home/sergio/SARTA_CLOUDY/Bin/sarta_dec05_iceaggr_waterdrop_volz_1_Mar08'; 

%% single cloud code given in late Nov 2008
sarta1x=...
  '/home/sergio/SARTA_CLOUDY/Bin/testPY1sarta_apr08_m130_iceaggr_waterdrop_volz_1_slabcloud_hg3';

%% two cloud code given in late Jan 2009
sarta2x=...
  '/home/sergio/SARTA_CLOUDY/Bin/testPY2sarta_apr08_m130_iceaggr_waterdrop_volz_1_slabcloud_hg3';

gprofer = ['gprof ' sarta2x ' gmon.out >& ugh; more ugh'];  %% <<<<-------

fip='/home/sergio/MATLABCODE/DCC_SIMPLE/pin_feb2002_sea_airsnadir_ip.so2.rtp';

for dd = 1 : length(deff)
  [h,ha,p,pa] = rtpread(fip);
  p.rcalc = zeros(2378,49);
  %hot = find(p.stemp >= 290);
  hot = 1;   %% use TROPICAL profile
  h.pfields = 5;
  [h,phot] = subset_rtp(h,p,h.glist,1:2378,hot);

  iCnt = 0;
  ss = 2.5;  
  ss = 22.5;  
  pp = 1;
  [h,p] = subset_rtp(h,phot,h.glist,1:2378,pp);
  [index_trop] = tropopause_rtp(h,p);
  trop_press = p.plevs(index_trop);
  fprintf(1,'tropopause is at %9.4f \n',trop_press);

  p.ctype  = 201;   %%% cirrus <-- Scott's new code 
  p.cpsize = deff(dd);
  p.cngwat = cngwat;
  p.cfrac  = cfrac;
  p.scanang = ss;
  p.satzen  = ss;

  p.cprtop = trop_press;
  p.cprtop = 80;
  p.cprbot = p.cprtop + 50;  %%% <<<<--- changed this from DCC 150 mb

  % from /asl/packages/sartaV105/yukyung_readme.txt 
  %      prof.udef(11,:) = cngwat 
  %      prof.udef(12,:) = cpsize 
  %      prof.udef(13,:) = cprtop 
  %      prof.udef(14,:) = cprbot 
  %      prof.udef(15,:) = cfrac
  %      prof.udef(16,:) = "cfrac12", fraction of FOV containing both clouds 
  %      prof.udef(17,:) = ctype {currently not used} 
  %      prof.udef(18,:) = cemis for cloud 2 
  p.udef(11) = 1.5*4000; %%% see Zimmer2008 paper 1.5g/m3, cloud from 1 to 6 km
  p.udef(12) = weff;
  p.udef(13) = p.cprbot + 50;
  p.udef(14) = h2p(2000);
  p.udef(15) = 1.0;
  p.udef(16) = 1.0;
  p.udef(17) = 101;      %%% water <--- Scott's new code

  cngwat = 0 : 0.2 : 20;
  cngwat = (0 : 0.1 : 2)*10; cngwat = [cngwat 300];  %% for DCC
  cngwat = (0 : 1 : 2)*10; cngwat = [cngwat 300];  %% for DCC
  %p.emis = ones(size(p.emis));
  %p.rho = zeros(size(p.rho));
  [headzz,profzz] = replicate_rtp_headprof(h,p,1,length(cngwat));
  profzz.cngwat = cngwat;
  profzz.cprtop = ones(size(profzz.cprtop)) * 200;
  profzz.cprbot = profzz.cprtop + 100;

  profzz.udef(11,:) = 1.5*4000; %%% see Zimmer2008 paper 1.5g/m3, cloud from 1 to 6 km
  profzz.udef(12,:) = weff;
  profzz.udef(13,:) = profzz.cprbot + 50;
  profzz.udef(14,:) = h2p(2000);
  profzz.udef(15,:) = 1.0;
  profzz.udef(16,:) = 1.0;
  profzz.udef(17,:) = 101;      %%% water <--- Scott's new code

  %************************************************************************

  %% now mess things up  ... put 4 um dust
  profzz.udef(11,:) = 0.0;
  profzz.udef(15,:) = 0.0;
  profzz.udef(16,:) = 0.0;
  profzz.ctype = 301 * ones(size(profzz.ctype));;      %%% dust <--- 
  profzz.cprtop = 800 * ones(size(profzz.ctype));;      
  profzz.cprbot = 850 * ones(size(profzz.ctype));;      
  profzz.cpsize = 4.0 * ones(size(profzz.ctype));;      

  %% now mess things up  ... put water
  profzz.udef(11,:) = 0.0;
  profzz.udef(15,:) = 0.0;
  profzz.udef(16,:) = 0.0;
  profzz.ctype = 101 * ones(size(profzz.ctype));;      %%% water <--- 
  profzz.cprtop = 800 * ones(size(profzz.ctype));;      
  profzz.cprbot = 850 * ones(size(profzz.ctype));;      
  profzz.cpsize = 20.0 * ones(size(profzz.ctype));;     
  PX = 'W';  GammaOrLogX = 1; 
  TTX = 222;  rhoX = 1.00; gammaX = 6;   %% for water

  %% now mess things up  ... put ice
  profzz.udef(11,:) = 0.0;
  profzz.udef(15,:) = 0.0;
  profzz.udef(16,:) = 0.0;
  profzz.ctype = 201 * ones(size(profzz.ctype));;      %%% ice <--- 
  profzz.cprtop = 350 * ones(size(profzz.ctype));;      
  profzz.cprbot = 370 * ones(size(profzz.ctype));;      
  profzz.cpsize = 75.0 * ones(size(profzz.ctype));;    
  PX = 'I';  GammaOrLogX = 1; 
  TTX = 2;  rhoX = 0.97; gammaX = 2;   %% for ice 

  profzz.cfrac(1) = 0;
  rtpwrite('junk.ip.rtp',headzz,ha,profzz,pa);
  klayerer = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp '];
  eval(klayerer)

  sartaer = ['!time ' sarta0 ' fin=junk.op.rtp fout=junk.rad.rtp '];
  eval(sartaer)
  [hh,hha,p0,ppa] = rtpread('junk.rad.rtp');

  iDeltaUnScale0 = +1; 
  raTau900Upper(dd,:) = ...
    loadmie_tau(900,profzz.cpsize(1),PX,TTX,GammaOrLogX,gammaX,1,iDeltaUnScale0);
  tau = raTau900Upper(1)*profzz.cngwat(2);

  %% so roughly 1g/m2 gives raTau900Upper(1)  in IR OD, which is about x2 smaller than VIS OD
  %p0.cngwat = p0.cngwat * 2.00 * raTau900Upper(1);  %%should roughly be what VIS OD we need
  if profzz.ctype(1) < 200
    %% adjust water so CNGWAT --> OD else keep things the same for dust
    %% should roughly be what VIS OD we need for 50 um waterparticles
    %% p0.cngwat = p0.cngwat * 0.75 * raTau900Upper(1);  
    %% should roughly be what VIS OD we need for 20 um waterparticles
    p0.cngwat = p0.cngwat * 1.0 * raTau900Upper(1);  
    p0.cfrac(1) = 1.0;
    rtpwrite('junk.op.rtp',hh,hha,p0,pa);
  elseif profzz.ctype(1) < 300
    %% adjust ice so CNGWAT --> OD else keep things the same for dust
    %% should roughly be what VIS OD we need for 50 um ice particles
     p0.cpsize = p0.cpsize * 0.75;
     p0.cngwat = p0.cngwat * 1.15 * raTau900Upper(1);  
    %% should roughly be what VIS OD we need for 65 um ice particles at 500 mb
    %% p0.cngwat = p0.cngwat * 1.0 * raTau900Upper(1);  
    p0.cfrac(1) = 1.0;
    rtpwrite('junk.op.rtp',hh,hha,p0,pa);
    end

  %p0 =  pclsam
  %p1x = PY1
  %p2x = PY2

  sartaer = ['!time ' sarta1x ' fin=junk.op.rtp fout=junk.rad.rtp '];
  eval(sartaer)
  [hh,hha,p1x,ppa] = rtpread('junk.rad.rtp');

  sartaer = ['!time ' sarta2x ' fin=junk.op.rtp fout=junk.rad.rtp '];
  %sartaer = ['!time gprof ' sarta2x ' fin=junk.op.rtp fout=junk.rad.rtp '];
  eval(sartaer)
  [hh,hha,p2x,ppa] = rtpread('junk.rad.rtp');

  figure(1); clf
  plot_wvn_lambda(hh.vchan,rad2bt(hh.vchan,p0.rcalc))
  hl = legend(num2str(cprtop'));
  set(hl,'FontSize',12)
  title('pclsam')
  axis([600 1650 190 300])

  figure(2); clf
  plot_wvn_lambda(hh.vchan,rad2bt(hh.vchan,p2x.rcalc))
  hl = legend(num2str(cprtop'));
  set(hl,'FontSize',12)
  title('ping yang 2cld version')
  axis([600 1650 190 300])

  figure(3); clf
  plot(hh.vchan,rad2bt(hh.vchan,p0.rcalc)-rad2bt(hh.vchan,p1x.rcalc)); grid on
  ylabel('\delta BT(K) = PCLSAM - PY1')

  figure(4); clf
  plot(hh.vchan,rad2bt(hh.vchan,p1x.rcalc)-rad2bt(hh.vchan,p2x.rcalc)); grid on
  ylabel('\delta BT(K) = PY1 - PY2')

  figure(5); clf
  plot(hh.vchan,...
    rad2bt(hh.vchan,p1x.rcalc(:,1))-rad2bt(hh.vchan,p0.rcalc(:,1)),...
       hh.vchan,...
    rad2bt(hh.vchan,p2x.rcalc(:,1))-rad2bt(hh.vchan,p0.rcalc(:,1)),...
       hh.vchan,...
    rad2bt(hh.vchan,p1x.rcalc(:,1))-rad2bt(hh.vchan,p2x.rcalc(:,1)))
  hl = legend('py1-pclasm','py2-pclsam','py1-py2','location','southwest')
  set(hl,'fontsize',12);
  axis([600 2800 -2 1]); grid on
  end
