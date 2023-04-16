%% this is the cngwat(SERGIO) --> cngwat(PY) conversion
dmeXX    = [10.00  20.00  30.00  40.00  50.00  75.00  100.0   150.0  200.0];
factorXX = [0.170  0.125  0.080  0.060  0.040  0.030  0.0225  0.015  0.015];

%% cp ~sergio/SARTA_CLOUDY/v108_Src_pclsam_ctype_HG3_PYang_1layer/make_rtp.m

dogoodchan; freqsNchans;

weff = [20]; %% typical water particles
%deff = input('enter dme : ');
%cprtop = input('enter cprtop : ');
%cngwat = input('enter cngwat : ');

xx = input('enter [deff cprtop cngwat] : ');
deff = xx(1);
cprtop = xx(2);
cngwat = xx(3);

% filename of klayers executable 
klayers = '/asl/packages/klayers/Bin/klayers_airs'; 

%% these should have same database; one for 2378 chans, other for 2834 chans
sarta0=...
  '/home/sergio/SARTA_CLOUDY/Bin/testPYsarta_apr08_m130_iceaggr_waterdrop_volz_1_slabcloud_hg3';
sarta0= ... 
  '/home/sergio/SARTA_CLOUDY/Bin/sarta_dec05_iceaggr_waterdrop_volz_1_Mar08'; 

sartax=...
  '/home/sergio/SARTA_CLOUDY/Bin/testPY1sarta_apr08_m130_iceaggr_waterdrop_volz_1_slabcloud_hg3';

fin='/home/sergio/MATLABCODE/DCC_SIMPLE/pin_feb2002_sea_airsnadir_ip.so2.rtp';
klayerer = ['!' klayers ' fin=' fin ' fout=junk.op.rtp '];
eval(klayerer)
[h,ha,p,pa] = rtpread('junk.op.rtp');
mm = mmwater_rtp(h,p);

[h,ha,p,pa] = rtpread(fin);
p.rcalc = zeros(2378,49);
hot = find(p.stemp >= p.stemp(1) & mm >= mm(1))
h.pfields = 5;

wmult0   = 1.0;
poffset0 = 0;
wmult   =  1.0;
poffset = poffset0;
dme = 10 : 20 : 150;
numtimes = length(wmult) * length(poffset) * length(dme);

%for pp = 1 : length(hot)
for pp = 1 : 1
  [h,phot] = subset_rtp(h,p,h.glist,1:2378,hot(pp));

  index_trop = tropopause_rtp(h,phot);
  trop_press = phot.plevs(index_trop);
  fprintf(1,'tropopause is at %9.4f \n',trop_press);

  [headyy,profyy] = replicate_rtp_headprof(h,phot,1,numtimes);
  profyy.scanang = 22.5 * ones(size(profyy.scanang));
  profyy.satzen = 22.5 * ones(size(profyy.scanang));

  profzz = profyy;

  profyy.ctype   = 201 * ones(size(profyy.stemp));   %%% cirrus <-- Scott
  profyy.cpsize  = deff * ones(size(profyy.stemp));   
  profyy.cngwat  = cngwat * ones(size(profyy.stemp));  
  profyy.cfrac   = 1 * ones(size(profyy.stemp));  

  iCnt = 0;
  for dd = 1 : length(dme)
    deffx = dme(dd);
    for oo = 1 : length(poffset)
      for ww = 1 : length(wmult)
        iCnt = iCnt + 1;
        profyy.cpsize(iCnt)  = deffx;
        profyy.cprtop(iCnt) = trop_press + poffset(oo);
        profyy.cprbot(iCnt) = trop_press + poffset(oo) + 100;
        profyy.cprtop(iCnt) = cprtop;
        profyy.cprbot(iCnt) = cprtop + 100;
        doink = profyy.plevs(:,iCnt); 
        boink = find(doink <= profyy.cprbot(iCnt)); 
        profyy.gas_1(boink,iCnt) =  profyy.gas_1(boink,iCnt) * wmult(ww);
        end
      end
    end

  rtpwrite('junk.ip.rtp',h,ha,profyy,pa);
  klayerer = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp '];
  eval(klayerer)

  sartaer = ['!time ' sarta0 ' fin=junk.op.rtp fout=junk.rad.rtp '];
  eval(sartaer)
  [hx,hxa,px0,pxa] = rtpread('junk.rad.rtp');
  tcal0 = rad2bt(hx.vchan,px0.rcalc);

  for dd = 1 : iCnt
    deffx = px0.cpsize(dd);
    xconv = interp1(dmeXX,factorXX,deffx);
     %% rough cngwat --> OD conversion
    px0.cngwat(dd) = px0.cngwat(dd) *xconv;
    end    
  rtpwrite('junk.op.rtp',hx,ha,px0,pa);
  sartaer = ['!time ' sartax ' fin=junk.op.rtp fout=junk.rad.rtp '];
  eval(sartaer)
  [hx,hxa,pxx,pxa] = rtpread('junk.rad.rtp');
  tcalX = rad2bt(hx.vchan,pxx.rcalc);

  poink = find(hx.vchan >= 1419.1); poink = intersect(poink,g); 
  ii1419 = poink(1);

  figure(1); clf
  plot(tcal0(cind1(5),:),tcal0(cind1(5),:)-tcal0(ii1419,:),'bo',...
       tcalX(cind1(5),:),tcalX(cind1(5),:)-tcalX(ii1419,:),'ro','Linewidth',2)
  title(num2str(hot(pp))); xlabel('BT1231'); ylabel('BT1231-BT1419')

  figure(2); clf
  plot(hx.vchan,tcal0,'b','Linewidth',2); hold on
  plot(hx.vchan,tcalX,'r--','linewidth',2); axis([650 1650 180 250]); grid on
  title(num2str(hot(pp))); xlabel('wavenumber'); ylabel('BT')

  figure(3); 
    nn = 1 : profyy.nlevs(1); 
    zonkA = interp1(log10(profyy.plevs(nn,1)),profyy.ptemp(nn,1),log10(profyy.cprtop(1)));
    zonkB = interp1(log10(profyy.plevs(nn,1)),profyy.ptemp(nn,1),log10(profyy.cprbot(1)));
    subplot(121);  semilogy(profyy.ptemp(nn,1),profyy.plevs(nn,1));
      set(gca,'ydir','reverse')
      axis([180 300 1e-1 1000]);
%    line([zonkA zonkA]);
%    line([zonkB zonkB]);
    fprintf(1,'T(cprtop) = %8.5f T(cprbot) = %8.5f \n',zonkA,zonkB); 
    subplot(122);  loglog(profyy.gas_1(nn,:),profyy.plevs(nn,1));
    zonkA = interp1(log10(profyy.plevs(nn,1)),profyy.ptemp(nn,1),log10(profyy.cprtop(1)));
    zonkB = interp1(log10(profyy.plevs(nn,1)),profyy.ptemp(nn,1),log10(profyy.cprbot(1)));
      set(gca,'ydir','reverse')
      axis([1e-6 1e-3 1e-1 1000]);
  title(num2str(hot(pp))); 

  end

