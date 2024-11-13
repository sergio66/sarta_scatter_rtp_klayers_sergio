addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE/PLOTTER

disp('see /asl/packages/klayersV205/Data/glatm_16Aug2010.dat')
disp('see /home/sergio/KCARTA/NONLTE2/sergio/VT_48PROFILES_FOR_MANUEL/AtomicO_SABER_Mlynczak/add_othergases_arb_pressures.m')

mult = 0;
mult = 1;

%load /home/chepplew/data/sarta/SAF/saf_sub300_reg49_plv147.mat
%load saf_sub300_reg49_plv147.mat
%load saf_sub300_reg49_plv147_400ppm_v3
%load saf_sub600_reg49_plv147_400ppm
load saf_sub300_reg49_plv147_400ppm_v3.mat

if pcat.plevs(1,1) > pcat.plevs(137,1)
  disp('bloody heck have to flip every matrix')
  pcat.cc     = flipud(pcat.cc);
  pcat.ciwc   = flipud(pcat.ciwc);
  pcat.clwc   = flipud(pcat.clwc);
  pcat.plevs  = flipud(pcat.plevs);
  pcat.ptemp  = flipud(pcat.ptemp);
  pcat.gas_1  = flipud(pcat.gas_1);
  pcat.gas_2  = flipud(pcat.gas_2);
  pcat.gas_3  = flipud(pcat.gas_3);
  pcat.gas_5  = flipud(pcat.gas_5);
  pcat.gas_6  = flipud(pcat.gas_6);
  pcat.gas_9  = flipud(pcat.gas_9);
  pcat.gas_34 = flipud(pcat.gas_34);
end

scatter_coast(pcat.rlon,pcat.rlat,40,pcat.spres)

bad = find(pcat.spres < 100);
if length(bad) > 0
  fprintf(1,'found %3i profiles with spres < 100 \n',length(bad))
end

[mm,nn] = find(isnan(pcat.plevs));
if length(mm) ~= 0
  disp('pcat.plevs has NaN elements')
end

[mm,nn] = find(isnan(pcat.ptemp));
if length(mm) ~= 0
  disp('pcat.ptemp has NaN elements')
end

[mm,nn] = find(isnan(pcat.gas_1));
if length(mm) ~= 0
  disp('pcat.gas_1 has NaN elements')
end

[mm,nn] = find(isnan(pcat.gas_3));
if length(mm) ~= 0
  disp('pcat.gas_3 has NaN elements')
end

if ~isfield(pcat,'spres')
  disp('setting spres = plevs(147,:)')
  pcat.spres = pcat.plevs(147,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hcatx,pcatx] = subset_rtp(hcat,pcat,[1 2 3 5 6 9],[],[]);
pcatx.gas_1_press = pcatx.plevs * 100 .* pcatx.gas_1/1e6;   %% gas units are ppmv, so change to pressures
pcatx.altitude = hypsometric_p2h(pcatx.salti,pcatx.spres,pcatx.plevs,pcatx.ptemp,mult*pcatx.gas_1_press);
rtpwrite('saf_sub300_reg49_plv147.ip.rtp',hcatx,[],pcatx,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pcatx080 = pcatx;

for ii = 1 : 349
  boo = pcatx.plevs(:,ii);
  boo = abs(boo - 0.005);
  i005(ii) = find(boo == min(boo),1);
  boo = abs(boo - 1.0);
  i1(ii) = find(boo == min(boo),1);
end
%i005 = i005 - 1;

iuse = 28:147; whos iuse
iuse = [i005:2:48 49:147]; whos iuse

pcatx080.plevs = pcat.plevs(iuse,:);
pcatx080.ptemp = pcat.ptemp(iuse,:);
pcatx080.gas_1 = pcat.gas_1(iuse,:);
pcatx080.gas_2 = pcat.gas_2(iuse,:);
pcatx080.gas_3 = pcat.gas_3(iuse,:);
pcatx080.gas_5 = pcat.gas_5(iuse,:);
pcatx080.gas_6 = pcat.gas_6(iuse,:);
pcatx080.gas_9 = pcat.gas_9(iuse,:);
pcatx080.nlevs = 120 * ones(size(pcatx080.nlevs));
pcatx080.gas_1_press = pcatx080.plevs * 100 .* pcatx080.gas_1/1e6;   %% gas units are ppmv, so change to pressures
pcatx080.altitude = hypsometric_p2h(pcatx080.salti,pcatx080.spres,pcatx080.plevs,pcatx080.ptemp,mult*pcatx080.gas_1_press);
rtpwrite('saf_sub300_reg49_plv120.ip.rtp',hcatx,[],pcatx080,[]);

rtpdumper = ['! /home/sergio/git/rtp/rtpV221_150levs/utils/rtpdump -p -n 1 saf_sub300_reg49_plv147.ip.rtp >& ugh'];
eval(rtpdumper)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runner = ['!test_300_49.sc']; eval(runner)

[h120,ha,p120,pa]         = rtpread('saf_sub300_reg49_plv147_120km.rtp');   %% sergio klayers  0-120 km, 150 levels
[h080,ha,p080,pa]         = rtpread('saf_sub300_reg49_plv120_080km.rtp');   %% usual  klayers  0-080 km, 120 levels
[h080,ha,p080_150levs,pa] = rtpread('saf_sub300_reg49_plv147_080km.rtp');   %% sergio klayersx 0-080 km, 150 levels

mm080_150levs = mmwater_rtp(h080,p080_150levs);
mm080_120levs = mmwater_rtp(h080,p080); mm080_usual = mm080_120levs;
mm120 = mmwater_rtp(h120,p120);

iNumProf = length(mm080_usual);

plot(1:iNumProf,mm080_usual,1:iNumProf,mm120,1:iNumProf,mm080_150levs)
plot(1:iNumProf,mm080_usual-mm120,1:iNumProf,mm080_usual-mm080_150levs)

iN = 100;

semilogy(p080.ptemp(:,iN),p080.plevs(:,iN),p080_150levs.ptemp(:,iN),p080_150levs.plevs(:,iN),'linewidth',2); set(gca,'ydir','reverse'); xlim([170 320]); ylim([0.005 1013])
  ylabel('P [mb]'); xlabel('T(p)'); legend('usual','new','location','best','fontsize',10);
loglog(p080.gas_1(:,iN),p080.plevs(:,iN),p080_150levs.gas_1(:,iN),p080_150levs.plevs(:,iN),'linewidth',2); set(gca,'ydir','reverse')
  ylabel('P [mb]'); xlabel('WV(p)'); legend('usual','new','location','best','fontsize',10);
semilogy(p080.palts(:,iN)/1000,p080.plevs(:,iN),'bx',p080_150levs.palts(:,iN)/1000,p080_150levs.plevs(:,iN),'linewidth',2); set(gca,'ydir','reverse')
  ylabel('P [mb]'); xlabel('H(z)'); legend('usual','new','location','best','fontsize',10);
semilogy(abs(diff(p080.palts(:,iN)/1000)),meanvaluebin(p080.plevs(:,iN)),abs(diff(p080_150levs.palts(:,iN)/1000)),meanvaluebin(p080_150levs.plevs(:,iN)),'linewidth',2); 
  set(gca,'ydir','reverse');
  ylabel('P [mb]'); xlabel('dH(z)'); legend('usual','new','location','best','fontsize',10);

semilogy(p080.ptemp(:,iN),p080.plevs(:,iN),p120.ptemp(:,iN),p120.plevs(:,iN),'linewidth',2); set(gca,'ydir','reverse'); xlim([170 320]); ylim([2e-5 1013])
  ylabel('P [mb]'); xlabel('T(p)'); legend('usual','new','location','best','fontsize',10);
loglog(p080.gas_1(:,iN),p080.plevs(:,iN),p120.gas_1(:,iN),p120.plevs(:,iN),'linewidth',2); set(gca,'ydir','reverse')
  ylabel('P [mb]'); xlabel('WV(p)'); legend('usual','new','location','best','fontsize',10);
semilogy(p080.palts(:,iN)/1000,p080.plevs(:,iN),'bx',p120.palts(:,iN)/1000,p120.plevs(:,iN),'linewidth',2); set(gca,'ydir','reverse')
  ylabel('P [mb]'); xlabel('H(z)'); legend('usual','new','location','best','fontsize',10);
semilogy(abs(diff(p080.palts(:,iN)/1000)),meanvaluebin(p080.plevs(:,iN)),abs(diff(p120.palts(:,iN)/1000)),meanvaluebin(p120.plevs(:,iN)),'linewidth',2); 
  set(gca,'ydir','reverse');
  ylabel('P [mb]'); xlabel('dH(z)'); legend('usual','new','location','best','fontsize',10);

plot(1:iNumProf,p120.palts(1,:)/1000,'b',1:iNumProf,p080.palts(1,:)/1000,'k',...
     1:iNumProf,pcatx.altitude(1,:)/1000,'c',1:iNumProf,pcatx080.altitude(1,:)/1000,'g')
ylim([110 125])
xlabel('Profile'); ylabel('HGT at 2.5e-5 mb TOA [km]')
set(gca,'fontsize',12)

ind = 1 : iNumProf;
ind = 150: 290;
x = [nanmean(p120.palts(1,ind)/1000-pcatx.altitude(1,ind)/1000) nanstd(p120.palts(1,ind)/1000-pcatx.altitude(1,ind)/1000)];
fprintf(1,'mult = %5.3f : 120 km (150 levels) TOA : klayers-simple hyyposometric = %8.4f +/- %8.4f km \n',mult,x);

dh = -2:0.1:+2; plot(dh,histc((p120.palts(1,ind)-pcatx.altitude(1,ind))/1000,dh)); grid
title('TOA : klayers palts - simple hyposometic')
set(gca,'fontsize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iY = input('change to 0:1:120 Z grid for Manuel (-1 default/+1) : ');
if length(iY) == 0
  iY = -1;
end
if iY > 0
  nlte_120kmgrid
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iY = input('compare to msis? (-1 default/+1) : ');
if length(iY) == 0
  iY = -1;
end
if iY > 0
  t_msis
end

iY = input('loop msis? (-1 default/+1) : ');
if length(iY) == 0
  iY = -1;
end
if iY > 0
  loop_t_msis
end
