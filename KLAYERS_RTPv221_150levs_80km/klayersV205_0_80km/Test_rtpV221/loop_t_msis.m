addpath /home/sergio/MATLABCODE/TIME
[yy,mm,dd,hh] = tai2utcSergio(pcatx.rtime);
doy = change2days(yy,mm,dd,2002);
plot(doy)
doy = change2days(2002*ones(size(yy)),mm,dd,2002);

pnlte2 = pnlte;

for ii = 50 : length(pnlte.stemp)
  if mod(ii,100) == 0
    fprintf(1,'+')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end

  sedder = ['!sed -e "s/RLAT/' num2str(pnlte.rlat(ii)) '/g"  -e "s/RLON/' num2str(pnlte.rlon(ii)) '/g" -e "s/DOY/' num2str(doy(ii)) '/g" '];
  sedder = [sedder ' /home/sergio/IR_NIR_VIS_UV_RTcodes/NRL_MSIS2.1/msis2.1_test_in_sed.txt > /home/sergio/IR_NIR_VIS_UV_RTcodes/NRL_MSIS2.1/msis2.1_test_in_x.txt'];
  eval(sedder)

  rmer = ['!/bin/rm /home/sergio/IR_NIR_VIS_UV_RTcodes/NRL_MSIS2.1/msis2.1_test_in.txt']; eval(rmer)
  lner = ['!ln -s /home/sergio/IR_NIR_VIS_UV_RTcodes/NRL_MSIS2.1/msis2.1_test_in_x.txt  /home/sergio/IR_NIR_VIS_UV_RTcodes/NRL_MSIS2.1/msis2.1_test_in.txt']; eval(lner)
  
  cd  /home/sergio/IR_NIR_VIS_UV_RTcodes/NRL_MSIS2.1/
  runner = ['!msis2.1_test.exe'];
  eval(runner)
  cd /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_120km/klayersV205_0_120km/Test_rtpV221/

  boo = load('/home/sergio/IR_NIR_VIS_UV_RTcodes/NRL_MSIS2.1/msis2.1_test_out.txt');

  glat = boo(1,4);
  glon = boo(1,5);
  woo = (glat-pcatx.rlat).^2 + (glon-pcatx.rlon).^2;
  ic = find(woo == min(woo));

  toff(ii) = boo(5,20) - pnlte.ptemp(5,ic);
  plot(boo(:,20),boo(:,3),pnlte.ptemp(:,ic)+toff(ii),pnlte.palts(:,ic)/1000,pcatx.ptemp(:,ic)+toff(ii),pcatx.palts(:,ic)/1000); ylim([0 120]); pause(0.1)
  pnlte2.ptemp(:,ii) = boo(1:110,20);
end
fprintf(1,'\n')
fprintf(1,'toffset to match at 5 km = %8.6f +/- %8.6f K \n',mean(toff),std(toff))

plot(pnlte2.ptemp - pnlte.ptemp,pnlte.palts/1000,'c',...
     nanmean(pnlte2.ptemp - pnlte.ptemp,2),pnlte.palts/1000,'b',nanstd(pnlte2.ptemp - pnlte.ptemp,[],2),pnlte.palts/1000,'b--')
plotaxis2;
xlabel('\delta T = MSIS - KLAYERS [K]'); ylabel('H [km]')

addpath /home/sergio/MATLABCODE/PLOTMISC
plot(nanmean(pnlte2.ptemp,2),pnlte.palts/1000,'r',nanmean(pnlte.ptemp,2),pnlte.palts/1000,'b','linewidth',2)
hold on
shadedErrorBarY(nanmean(pnlte2.ptemp,2),nanmean(pnlte2.palts,2)/1000,nanstd(pnlte2.ptemp,[],2),'r',0.25);
shadedErrorBarY(nanmean(pnlte.ptemp,2), nanmean(pnlte.palts,2)/1000, nanstd(pnlte.ptemp,[],2),'b',0.25);
plotaxis2;
xlabel('(r) MSIS  (b) KLAYERS [K]'); ylabel('H [km]')
hold off

i60 = find(pnlte.palts(:,1)/1000 < 20);
i60 = find(pnlte.palts(:,1)/1000 > 60);
i60 = find(pnlte.palts(:,1)/1000 > 90);
i60 = find(pnlte.palts(:,1)/1000 > 40);
Lklayers = driver_eof1(pnlte.ptemp(i60,:));
Lmsis  = driver_eof1(pnlte2.ptemp(i60,:));
LLklayers = Lklayers/sum(Lklayers);
LLmsis  = Lmsis/sum(Lmsis);
plot(1:length(Lklayers),cumsum(LLklayers),1:length(Lklayers),cumsum(LLmsis))
plot(1:length(Lklayers),cumsum(LLklayers),1:length(Lklayers),cumsum(LLmsis),'linewidth',2)
axis([0 10 0 1])
axis([0 10 0.8 1])
title('(r) MSIS  (b) KLAYERS')

%{
pnlte_msis = pnlte2;
comment = 'see /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_120km/klayersV205_0_120km/Test_rtpV221/test_300_49.m';
comment = 'see /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_120km/klayersV205_0_120km/Test_rtpV221/driver_test_300_49_make_nlte_Z.m';
save pZ_0_120km.mat hnlte pnlte pnlte_msis comment
%}
