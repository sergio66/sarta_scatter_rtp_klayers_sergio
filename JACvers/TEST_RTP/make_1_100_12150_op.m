addpath /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/matlibSergio/matlib2025/h4tools
addpath /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/matlibSergio/matlib2025/rtptools
addpath /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/matlibSergio/matlib2025/aslutil
addpath /home/sergio/git/sergio_matlib/matlib/clouds/sarta/
addpath /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode/PLOTTER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fin = 'your_favorite_rtp  SO you can easily re-generate if needed';
fin = '/home/sergio/nogit/sergio_temp_rtp_files/rtp_airicrad_v6/2018/06/29/cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.186_cumsum_-1.rtp';

iSkip = 10;
iSkip = 100;  %% for saving small rtp files if needed)

[h0,ha,p0,pa] = rtpread(fin);
[h,p] = subset_rtp_allcloudfields(h0,p0,[],[],1:iSkip:12150);

fip = 'newdayx_1_100_12150.ip.rtp';
fop = 'newdayx_1_100_12150.op.rtp';

rtpwrite(fip,h,ha,p,pa);

klayers = '/home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_80km/klayersV205_0_80km/BinV221/klayers_airs';

klayerser = ['!time ' klayers ' fin=' fip ' fout=' fop];
eval(klayerser)
disp('did klayers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sarta = '../bin/jac_airs_l1c_2834_cloudy_jan25_H2020';
frp20 = 'newdayx_1_100_12150.H2020.rp.rtp';

sartaer = ['!time ' sarta ' fin=' fop ' fout=' frp20];
eval(sartaer)
disp('did sarta H2020')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sarta = '../bin/jac_airs_l1c_2834_cloudy_apr26_H2024';
frp24 = 'newdayx_1_100_12150.H2024.rp.rtp';

sartaer = ['!time ' sarta ' fin=' fop ' fout=' frp24];
eval(sartaer)
disp('did sarta H2024')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hh,ha,p20,pa] = rtpread(frp20); tobs = rad2bt(hh.vchan,p20.robs1);  tcal20 = rad2bt(hh.vchan,p20.rcalc);
[hh,ha,p24,pa] = rtpread(frp24); tobs = rad2bt(hh.vchan,p24.robs1);  tcal24 = rad2bt(hh.vchan,p24.rcalc);

[Y,I] = sort(hh.vchan);
figure(1); clf; plot(hh.vchan(I),nanmean(tobs(I,:)'),'k',hh.vchan(I),nanmean(tcal20(I,:)'),'b',hh.vchan(I),nanmean(tcal24(I,:)'),'r');
  legend('obs','H2020','H2024','location','best');
figure(2); clf; plot(hh.vchan(I),nanmean(tobs(I,:)'-tcal20(I,:)'),'b',hh.vchan(I),nanmean(tobs(I,:)'-tcal24(I,:)'),'r');    legend('bias H2020','bias H2024','location','best');
figure(3); clf; plot(hh.vchan(I),nanstd(tobs(I,:)'-tcal20(I,:)'),'b',hh.vchan(I),nanstd(tobs(I,:)'-tcal24(I,:)'),'r');      legend('std H2020','std H2024','location','best');
figure(4); clf; plot(hh.vchan(I),nanmean(tcal20(I,:)'-tcal24(I,:)'),'b',hh.vchan(I),nanstd(tcal20(I,:)'-tcal24(I,:)'),'r'); legend('bias H2020-H2024','std H2020-H2024','location','best');
figure(5); clf; scatter_coast(p20.rlon,p20.rlat,10,tobs(1520,:)); title('BT 1231 Obs'); colormap jet
