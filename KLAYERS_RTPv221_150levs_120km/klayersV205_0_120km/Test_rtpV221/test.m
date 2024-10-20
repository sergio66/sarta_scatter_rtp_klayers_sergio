%{
[h,ha,p,pa] = rtpread('/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_120km/klayersV205_0_120km/Test_rtpV221/pin_feb2002_raw_ip.rtp');
for ii = 1 : 49
  nlevs = p.nlevs(ii);
  p.spres(ii) = p.plevs(nlevs,ii);
  p.spres(ii) = p.plevs(1,ii);
end
rtpwrite('/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_120km/klayersV205_0_120km/Test_rtpV221/pin_feb2002_raw_ip.rtp',h,ha,p,pa);
%}

[h120,ha,p120,pa] = rtpread('/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_120km/klayersV205_0_120km/Test_rtpV221/pin_feb2002_raw_op_0_120km.rtp');
[h080,ha,p080,pa] = rtpread('/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_RTPv221_150levs_120km/klayersV205_0_120km/Test_rtpV221/pin_feb2002_raw_op_0_080km.rtp');

mm080 = mmwater_rtp(h080,p080);
mm120 = mmwater_rtp(h120,p120);

plot(1:49,mm080,1:49,mm120)
plot(1:49,mm080-mm120)

semilogy(p080.ptemp(:,1),p080.plevs(:,1),p120.ptemp(:,1),p120.plevs(:,1),'linewidth',2); set(gca,'ydir','reverse')
loglog(p080.gas_1(:,1),p080.plevs(:,1),p120.gas_1(:,1),p120.plevs(:,1),'linewidth',2); set(gca,'ydir','reverse')
semilogy(p080.palts(:,1)/1000,p080.plevs(:,1),'bx',p120.palts(:,1)/1000,p120.plevs(:,1),'linewidth',2); set(gca,'ydir','reverse')

plot(p120.palts(1,:)/1000)
