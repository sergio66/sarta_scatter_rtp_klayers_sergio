boo = load('/home/sergio/IR_NIR_VIS_UV_RTcodes/NRL_MSIS2.1/msis2.1_test_out.txt');

glat = boo(1,4);
glon = boo(1,5);
woo = (glat-pcatx.rlat).^2 + (glon-pcatx.rlon).^2;
ic = find(woo == min(woo));

toff = boo(5,20) - pnlte.ptemp(5,ic);
plot(boo(:,20),boo(:,3),pnlte.ptemp(:,ic)+toff,pnlte.palts(:,ic)/1000,pcatx.ptemp(:,ic)+toff,pcatx.palts(:,ic)/1000); ylim([0 120])

