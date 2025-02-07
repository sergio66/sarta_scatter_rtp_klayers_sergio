klayers = toptsSARTA.klayers;
sarta   = toptsSARTA.sarta_cris;

sartaer   = ['!time ' sarta   ' fin=junkA.op.rtp fout=junkAX.rp.rtp listp=1 listj=300,200,100,1,2,3 > ugh3A'];
eval(sartaer);
[hA,ha,pAX,pa] = rtpread('junkAX.rp.rtp');

fname = 'junkAX.rp.rtp_jacCLD'; [w,dA_cld] = readsarta_jacV2(fname,300);
fname = 'junkAX.rp.rtp_WGTFCN'; [w,dA_wgt] = readsarta_jacV2(fname,200);
fname = 'junkAX.rp.rtp_jacTZ' ; [w,dA_TZZ] = readsarta_jacV2(fname,100);
fname = 'junkAX.rp.rtp_jacG1' ; [w,dA_g1]  = readsarta_jacV2(fname,001);
fname = 'junkAX.rp.rtp_jacG2' ; [w,dA_g2]  = readsarta_jacV2(fname,002);
fname = 'junkAX.rp.rtp_jacG3' ; [w,dA_g3]  = readsarta_jacV2(fname,003);

%%%%%%%%%%%%%%%%%%%%%%%%%

klayers = toptsSARTA.PBL.klayers;
sarta   = toptsSARTA.PBL.sarta_cris;

sartaer   = ['!time ' sarta   ' fin=junkB.op.rtp fout=junkBX.rp.rtp listp=1 listj=300,200,100,1,2,3 > ugh3B'];
eval(sartaer);
[hB,ha,pBX,pa] = rtpread('junkBX.rp.rtp');

fname = 'junkBX.rp.rtp_jacCLD'; [w,dB_cld] = readsarta_jacV2(fname,300);
fname = 'junkBX.rp.rtp_WGTFCN'; [w,dB_wgt] = readsarta_jacV2(fname,200);
fname = 'junkBX.rp.rtp_jacTZ' ; [w,dB_TZZ] = readsarta_jacV2(fname,100);
fname = 'junkBX.rp.rtp_jacG1' ; [w,dB_g1]  = readsarta_jacV2(fname,001);
fname = 'junkBX.rp.rtp_jacG2' ; [w,dB_g2]  = readsarta_jacV2(fname,002);
fname = 'junkBX.rp.rtp_jacG3' ; [w,dB_g3]  = readsarta_jacV2(fname,003);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hA,pAX] = proxy_box_to_ham_fsr_CHepplewhite(hA,pAX,0);
[hB,pBX] = proxy_box_to_ham_fsr_CHepplewhite(hA,pBX,0);

tobsx = rad2bt(hA.vchan,pAX.robs1);
tcalAx = rad2bt(hA.vchan,pAX.rcalc);
tcalBx = rad2bt(hA.vchan,pBX.rcalc);

figure(4); plot(hA.vchan,tobsx-tcalAx,'b',hA.vchan,tobsx-tcalBx,'r')
axis([625 2750 -10 +10])

nA = pAX.nlevs-1;
pA = pAX.plevs(1:nA);
nB = pBX.nlevs-1;
pB = pBX.plevs(1:nB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); colormap jet; pcolor(hA.vchan,pA,dA_wgt(1:nA,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([0 1]*0.1)
ylim([10 1000])

figure(6); colormap jet; pcolor(hB.vchan,pB,dB_wgt(1:nB,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([0 1]*0.1)
ylim([10 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); colormap jet; pcolor(hA.vchan,[pA; pA(end)*1.01],dA_TZZ(1:nA+1,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([0 1]*0.1)
ylim([10 1000])

figure(6); colormap jet; pcolor(hB.vchan,[pB; pB(end)*1.01],dB_TZZ(1:nB+1,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([0 1]*0.1)
ylim([10 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); colormap jet; pcolor(hA.vchan,pA,dA_g1(1:nA,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([-1 0.1]*0.25)
ylim([100 1000])

figure(6); colormap jet; pcolor(hB.vchan,pB,dB_g1(1:nB,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([-1 0.1]*0.25)
ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); colormap jet; pcolor(hA.vchan,pA,dA_g2(1:nA,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([-1 +1]*0.25)
ylim([10 1000])

figure(6); colormap jet; pcolor(hB.vchan,pB,dB_g2(1:nB,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([-1 +1]*0.25)
ylim([10 1000])


%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); colormap jet; pcolor(hA.vchan,pA,dA_g3(1:nA,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([-1 +1]*0.25/3)
ylim([10 1000])

figure(6); colormap jet; pcolor(hB.vchan,pB,dB_g3(1:nB,:)); shading interp; colorbar
  set(gca,'ydir','reverse');   set(gca,'yscale','log')
caxis([-1 +1]*0.25/3)
ylim([10 1000])
