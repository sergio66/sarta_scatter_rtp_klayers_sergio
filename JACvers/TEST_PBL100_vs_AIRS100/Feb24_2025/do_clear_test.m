[h,ha,p,pa] = rtpread(fip);
p.cngwat  = 0 * p.cngwat;
p.cngwat2 = 0 * p.cngwat2;
p.cfrac  = 0 * p.cfrac;
p.cfrac2 = 0 * p.cfrac2;
p.cfrac12 = 0 * p.cfrac12;
rtpwrite('clear_airs_l1c_ecm_sarta_baum_ice.2025.01.08.099.rtp',h,ha,p,pa);
fip = 'clear_airs_l1c_ecm_sarta_baum_ice.2025.01.08.099.rtp';

klayers = toptsSARTA.klayers;
sarta   = toptsSARTA.sarta_cris;
klayerser = ['!time ' klayers ' fin=' fip '      fout=junkA.op.rtp > ugh1A'];
sartaer   = ['!time ' sarta   ' fin=junkA.op.rtp fout=junkA.rp.rtp > ugh2A'];
eval(klayerser);
eval(sartaer);
[hA,ha,pA,pa] = rtpread('junkA.rp.rtp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

klayers = toptsSARTA.PBL.klayers;
sarta   = toptsSARTA.PBL.sarta_cris;
klayerser = ['!time ' klayers ' fin=' fip '      fout=junkB.op.rtp > ugh1B'];
sartaer   = ['!time ' sarta   ' fin=junkB.op.rtp fout=junkB.rp.rtp > ugh2B'];
eval(klayerser);
eval(sartaer);
[hB,ha,pB,pa] = rtpread('junkB.rp.rtp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hA,pA] = proxy_box_to_ham_fsr_CHepplewhite(hA,pA,0);
[hB,pB] = proxy_box_to_ham_fsr_CHepplewhite(hA,pB,0);

compare_two_structures(pA,pB)

tobs   = rad2bt(hA.vchan,pA.robs1);
tcalcA = rad2bt(hA.vchan,pA.rcalc);
tcalcB = rad2bt(hB.vchan,pB.rcalc);

figure(1); clf; 
plot(hA.vchan,nanmean(tobs'-tcalcA'),'b.-',hA.vchan,nanmean(tobs'-tcalcB'),'r',hA.vchan,nanstd(tobs'-tcalcA'),'c.-',hA.vchan,nanstd(tobs'-tcalcB'),'m')
plotaxis2;
legend('mean obs-100A','mean obs-PBL','std obs-100A','std obs-PBL','location','best','fontsize',10);
axis([625 2750 -20 +20])

figure(2); clf; 
plot(hA.vchan,nanmean(tcalcB'-tcalcA'),'b',hA.vchan,nanstd(tcalcB'-tcalcA'),'c')
plotaxis2;
legend('mean PBL-100A','std PBL-100A','location','best','fontsize',10);
title('this is CLEAR SARTA CrIS : PBL 100 - L2 100')
axis([625 2750 -10 +10])
figure(2); ylim([-1 +1])
