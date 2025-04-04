addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/git/CRODGERS_FAST_CLOUD
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p,pa] = rtpread('mktemp_4PjEkGC2_cris_hr_pbl_op.rtp');
fprintf(1,'there are %6i profiles \n',length(p.stemp))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ../../bin/jac_crisg4_hires_feb25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new_PBL    fin=mktemp_4PjEkGC2_cris_hr_pbl_op.rtp fout=testall.rp.rtp

% ../../bin/jac_crisg4_hires_jan25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_newFeb2025 fin=mktemp_4PjEkGC2_cris_hr_pbl_op.rtp fout=test.rp.rtp listp=1,6,1749,5274,7128,7740,8609,9276,12313,16712,18118,18127,19051  OOOOOPS
% ../../bin/jac_crisg4_hires_feb25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new_PBL    fin=mktemp_4PjEkGC2_cris_hr_pbl_op.rtp fout=test.rp.rtp listp=1,6,1749,5274,7128,7740,8609,9276,12313,16712,18118,18127,19051 

% ../../bin/jac_crisg4_hires_feb25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new_PBL    fin=mktemp_4PjEkGC2_cris_hr_pbl_op.rtp fout=test1749.rp.rtp listp=1749
% ../../bin/jac_crisg4_hires_feb25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new_PBL    fin=mktemp_4PjEkGC2_cris_hr_pbl_op.rtp fout=test1749_900cm.rp.rtp listp=1749 listc=403
thelist = [1 6 1749 5274 7128 7740 8609 9276 12313 16712 18118 18127 19051];

[hall,pa,pall,pa] = rtpread('testall.rp.rtp');
  [h1,p1] = proxy_box_to_ham_fsr(hall,pall,2223);
  [h2,p2] = proxy_box_to_ham_fsr_CHepplewhite(hall,pall,0);
  pall = p2;

i900 = find(hall.vchan >= 900,1);
tobsall = rad2bt(hall.vchan,pall.robs1);
tcalall = rad2bt(hall.vchan,pall.rcalc);

plot(pall.stemp,rad2bt(900,pall.robs1(i900,:)),'.',pall.stemp,rad2bt(900,pall.rcalc(i900,:)),'.')
xlabel('SKT'); ylabel('BT 900'); legend('obs','cal')
plot(pall.stemp,rad2bt(900,pall.robs1(i900,:)) - rad2bt(900,pall.rcalc(i900,:)),'.')
xlabel('SKT'); ylabel('BT 900'); legend('obs-cal')

plot(hall.vchan,mean(tobsall'-tcalall'),hall.vchan,std(tobsall'-tcalall'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hx,ha,px,pa] = rtpread('test.rp.rtp');
px.zobs = 805000 * ones(size(px.stemp));
px.scanang = saconv(px.satzen,px.zobs);
rtpwrite('test.rp.rtp',hx,ha,px,pa)


iProf = 3;
[dkc,fkc] = readkcstd(['profile' num2str(iProf) '.dat']);
[fc,qc] = quickconvolve(fkc,dkc,1,1);
tc = rad2bt(fc,qc);
[mm,nn] = size(tc);

[p.stemp(thelist(iProf)) px.stemp(iProf)]


[h1749,~,p1749,~] = rtpread('test1749.rp.rtp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see /home/sergio/git/CRODGERS_FAST_CLOUD/driver_stagel_getparam_setoutput.m
  disp('need to change CrIS FSR obs from sinc to hamming (SARTA calc are hamming)')
  tobsBeforeHamm = rad2bt(h.vchan,p.robs1);
  [h1,p1] = proxy_box_to_ham_fsr(hx,px,2223);
  [h2,p2] = proxy_box_to_ham_fsr_CHepplewhite(hx,px,0);
  px = p2;
  clear p1 p2 h2 h1
  junk = find(h.ichan <= 2223-12);

tobs = rad2bt(hx.vchan,px.robs1);
tcal = rad2bt(hx.vchan,px.rcalc);
plot(1:length(px.stemp),tobs(i900,:),1:length(px.stemp),tcal(i900,:))

plot(hx.vchan,tobs(:,iProf),hx.vchan,tcal(:,iProf),fc,tc(:,nn),h1749.vchan,rad2bt(h1749.vchan,p1749.rcalc),'k'); 
legend('Obs','SARTA','KCARTA','location','best'); ylabel('BT'); xlabel('Wavenumber cm-1'); set(gca,'fontsize',12)

plot(hx.vchan,nanmean(tobs,2))
plot(hx.vchan,nanmean(tobs-tcal,2))

plevs101 = flipud(load('/home/sergio/MATLABCODE/airslevels.dat'));
[px.plevs(:,1) plevs101]
[px.plevs(:,1) ./ plevs101]

