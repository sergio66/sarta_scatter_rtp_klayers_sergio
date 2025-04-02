addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/git/CRODGERS_FAST_CLOUD
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE

[h,ha,p,pa] = rtpread('mktemp_4PjEkGC2_cris_hr_pbl_op.rtp');
fprintf(1,'there are %6i profiles \n',length(p.stemp))

[hx,pa,px,pa] = rtpread('test.rp.rtp');
px.zobs = 805000 * ones(size(px.stemp));
px.scanang = saconv(px.satzen,px.zobs);
rtpwrite('test.rp.rtp',hx,ha,px,pa)

iProf = 3;
[dkc,fkc] = readkcstd(['profile' num2str(iProf) '.dat']);
[fc,qc] = quickconvolve(fkc,dkc,1,1);
tc = rad2bt(fc,qc);
[mm,nn] = size(tc);


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
i900 = find(hx.vchan >= 900,1);
plot(1:length(px.stemp),tobs(i900,:),1:length(px.stemp),tcal(i900,:))

plot(hx.vchan,tobs(:,iProf),hx.vchan,tcal(:,iProf),fc,tc(:,nn)); 
legend('Obs','SARTA','KCARTA','location','best'); ylabel('BT'); xlabel('Wavenumber cm-1'); set(gca,'fontsize',12)

plot(hx.vchan,nanmean(tobs,2))
plot(hx.vchan,nanmean(tobs-tcal,2))

plevs101 = flipud(load('/home/sergio/MATLABCODE/airslevels.dat'));
[px.plevs(:,1) plevs101]
[px.plevs(:,1) ./ plevs101]

