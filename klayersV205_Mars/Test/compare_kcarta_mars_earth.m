addpath /asl/matlib/h4tools
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE

[dMars,w] = readkcstd('/home/sergio/OTHERSTUFF/klayersV205_Mars/Data/mars_kcarta.dat');
[d121,w] = readkcstd('/home/sergio/OTHERSTUFF/klayersV205_Mars/Data/earth_kcarta_121.dat');
[d122,w] = readkcstd('/home/sergio/OTHERSTUFF/klayersV205_Mars/Data/earth_kcarta_122.dat');

disp(' ')
fprintf(1,'comparing Earth kCARTA : v121-v122 = %10.6f \n',sum(d121-d122))
disp(' ')

figure(1); plot(w,rad2bt(w,d121),w,rad2bt(w,dMars))
[fc,qc] = quickconvolve(w,[d122 dMars],0.5,0.5); tc = rad2bt(fc,qc);

figure(1); plot(fc,tc(:,1),fc,tc(:,2))
  title('kCARTA --> 0.5 cm^{-1}'); hl = legend('Earth','Mars','location','best','fontsize',10);

[hm,~,pm,~] = rtpread('mars_sergio_rtp.op.rtp');
[he,~,pe,~] = rtpread('latbin1_40.op.rtp');

figure(2); 
  subplot(121); semilogy(pe.ptemp(1:97,21),pe.plevs(1:97,21),'b',pm.ptemp(1:97,1),pm.plevs(1:97,1),'r')
    set(gca,'ydir','reverse'); ylabel('p(mb)'); xlabel('T(K)')
    hl = legend('Earth','Mars','location','best','fontsize',10);
  subplot(122); loglog(pe.gas_1(1:97,21),pe.plevs(1:97,21),'b',pm.gas_1(1:97,1),pm.plevs(1:97,1),'r')
    set(gca,'ydir','reverse'); xlabel('WV molecules/cm-2')
    hl = legend('Earth','Mars','location','best','fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iSave = input('save data for Guiliano? ');
if length(iSave) == 0
  iSave = -1;
end
if iSave > 0
  save_for_giuliano
end
