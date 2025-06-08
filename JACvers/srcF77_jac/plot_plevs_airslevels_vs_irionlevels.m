a = load('/home/sergio/MATLABCODE/airslevels.dat');
b = load('bill_pbl_levs.txt');

figure(1); clf
hold on
for ii = 1 : length(a)
  line([1 2],[a(ii) a(ii)],'color','b');
  line([3 4],[b(ii) b(ii)],'color','r');
end
set(gca,'ydir','reverse');
ylim([0 1100])
set(gca,'yscale','log')
ylabel('level P [mb]')
ylim([600 1100])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE
%% p2h is for US STD

ah = p2h(a);            ah = diff(ah);
alays = plevs2plays(a); alays = alays(1:end-1);

bh = p2h(b);            bh = -diff(bh);
blays = plevs2plays(b); blays = blays(1:end-1);

figure(2); clf
plot(ah,alays,'bx-',bh,blays,'ro-');
set(gca,'ydir','reverse');
ylim([0 1100])
set(gca,'yscale','log')

ylabel('layer avg P [mb]')
xlabel('Thickness [m]')
ylim([600 1100])
