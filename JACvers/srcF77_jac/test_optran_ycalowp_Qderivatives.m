%% checking the SEC(theta)*Q(L)  weighted averages and jacs
if ~exist('p')
  [h,ha,p,pa] = rtpread('newdayx_1_100_12150.op.rtp');
  xavg = p.gas_1/6.023e26;       xavg = xavg(1:97,1);    dx = xavg*0.1;
  Tavg = p.ptemp;                Tavg = Tavg(1:97,1);
  pavg = plevs2plays(p.plevs);   pavg = pavg(1:97,1)/1013;

  %% line 236 in ycalowp.f
  %% write(6,'(A,I3,X,8(E12.4))') 'Y1',L,P(L),T(L),WAMNT(L),WAANG(L),WAZ(L),PZ(L),TZ(L),WAZSUM, WAANG_1(L),WAZ_1(L),PZ_1(L),TZ_1(L)
  %%                                   1  2    3     4       5        6      7     8      9        10       11      12       13
  ugh0 = load('ugh0');  %% raw WV
  ugh1 = load('ugh1');  %% WV * 1.1
end

figure(1); clf; ii = 2; plot(ugh0(:,ii),1:97,'b.-',pavg,1:97,ugh0(:,ii+5),1:97,'k'); set(gca,'ydir','reverse'); title('PAVG')
figure(2); clf; ii = 3; plot(ugh0(:,ii),1:97,'b.-',Tavg,1:97,ugh0(:,ii+5),1:97,'k'); set(gca,'ydir','reverse'); title('TAVG')
figure(3); clf; ii = 4; plot(ugh0(:,ii),1:97,'b.-',xavg,1:97,ugh0(:,ii+1),1:97,'k'); set(gca,'ydir','reverse'); title('QAVG')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('THESE RESULTS SUGGEST PZ_1 and TZ_1 should be 0')
disp('THESE RESULTS SUGGEST PZ_1 and TZ_1 should be 0')
disp('THESE RESULTS SUGGEST PZ_1 and TZ_1 should be 0')
ii = input('Enter index WZ=WAANG WAZ, PZ, TZ  (5:8) : ');
ugh2 = zeros(size(ugh0)); %% goinna guesstimate this using dy/dx ~ delta y /delta x ... so delta y = dy/dx delta x
ugh2(:,ii) = ugh0(:,ii) + dx.*ugh0(:,ii+5); 

figure(1); clf; plot(ugh0(:,ii),1:97,ugh1(:,ii),1:97); set(gca,'ydir','reverse'); title('Raw(b) Perturbed(r) Variable')
figure(2); clf; plot(ugh1(:,ii)-ugh0(:,ii),1:97); set(gca,'ydir','reverse'); title('Raw-Perturbed')
figure(3); clf; plot(ugh0(:,ii+5),1:97); set(gca,'ydir','reverse'); title('Jacobian')
figure(4); clf; plot(ugh1(:,ii)-ugh0(:,ii),1:97,'bx-',ugh2(:,ii)-ugh0(:,ii),1:97,'r'); set(gca,'ydir','reverse'); 
  title('Did we fix it???? \newline red should lie on blue')
