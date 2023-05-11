%% checking the SEC(theta)*Q(L)  weighted averages and jacs
%{
[h,ha,p,pa] = rtpread('newdayx_1_100_12150.op.rtp'); p.gas_1 = p.gas_1*1.1; rtpwrite('newdayx_1_100_12150_g1_x1.1.op.rtp',h,ha,p,pa);
[h,ha,p,pa] = rtpread('newdayx_1_100_12150.op.rtp'); p.ptemp = p.ptemp+1;   rtpwrite('newdayx_1_100_12150_Tz_add1.op.rtp',h,ha,p,pa);

date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1614 >& ugh01,ugh0T
date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150_g1_x1.1.op.rtp  fout=newdayx_1_100_12150.rad.rtp listp=1  listj=1  listc=1614 >& ugh1
date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150_Tz_add1.op.rtp  fout=newdayx_1_100_12150.rad.rtp listp=1  listj=1  listc=1614 >& ughT
%} 


if ~exist('p')
  [h,ha,p,pa] = rtpread('newdayx_1_100_12150.op.rtp');
  xavg = p.gas_1/6.023e26;       xavg = xavg(1:97,1);    dx = xavg*0.1;
  Tavg = p.ptemp;                Tavg = Tavg(1:97,1);
  pavg = plevs2plays(p.plevs);   pavg = pavg(1:97,1)/1013;

  %% line 236 in ycalowp.f
  %% write(6,'(A,I3,X,8(E12.4))') 'Y1 001',L,P(L),T(L),WAMNT(L),WAANG(L),WAZ(L),PZ(L),TZ(L),WAZSUM, WAANG_1(L),WAZ_1(L),PZ_1(L),TZ_1(L)
  %%                                       1  2    3     4       5        6      7     8      9        10       11      12       13
  %% write(6,'(A,I3,X,8(E12.4))') 'YT 100',L,P(L),T(L),WAMNT(L),WAANG(L),WAZ(L),PZ(L),TZ(L),WAZSUM, WAANG_T(L),WAZ_T(L),PZ_T(L),TZ_T(L)
  %%                                       1  2    3     4       5        6      7     8      9        10       11      12       13
  ugh0 = load('UGH_TEST_OPTRAN/V0_May3_2023/ugh01');  %% raw WV
  ugh0 = load('ugh01');  %% raw WV
  boo = ugh0(:,1); 
  if mean(boo) ~= 1; 
    error('wxpeting 001 in first column, as this is G1');
  else
    ugh0 = ugh0(:,2:14);
  end

  ugh1 = load('UGH_TEST_OPTRAN/V0_May3_2023/ugh1');  %% WV * 1.1
  ugh1 = load('ugh1');  %% WV * 1.1
  boo = ugh1(:,1); 
  if mean(boo) ~= 1; 
    error('wxpeting 001 in first column, as this is G1');
  else
    ugh1 = ugh1(:,2:14);
  end

end

figure(1); clf; ii = 2; plot(ugh0(:,ii),1:97,'b.-',pavg,1:97,ugh0(:,ii+5),1:97,'k'); set(gca,'ydir','reverse'); title('PAVG')
figure(2); clf; ii = 3; plot(ugh0(:,ii),1:97,'b.-',Tavg,1:97,ugh0(:,ii+5),1:97,'k'); set(gca,'ydir','reverse'); title('TAVG')
figure(3); clf; ii = 4; plot(ugh0(:,ii),1:97,'b.-',xavg,1:97,ugh0(:,ii+1),1:97,'k'); set(gca,'ydir','reverse'); title('QAVG')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('THESE RESULTS SUGGEST PZ_1 and TZ_1 should be 0')
disp('THESE RESULTS SUGGEST PZ_1 and TZ_1 should be 0')
disp('THESE RESULTS SUGGEST PZ_1 and TZ_1 should be 0')
ii = input('Enter index WZ=WAANG WAZ, PZ, TZ  (5:8) : ');
if ii < 5 | ii > 8
  ii = 5;
  disp('reset ii to 5')
end
ugh2 = zeros(size(ugh0)); %% goinna guesstimate this using dy/dx ~ delta y /delta x ... so delta y = dy/dx delta x
ugh2(:,ii) = ugh0(:,ii) + dx.*ugh0(:,ii+5); 

homejac = (ugh1(:,ii)-ugh0(:,ii))./dx;
figure(1); clf; plot(ugh0(:,ii),1:97,ugh1(:,ii),1:97);         set(gca,'ydir','reverse'); title('Raw(b) Perturbed(r) Variable')
figure(2); clf; plot(ugh1(:,ii)-ugh0(:,ii),1:97);              set(gca,'ydir','reverse'); title('Raw-Perturbed')
figure(2); clf; plot(homejac,1:97);                            set(gca,'ydir','reverse'); title('(Raw-Perturbed)/dx')
figure(3); clf; plot(ugh0(:,ii+5),1:97,'b',homejac,1:97,'r');  set(gca,'ydir','reverse'); title('Jacobian')
if ii == 6
  moo = 1./xavg; plot(cumsum(moo),1:97);                 figure(3); clf; plot(ugh0(:,ii+5),1:97,'b',cumsum(moo)/1e11,1:97);               set(gca,'ydir','reverse'); title('Jacobian')
  moo = xavg; plot(cumsum(moo),1:97);                    figure(3); clf; plot(ugh0(:,ii+5),1:97,'b',xavg./cumsum(xavg),1:97);             set(gca,'ydir','reverse'); title('Jacobian')
  boo = 0.5*[1; cumsum(xavg(1:length(xavg)-1))];         figure(3); clf; plot(ugh0(:,ii+5),1:97,'b',homejac,1:97,'r',xavg./boo,1:97,'k'); set(gca,'ydir','reverse'); title('Jacobian')
  boo = 0.5 + [0; cumsum(xavg(1:length(xavg)-1))]./xavg; figure(3); clf; plot(ugh0(:,ii+5),1:97,'b',homejac,1:97,'r',sec(p.satzen(1)*pi/180)*boo,1:97,'k');       set(gca,'ydir','reverse'); title('Jacobian')
end
figure(4); clf; plot(ugh1(:,ii)-ugh0(:,ii),1:97,'bx-',ugh2(:,ii)-ugh0(:,ii),1:97,'r'); set(gca,'ydir','reverse'); 
  title('Did we fix it???? \newline red should lie on blue')
