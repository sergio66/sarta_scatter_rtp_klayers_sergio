%% checking the SEC(theta)*Q(L)  weighted averages and jacs
%{
[h,ha,p,pa] = rtpread('newdayx_1_100_12150.op.rtp'); p.gas_1 = p.gas_1*1.1; rtpwrite('newdayx_1_100_12150_g1_x1.1.op.rtp',h,ha,p,pa);
[h,ha,p,pa] = rtpread('newdayx_1_100_12150.op.rtp'); p.ptemp = p.ptemp+1;   rtpwrite('newdayx_1_100_12150_Tz_add1.op.rtp',h,ha,p,pa);

date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1614 >& ughLOP_10
date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150_g1_x1.1.op.rtp  fout=newdayx_1_100_12150.rad.rtp listp=1  listj=1  listc=1614 >& ughLOP_1

date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1603 >& ughLOP_10
date; time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150_g1_x1.1.op.rtp  fout=newdayx_1_100_12150.rad.rtp listp=1  listj=1  listc=1603 >& ughLOP_1

%} 


if ~exist('p')
  [h,ha,p,pa] = rtpread('newdayx_1_100_12150.op.rtp');
  xavg = p.gas_1/6.023e26;       xavg = xavg(1:97,1);    dx = 0*xavg*0.1;
  Tavg = p.ptemp;                Tavg = Tavg(1:97,1);    dx = ones(size(dx));
  pavg = plevs2plays(p.plevs);   pavg = pavg(1:97,1)/1013;

  %% line 236 in ycalowp.f
  %%                                           1    2     3        4 5   6     7    8    9    10     11    12  13    14      15
  %%                                                              [1 2   3    4     5    6]
  %%       write(6,'(A,I4,X,12(E12.4))') 'LOP 001',LOP,WAZOP(LOP),DA,POP,TOP,PZOP,TZOP,ANGOP,DA_1,POP_1,TOP_1,PZOP_1,TZOP_1,ANGOP_1
  %%       write(6,'(A,I4,X,12(E12.4))') 'LOP 100',LOP,WAZOP(LOP),DA,POP,TOP,PZOP,TZOP,ANGOP,DA_T,POP_T,TOP_T,PZOP_T,TZOP_T,ANGOP_T
  ugh0 = load('UGH_TEST_OPTRAN/V0_May3_2023/ughLOP_10');  %% raw WV
  ugh0 = load('UGH_TEST_OPTRAN/V2_May5_2023/ughLOP_10');  %% raw WV
  ugh0 = load('ughLOP_10');  %% raw WV
  boo = ugh0(:,1); 
  if mean(boo) ~= 1; 
    error('wxpeting 001 in first column, as this is G1');
  else
    ugh0 = ugh0(:,2:15);
  end

  ugh1 = load('UGH_TEST_OPTRAN/V0_May3_2023/ughLOP_1');  %% WV * 1.1
  ugh1 = load('UGH_TEST_OPTRAN/V2_May5_2023/ughLOP_1');  %% WV * 1.1
  ugh1 = load('ughLOP_1');  %% WV * 1.1
  boo = ugh1(:,1); 
  if mean(boo) ~= 001; 
    error('wxpeting 001 in first column, as this is G1');
  else
    ugh1 = ugh1(:,2:15);
  end

end

[len,~] = size(ugh0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('THESE RESULTS SUGGEST PZ_T and TZ_T should be 0,1')
disp('THESE RESULTS SUGGEST PZ_T and TZ_T should be 0,1')
disp('THESE RESULTS SUGGEST PZ_T and TZ_T should be 0,1')
ii = input('Enter index DA,POP,TOP,PZOP,TZOP,ANGOP (3:8) : ');
if ii < 3 | ii > 8
  ii = 3;
  disp('reset ii to 3')
end
ugh2 = zeros(size(ugh0)); %% goinn guesstimate this using dy/dx ~ delta y /delta x ... so delta y = dy/dx delta x
dx = ones(len,1);
dx = ugh0(:,2)*0.1;
ugh2(:,ii) = ugh0(:,ii) + dx.*ugh0(:,ii+6); 

homejac = (ugh1(:,ii)-ugh0(:,ii))./dx;
figure(1); clf; plot(ugh0(:,ii),1:len,ugh1(:,ii),1:len);         set(gca,'ydir','reverse'); title('Raw(b) Perturbed(r) Variable')
figure(2); clf; plot(ugh1(:,ii)-ugh0(:,ii),1:len);              set(gca,'ydir','reverse'); title('Raw-Perturbed')
figure(2); clf; plot(homejac,1:len);                            set(gca,'ydir','reverse'); title('(Raw-Perturbed)/dx')
figure(3); clf; plot(ugh0(:,ii+6),1:len,'b.-',homejac,1:len,'r');  set(gca,'ydir','reverse'); title('Jacobian')
%if ii == 6
%  moo = 1./xavg; plot(cumsum(moo),1:len);
%  figure(3); clf; plot(ugh0(:,ii+6),1:len,cumsum(moo)/1e11,1:len); set(gca,'ydir','reverse'); title('Jacobian')
%end
figure(4); clf; plot(ugh1(:,ii)-ugh0(:,ii),1:len,'bx-',ugh2(:,ii)-ugh0(:,ii),1:len,'r'); set(gca,'ydir','reverse'); 
  title('Did we fix it???? \newline red should lie on blue')
