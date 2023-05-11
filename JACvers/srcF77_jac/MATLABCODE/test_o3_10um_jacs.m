%{
% see lines 459 and 699 in ycalt2_od.f
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_nocldfields.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=3  listc=1092 >& ugh0_ozone
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_nocldfields_pertO3x1.1.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=3  listc=1092 >& ughpert_ozone
%}

[h,ha,p0,pa] = rtpread('newdayx_nocldfields.op.rtp');
[h,ha,p3,pa] = rtpread('newdayx_nocldfields_pertO3x1.1.op.rtp');
nlays = p0.nlevs(1)-1;
x = p0.gas_3(:,1); x = x(1:nlays); dx = x(1:nlays)/6.02214199E+26 * 0.1;   %% see rdrtp.f
t = p0.ptemp(:,1); t = t(1:nlays); st = p0.stemp(1); mu = cos(p0.satzen(1)*pi/180);

junk = load('ugh0_ozone');    len = length(junk); rawrad  = junk(len,3); junk = junk(1:len-1,3); rad_ODjac  = reshape(junk,nlays,2); %% OZ only, but has correct jas in last 100
junk = load('ughpert_ozone'); len = length(junk); pertrad = junk(len,3); junk = junk(1:len-1,3); pert_ODjac = reshape(junk,nlays,2); %% OZ only, but has correct jac in last 100

figure(1); clf;  plot(rad_ODjac(:,1),1:nlays,pert_ODjac(:,1),1:nlays); title('ODs raw(b) and pert (r)'); set(gca,'ydir','reverse');
figure(2); clf; jac_here = (pert_ODjac(:,1)-rad_ODjac(:,1))./dx; plot(jac_here,1:nlays,rad_ODjac(:,2),1:nlays); set(gca,'ydir','reverse'); hl = legend('finite difference','from SARTA analytic jac','location','best','fontsize',10);

%junk = load('ugh0_ozone');    len = length(junk); rawrad  = junk(len,4); junk = junk(1:len-1,3); rad_ODjac  = reshape(junk,nlays,2);  %% OZ + other gases is more correct, but has 0 in last 100
%junk = load('ughpert_ozone'); len = length(junk); pertrad = junk(len,4); junk = junk(1:len-1,3); pert_ODjac = reshape(junk,nlays,2);  %% OZ + other gases is more correct, but has 0 in last 100

%{
miaow = load('/home/sergio/KCARTA/L2SComparisons/l2s_kc122_H16_605_2830.mat');
[fc,qc] = convolve_airs(miaow.w,miaow.d(:,3),1:2378); figure(3); plot(fc,qc); 
booc = find(fc >= 1041,1); [qc(booc) sum(rad_ODjac(:,1))]
%}

freq = 1042; nemis = p0.nemis(1); emis = interp1(p0.efreq(1:nemis,1),p0.emis(1:nemis,1),freq);

r0 = emis * ttorad(freq,st);
r1 = emis * ttorad(freq,st);
for ii = 97 : -1 : 1
  fprintf(1,'[ii t(ii)] = %3i %8.6f \n',[ii t(ii)])
  R = ttorad(freq,t(ii));
  q0 = rad_ODjac(ii,1);                         r0 = r0 *exp(-q0) + R * (1-exp(-q0));
  q1 = pert_ODjac(ii,1);                        r1 = r1 *exp(-q1) + R * (1-exp(-q1));
end  
[rad2bt(1040,[rawrad*1000 r0 pertrad*1000 r1]) rad2bt(1040,pertrad*1000)-rad2bt(1040,rawrad*1000) rad2bt(1040,r1)-rad2bt(1040,r0)]
r0save = r0;
r1save = r1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now do jac calcs
rpert = emis * ttorad(freq,st);
rpert = rpert * ones(1,97);
rest  = rpert;
for jacloop = 1 : 97
  r0 = rpert(jacloop);
  r1 = rpert(jacloop);
  for ii = 97 : -1 : 1
    R = ttorad(freq,t(ii));
    q0 = rad_ODjac(ii,1);
    q1 = rad_ODjac(ii,1);
    if jacloop == ii
      q0 = pert_ODjac(ii,1);
      q1 = rad_ODjac(ii,1) + dx(ii)*jac_here(ii);     %% makes things PERFECT
      q1 = rad_ODjac(ii,1) + dx(ii)*rad_ODjac(ii,2);  %% SARTA analytic jacs TOO BIG
    end
    r0 = r0 *exp(-q0) + R * (1-exp(-q0));
    r1 = r1 *exp(-q1) + R * (1-exp(-q1));
  end
  rpert(jacloop) = r0;
  rest(jacloop)  = r1;
end  
figure(3); plot(rpert,1:97,rest,1:97); set(gca,'ydir','reverse'); title('Comparing perturbed calcs'); hl = legend('from ODs','from Jac OD x dx','fontsize',10,'location','best');

jac   = (rpert - r0save)./dx';  
jace  = (rest - r0save)./dx';  
figure(4); plot(jac,1:97,jace,1:97); set(gca,'ydir','reverse'); title('Jac in RAD space : jace too large!!!'); hl = legend('from ODs','from Jac OD x dx','fontsize',10,'location','best');

jacBT  = rad2bt(freq,rpert) - rad2bt(freq,r0save); 
jacBTe = rad2bt(freq,rest) - rad2bt(freq,r0save); 
%% recall log(1+d) = log(1) + d 1/1 + d2(-1/1**2) etc ~ d
%% thechnically now we do delta BT      =        deltaBT              = delta BT        delta BT             delta(BT)  * q      delta BT  * X     delta BT    * X    delta BT
%%                        ----------          --------------------      --------    ~~  --------     ======  ---------       =   --------       =  --------        = ---------
%%                        delta log(q)         log(X(1+d))-log(X)        log(1+d)         d                    delta(q)          X+dx - X          X(1+0.1)-X            0.1
  jacBT  = jacBT/log(1.1);  
  jacBTe = jacBTe/log(1.1);  
figure(5); clf; 
plot(jacBT,1:97,jacBTe,1:97); set(gca,'ydir','reverse'); title('Jac in BT space : jace too large!!!'); hl = legend('from ODs','from Jac OD x dx','fontsize',10,'location','best');
fprintf(1,'COLUMN SUM JAC BT = %8.6f %8.6f\n',sum(jacBT),sum(jacBTe))

%{
save o3_jac_column_sum.mat jacBT
%}
