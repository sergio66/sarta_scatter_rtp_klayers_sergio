%{
[h,ha,p0,pa] = rtpread('newdayx_1_100_12150.op.rtp');
dQ = 0.1;
dQ = 0.01;
p = p0; p.cngwat = p.cngwat*(1+dQ);   rtpwrite('TEST_RTP/newdayx_1_100_12150_pert_cngwat1.op.rtp',h,ha,p,pa);
p = p0; p.cpsize = p.cpsize*(1+dQ);   rtpwrite('TEST_RTP/newdayx_1_100_12150_pert_cpsize1.op.rtp',h,ha,p,pa);
p = p0; p.cprtop = p.cprtop*(1+dQ);   rtpwrite('TEST_RTP/newdayx_1_100_12150_pert_cprtop1.op.rtp',h,ha,p,pa);
p = p0; p.cprbot = p.cprbot*(1+dQ);   rtpwrite('TEST_RTP/newdayx_1_100_12150_pert_cprbot1.op.rtp',h,ha,p,pa);

p = p0; p.cngwat2 = p.cngwat2*(1+dQ); rtpwrite('TEST_RTP/newdayx_1_100_12150_pert_cngwat2.op.rtp',h,ha,p,pa);
p = p0; p.cpsize2 = p.cpsize2*(1+dQ); rtpwrite('TEST_RTP/newdayx_1_100_12150_pert_cpsize2.op.rtp',h,ha,p,pa);
p = p0; p.cprtop2 = p.cprtop2*(1+dQ); rtpwrite('TEST_RTP/newdayx_1_100_12150_pert_cprtop2.op.rtp',h,ha,p,pa);
p = p0; p.cprbot2 = p.cprbot2*(1+dQ); rtpwrite('TEST_RTP/newdayx_1_100_12150_pert_cprbot2.op.rtp',h,ha,p,pa);

% see ccprep_slab.f
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=TEST_RTP/newdayx_1_100_12150.op.rtp               fout=TEST_RTP/newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 >& ugh0_cloudjac_cldfracsL
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=TEST_RTP/newdayx_1_100_12150_pert_cprtop1.op.rtp  fout=TEST_RTP/newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cprtop1_cloudjac_cldfracsL
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=TEST_RTP/newdayx_1_100_12150_pert_cprbot1.op.rtp  fout=TEST_RTP/newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cprbot1_cloudjac_cldfracsL
%}

[h,ha,p0,pa] = rtpread('TEST_RTP/newdayx_1_100_12150.op.rtp');
[h,p0] = subset_rtp_allcloudfields(h,p0,[],[],1);

% ccprep_slab.f
%               write(*,'(I5,6(F12.5))') L,CFRCL(L),CPRBOT,CPRTOP,PLEV(L),JACTOP_CFRCL(L),JACBOT_CFRCL(L)
%                                        1     2      3      4      5          6              7

raw_cldfracsL = load('UGH_TEST_OPTRAN/V6_May11_2023_cldjacs/ugh0_cloudjac_cldfracsL');
pert_cldfracsT1 = load('UGH_TEST_OPTRAN/V6_May11_2023_cldjacs/ugh_cprtop1_cloudjac_cldfracsL');
pert_cldfracsB1 = load('UGH_TEST_OPTRAN/V6_May11_2023_cldjacs/ugh_cprbot1_cloudjac_cldfracsL');

%% check that L is the same
ii= 1; [sum(raw_cldfracsL(:,ii)-pert_cldfracsT1(:,ii))] %sum(raw_cldfracsL(:,ii)-pert_cldfracsB1(:,ii))]
ii= 5; [sum(raw_cldfracsL(:,ii)-pert_cldfracsT1(:,ii))] %sum(raw_cldfracsL(:,ii)-pert_cldfracsB1(:,ii))]

%% the raw and T1 plevs are same
[raw_cldfracsL(:,5) pert_cldfracsT1(:,5)]

%% the raw and B1 plevs are different
[[0; raw_cldfracsL(:,5)] pert_cldfracsB1(:,5)]


dQ = 0.01;

ctop0 = p0.cprtop * 1;
ctopX = p0.cprtop * (1+dQ);

cbot0 = p0.cprbot * 1;
cbotX = p0.cprbot * (1+dQ);

[Lmax,~] = size(raw_cldfracsL);

%% yestimate = y0 + dy/dx deltax
for ind = 1 : Lmax
 %%%            L0                   y0                  y0 + dy/dx  dx                                              L1             y1
 %%%             1                   2                        3                                                       4              5
 junk = [raw_cldfracsL(ind,1) raw_cldfracsL(ind,2) raw_cldfracsL(ind,2)+raw_cldfracsL(ind,6)*dQ*ctop0 pert_cldfracsT1(ind,1) pert_cldfracsT1(ind,2)]; 
 junk = [junk (junk(5)-junk(3))*100/junk(3)];
 fprintf(1,'cldtop : ind L0 y0 yestimate y1 L1 percentdiff = %3i %3i %8.6f %8.6f %3i %8.6f %8.6f \n',ind,junk);
end

disp(' ')
for ind = 1 : Lmax
 junk = [raw_cldfracsL(ind,1) raw_cldfracsL(ind,2) raw_cldfracsL(ind,2)+raw_cldfracsL(ind,7)*dQ*cbot0 pert_cldfracsB1(ind,1) pert_cldfracsB1(ind,2)]; 
 junk = [junk (junk(5)-junk(3))*100/junk(3)];
 fprintf(1,'cldbot : ind L0 y0 yestimate y1 L1 percentdiff = %3i %3i %8.6f %8.6f %3i %8.6f %8.6f \n',ind,junk);
end
