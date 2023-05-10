%{
[h,ha,p0,pa] = rtpread('newdayx_1_100_12150.op.rtp');
dQ = 0.1;
dQ = 0.01;
p = p0; p.cngwat = p.cngwat*(1+dQ);   rtpwrite('newdayx_1_100_12150_pert_cngwat1.op.rtp',h,ha,p,pa);
p = p0; p.cpsize = p.cpsize*(1+dQ);   rtpwrite('newdayx_1_100_12150_pert_cpsize1.op.rtp',h,ha,p,pa);
p = p0; p.cprtop = p.cprtop*(1+dQ);   rtpwrite('newdayx_1_100_12150_pert_cprtop1.op.rtp',h,ha,p,pa);
p = p0; p.cprbot = p.cprbot*(1+dQ);   rtpwrite('newdayx_1_100_12150_pert_cprbot1.op.rtp',h,ha,p,pa);

p = p0; p.cngwat2 = p.cngwat2*(1+dQ); rtpwrite('newdayx_1_100_12150_pert_cngwat2.op.rtp',h,ha,p,pa);
p = p0; p.cpsize2 = p.cpsize2*(1+dQ); rtpwrite('newdayx_1_100_12150_pert_cpsize2.op.rtp',h,ha,p,pa);
p = p0; p.cprtop2 = p.cprtop2*(1+dQ); rtpwrite('newdayx_1_100_12150_pert_cprtop2.op.rtp',h,ha,p,pa);
p = p0; p.cprbot2 = p.cprbot2*(1+dQ); rtpwrite('newdayx_1_100_12150_pert_cprbot2.op.rtp',h,ha,p,pa);

% see ccprep_slab.f
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150.op.rtp               fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 >& ugh0_cloudjac
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150_pert_cpsize1.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cpsize1_cloudjac
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150_pert_cpsize2.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cpsize2_cloudjac
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150_pert_cprtop1.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cprtop1_cloudjac
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150_pert_cprbot1.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cprbot1_cloudjac
%}

[h,ha,p0,pa] = rtpread('newdayx_1_100_12150.op.rtp');
[h,p0] = subset_rtp_allcloudfields(h,p0,[],[],1);

% docloudyTwoSlab_RT.f
%        write(*,'(A,I4,9(F12.4))') 'check taucldOD and radpure',IJUNK,TAU4(1,IJUNK,I),RAACLDOD4(1:4,IJUNK,I),PURE_RAD4(1:4,IJUNK,I), PLANK_RAD4(1:4,IJUNK,I)
%                                                                  1       2                  3-6                7-10                   11-14

raw_ODs = load('ugh0_cloudjac'); %% see test_cld_jacs but only sze,cng jacs SO A WASTE OF TIME HERE

raw_rads = load('ugh0_cloudjac_rad');
pert_rads1 = load('ugh_cprtop1_cloudjac');
pert_rads2 = load('ugh_cprbot1_cloudjac');

dQ = 0.01;

dme0 = p0.cpsize * 1;
dmeX = p0.cpsize * (1+dQ);

dme20 = p0.cpsize2 * 1;
dme2X = p0.cpsize2 * (1+dQ);

clf
plot(pert_rads1(:,7) -raw_rads(:,7), 1:97, pert_rads2(:,7) -raw_rads(:,7), 1:97); set(gca,'ydir','reverse')  %% clr rad
plot(pert_rads1(:,8) -raw_rads(:,8), 1:97, pert_rads2(:,8) -raw_rads(:,8), 1:97); set(gca,'ydir','reverse')  %% cld1  rad   ***
plot(pert_rads1(:,9) -raw_rads(:,9), 1:97, pert_rads2(:,9) -raw_rads(:,9), 1:97); set(gca,'ydir','reverse')  %% cld2  rad
plot(pert_rads1(:,10)-raw_rads(:,10),1:97, pert_rads2(:,10) -raw_rads(:,10),1:97); set(gca,'ydir','reverse')  %% cld12 rad

%% dr/dX   = dr/dtau dtau/dX
%% dr/dCNG = dr/dtau dtau/dCNG
%%   1/mu dtau/dCNG = raw_ODs(1,13)  from test_cld_jacs.m
%%   mu dr/dtau     = -r + B                
dtau_dCNG004 =  zeros(97,1); dtau_dCNG004(57:61) = raw_ODs(1,13);
dtau_dCNGALL =  zeros(97,1); dtau_dCNGALL(01:61) = raw_ODs(1,13);
dr_dtau   =  raw_rads(:,11); 
dr_dCNG004 = dr_dtau .* dtau_dCNG004 * 0.25;  %% remember we have 4 cloud layers so CNG is divided into 4
dr_dCNGALL = dr_dtau .* dtau_dCNGALL * 0.25;  %% remember we have 4 cloud layers so CNG is divided into 4

plot(pert_rads1(:,8)-raw_rads(:,8), 1:97, pert_rads2(:,8) -raw_rads(:,8), 1:97, dr_dCNG004, 1:97, 'kx-', dr_dCNGALL, 1:97, 'go-'); set(gca,'ydir','reverse')  %% cld1  rad   ***
plot(rad2bt(1231,1000*pert_rads1(:,8))-rad2bt(1231,1000*raw_rads(:,8)), 1:97); set(gca,'ydir','reverse')       %% cld1  rad   ***  sarta does things in W/m2/sr-1/cm-1 so change to mW
plot(rad2bt(1231,1000*pert_rads1(:,8))-rad2bt(1231,1000*raw_rads(:,8)), 1:97,rad2bt(1231,1000*pert_rads2(:,8))-rad2bt(1231,1000*raw_rads(:,8)), 1:97); set(gca,'ydir','reverse'); plotaxis2
plot(rad2bt(1231,1000*pert_rads1(:,8)),1:97,rad2bt(1231,1000*raw_rads(:,8)), 1:97); set(gca,'ydir','reverse')  %% cld1  rad   ***  sarta does things in W/m2/sr-1/cm-1 so change to mW
plot(rad2bt(1231,1000*pert_rads1(:,8)),1:97,rad2bt(1231,1000*pert_rads2(:,8)),1:97,rad2bt(1231,1000*raw_rads(:,8)), 1:97); set(gca,'ydir','reverse')  %% cld1  rad   ***  sarta does things in W/m2/sr-1/cm-1 so change to mW
