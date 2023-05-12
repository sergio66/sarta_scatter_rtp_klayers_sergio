%{
[h,ha,p0,pa] = rtpread('TEST_RTP/newdayx_1_100_12150.op.rtp');
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
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=TEST_RTP/newdayx_1_100_12150.op.rtp               fout=TEST_RTP/newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 >& ugh0_cloudjac
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=TEST_RTP/newdayx_1_100_12150_pert_cpsize1.op.rtp  fout=TEST_RTP/newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cpsize1_cloudjac
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=TEST_RTP/newdayx_1_100_12150_pert_cpsize2.op.rtp  fout=TEST_RTP/newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cpsize2_cloudjac
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=TEST_RTP/newdayx_1_100_12150_pert_cprtop1.op.rtp  fout=TEST_RTP/newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cprtop1_cloudjac
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=TEST_RTP/newdayx_1_100_12150_pert_cprbot1.op.rtp  fout=TEST_RTP/newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1291 > ugh_cprbot1_cloudjac
%}

addpath MATLABCODE

[h,ha,p0,pa] = rtpread('TEST_RTP/newdayx_1_100_12150.op.rtp');
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

figure(1); clf

% plot(pert_rads1(:,7) -raw_rads(:,7), 1:97,'b.-',pert_rads2(:,7) -raw_rads(:,7), 1:97,'r');  plotaxis2; hl=legend('cprtop pert','cprbot pert','location','best','fontsize',10); title('Pert-Raw'); set(gca,'ydir','reverse')  %% clr rad
plot(pert_rads1(:,8) -raw_rads(:,8), 1:97,'b.-',pert_rads2(:,8) -raw_rads(:,8), 1:97,'r');  plotaxis2; hl=legend('cprtop pert','cprbot pert','location','best','fontsize',10); title('Pert-Raw'); set(gca,'ydir','reverse')  %% cld1  rad  ******
% plot(pert_rads1(:,9) -raw_rads(:,9), 1:97,'b.-',pert_rads2(:,9) -raw_rads(:,9), 1:97,'r');  plotaxis2; hl=legend('cprtop pert','cprbot pert','location','best','fontsize',10); title('Pert-Raw'); set(gca,'ydir','reverse')  %% cld2  rad
% plot(pert_rads1(:,10)-raw_rads(:,10),1:97,'b.-',pert_rads2(:,10) -raw_rads(:,10),1:97,'r'); plotaxis2; hl=legend('cprtop pert','cprbot pert','location','best','fontsize',10); title('Pert-Raw'); set(gca,'ydir','reverse')  %% cld12 rad **

%% write(*,'(I5,6(F12.5))') L,CFRCL(L),CPRBOT,CPRTOP,PLEV(L),JACTOP_CFRCL(L),JACBOT_CFRCL(L)  .. first half (rows 1:N/2) are cloud 1, second half (N/2+1:N) are cloud two 
frac_jacs = load('frac_jacs.txt'); %% ccpreop.f
disp('showing cprtop fracjacs')
[frac_jacs(1:5,1) frac_jacs(1:5,6)]

%% dr/dX   = dr/dtau dtau/dX
%% dr/dCNG = dr/dtau dtau/dCNG
%%   1/mu dtau/dCNG = raw_ODs(1,13)  from test_cld_jacs.m
%%   mu dr/dtau     = -r + B                

%% who knows what this is?????
%% dtau_dCNG004 =  zeros(97,1); dtau_dCNG004(57:61) = raw_ODs(1,13);
%% dtau_dCNGALL =  zeros(97,1); dtau_dCNGALL(01:61) = raw_ODs(1,13);

%% remember cldOD(cng) = cng * cldOD(1) ==> dtau/dCNG = cldOD(1) = raw_rads(4)/cngwat (for cld1) and raw_rads(5)/cngwat2 (for cld2)
dtau_dCNG004 =  zeros(97,1); dtau_dCNG004 = raw_rads(:,4)/p0.cngwat;  
dtau_dCNGALL =  zeros(97,1); dtau_dCNGALL = raw_rads(:,4)*mean(abs(frac_jacs(:,6)));

dr_dtau   =  raw_rads(:,11);                   %% - r + B
%dr_dCNG004 = dr_dtau .* dtau_dCNG004 * 0.25;  %% remember we have 4 cloud layers so CNG is divided into 4
%dr_dCNGALL = dr_dtau .* dtau_dCNGALL * 0.25;  %% remember we have 4 cloud layers so CNG is divided into 4
dr_dCNG004 = dr_dtau .* dtau_dCNG004;
dr_dCNGALL = dr_dtau .* dtau_dCNGALL;

figure(2); clf
plot(pert_rads1(:,8)-raw_rads(:,8), 1:97, 'bx-', pert_rads2(:,8) -raw_rads(:,8), 1:97, 'ro-', dr_dCNG004, 1:97, 'kx-', dr_dCNGALL, 1:97, 'go-', 'linewidth',2); plotaxis2; set(gca,'ydir','reverse');  %% cld1  rad   ***
hl = legend('cprtop pert-raw','cprbot pert-raw','guess 1','guess 2','location','best','fontsize',10); title('Cld 1 rad differences')'

figure(3); clf
plot(rad2bt(1231,1000*pert_rads1(:,8))-rad2bt(1231,1000*raw_rads(:,8)), 1:97); plotaxis2; set(gca,'ydir','reverse')       %% cld1  rad   ***  sarta does things in W/m2/sr-1/cm-1 so change to mW
plot(rad2bt(1231,1000*pert_rads1(:,8))-rad2bt(1231,1000*raw_rads(:,8)), 1:97,rad2bt(1231,1000*pert_rads2(:,8))-rad2bt(1231,1000*raw_rads(:,8)), 1:97); plotaxis2; set(gca,'ydir','reverse'); plotaxis2;
plot(rad2bt(1231,1000*pert_rads1(:,8)),1:97,rad2bt(1231,1000*raw_rads(:,8)), 1:97); plotaxis2; set(gca,'ydir','reverse')  %% cld1  rad   ***  sarta does things in W/m2/sr-1/cm-1 so change to mW
plot(rad2bt(1231,1000*pert_rads1(:,8)),1:97,rad2bt(1231,1000*pert_rads2(:,8)),1:97,rad2bt(1231,1000*raw_rads(:,8)), 1:97); plotaxis2; set(gca,'ydir','reverse')  %% cld1  rad   ***  sarta does things in W/m2/sr-1/cm-1 so change to mW

%% see calrad1.f
%%          print *,'calrad1',L,RADUP,TAULX(L),RPLNCK(L),RADUP*TAULX(L) + RPLNCK(L)*(1.0 - TAULX(L)) + RSUNSC
x0 = load('ugh0_cloudjac_radX');
xT = load('ugh_cprtop1_cloudjacX');
xB = load('ugh_cprbot1_cloudjacX');

iaL = xT(1:97,1); 
%% remember we call calrad1 for cloud2 before cloud 1 ...... and we only perturbed cloud 1. So use (1:97)+97
ii=3; figure(4); plot(xT((01:97)+97,ii)-x0((01:97)+97,ii),iaL,'bx-',xB((01:97)+97,ii)-x0((01:97)+97,ii),iaL,'r.-'); title('pert cld 1 TAUX(L)-TAUX0(L)');   plotaxis2; hl = legend('cprtop','cprbot','location','best'); set(gca,'ydir','reverse');
ii=4; figure(4); plot(xT((01:97)+97,ii)-x0((01:97)+97,ii),iaL,'bx-',xB((01:97)+97,ii)-x0((01:97)+97,ii),iaL,'r.-'); title('pert cld 1 PLNCK(L)-PLNCK0(L)'); plotaxis2; hl = legend('cprtop','cprbot','location','best'); set(gca,'ydir','reverse');

ii=5; figure(4); plot(xT((01:97)+97,ii)-x0((01:97)+97,ii),iaL,'bx-',xB((01:97)+97,ii)-x0((01:97)+97,ii),iaL,'r.-'); title('pert cld 1 RAD(L)-RAD0(L)');     plotaxis2; 
[Y,I] = sort(frac_jacs(1:5,1)); figure(4); hold on; plot(frac_jacs(I(1:5),6)/100,frac_jacs(I(1:5),1),'k+-'); hold off
[Y,I] = sort(frac_jacs(1:5,1)); figure(4); hold on; plot(frac_jacs(I(1:5),7)/100,frac_jacs(I(1:5),1),'go-'); hold off
hl = legend('cprtop','cprbot','location','best'); set(gca,'ydir','reverse'); 
ylim([55 65])
