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

% ccprep_slab.f, will have two rows, one for cld1, one for cld2
%       I = 1 
%       write (*,'(A,I4,12(F12.4))') 'cld jacs',I,  NEXTOD(I),NSCAOD(I),G_ASYM(I),  FINALOD(I)     JACA_NEXTOD(I),JACA_NSCAOD(I),JACA_G_ASYM(I),JACA_FINAL(I),    JACS_NEXTOD(I),JACS_NSCAOD(I),JACS_G_ASYM(I),JACS_FINAL(I)
%                                               1     2         3         4             5              6              7              8                  9            10             11             12             13
raw_ODs = load('ugh0_cloudjac');
pert_ODs1 = load('ugh_cpsize1_cloudjac');
pert_ODs2 = load('ugh_cpsize2_cloudjac');

dQ = 0.01;

dme0 = p0.cpsize * 1;
dmeX = p0.cpsize * (1+dQ);

dme20 = p0.cpsize2 * 1;
dme2X = p0.cpsize2 * (1+dQ);

% [y0      y0+deltay=y0+dy/dx deltax   y1actual        percent difference (y1est-y1actual)*100/y1actual   %%
junk = [raw_ODs(1,2) raw_ODs(1,2)+raw_ODs(1,10)*dme0*dQ  pert_ODs1(1,2) (raw_ODs(1,2)+raw_ODs(1,10)*dme0*dQ - pert_ODs1(1,2))*100/pert_ODs1(1,2)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NEXTOD01  NEXT1_ESTMATED     NEXTOD1
junk = [raw_ODs(1,3) raw_ODs(1,3)+raw_ODs(1,11)*dme0*dQ  pert_ODs1(1,3) (raw_ODs(1,3)+raw_ODs(1,11)*dme0*dQ - pert_ODs1(1,3))*100/pert_ODs1(1,3)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NSXTOD01  NSXT1_ESTMATED     NSXTOD1
junk = [raw_ODs(1,4) raw_ODs(1,4)+raw_ODs(1,12)*dme0*dQ  pert_ODs1(1,4) (raw_ODs(1,4)+raw_ODs(1,12)*dme0*dQ - pert_ODs1(1,4))*100/pert_ODs1(1,4)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NAXTOD01  NAXT1_ESTMATED     NAXTOD1
junk = [raw_ODs(1,5) raw_ODs(1,5)+raw_ODs(1,13)*dme0*dQ  pert_ODs1(1,5) (raw_ODs(1,5)+raw_ODs(1,13)*dme0*dQ - pert_ODs1(1,5))*100/pert_ODs1(1,5)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NFXTOD01  NFXT1_ESTMATED     NEXTOD1
disp(' ')

junk = [raw_ODs(2,2) raw_ODs(2,2)+raw_ODs(2,10)*dme20*dQ  pert_ODs2(2,2) (raw_ODs(2,2)+raw_ODs(2,10)*dme0*dQ - pert_ODs1(2,2))*100/pert_ODs1(2,2)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NEXTOD02  NEXT2_ESTMATED     NEXTOD2
junk = [raw_ODs(2,3) raw_ODs(2,3)+raw_ODs(2,11)*dme20*dQ  pert_ODs2(2,3) (raw_ODs(2,3)+raw_ODs(2,11)*dme0*dQ - pert_ODs1(2,3))*100/pert_ODs1(2,3)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NSXTOD02  NSXT2_ESTMATED     NSXTOD2
junk = [raw_ODs(2,4) raw_ODs(2,4)+raw_ODs(2,12)*dme20*dQ  pert_ODs2(2,4) (raw_ODs(2,4)+raw_ODs(2,12)*dme0*dQ - pert_ODs1(2,4))*100/pert_ODs1(2,4)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NAXTOD02  NAXT2_ESTMATED     NAXTOD2
junk = [raw_ODs(2,5) raw_ODs(2,5)+raw_ODs(2,13)*dme20*dQ  pert_ODs2(2,5) (raw_ODs(2,5)+raw_ODs(2,13)*dme0*dQ - pert_ODs1(2,5))*100/pert_ODs1(2,5)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NFXTOD02  NFXT2_ESTMATED     NFXTOD2
disp(' ')

%{
time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_1_100_12150.op.rtp               fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=300 >& ugh
[w,d] = readsarta_jac('newdayx_1_100_12150.rad.rtp_jacCLD',300); d = squeeze(d(1,:,:));
jacx = quicksartajac(h,p0,1,1); 

clf

%% cfrac
 ii=1; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));
 ii=6; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));
 ii=11; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));

%% cngwat
 ii=2; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));
 ii=7; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));

%% cpsize
 ii=3; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));
 ii=8; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));
 ii=3; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii)); [raw_ODs(1,5) raw_ODs(1,5)+raw_ODs(1,13)*dme0*dQ   pert_ODs1(1,5)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NEXTOD0  NEXT1_ESTMATED     NEXTOD1
 ii=8; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii)); [raw_ODs(2,5) raw_ODs(2,5)+raw_ODs(2,13)*dme20*dQ  pert_ODs2(2,5)]; fprintf(1,'%8.5f %8.6f %8.6f %8.6f \n',junk);      %%  NEXTOD0  NEXT1_ESTMATED     NEXTOD1

%% cprtop
 ii=4; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));
 ii=9; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));

%% cprbot
 ii=05; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));
 ii=10; plot(w,jacx.cldjac(:,ii),'b.-',w,d(:,ii));


%}
