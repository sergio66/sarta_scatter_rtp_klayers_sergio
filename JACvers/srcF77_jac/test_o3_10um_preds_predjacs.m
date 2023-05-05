%{
% see lines 459 and 699 in ycalt2_od.f
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_nocldfields_pertO3x1.1.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=3  listc=1092 >& ughtestozonepredsA_pert
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_nocldfields.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=3  listc=1092 >& ughtestozonepredsA_0
%}

if ~exist('p')
  [h,ha,p0,pa] = rtpread('newdayx_nocldfields.op.rtp');
  [h,ha,p3,pa] = rtpread('newdayx_nocldfields_pertO3x1.1.op.rtp');
  nlays = p0.nlevs(1)-1;
  x = p0.gas_3(:,1); x = x(1:nlays); dx = x(1:nlays)/6.02214199E+26 * 0.1;   %% see rdrtp.f
  t = p0.ptemp(:,1); t = t(1:nlays); st = p0.stemp(1); mu = cos(p0.satzen(1)*pi/180);
end

predA_0    = load('ughtestozonepredsA_0');
predA_pert = load('ughtestozonepredsA_pert');

disp('XZ_O is not working : see ycalpar_jacINIT.f, line 212')
% CjacX      XZREF = XZREF + ROAMNT(L) 
% CjacX      XZ = XZ + POAMNT(L)       
% c      XZ_O = XZ/XZREF               

ii = 2; miaow3     = cumsum(predA_0(:,ii));
ii = 3; miaow3Ref  = cumsum(predA_0(:,ii));
ii = 7; junk = miaow3./miaow3Ref;   plot(predA_0(:,ii),1:97,'b.-',junk,1:97,'r'); set(gca,'ydir','reverse')
ii = 2; miaow3x    = cumsum(predA_pert(:,ii));
ii = 3; miaow3Refx = cumsum(predA_pert(:,ii));
ii = 7; junk = miaow3x./miaow3Refx; plot(predA_pert(:,ii),1:97,'b.-',junk,1:97,'r'); set(gca,'ydir','reverse')
ii = 7; junk = miaow3x./miaow3Refx; semilogx(predA_0(:,ii+1),1:97,'b.-',1./miaow3Refx,1:97,'ro-'); hl = legend('I predict jac 1/cumsum(Ozref)','this is hand 1/cumsum(Ozref)'); title('BUT THIS JAC IS WRONG .. WHY???')
ii = 7; junk = miaow3x./miaow3Refx; plot(predA_pert(:,ii),1:97,'b.-',junk,1:97,'r',predA_0(:,ii) + predA_0(:,ii+1).*dx,1:97,'k'); set(gca,'ydir','reverse')
ii = 7; semilogx((predA_pert(:,ii) - predA_0(:,ii))./dx,1:97,'b.-',predA_0(:,ii+1),1:97,'rx-',1./miaow3Refx,1:97,'k'); 
set(gca,'ydir','reverse'); title('XZ\_0 jac'); hl = legend('finite diff','predicted','1/cumsum(O3ref','fontsize',10,'location','best');

%% do it a little more carefully, only perturb layer IX and see how jac changes
ii = 2;
%% this is correct, and looks like MY derivative
%% saying derivative of XZ_0 = XZ/XZREF = 1/sum(qref) since only ONE layer perturbed
clear *oo*
for ix = 1 : 97
  oo0 = predA_0(:,ii);
  oo1 = oo0;
    oo1(ix) = predA_pert(ix,ii);
  moo = predA_0(:,ii);
  oo0  = oo0(1:ix);
  oo1 = oo1(1:ix);
  moo = moo(1:ix);
  a1 = sum(oo1)/sum(moo);
  a0 = sum(oo0)/sum(moo);
  ratio(ix) = (a1-a0)/dx(ix);
end
ii = 7; semilogx((predA_pert(:,ii) - predA_0(:,ii))./dx,1:97,'b.-',predA_0(:,ii+1),1:97,'rx-',1./miaow3Refx,1:97,'k',ratio,1:97,'gx-'); 
set(gca,'ydir','reverse'); title('XZ\_0 jac'); hl = legend('finite diff','predicted','1/cumsum(O3ref','haha is this right','fontsize',10,'location','best');

%% this is wrong, but looks like numerical derivative!!!! other than factor of 11 too large
clear *oo*
for ix = 1 : 97
  oo = predA_0(:,ii);
  oo(ix) = predA_pert(ix,ii);
  moo = predA_0(:,ii);
  oo = oo(1:ix);
  moo = moo(1:ix);
  ratio(ix) = sum(oo)/sum(moo)/dx(ix);
end
ii = 7; semilogx((predA_pert(:,ii) - predA_0(:,ii))./dx,1:97,'b.-',predA_0(:,ii+1),1:97,'rx-',1./miaow3Refx,1:97,'k',ratio/11,1:97,'gx-'); 
set(gca,'ydir','reverse'); title('XZ\_0 jac'); hl = legend('finite diff','predicted','1/cumsum(O3ref','haha is this right','fontsize',10,'location','best');

%% how about this, since we perturb EVERY layer above L!!!! 
%% saying derivative of XZ_0 = XZ/XZREF = sum(q)/sum(qref)/q(L)
clear *oo*
for ix = 1 : 97
  oo0 = predA_0(:,ii);
  oo1 = predA_pert(:,ii);
  moo = predA_0(:,ii);
  oo0  = oo0(1:ix);
  oo1 = oo1(1:ix);
  moo = moo(1:ix);
  ratio(ix) = (sum(oo1)-sum(oo0))/sum(moo)/dx(ix);
end
ii = 7; semilogx((predA_pert(:,ii) - predA_0(:,ii))./dx,1:97,'b.-',predA_0(:,ii+1),1:97,'rx-',1./miaow3Refx,1:97,'k',ratio,1:97,'gx-'); 
set(gca,'ydir','reverse'); title('XZ\_0 jac'); hl = legend('finite diff','predicted','1/cumsum(O3ref','haha is this right','fontsize',10,'location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii =  5; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) '  A\_0']);  %% yay
disp('ret to continue'); pause

ii =  7; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' XZ\_0']);  %% bah
disp('ret to continue'); pause

ii =  9; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' OJUNKA']); %% yay
disp('ret to continue'); pause

ii = 11; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' OJUNKR']); %% yay
disp('ret to continue'); pause

ii = 13; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' OJUNKZ']); %% bah
disp('ret to continue'); pause

ii = 15; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' OJUNKX']); %% bah
disp('ret to continue'); pause
