%{
% see ycalpar_WVjacpredsINC.f
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_nocldfields_pertWVx1.1.op.rtp  fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1603 > ughtestwaterpredsA_pert
 time ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug fin=newdayx_nocldfields.op.rtp             fout=newdayx_1_100_12150.rad.rtp         listp=1  listj=1  listc=1603 > ughtestwaterpredsA_0
%}

if ~exist('p')
  [h,ha,p0,pa] = rtpread('newdayx_nocldfields.op.rtp');
  [h,ha,p3,pa] = rtpread('newdayx_nocldfields_pertWVx1.1.op.rtp');
  nlays = p0.nlevs(1)-1;
  x = p0.gas_1(:,1); x = x(1:nlays); dx = x(1:nlays)/6.02214199E+26 * 0.1;   %% see rdrtp.f
  t = p0.ptemp(:,1); t = t(1:nlays); st = p0.stemp(1); mu = cos(p0.satzen(1)*pi/180);
end

predA_0    = load('ughtestwaterpredsA_0');
predA_pert = load('ughtestwaterpredsA_pert');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii =  5; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) '  A\_W']);  %% yay
disp('ret to continue'); pause

ii =  7; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' AZ\_W']);  %% bah
disp('ret to continue'); pause

ii =  9; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' WJUNKA']); %% yay
disp('ret to continue'); pause

ii = 11; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' WJUNKR']); %% yay
disp('ret to continue'); pause

ii = 13; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' WJUNKS']); %% bah
disp('ret to continue'); pause

ii = 15; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' WJUNKZ']); %% bah
disp('ret to continue'); pause

ii = 17; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' WJUNK4']); %% bah
disp('ret to continue'); pause

ii = 19; plot(predA_0(:,ii),1:97,'b',predA_pert(:,ii),1:97,'r.-',predA_0(:,ii) + predA_0(:,ii+1) .* dx,1:97,'ko-'); 
set(gca,'ydir','reverse'); legend('raw','pert','raw + deriv*dq','location','best','fontsize',10); title([num2str(ii) ' MJUNKZ']); %% bah
disp('ret to continue'); pause
