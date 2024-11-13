function [profnew,proftop] = add_othergases_arb_pressures(headin,profin,more_p);

% used "add_othergases_arb_pressures.m" in 
%   /home/sergio/MATLABCODE/REGR_PROFILES_SARTA/ECMWF_SAF_137Profiles/ 
%   ECMWF_RTTOV_91_25000_and_7%04_Profiles/SAF704_2024/driver_look_at_Vers137_2024_rtpX.m

profnew = profin;
for ii = 1 : headin.ngas
  fprintf(1,'gasID = %2i gasunit = %2i \n',headin.glist(ii),headin.gunit(ii))
  if headin.gunit(ii) ~= 10 
    error('need gunit = 10')
  end
end

more_p = sort(more_p);
[mmx,nnx] = size(profin.plevs);

profnew.plevs                                       = zeros(mmx+length(more_p),nnx);
profnew.plevs(1:length(more_p),:)                   = more_p * ones(1,nnx);
profnew.plevs(length(more_p)+1:mmx+length(more_p),:) = profin.plevs;
profnew.nlevs = profin.nlevs + length(more_p);

[yy,mm,dd,hh] = tai2utcSergio(profin.rtime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iaTRP = find(abs(profin.rlat) < 30);

iaMLS_N = find(profin.rlat < +60  & profin.rlat >= +30 & mm >= 6 & mm <= 9);
iaMLS_S = find(profin.rlat > -60  & profin.rlat <= -30 & mm >= 1 & mm <= 4);
iaMLS = union(iaMLS_N,iaMLS_S);

iaMLW_N = find(profin.rlat < +60  & profin.rlat >= +30 & mm >= 1 & mm <= 4);
iaMLW_S = find(profin.rlat > -60  & profin.rlat <= -30 & mm >= 6 & mm <= 9);
iaMLW = union(iaMLW_N,iaMLW_S);

iaSTD = find(abs(profin.rlat) < 60 & abs(profin.rlat) >= 30);
iaSTD = setdiff(setdiff(iaSTD,iaMLS),iaMLW);

iaSAS_N = find(profin.rlat >= +60 & mm >= 4 & mm <= 9);
iaSAS_S = find(profin.rlat <= -60 & ((mm >= 10 & mm <= 12) | (mm >= 1 & mm <= 3)));
iaSAS = union(iaSAS_N,iaSAS_S);

iaSAW_N = find(profin.rlat >= +60 & ((mm >= 10 & mm <= 12) | (mm >= 1 & mm <= 3)));
iaSAW_S = find(profin.rlat <= -60 & mm >= 4 & mm <= 9);
iaSAW = union(iaSAW_N,iaSAW_S);

%length(iaTRP) + length(iaSTD) + length(iaMLS) + length(iaMLW) + length(iaSAW) + length(iaSAS)
wah = setdiff(1:length(profin.stemp),[iaTRP iaMLS iaMLW iaSAS iaSAW iaSTD]); 
if length(wah) > 0
  error('oopsy left out some FOVS')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see add_afgl_g34.m
for gg = 1 : headin.ngas
  gasID = headin.glist(gg);
  fprintf(1,'processing gas %2i \n',gasID)

  clear miaow

  profnew.ptemp                                        = zeros(mmx+length(more_p),nnx);
  profnew.ptemp(length(more_p)+1:mmx+length(more_p),:) = profin.ptemp;

  if gasID == 1
    profnew.gas_1                                        = zeros(mmx+length(more_p),nnx);
    profnew.gas_1(length(more_p)+1:mmx+length(more_p),:) = profin.gas_1;
  elseif gasID == 2
    profnew.gas_2                                        = zeros(mmx+length(more_p),nnx);
    profnew.gas_2(length(more_p)+1:mmx+length(more_p),:) = profin.gas_2;
  elseif gasID == 3
    profnew.gas_3                                        = zeros(mmx+length(more_p),nnx);
    profnew.gas_3(length(more_p)+1:mmx+length(more_p),:) = profin.gas_3;
  elseif gasID == 4
    profnew.gas_4                                        = zeros(mmx+length(more_p),nnx);
    profnew.gas_4(length(more_p)+1:mmx+length(more_p),:) = profin.gas_4;
  elseif gasID == 5
    profnew.gas_5                                        = zeros(mmx+length(more_p),nnx);
    profnew.gas_5(length(more_p)+1:mmx+length(more_p),:) = profin.gas_5;
  elseif gasID == 6
    profnew.gas_6                                        = zeros(mmx+length(more_p),nnx);
    profnew.gas_6(length(more_p)+1:mmx+length(more_p),:) = profin.gas_6;
  elseif gasID == 9
    profnew.gas_9                                        = zeros(mmx+length(more_p),nnx);
    profnew.gas_9(length(more_p)+1:mmx+length(more_p),:) = profin.gas_9;
  elseif gasID == 34
    profnew.gas_34                                       = zeros(mmx+length(more_p),nnx);
    profnew.gas_34(length(more_p)+1:mmx+length(more_p),:) = profin.gas_34;
  end

  ahah = ['GN = profin.gas_' num2str(gasID) ';']; eval(ahah);
  if gg == 1
    ahah = ['TN = profin.ptemp;']; eval(ahah);
  end

  for iAtm = 1 : 6
    if iAtm == 1 
      iaX = iaTRP;
    elseif iAtm == 2
      iaX = iaMLS;
    elseif iAtm == 3
      iaX = iaMLW;
    elseif iAtm == 4
      iaX = iaSAS;
    elseif iAtm == 5
      iaX = iaSAW;
    elseif iAtm == 6
      iaX = iaSTD;
    end

    afglX = quick_read_afgl(gasID,iAtm);
    if gasID <= 7
      pX = afglX.piAtm;
      qX = afglX.qiAtm;
      tX = afglX.tiAtm;
    else
      pX = afglX.pstd;
      qX = afglX.qstd;
      tX = afglX.tstd;
    end

%if gasID == 3
%  figure(6); clf; loglog(qX,pX); set(gca,'ydir','reverse')
%  figure(1)
%end

    x2 = interp1(log(pX),log(qX),log(more_p),[],'extrap'); x2 = exp(x2); 
    xx = interp1(log(pX),log(qX),log(min(profin.plevs(:))),[],'extrap'); xx = exp(xx); 
    for pp = 1 : length(iaX)
      offset =  GN(1,iaX(pp)) - xx;
      miaow(:,iaX(pp)) = x2 + offset;
    end

    if gg == 1
      x2 = interp1(log(pX),tX,log(more_p),[],'extrap');
      xx = interp1(log(pX),tX,log(min(profin.plevs(:))),[],'extrap');
      for pp = 1 : length(iaX)
        offset =  TN(1,iaX(pp)) - xx;
        miaowT(:,iaX(pp)) = x2 + offset;
      end
    end

  end
  ahah = ['proftop.gas_' num2str(gasID) ' = miaow;'];  eval(ahah);
  if gg == 1
    ahah = ['proftop.ptemp                = miaowT;']; eval(ahah);
  end

end
proftop.plevs = more_p * ones(1,length(profin.rlon));
profnew.plevs(1:length(more_p),:)  = proftop.plevs;
profnew.ptemp(1:length(more_p),:)  = proftop.ptemp;
profnew.gas_1(1:length(more_p),:)  = proftop.gas_1;
profnew.gas_2(1:length(more_p),:)  = proftop.gas_2;
profnew.gas_3(1:length(more_p),:)  = proftop.gas_3;
profnew.gas_5(1:length(more_p),:)  = proftop.gas_5;
profnew.gas_6(1:length(more_p),:)  = proftop.gas_6;
profnew.gas_9(1:length(more_p),:)  = proftop.gas_9;
%profnew.gas_34(1:length(more_p),:) = proftop.gas_34;

plot_add_othergases_arb_pressures

% figure(1); clf; 
%   semilogy(nanmean(profin.ptemp,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.ptemp,2),nanmean(proftop.plevs,2),'r.-',...
%            nanmean(profnew.ptemp,2),nanmean(profnew.plevs,2),'k'); 
%   set(gca,'ydir','reverse')
%   title('T(z)')
% 
% figure(2); clf; 
%   loglog(nanmean(profin.gas_1,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_1,2),nanmean(proftop.plevs,2),'r.-',...
%         nanmean(profnew.gas_1,2),nanmean(profnew.plevs,2),'k'); 
%   set(gca,'ydir','reverse')
%   title('WV(z)')
% 
% figure(3); clf; 
%   loglog(nanmean(profin.gas_2,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_2,2),nanmean(proftop.plevs,2),'r.-',...
%          nanmean(profnew.gas_2,2),nanmean(profnew.plevs,2),'k'); 
%   set(gca,'ydir','reverse')
%   title('CO2(z)')
% 
% figure(4); clf; 
%   loglog(nanmean(profin.gas_3,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_3,2),nanmean(proftop.plevs,2),'r.-',...
%          nanmean(profnew.gas_3,2),nanmean(profnew.plevs,2),'k'); 
%   set(gca,'ydir','reverse')
%   title('O3(z)')
% 
% figure(5); clf; 
%   loglog(nanmean(profin.gas_34,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_34,2),nanmean(proftop.plevs,2),'r.-',...
%          nanmean(profnew.gas_34,2),nanmean(profnew.plevs,2),'k'); 
%   set(gca,'ydir','reverse')
%   title('O(z)')
% 
