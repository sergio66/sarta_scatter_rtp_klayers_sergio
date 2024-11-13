function [profnew,proftop] = add_othergases_arb_pressures(headin,profin,more_p,iVaryNMANGL_or_STD);

% used "add_othergases_arb_pressures.m" in 
%   /home/sergio/MATLABCODE/REGR_PROFILES_SARTA/ECMWF_SAF_137Profiles/ 
%   ECMWF_RTTOV_91_25000_and_7%04_Profiles/SAF704_2024/driver_look_at_Vers137_2024_rtpX.m

if nargin < 3
  error('need h,p,levs_additional')
end
if nargin == 3
  iVaryNMANGL_or_STD = +1;
end


more_p = reshape(more_p,length(more_p),1);

for ii = 1 : headin.ngas
  fprintf(1,'gasID = %2i gasunit = %2i \n',headin.glist(ii),headin.gunit(ii))
  if headin.gunit(ii) ~= 10 
    error('need gunit = 10')
  end
end

more_p = sort(more_p);
more_p = sort(more_p,'asc');
[mmx,nnx] = size(profin.plevs);

for ii = 1 : length(profin.stemp)
  nlevs = profin.nlevs(ii);
  plevs = profin.plevs(1:nlevs,ii);

  [Y,I] = sort(plevs,'asc');
  profin.plevs(1:nlevs,ii) = Y;

  ptemp = profin.ptemp(1:nlevs,ii);
  profin.ptemp(1:nlevs,ii) = ptemp(I);
  for gg = 1 : headin.ngas
    gasID = headin.glist(gg);
    str = ['gas = profin.gas_' num2str(gasID) '(1:nlevs,ii);'];
    eval(str);
    gas = gas(I);
    str = ['profin.gas_' num2str(gasID) '(1:nlevs,ii) = gas;'];
    eval(str);
  end
end

profnew = profin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% profnew.plevs                                        = zeros(mmx+length(more_p),nnx);
% profnew.plevs(1:length(more_p),:)                    = more_p * ones(1,nnx);
% profnew.plevs(length(more_p)+1:mmx+length(more_p),:) = profin.plevs;
% profnew.nlevs = profin.nlevs + length(more_p);
  
profnew.plevs = zeros(mmx+length(more_p),nnx);
for ii = 1 : length(profin.stemp)
  nlevs = profin.nlevs(ii);
  plevs = profin.plevs(1:nlevs,ii);
  pnewlevs = [more_p; plevs];
  profnew.plevs(1:nlevs + length(more_p),ii) = pnewlevs;
  profnew.nlevs(ii) = profin.nlevs(ii) + length(more_p);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if iVaryNMANGL_or_STD == -1
  disp('tacking on US Std to ALL profiles')
  iaSTD = 1:length(profin.stemp);
  iaTRP = [];
  iaMLS = [];
  iaMLW = [];
  iaSAS = [];
  iaSAW = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see add_afgl_g34.m
for gg = 1 : headin.ngas
  gasID = headin.glist(gg);
  fprintf(1,'processing gas %2i \n',gasID)

  clear miaow

  profnew.ptemp                                        = zeros(mmx+length(more_p),nnx);
  profnew.ptemp(length(more_p)+1:mmx+length(more_p),:) = profin.ptemp;

%  if gasID == 1
%    profnew.gas_1                                        = zeros(mmx+length(more_p),nnx);
%    profnew.gas_1(length(more_p)+1:mmx+length(more_p),:) = profin.gas_1;
%  elseif gasID == 2
%    profnew.gas_2                                        = zeros(mmx+length(more_p),nnx);
%    profnew.gas_2(length(more_p)+1:mmx+length(more_p),:) = profin.gas_2;
%  elseif gasID == 3
%    profnew.gas_3                                        = zeros(mmx+length(more_p),nnx);
%    profnew.gas_3(length(more_p)+1:mmx+length(more_p),:) = profin.gas_3;
%  elseif gasID == 6
%    profnew.gas_6                                        = zeros(mmx+length(more_p),nnx);
%    profnew.gas_6(length(more_p)+1:mmx+length(more_p),:) = profin.gas_6;
%  elseif gasID == 34
%    profnew.gas_34                                       = zeros(mmx+length(more_p),nnx);
%    profnew.gas_34(length(more_p)+1:mmx+length(more_p),:) = profin.gas_34;
%  end
  
  str = ['profnew.gas_' num2str(gasID) ' = zeros(mmx+length(more_p),nnx);'];                                          eval(str);
  str = ['profnew.gas_' num2str(gasID) ' (length(more_p)+1:mmx+length(more_p),:) = profin.gas_'  num2str(gasID) ';']; eval(str)

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

    if length(iaX) > 0
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
  
      xg = interp1(log(pX),log(qX),log(more_p),[],'extrap'); xg = exp(xg); 
      xt = interp1(log(pX),tX,log(more_p),[],'extrap');
      for pp = 1 : length(iaX)
        nlevs = profin.nlevs(iaX(pp));
        plevs = profin.plevs(1:nlevs,iaX(pp));
  
        xx = interp1(log(pX),log(qX),log(min(plevs)),[],'extrap'); xx = exp(xx); 
        offset =  GN(1,iaX(pp)) - xx;
        miaow(:,iaX(pp)) = xg + offset;
  
        if gg == 1
          xx = interp1(log(pX),tX,log(min(plevs)),[],'extrap');
          offset =  TN(1,iaX(pp)) - xx;
          miaowT(:,iaX(pp)) = xt + offset;
        end
      end
    else
      disp('mnafgl = 1,2,3,4,5,6 = RTP/MLS/MLW/SAS/SAW/STD')
      fprintf(1,'length(iaX) = 0 for iAtm/mnafgl = %2i \n',iAtm)
    end      %% if length(iaX) > 0
  end        %%  for iAtm = 1 : 6

  ahah = ['proftop.gas_' num2str(gasID) ' = miaow;'];  eval(ahah);
  if gg == 1
    ahah = ['proftop.ptemp                = miaowT;']; eval(ahah);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

proftop.plevs = sort(more_p,'asc') * ones(1,length(profin.rlon));
profnew.plevs(1:length(more_p),:)  = proftop.plevs;
profnew.ptemp(1:length(more_p),:)  = proftop.ptemp;

for pp = 1 : length(profin.stemp)
  for gg = 1 : headin.ngas
    gasID = headin.glist(gg);
    str = ['junk = proftop.gas_' num2str(gasID) '(:,pp);'];
    eval(str);
    good = find(junk > 0);
    bad  = find(junk <= 0);
    if length(bad) > 0
      minn = min(junk(good));
      junk(bad) = minn;
    end
    str = ['proftop.gas_' num2str(gasID) '(:,pp) = junk;'];
    eval(str);
  end
end

%profnew.gas_1(1:length(more_p),:)  = proftop.gas_1;
%profnew.gas_2(1:length(more_p),:)  = proftop.gas_2;
%profnew.gas_3(1:length(more_p),:)  = proftop.gas_3;
%profnew.gas_34(1:length(more_p),:) = proftop.gas_34;
for gg = 1 : headin.ngas
  gasID = headin.glist(gg);
  doer = ['profnew.gas_' num2str(gasID) '(1:length(more_p),:)  = eps + max(0,proftop.gas_' num2str(gasID) ');'];
  eval(doer);
end

plot_add_othergases_arb_pressures
