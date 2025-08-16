addpath /home/sergio/git/sarta_utilities/
addpath /home/sergio/MATLABCODE/COLORMAP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [ichan, fchan, coef, info] = rdcoef(set, lmerged, fname);
%
% Reads in binary FORTRAN coefficient data file
%
% Input:
%    set     = {integer} set number [1:7} and 8=CO2/SO2/HNO3/NH3 4term,
%              9=OPTRAN, 10=N2O 7term, 11=therm/nonLTE 12=CO2/SO2/HNO3 5term
%              13=HDO 11term
%    lmerged = {integer} 0=raw (not merged), 1=merged with water continuum
%    fname   = {string} name of binary FORTRAN data file to read
%
% Output:
%    ichan   = [nchan x 1] channel ID (JPL numbering)
%    fchan   = [nchan x 1] center frequency (wavenumber)
%    coef    = [nchan x nlay x ncoef] coeffcients
%    info    = {structure} various info
%
% Warning: if set=8 then info.gasid=2 even if SO2(9) or HNO3(12)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coef_vers0 = '../SymbolicLinks/L2_100/AIRS/2022/';
coef_vers1 = '../SymbolicLinks/L2_100/AIRS/2025/ECM83/';
coef_vers1 = '../SymbolicLinks/L2_100/AIRS/2025/SAF703/';

dir0 = dir(coef_vers0);
dir1 = dir(coef_vers1);

fprintf(1,' dir0 has %2i files \n dir1 has %2i files \n',length(dir0),length(dir1))
disp(' ')
iCnt = 0;
for ii = 1 : length(dir0)
  disp('  ')
  zname0 = dir0(ii).name;
  
  if ~strcmp(zname0,'.') & ~strcmp(zname0,'..') & ~strcmp(zname0,'refprof_trace400')
    iFound = -1;
    jj = 1;
    while iFound < 0 & jj <= length(dir1)
      zname1 = dir1(jj).name;
      if strcmp(zname0,zname1)
        iCnt = iCnt + 1;
        iFound = +1;
        fprintf(1,'  %s %s \n',zname0,zname1)
	
	if strfind(zname0,'set')
          setnum = str2num(zname0(4:4));
	  lmerged = 1;
          [ichan0, fchan0, coef0, info0] = rdcoef(setnum, lmerged, [coef_vers0 '/' zname0]);
          [ichan1, fchan1, coef1, info1] = rdcoef(setnum, lmerged, [coef_vers1 '/' zname1]);
	  
	%elseif strfind(zname0,'co2') | strfind(zname0,'so2') | strfind(zname0,'hno3') | strfind(zname0,'nh3')
	elseif strfind(zname0,'nh3')	
          setnum = 8;
	  lmerged = 1;
          [ichan0, fchan0, coef0, info0] = rdcoef(setnum, lmerged, [coef_vers0 '/' zname0]);
          [ichan1, fchan1, coef1, info1] = rdcoef(setnum, lmerged, [coef_vers1 '/' zname1]);

	elseif strfind(zname0,'optran')
          setnum = 9;
	  lmerged = 1;
          [ichan0, fchan0, coef0, info0] = rdcoef(setnum, lmerged, [coef_vers0 '/' zname0]);
          [ichan1, fchan1, coef1, info1] = rdcoef(setnum, lmerged, [coef_vers1 '/' zname1]);

	elseif strfind(zname0,'n2o')
          setnum = 10;
	  lmerged = 1;
          [ichan0, fchan0, coef0, info0] = rdcoef(setnum, lmerged, [coef_vers0 '/' zname0]);
          [ichan1, fchan1, coef1, info1] = rdcoef(setnum, lmerged, [coef_vers1 '/' zname1]);

	elseif strfind(zname0,'refl') | strfind(zname0,'nte')
          setnum = 11;
	  lmerged = 1;
          [ichan0, fchan0, coef0, info0] = rdcoef(setnum, lmerged, [coef_vers0 '/' zname0]);
          [ichan1, fchan1, coef1, info1] = rdcoef(setnum, lmerged, [coef_vers1 '/' zname1]);

	elseif strfind(zname0,'co2') | strfind(zname0,'so2') | strfind(zname0,'hno3')
          setnum = 12;
	  lmerged = 1;
          [ichan0, fchan0, coef0, info0] = rdcoef(setnum, lmerged, [coef_vers0 '/' zname0]);
          [ichan1, fchan1, coef1, info1] = rdcoef(setnum, lmerged, [coef_vers1 '/' zname1]);

        elseif strfind(zname0,'hdo')
          setnum = 13;
	  lmerged = 1;
          [ichan0, fchan0, coef0, info0] = rdcoef(setnum, lmerged, [coef_vers0 '/' zname0]);
          [ichan1, fchan1, coef1, info1] = rdcoef(setnum, lmerged, [coef_vers1 '/' zname1]);
        end
	
        err = 0;
        x = size(ichan0); y = size(ichan1);
          if x ~= y
	    err = err + 1;
	    disp('ichan not same size');
	  end
        x = size(fchan0); y = size(fchan1);
	  if x ~= y
	    err = err + 1;
	    disp('fchan not same size');
	  end
        x = size(coef0); y = size(coef1);
	  if x ~= y
	    err = err + 1;
	    disp('coef not same size');
	  end
	  
	%boo = compare_two_structures(info0,info1);

        nlay = 100;
	if strfind(zname0,'optran')
	  nlay = 300;
	end
	
        figure(1); pcolor(fchan0,1:nlay,squeeze(coef0(:,:,1))');                set(gca,'ydir','reverse'); shading interp; colorbar; title('A');   colormap jet
        figure(2); pcolor(fchan1,1:nlay,squeeze(coef1(:,:,1))');                set(gca,'ydir','reverse'); shading interp; colorbar; title('B');   colormap jet
        figure(3); pcolor(fchan0,1:nlay,squeeze(coef0(:,:,1)'./coef1(:,:,1)')); set(gca,'ydir','reverse'); shading interp; colorbar; title('A/B'); colormap(usa2); caxis([-1 +1]*10)
        disp('ret to continue'); pause
        disp(' ')
	
      else
        jj = jj + 1;
      end    %% if
    end      %% while iFound < 0
    if iFound < 0
      fprintf(1,'  %s is not in dir1 \n',zname0);
    end
    
  end        %% if loop . and ..
end          %% loop over ii

disp(' ')
fprintf(1,'dir0 has %2i files while dir1 has %2i files ... found %2i names in common ... recall 3 are . and .. and refprof_trace400 \n',length(dir0),length(dir1),iCnt)
