function [ierr] = wrtcoef_nte(fname, id, freq, coef);

% function [ierr] = wrtcoef_nte(fname, ichan, freq, coef);
%
% Write out nonLTE coefficients
%
% Input:
%    fname : {string} name of coef file to create
%    id    : [n x 1] channel ID
%    freq  : [n x 1] channel frequency
%    coef  : [n x m] coefficients
%
% Output
%    ierr  : [1 x 1] error flag, 0=OK, 1=error
%

% Created: 03 October 2008, Scott Hannon - based on wrt_nte.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize ierr to bad
ierr = 1;

% Check input
%
if (nargin ~= 4)
   disp('insufficient arguments')
   return
end
%
d = size(id);
if (length(d) ~= 2 | min(d) ~= 1)
   disp('id must be a 1-D vector');
   return
end
nchan = max(d);
%
d = size(freq);
if (length(d) ~= 2 | min(d) ~= 1)
   disp('freq must be a 1-D vector');
   return
end
if (max(d) ~= nchan)
   disp('id and freq must be the same length')
   return
end
%
d = size(coef);
if (length(d) ~= 2 | d(1) ~= nchan)
   disp('coef must be a 2-D matrix with lead dimension n same as id & freq')
   return
end
ncoef = d(2);


% Open output file
fprintf(1,'writing to %s \n',fname);
fid = fopen(fname,'w','ieee-be');

% FORTRAN record marker info
ifm = 4*(1 + 1 + ncoef); % 4 bytes each * (id, freq, coef)

% write out data to fortran file
for ic=1:nchan
   fwrite(fid,ifm,'integer*4');
   fwrite(fid,id(ic),'integer*4');
   fwrite(fid,freq(ic),'real*4');
   fwrite(fid,coef(ic,:),'real*4');
   fwrite(fid,ifm,'integer*4');
end
fclose(fid);

disp(['finished writing data for ' int2str(nchan) ' channels to file ' fname])
ierr = 0;

%%% end of program %%%
