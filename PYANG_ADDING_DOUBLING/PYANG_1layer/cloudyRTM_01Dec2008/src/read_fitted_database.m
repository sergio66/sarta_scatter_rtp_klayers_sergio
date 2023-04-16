%% based on read_fitted_database.f
fname = '../cloudcoeffs/fitedcoerc.dat';
fid = fopen(fname,'r','ieee-le');

%% do not worry about record info
version = fread(fid, 50, 'char');
version = setstr(version')
%str = fread(fid,50,'uint8=>char'); str'  %% junk str

Nwl = 201; %% number of wavelength pts
Ntcld=50;
MM=4;

for ii = 1 : Nwl
  C = textscan(fid,'%f %s'); C{1}
  wav(ii) = C{1}; str = C{2};
  for jj = 1 : Ntcld
    C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f')
    end
  end
fclose(fid)
