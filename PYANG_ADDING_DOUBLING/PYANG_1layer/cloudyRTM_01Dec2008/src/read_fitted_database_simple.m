%% based on read_fitted_database.f
fname = '../cloudcoeffs/fitedcoerc_simple.dat';
fid = fopen(fname,'r');

Nwl = 201; %% number of wavelength pts
Ntcld=50;
MM=4;

for ii = 1 : Nwl
  fread(fid,1,'integer*4');
  wav(ii) = fread(fid, 1, 'real*4');
  fread(fid,1,'integer*4');
  for jj = 1 : Ntcld
    fread(fid,1,'integer*4');
    tau(jj) = fread(fid, 1, 'real*4');
    fread(fid,1,'integer*4');

    fread(fid,1,'integer*4');
    [blah,count] = fread(fid, [1,MM*MM] ,       'real*4');
    fread(fid,1,'integer*4');

    end
  end
fclose(fid);

plot(wav)