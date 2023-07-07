%{
-rw-r--r-- 1 sergio pi_strow  2375 Jul  1  2015 ../../getcoef.f
-rw-r--r-- 1 sergio pi_strow 12689 Mar 30  2006 ../../opac2.f
-rwxr-xr-x 1 sergio pi_strow  5280 Mar 29  2006 ../../mwtran.f
-rw-rw-r-- 1 sergio pi_strow  2124 May 17  2001 ../../tb11.for
-rw-rw-r-- 1 sergio pi_strow  2429 May 16  2001 ../../tb10.for
-rwxr-xr-x 1 sergio pi_strow  7378 Apr 25  1997 ../../shval3.for
-rwxr-xr-x 1 sergio pi_strow  2018 Apr 24  1997 ../../bfield.for
-rwxr-xr-x 1 sergio pi_strow  1532 Nov  1  1995 ../../vlint.for
%}

thedir = dir('../../*.f*');
for ii = 1 : length(thedir)

  fname0 = thedir(ii).name;
  fprintf(1,'diffing %s \n',fname0)
  
  f1 = ['../../' fname0];

  if fname0(end) == 'r'
    f2 = [fname0(1:end-2)];
  else
    f2 = fname0;
  end
  
  differ = ['!diff ' f1 ' ' f2];
  eval(differ);
  disp('ret to cont');
  pause
end
