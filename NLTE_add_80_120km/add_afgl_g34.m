function [y,yp] = add_afgl_g34(x,p,afgl34);

%% this adds on AFGL pressures to SABER (x(p),p) for gasID using iAtm = 1,2,3,4,5,6 for TRP,MLS,MLW,SAS,SAW,STD

[Y,I] = sort(p);
p = p(I);
x = x(I);

iAtm = 6;
gasID = 34;
%afgl34 = quick_read_afgl(gasID,iAtm);

if gasID > 7
  p0 = afgl34.pstd;
  x0 = afgl34.qstd;
else
  p0 = afgl34.piAtm;
  x0 = afgl34.qiAtm;
end
[Y,I] = sort(p0);
p0 = p0(I);
x0 = x0(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REGION 1 : from GND to 100 mb
moo = find(p0 > 100);
y(moo) = 1e-6;

%% REGION 2 : from 100 mb to max(SABER p))
moo = find(p0 > max(p) & p0 < 100);
x1 = interp1(log(p0),log(x0),log(max(p))); x1 = exp(x1); 
if length(moo) > 0
  offset  = x(end)/x1;
  slope = (offset - 1)/(100 - max(p));
  scale(moo) = 1 + (100 - p0(moo)).*slope;
  y(moo) = x0(moo)' .*scale(moo);
end

%% REGION 3 : SABER DATA
moo = find(p0 >= min(p) & p0 <= max(p)); 
y(moo) = interp1(log(p),log(x),log(p0(moo)));
y(moo) = exp(y(moo));

%% REGION 4 : from top of SABER data to 0.000025 mb TOA
moo = find(p0 < min(p)); 
x2 = interp1(log(p0),log(x0),log(min(p))); x2 = exp(x2); 
if length(moo) > 0
  offset  = x(1) - x2;
  y(moo) = x0(moo) + offset; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yp = p0;
loglog(x,p,'bx-',y,yp,'r'); set(gca,'ydir','reverse'); grid
