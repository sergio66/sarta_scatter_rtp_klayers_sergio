function [plev]=xaerilayers(pstart, nlay);

% function [plev]=xaerilayers(pstart, nlay);
%
% Calculate AERI pressure layer boundary levels
%
% Input:
%    pstart = start pressure (ie pressure at AERI location) (mb)
%    nlay   = desired number of layers (integer)
%
% Output:
%    plev = [(nlay+1) x 1] AERI layer pressure boundary levels (mb)
%

% Created by Scott Hannon, 6 July 2001 based on aerilayers.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pend=0.005;

nx=nlay/100;

ix=round( [ 1, nx*33, nx*85, nlay] );  % layer index
iy=[ 6, 280,  1300, 6500]/nx;  % approx desired layer thickness (meters)

% Do linear interpolation in log(thickness)
xf=1:nlay;
liy=log(iy);
lyf=interp1(ix,liy,xf,'linear');
yf=exp(lyf);

%semilogy(xf,yf);

% Approximaton P = P_init * exp( -Z / H )
% where
%    P is pressure in millibars
%    Z is altitude in meters
%    H is the scale height in meters
% This can be solved for layer thickness,
%    dz = H * ln( P1 / P2 )

dzsum=cumsum(yf);
dztot=dzsum(nlay);
h=dztot/log(pstart/pend);
%disp(['scale height h=' num2str(h)])

plev=[pstart, pstart*exp(-dzsum/h)]';

%plot(plevs,'bo-'),grid

%%% end of function %%%
