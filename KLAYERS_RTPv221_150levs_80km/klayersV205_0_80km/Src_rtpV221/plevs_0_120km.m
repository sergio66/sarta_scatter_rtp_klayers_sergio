%% IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, 
%% VOL. 41, NO. 2, FEBRUARY 2003 303A
%% An Overview of the AIRS Radiative Transfer Model
%% L. Larrabee Strow, Scott E. Hannon, Sergio De Souza-Machado, Howard E. Motteler, Member, IEEE, andDavid Tobin 

%C. Atmospheric Layering

% The AIRS fast model uses two different layering methods tomodel the
% optical depths. The most intuitive is a grid of verticalslabs with
% constant pressure. In this case, each layer is definedby the two
% bounding grid pressure levels (so layers require alevel grid). The
% discretized version of the radiative transferequation discussed in
% the previous section uses this 100-layerpressure grid and ultimately
% all s must be available on thisgrid. All component gas
% transmittances, with the exception ofwater vapor for 600 channels,
% are parameterized on this gridwhich spans the range 1100–0.005
% hPa. The 101 atmosphericpressure levels which divide the atmosphere
% into 100 layers aredefined as(19) 
%     P(i) = (Ai^2 + Bi +C )^((7/2)
% where is level number, and A,B,C and are constants. By fixing P(i) =
% 1100, P(38)= 300, P(101) = 0.005 Pa, we can then solve for these
% three constants. This re-lation gives us smoothly varying layers and
% is fine enough tonot limit the accuracy of the radiative transfer
% equation.

ii = 1 : 101;
P(1)   = 1100;
P(38)  = 300;
P(101) = 0.005;

%[1 1 1; 38^2 38 1; 101^2 101 1](x) = [1100 300 0.005].^(2/7);
x = inv([1 1 1; 38^2 38 1; 101^2 101 1]) * [1100 300 0.005]'.^(2/7)

A = x(1); B = x(2); C = x(3); 
% x = -1.550789414500297e-04     -5.593654380586063e-02      7.451622227151780e+00

P = (x(1)*ii.*ii + x(2) * ii + x(3)).^(7/2); plot(ii,P)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = 1 : 121;
P(1)   = 1100;
P(30)  = 300;
P(121) = 2.5e-5;

%[1 1 1; 38^2 38 1; 121^2 121 1](x) = [1100 300 2.5e-5].^(2/7);
x = inv([1 1 1; 30^2 30 1; 121^2 121 1]) * [1100 300 2.5e-5]'.^(2/7)

A = x(1); B = x(2); C = x(3); 
% x =      9.142103672821475e-06    -6.234116456565193e-02     7.457862626866450e+00

P = (x(1)*ii.*ii + x(2) * ii + x(3)).^(7/2); plot(ii,P)
semilogy(ii,P,'.-'); line([1 150],[1 1]*5e-3,'color','k','linewidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = 1 : 101;
P(1)   = 1100;
P(20)  = 300;
P(101) = 2.25e-5;

%[1 1 1; 38^2 38 1; 121^2 121 1](x) = [1100 300 2.25e-5].^(2/7);
x = inv([1 1 1; 20^2 20 1; 101^2 101 1]) * [1100 300 2.25e-5]'.^(2/7)

A = x(1); B = x(2); C = x(3); 
% x =    5.831592961719922e-04    -1.329532417106244e-01     7.527900686818923e+00

P = (x(1)*ii.*ii + x(2) * ii + x(3)).^(7/2); plot(ii,P)
semilogy(ii,P,'.-'); line([1 150],[1 1]*5e-3,'color','k','linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = 1;
fprintf(fid,'       DATA  ( PLEV(I), I = 101, 52, -1 ) \n');
  ind = 101 : -1 : 97;  fprintf(fid,'     $         / %10.6f,%10.6f,%10.6f,%10.6f,%10.6f,\n',P(ind));
for jj = 1 : 8
  ind = ind  - 5;  fprintf(fid,'     $           %10.6f,%10.6f,%10.6f,%10.6f,%10.6f,\n',P(ind));
end
  ind = ind  - 5;  fprintf(fid,'     $           %10.6f,%10.6f,%10.6f,%10.6f,%10.6f/ \n',P(ind));

fprintf(fid,'       DATA  ( PLEV(I), I = 51, 1, -1 ) \n');
  ind = 051 : -1 : 47;  fprintf(fid,'     $         / %10.4f,%10.4f,%10.4f,%10.4f,%10.4f,\n',P(ind));
for jj = 1 : 9
  ind = ind  - 5;  fprintf(fid,'     $           %10.4f,%10.4f,%10.4f,%10.4f,%10.4f,\n',P(ind));
end
  ind = 001;       fprintf(fid,'     $           %10.4f /\n',P(ind));
fprintf(fid,'       END \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('../Grid/xcbplev_airs120km.f','w');
fprintf(fid,'       DATA  ( PLEV(I), I = 101, 52, -1 ) \n');
  ind = 101 : -1 : 97;  fprintf(fid,'     $         / %10.6f,%10.6f,%10.6f,%10.6f,%10.6f,\n',P(ind));
for jj = 1 : 8
  ind = ind  - 5;  fprintf(fid,'     $           %10.6f,%10.6f,%10.6f,%10.6f,%10.6f,\n',P(ind));
end
  ind = ind  - 5;  fprintf(fid,'     $           %10.6f,%10.6f,%10.6f,%10.6f,%10.6f/ \n',P(ind));

fprintf(fid,'       DATA  ( PLEV(I), I = 51, 1, -1 ) \n');
  ind = 051 : -1 : 47;  fprintf(fid,'     $         / %10.4f,%10.4f,%10.4f,%10.4f,%10.4f,\n',P(ind));
for jj = 1 : 9
  ind = ind  - 5;  fprintf(fid,'     $           %10.4f,%10.4f,%10.4f,%10.4f,%10.4f,\n',P(ind));
end
  ind = 001;       fprintf(fid,'     $           %10.4f /\n',P(ind));
fprintf(fid,'       END \n');
fclose(fid);
