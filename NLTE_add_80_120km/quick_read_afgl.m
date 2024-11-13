function afgl = quick_read_afgl(gasID,iAtm)

%% gasID = 1:42 (main) or 51-81 (xsc)
%% iAtm = 1,2,3,4,5,6 for TRP,MLS,MLW,SAS,SAW,STD

if nargin == 1
  iAtm = 6;
end

if (gasID > 42 & gasID < 51) | gasID > 81
  error('invalid GASID, gasID = 1:42 (main) or 51-81 (xsc)')
end
if iAtm > 6 | iAtm < 1
  error('invalid iAtm : 1,2,3,4,5,6 for TRP,MLS,MLW,SAS,SAW,STD')
end

gid = 1:81;
%%         1     2    3     4    5                      6    7     8     9    10
gname = {'H2O','CO2','O3','N2O','CO',                  'CH4','O2','NO','SO2','NO2',...
         'NH3','HNO3','OH','HF','HCl',                 'HBr','HI','ClO','OCS','H2CO',...
         'HOCl','N2','NCN','CH3Cl','H2O2',             'C2H2','C2H6','PH3','COF2','SF6',...
         'H2S','HCOOH','HO2','O','ClONO2',             'NO+','HOBr','C2H4','CH3OH','CH3Br',...
         'CH3CN','CF4','X','X','X',                    'X','X','X','X','X',...
         'CFC11','CFC12','CFC13','CFC14','CFC21',      'CFC22','CFC113','CFC114','CFC115','CCl4',...
         'ClONO2','N2O5','HNO4','CFC116','CFC123',     'CFC124','CFC141b','CFC142b','CFC225a','CFC225b',...
         'HFC32','HFC134a','HFC143a','HFC152a','C6H6', 'HFC125','HFC134','SF5CF3','PAN','CH3CN',...
         'SF6'};

iQuiet = +1;
if gasID <= 7 & iQuiet < 0
  fprintf(1,'finding US Std and iAtm AFGL profile for gasID = %2i = %s \n',gasID,gname{gasID})
elseif gasID > 7 & iQuiet < 0
  fprintf(1,'finding US Std               profile for gasID = %2i = %s \n',gasID,gname{gasID})
end

iAtm6 = 6; %% US Std

x0 = load('glatm_16Aug2010.dat');

x = x0(2:871,:);
%for lljunk = 1 : 6
%  for cc = 1 : 11
%    ind = (1:10) + (cc-1)*10 + (lljunk-1)*110
%    junk = x(ind,:);
%    junk = junk';
%    junk = junk(:);
%    matr(lljunk,cc,:) = junk;  %% alt(km), p(mb),T(K),den(cm-3),1,2,3,4,5,6,7
%  end
%end

for lljunk = 1 : 6
  for cc = 1 : 11 %% [alt(km), p(mb), T(K), den(cm-3)], [G1,G2,G3,G4,G5,G6,G7]
    ind = (1:10) + (cc-1)*10 + (lljunk-1)*110;
    junk = x(ind,:);
    junk = junk';
    junk = junk(:);
    matr(lljunk,cc,:) = junk;  
%    fprintf(1,'read iAtm = %2i gas %2i \n',lljunk,cc)
  end
end

%[alt(km), p(mb), T(K), den(cm-3)] takes 01-04
%[G1,G2,G3,G4,G5,G6,G7]            takes 05-11
%[G8-G28]                          takes 12-32
lljunk = 6;
for cc = 8+4 : 28+4
  ind = (1:10) + (cc-1)*10 + (lljunk-1)*110;
  junk = x(ind,:);
  junk = junk';
  junk = junk(:);
  matr(lljunk,cc,:) = junk;  %% alt(km), p(mb),T(K),den(cm-3),1,2,3,4,5,6,7
%  fprintf(1,'read iAtm = %2i gas %2i \n',lljunk,cc-4)
end
  
%% main stuff
%>> plot(squeeze(matr(:,1,:))',1:50) %% alt
%>> plot(squeeze(matr(:,2,:))',1:50) %% p
%>> plot(squeeze(matr(:,3,:))',1:50) %% T
%>> plot(squeeze(matr(:,4,:))',1:50) %% density
%% profiles
%>> plot(squeeze(matr(:,5,:))',1:50) %% WV
%>> plot(squeeze(matr(:,6,:))',1:50) %% CO2
%>> plot(squeeze(matr(:,7,:))',1:50) %% O3
  
afgl.pstd   = squeeze(matr(iAtm6,2,:));
afgl.tstd   = squeeze(matr(iAtm6,3,:));
afgl.dstd   = squeeze(matr(iAtm6,4,:));
if gasID <= 28
  afgl.qstd   = squeeze(matr(iAtm6,4+gasID,:));
end

if gasID <= 7
  afgl.piAtm  = squeeze(matr(iAtm,2,:));
  afgl.tiAtm  = squeeze(matr(iAtm,3,:));
  afgl.diAtm  = squeeze(matr(iAtm,4,:));
  afgl.qiAtm  = squeeze(matr(iAtm,4+gasID,:));
  
  %% before "gasstd"   now "qstd"
  loglog(afgl.qstd,afgl.pstd,'bx-',afgl.qiAtm,afgl.piAtm,'ro-'); set(gca,'ydir','reverse'); ylim([1e-6 1000]); grid
   hl = legend('US Std','iAtm','location','best','fontsize',10);

elseif gasID >= 8 & gasID <= 28
  loglog(afgl.qstd,afgl.pstd,'bx-'); set(gca,'ydir','reverse'); ylim([1e-6 1000]); grid
   hl = legend('US Std','location','best','fontsize',10);

elseif gasID >= 28 & gasID <= 81
  y = load('glatm_16Aug2010_minorgases.dat');
  lljunk = 6;

  for cc = 29 : 42
    ind = (1:10) + (cc-29)*10;
    junk = y(ind,:);
    junk = junk';
    junk = junk(:);
    matr29_81(lljunk,cc-29+1,:) = junk;  
    gas29_81(cc-29+1) = cc;
    %fprintf(1,'read iAtm = %2i gas %2i \n',lljunk,cc)
  end
  %whos matr matr29_81 gas29_81

  for cc = 51:81
    ind = (1:10) + (42-29)*10 + (cc-51)*10;
    junk = y(ind,:);
    junk = junk';
    junk = junk(:);
    matr29_81(lljunk,cc-29+1,:) = junk;  
    gas29_81(cc-29+1) = cc;
    %fprintf(1,'read iAtm = %2i gas %2i \n',lljunk,cc)
  end
  %whos matr matr29_81 gas29_81

  [~,iX] = intersect(gas29_81,gasID); 
  afgl.qstd   = squeeze(matr29_81(iAtm6,iX,:));
  %printarray(gas29_81)

  loglog(afgl.qstd,afgl.pstd,'bx-'); set(gca,'ydir','reverse'); ylim([1e-6 1000]); grid
   hl = legend('US Std','location','best','fontsize',10);

end
