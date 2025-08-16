%% see ~/git/MATLABCODE_Git/AveragingKernel_DOF/driver_vertical_resolution.m

%for ii=1:9
%  figure(ii); clf; colormap jet
%end
  
iP = input('Enter Profile : ');
if length(iP) == 0
  iP = 1;
end
  
iAIRS_or_CrIS = +1;

clear JACS Z
i100_L2_or_PBL = +1;
[Z1,P1,JACS1,p1,fout1,rout1,dz1,dp1] = comparejacs_AIRSL2_vs_PBL(iAIRS_or_CrIS,i100_L2_or_PBL,iP);

clear JACS Z
i100_L2_or_PBL = -1;
[Z2,P2,JACS2,p2,fout2,rout2,dz2,dp2] = comparejacs_AIRSL2_vs_PBL(iAIRS_or_CrIS,i100_L2_or_PBL,iP);

t1 = rad2bt(fout1,rout1);
t2 = rad2bt(fout2,rout2);

ii = 100-length(P1)+1:100; Z1 = Z1(ii); dz1 = dz1(ii);
ii = 100-length(P2)+1:100; Z2 = Z2(ii); dz2 = dz2(ii);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fout1;
figure(1); clf
  haha = tiledlayout(2,1);
  h1 = nexttile; plot(fout,t1,fout,t2);
  h2 = nexttile; plot(fout,t1-t2);

figure(1); clf
    yyaxis left;
    ax = gca; % Get current axes handle
    ax.YAxis(1).Color = 'r'; % Set left Y-axis color to red
    yyaxis right;
    ax.YAxis(2).Color = 'b'; % Set right Y-axis color to blue
    yyaxis left;
    plot(fout,t1,fout,t2); % Plots will use the first color ('r')
    ylabel('BT [K]')
    yyaxis right;
    plot(fout,t1-t2); % Plots will use the second color ('b')
    ylabel('BT [K]')
    grid on;
    xlim([645 1625]); xlabel('Wavenumber [cm-1]')
    
figure(1); clf
colororder({'b', 'r'}); % Set default colors for left and right axes
    yyaxis left;
    plot(fout,t1,fout,t2); % Plots will use the first color ('r')
    ylabel('BT [K]')    
    yyaxis right;
    plot(fout,t1-t2); % Plots will use the second color ('b')
    ylabel('dBT [K]')
    grid on;
    xlim([645 1625]); xlabel('Wavenumber [cm-1]')    
    
figure(2); clf; plot(1:length(Z1),Z1,1:length(Z2),Z2)
  ylabel('Z [km]')
  xlabel('index')
  
figure(3); clf; semilogy(Z1(1:length(P1)),P1,Z2(1:length(P2)),P2); set(gca,'ydir','reverse');
  ylim([0.1 1000]); xlabel('Z(km)');  ylabel('P(mb)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%JAC = [jacWV jacO3 jacTZ jacWGT jacST];  %% WV O3 T WGT SURF
nlays1 = length(P1);
nlays2 = length(P2);

ii = 1; ind = (1:nlays1) + (ii-1)*nlays1; jWV1 = JACS1(:,ind);
ii = 2; ind = (1:nlays1) + (ii-1)*nlays1; jOZ1 = JACS1(:,ind);
ii = 3; ind = (1:nlays1) + (ii-1)*nlays1; jTZ1 = JACS1(:,ind);
ii = 4; ind = (1:nlays1) + (ii-1)*nlays1; jWG1 = JACS1(:,ind);

ii = 1; ind = (1:nlays2) + (ii-1)*nlays2; jWV2 = JACS2(:,ind);
ii = 2; ind = (1:nlays2) + (ii-1)*nlays2; jOZ2 = JACS2(:,ind);
ii = 3; ind = (1:nlays2) + (ii-1)*nlays2; jTZ2 = JACS2(:,ind);
ii = 4; ind = (1:nlays2) + (ii-1)*nlays2; jWG2 = JACS2(:,ind);

dzz1 = ones(2378,1) * dz1';
zjWV1 = jWV1 ./ dzz1;
zjOZ1 = jOZ1 ./ dzz1;
zjTZ1 = jTZ1 ./ dzz1;
zjWG1 = jWG1 ./ dzz1;

dzz2 = ones(2378,1) * dz2';
zjWV2 = jWV2 ./ dzz2;
zjOZ2 = jOZ2 ./ dzz2;
zjTZ2 = jTZ2 ./ dzz2;
zjWG2 = jWG2 ./ dzz2;

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); clf; pcolor(fout,P1,jWV1'); title('jac wv1 AIRSL2'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');
figure(5); clf; pcolor(fout,P2,jWV2'); title('jac wv2 PBL');    colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');

figure(4); clf; pcolor(fout,Z1,jWV1'); title('jac wv1 AIRSL2'); colorbar; shading interp
figure(5); clf; pcolor(fout,Z2,jWV2'); title('jac wv2 PBL'); colorbar; shading interp

figure(6); clf; %plot(fout,sum(jWV1,2),fout,sum(jWV2,2))
colororder({'b', 'r'}); % Set default colors for left and right axes
    yyaxis left;
    plot(fout,sum(jWV1,2),fout,sum(jWV2,2)); % Plots will use the first color ('r')
    ylabel('sum jacWV')
    yyaxis right;
    plot(fout,sum(jWV1,2)./sum(jWV2,2)-1); % Plots will use the second color ('b')
    ylim([-1 +1]*5)    
    ylabel('ratio jacWV1/jacWV2')
    grid on;
    xlim([645 1625]); xlabel('Wavenumber [cm-1]')

%%%

figure(7); clf; pcolor(fout,P1,zjWV1'); title('zjac wv1 AIRSL2'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');
figure(8); clf; pcolor(fout,P2,zjWV2'); title('zjac wv2 PBL');    colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');

figure(7); clf; pcolor(fout,Z1,zjWV1'); title('zjac wv1 AIRSL2'); colorbar; shading interp
figure(8); clf; pcolor(fout,Z2,zjWV2'); title('zjac wv2 PBL'); colorbar; shading interp

figure(9); clf; %plot(fout,sum(zjWV1,2),fout,sum(zjWV2,2))
colororder({'b', 'r'}); % Set default colors for left and right axes
    yyaxis left;
    plot(fout,sum(zjWV1,2),fout,sum(zjWV2,2)); % Plots will use the first color ('r')
    ylabel('sum zjacWV')
    yyaxis right;
    plot(fout,sum(zjWV1,2)./sum(zjWV2,2)-1); % Plots will use the second color ('b')
    ylim([-1 +1]*5)
    ylabel('ratio zjacWV1/zjacWV2')
    grid on;
    xlim([645 1625]); xlabel('Wavenumber [cm-1]')

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); clf; pcolor(fout,P1,jOZ1'); title('jac oz1 AIRSL2'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');
figure(5); clf; pcolor(fout,P2,jOZ2'); title('jac oz2 PBL');    colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');

figure(4); clf; pcolor(fout,Z1,jOZ1'); title('jac oz1 AIRSL2'); colorbar; shading interp
figure(5); clf; pcolor(fout,Z2,jOZ2'); title('jac oz2 PBL'); colorbar; shading interp

figure(6); clf; %plot(fout,sum(jOZ1,2),fout,sum(jOZ2,2))
colororder({'b', 'r'}); % Set default colors for left and right axes
    yyaxis left;
    plot(fout,sum(jOZ1,2),fout,sum(jOZ2,2)); % Plots will use the first color ('r')
    ylabel('sum jacOZ')
    yyaxis right;
    plot(fout,sum(jOZ1,2)./sum(jOZ2,2)-1); % Plots will use the second color ('b')
    ylim([-1 +1]*5)    
    ylabel('ratio jacOZ1/jacOZ2-1')
    grid on;
    xlim([645 1625]); xlabel('Wavenumber [cm-1]')

%%%

figure(7); clf; pcolor(fout,P1,zjOZ1'); title('zjac oz1 AIRSL2'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');
figure(8); clf; pcolor(fout,P2,zjOZ2'); title('zjac oz2 PBL');    colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');

figure(7); clf; pcolor(fout,Z1,zjOZ1'); title('zjac oz1 AIRSL2'); colorbar; shading interp
figure(8); clf; pcolor(fout,Z2,zjOZ2'); title('zjac oz2 PBL'); colorbar; shading interp

figure(9); clf; %plot(fout,sum(zjOZ1,2),fout,sum(zjOZ2,2))
colororder({'b', 'r'}); % Set default colors for left and right axes
    yyaxis left;
    plot(fout,sum(zjOZ1,2),fout,sum(zjOZ2,2)); % Plots will use the first color ('r')
    ylabel('sum zjacOZ')
    yyaxis right;
    plot(fout,sum(zjOZ1,2)./sum(zjOZ2,2)-1); % Plots will use the second color ('b')
    ylim([-1 +1]*5)
    ylabel('ratio zjacOZ1/zjacOZ2-1')
    grid on;
    xlim([645 1625]); xlabel('Wavenumber [cm-1]')

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); clf; pcolor(fout,P1,jTZ1'); title('jac tz1 AIRSL2'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');
figure(5); clf; pcolor(fout,P2,jTZ2'); title('jac tz2 PBL');    colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');

figure(4); clf; pcolor(fout,Z1,jTZ1'); title('jac tz1 AIRSL2'); colorbar; shading interp
figure(5); clf; pcolor(fout,Z2,jTZ2'); title('jac tz2 PBL'); colorbar; shading interp

figure(6); clf; %plot(fout,sum(jTZ1,2),fout,sum(jTZ2,2))
colororder({'b', 'r'}); % Set default colors for left and right axes
    yyaxis left;
    plot(fout,sum(jTZ1,2),fout,sum(jTZ2,2)); % Plots will use the first color ('r')
    ylabel('sum jacTZ')
    yyaxis right;
    plot(fout,sum(jTZ1,2)./sum(jTZ2,2)-1); % Plots will use the second color ('b')
    ylim([-1 +1]*5)    
    ylabel('ratio jacTZ1/jacTZ2')
    grid on;
    xlim([645 1625]); xlabel('Wavenumber [cm-1]')

%%%

figure(7); clf; pcolor(fout,P1,zjTZ1'); title('zjac tz1 AIRSL2'); colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');
figure(8); clf; pcolor(fout,P2,zjTZ2'); title('zjac tz2 PBL');    colorbar; shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log');

figure(7); clf; pcolor(fout,Z1,zjTZ1'); title('zjac tz1 AIRSL2'); colorbar; shading interp
figure(8); clf; pcolor(fout,Z2,zjTZ2'); title('zjac tz2 PBL'); colorbar; shading interp

figure(9); clf; %plot(fout,sum(zjTZ1,2),fout,sum(zjTZ2,2))
colororder({'b', 'r'}); % Set default colors for left and right axes
    yyaxis left;
    plot(fout,sum(zjTZ1,2),fout,sum(zjTZ2,2)); % Plots will use the first color ('r')
    ylabel('sum zjacTZ')
    yyaxis right;
    plot(fout,sum(zjTZ1,2)./sum(zjTZ2,2)-1); % Plots will use the second color ('b')
    ylim([-1 +1]*5)
    ylabel('ratio zjacTZ1/zjacTZ2-1')
    grid on;
    xlim([645 1625]); xlabel('Wavenumber [cm-1]')

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

