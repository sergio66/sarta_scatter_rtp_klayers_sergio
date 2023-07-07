figure(9); clf
%% daytime : saved sarta data has 47 profs x 12 solarangles * 6 viewangles
%% dayime : my fits           has 48 profs x 12 solarangles * 6 viewangles

dx = 1 : 0.05 : 2;    plot(dx,histc(weightD,dx)/length(weightD),'bo-',dx,histc(sarta_debug.weight,dx)/length(sarta_debug.weight),'r.-');                  title('weight');    hl = legend('sergio','scott','location','best'); 
  disp('ret to continue'); pause
dx = 0 : 0.01 : 0.15; plot(dx,histc(deltaradD,dx)/length(weightD),'bo-',dx,histc(sarta_debug.f./sarta_debug.weight,dx)/length(sarta_debug.weight),'r.-'); title('delta rad'); hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 1 : 0.05 : 2;    plot(dx,histc(weightD.*pred1D,dx)/length(weightD),'bo-',dx,histc(sarta_debug.term1,dx)/length(sarta_debug.weight),'r.-');           title('Pred1');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 0.05 : 2;    plot(dx,histc(weightD.*pred2D,dx)/length(weightD),'bo-',dx,histc(sarta_debug.term2,dx)/length(sarta_debug.weight),'r.-');           title('Pred2');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 0.05 : 2;    plot(dx,histc(weightD.*pred3D,dx)/length(weightD),'bo-',dx,histc(sarta_debug.term3,dx)/length(sarta_debug.weight),'r.-');           title('Pred3');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 0.05 : 4;    plot(dx,histc(weightD.*pred4D,dx)/length(weightD),'bo-',dx,histc(sarta_debug.term4,dx)/length(sarta_debug.weight),'r.-');           title('Pred4');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 5 : 500;     plot(dx,histc(weightD.*pred5D,dx)/length(weightD),'bo-',dx,histc(sarta_debug.term5,dx)/length(sarta_debug.weight),'r.-');           title('Pred5');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 0.05 : 2;    plot(dx,histc(weightD.*pred6D,dx)/length(weightD),'bo-',dx,histc(sarta_debug.term6,dx)/length(sarta_debug.weight),'r.-');           title('Pred6');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 200 : 2 : 240;   plot(dx,histc(T_highD,dx)/length(weightD),'bo-',dx,histc(chris2.thbi,dx)/length(sarta_debug.weight),'r.-');                         title('Thigh');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause

len = length(deltaradD);
plot(1:len,deltaradD' - sarta_debug.f./sarta_debug.weight); title('\delta deltarad'); disp('ret to continue'); pause
plot(1:len,weightD' -   sarta_debug.weight);                title('\delta weight');   disp('ret to continue'); pause
plot(1:len,weightD'.*pred1D' -   sarta_debug.term1);        title('\delta term1');    disp('ret to continue'); pause
plot(1:len,weightD'.*pred2D' -   sarta_debug.term2);        title('\delta term2');    disp('ret to continue'); pause
plot(1:len,weightD'.*pred3D' -   sarta_debug.term3);        title('\delta term3');    disp('ret to continue'); pause
plot(1:len,weightD'.*pred4D' -   sarta_debug.term4);        title('\delta term4');    disp('ret to continue'); pause
plot(1:len,weightD'.*pred5D' -   sarta_debug.term5);        title('\delta term5');    disp('ret to continue'); pause
plot(1:len,weightD'.*pred6D' -   sarta_debug.term6);        title('\delta term6');    disp('ret to continue'); pause
plot(1:len,T_highD'          -   chris2.thbi);              title('\delta Thigh');    disp('ret to continue'); pause

junkpreds = [sarta_debug.term1; sarta_debug.term2; sarta_debug.term3; sarta_debug.term4; sarta_debug.term5; sarta_debug.term6]';
junkdeltarads = sarta_debug.f';
junkcoef = junkpreds \ junkdeltarads;

predsDX = [weightD.*pred1D weightD.*pred2D weightD.*pred3D weightD.*pred4D weightD.*pred5D weightD.*pred6D];
coefDX0(:,ii) = predsDX \ (weightD.*deltaradD);     %%% predsDX,junkpreds are so close yet computation is so far
coefDX1(:,ii) = junkpreds \ (weightD.*deltaradD);   %%% perfecto
coefDX2(:,ii) = linsolve(predsDX, weightD.*deltaradD);
coefDX3(:,ii) = pinv(predsDX' * predsDX) * (predsDX' * (weightD.*deltaradD));
coefDX4(:,ii) = inv(predsDX' * predsDX) * (predsDX' * (weightD.*deltaradD));

isarta = find(sarta.vchan >= fairs(ii)-0.1,1);
wah = [coefD(:,ii) coefDX0(:,ii) coefDX1(:,ii) coefDX2(:,ii) coefDX3(:,ii) coefDX4(:,ii) sarta.cofn(isarta,1:6)' junkcoef];
disp('Naive WA\Wb  Naive2Wa\Wb   NaiveC\Wb    linsolve      pinv          inv          SARTA    JUNKSARTA')
fprintf(1,'%8.5f     %8.5f     %8.5f     %8.5f     %8.5f     %8.5f     %8.5f    %8.5f \n',wah');
disp(' ')

