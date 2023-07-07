figure(9); clf
%% daytime : saved sarta data has 47 profs x 12 solarangles * 6 viewangles
%% dayime : my fits           has 48 profs x 12 solarangles * 6 viewangles

dx = 1 : 0.05 : 2;    plot(dx,histc(weightN,dx)/length(weightN),'bo-',dx,histc(sarta_debug.weight,dx)/length(sarta_debug.weight),'r.-');                  title('weight');    hl = legend('sergio','scott','location','best'); 
  disp('ret to continue'); pause
dx = 0 : 0.01 : 0.15; plot(dx,histc(deltaradN,dx)/length(weightN),'bo-',dx,histc(sarta_debug.f./sarta_debug.weight,dx)/length(sarta_debug.weight),'r.-'); title('delta rad'); hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 1 : 0.05 : 2;    plot(dx,histc(weightN.*pred1N,dx)/length(weightN),'bo-',dx,histc(sarta_debug.term1,dx)/length(sarta_debug.weight),'r.-');           title('Pred1');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 0.05 : 2;    plot(dx,histc(weightN.*pred2N,dx)/length(weightN),'bo-',dx,histc(sarta_debug.term2,dx)/length(sarta_debug.weight),'r.-');           title('Pred2');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 0.05 : 2;    plot(dx,histc(weightN.*pred3N,dx)/length(weightN),'bo-',dx,histc(sarta_debug.term3,dx)/length(sarta_debug.weight),'r.-');           title('Pred3');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 0.05 : 4;    plot(dx,histc(weightN.*pred4N,dx)/length(weightN),'bo-',dx,histc(sarta_debug.term4,dx)/length(sarta_debug.weight),'r.-');           title('Pred4');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 5 : 500;     plot(dx,histc(weightN.*pred5N,dx)/length(weightN),'bo-',dx,histc(sarta_debug.term5,dx)/length(sarta_debug.weight),'r.-');           title('Pred5');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 0 : 0.05 : 2;    plot(dx,histc(weightN.*pred6N,dx)/length(weightN),'bo-',dx,histc(sarta_debug.term6,dx)/length(sarta_debug.weight),'r.-');           title('Pred6');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause
dx = 200 : 2 : 240;   plot(dx,histc(T_highN,dx)/length(weightN),'bo-',dx,histc(chris2.thbi,dx)/length(sarta_debug.weight),'r.-');                         title('Thigh');     hl = legend('sergio','scott','location','best');
  disp('ret to continue'); pause

len = length(deltaradN);
plot(1:len,deltaradN' - sarta_debug.f./sarta_debug.weight); title('\delta deltarad'); disp('ret to continue'); pause
plot(1:len,weightN' -   sarta_debug.weight);                title('\delta weight');   disp('ret to continue'); pause
plot(1:len,weightN'.*pred1N' -   sarta_debug.term1);        title('\delta term1');    disp('ret to continue'); pause
plot(1:len,weightN'.*pred2N' -   sarta_debug.term2);        title('\delta term2');    disp('ret to continue'); pause
plot(1:len,weightN'.*pred3N' -   sarta_debug.term3);        title('\delta term3');    disp('ret to continue'); pause
plot(1:len,weightN'.*pred4N' -   sarta_debug.term4);        title('\delta term4');    disp('ret to continue'); pause
plot(1:len,weightN'.*pred5N' -   sarta_debug.term5);        title('\delta term5');    disp('ret to continue'); pause
plot(1:len,weightN'.*pred6N' -   sarta_debug.term6);        title('\delta term6');    disp('ret to continue'); pause
plot(1:len,T_highN'          -   chris2.thbi);              title('\delta Thigh');    disp('ret to continue'); pause

junkpreds = [sarta_debug.term1; sarta_debug.term2; sarta_debug.term3; sarta_debug.term4; sarta_debug.term5; sarta_debug.term6]';
junkdeltarads = sarta_debug.f';
junkcoef = junkpreds \ junkdeltarads;

predsNX = [weightN.*pred1N weightN.*pred2N weightN.*pred3N weightN.*pred4N weightN.*pred5N weightN.*pred6N];
coefNX0(:,ii) = predsNX \ (weightN.*deltaradN);     %%% predsNX,junkpreds are so close yet computation is so far
coefNX1(:,ii) = junkpreds \ (weightN.*deltaradN);   %%% perfecto
coefNX2(:,ii) = linsolve(predsNX, weightN.*deltaradN);
coefNX3(:,ii) = pinv(predsNX' * predsNX) * (predsNX' * (weightN.*deltaradN));
coefNX4(:,ii) = inv(predsNX' * predsNX) * (predsNX' * (weightN.*deltaradN));

isarta = find(sarta.vchan >= fairs(ii)-0.1,1);
wah = [coefN(:,ii) coefNX0(:,ii) coefNX1(:,ii) coefNX2(:,ii) coefNX3(:,ii) coefNX4(:,ii) sarta.cofn(isarta,1:6)' junkcoef];
disp('Naive WA\Wb  Naive2Wa\Wb   NaiveC\Wb    linsolve      pinv          inv          SARTA    JUNKSARTA')
fprintf(1,'%8.5f     %8.5f     %8.5f     %8.5f     %8.5f     %8.5f     %8.5f    %8.5f \n',wah');
disp(' ')

