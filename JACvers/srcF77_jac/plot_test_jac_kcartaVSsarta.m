usechans = h.ichan;
figure(20); plot(h.vchan - kcrads.fKc(usechans))
figure(20); plot(h.vchan,rad2bt(h.vchan,porig.rcalc),h.vchan,rad2bt(h.vchan,kcrads.rKc(usechans)));
figure(20); plot(h.vchan,rad2bt(h.vchan,porig.rcalc) - rad2bt(h.vchan,kcrads.rKc(usechans))); title('BT diff SARTA-KCARTA \newline (clds are slightly different) ')
disp('ret to continue'); pause

%figure(20); clf; plot(h.vchan,rad2bt(h.vchan,porig.rcalc)-rad2bt(h.vchan,pnew.rcalc),'kx-',h.vchan,jacx.jac(:,7),'b',h.vchan,stempjac_sarta_fast,'r',h.vchan,kcjacs.rKc(usechans,389),'g'); xlim([640 1640]); 
%  hl = legend('BT DIFF orig-new','STEMP JAC finite diff','STEMP JAC analytic','STEMP JAC kcarta','location','best','fontsize',10);
figure(20); clf; plot(h.vchan,jacx.jac(:,7),'b',h.vchan,stempjac_sarta_fast,'r',h.vchan,kcjacs.rKc(usechans,389),'k'); xlim([640 1640]); 
  hl = legend('SARTA finite diff','SARTA analytic','KCARTA','location','best','fontsize',10); title('Comparing STEMP jacs')
xlim([840 1440])
disp('ret to continue'); pause

figure(20); clf; ind = 3; ind = (1:nlays) + (ind-1)*nlays; plot(h.vchan,sum(jacx.tjac,2),'b',h.vchan,sum(ptempjac_sarta_fast,2),'r',h.vchan,sum(kcjacs.rKc(usechans,ind),2),'k');
  title('Comparing column Tjacs'); hl = legend('SARTA Finite','SARTA analytic','KCARTA','location','best','fontsize',10); xlim([640 1640])
xlim([640 840])
disp('ret to continue'); pause

figure(20); clf; ind = 1; ind = (1:nlays) + (ind-1)*nlays; plot(h.vchan,sum(jacx.wvjac,2),'b',h.vchan,sum(g1jac_sarta_fast,2),'r',h.vchan,sum(kcjacs.rKc(usechans,ind),2),'k');
  title('Comparing column WVjacs'); hl = legend('SARTA Finite','SARTA analytic','KCARTA','location','best','fontsize',10); xlim([640 1640]) 
xlim([1240 1440])
disp('ret to continue'); pause

figure(20); clf; ind = 2; ind = (1:nlays) + (ind-1)*nlays; plot(h.vchan,sum(jacx.o3jac,2),'b',h.vchan,sum(g3jac_sarta_fast,2),'r',h.vchan,sum(kcjacs.rKc(usechans,ind),2),'k');
  title('Comparing column O3jacs'); hl = legend('SARTA Finite','SARTA analytic','KCARTA','location','best','fontsize',10); xlim([640 1640])
xlim([640 840])
xlim([940 1140])
disp('ret to continue'); pause

figure(20); clf; ind = 3; ind = (1:nlays) + (ind-1)*nlays; plot(h.vchan,0*sum(jacx.o3jac,2),'b',h.vchan,sum(wgtfcn_sarta_fast,2),'r',h.vchan,sum(kcjacs.rKc(usechans,ind),2),'k');
  title('Comparing WGT FCNs'); hl = legend('SARTA Finite','SARTA analytic','KCARTA','location','best','fontsize',10); xlim([640 1640])
xlim([640 840])
disp('ret to continue'); pause
