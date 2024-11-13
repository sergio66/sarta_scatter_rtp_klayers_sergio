
if std(double(profin.nlevs)) == 0
  figure(1); clf; 
    semilogy(nanmean(profin.ptemp,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.ptemp,2),nanmean(proftop.plevs,2),'r.-',...
             nanmean(profnew.ptemp,2),nanmean(profnew.plevs,2),'k'); 
    set(gca,'ydir','reverse')
    title('T(z)')
  
  figure(2); clf; 
    loglog(nanmean(profin.gas_1,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_1,2),nanmean(proftop.plevs,2),'r.-',...
          nanmean(profnew.gas_1,2),nanmean(profnew.plevs,2),'k'); 
    set(gca,'ydir','reverse')
    title('WV(z)')

  if isfield(proftop,'gas_2')  
    figure(3); clf; 
    loglog(nanmean(profin.gas_2,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_2,2),nanmean(proftop.plevs,2),'r.-',...
           nanmean(profnew.gas_2,2),nanmean(profnew.plevs,2),'k'); 
    set(gca,'ydir','reverse')
    title('CO2(z)')
  end

  if isfield(proftop,'gas_3')  
    figure(4); clf; 
    loglog(nanmean(profin.gas_3,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_3,2),nanmean(proftop.plevs,2),'r.-',...
           nanmean(profnew.gas_3,2),nanmean(profnew.plevs,2),'k'); 
    set(gca,'ydir','reverse')
    title('O3(z)')
  end

  if isfield(proftop,'gas_5')  
    figure(5); clf; 
    loglog(nanmean(profin.gas_5,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_5,2),nanmean(proftop.plevs,2),'r.-',...
           nanmean(profnew.gas_5,2),nanmean(profnew.plevs,2),'k'); 
    set(gca,'ydir','reverse')
    title('CO(z)')
  end

  if isfield(proftop,'gas_6')  
    figure(6); clf; 
    loglog(nanmean(profin.gas_6,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_6,2),nanmean(proftop.plevs,2),'r.-',...
           nanmean(profnew.gas_6,2),nanmean(profnew.plevs,2),'k'); 
    set(gca,'ydir','reverse')
    title('CH4(z)')
  end

  %figure(5); clf; 
  %  loglog(nanmean(profin.gas_34,2),nanmean(profin.plevs,2),'b.-',nanmean(proftop.gas_34,2),nanmean(proftop.plevs,2),'r.-',...
  %         nanmean(profnew.gas_34,2),nanmean(profnew.plevs,2),'k'); 
  %  set(gca,'ydir','reverse')
  %  title('O(z)')

else

  figure(1); clf; 
    for pp = 1 : length(profin.stemp)
      nlevs = profin.nlevs(pp);
      nn =  1:nlevs;
      nnx = 1 : nlevs + length(more_p);
      semilogy(profin.ptemp(nn,pp),profin.plevs(nn,pp),'b.-',proftop.ptemp(:,pp),proftop.plevs(:,pp),'r.-',...
               profnew.ptemp(nnx,pp),profnew.plevs(nnx,pp),'k'); 
      hold on
    end
    set(gca,'ydir','reverse')
    title('T(z)')
    hold off

  figure(2); clf; 
    for pp = 1 : length(profin.stemp)
      nlevs = profin.nlevs(pp);
      nn =  1:nlevs;
      nnx = 1 : nlevs + length(more_p);
      loglog(profin.gas_1(nn,pp),profin.plevs(nn,pp),'b.-',proftop.gas_1(:,pp),proftop.plevs(:,pp),'r.-',...
               profnew.gas_1(nnx,pp),profnew.plevs(nnx,pp),'k'); 
      hold on
    end
    set(gca,'ydir','reverse')
    title('WV(z)')
    hold off

  if isfield(proftop,'gas_2')  
    figure(3); clf; 
    for pp = 1 : length(profin.stemp)
      nlevs = profin.nlevs(pp);
      nn =  1:nlevs;
      nnx = 1 : nlevs + length(more_p);
      loglog(profin.gas_2(nn,pp),profin.plevs(nn,pp),'b.-',proftop.gas_2(:,pp),proftop.plevs(:,pp),'r.-',...
               profnew.gas_2(nnx,pp),profnew.plevs(nnx,pp),'k'); 
      hold on
    end
    set(gca,'ydir','reverse')
    title('CO2(z)')
    hold off
  end

  if isfield(proftop,'gas_3')  
    figure(4); clf; 
    for pp = 1 : length(profin.stemp)
      nlevs = profin.nlevs(pp);
      nn =  1:nlevs;
      nnx = 1 : nlevs + length(more_p);
      loglog(profin.gas_3(nn,pp),profin.plevs(nn,pp),'b.-',proftop.gas_3(:,pp),proftop.plevs(:,pp),'r.-',...
               profnew.gas_3(nnx,pp),profnew.plevs(nnx,pp),'k'); 
      hold on
    end
    set(gca,'ydir','reverse')
    title('O3(z)')
    hold off
  end

  if isfield(proftop,'gas_5')  
    figure(5); clf; 
    for pp = 1 : length(profin.stemp)
      nlevs = profin.nlevs(pp);
      nn =  1:nlevs;
      nnx = 1 : nlevs + length(more_p);
      loglog(profin.gas_5(nn,pp),profin.plevs(nn,pp),'b.-',proftop.gas_5(:,pp),proftop.plevs(:,pp),'r.-',...
               profnew.gas_5(nnx,pp),profnew.plevs(nnx,pp),'k'); 
      hold on
    end
    set(gca,'ydir','reverse')
    title('CO(z)')
    hold off
  end

  if isfield(proftop,'gas_6')  
    figure(6); clf; 
    for pp = 1 : length(profin.stemp)
      nlevs = profin.nlevs(pp);
      nn =  1:nlevs;
      nnx = 1 : nlevs + length(more_p);
      loglog(profin.gas_6(nn,pp),profin.plevs(nn,pp),'b.-',proftop.gas_6(:,pp),proftop.plevs(:,pp),'r.-',...
               profnew.gas_6(nnx,pp),profnew.plevs(nnx,pp),'k'); 
      hold on
    end
    set(gca,'ydir','reverse')
    title('CH4(z)')
    hold off
  end
end  
