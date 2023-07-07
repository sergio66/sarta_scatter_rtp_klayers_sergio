addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE

[rad,w] = readkcstd('mars_kcartaX.dat');
[jac,w] = readkcjac('mars_kcartaX.jac');

figure(1); [fc,qc] = quickconvolve(w,rad,0.5,0.5); tc = rad2bt(fc,qc); plot(w,rad2bt(w,rad),'b',fc,tc,'c')

wv  = jac(:,001:100) + jac(:,201:300) + jac(:,301:400) +  jac(:,401:500);
co2 = jac(:,101:200);
T   = jac(:,501:600);
wgt = jac(:,601:700);
st  = jac(:,701);

[fc,jacWV]  = quickconvolve(w,wv,0.5,0.5);
[fc,jacCO2] = quickconvolve(w,co2,0.5,0.5);
[fc,jacT]   = quickconvolve(w,T,0.5,0.5);
[fc,jacWGT] = quickconvolve(w,wgt,0.5,0.5);
[fc,jacST]  = quickconvolve(w,st,0.5,0.5);

%% look at /home/sergio/KCARTA/INCLUDE/MARS_database_params_2021/airslevels.param_mars
marslevels
playsN = plevs(1:100)-plevs(2:101);
playsD = log(plevs(1:100)./plevs(2:101));
plays = playsN ./ playsD;
plays = fliplr(plays);

figure(1); pcolor(fc,plays,jacWV');  set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; colormap jet; xlabel('wavenumber cm-1'); ylabel('p(mb)'); colorbar; title('WV jac');
figure(2); pcolor(fc,plays,jacCO2'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; colormap jet; xlabel('wavenumber cm-1'); ylabel('p(mb)'); colorbar; title('CO2 jac');
figure(3); pcolor(fc,plays,jacT');   set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; colormap jet; xlabel('wavenumber cm-1'); ylabel('p(mb)'); colorbar; title('T jac');
figure(4); pcolor(fc,plays,jacWGT'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; colormap jet; xlabel('wavenumber cm-1'); ylabel('p(mb)'); colorbar; title('wgt fcn');
