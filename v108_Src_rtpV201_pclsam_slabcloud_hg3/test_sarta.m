%% copied from /home/sergio/MATLABCODE/JIGOU/CIRRUS_BRYANBAUM
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /asl/matlib/science
addpath /asl/matlib/h4tools

clear all

run_sarta.klayers_code = '/asl/packages/klayersV205/BinV201/klayers_airs';

%% this is default Baran ICE AGGR
run_sarta.sartacloud_code = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
iDebug = +1;

%% this is new Baum scattering General Habit Model
run_sarta.sartacloud_code = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
iDebug = -1;

fip = mktemp('temp.ip.rtp');
fop = mktemp('temp.op.rtp');
frp = mktemp('temp.rp.rtp');

for yy = 2012 : 2012
  for mm = 06 : 06
    for dd = 01 : 01
      fDIR = ['/asl/data/rtprod_airs/' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d') '/'];
      for hh = 00 : 23
        fCLR  = [fDIR 'clr_era.airs_ctr.' num2str(yy) '.' num2str(mm,'%02d') '.' num2str(dd,'%02d') '.'];
        fCLR  = [fCLR num2str(hh,'%02d') '.rtp'];
        fCLD  = [fDIR 'cld_era.airs_ctr.' num2str(yy) '.' num2str(mm,'%02d') '.' num2str(dd,'%02d') '.'];
        fCLD  = [fCLD num2str(hh,'%02d') '.rtp'];
        fOUT  = [fDIR 'baumGHM_era.airs_ctr.' num2str(yy) '.' num2str(mm,'%02d') '.' num2str(dd,'%02d') '.'];
        fOUT  = [fOUT num2str(hh,'%02d') '.rtp'];
        eeR = exist(fCLR);
        eeD = exist(fCLD);
        eeO = exist(fOUT);
        if eeR > 0 & eeD > 0 & eeO == 0
          fprintf(1,'processing %4i/%2i/%2i : H %2i \n',yy,mm,dd,hh)

          [hc,hac,pc,pac] = rtpread(fCLD);
          [hc,hac,pc,pac] = rtpgrow(hc,hac,pc,pac);

          [h,ha,p,pa] = rtpread(fCLR);
          [h,ha,p,pa] = rtpgrow(h,ha,p,pa);

          rtpwrite(fip,h,ha,p,pa);
          klayerser = ['!' run_sarta.klayers_code  ' fin=' fip ' fout=' fop];
          sartaer = ['!' run_sarta.sartacloud_code ' fin=' fop ' fout=' frp];
          eval(klayerser);
          eval(sartaer);
          [h1,ha1,p1,pa1] = rtpread(frp);
          p1x.rcalc = p1.rcalc;
          p1x.rtime = p1.rtime;
          p1x.rlat  = p1.rlat;
          p1x.rlon  = p1.rlon;
          p1x.atrack = p1.atrack;
          p1x.xtrack = p1.xtrack;

          if iDebug > 0
            tclr = rad2bt(h.vchan,p.rcalc);
            tcld0 = rad2bt(h.vchan,pc.rcalc);
            tcldX = rad2bt(h.vchan,p1x.rcalc);

            plot(h.vchan,tclr - tcld0)
            plot(h.vchan,tcldX - tcld0)
            error('oo')
            pause(0.1)
          end  %% if iDebug

          rtpwrite(fOUT,h1,ha1,p1x,pa1);

        end    %% if exist
      end      %% loop hh
    end        %% loop dd
  end          %% loop mm
end            %% loop yy

rmer = ['!/bin/rm ' fip ' ' fop ' ' frp];
eval(rmer);