make clean; make -f Makefile pbl_wetwater_test
make clean; make -f Makefile pbl_wetwater_g80_test
make clean; make -f Makefile airs_wetwater_test
make clean; make -f Makefile airs_wetwater
make clean; make -f Makefile g80_wetwater

#make clean;  make -f make_klayers_pbl_g80; mv a.out ../BinV201/klayers_pbl_wetwater_g80_test
#make clean;  make -f make_klayers_g80    ; mv a.out ../BinV201/klayers_airs_wetwater_g80_test
#make clean;  make -f make_klayers_pbl    ; mv a.out ../BinV201/klayers_pbl_wetwater_test
#make clean;  make -f make_klayers        ; mv a.out ../BinV201/klayers_airs_wetwater_test

../BinV201/klayers_pbl_wetwater_g80_test  fin=../Data/adafgl_16Aug2010_ip.rtp fout=../Data/adafgl_03Jan2025_pbl_g80.op.rtp   nwantp=-1
../BinV201/klayers_airs_wetwater_g80_test fin=../Data/adafgl_16Aug2010_ip.rtp fout=../Data/adafgl_03Jan2025_airs_g80.op.rtp  nwantp=-1
../BinV201/klayers_pbl_wetwater_test      fin=../Data/adafgl_16Aug2010_ip.rtp fout=../Data/adafgl_03Jan2025_pbl.op.rtp
../BinV201/klayers_airs_wetwater_test     fin=../Data/adafgl_16Aug2010_ip.rtp fout=../Data/adafgl_03Jan2025_airs.op.rtp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/
/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_PBL/Jan02_2025/klayersV205/BinV201/klayers_pbl_wetwater_g80_test  fin=sub6000clr.ip.rtp fout=sergio_sub6000clr_pbl_g80.op.rtp   nwantp=-1 >& ugh
/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_PBL/Jan02_2025/klayersV205/BinV201/klayers_airs_wetwater_g80_test fin=sub6000clr.ip.rtp fout=sergio_sub6000clr_airs_g80.op.rtp  nwantp=-1 >& ugh
/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_PBL/Jan02_2025/klayersV205/BinV201/klayers_pbl_wetwater_test      fin=sub6000clr.ip.rtp fout=sergio_sub6000clr_pbl.op.rtp        >& ugh
/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_PBL/Jan02_2025/klayersV205/BinV201/klayers_airs_wetwater_test     fin=sub6000clr.ip.rtp fout=sergio_sub6000clr_airs.op.rtp       >& ugh
cd /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/KLAYERS_PBL/Jan02_2025/klayersV205/Src_rtpV201/

