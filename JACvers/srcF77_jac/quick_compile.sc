ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include calpar.f
ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include sunpar.f
ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include ycalt1_od.f
ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include ycalt2_od.f
ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include ycalt3_od.f
ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include ycalt4_od.f
ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include ycalt5_od.f
ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include ycalt6_od.f
ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include ycalt7_od.f
ifort -c -O2 -convert big_endian -extend-source 132 -check bounds -g -traceback -debug -I/asl/packages/rtp/include sarta_pclsam.f

ifort rdprof.o  vaconv.o  calpar.o  qikexp.o  hg3.o  ysetems_pclsam.o  getbot.o ycalt1_od.o  ycalt2_od.o  ycalt3_od.o  ycalt4_od.o  ycalt5_od.o ycalt6_od.o  ycalt7_od.o  bkprep.o  getcld_slab.o faketz_pclsam.o  sunpar.o   rdsun.o  saconv.o  calowp.o  calokw.o  tunmlt.o util.o  opnrtp_pclsam.o  rdrtp_pclsam.o  wrtrtp.o  mean_t.o  calnte.o rdinfo.o  cbplev.o  calrad0.o  calrad1.o  calrad2.o rdcldt.o getmie.o  fnmie.o  rdcoef.o  ccprep_slab.o  intersect.o sarta_pclsam.o -L/asl/packages/rtp/lib -lrtp -L/asl/packages/external/hdf/hdf4/lib -ldf -L/asl/packages/external/jpeg -ljpeg -L/asl/packages/external/zlib -lz

mv a.out ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug
