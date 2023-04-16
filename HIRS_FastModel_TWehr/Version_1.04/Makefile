# MOLOTOV
F77FLAGS  = -c -convert big_endian -check overflow -check_bounds -O0 -tune host
F77FLAGS2 =    -convert big_endian -check overflow -check_bounds -O0 -tune host
FCOMPILER = f90

# LINUX
F77FLAGS  = -c -O
F77FLAGS2 =    -O
FCOMPILER = g77

# SGI
F77FLAGS  = -c -O
F77FLAGS2 =    -O
FCOMPILER = f77



hffp.x:		hffp_main.o hffp_tau.o hffp_util.o
		$(FCOMPILER) hffp_main.o hffp_tau.o hffp_util.o $(F77FLAGS2) -o hffp.x

hffp_main.o:	hffp_main.f hffp_aux_dec.f hffp_glob_dec.f hffp_init_dec.f hffp_ref_dec.f
		$(FCOMPILER) hffp_main.f $(F77FLAGS) -o hffp_main.o

hffp_tau.o:	hffp_tau.f hffp_aux_dec.f hffp_glob_dec.f hffp_ref_dec.f
		$(FCOMPILER) hffp_tau.f $(F77FLAGS) -o hffp_tau.o

hffp_util.o:	hffp_util.f hffp_aux_dec.f hffp_glob_dec.f hffp_ref_dec.f hffp_rt_dec.f
		$(FCOMPILER) hffp_util.f $(F77FLAGS) -o hffp_util.o

