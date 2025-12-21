echo "making analystic jacs : H2020"

echo "making AIRS analytic jac SARTA"
# make clean; make -f Makefile_F77 airs_jul22_dev
make clean; make -f Makefile airs_cloudy_jan25_H2020

echo "making CRIS analytic jac SARTA"
# make clean; make -f Makefile_F77 cris_hrg4_p2022_jul22
make clean; make -f Makefile crisg4_hires_jan25_H2020_icebaumGHM_waterdrop_desertdust

echo "making IASI analytic jac SARTA"
# make clean; make -f Makefile_F77 iasi_p2022mar23_dev
make clean; make -f Makefile iasi_jan25_H2020_iceaggr_waterdrop_desertdust_wcon_nte_swch4

echo "making CHIRP analytic jac SARTA"
# make clean; make -f Makefile_F77 chirp_feb20
make clean; make -f Makefile chirp_jul22_dev
