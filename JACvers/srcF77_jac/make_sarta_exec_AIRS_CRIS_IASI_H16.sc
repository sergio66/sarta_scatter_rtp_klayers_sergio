echo "making analystic jacs : H2016"

echo "making AIRS analytic jac SARTA"
make clean; make -f Makefile airs_cloudy_may19

echo "making CRIS analytic jac SARTA"
make clean; make -f Makefile crisg4_hires_dec17_icebaumGHM_waterdrop_desertdust

echo "making IASI analytic jac SARTA"
make clean; make -f Makefile iasi_may09_iceaggr_waterdrop_desertdust_wcon_nte_swch4

