echo "making analystic jacs : H2020"

echo "making PBL AIRS analytic jac SARTA"
# make clean; make -f Makefile_F77 airs_jul22_dev
# make clean; make -f Makefile airs_cloudy_jan25_H2020
make clean; make -f Makefile airs_cloudy_feb25_H2020_PBL

echo "making PBL CRIS analytic jac SARTA"
# make clean; make -f Makefile_F77 cris_hrg4_p2022_jul22
# make clean; make -f Makefile crisg4_hires_jan25_H2020_icebaumGHM_waterdrop_desertdust
make clean; make -f Makefile crisg4_hires_feb25_H2020_icebaumGHM_waterdrop_desertdust_PBL

# echo "making PBL IASI analytic jac SARTA"
# make clean; make -f Makefile_F77 iasi_p2022mar23_dev
# make clean; make -f Makefile iasi_jan25_H2020_iceaggr_waterdrop_desertdust_wcon_nte_swch4
