# Makefile for SARTA.  Most of the work is done by make_sarta.

#
# Makefile for SARTA v1.08.  Most of the work is either by
# make_sarta_pclsam (sarta_cloudy)
# or by
# make_sarta (sarta clear only)
#
# Default is sarta_cloudy named m140_icebaumGHM_waterdrop_desertdust

#----------------------------------------------------------------------------------------
# Builds cloudy
#----------------------------------------------------------------------------------------
# AIRS for L1C Jan 2025, H2020
# YAYAYAYAYAYAYAYAYA
airs_cloudy_jan25_H2020:
	rm -f *.o incFTC.f fnmie.f cbplev.f 
	ln -s incFTC_airs_jan25_H2020_pclsam.f incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	ln -s cbplev_L2100.f cbplev.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_airs_l1c_2834_cloudy_jan25_H2020

# AIRS for L1C Feb 2025, H2020, PBL
# YAYAYAYAYAYAYAYAYA
airs_cloudy_feb25_H2020_PBL:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	ln -s incFTC_airs_feb25_H2020_PBL_pclsam.f incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	ln -s cbplev_PBL100.f cbplev.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_airs_l1c_2834_cloudy_feb25_H2020_PBL

# AIRS for L1C May 2019
# YAYAYAYAYAYAYAYAYA
airs_cloudy_may19:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	ln -s incFTC_airs_may19_pclsam.f incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	ln -s cbplev_L2100.f cbplev.f
	make -f make_sarta_pclsam
#	mv a.out ../bin/jac_airs_l1c_2834_cloudy_may19_prod_debug
	mv a.out ../bin/jac_airs_l1c_2834_cloudy_may19_prod

# AIRS for L1C May 2019 w/tuning
airs_cloudy_may19_tun:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	ln -s incFTC_airs_may19_pclsam_tun_fix.f incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	ln -s cbplev_L2100.f cbplev.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_airs_l1c_2834_cloudy_may19_tun_fix_v1

# AIRS
m140_icebaumGHM_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	ln -s cbplev_L2100.f cbplev.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_sarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3X

########################################################################
# CrisHiRes YAYAYAYAYAYAYAYAYA
crisg4_hires_jan25_H2020_icebaumGHM_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	ln -s incFTC_crisg4_jan25_H2020_new.f incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust_cris_hires_new.f fnmie.f
	ln -s cbplev_L2100.f cbplev.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_crisg4_hires_jan25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new

# CrisHiRes YAYAYAYAYAYAYAYAYA
crisg4_hires_feb25_H2020_icebaumGHM_waterdrop_desertdust_PBL:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	ln -s incFTC_crisg4_feb25_H2020_PBL_new.f incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust_cris_hires_new.f fnmie.f
	ln -s cbplev_PBL100.f cbplev.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_crisg4_hires_feb25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new_PBL

crisg4_hires_dec17_icebaumGHM_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	ln -s incFTC_crisg4_oct16_new.f incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust_cris_hires_new.f fnmie.f
	ln -s cbplev_L2100.f cbplev.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new

crisg4_hires_dec17_icebaumGHM_waterdrop_biosmoke:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	ln -s incFTC_crisg4_oct16_new.f incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust_cris_hires_new.f fnmie.f
	ln -s cbplev_L2100.f cbplev.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_cris_fsr_dec17_iceGHMbaum_wdrop_biosmoke

#----------------------------------------------------------------------------------------
# Builds clear
#----------------------------------------------------------------------------------------

########## AIRS ################
apr08_m140_wcon_nte_2016tuning:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	sed 's/XCSNAMEX/m140/g' incFTC_airs_apr08_template_m130_m140_m150_wcon_nte.f > incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_sarta_apr08_m140_wcon_nte_2016tuning

apr08_m140_pge_v6_tunmlt:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	sed 's/XCSNAMEX/m140/g' incFTC_airs_apr08_template_m130_m140_m150_pge_v6_tunmlt.f > incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_apr08_m140_pge_v6_tunmlt

# AIRS L1C
# this for testing new build with jpl tuning
airs_may19_jpl_tunmlt:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_airs_may19_jpl_tunmlt.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_airs_l1c_2834_may19_jpl_tunmlt

#  build with ASL tuning
airs_may19_asl_tunmlt:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_airs_may19_asl_tun_wlin.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_airs_l1c_2834_may19_tun_wlin_0524a

# -------------------------------------------------------- #
# this one for production testing
airs_may19:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_airs_may19_prod_tra.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_airs_l1c_2834_may19_prod_v2

airs_may19_tune:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_airs_may19_prod_tune.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_airs_l1c_2834_may19_prod_tune

# test and dev versions
airs_basic:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_airs_prod2019_dec2018_basic_debug.f  incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_airs_2834_mar19_basic_optr_co2_debug

airs_thrm:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_airs_prod2019_dec2018_basic_optr_tra_nte_thrm.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_airs_2834_mar19_basic_optr_tra_nte_thrm

airs_rd:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_airs_prod2019_dec2018_basic_optr_tra_nte_thrm.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_airs_2834_mar19_rdown

#  ------------  CrIS   ------------ #######
# CrIS (low res) standard resolution
crisg4_nov09_wcon_nte:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_crisg4_nov09_wcon_nte_new.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_crisg4_nov09_wcon_nte_new

# CrIS (low res) standard resolution
cris_lrg4_p2021_dev:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_cris_lrg4_p2021_may21_dev.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_cris_lrg4_p2021_dev

# 
# CrIS extended resolution (high res)
crisg4_oct16:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_crisg4_oct16_new.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_crisg4_oct16_aug20

# CrIS  (high res) basic coef set[1-7] and wcon.
crisg4_mar18:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_cris_hr_g4_mar18_basic_optr_co2_so2_nte.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_crisg4_mar18_basic_so2_x

# CrIS (high Res) with var HDO
crisg4_may18_hdo:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_cris_hr_g4_may18_tra_hdo.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_crisg4_may18_hdo_d6

# CrIS FSR v2.01 (prod2019 dec2018)
cris_hrg4_p2019_dec18:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_cris_hrg4_p2019_dec2018_prod.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_cris_hrg4_p2019dec18_set3

######  -------  IASI  -----------
# this is what I kinda used in the SingleFootprint retrievals
# copied from ../../v108_Src_rtpV201_pclsam_slabcloud_hg3/Makefile
# YAYAYAYAYAYAYAYAYA
iasi_jan25_H2020_iceaggr_waterdrop_desertdust_wcon_nte_swch4:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_iasi_jan25_H2020_wcon_nte_swch4.f incFTC.f
	ln -s fnmie_iasi_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_sarta_iasi_jan25_H2020_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4

iasi_may09_iceaggr_waterdrop_desertdust_wcon_nte_swch4:
	rm -f *.o incFTC.f fnmie.f cbplev.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_iasi_may09_wcon_nte_swch4.f incFTC.f
	ln -s fnmie_iasi_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../bin/jac_sarta_iasi_may09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4

# Legacy
iasi_may09:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_iasi_may09_wcon_nte_new.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_sarta_iasi_may09_wcon_nte_new

# IASI basic 7-set w/continuum
iasi_mar2019:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_iasi_dec2018_all_clh.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_iasi_mar2019_test

iasi_jun19_test:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_iasi_jun19_test.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_iasi_jun19_test

# IASI w/trace gases and updated refl thermal (HDO not computed)
iasi_jun19:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_iasi_jun19_xnte.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_iasi_jun19_xnte

# --------------------- CHIRP -------------------
# CHIRP Bare bones
chirp_feb20:
	rm -f *.o incFTC.f
	ln -s cbplev_L2100.f cbplev.f
	ln -s incFTC_chirp_base_tra_thrm_nte.f incFTC.f
	make -f make_sarta
	mv a.out ../bin/jac_chirp_feb20_base_tra_thrm_nte

#----------------------------------------------------------------------------------------
# File management
#----------------------------------------------------------------------------------------

clean:
	rm -f *.o
	rm -f *_genmod.f90 *_genmod.mod
