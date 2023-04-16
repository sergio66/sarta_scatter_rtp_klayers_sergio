#
# Makefile for SARTA v1.08 PCLSAM using rtpV201.  Most of the work is
# done by make_sarta_pclsam.
#

# ------
# Builds
# ------

### AIRS pre-Nov2003 ###

m140x_icebaumGHM_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140x/' incFTC_airs_apr08x_370.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_100layerNEW

m140x_iceaggr_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140x/' incFTC_airs_apr08x_370.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140x_iceaggr_waterdrop_desertdust_slabcloud_hg3_100layerNEW

m140x_iceaggr_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140x/' incFTC_airs_apr08x_370_wcon_nte.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140x_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_100layerNEW


### AIRS post-Nov2003 ###

m140_iceaggr_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_100layerNEW

m140_iceaggr_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08_wcon_nte.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_100layerNEW


### IASI ###

iasi_may09_iceaggr_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_iasi_may09_wcon_nte.f incFTC.f
	ln -s fnmie_iasi_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_iasi_may09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_100layerNEW

iasi_may09_iceaggr_waterdrop_desertdust_wcon_nte_swch4:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_iasi_may09_wcon_nte_swch4.f incFTC.f
	ln -s fnmie_iasi_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_swch4_pclsam
	mv a.out ../BinV201/sarta_iasi_may09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4_100layerNEW

iasi_may09_iceaggr_waterdrop_biom_wcon_nte_swch4:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_iasi_may09_wcon_nte_swch4.f incFTC.f
	ln -s fnmie_iasi_iceaggr_waterdrop_biom.f fnmie.f
	make -f make_sarta_swch4_pclsam
	mv a.out ../BinV201/sarta_iasi_may09_iceaggr_waterdrop_biom_slabcloud_hg3_wcon_nte_swch4_100layerNEW

iasi_may09_AmNit_AmSul_biom_wcon_nte_swch4:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_iasi_may09_wcon_nte_swch4.f incFTC.f
	ln -s fnmie_iasi_AmNit_AmSul_biom.f fnmie.f
	make -f make_sarta_swch4_pclsam
	mv a.out ../BinV201/sarta_iasi_may09_AmNit_AmSul_biom_slabcloud_hg3_wcon_nte_swch4_100layerNEW

### CrIS g4 ###

crisg4_nov09_iceaggr_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_crisg4_nov09_wcon_nte.f incFTC.f
	ln -s fnmie_crisg4_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_crisg4_nov09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_100layerNEW

# ---------------
# File management
# ---------------

clean:
	rm -f *.o
