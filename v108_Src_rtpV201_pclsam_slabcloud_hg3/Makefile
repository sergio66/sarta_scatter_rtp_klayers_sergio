#
# Makefile for SARTA v1.08 PCLSAM using rtpV201.  Most of the work is
# done by make_sarta_pclsam.
#

# ------
# Builds
# ------
###   1     2         3      4     5       6       7        8      9
###  ID   freq     fixed  wlines  wcon    ozone    CO      CH4   nonLTE
###
### hmm what is this    XCSNAMEX = m130x,m140x,m150x
### incFTC_airs_apr08x_370.f has FNMLT = /asl/rta/sarta_database/Data_AIRS_apr08/Coef/tunmlt_XCSNAMEX.txt
### and you can see those have tuning!!!!!!
###
### incFTC_airs_apr08.f has FNTMLT = /asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_XCSNAMEX.txt eg
###
### incFTC_airs_apr08.f has FNTMLT = /asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_m140.txt
###
### vs
###
### incFTC_airs_apr08.f has FNTMLT = /asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_m140x.txt
###
### diff /asl/rta/sarta_database/Data_AIRS_apr08/Coef/tunmlt_m140.txt /asl/rta/sarta_database/Data_AIRS_apr08/Coef/tunmlt_m140x.txt is BLAMK ZILCH NADA

### basic only has wcon and nlte tuning (columns 5,9)
### incFTC_airs_apr08_wcon_nte.f has FNTMLT = /asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_wcon_nte.txt
### incFTC_airs_apr08_wcon_nte.f has FNTMLT = /asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_wcon_nte.txt
### incFTC_airs_apr08_wcon_nte.f has FNTMLT = /asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_wcon_nte.txt

########################################################################
########################################################################
########################################################################
### AIRS post-Nov2003 , with tuning ###

m140_icebaumGHM_waterdrop_biom:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_biom.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/xsarta_apr08_m140_iceGHMbaum_waterdrop_biom_slabcloud_hg3

## "volz" is volz
m140_icebaumGHM_waterdropMODIS_volzdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdropMODIS_volz.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/xsarta_apr08_m140_iceGHMbaum_waterdropMODIS_volzdust_slabcloud_hg3

## "desertdust" maybe OPAC, not sure!!!
m140_icebaumGHM_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/xsarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3

## "desertdust" maybe OPAC, not sure!!!
m140_icebaumGHM_waterdropMODIS_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdropMODIS_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceGHMbaum_waterdropMODIS_desertdust_slabcloud_hg3

m140_icebaumGHM_waterdrop_desertdust_newdirs:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08_newdirs.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust_newdirs.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_newdirs

m140_iceaggr_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3

m140_iceaggr_waterdrop_andesite:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_andesite.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceaggr_waterdrop_andesite_slabcloud_hg3

m140_iceaggr_waterdrop_basalt:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_basalt.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceaggr_waterdrop_basalt_slabcloud_hg3

########################################################################
########################################################################
########################################################################
### AIRS post-Nov2003 , with NO tuning ###
m140_icebaumGHM_waterdropMODIS_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08_wcon_nte.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdropMODIS_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceGHMbaum_waterdropMODIS_desertdust_slabcloud_hg3_wcon_nte

m140_icebaumGHM_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08_wcon_nte.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_wcon_nte

m140_iceaggr_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08_wcon_nte.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte

m140_iceaggr_waterdrop_desertdust_wcon_nte_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08_wcon_nte_wcon_nte.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_wcon_nte

m140_iceaggr_waterdrop_andesite_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08_wcon_nte.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_andesite.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceaggr_waterdrop_andesite_slabcloud_hg3_wcon_nte

m140_iceaggr_waterdrop_basalt_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_airs_apr08_wcon_nte.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_basalt.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140_iceaggr_waterdrop_basalt_slabcloud_hg3_wcon_nte

########################################################################
########################################################################
########################################################################
### AIRS pre-Nov2003 ###
m140x_iceGHMbaum_waterdropMODIS_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140x/' incFTC_airs_apr08x_370.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdropMODIS_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140x_iceGHMbaum_waterdropMODIS_desertdust_slabcloud_hg3

## this was just testing the above code for blackclouds
m140x_blackcld:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140x/' incFTC_airs_apr08x_370.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdropMODIS_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140x_blackcld

m140x_icebaumGHM_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140x/' incFTC_airs_apr08x_370.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3

m140x_iceaggr_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140x/' incFTC_airs_apr08x_370.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140x_iceaggr_waterdrop_desertdust_slabcloud_hg3

m140x_iceaggr_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140x/' incFTC_airs_apr08x_370_wcon_nte.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m140x_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte

## for JPL EOF, used the pre_Nov2003 SRFs plus tuning
m130_icebaumGHM_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m130/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m130_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3

## this would be pre_Nov2003 SRFs with no tuning
m130_icebaumGHM_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m130/' incFTC_airs_apr08_wcon_nte.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m130_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_wcon_nte

## anticipates even more shifts, and has tuning
m150_icebaumGHM_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m150/' incFTC_airs_apr08.f > incFTC.f
	ln -s fnmie_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_apr08_m150_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3

########################################################################
########################################################################
########################################################################
### IASI ###

iasi_may09_iceaggr_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_iasi_may09_wcon_nte.f incFTC.f
	ln -s fnmie_iasi_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_iasi_may09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte

iasi_may09_iceaggr_waterdrop_desertdust_wcon_nte_swch4:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_iasi_may09_wcon_nte_swch4.f incFTC.f
	ln -s fnmie_iasi_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_swch4_pclsam
	mv a.out ../BinV201/sarta_iasi_may09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte_swch4

iasi_may09_iceaggr_waterdrop_biom_wcon_nte_swch4:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_iasi_may09_wcon_nte_swch4.f incFTC.f
	ln -s fnmie_iasi_iceaggr_waterdrop_biom.f fnmie.f
	make -f make_sarta_swch4_pclsam
	mv a.out ../BinV201/sarta_iasi_may09_iceaggr_waterdrop_biom_slabcloud_hg3_wcon_nte_swch4

iasi_may09_AmNit_AmSul_biom_wcon_nte_swch4:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_iasi_may09_wcon_nte_swch4.f incFTC.f
	ln -s fnmie_iasi_AmNit_AmSul_biom.f fnmie.f
	make -f make_sarta_swch4_pclsam
	mv a.out ../BinV201/sarta_iasi_may09_AmNit_AmSul_biom_slabcloud_hg3_wcon_nte_swch4

########################################################################
########################################################################
########################################################################
### CrIS g4 ###

crisg4_nov09_iceGHMbaum_waterdropMODIS_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_crisg4_nov09_wcon_nte.f incFTC.f
	ln -s fnmie_crisg4_iceGHMbaum_waterdropMODIS_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_crisg4_nov09_iceGHMbaum_waterdropMODIS_desertdust_slabcloud_hg3_wcon_nte

crisg4_nov09_iceaggr_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_crisg4_nov09_wcon_nte.f incFTC.f
	ln -s fnmie_crisg4_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_crisg4_nov09_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte

crisg4_nov09_iceGHMbaum_waterdrop_desertdust_wcon_nte:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_crisg4_nov09_wcon_nte.f incFTC.f
	ln -s fnmie_crisg4_iceGHMbaum_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../BinV201/sarta_crisg4_nov09_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3_wcon_nte


# ---------------
# File management
# ---------------

clean:
	rm -f *.o
