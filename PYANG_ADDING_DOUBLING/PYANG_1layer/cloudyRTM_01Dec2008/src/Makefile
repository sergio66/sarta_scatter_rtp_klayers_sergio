#
# Makefile for SARTA v1.08 PCLSAM.  Most of the work is done by make_sarta.
#

# ------
# Builds
# ------

m130x_iceaggr_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m130x/' incFTC_apr08x.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam_4500
	mv a.out ../Bin/sarta_apr08_m130x_iceaggr_waterdrop_desertdust_slabcloud_hg3

m140x_iceaggr_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140x/' incFTC_apr08x.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam_4500
	mv a.out ../Bin/sarta_apr08_m140x_iceaggr_waterdrop_desertdust_slabcloud_hg3

m130_iceaggr_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m130/' incFTC_apr08.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam_4500
	mv a.out ../Bin/sarta_apr08_m130_iceaggr_waterdrop_desertdust_slabcloud_hg3

m140_iceaggr_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_apr08.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam_4500
	mv a.out ../Bin/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3

m140_iceaggr_waterdrop_cirrus:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m140/' incFTC_apr08.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_cirrus.f fnmie.f
	make -f make_sarta_pclsam_4500
	mv a.out ../Bin/sarta_apr08_m140_iceaggr_waterdrop_cirrus_slabcloud_hg3

m150_iceaggr_waterdrop_desertdust:
	rm -f *.o incFTC.f fnmie.f
	sed 's/XCSNAMEX/m150/' incFTC_apr08.f > incFTC.f
	ln -s fnmie_iceaggr_waterdrop_desertdust.f fnmie.f
	make -f make_sarta_pclsam_4500
	mv a.out ../Bin/sarta_apr08_m150_iceaggr_waterdrop_desertdust_slabcloud_hg3


# ---------------
# File management
# ---------------

clean:
	rm -f *.o
