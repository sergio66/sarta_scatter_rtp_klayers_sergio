#
# Makefile for SARTA v1.07.  Most of the work is done by make_sarta.
#

# ------
# Builds
# ------

iceaggr_waterdrop_volz_1_May07_100layer:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_dec05.f incFTC.f
	ln -s fnmie_iceaggr_waterdrop_volz_1_May07.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../Bin/sarta_dec05_iceaggr_waterdrop_volz_1_May07_100layer_testme

iceaggr_waterdrop_volzMay07_100layer:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_dec05.f incFTC.f
	ln -s fnmie_iceaggr_waterdrop_volzMay07.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../Bin/sarta_dec05_iceaggr_waterdrop_volzMay07_100layer_testme

iceaggr_waterdrop_andesite_100layer:
	rm -f *.o incFTC.f fnmie.f
	ln -s incFTC_dec05.f incFTC.f
	ln -s fnmie_iceaggr_waterdrop_andesite.f fnmie.f
	make -f make_sarta_pclsam
	mv a.out ../Bin/sarta_dec05_iceaggr_waterdrop_andesite_100layer_testme


# ---------------
# File management
# ---------------

clean:
	rm -f *.o
