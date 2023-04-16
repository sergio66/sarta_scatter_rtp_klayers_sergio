#
# Makefile for SARTA.  Most of the work is done by make_sarta.
#

# ------
# Builds
# ------

modis:
	rm -f *.o incFTC.f
	ln -s incFTC_modis.f incFTC.f
	make -f make_sarta_modis
	mv a.out ../Bin/sarta_modis

wcononly:
	rm -f *.o incFTC.f
	ln -s incFTC_modis_wcononly.f incFTC.f
	make -f make_sarta_modis
	mv a.out ../Bin/sarta_modis_wcononly

notuning:
	rm -f *.o incFTC.f
	ln -s incFTC_modis_notuning.f incFTC.f
	make -f make_sarta_modis
	mv a.out ../Bin/sarta_modis_notuning

wcononly_test:
	rm -f *.o incFTC.f
	ln -s incFTC_modis_wcononly.f incFTC.f
	make -f make_sarta_modis_test
	mv a.out ../Bin/sarta_modis_wcononly_test

# ---------------
# File management
# ---------------

clean:
	rm -f *.o
