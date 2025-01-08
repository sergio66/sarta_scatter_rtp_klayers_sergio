#
# Makefile for KLAYERS.  Most of the actual work is done by make_klayers.
#

# ------
# Builds
# ------
pbl_wetwater_test:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs_oco2_pbl_v2.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers_pbl
#	mv a.out ../BinV201/klayers_pbl_wetwater_test

pbl_wetwater_g80_test:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs_oco2_pbl_v2.f cbplev.f
	ln -s incLAY_wetwater_g80.f incLAY.f
	make -f make_klayers_pbl_g80
	mv a.out ../BinV201/klayers_pbl_wetwater_g80_test

airs_wetwater_test:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV201/klayers_airs_wetwater_test

g80_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_wetwater_g80.f incLAY.f
	make -f make_klayers_g80
	mv a.out ../BinV201/klayers_g80_airs_wetwater

g80_samewater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers_g80
	mv a.out ../BinV201/klayers_g80_samewater

airs_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV201/klayers_airs_wetwater

airs_samewater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV201/klayers_airs_samewater

oss_samewater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_oss.f cbplev.f
	ln -s incLAY_samewater_oss.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV201/klayers_oss_Samewater

clean:
	rm -f *.o
