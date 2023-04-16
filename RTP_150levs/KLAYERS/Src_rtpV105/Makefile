#
# Makefile for KLAYERS.  Most of the actual work is done by make_klayers.
#

# ------
# Builds
# ------

g80_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers_g80
	mv a.out ../Bin/klayers_g80_wetwater

g80_samewater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers_g80
	mv a.out ../Bin/klayers_g80_samewater

airs_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_airs_wetwater

airs_samewater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_airs_samewater

oss_samewater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_oss.f cbplev.f
	ln -s incLAY_samewater_oss.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_oss_samewater

clean:
	rm -f *.o
