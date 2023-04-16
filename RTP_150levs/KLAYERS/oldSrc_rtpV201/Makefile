#
# Makefile for KLAYERS.  Most of the actual work is done by make_klayers.
#

# ------
# Builds
# ------

airs:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV201/klayers_airs_samewater

airs_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV201/klayers_airs_wetwater

oss:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_oss.f cbplev.f
	ln -s incLAY_samewater_oss.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV201/klayers_oss

clean:
	rm -f *.o
