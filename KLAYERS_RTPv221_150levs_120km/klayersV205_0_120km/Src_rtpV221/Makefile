#
# Makefile for KLAYERS.  Most of the actual work is done by make_klayers.
#

# ------
# Builds
# ------

airs_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs_0_120km.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV221/klayers_airs_wetwater140levs_120km

airs_wetwater_test:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs_0_120km.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV221/klayers_airs_wetwater_test140levs_120km

g80_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs_0_120km.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers_g80
	mv a.out ../BinV221/klayers_g80_wetwater140levs_120km

g80_samewater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs_0_120km.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers_g80
	mv a.out ../BinV221/klayers_g80_samewater140levs_120km

airs_samewater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs_0_120km.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV221/klayers_airs_samewater140levs_120km

oss_samewater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_oss.f cbplev.f
	ln -s incLAY_samewater_oss.f incLAY.f
	make -f make_klayers
	mv a.out ../BinV221/klayers_oss_Samewater140levs_120km

clean:
	rm -f *.o
