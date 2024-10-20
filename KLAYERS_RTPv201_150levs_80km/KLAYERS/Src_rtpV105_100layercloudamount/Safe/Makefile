#
# Makefile for KLAYERS.  Most of the actual work is done by make_klayers.
#

# ------
# Builds
# ------

airs_v5:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v5
	mv a.out ../Bin/klayers_airs_v5_testme

airs_v4:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_airs_v4_testme

oss:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_oss.f cbplev.f
	ln -s incLAY_oss.f incLAY.f
	make -f make_klayers_v5
	mv a.out ../Bin/klayers_oss

aeri1033:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri1033.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_aeri1033

aeri1013:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri1013.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_aeri1013

aeri999:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri999.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_aeri999

aeri993:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri993.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_aeri993

aeri980:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri980.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_aeri980

aeri961:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri961.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_aeri961

aeri693:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri693.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_aeri693

aeri680:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri680.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_aeri680

aeri667:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri667.f cbplev.f
	ln -s incLAY_airs.f incLAY.f
	make -f make_klayers_v4
	mv a.out ../Bin/klayers_aeri667

clean:
	rm -f *.o
