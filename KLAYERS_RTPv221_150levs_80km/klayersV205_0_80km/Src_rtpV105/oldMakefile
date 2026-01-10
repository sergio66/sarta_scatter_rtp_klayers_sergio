#
# Makefile for KLAYERS.  Most of the actual work is done by make_klayers.
#

# ------
# Builds
# ------

airs_trace:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers_trace
	mv a.out ../Bin/klayers_airs_trace

airs_trace_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers_trace
	mv a.out ../Bin/klayers_airs_trace_wetwater

910to1030_trace_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_910to1030.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers_trace
	mv a.out ../Bin/klayers_910to1030_trace_wetwater

airs:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_airs

oss_trace:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_oss.f cbplev.f
	ln -s incLAY_samewater_oss.f incLAY.f
	make -f make_klayers_trace
	mv a.out ../Bin/klayers_oss_trace

airs_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_airs_wetwater

airs_old:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s old_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_airs_old

aeri1033:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri1033.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_aeri1033

aeri1013:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri1013.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_aeri1013

aeri999:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri999.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_aeri999

aeri993:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri993.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_aeri993

aeri980:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri980.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_aeri980

aeri961:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri961.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_aeri961

aeri693:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri693.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_aeri693

aeri680:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri680.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_aeri680

aeri667:
	rm -f cbplev.f wexsvp.f incLAY.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_aeri667.f cbplev.f
	ln -s incLAY_samewater.f incLAY.f
	make -f make_klayers
	mv a.out ../Bin/klayers_aeri667

clean:
	rm -f *.o
