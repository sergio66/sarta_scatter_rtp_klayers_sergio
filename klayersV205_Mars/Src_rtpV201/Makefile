#
# Makefile for KLAYERS.  Most of the actual work is done by make_klayers.
#

# ------
# Builds
# ------

## MARS default
## in BinV201      klayers_mars -> klayers_mars_wetwater
mars_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f grav.f cbgids.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_mars.f cbplev.f
	ln -s incLAY_wetwater_mars.f incLAY.f
	ln -s grav_mars.f grav.f
	ln -s cbgids_mars.f cbgids.f
	make -f make_klayers
	mv a.out ../BinV201/klayers_mars_wetwater

## AIRS default
## in BinV201      klayers_airs -> klayers_airs_wetwater
airs_wetwater:
	rm -f cbplev.f wexsvp.f incLAY.f grav.f cbgids.f
	ln -s new_wexsvp.f wexsvp.f
	ln -s ../Grid/cbplev_airs.f cbplev.f
	ln -s incLAY_wetwater.f incLAY.f
	ln -s grav_airs.f grav.f
	ln -s cbgids_airs.f cbgids.f
	make -f make_klayers
	mv a.out ../BinV201/klayers_airs_wetwater

clean:
	rm -f *.o
