#
# Makefile for SARTA MW.  Most of the work is done by makefile
# run with eg       make -f Makefile mwrta_atms

# ------
# Builds
# ------

mwrta_amsua:
	rm -f *.o
	rm m_rta.f
	ln -s m_rta1.f m_rta.f
	make -f makefile
	mv ../BinV201/mrta ../BinV201/mwsarta_amsua

mwrta_amsub:
	rm -f *.o
	rm m_rta.f	
	ln -s m_rta2.f m_rta.f	
	make -f makefile
	mv ../BinV201/mrta ../BinV201/mwsarta_amsub

mwrta_atms:
	rm -f *.o
	rm m_rta.f	
	ln -s m_rta3.f m_rta.f	
	make -f makefile
	mv ../BinV201/mrta ../BinV201/mwsarta_atms

