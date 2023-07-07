#   This is the makefile for program "mrta"
############################################
### Object files ###
OBJS = bfield.o  dcli.o  dilec4.o  fastem_mod.o  foam.o  fresnel.o \
       getcoef.o  largesc.o  mwtran.o amsutau.o  opac.o  opac2.o  \
       shval3.o  smallsc.o mw_bt.o tb11.o  vlint.o lsurface.o  \
       m_rta.o
# ---------
# HDF flags
# ---------
# HDF_INSTALL : set this to point to the HDF location
# HDFLIB : hdf lib, relative to install path above
#
# linux
HDF_INSTALL = /home/motteler/asl/hdf/HDF4.1r3/NewHDF
HDFLIB   = -L$(HDF_INSTALL)/lib -lmfhdf -ldf -ljpeg -lz
#
# sgi
# HDF_INSTALL = /salsify/users/motteler/hdftest
# HDFLIB   = -L$(HDF_INSTALL)/lib -lmfhdf -ldf -ljpeg -lz

 
# ---------
# RTP flags
# ---------
#
RTPLIB = -L/asl/packages/rtpV105/lib -lrtp

#-----
#Barnet library
#------
CBLib = -L/home/cbarnet/code/lib -lmw_v4 -lairs -lgsfc

### Compiler name ###
# Linux absoft:
F77 = /usr/local/absoft/bin/f77
f77lib = -L/usr/local/absoft/lib -lU77
#
# Generic:
#F77 = f77

### Compiler options ###
# Linux absoft:
#    -f    fold variables to lower case
#    -C    array bounds checking
#    -O    moderate optimization
#    -s    static allocation
#    -g    include symbol table info
#    -N26  read big-endian binary files
#    -N3   use embedded record lengths in binary files 
#OPT = -s -O -N3 -N26 -N109 -f
OPT = -s -O -N3 -N26 -N109
#
# SGI:
#    -O       turn on default optimization (same as -O2)
#    -O0      turn off all optimization
#    -static
#    -check_bounds
# OPT = -O

### Build executable ###
mrta: $(OBJS)

#sgi	$(F77) $(OBJS) -o ../Bin/mrta  $(RTPLIB) $(HDFLIB)
#linux	$(F77) $(OBJS) -o ../Bin/mrta  $(RTPLIB) $(HDFLIB) -lU77
#	$(F77) $(OBJS) -o ../Bin/mrta  $(RTPLIB) $(HDFLIB) $(CBLib) -lU77
	$(F77) $(OBJS) -o ../Bin/mrta  $(RTPLIB) $(HDFLIB) -lU77
	
### Compile, create object files ###
.f.o:
	$(F77) -c $(OPT) $*.f


### Set dependency on include files ###
$(OBJS): ./rtpdefs.f  ./emismw.h  ./amsurta.com ./paramet.com 

### Remove object files ###
clean:
	rm $(OBJS)

### End of makefile for MRTA ###
