#   This is the makefile for program "mrta"
############################################
### Object files ###
OBJS = bfield.o  dcli.o  dilec4.o  fastem_mod.o  foam.o  fresnel.o \
       getcoef.o  largesc.o  mwtran.o amsutau.o opac.o opac2.o  \
       shval3.o  smallsc.o mw_bt.o tb11.o  vlint.o lsurface.o  \
       wrtrtp.o m_rta.o

OBJS = bfield.o  dcli.o  dilec4.o  fastem_mod.o  foam.o  fresnel.o \
       getcoef.o  largesc.o  mwtran.o opac2.o  \
       shval3.o  smallsc.o mw_bt.o tb11.o  vlint.o lsurface.o  \
       wrtrtp.o m_rta.o

# ---------
# HDF flags
# ---------
HDF_INSTALL1 = /usr/cluster/contrib/hdf4/4.2.11
HDFHOME1 = -L$(HDF_INSTALL1)/lib -ldf -ljpeg -lz ## -lmfhdf

# more libs
HDF_INSTALL2 = /usr/cluster/contrib/szip2/szip2.1
HDFHOME2   = -L$(HDF_INSTALL2)/lib -lsz -L/usr/lib64/libsz.so.2

HDFLIB = $(HDFHOME1) $(HDFHOME2)
#HDFHOME = $(HDFHOME1)
 
# ---------
# RTP flags
# ---------
#
# set RTPHOME to point to the local RTP installation
#RTPHOME = ../../rtpV201
RTPHOME = /home/sergio/git/rtp/rtpV201

# RTP library
#RTPLIB = -L$(RTPHOME)/lib -lrtp
RTPLIB = -L$(RTPHOME)/lib -lrtp

# RTP include
RTPINC = -I$(RTPHOME)/include

#-----
#compilers
#-----

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

#---
# INTEL
# ifort F77
# ------------
# -u  : turn off implicit typing of variables
# -g  : generate debugging information (turns off optimization)
# -C  : do run time subscript range checking
# -w0 : inform about unused variables
# -O3 : heavy optimization
# -64 : 64-bit objects (libraries must match)
# -us : puts trailing underscore
#
F77 = /cm/shared/apps/intel/composer_xe/2015.3.187/bin/intel64/ifort
F77 = /cm/shared/apps/intel/composer_xe/2015.5.223/bin/intel64/ifort
OPT = -O -convert big_endian -extend-source 132
OPT = -O -check bounds -heap-arrays -g -traceback -debug inline-debug-info -check all -fp-stack-check -convert big_endian -extend-source 132
FLIB = -lblas -mkl -lU77

########################################################################
### Build executable ###
mrta: $(OBJS)
#sgi	$(F77) $(OBJS) -o ../BinV201/mrta  $(RTPLIB) $(HDFLIB)
#linux	$(F77) $(OBJS) -o ../BinV201/mrta  $(RTPLIB) $(HDFLIB) -lU77
	$(F77) $(OBJS) -o ../BinV201/mrta  $(RTPLIB) $(HDFLIB) 
	
### Compile, create object files ###
.f.o:
	$(F77) -c $(OPT) $*.f

### Set dependency on include files ###
$(OBJS): /home/sergio/git/rtp/rtpV201/include/rtpdefs.f  ./emismw.h  ./amsurta.com ./paramet.com 

### Remove object files ###
clean:
	rm $(OBJS)

### End of makefile for MRTA ###
