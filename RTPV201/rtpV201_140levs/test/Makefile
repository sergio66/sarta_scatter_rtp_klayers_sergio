#
# RTP test Makefile
#
# this file builds assorted C and Fortran tests
#
# NOTE: these tests are intended for debugging and quick checks,
# and as simple examples of using the RTP libraries.  They are not
# validation tests; that was done by testing the interoperation of
# these routines with the Matlab RTP implementation.
#
# the site-specific parameters HDFHOME, RTPHOME, and the C and
# Fortran compiler options can be set in the top level makefile

# --------------
# HDF parameters
# --------------

## orig
# set HDFHOME to point to the local HDF installation
#HDFHOME = /asl/opt
# HDF definitions, relative to HDFHOME
#HDFINC   = -I$(HDFHOME)/include
#HDFLIB   = -L$(HDFHOME)/lib -ldf -ljpeg -lz

################## more recent Jan 2016
HDF_INSTALL1 = /usr/cluster/contrib/hdf4/4.2.11
HDFLIB1   = -L$(HDF_INSTALL1)/lib -ldf -ljpeg -lz       ## -lmfhdf

# more libs
HDF_INSTALL2 = /usr/cluster/contrib/szip2/szip2.1
HDF_INSTALL2 = /usr/ebuild/software/Szip/2.1.1-GCCcore-6.4.0
HDFLIB2   = -L$(HDF_INSTALL2)/lib -lsz -L/usr/lib64/libsz.so.2

HDFINC   = -I$(HDF_INSTALL1)/include
HDFLIB = $(HDFLIB1) $(HDFLIB2)

# --------------
# RTP parameters
# --------------

# set RTPHOME to point to the local RTP installation
RTPHOME = ..

# RTP definitions, relative to RTPHOME
RTPINC  = -I$(RTPHOME)/include		# RTP C and F77 includes
RTPLIB  = -L$(RTPHOME)/lib -lrtp	# RTP C and F77 library
SRCINC  = -I$(RTPHOME)/src		# low-level test includes

# ------------------
# C compiler options
# -------------------

CFLAGS = -O
CC = cc

# ------------------------
# Fortran compiler options
# ------------------------

#FFLAGS = -O
#FLIB = -lU77
#F77 = /asl/opt/absoft/absoft10.0/bin/af77

# INTEL ifort F77
# ------------------------------------------------------------------------
# -u  : turn off implicit typing of variables
# -g  : generate debugging information (turns off optimization)
# -C  : do run time subscript range checking
# -w0 : inform about unused variables
# -O3 : heavy optimization
# -64 : 64-bit objects (libraries must match)
# -us : puts trailing underscore
F77 = ifort
#FLIB =
OPT = -O -convert big_endian -extend-source 132

# gfortran
# ------------------------------------------------------------------------
#    -g           include symbol table info
#    -O           moderate optimization
F77 = /usr/ebuild/software/GCCcore/11.2.0/bin/gfortran
# FLIB =
OPT = -O2 -fconvert=big-endian -ffixed-line-length-none -fdec-structure

# ----------------
# Main target list
# ----------------

all: pvtest rtptest1 rtptest2 rtpsize ftest1 ftest2

# -------------------------
# rtpdump, pv and RTP tests
# -------------------------

# "pv" test (low level test)
pvtest: pvtest.o
	$(CC) $(CFLAGS) -o pvtest pvtest.o $(RTPLIB) $(HDFLIB)
pvtest.o: pvtest.c
	$(CC) $(CFLAGS) -c pvtest.c $(SRCINC) $(HDFINC)

# RTP main test 1
rtptest1: rtptest1.o
	$(CC) $(CFLAGS) -o rtptest1 rtptest1.o $(RTPLIB) $(HDFLIB)
rtptest1.o: rtptest1.c
	$(CC) $(CFLAGS) -c rtptest1.c $(RTPINC) $(HDFINC)

# RTP main test 2
rtptest2: rtptest2.o
	$(CC) $(CFLAGS) -o rtptest2 rtptest2.o $(RTPLIB) $(HDFLIB)
rtptest2.o: rtptest2.c
	$(CC) $(CFLAGS) -c rtptest2.c $(RTPINC) $(HDFINC)

# rtp size checks
rtpsize: rtpsize.o 
	$(CC) $(CFLAGS) -o rtpsize rtpsize.o $(RTPLIB) $(HDFLIB)
rtpsize.o: rtpsize.c
	$(CC) $(CFLAGS) -c rtpsize.c $(RTPINC) $(HDFINC)

# -------------
# Fortran tests
# -------------

ftest1: ftest1.o
	$(F77) $(FFLAGS) -o ftest1 ftest1.o $(RTPLIB) $(HDFLIB) 
ftest1.o: ftest1.f
	$(F77) $(FFLAGS) -c ftest1.f $(RTPINC)

ftest2: ftest2.o
	$(F77) $(FFLAGS) -o ftest2 ftest2.o $(RTPLIB) $(HDFLIB) 
ftest2.o: ftest2.f
	$(F77) $(FFLAGS) -c ftest2.f $(RTPINC)

# -------
# cleanup
# -------

clean:
	rm *.o rtpdump rtpsize \
		pvtest pvtest.hdf \
		rtptest? rtptest?.hdf \
		ftest? ftest?.hdf \
		2> /dev/null || true
