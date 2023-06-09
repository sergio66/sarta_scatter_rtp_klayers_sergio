#
# RTP library makefile
#
# this makefile builds the library librtp.a for the RTP C and 
# Fortran application interfaces
#
# the site-specific parameters HDFHOME, RTPHOME, and the C and 
# Fortran compiler options can be set in the top level makefile

# --------------
# HDF parameters
# --------------

# orig
# set HDFHOME to point to the local HDF installation
#HDFHOME = /asl/opt
# HDF definitions, relative to HDFHOME
#HDFINC   = -I$(HDFHOME)/include
#HDFLIB   = -L$(HDFHOME)/lib -ldf -ljpeg -lz

# new
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
# RTPHOME = /usr/local/rtp
RTPHOME = ..

# RTP definitions, relative to RTPHOME
RTPINC  = -I$(RTPHOME)/include		# RTP C and F77 includes
RTPLIB  = -L$(RTPHOME)/lib -lrtp	# RTP C and F77 library

# ------------------
# C compiler options
# -------------------

CFLAGS = -O $(RTPINC)
CC = cc

# ------------
# make targets
# -------------

all: prep librtp.a install

prep:
	cp -p ../include/rtp.h .
	cp -p ../include/pvdefs.h .

# --------------------------
# vdata interface procedures
# --------------------------

PVSRC = pvopen.c pvclose.c \
	pvread1.c pvwrite1.c pvread2.c pvwrite2.c pvread3.c pvwrite3.c 

PVOBS = pvopen.o pvclose.o \
	pvread1.o pvwrite1.o pvread2.o pvwrite2.o pvread3.o pvwrite3.o 

$(PVOBS): $(PVSRC) rtp.h pvdefs.h
	$(CC) $(CFLAGS) -c $(PVSRC) $(HDFINC)

# ----------------
# RTP C procedures
# ----------------

RTPSRC = rtpwrite1.c rtpwrite2.c rtpread1.c rtpread2.c getchan.c \
		rtpclose1.c gsplit.c rtphcheck.c headinit.c profinit.c \
		fmatch.c hsize.c
RTPOBS = rtpwrite1.o rtpwrite2.o rtpread1.o rtpread2.o getchan.o \
		rtpclose1.o gsplit.o rtphcheck.o headinit.o profinit.o \
		fmatch.o hsize.o

$(RTPOBS): $(RTPSRC) rtp.h pvdefs.h
	$(CC) $(CFLAGS) -c $(RTPSRC) $(HDFINC)

# -----------------------
# test support procedures
# -----------------------

TSTSRC = dump_flist.c dump_vbuf.c dump_pstr.c dump_hstr.c \
		dump_attrs.c dump_chan.c

TSTOBS = dump_flist.o dump_vbuf.o dump_pstr.o dump_hstr.o \
		dump_attrs.o dump_chan.o

$(TSTOBS): $(TSTSRC) rtp.h pvdefs.h
	$(CC) $(CFLAGS) -c $(TSTSRC) $(HDFINC)

# -------------
# C/Fortran API
# -------------

FORSRC = rtpopen.c rtpread.c rtpwrite.c rtpclose.c rtpvers.c rtpinit.c
FOROBS = rtpopen.o rtpread.o rtpwrite.o rtpclose.o rtpvers.o rtpinit.o

$(FOROBS): $(FORSRC) rtp.h pvdefs.h
	$(CC) $(CFLAGS) -c $(FORSRC) $(HDFINC)

# -----------
# RTP library
# -----------

librtp.a: $(PVOBS) $(RTPOBS) $(FOROBS) $(TSTOBS)
	ar -rcs librtp.a $(PVOBS) $(RTPOBS) $(FOROBS) $(TSTOBS)

# -------
# install
# -------

install:librtp.a
	cp -a librtp.a $(RTPHOME)/lib

# -------
# cleanup
# -------

clean:
	rm *.o librtp.a 2> /dev/null || true

