#
# RTP top level Makefile for rtpV201
#
# this file sets some common parameters and then calls makefiles
# in the src, test, and utils subdirectories
#
# parameters should normally be set here in preference to setting
# them in the subdirectories.  In most cases, all you need to set
# are HDFHOME, RTPHOME, and your local compiler options
#

# --------------
# HDF parameters
# --------------

#V1
# set HDFHOME to point to the local HDF installation
# set HDFHOME to point to the local HDF installation
#HDFHOME = /asl/opt
# HDF libraries
#HDFLIB   = -L$(HDFHOME)/lib -lmfhdf -ldf -ljpeg -lz

#########################
#V2
# set HDFHOME to point to the local HDF installation
##HDFHOME1 = /asl/opt
#HDFHOME1 = /asl/opt/hdf/4.2r4
## HDF libraries
#HDFLIB1   = -L$(HDFHOME1)/lib -lmfhdf -ldf -ljpeg -lz
#
## more libs
#HDFHOME2 = /asl/opt/szip/2.1
##HDFLIB2   = -L$(HDFHOME2)/lib -lsz
##HDFLIB2   = -L$(HDFHOME2)/lib -lsz -libsz.so.2
#HDFLIB2   = -L$(HDFHOME2)/lib -lsz -L/usr/lib64/libsz.so.2
#
#HDFLIB = $(HDFLIB1) $(HDFLIB2)

################## more recent Jan 2016
HDF_INSTALL1 = /usr/cluster/contrib/hdf4/4.2.11
HDFLIB1   = -L$(HDF_INSTALL1)/lib -ldf -ljpeg -lz       ## -lmfhdf

# more libs
HDF_INSTALL2 = /usr/cluster/contrib/szip2/szip2.1
HDF_INSTALL2 = /usr/ebuild/software/Szip/2.1.1-GCCcore-6.4.0
HDFLIB2   = -L$(HDF_INSTALL2)/lib -lsz -L/usr/lib64/libsz.so.2

HDFLIB = $(HDFLIB1) $(HDFLIB2)

# --------------
# RTP parameters
# --------------

# set RTPHOME to point to the local RTP installation
RTPHOME = ..

# ------------------
# C compiler options
# -------------------

# -64 for 64-bit IRIX
# CFLAGS = -g -64
# CFLAGS = -g
CFLAGS = -O
CC = cc

# ------------------------
# Fortran compiler options
# ------------------------

# Absoft Fortran
# --------------
# -N109  fold all names to upper case
# -C     check array bounds
# -O     some optimizations
# -N3    add record info to unformatted files
# -s     static allocation
# FFLAGS = -N109 -C -O -N3 -s
# FFLAGS = -g
#FFLAGS = -O
#FLIB = -lU77	# Linux Absoft needs -lU77
#F77 = /asl/opt/absoft/absoft10.0/bin/af77

# SunOS options
# -------------
# FFLAGS = -e -fast -w     
# F77 = f77

# SGI options
# -----------
# include -cpp option if not default 
# FFLAGS = -O
# -64 for 64-bit IRIX
# FFLAGS  = -O -64
# F77 = f77

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
FFLAGS = -O -convert big_endian -extend-source 132

# gfortran
# ------------------------------------------------------------------------
#    -C           array bounds checking
#    -g           include symbol table info
#    -O           moderate optimization
F77 = /usr/ebuild/software/GCCcore/11.2.0/bin/gfortran
# FLIB =
FFLAGS = -O -fconvert=big-endian -ffixed-line-length-none -fdec-structure

########################################################################
########################################################################
########################################################################

# pass the variables set above to the subdirectory makefiles
#
EXPORTS = HDFHOME="$(HDFHOME)" RTPHOME="$(RTPHOME)" \
	CC="$(CC)" CFLAGS="$(CFLAGS)" CLIB="$(CLIB)" \
	F77="$(F77)" FFLAGS="$(FFLAGS)" FLIB="$(FLIB)"

# -------------
# Make Targets
# -------------

all: SRC TEST UTILS

SRC:
	cd src && make $(EXPORTS)
TEST:
	cd test && make $(EXPORTS)
UTILS:
	cd utils && make $(EXPORTS)

clean:
	cd src && make clean
	cd test && make clean
	cd utils && make clean

# -------------------------
# make an RTP distribution
# -------------------------
#
# "make dist" makes a distribution named by the current working
# directory.  For example, if we are in the subdirectory rtpV201
# "make dist" will clean things up and then create an rtpV201.tar 
# in the parent directory that unpacks to rtpV201/<etc>.
# 
dist: clean
	rm rtp.tar bin/* lib/* 2> /dev/null || true
	rbase=`/bin/pwd`                && \
	    rbase=`basename $${rbase}`  && \
	    cd ..                       && \
	    tar -cf $${rbase}.tar          \
		$${rbase}/bin              \
		$${rbase}/doc              \
		$${rbase}/include          \
		$${rbase}/lib              \
		$${rbase}/Makefile         \
		$${rbase}/README           \
		$${rbase}/src              \
		$${rbase}/test             \
		$${rbase}/utils
	@echo created `/bin/pwd`.tar
	@echo "\"make all\" to rebuild the local distribution"

