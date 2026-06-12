
   ===================================
   README file for the KLAYERS package
   ===================================
   KLAYERS programmer: originally Scott Hannon.
           followed by: Chris Hepplewhite (email: chepplew@umbc.edu)
   UMBC Atmospheric Spectroscopy Laboratory (ASL)
   Previous update: 15 March 2002
   Present update:   2 Nov 2024.

-------------------------------------------------------------------------------
* Introduction:
---------------

This package contains the KLAYERS program; you should also have
the RTP package.  You will not be able to compile KLAYERS until
you have built the RTP library, so do that first.

KLAYERS is a program to convert "level" point profiles into
"layer" integrated slab profiles on a particular layering grid.
The makefile is set up to compile using the 100 AIRS layers grid
unless another grid is specified.  See "Doc/description.txt" for
details on using the program.

KLAYERS has two build versions: the standard and the extended gas 
version referred to as 'g80'. Each version requires it's own RTP
references (refer to the build and make files) and are built 
separately.

The standard version is for use with regular processing of hyperspectral
models with the full spectral bandwidth and limited gases. Approximately
4800 spectral channels and 12 gases can be used in combination. Note that
calling KLAYERS without a wanted gas list, 8 default gases are used.  

The 'g80' extended gas version can accept a maximum of 50 spectral 
channels and the full list of 80 gases (actual number depends on the 
reference atmospheric state vector available). 

For each of the two version types, different atmospheric layering
can be used, the long-standing version used in most models uses the
AIRS pressure layering, and a new layering for planetary boundary
layer focussed research called 'PBL' layering. See associated
definitions in the Grid/ directory. Both use 100 layers.

* How to Build KLAYERS:
------------------------
First decide to build the standard or 'g80' version and update the
symbolic flie links in the RTP library
  ../rtpV221/include/rtpdefs.h, rtp.h.
  ../rtpV221/lib/librtp.a

Be aware that the tests outlined below are designed to work with the 
standard build.

Second: update the Makefile so that the correct Grid/ and
make_klayers_<> are selected. Ensure that the Bin/<exec> file
is not overwritten unless intended.

Third: check and update the appropriate  make_klayers_<> file.

Fourth: Build the executable. When running KLAYERS remember to
use the correct input parameters for the purpose intended.

* Overview of the package Contents:
-----------------------------------
Data/glatm.dat: model profile database (ascii text).

Doc: documentation for KLAYERS.

Grid: layer grids.  The usual grid is the 100 AIRS layers.

Src: FORTRAN source code for the KLAYERS program.

Test: contains the "run_test" script and example RTP files used to
   test out proper working of the package.


-------------------------------------------------------------------------------
It is recommended you follow the procedure detailed below to test out
the code before doing any other work with this package.
-------------------------------------------------------------------------------
We have made the following assumptions:
1) The directory structure and all files mentioned above are as supplied.
2) Your computer uses some flavor of the UNIX operating system, and has
   a FORTRAN 77 (or FORTRAN 90) compiler that allows the non-standard
   "STRUCTURE" variable type.  (This is a common extension to F77.)
3) You need the HDF4 library.  You will need to link to it.
4) You need our RTP package.  You need the library and include files.


How to test KLAYERS:
-------------------

step 1) Verify you have the HDF and RTP packages installed.  Note the
   location of each installion, you'll need this info in step 2.


step 2) Compile the KLAYERS program in "Src".  Edit "make_klayers" as
   needed for your compiler and the location of the RTP and HDF
   libraries and include files.  Edit the include file "incLAY.f"
   for the location of the AFGL models (variable DFAFGL).  When ready
   to compile, type "make".


step 3) Run the test script in the "Test" directory: "sh run_test".
   This will run KLAYERS using the "test_in.rtp" example RTP file
   as input.  If successful, KLAYERS will create a new output file
   called "test_out.rtp".


step 4) Compare the output RTP file you just created in step 3 to
   the "comparison.rtp" file; they should be virtually identical.

   One way to check the output is to use the "rtpdump" utility
   program to create an ASCII dump of some portion of the file
   contents, and then "diff" the ASCII dumps from the two files.
   For example, using the "-p" flag (profile data)
       "rtpdump -p test_out.rtp   > test.p"
       "rtpdump -p comparison.rtp > comp.p"
       "diff test.p comp.p"
   There are likely to be minor differences in the last digits of some
   of the floating point numbers due to machine/compiler differences.
   If "diff" returns too much for you to look, do some spot checks by
   hand to confirm the two output files are essentially the same.


step 5) Read the "description.txt" file in the "Doc" directory.
   This will give you an overview of how the program works and
   how to run it.


---end of file---
