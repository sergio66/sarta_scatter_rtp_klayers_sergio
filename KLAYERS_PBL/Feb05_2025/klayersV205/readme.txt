
   ===================================
   README file for the KLAYERS package
   ===================================
   KLAYERS programmer: Scott Hannon (email: hannon@umbc.edu)
   UMBC Atmospheric Spectroscopy Laboratory (ASL)
   Last updated: 15 March 2002

-------------------------------------------------------------------------------
Introduction:
------------

This package contains the KLAYERS program; you should also have
the RTP package.  You will not be able to compile KLAYERS until
you have built the RTP library, so do that first.

KLAYERS is a program to convert "level" point profiles into
"layer" integrated slab profiles on a particular layering grid.
The makefile is set up to compile using the 100 AIRS layers grid
unless another grid is specified.  See "Doc/description.txt" for
details on using the program.


Overview of the package Contents:
--------------------------------
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
