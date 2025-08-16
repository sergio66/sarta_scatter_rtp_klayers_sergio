make clean; make -f Makefile airs_wetwater
make clean; make -f Makefile airs_samewater

echo "now going to ../BinV201_chip;  ln -s klayers_airs_wetwater klayers_airs"
cd ../BinV201_chip;  rm klayers_airs; ln -s klayers_airs_wetwater klayers_airs
