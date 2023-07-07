## obviously have to edit the date ....
## /bin/cp -a *.f *.f90 *.m *.nml *.sbatch *.sc *eadme* WORKS_DecXYZ_2018

echo "do you need to copy over any special files eg Makefiles? had to add this to the SARTA F90 code  Makefile  Makefile_pclsam  make_sarta_pclsam"

/bin/cp -a *.f *.f90 *.m *.nml *.sbatch *.sc *eadme* $1 
