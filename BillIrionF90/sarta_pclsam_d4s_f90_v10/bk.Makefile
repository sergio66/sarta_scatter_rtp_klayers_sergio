F90OBJS = \
	airs_sarta_pclsam.o \
	airs_sarta_variables.o \
	incftc.o \
	rdprof.o \
	fnmie_iceaggr_waterdrop_desertdust.o \
	fnmie_iceGHMbaum_waterdrop_desertdust.o \
	rdcldt_netcdf.o \
	rdcoef.o \
	tunmlt.o \
	rdsun.o \
	getbot.o \
	mean_t.o \
	vaconv.o \
	calpar.o \
	calowp.o \
	calokw.o \
	calt1_od.o \
	calt2_od.o \
	calt3_od.o \
	calt4_od.o \
	calt5_od.o \
	calt6_od.o \
	calt7_od.o \
	faketz.o \
	bkprep.o \
	getmie.o \
	qikexp.o \
	calrad0.o \
	calrad1.o \
	calrad2.o \
	calnte.o \
	hg3.o \
	sunpar.o \
	ccprep_slab.o \
	lennb.o \
	saconv.o \
	test_mie_netcdf.o 
#	setems.f90 \
#	util.f90 \
#	getcld_slab.o \


MODS = incftc.mod airs_sarta_variables.mod 

F90 = gfortran
LD = gfortran


#    Compiler/loader options.f90or abs.f90t.f9090 compiler
#               -N26.f90or reading Big-Endian binary.f90iles
#               -W 132.f90or 132 wide.f90ixed width .f90.f90iles
#               -YEXT_NAMES=LCS    external names are in lower case
#               -YEXT_SFX=_ external names are appended by "_"
#               .fPIC generate position independent code
#FFLAGS = -O -ffree-form -fPIC

#FFLAGS = -fPIC -c -O3 -fbounds-check 
FFLAGS = -fPIC -c -O 

LIB_FLAGS = -L/opt/packages/netcdf-4.2f_64/lib \
     -L/opt/packages/netcdf-4.2.1_64/lib

WLIB_FLAGS = -Wl,-rpath=/opt/packages/netcdf-4.2f_64/lib \
     -Wl,-rpath=/opt/packages/netcdf-4.2.1_64/lib

LD_POST=  -fPIC -lnetcdf -lnetcdff -lm -lc -shared


INC = -I/opt/packages/netcdf-4.2f_64/include

all : lib

lib : libsarta_forward.so

libsarta_forward.so : ${MODS} ${F90OBJS}
	${LD} $(INC) -o libsarta_forward.so ${F90OBJS}  $(LIB_FLAGS)   ${LD_POST}

clean:
	@rm *.mod *.o *.so

test_mie_netcdf.o : test_mie_netcdf.f90
	${F90} ${INC} ${FFLAGS} test_mie_netcdf.f90 

incftc.mod : incftc.f90
	${F90} ${INC} ${FFLAGS} incftc.f90 

airs_sarta_variables.mod : airs_sarta_variables.f90 incftc.mod
	${F90} ${INC} ${FFLAGS} airs_sarta_variables.f90 

airs_sarta_pclsam.o : airs_sarta_pclsam.f90 airs_sarta_variables.mod
	${F90} ${INC} ${FFLAGS} airs_sarta_pclsam.f90 

rdcoef.o: rdcoef.f90 incftc.mod
	${F90} ${INC} ${FFLAGS} rdcoef.f90 

rdprof.o: rdprof.f90 incftc.mod
	${F90} ${INC} ${FFLAGS} rdprof.f90 

fnmie_iceaggr_waterdrop_desertdust.o : fnmie_iceaggr_waterdrop_desertdust.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} fnmie_iceaggr_waterdrop_desertdust.f90 

fnmie_iceGHMbaum_waterdrop_desertdust.o : fnmie_iceGHMbaum_waterdrop_desertdust.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} fnmie_iceGHMbaum_waterdrop_desertdust.f90 
	
rdcldt_netcdf.o : rdcldt_netcdf.f90 incftc.mod
	${F90} ${LD_FLAGS} ${INC} ${FFLAGS} rdcldt_netcdf.f90 
	
tunmlt.o : tunmlt.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} tunmlt.f90 

rdsun.o : rdsun.f90 incftc.mod
	${F90} ${INC} ${FFLAGS} rdsun.f90 

getbot.o : getbot.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} getbot.f90 

mean_t.o : mean_t.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} mean_t.f90 

vaconv.o : vaconv.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} vaconv.f90 

calpar.o : calpar.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calpar.f90 

calowp.o : calowp.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calowp.f90 

calokw.o : calokw.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calokw.f90 

calt1_od.o : calt1_od.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt1_od.f90 

calt2_od.o : calt2_od.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt2_od.f90 

calt3_od.o : calt3_od.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt3_od.f90 

calt4_od.o : calt4_od.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt4_od.f90 

calt5_od.o : calt5_od.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt5_od.f90 

calt6_od.o : calt6_od.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt6_od.f90 

calt7_od.o : calt7_od.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt7_od.f90 

faketz.o : faketz.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} faketz.f90 

#getcld_slab.o : getcld_slab.f90 incftc.mod 
#	${F90} ${INC} ${FFLAGS} getcld_slab.f90 

bkprep.o : bkprep.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} bkprep.f90 

getmie.o : getmie.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} getmie.f90 

qikexp.o : qikexp.f90 
	${F90} ${INC} ${FFLAGS} qikexp.f90 

calrad0.o : calrad0.f90 
	${F90} ${INC} ${FFLAGS} calrad0.f90 

calrad1.o : calrad1.f90 
	${F90} ${INC} ${FFLAGS} calrad1.f90 

calrad2.o : calrad2.f90 
	${F90} ${INC} ${FFLAGS} calrad2.f90 

calnte.o : calnte.f90 
	${F90} ${INC} ${FFLAGS} calnte.f90 

hg3.o : hg3.f90 
	${F90} ${INC} ${FFLAGS} hg3.f90 

sunpar.o : sunpar.f90 
	${F90} ${INC} ${FFLAGS} sunpar.f90 
	
ccprep_slab.o : ccprep_slab.f90 
	${F90} ${INC} ${FFLAGS} ccprep_slab.f90 
	
lennb.o : lennb.f90 
	${F90} ${INC} ${FFLAGS} lennb.f90 
	
saconv.o : saconv.f90 
	${F90} ${INC} ${FFLAGS} saconv.f90 
	



