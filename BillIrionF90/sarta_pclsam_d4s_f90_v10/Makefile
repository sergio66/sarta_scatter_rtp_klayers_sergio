F90OBJS = \
	airs_sarta_pclsam_d4s.o \
	airs_sarta_variables.o \
	incftc.o \
	rdprof.o \
	fnmie_iceaggr_waterdrop_desertdust.o \
	fnmie_iceGHMbaum_waterdrop_desertdust.o \
	rdcldt.o \
	rdcoef.o \
	tunmlt.o \
	rdsun.o \
	getbot.o \
	mean_t.o \
	vaconv.o \
	calpar.o \
	calowp.o \
	calokw.o \
	calt1.o \
	calt2.o \
	calt3.o \
	calt4.o \
	calt4_solar.o \
	calt5.o \
	calt5_solar.o \
	calt6.o \
	calt6_solar.o \
	calt7.o \
	calt7_solar.o \
	faketz.o \
	bkprep.o \
	getmie.o \
	qikexp.o \
	calrad0.o \
	calrad1.o \
	calrad2.o \
	calrad_D4S.o \
	calnte.o \
	hg3.o \
	sunpar.o \
	ccprep_slab.o \
	lennb.o \
	saconv.o \
	ref_profiles.o \
	select_species.o \
	secang_module.o
#	setems.f90 \
#	util.f90 \
#	getcld_slab.o \


MODS = incftc.mod airs_sarta_variables.mod ref_profiles.mod calpar.mod sunpar.mod \
	getbot.mod select_species.mod secang_module.mod
CALMODS = calt1.mod calt2.mod calt3.mod calt4.mod calt4_solar.mod calt5.mod calt5_solar.mod \
	calt6.mod calt6_solar.mod calt7.mod calt7_solar.mod

F90 = gfortran
LD = gfortran

FFLAGS = -fPIC -c -O3 

LD_POST=  -fPIC -lm -lc -shared

INC = 

all : lib

lib : libsarta_forward.so

libsarta_forward.so : ${MODS} ${F90OBJS}
	${LD} -o libsarta_forward.so ${F90OBJS} ${LD_POST}

clean:
	@rm *.mod *.o *.so

incftc.mod : incftc.f90
	${F90} ${INC} ${FFLAGS} incftc.f90 

ref_profiles.mod : ref_profiles.f90 incftc.mod
	${F90} ${INC} ${FFLAGS} ref_profiles.f90 

airs_sarta_variables.mod : airs_sarta_variables.f90 incftc.mod ref_profiles.mod 
	${F90} ${INC} ${FFLAGS} airs_sarta_variables.f90 

calpar.mod : calpar.f90 incftc.mod ref_profiles.mod airs_sarta_variables.mod select_species.mod
	${F90} ${INC} ${FFLAGS} calpar.f90 

sunpar.mod : sunpar.f90 incftc.mod ref_profiles.mod airs_sarta_variables.mod select_species.mod
	${F90} ${INC} ${FFLAGS} sunpar.f90 

airs_sarta_pclsam_d4s.o : airs_sarta_pclsam_d4s.f90 airs_sarta_variables.mod calpar.mod sunpar.mod ${CALMODS} getbot.mod
	${F90} ${INC} ${FFLAGS} airs_sarta_pclsam_d4s.f90 

rdcoef.o: rdcoef.f90 incftc.mod
	${F90} ${INC} ${FFLAGS} rdcoef.f90 

rdprof.o: rdprof.f90 incftc.mod
	${F90} ${INC} ${FFLAGS} rdprof.f90 

secang_module.mod : incftc.mod saconv.o vaconv.o
	${F90} ${INC} ${FFLAGS} secang_module.f90 

fnmie_iceaggr_waterdrop_desertdust.o : fnmie_iceaggr_waterdrop_desertdust.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} fnmie_iceaggr_waterdrop_desertdust.f90 

fnmie_iceGHMbaum_waterdrop_desertdust.o : fnmie_iceGHMbaum_waterdrop_desertdust.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} fnmie_iceGHMbaum_waterdrop_desertdust.f90 
	
rdcldt.o : rdcldt.f90 incftc.mod
	${F90} ${INC} ${FFLAGS} rdcldt.f90 
	
tunmlt.o : tunmlt.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} tunmlt.f90 

rdsun.o : rdsun.f90 incftc.mod
	${F90} ${INC} ${FFLAGS} rdsun.f90 

getbot.mod : getbot.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} getbot.f90 

mean_t.o : mean_t.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} mean_t.f90 

vaconv.o : vaconv.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} vaconv.f90 

calowp.o : calowp.f90 incftc.mod ref_profiles.mod
	${F90} ${INC} ${FFLAGS} calowp.f90 

calokw.o : calokw.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calokw.f90 

calt1.mod : calt1.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt1.f90 

calt2.mod : calt2.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt2.f90 

calt3.mod : calt3.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt3.f90 

calt4.mod : calt4.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt4.f90 

calt4_solar.mod : calt4_solar.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt4_solar.f90 

calt5.mod : calt5.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt5.f90 

calt5_solar.mod : calt5_solar.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt5_solar.f90 

calt6.mod : calt6.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt6.f90 

calt6_solar.mod : calt6_solar.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt6_solar.f90 


calt7.mod : calt7.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt7.f90 

calt7_solar.mod : calt7_solar.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calt7_solar.f90 

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

calrad0.o : calrad0.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calrad0.f90 

calrad1.o : calrad1.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calrad1.f90 

calrad2.o : calrad2.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calrad2.f90 

calrad_D4S.o : calrad_D4S.f90 incftc.mod 
	${F90} ${INC} ${FFLAGS} calrad_D4S.f90 

calnte.o : calnte.f90 
	${F90} ${INC} ${FFLAGS} calnte.f90 

hg3.o : hg3.f90 
	${F90} ${INC} ${FFLAGS} hg3.f90 

ccprep_slab.o : ccprep_slab.f90 
	${F90} ${INC} ${FFLAGS} ccprep_slab.f90 
	
lennb.o : lennb.f90 
	${F90} ${INC} ${FFLAGS} lennb.f90 
	
saconv.o : saconv.f90 
	${F90} ${INC} ${FFLAGS} saconv.f90 

select_species.mod : select_species.f90
	${F90} ${INC} ${FFLAGS} select_species.f90 
	



