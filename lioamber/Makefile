# sources and objects
MOD_SRCS=$(wildcard mods/*.f)
# COMMON_SRCS=$(wildcard *.f)
GARCHA_SRCS=$(wildcard *.f)

MOD_OBJ=$(MOD_SRCS:%.f=%.o) 
GARCHA_OBJ=$(MOD_OBJ) $(COMMON_SRCS:%.f=%.o) $(GARCHA_SRCS:%.f=%.o)

# flags
ifeq ($(profile),1)
  PROFILE  = -pg
else
  PROFILE  =  
endif

ifeq ($(cpu_recompute),0)
  CXXFLAGS += -DCPU_RECOMPUTE=0
else
  CXXFLAGS += -DCPU_RECOMPUTE=1
endif

ifeq ($(full_double),1)
  CXXFLAGS += -DFULL_DOUBLE=1
  NVCCFLAGS += -DFULL_DOUBLE=1
else
  CXXFLAGS += -DFULL_DOUBLE=0
  NVCCFLAGS += -DFULL_DOUBLE=0
endif

DEFINE   = -Dpack -fPIC   
CFLAGS   = -g -mp1 -ip -O3 -c -cpp  $(DEFINE) $(PROFILE)
CFLAGS2  = -g -O1 -c -cpp  $(DEFINE) $(PROFILE)
CFLAGS3  = -g -O3 -c -cpp -parallel   $(DEFINE) $(PROFILE)
CPPFLAGS = -P -traditional-cpp 
LIBS=-L/usr/lib -L/usr/lib64 -L../g2g -lg2g -lstdc++ -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include $(PROFILE)  -lmkl_lapack95_lp64 -Wl,--start-group -lmkl_intel_lp64\
     -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread -lm -Wl,-rpath='$$ORIGIN/' -Wl,-rpath='$$ORIGIN/../g2g' 
#eliminada   -lmkl_solver_lp64 

# directories
G2G_CFLAGS=-DG2G

ifeq ($(ultima_garcha),1)
	G2G_CFLAGS += -DULTIMA_GARCHA
endif

ifeq ($(convergencia_garcha),1)
	G2G_CFLAGS += -DINT3N_GARCHA
endif

ifeq ($(ultima_g2g),1)
	CFLAGS += -DULTIMA_G2G
	CFLAGS2 += -DULTIMA_G2G
endif

## Targets: make ##

all: dirs liblio-g2g.so 

dirs:
	mkdir -p obj/garcha_g2g_obj/mods obj/garcha_g2g_obj/

liblio-g2g.so: dirs $(GARCHA_OBJ:%.o=obj/garcha_g2g_obj/%.o)
	ifort -shared -fPIC -o $@ $(GARCHA_OBJ:%.o=obj/garcha_g2g_obj/%.o) $(LIBS)


## Targets: clean ##
clean:
	rm -rf obj
	rm -rf libgarcha-g2g.so liblio-g2g.so


# Different flags for some .o
obj/garcha_g2g_obj/dip2.o: dip2.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/dipG.o: dipG.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<
obj/garcha_g2g_obj/dip.o: dip.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/int1G.o: int1G.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/int3G.o: int3G.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/intsolG.o: intsolG.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS)  $<

obj/garcha_g2g_obj/init.o: init.f
	ifort -o $@ $(CFLAGS) $(G2G_CFLAGS)  $<

obj/garcha_g2g_obj/SCFop.o: SCFop.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/intsolGs.o: intsolGs.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/matmulnano.o: matmulnano.f
	ifort -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/conmut.o: conmut.f
	ifort -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<
obj/garcha_g2g_obj/matmuldiag.o: matmuldiag.f
	ifort -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<
obj/garcha_g2g_obj/int3lu.o: int3lu.f
	ifort -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/intsol.o: intsol.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS)  $<


### Objects


obj/garcha_g2g_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $(G2G_CFLAGS) $<

