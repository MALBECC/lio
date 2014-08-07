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

ifeq ($(dbg),1)
  DEBUG  = -g
  DEBUG2 = -g
  DEBUG3 = -g
else
  DEBUG  = -O3
  DEBUG2 = -O1
  DEBUG3 = -O3
endif

ifeq ($(gcc),1)
  CXX = gfortran
  PACK_LIBS = -L$(PACK_LIB_DIR) -llapack -lblas
else
  CXX = ifort
  ifneq ($(dbg),1)
    DEBUG  += -mp1 -ip
    DEBUG3 += -parallel
  endif
  PACK_LIBS = -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include -lmkl_lapack95_lp64 -Wl,--start-group -lmkl_intel_lp64\
     -lmkl_intel_thread -lmkl_core -Wl,--end-group -fopenmp
endif

DEFINE   = -Dpack -fPIC   
CFLAGS   = $(DEBUG)  -c -cpp $(DEFINE) $(PROFILE)
CFLAGS2  = $(DEBUG2) -c -cpp $(DEFINE) $(PROFILE)
CFLAGS3  = $(DEBUG3) -c -cpp $(DEFINE) $(PROFILE)
CPPFLAGS = -P -traditional-cpp 
LIBS=-L/usr/lib -L/usr/lib64 -L../g2g -lg2g -lstdc++ $(PACK_LIBS) -lm -Wl,-rpath='$$ORIGIN/' -Wl,-rpath='$$ORIGIN/../g2g' 

ifeq ($(magma),1)
	DEFINE   += -Dmagma
	LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif

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
	$(CXX) -shared -fPIC -o $@ $(GARCHA_OBJ:%.o=obj/garcha_g2g_obj/%.o) $(LIBS)


## Targets: clean ##
clean:
	rm -rf obj
	rm -rf libgarcha-g2g.so liblio-g2g.so
	rm -f garcha_mod.mod


# Different flags for some .o
obj/garcha_g2g_obj/dip2.o: dip2.f
	$(CXX) -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/dipG.o: dipG.f
	$(CXX) -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<
obj/garcha_g2g_obj/dip.o: dip.f
	$(CXX) -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/int1G.o: int1G.f
	$(CXX) -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/int3G.o: int3G.f
	$(CXX) -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/intsolG.o: intsolG.f
	$(CXX) -o $@ $(CFLAGS2) $(G2G_CFLAGS)  $<

obj/garcha_g2g_obj/init.o: init.f
	$(CXX) -o $@ $(CFLAGS) $(G2G_CFLAGS)  $<

obj/garcha_g2g_obj/SCFop.o: SCFop.f
	$(CXX) -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/intsolGs.o: intsolGs.f
	$(CXX) -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/matmulnano.o: matmulnano.f
	$(CXX) -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/conmut.o: conmut.f
	$(CXX) -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<
obj/garcha_g2g_obj/matmuldiag.o: matmuldiag.f
	$(CXX) -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<
obj/garcha_g2g_obj/int3lu.o: int3lu.f
	$(CXX) -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<

obj/garcha_g2g_obj/intsol.o: intsol.f
	$(CXX) -o $@ $(CFLAGS2) $(G2G_CFLAGS)  $<


### Objects


obj/garcha_g2g_obj/%.o: %.f
	$(CXX) -o $@ $(CFLAGS) $(G2G_CFLAGS) $<

