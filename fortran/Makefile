#-*- mode: make -*-

# sources and objects
COMMON_SRCS=$(wildcard *.f)
GARCHA_SRCS=$(wildcard garcha_src/*.f)
HIBRIDO_SRCS=$(wildcard hibrido_src/*.f)

GARCHA_OBJ=$(COMMON_SRCS:%.f=%.o) $(GARCHA_SRCS:%.f=%.o)
HIBRIDO_OBJ=$(COMMON_SRCS:%.f=%.o) $(HIBRIDO_SRCS:%.f=%.o)

# flags
ifeq ($(profile),1)
  PROFILE  = -pg
else
  PROFILE  =
endif

DEFINE   = -Dpack
CFLAGS   = -mp1 -ip -O3 -c -cpp $(DEFINE) $(PROFILE)
CFLAGS2  = -O3 -c -cpp $(DEFINE) $(PROFILE)
CFLAGS3  = -O1 -c -cpp $(DEFINE) $(PROFILE)
CPPFLAGS = -P -traditional-cpp
LIBS=-L/usr/lib $(CUDA_LIBS) -L../g2g -lg2g -lstdc++ -L$(MKL_LIBRARY_PATH) $(PROFILE) -lmkl_solver_ilp64_sequential -lblas -lmkl_core -lmkl_sequential -lguide -lpthread -lm -llapack

# directories
CUDA_COMMONDIR := $(CUDA_INSTALL_PATH)/../sdk/common
CUDA_LIBS=-L$(CUDA_INSTALL_PATH)/lib64 -L$(CUDA_COMMONDIR)/lib64 -lcuda -lcudart
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

all: dirs garcha-g2g hibrido-g2g garcha hibrido

dirs:
	mkdir -p obj/garcha_g2g_obj/garcha_src
	mkdir -p obj/hibrido_g2g_obj/hibrido_src
	mkdir -p obj/garcha_obj/garcha_src
	mkdir -p obj/hibrido_obj/hibrido_src
	
garcha-g2g: $(GARCHA_OBJ:%.o=obj/garcha_g2g_obj/%.o)
	ifort -o $@ $(GARCHA_OBJ:%.o=obj/garcha_g2g_obj/%.o) $(LIBS) $(CUDA_LIBS)

hibrido-g2g: $(HIBRIDO_OBJ:%.o=obj/hibrido_g2g_obj/%.o)	
	ifort -o $@ $(HIBRIDO_OBJ:%.o=obj/hibrido_g2g_obj/%.o) $(LIBS) $(CUDA_LIBS)
	
garcha: $(GARCHA_OBJ:%.o=obj/garcha_obj/%.o)
	ifort -o $@ $(GARCHA_OBJ:%.o=obj/garcha_obj/%.o) $(LIBS)
	
hibrido: $(HIBRIDO_OBJ:%.o=obj/hibrido_obj/%.o)
	ifort -o $@ $(HIBRIDO_OBJ:%.o=obj/hibrido_obj/%.o) $(LIBS)
	
## Targets: clean ##
clean:
	rm -rf obj
	rm -rf garcha garcha-g2g hibrido hibrido-g2g


# Different flags for some .o
obj/garcha_g2g_obj/dip2.o: dip2.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/garcha_obj/dip2.o: dip2.f
	ifort -o $@ $(CFLAGS3) $<
obj/hibrido_g2g_obj/dip2.o: dip2.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/hibrido_obj/dip2.o: dip2.f
	ifort -o $@ $(CFLAGS3) $<

obj/garcha_g2g_obj/dipG.o: dipG.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/garcha_obj/dipG.o: dipG.f
	ifort -o $@ $(CFLAGS3) $<
obj/hibrido_g2g_obj/dipG.o: dipG.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/hibrido_obj/dipG.o: dipG.f
	ifort -o $@ $(CFLAGS3) $<

obj/garcha_g2g_obj/garcha_src/dip.o: garcha_src/dip.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/garcha_obj/garcha_src/dip.o: garcha_src/dip.f
	ifort -o $@ $(CFLAGS3) $<
obj/hibrido_g2g_obj/hibrido_src/dip.o: hibrido_src/dip.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/hibrido_obj/hibrido_src/dip.o: hibrido_src/dip.f
	ifort -o $@ $(CFLAGS3) $<

obj/garcha_g2g_obj/garcha_src/int1G.o: garcha_src/int1G.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/garcha_obj/garcha_src/int1G.o: garcha_src/int1G.f
	ifort -o $@ $(CFLAGS3) $<
obj/hibrido_g2g_obj/hibrido_src/int1G.o: hibrido_src/int1G.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/hibrido_obj/hibrido_src/int1G.o: hibrido_src/int1G.f
	ifort -o $@ $(CFLAGS3) $<

obj/garcha_g2g_obj/garcha_src/int3G.o: garcha_src/int3G.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/garcha_obj/garcha_src/int3G.o: garcha_src/int3G.f
	ifort -o $@ $(CFLAGS3) $<
obj/hibrido_g2g_obj/hibrido_src/int3G.o: hibrido_src/int3G.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/hibrido_obj/hibrido_src/int3G.o: hibrido_src/int3G.f
	ifort -o $@ $(CFLAGS3) $<

obj/garcha_g2g_obj/garcha_src/intsolG.o: garcha_src/intsolG.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/garcha_obj/garcha_src/intsolG.o: garcha_src/intsolG.f
	ifort -o $@ $(CFLAGS3) $<
obj/hibrido_g2g_obj/hibrido_src/intsolG.o: hibrido_src/intsolG.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
obj/hibrido_obj/hibrido_src/intsolG.o: hibrido_src/intsolG.f
	ifort -o $@ $(CFLAGS3) $<

### Objects
obj/garcha_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $<

obj/hibrido_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $<
	
obj/garcha_g2g_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $(G2G_CFLAGS) $<

obj/hibrido_g2g_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $(G2G_CFLAGS) $<
