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
LIBS=-L/usr/lib $(CUDA_LIBS) -L../g2g -lg2g -lstdc++ -L$(MKL_LIBRARY_PATH) $(PROFILE) -lmkl_lapack64 -lmkl_em64t -lguide -lpthread -lm

# directories
CUDA_COMMONDIR := $(CUDA_INSTALL_PATH)/../sdk/common
CUDA_LIBS=-L$(CUDA_INSTALL_PATH)/lib64 -L$(CUDA_COMMONDIR)/lib64 -lcudart
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
	echo $(GARCHA_OBJ)
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
%/dipG.o: dipG.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
%/dip.o: garcha_src/dip.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
%/int1G.o: garcha_src/int1G.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<
%/intsolG.o: garcha_src/intsolG.f
	ifort -o $@ $(CFLAGS3) $(G2G_CFLAGS) $<



### Objects
obj/garcha_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $<

obj/hibrido_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $<
	
obj/garcha_g2g_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $(G2G_CFLAGS) $<

obj/hibrido_g2g_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $(G2G_CFLAGS) $<

#obj/init.o    : init.f
#	ifort -o $@ -cpp $(COMPFLG) init.f
#obj/int.o    : int.f param
#	ifort -o $@ $(COMPFLG) int.f 
#obj/int1G.o    : int1G.f param
#	ifort -o $@ $(COMPFLG3) int1G.f 
#obj/drive.o : drive.f param 
#	ifort -o $@ $(COMPFLG2) drive.f 
#obj/func.o : func.f
#	ifort -o $@ $(COMPFLG) func.f 
#obj/int2.o    : int2.f param
#	ifort -o $@ $(COMPFLG) int2.f 
#obj/int2G.o    : int2G.f param
#	ifort -o $@ $(COMPFLG) int2G.f 
#obj/int3G.o    : int3G.f param
#	ifort -o $@ $(COMPFLG2) int3G.f
#obj/int3.o    : int3.f param
#	ifort -o $@ $(COMPFLG2) int3.f
#obj/int3N.o    : int3N.f param
#	ifort -o $@ $(COMPFLG2) int3N.f
#obj/int3lu.o    : int3lu.f param
#	ifort -o $@ $(COMPFLG2) int3lu.f
#obj/int3mem.o    : int3mem.f param
#	ifort -o $@ $(COMPFLG2) int3mem.f
#obj/exch.o    : exch.f param
#	ifort -o $@ $(COMPFLG) exch.f
#obj/exch2.o    : exch2.f param
#	ifort -o $@ $(COMPFLG) exch2.f
#obj/pot.o    : pot.f 
#	ifort -o $@ $(COMPFLG) pot.f 
#obj/SCF.o    : SCF.f param
#	ifort -o $@ $(COMPFLG) SCF.f
#obj/MD2.o     : MD2.f param
#	ifort -o $@ $(COMPFLG) MD2.f 
#obj/geom.o     : geom.f param
#	ifort -o $@ $(COMPFLG) geom.f
#obj/intSG.o     : intSG.f param
#	ifort -o $@ $(COMPFLG) intSG.f 
#obj/dfp2.o     : dfp2.f param
#	ifort -o $@ $(COMPFLG) dfp2.f 
#obj/lsearch.o     : lsearch.f param
#	ifort -o $@ $(COMPFLG) lsearch.f 
#obj/write.o     : write.f param
#	ifort -o $@ $(COMPFLG) write.f
#obj/nwrite.o     : nwrite.f param
#	ifort -o $@ $(COMPFLG) nwrite.f
#obj/dns.o     : dns.f param
#	ifort -o $@ -cpp $(COMPFLG) dns.f 
#obj/dns2.o     : dns2.f param
#	ifort -o $@ $(COMPFLG) dns2.f 
#obj/dnsg.o     : dnsg.f param
#	ifort -o $@ $(COMPFLG) dnsg.f 
#obj/densg.o     : densg.f param
#	ifort -o $@ $(COMPFLG) densg.f 
#obj/potg.o     : potg.f 
#	ifort -o $@ $(COMPFLG) potg.f 
#obj/dip.o     : dip.f param
#	ifort -o $@ $(COMPFLG) dip.f 
#obj/dipG.o     : dipG.f param
#	ifort -o $@ $(COMPFLG3) dipG.f 
#obj/dip2.o     : dip2.f param
#	ifort -o $@ $(COMPFLG) dip2.f 
#obj/vol.o     : vol.f param
#	ifort -o $@ $(COMPFLG) vol.f 
#obj/exchnum.o     : exchnum.f param
#	ifort -o $@ $(COMPFLG) exchnum.f 
#obj/grid.o     : grid.f
#	ifort -o $@ $(COMPFLG) grid.f 
#obj/exchnum2.o     : exchnum2.f param
#	ifort -o $@ $(COMPFLG) exchnum2.f 
#obj/exchfock.o     : exchfock.f param
#	ifort -o $@ $(COMPFLG) exchfock.f 
#obj/SCFop.o    : SCFop.f param
#	ifort -o $@ $(COMPFLG) SCFop.f
#obj/densgop.o     : densgop.f param
#	ifort -o $@ $(COMPFLG) densgop.f
#obj/dns2op.o     : dns2op.f param
#	ifort -o $@ $(COMPFLG) dns2op.f
#obj/dnsop.o     : dnsop.f param
#	ifort -o $@ $(COMPFLG) dnsop.f
#obj/exchnum2op.o     : exchnum2op.f param
#	ifort -o $@ $(COMPFLG) exchnum2op.f
#obj/exchnumop.o     : exchnumop.f param
#	ifort -o $@ $(COMPFLG) exchnumop.f
#obj/potop.o     : potop.f param
#	ifort -o $@ $(COMPFLG) potop.f
#obj/potgop.o     : potgop.f 
#	ifort -o $@ $(COMPFLG) potgop.f
#obj/dnsgop.o     : dnsgop.f param
#	ifort -o $@ $(COMPFLG) dnsgop.f
#obj/alg.o     : alg.f 
#	ifort -o $@ $(COMPFLG) alg.f
#obj/eig.o     : eig.f
#	ifort -o $@ $(COMPFLG2) eig.f
#obj/svd.o     : svd.f
#	ifort -o $@ $(COMPFLG) svd.f
#obj/elec.o     : elec.f param
#	ifort -o $@ $(COMPFLG) elec.f 
#obj/charge.o   : charge.f param
#	ifort -o $@ $(COMPFLG) charge.f
#obj/intsol.o    : intsol.f param
#	ifort -o $@ $(COMPFLG) intsol.f
#obj/mmsol.o    : mmsol.f param
#	ifort -o $@ $(COMPFLG) mmsol.f 
#obj/intsolG.o    : intsolG.f param
#	ifort -o $@ $(COMPFLG2) intsolG.f 
#obj/mmsolG.o    : mmsolG.f param
#	ifort -o $@ $(COMPFLG3) mmsolG.f 
#obj/resp.o    : resp.f param
#	ifort -o $@ $(COMPFLG) resp.f 
#obj/efield.o    : efield.f param
#	ifort -o $@ $(COMPFLG) efield.f 
#obj/lalg.o    : lalg.f param
#	ifort -o $@ $(COMPFLG) lalg.f 
#obj/popu.o    : popu.f param
#	ifort -o $@ $(COMPFLG) popu.f 
