#-*- mode: make -*-
SRCS = $(wildcard *.f)
OBJ=$(SRCS:%.f=%.o)
DEFINE   = -Dpack
COMPFLG  = -mp1 -ip -O3 -c $(CFLAGS) -cpp $(DEFINE)
COMPFLG2 = -O3 -c $(CFLAGS) -cpp $(DEFINE)
CPPFLAGS = -P -traditional-cpp
LIBS= -L/usr/lib -lblas -L../gpu -lg2g -lstdc++ $(CUDA_LIBS) -L$(MKL_LIBRARY_PATH) -lmkl_lapack64 -lmkl_em64t -lguide \
 	-lpthread -lm

CUDA_COMMONDIR := $(CUDA_INSTALL_PATH)/../sdk/common
CUDA_LIBS=-L$(CUDA_INSTALL_PATH)/lib64 -L$(CUDA_COMMONDIR)/lib64 -lcuda -lcudart
G2G_CFLAGS=-DG2G

ifeq ($(profile),1)
	COMPFLG += -pg
	COMPFLG2 += -pg
	LIBS += -pg
endif

ifeq ($(ultima_garcha),1)
	G2G_CFLAGS += -DULTIMA_CPU
endif

ifeq ($(convergencia_garcha),1)
	G2G_CFLAGS += -DINT3N_CPU
endif

ifeq ($(ultima_g2g),1)
	COMPFLG += -DULTIMA_G2G
	COMPFLG2 += -DULTIMA_G2G
endif

all: garcha-g2g garcha

clean-original:
	rm -f obj/* garcha
	
clean-g2g:
	rm -f g2g_obj/* garcha-g2g	
	
clean: clean-original clean-g2g

garcha-g2g: $(OBJ:%.o=g2g_obj/%.o)
	ifort -o $@ $(OBJ:%.o=g2g_obj/%.o) $(LIBS) $(EXTRA_FLAGS)
	
garcha: $(OBJ:%.o=obj/%.o)
	ifort -o $@ $(OBJ:%.o=obj/%.o) $(LIBS) $(EXTRA_FLAGS)

# Original code object files
obj/init.o    : init.f
	ifort -o $@ -cpp $(COMPFLG) init.f
obj/int.o    : int.f param
	ifort -o $@ $(COMPFLG) int.f 
obj/int1G.o    : int1G.f param
	ifort -o $@ $(COMPFLG2) int1G.f 
obj/drive.o : drive.f param 
	ifort -o $@ $(COMPFLG2) drive.f 
obj/func.o : func.f
	ifort -o $@ $(COMPFLG) func.f 
obj/int2.o    : int2.f param
	ifort -o $@ $(COMPFLG) int2.f 
obj/int2G.o    : int2G.f param
	ifort -o $@ $(COMPFLG) int2G.f 
obj/int3G.o    : int3G.f param
	ifort -o $@ $(COMPFLG2) int3G.f
obj/int3.o    : int3.f param
	ifort -o $@ $(COMPFLG2) int3.f
obj/int3N.o    : int3N.f param
	ifort -o $@ $(COMPFLG2) int3N.f
obj/int3lu.o    : int3lu.f param
	ifort -o $@ $(COMPFLG2) int3lu.f
obj/int3mem.o    : int3mem.f param
	ifort -o $@ $(COMPFLG2) int3mem.f
obj/exch.o    : exch.f param
	ifort -o $@ $(COMPFLG) exch.f
obj/exch2.o    : exch2.f param
	ifort -o $@ $(COMPFLG) exch2.f
obj/pot.o    : pot.f 
	ifort -o $@ $(COMPFLG) pot.f 
obj/SCF.o    : SCF.f param
	ifort -o $@ $(COMPFLG) SCF.f
obj/MD2.o     : MD2.f param
	ifort -o $@ $(COMPFLG) MD2.f 
obj/geom.o     : geom.f param
	ifort -o $@ $(COMPFLG) geom.f
obj/intSG.o     : intSG.f param
	ifort -o $@ $(COMPFLG) intSG.f 
obj/dfp2.o     : dfp2.f param
	ifort -o $@ $(COMPFLG) dfp2.f 
obj/lsearch.o     : lsearch.f param
	ifort -o $@ $(COMPFLG) lsearch.f 
obj/write.o     : write.f param
	ifort -o $@ $(COMPFLG) write.f
obj/nwrite.o     : nwrite.f param
	ifort -o $@ $(COMPFLG) nwrite.f
obj/dns.o     : dns.f param
	ifort -o $@ -cpp $(COMPFLG) dns.f 
obj/dns2.o     : dns2.f param
	ifort -o $@ $(COMPFLG) dns2.f 
obj/dnsg.o     : dnsg.f param
	ifort -o $@ $(COMPFLG) dnsg.f 
obj/densg.o     : densg.f param
	ifort -o $@ $(COMPFLG) densg.f 
obj/potg.o     : potg.f 
	ifort -o $@ $(COMPFLG) potg.f 
obj/dip.o     : dip.f param
	ifort -o $@ $(COMPFLG) dip.f 
obj/dipG.o     : dipG.f param
	ifort -o $@ $(COMPFLG) dipG.f 
obj/dip2.o     : dip2.f param
	ifort -o $@ $(COMPFLG) dip2.f 
obj/vol.o     : vol.f param
	ifort -o $@ $(COMPFLG) vol.f 
obj/exchnum.o     : exchnum.f param
	ifort -o $@ $(COMPFLG) exchnum.f 
obj/grid.o     : grid.f
	ifort -o $@ $(COMPFLG) grid.f 
obj/exchnum2.o     : exchnum2.f param
	ifort -o $@ $(COMPFLG) exchnum2.f 
obj/exchfock.o     : exchfock.f param
	ifort -o $@ $(COMPFLG) exchfock.f 
obj/SCFop.o    : SCFop.f param
	ifort -o $@ $(COMPFLG) SCFop.f
obj/densgop.o     : densgop.f param
	ifort -o $@ $(COMPFLG) densgop.f
obj/dns2op.o     : dns2op.f param
	ifort -o $@ $(COMPFLG) dns2op.f
obj/dnsop.o     : dnsop.f param
	ifort -o $@ $(COMPFLG) dnsop.f
obj/exchnum2op.o     : exchnum2op.f param
	ifort -o $@ $(COMPFLG) exchnum2op.f
obj/exchnumop.o     : exchnumop.f param
	ifort -o $@ $(COMPFLG) exchnumop.f
obj/potop.o     : potop.f param
	ifort -o $@ $(COMPFLG) potop.f
obj/potgop.o     : potgop.f 
	ifort -o $@ $(COMPFLG) potgop.f
obj/dnsgop.o     : dnsgop.f param
	ifort -o $@ $(COMPFLG) dnsgop.f
obj/alg.o     : alg.f 
	ifort -o $@ $(COMPFLG) alg.f
obj/eig.o     : eig.f
	ifort -o $@ $(COMPFLG2) eig.f
obj/svd.o     : svd.f
	ifort -o $@ $(COMPFLG) svd.f
obj/elec.o     : elec.f param
	ifort -o $@ $(COMPFLG) elec.f 
obj/charge.o   : charge.f param
	ifort -o $@ $(COMPFLG) charge.f
obj/intsol.o    : intsol.f param
	ifort -o $@ $(COMPFLG) intsol.f 
obj/mmsol.o    : mmsol.f param
	ifort -o $@ $(COMPFLG) mmsol.f 
obj/intsolG.o    : intsolG.f param
	ifort -o $@ $(COMPFLG2) intsolG.f 
obj/mmsolG.o    : mmsolG.f param
	ifort -o $@ $(COMPFLG) mmsolG.f 
obj/resp.o    : resp.f param
	ifort -o $@ $(COMPFLG) resp.f 
obj/efield.o    : efield.f param
	ifort -o $@ $(COMPFLG) efield.f 
obj/lalg.o    : lalg.f param
	ifort -o $@ $(COMPFLG) lalg.f 
obj/popu.o    : popu.f param
	ifort -o $@ $(COMPFLG) popu.f 

# G2G Object files	
g2g_obj/init.o    : init.f
	ifort -o $@ -cpp $(COMPFLG) $(G2G_CFLAGS) init.f 
g2g_obj/int.o    : int.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) int.f 
g2g_obj/int1G.o    : int1G.f param
	ifort -o $@ $(COMPFLG2) $(G2G_CFLAGS) int1G.f 
g2g_obj/drive.o : drive.f param 
	ifort -o $@ $(COMPFLG2) $(G2G_CFLAGS) drive.f 
g2g_obj/func.o : func.f
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) func.f 
g2g_obj/int2.o    : int2.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) int2.f 
g2g_obj/int2G.o    : int2G.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) int2G.f 
g2g_obj/int3G.o    : int3G.f param
	ifort -o $@ $(COMPFLG2) $(G2G_CFLAGS) int3G.f
g2g_obj/int3.o    : int3.f param
	ifort -o $@ $(COMPFLG2) $(G2G_CFLAGS) int3.f
g2g_obj/int3N.o    : int3N.f param
	ifort -o $@ $(COMPFLG2) $(G2G_CFLAGS) int3N.f
g2g_obj/int3lu.o    : int3lu.f param
	ifort -o $@ $(COMPFLG2) $(G2G_CFLAGS) int3lu.f
g2g_obj/int3mem.o    : int3mem.f param
	ifort -o $@ $(COMPFLG2) $(G2G_CFLAGS) int3mem.f
g2g_obj/exch.o    : exch.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) exch.f
g2g_obj/exch2.o    : exch2.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) exch2.f
g2g_obj/pot.o    : pot.f 
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) pot.f 
g2g_obj/SCF.o    : SCF.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) SCF.f
g2g_obj/MD2.o     : MD2.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) MD2.f 
g2g_obj/geom.o     : geom.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) geom.f
g2g_obj/intSG.o     : intSG.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) intSG.f 
g2g_obj/dfp2.o     : dfp2.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) dfp2.f 
g2g_obj/lsearch.o     : lsearch.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) lsearch.f 
g2g_obj/write.o     : write.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) write.f
g2g_obj/nwrite.o     : nwrite.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) nwrite.f
g2g_obj/dns.o     : dns.f param
	ifort -o $@ -cpp $(COMPFLG) $(G2G_CFLAGS) dns.f 
g2g_obj/dns2.o     : dns2.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) dns2.f 
g2g_obj/dnsg.o     : dnsg.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) dnsg.f 
g2g_obj/densg.o     : densg.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) densg.f 
g2g_obj/potg.o     : potg.f 
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) potg.f 
g2g_obj/dip.o     : dip.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) dip.f 
g2g_obj/dipG.o     : dipG.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) dipG.f 
g2g_obj/dip2.o     : dip2.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) dip2.f 
g2g_obj/vol.o     : vol.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) vol.f 
g2g_obj/exchnum.o     : exchnum.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) exchnum.f 
g2g_obj/grid.o     : grid.f
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) grid.f 
g2g_obj/exchnum2.o     : exchnum2.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) exchnum2.f 
g2g_obj/exchfock.o     : exchfock.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) exchfock.f 
g2g_obj/SCFop.o    : SCFop.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) SCFop.f
g2g_obj/densgop.o     : densgop.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) densgop.f
g2g_obj/dns2op.o     : dns2op.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) dns2op.f
g2g_obj/dnsop.o     : dnsop.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) dnsop.f
g2g_obj/exchnum2op.o     : exchnum2op.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) exchnum2op.f
g2g_obj/exchnumop.o     : exchnumop.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) exchnumop.f
g2g_obj/potop.o     : potop.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) potop.f
g2g_obj/potgop.o     : potgop.f 
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) potgop.f
g2g_obj/dnsgop.o     : dnsgop.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) dnsgop.f
g2g_obj/alg.o     : alg.f 
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) alg.f
g2g_obj/eig.o     : eig.f
	ifort -o $@ $(COMPFLG2) $(G2G_CFLAGS) eig.f
g2g_obj/svd.o     : svd.f
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) svd.f
g2g_obj/elec.o     : elec.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) elec.f 
g2g_obj/charge.o   : charge.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) charge.f
g2g_obj/intsol.o    : intsol.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) intsol.f 
g2g_obj/mmsol.o    : mmsol.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) mmsol.f 
g2g_obj/intsolG.o    : intsolG.f param
	ifort -o $@ $(COMPFLG2) $(G2G_CFLAGS) intsolG.f 
g2g_obj/mmsolG.o    : mmsolG.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) mmsolG.f 
g2g_obj/resp.o    : resp.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) resp.f 
g2g_obj/efield.o    : efield.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) efield.f 
g2g_obj/lalg.o    : lalg.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) lalg.f 
g2g_obj/popu.o    : popu.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) popu.f 
