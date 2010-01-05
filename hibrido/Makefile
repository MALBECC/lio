#-*- mode: make -*-

SRCS=$(wildcard *.f)
OBJ =$(SRCS:%.f=%.o)
DEFINE   = -Dpack
COMPFLG  = -mp1 -ip -O3 -c $(CFLAGS) -cpp $(DEFINE)
COMPFLG2 = -O3     -c $(CFLAGS) -cpp $(DEFINE)
#more or less general
CPPFLAGS =  -P -traditional-cpp
LIBS=-L/usr/lib -lblas -lstdc++ -L/opt/intel/mkl/9.0/lib/em64t -lmkl_em64t -lguide -lpthread  -L../gpu -lg2g -lm

CUDA_COMMONDIR := $(CUDA_INSTALL_PATH)/../sdk/common
CUDA_LIBS=-L$(CUDA_INSTALL_PATH)/lib64 -L$(CUDA_COMMONDIR)/lib64 -lcuda -lcudart
CPU_LIBS=$(LIBS)
G2G_LIBS=$(LIBS) $(CUDA_LIBS) 
G2G_CFLAGS=-DG2G

ifeq ($(profile),1)
	COMPFLG += -pg
	COMPFLG2 += -pg	
	LIBS += -pg
endif

ifeq ($(ultima_cpu),1)
	G2G_CFLAGS += -DULTIMA_GARCHA
endif

ifeq ($(convergencia_cpu),1)
	G2G_CFLAGS += -DINT3N_GARCHA
endif

ifeq ($(ultima_g2g),1)
	COMPFLG += -DULTIMA_G2G
	COMPFLG2 += -DULTIMA_G2G
endif

all: hibrido-g2g hibrido

clean-original:
	rm -f obj/* hibrido
	
clean-g2g:
	rm -f g2g_obj/* hibrido-g2g
	
clean: clean-original clean-g2g

hibrido-g2g: $(OBJ:%.o=g2g_obj/%.o)
	ifort -o $@ $(OBJ:%.o=g2g_obj/%.o) $(G2G_LIBS) $(EXTRA_FLAGS)
	
hibrido: $(OBJ:%.o=obj/%.o)
	ifort -o $@ $(OBJ:%.o=obj/%.o) $(CPU_LIBS) $(EXTRA_FLAGS)

# CPU OBJs
obj/config.o    : config.f COMM
	ifort -o $@  $(COMPFLG) config.f
obj/corecte.o   : corecte.f COMM
	ifort -o $@  $(COMPFLG) corecte.f
obj/finte.o     : finte.f COMM
	ifort -o $@  $(COMPFLG) finte.f
obj/cvemot.o    : cvemot.f COMM
	ifort -o $@  $(COMPFLG) cvemot.f 
obj/cvemotn.o   : cvemotn.f COMM
	ifort -o $@  $(COMPFLG) cvemotn.f 
obj/cvemotnn.o  : cvemotnn.f COMM
	ifort -o $@  $(COMPFLG) cvemotnn.f
obj/cvemotnnn.o  : cvemotnnn.f COMM
	ifort -o $@  $(COMPFLG) cvemotnnn.f
obj/cvemot1.o   : cvemot1.f COMM
	ifort -o $@  $(COMPFLG) cvemot1.f 
obj/cvemotn1.o   : cvemotn1.f COMM
	ifort -o $@  $(COMPFLG) cvemotn1.f 
obj/gamma.o     : gamma.f COMM
	ifort -o $@  $(COMPFLG) gamma.f
obj/gauss.o    : gauss.f COMM
	ifort -o $@  $(COMPFLG) gauss.f 
obj/newxyz.o    : newxyz.f COMM
	ifort -o $@  $(COMPFLG) newxyz.f 
obj/lath.o     : lath.f COMM
	ifort -o $@  $(COMPFLG) lath.f 
obj/main.o    : main.f COMM param
	ifort -o $@  $(COMPFLG) main.f 
obj/doda.o    : doda.f COMM 
	ifort -o $@  $(COMPFLG) doda.f 
obj/dado.o    : dado.f COMM param
	ifort -o $@  $(COMPFLG) dado.f 
obj/tepchr.o     : tepchr.f COMM
	ifort -o $@  $(COMPFLG) tepchr.f
obj/tepchr2.o     : tepchr2.f COMM
	ifort -o $@  $(COMPFLG) tepchr2.f
obj/corect.o    : corect.f COMM
	ifort -o $@  $(COMPFLG) corect.f 
obj/corect2.o    : corect2.f COMM
	ifort -o $@  $(COMPFLG) corect2.f
obj/fspher.o     : fspher.f COMM
	ifort -o $@  $(COMPFLG) fspher.f
obj/fspher2.o     : fspher2.f COMM
	ifort -o $@  $(COMPFLG) fspher2.f
obj/furim.o    : furim.f COMM
	ifort -o $@  $(COMPFLG) furim.f
obj/printer.o     : printer.f COMM
	ifort -o $@  $(COMPFLG) printer.f
obj/dipolos.o     : dipolos.f COMM
	ifort -o $@  $(COMPFLG) dipolos.f
obj/qemot.o      : qemot.f COMM
	ifort -o $@  $(COMPFLG) qemot.f
obj/qtmotb.o      : qtmotb.f COMM
	ifort -o $@  $(COMPFLG) qtmotb.f
obj/tepq.o      : tepq.f COMM
	ifort -o $@  $(COMPFLG) tepq.f
obj/inicio.o      : inicio.f COMM
	ifort -o $@  $(COMPFLG) inicio.f
obj/fint.o      : fint.f COMM
	ifort -o $@  $(COMPFLG) fint.f
obj/fint2.o      : fint2.f COMM
	ifort -o $@  $(COMPFLG) fint2.f
obj/grilla.o      : grilla.f COMM
	ifort -o $@  $(COMPFLG) grilla.f
obj/optim.o      : optim.f COMM
	ifort -o $@  $(COMPFLG) optim.f
obj/optimq.o      : optimq.f COMM
	ifort -o $@  $(COMPFLG) optimq.f
obj/proton.o      : proton.f COMM
	ifort -o $@  $(COMPFLG) proton.f
obj/int.o    : int.f param
	ifort -o $@  $(COMPFLG) int.f 
obj/int1G.o    : int1G.f param
	ifort -o $@  $(COMPFLG) int1G.f 
obj/drive.o : drive.f COMM param 
	ifort -o $@  $(COMPFLG2) drive.f 
obj/func.o : func.f
	ifort -o $@  $(COMPFLG) func.f 
obj/int2.o    : int2.f param
	ifort -o $@ $(COMPFLG) int2.f
obj/int2G.o    : int2G.f param
	ifort -o $@  $(COMPFLG) int2G.f 
obj/int3G.o    : int3G.f param
	ifort -o $@  $(COMPFLG) $(COMPFLG) int3G.f 
obj/int3.o    : int3.f param
	ifort -o $@ $(COMPFLG) int3.f
obj/int3N.o    : int3N.f param
	ifort -o $@ $(COMPFLG) int3N.f
obj/int3lu.o    : int3lu.f param
	ifort -o $@ $(COMPFLG) int3lu.f
obj/int3mem.o    : int3mem.f param
	ifort -o $@  $(COMPFLG) int3mem.f
obj/exch.o    : exch.f param
	ifort -o $@  $(COMPFLG) exch.f
obj/exch2.o    : exch2.f param
	ifort -o $@  $(COMPFLG) exch2.f
obj/pot.o    : pot.f 
	ifort -o $@  $(COMPFLG) pot.f 
obj/SCF.o    : SCF.f param 
	ifort -o $@  $(COMPFLG) SCF.f
obj/intSG.o     : intSG.f param
	ifort -o $@  $(COMPFLG) intSG.f 
obj/write.o     : write.f param
	ifort -o $@  $(COMPFLG) write.f
obj/nwrite.o     : nwrite.f param
	ifort -o $@  $(COMPFLG) nwrite.f
obj/dns.o     : dns.f param
	ifort -o $@  $(COMPFLG) dns.f 
obj/dns2.o     : dns2.f param
	ifort -o $@  $(COMPFLG) dns2.f 
obj/dnsg.o     : dnsg.f param
	ifort -o $@  $(COMPFLG) dnsg.f 
obj/densg.o     : densg.f param
	ifort -o $@  $(COMPFLG) densg.f 
obj/potg.o     : potg.f 
	ifort -o $@  $(COMPFLG) potg.f 
obj/dip.o     : dip.f param
	ifort -o $@  $(COMPFLG) dip.f 
obj/dipG.o     : dipG.f param
	ifort -o $@  $(COMPFLG) dipG.f 
obj/dip2.o     : dip2.f param
	ifort -o $@  $(COMPFLG) dip2.f 
obj/exchnum.o     : exchnum.f param
	ifort -o $@  $(COMPFLG) exchnum.f 
obj/grid.o     : grid.f
	ifort -o $@  $(COMPFLG) grid.f 
obj/exchnum2.o     : exchnum2.f param
	ifort -o $@  $(COMPFLG) exchnum2.f 
obj/exchfock.o     : exchfock.f param
	ifort -o $@  $(COMPFLG) exchfock.f 
obj/SCFop.o    : SCFop.f param COMM
	ifort -o $@  $(COMPFLG) SCFop.f
obj/densgop.o     : densgop.f param
	ifort -o $@  $(COMPFLG) densgop.f
obj/dns2op.o     : dns2op.f param
	ifort -o $@  $(COMPFLG) dns2op.f
obj/dnsop.o     : dnsop.f param
	ifort -o $@  $(COMPFLG) dnsop.f
obj/exchnum2op.o     : exchnum2op.f param
	ifort -o $@  $(COMPFLG) exchnum2op.f
obj/exchnumop.o     : exchnumop.f param
	ifort -o $@  $(COMPFLG) exchnumop.f
obj/potop.o     : potop.f param
	ifort -o $@  $(COMPFLG) potop.f
obj/potgop.o     : potgop.f 
	ifort -o $@  $(COMPFLG) potgop.f
obj/dnsgop.o     : dnsgop.f param
	ifort -o $@  $(COMPFLG) dnsgop.f
obj/alg.o     : alg.f 
	ifort -o $@  $(COMPFLG) alg.f
obj/eig.o     : eig.f
	ifort -o $@  $(COMPFLG) eig.f
obj/svd.o     : svd.f
	ifort -o $@  $(COMPFLG) svd.f
obj/intsol.o    : intsol.f param
	ifort -o $@  $(COMPFLG) intsol.f 
obj/intsolG.o    : intsolG.f param
	ifort -o $@  $(COMPFLG) intsolG.f 
obj/lalg.o    : lalg.f param
	ifort -o $@  $(COMPFLG) lalg.f 
obj/ungles.o    : ungles.f COMM param
	ifort -o $@  $(COMPFLG) ungles.f 
obj/dipcor.o    : dipcor.f COMM 
	ifort -o $@  $(COMPFLG) dipcor.f
	
# G2G OBJs
g2g_obj/config.o    : config.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G config.f
g2g_obj/corecte.o   : corecte.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G corecte.f
g2g_obj/finte.o     : finte.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G finte.f
g2g_obj/cvemot.o    : cvemot.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G cvemot.f 
g2g_obj/cvemotn.o   : cvemotn.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G cvemotn.f 
g2g_obj/cvemotnn.o  : cvemotnn.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G cvemotnn.f
g2g_obj/cvemotnnn.o  : cvemotnnn.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G cvemotnnn.f
g2g_obj/cvemot1.o   : cvemot1.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G cvemot1.f 
g2g_obj/cvemotn1.o   : cvemotn1.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G cvemotn1.f 
g2g_obj/gamma.o     : gamma.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G gamma.f
g2g_obj/gauss.o    : gauss.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G gauss.f 
g2g_obj/newxyz.o    : newxyz.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G newxyz.f 
g2g_obj/lath.o     : lath.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G lath.f 
g2g_obj/main.o    : main.f COMM param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G main.f 
g2g_obj/doda.o    : doda.f COMM 
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G doda.f 
g2g_obj/dado.o    : dado.f COMM param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dado.f 
g2g_obj/tepchr.o     : tepchr.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G tepchr.f
g2g_obj/tepchr2.o     : tepchr2.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G tepchr2.f
g2g_obj/corect.o    : corect.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G corect.f 
g2g_obj/corect2.o    : corect2.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G corect2.f
g2g_obj/fspher.o     : fspher.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G fspher.f
g2g_obj/fspher2.o     : fspher2.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G fspher2.f
g2g_obj/furim.o    : furim.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G furim.f
g2g_obj/printer.o     : printer.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G printer.f
g2g_obj/dipolos.o     : dipolos.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dipolos.f
g2g_obj/qemot.o      : qemot.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G qemot.f
g2g_obj/qtmotb.o      : qtmotb.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G qtmotb.f
g2g_obj/tepq.o      : tepq.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G tepq.f
g2g_obj/inicio.o      : inicio.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G inicio.f
g2g_obj/fint.o      : fint.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G fint.f
g2g_obj/fint2.o      : fint2.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G fint2.f
g2g_obj/grilla.o      : grilla.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G grilla.f
g2g_obj/optim.o      : optim.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G optim.f
g2g_obj/optimq.o      : optimq.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G optimq.f
g2g_obj/proton.o      : proton.f COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G proton.f
g2g_obj/int.o    : int.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G int.f 
g2g_obj/int1G.o    : int1G.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G int1G.f 
g2g_obj/drive.o : drive.f COMM param 
	ifort -o $@  $(COMPFLG2) $(G2G_CFLAGS) drive.f
g2g_obj/func.o : func.f
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G func.f 
g2g_obj/int2.o    : int2.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) -DG2G int2.f
g2g_obj/int2G.o    : int2G.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G int2G.f 
g2g_obj/int3G.o    : int3G.f param
	ifort -o $@ $(COMPFLG) $(COMPFLG) $(G2G_CFLAGS) -DG2G int3G.f 
g2g_obj/int3.o    : int3.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) -DG2G int3.f
g2g_obj/int3N.o    : int3N.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) -DG2G int3N.f
g2g_obj/int3lu.o    : int3lu.f param
	ifort -o $@ $(COMPFLG) $(G2G_CFLAGS) -DG2G int3lu.f
g2g_obj/int3mem.o    : int3mem.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G int3mem.f
g2g_obj/exch.o    : exch.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G exch.f
g2g_obj/exch2.o    : exch2.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G exch2.f
g2g_obj/pot.o    : pot.f 
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G pot.f 
g2g_obj/SCF.o    : SCF.f param 
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G SCF.f
g2g_obj/intSG.o     : intSG.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G intSG.f 
g2g_obj/write.o     : write.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G write.f
g2g_obj/nwrite.o     : nwrite.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G nwrite.f
g2g_obj/dns.o     : dns.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dns.f 
g2g_obj/dns2.o     : dns2.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dns2.f 
g2g_obj/dnsg.o     : dnsg.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dnsg.f 
g2g_obj/densg.o     : densg.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G densg.f 
g2g_obj/potg.o     : potg.f 
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G potg.f 
g2g_obj/dip.o     : dip.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dip.f 
g2g_obj/dipG.o     : dipG.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dipG.f 
g2g_obj/dip2.o     : dip2.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dip2.f 
g2g_obj/exchnum.o     : exchnum.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G exchnum.f 
g2g_obj/grid.o     : grid.f
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G grid.f 
g2g_obj/exchnum2.o     : exchnum2.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G exchnum2.f 
g2g_obj/exchfock.o     : exchfock.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G exchfock.f 
g2g_obj/SCFop.o    : SCFop.f param COMM
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G SCFop.f
g2g_obj/densgop.o     : densgop.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G densgop.f
g2g_obj/dns2op.o     : dns2op.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dns2op.f
g2g_obj/dnsop.o     : dnsop.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dnsop.f
g2g_obj/exchnum2op.o     : exchnum2op.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G exchnum2op.f
g2g_obj/exchnumop.o     : exchnumop.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G exchnumop.f
g2g_obj/potop.o     : potop.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G potop.f
g2g_obj/potgop.o     : potgop.f 
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G potgop.f
g2g_obj/dnsgop.o     : dnsgop.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dnsgop.f
g2g_obj/alg.o     : alg.f 
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G alg.f
g2g_obj/eig.o     : eig.f
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G eig.f
g2g_obj/svd.o     : svd.f
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G svd.f
g2g_obj/intsol.o    : intsol.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G intsol.f 
g2g_obj/intsolG.o    : intsolG.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G intsolG.f 
g2g_obj/lalg.o    : lalg.f param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G lalg.f 
g2g_obj/ungles.o    : ungles.f COMM param
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G ungles.f 
g2g_obj/dipcor.o    : dipcor.f COMM 
	ifort -o $@  $(COMPFLG) $(G2G_CFLAGS) -DG2G dipcor.f	
