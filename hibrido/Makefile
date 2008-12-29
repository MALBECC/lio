#-*- mode: make -*-

OBJ =     config.o cvemot.o  cvemotn.o  cvemotnn.o cvemotnnn.o\
	   gamma.o lath.o   main.o  doda.o dado.o \
	   gauss.o  tepchr.o  newxyz.o corect.o corect2.o\
           fspher.o fspher2.o  printer.o  dipolos.o \
           qemot.o qtmotb.o  tepq.o tepchr2.o \
           inicio.o  fint.o fint2.o optim.o optimq.o\
           int.o drive.o SCF.o  func.o     \
           int2.o  int3.o int3G.o exch.o    \
           exch2.o pot.o int2G.o int1G.o \
           intSG.o write.o finte.o corecte.o\
           dns.o dns2.o dnsg.o potg.o dip.o dipG.o  \
           int3N.o exchnum.o exchnum2.o grid.o\
           densg.o  exchfock.o cvemotn1.o\
           nwrite.o  SCFop.o\
           densgop.o dns2op.o dnsop.o exchnum2op.o\
           exchnumop.o potop.o alg.o eig.o svd.o \
           dnsgop.o potgop.o  dip2.o \
           intsolG.o intsol.o grilla.o\
           lalg.o int3mem.o int3lu.o ungles.o \
           cvemot1.o  proton.o dipcor.o  furim.o

#COMPFLG =   -O3 -pipe -march=i686 -fomit-frame-pointer -c\
#           -falign-loops=2 -fno-strength-reduce \
#           -ffast-math  -falign-functions=2 -static
DEFINE   = -Dpack
COMPFLG  = -mp1 -ip -O3 -c $(CFLAGS) -cpp $(DEFINE)
COMPFLG2 = -O3     -c $(CFLAGS) -cpp $(DEFINE)
#more or less general
LINKFLG = -lm
CPPFLAGS =  -P -traditional-cpp


CUDA_COMMONDIR := $(CUDA_INSTALL_PATH)/../sdk/common
CUDA_LIBS=-L$(CUDA_INSTALL_PATH)/lib -L$(CUDA_COMMONDIR)/lib -lcuda -lcudart
LIBS= -L/usr/lib -lblas -lstdc++ -L/opt/intel/mkl/9.0/lib/em64t -lmkl_em64t -lguide -lpthread  -L../gpu -lg2g
CPU_LIBS=$(LIBS)
GPU_LIBS=$(LIBS) $(CUDA_LIBS) 
GPU_CFLAGS=-DGPU

ifeq ($(profile),1)
	COMPFLG += -pg
	COMPFLG2 += -pg	
	LINKFLG += -pg
endif

ifeq ($(ultima_cpu),1)
	GPU_CFLAGS += -DULTIMA_CPU
endif

ifeq ($(convergencia_cpu),1)
	GPU_CFLAGS += -DINT3N_CPU
endif

ifeq ($(ultima_gpu),1)
	COMPFLG += -DULTIMA_GPU
	COMPFLG2 += -DULTIMA_GPU
endif

all: hibrido-gpu hibrido-cpu

clean-cpu:
	rm -f cpu_obj/* hibrido-cpu
	
clean-gpu:
	rm -f gpu_obj/* hibrido-gpu	
	
clean: clean-cpu clean-gpu

hibrido-gpu: $(OBJ:%.o=gpu_obj/%.o)
	ifort -o $@ $(OBJ:%.o=gpu_obj/%.o) $(LINKFLG) $(GPU_LIBS) $(EXTRA_FLAGS)
	
hibrido-cpu: $(OBJ:%.o=cpu_obj/%.o)
	ifort -o $@ $(OBJ:%.o=cpu_obj/%.o) $(LINKFLG) $(CPU_LIBS) $(EXTRA_FLAGS)

# CPU OBJs
cpu_obj/config.o    : config.f COMM
	ifort -o $@  $(COMPFLG) config.f
cpu_obj/corecte.o   : corecte.f COMM
	ifort -o $@  $(COMPFLG) corecte.f
cpu_obj/finte.o     : finte.f COMM
	ifort -o $@  $(COMPFLG) finte.f
cpu_obj/cvemot.o    : cvemot.f COMM
	ifort -o $@  $(COMPFLG) cvemot.f 
cpu_obj/cvemotn.o   : cvemotn.f COMM
	ifort -o $@  $(COMPFLG) cvemotn.f 
cpu_obj/cvemotnn.o  : cvemotnn.f COMM
	ifort -o $@  $(COMPFLG) cvemotnn.f
cpu_obj/cvemotnnn.o  : cvemotnnn.f COMM
	ifort -o $@  $(COMPFLG) cvemotnnn.f
cpu_obj/cvemot1.o   : cvemot1.f COMM
	ifort -o $@  $(COMPFLG) cvemot1.f 
cpu_obj/cvemotn1.o   : cvemotn1.f COMM
	ifort -o $@  $(COMPFLG) cvemotn1.f 
cpu_obj/gamma.o     : gamma.f COMM
	ifort -o $@  $(COMPFLG) gamma.f
cpu_obj/gauss.o    : gauss.f COMM
	ifort -o $@  $(COMPFLG) gauss.f 
cpu_obj/newxyz.o    : newxyz.f COMM
	ifort -o $@  $(COMPFLG) newxyz.f 
cpu_obj/lath.o     : lath.f COMM
	ifort -o $@  $(COMPFLG) lath.f 
cpu_obj/main.o    : main.f COMM param
	ifort -o $@  $(COMPFLG) main.f 
cpu_obj/doda.o    : doda.f COMM 
	ifort -o $@  $(COMPFLG) doda.f 
cpu_obj/dado.o    : dado.f COMM param
	ifort -o $@  $(COMPFLG) dado.f 
cpu_obj/tepchr.o     : tepchr.f COMM
	ifort -o $@  $(COMPFLG) tepchr.f
cpu_obj/tepchr2.o     : tepchr2.f COMM
	ifort -o $@  $(COMPFLG) tepchr2.f
cpu_obj/corect.o    : corect.f COMM
	ifort -o $@  $(COMPFLG) corect.f 
cpu_obj/corect2.o    : corect2.f COMM
	ifort -o $@  $(COMPFLG) corect2.f
cpu_obj/fspher.o     : fspher.f COMM
	ifort -o $@  $(COMPFLG) fspher.f
cpu_obj/fspher2.o     : fspher2.f COMM
	ifort -o $@  $(COMPFLG) fspher2.f
cpu_obj/furim.o    : furim.f COMM
	ifort -o $@  $(COMPFLG) furim.f
cpu_obj/printer.o     : printer.f COMM
	ifort -o $@  $(COMPFLG) printer.f
cpu_obj/dipolos.o     : dipolos.f COMM
	ifort -o $@  $(COMPFLG) dipolos.f
cpu_obj/qemot.o      : qemot.f COMM
	ifort -o $@  $(COMPFLG) qemot.f
cpu_obj/qtmotb.o      : qtmotb.f COMM
	ifort -o $@  $(COMPFLG) qtmotb.f
cpu_obj/tepq.o      : tepq.f COMM
	ifort -o $@  $(COMPFLG) tepq.f
cpu_obj/inicio.o      : inicio.f COMM
	ifort -o $@  $(COMPFLG) inicio.f
cpu_obj/fint.o      : fint.f COMM
	ifort -o $@  $(COMPFLG) fint.f
cpu_obj/fint2.o      : fint2.f COMM
	ifort -o $@  $(COMPFLG) fint2.f
cpu_obj/grilla.o      : grilla.f COMM
	ifort -o $@  $(COMPFLG) grilla.f
cpu_obj/optim.o      : optim.f COMM
	ifort -o $@  $(COMPFLG) optim.f
cpu_obj/optimq.o      : optimq.f COMM
	ifort -o $@  $(COMPFLG) optimq.f
cpu_obj/proton.o      : proton.f COMM
	ifort -o $@  $(COMPFLG) proton.f
cpu_obj/int.o    : int.f param
	ifort -o $@  $(COMPFLG) int.f 
cpu_obj/int1G.o    : int1G.f param
	ifort -o $@  $(COMPFLG) int1G.f 
cpu_obj/drive.o : drive.f COMM param 
	ifort -o $@  $(COMPFLG2) drive.f 
cpu_obj/func.o : func.f
	ifort -o $@  $(COMPFLG) func.f 
cpu_obj/int2.o    : int2.f param
	ifort -o $@ $(COMPFLG) int2.f
cpu_obj/int2G.o    : int2G.f param
	ifort -o $@  $(COMPFLG) int2G.f 
cpu_obj/int3G.o    : int3G.f param
	ifort -o $@  $(COMPFLG) $(COMPFLG) int3G.f 
cpu_obj/int3.o    : int3.f param
	ifort -o $@ $(COMPFLG) int3.f
cpu_obj/int3N.o    : int3N.f param
	ifort -o $@ $(COMPFLG) int3N.f
cpu_obj/int3lu.o    : int3lu.f param
	ifort -o $@ $(COMPFLG) int3lu.f
cpu_obj/int3mem.o    : int3mem.f param
	ifort -o $@  $(COMPFLG) int3mem.f
cpu_obj/exch.o    : exch.f param
	ifort -o $@  $(COMPFLG) exch.f
cpu_obj/exch2.o    : exch2.f param
	ifort -o $@  $(COMPFLG) exch2.f
cpu_obj/pot.o    : pot.f 
	ifort -o $@  $(COMPFLG) pot.f 
cpu_obj/SCF.o    : SCF.f param 
	ifort -o $@  $(COMPFLG) SCF.f
cpu_obj/intSG.o     : intSG.f param
	ifort -o $@  $(COMPFLG) intSG.f 
cpu_obj/write.o     : write.f param
	ifort -o $@  $(COMPFLG) write.f
cpu_obj/nwrite.o     : nwrite.f param
	ifort -o $@  $(COMPFLG) nwrite.f
cpu_obj/dns.o     : dns.f param
	ifort -o $@  $(COMPFLG) dns.f 
cpu_obj/dns2.o     : dns2.f param
	ifort -o $@  $(COMPFLG) dns2.f 
cpu_obj/dnsg.o     : dnsg.f param
	ifort -o $@  $(COMPFLG) dnsg.f 
cpu_obj/densg.o     : densg.f param
	ifort -o $@  $(COMPFLG) densg.f 
cpu_obj/potg.o     : potg.f 
	ifort -o $@  $(COMPFLG) potg.f 
cpu_obj/dip.o     : dip.f param
	ifort -o $@  $(COMPFLG) dip.f 
cpu_obj/dipG.o     : dipG.f param
	ifort -o $@  $(COMPFLG) dipG.f 
cpu_obj/dip2.o     : dip2.f param
	ifort -o $@  $(COMPFLG) dip2.f 
cpu_obj/exchnum.o     : exchnum.f param
	ifort -o $@  $(COMPFLG) exchnum.f 
cpu_obj/grid.o     : grid.f
	ifort -o $@  $(COMPFLG) grid.f 
cpu_obj/exchnum2.o     : exchnum2.f param
	ifort -o $@  $(COMPFLG) exchnum2.f 
cpu_obj/exchfock.o     : exchfock.f param
	ifort -o $@  $(COMPFLG) exchfock.f 
cpu_obj/SCFop.o    : SCFop.f param COMM
	ifort -o $@  $(COMPFLG) SCFop.f
cpu_obj/densgop.o     : densgop.f param
	ifort -o $@  $(COMPFLG) densgop.f
cpu_obj/dns2op.o     : dns2op.f param
	ifort -o $@  $(COMPFLG) dns2op.f
cpu_obj/dnsop.o     : dnsop.f param
	ifort -o $@  $(COMPFLG) dnsop.f
cpu_obj/exchnum2op.o     : exchnum2op.f param
	ifort -o $@  $(COMPFLG) exchnum2op.f
cpu_obj/exchnumop.o     : exchnumop.f param
	ifort -o $@  $(COMPFLG) exchnumop.f
cpu_obj/potop.o     : potop.f param
	ifort -o $@  $(COMPFLG) potop.f
cpu_obj/potgop.o     : potgop.f 
	ifort -o $@  $(COMPFLG) potgop.f
cpu_obj/dnsgop.o     : dnsgop.f param
	ifort -o $@  $(COMPFLG) dnsgop.f
cpu_obj/alg.o     : alg.f 
	ifort -o $@  $(COMPFLG) alg.f
cpu_obj/eig.o     : eig.f
	ifort -o $@  $(COMPFLG) eig.f
cpu_obj/svd.o     : svd.f
	ifort -o $@  $(COMPFLG) svd.f
cpu_obj/intsol.o    : intsol.f param
	ifort -o $@  $(COMPFLG) intsol.f 
cpu_obj/intsolG.o    : intsolG.f param
	ifort -o $@  $(COMPFLG) intsolG.f 
cpu_obj/lalg.o    : lalg.f param
	ifort -o $@  $(COMPFLG) lalg.f 
cpu_obj/ungles.o    : ungles.f COMM param
	ifort -o $@  $(COMPFLG) ungles.f 
cpu_obj/dipcor.o    : dipcor.f COMM 
	ifort -o $@  $(COMPFLG) dipcor.f
	
# GPU OBJs
gpu_obj/config.o    : config.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU config.f
gpu_obj/corecte.o   : corecte.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU corecte.f
gpu_obj/finte.o     : finte.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU finte.f
gpu_obj/cvemot.o    : cvemot.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU cvemot.f 
gpu_obj/cvemotn.o   : cvemotn.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU cvemotn.f 
gpu_obj/cvemotnn.o  : cvemotnn.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU cvemotnn.f
gpu_obj/cvemotnnn.o  : cvemotnnn.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU cvemotnnn.f
gpu_obj/cvemot1.o   : cvemot1.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU cvemot1.f 
gpu_obj/cvemotn1.o   : cvemotn1.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU cvemotn1.f 
gpu_obj/gamma.o     : gamma.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU gamma.f
gpu_obj/gauss.o    : gauss.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU gauss.f 
gpu_obj/newxyz.o    : newxyz.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU newxyz.f 
gpu_obj/lath.o     : lath.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU lath.f 
gpu_obj/main.o    : main.f COMM param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU main.f 
gpu_obj/doda.o    : doda.f COMM 
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU doda.f 
gpu_obj/dado.o    : dado.f COMM param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dado.f 
gpu_obj/tepchr.o     : tepchr.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU tepchr.f
gpu_obj/tepchr2.o     : tepchr2.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU tepchr2.f
gpu_obj/corect.o    : corect.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU corect.f 
gpu_obj/corect2.o    : corect2.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU corect2.f
gpu_obj/fspher.o     : fspher.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU fspher.f
gpu_obj/fspher2.o     : fspher2.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU fspher2.f
gpu_obj/furim.o    : furim.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU furim.f
gpu_obj/printer.o     : printer.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU printer.f
gpu_obj/dipolos.o     : dipolos.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dipolos.f
gpu_obj/qemot.o      : qemot.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU qemot.f
gpu_obj/qtmotb.o      : qtmotb.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU qtmotb.f
gpu_obj/tepq.o      : tepq.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU tepq.f
gpu_obj/inicio.o      : inicio.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU inicio.f
gpu_obj/fint.o      : fint.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU fint.f
gpu_obj/fint2.o      : fint2.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU fint2.f
gpu_obj/grilla.o      : grilla.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU grilla.f
gpu_obj/optim.o      : optim.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU optim.f
gpu_obj/optimq.o      : optimq.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU optimq.f
gpu_obj/proton.o      : proton.f COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU proton.f
gpu_obj/int.o    : int.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU int.f 
gpu_obj/int1G.o    : int1G.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU int1G.f 
gpu_obj/drive.o : drive.f COMM param 
	ifort -o $@  $(COMPFLG2) $(GPU_CFLAGS) drive.f
gpu_obj/func.o : func.f
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU func.f 
gpu_obj/int2.o    : int2.f param
	ifort -o $@ $(COMPFLG) $(GPU_CFLAGS) -DGPU int2.f
gpu_obj/int2G.o    : int2G.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU int2G.f 
gpu_obj/int3G.o    : int3G.f param
	ifort -o $@ $(COMPFLG) $(COMPFLG) $(GPU_CFLAGS) -DGPU int3G.f 
gpu_obj/int3.o    : int3.f param
	ifort -o $@ $(COMPFLG) $(GPU_CFLAGS) -DGPU int3.f
gpu_obj/int3N.o    : int3N.f param
	ifort -o $@ $(COMPFLG) $(GPU_CFLAGS) -DGPU int3N.f
gpu_obj/int3lu.o    : int3lu.f param
	ifort -o $@ $(COMPFLG) $(GPU_CFLAGS) -DGPU int3lu.f
gpu_obj/int3mem.o    : int3mem.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU int3mem.f
gpu_obj/exch.o    : exch.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU exch.f
gpu_obj/exch2.o    : exch2.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU exch2.f
gpu_obj/pot.o    : pot.f 
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU pot.f 
gpu_obj/SCF.o    : SCF.f param 
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU SCF.f
gpu_obj/intSG.o     : intSG.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU intSG.f 
gpu_obj/write.o     : write.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU write.f
gpu_obj/nwrite.o     : nwrite.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU nwrite.f
gpu_obj/dns.o     : dns.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dns.f 
gpu_obj/dns2.o     : dns2.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dns2.f 
gpu_obj/dnsg.o     : dnsg.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dnsg.f 
gpu_obj/densg.o     : densg.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU densg.f 
gpu_obj/potg.o     : potg.f 
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU potg.f 
gpu_obj/dip.o     : dip.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dip.f 
gpu_obj/dipG.o     : dipG.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dipG.f 
gpu_obj/dip2.o     : dip2.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dip2.f 
gpu_obj/exchnum.o     : exchnum.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU exchnum.f 
gpu_obj/grid.o     : grid.f
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU grid.f 
gpu_obj/exchnum2.o     : exchnum2.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU exchnum2.f 
gpu_obj/exchfock.o     : exchfock.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU exchfock.f 
gpu_obj/SCFop.o    : SCFop.f param COMM
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU SCFop.f
gpu_obj/densgop.o     : densgop.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU densgop.f
gpu_obj/dns2op.o     : dns2op.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dns2op.f
gpu_obj/dnsop.o     : dnsop.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dnsop.f
gpu_obj/exchnum2op.o     : exchnum2op.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU exchnum2op.f
gpu_obj/exchnumop.o     : exchnumop.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU exchnumop.f
gpu_obj/potop.o     : potop.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU potop.f
gpu_obj/potgop.o     : potgop.f 
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU potgop.f
gpu_obj/dnsgop.o     : dnsgop.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dnsgop.f
gpu_obj/alg.o     : alg.f 
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU alg.f
gpu_obj/eig.o     : eig.f
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU eig.f
gpu_obj/svd.o     : svd.f
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU svd.f
gpu_obj/intsol.o    : intsol.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU intsol.f 
gpu_obj/intsolG.o    : intsolG.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU intsolG.f 
gpu_obj/lalg.o    : lalg.f param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU lalg.f 
gpu_obj/ungles.o    : ungles.f COMM param
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU ungles.f 
gpu_obj/dipcor.o    : dipcor.f COMM 
	ifort -o $@  $(COMPFLG) $(GPU_CFLAGS) -DGPU dipcor.f	
