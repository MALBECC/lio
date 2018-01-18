# Makefile for HYBRID 
#

default: what hybrid
.SUFFIXES: .f .F .o .a  .f90 .F90
#
include ../Makefile.arch
#
#AGREGAR ESTA LINEA SI DA ERROR DE M2C
#%.o: %.mod
#------------------- 
what:
	@echo 
	@echo "Compilation options used: "
	@echo "Compiler : ${FC}"
	@echo "Analytics: ${analytics}"
	@echo


# Note that machine-specific files are now in top Src directory.
OBJS =  chkdim.o conjgr.o fixed.o iofa.o ofc.o functions.o \
	inver.o ioxv.o paste.o reclat.o \
	reinit.o timestamp.o volcel.o cgvc.o iocg.o \
	centerall.o ran3.o dot.o readall.o subconstr.o \
	writeall.o qmmm_lst_blk.o solv_assign.o \
	ljef.o solv_ene_fce.o linkatom.o recipes.o \
	assignation.o center_rotation.o hybrid.o \

# FDF Module
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) module)

# Routines using fdf calls.
ioxv.o iofa.o : $(FDF)
cgvc.o iocg.o ioxvconstr.o ofc.o : $(FDF)
fixed.o reinit.o readall.o writeall.o : $(FDF)
subconstr.o : $(FDF)
linkatom.o solv_assign.o : $(FDF)
asign.o : $(FDF)

# Rest of Modules
# Note : precision must be the first module
MOD_OBJS=precision.o \
	 ionew.o \
	 parallel.o \
	 sys.o \
	 fsiesta.o \
	 version.o \
	 neutralatom.o \
         parsing.o \


COM_OBJS=$(OBJS) $(SYSOBJ)
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)

# Routines using modules
solv_assign.o ioxv.o : precision.o
qmmm_lst_blk.o readall.o cgvc.o chkdim.o : precision.o
iofa.o fixed.o ioxvconstr.o assign.o hybrid.o : precision.o
ionew.o chkdim.o solv_assign.o hybrid.o readall.o qmmm_lst_blk.o fixed.o assign.o : sys.o
subconstr.o :sys.o
hybrid.o : fsiesta.o
hybrid.o : neutralatom.o
neutralatom.o hybrid.o solv_assign.o :ionew.o
reinit.o readall.o writeall.o : ionew.o
qmmm_lst_blk.o : ionew.o
iofa.o iocg.o ioxv.o ioxvconstr.o ofc.o : ionew.o
assign.o : parse.o


version:
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed  "s/SIESTA_ARCH/${SIESTA_ARCH}/g" version.f > temp0.f
	@sed  "s/FFLAGS/${FC} ${FFLAGS}/g" temp0.f > compinfo.F90
	@($(MAKE) -s "FC=$(FC)" "FFLAGS=$(FFLAGS)" compinfo.o)
	@rm -f temp0.f compinfo.F90
	@mv compinfo.o version.o
	@echo

# Lio objects
LIO_PATH = ../..
#cambiar esta direccion cuando hybrid este dentro de lio

LIO_H = $(wildcard $(LIO_PATH)/lioamber/obj/*.o) 
G2G = -lg2g
LFLAG= -L$(LIO_PATH)/g2g

hybrid: version $(FDF) $(ALL_OBJS)
	@mkdir -p ../bin
	@echo $(LIO_H)
	@echo $(LFLAG)
	$(FC) -o ../bin/hybrid \
		$(FFLAGS) $(LFLAG) $(G2G) $(ALL_OBJS) $(LIO_H) $(FDF) $(LIBS)

clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f *.o  *.a *.pcl *.pc fdf.Standard
	rm -f *.mod
	rm -f ../bin/hybrid
	(cd fdf ; $(MAKE) clean)

