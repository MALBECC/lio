################################################################################
# BASIC SETUPS
PROGRAM := inputgen.x
.DEFAULT_GOAL := inputgen.x

# COMPILATION SETUP
FC     = gfortran
FFLAGS =

# LINK SETUP
LC     = gfortran
LIBS   =
LFLAGS =

vpath %.f90 src

################################################################################
# LIST OF OBJECTS AND DEPENDENCIES
OBJECTS_L0 :=
OBJECTS_L1 :=
OBJECTS_L2 :=

OBJECTS_L0 += keywords.o
OBJECTS_L0 += auxmod_subs.o

OBJECTS_L1 += parser_gaussian.o
OBJECTS_L1 += parser_lio.o

OBJECTS_L2 += data_control.o

OBJECTS := $(OBJECTS_L2) $(OBJECTS_L1) $(OBJECTS_L0)



$(OBJECTS_L2) : $(OBJECTS_L1)
$(OBJECTS_L1) : $(OBJECTS_L0)

inputgen.o: $(OBJECTS)

################################################################################
# RULES SECTION

%.mod : %.o
	@touch $@

%.o :  %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o :  %.F90
	$(FC) $(FFLAGS) -c $< -o $@

inputgen.x : inputgen.o $(OBJECTS)
	$(LC) $^ $(LFLAGS) $(LIBS) -o $@

.PHONY: clean
clean:
	rm -rf inputgen.x *.o *.mod

################################################################################
