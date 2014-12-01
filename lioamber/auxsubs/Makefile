#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# COMPILER AND FLAGS
#--------------------------------------------------------------------#
FC     = ifort
FFLAGS = -g -O3 -c -cpp -Dpack -fPIC
#
#
# DEPENDENCIES
#--------------------------------------------------------------------#
INCLUDES  =
INCLUDES += gaussbell.f gaussbell_h.f
INCLUDES += commutate.f commutate_h.f
INCLUDES += matmulnanox.f matmulnanox_h.f
INCLUDES += transbase.f transbase_h.f
#
#
# COMMANDS
#--------------------------------------------------------------------#
all : auxsubs.o
	cp auxsubs.mod ../.
	cp auxsubs.o ../obj/garcha_g2g_obj/.

auxsubs.o : $(INCLUDES)

test : testprog.x

testprog.x : testprog.o auxsubs.o
	$(FC) $(FFLAGS) auxsubs.o testprog.o -o testprog.x

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: clean
clean :
	rm -f auxsubs.o auxsubs.mod testprog.o testprog.x

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
