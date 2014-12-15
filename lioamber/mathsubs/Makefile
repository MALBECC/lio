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
INCLUDES += transbase.f transbase_h.f
#
#
# MODULE COMPILATION COMMANDS
#--------------------------------------------------------------------#
all : mathsubs.o
	cp mathsubs.mod ../.
	cp mathsubs.o ../obj/.


mathsubs.o : $(INCLUDES)

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@
#
# TEST
#--------------------------------------------------------------------#
.PHONY: test
TOFLAGS = -g
TXFLAGS = -g
test : test.x

test.x : mathsubs.o test.o
	$(FC) $(TXFLAGS) $^ -o $@
test.o : test_mathsubs.f mathsubs.o
	$(FC) $(TOFLAGS) -c $< -o $@
#
# CLEAN
#--------------------------------------------------------------------#
.PHONY: clean
clean :
	rm -f mathsubs.o mathsubs.mod
	rm -f test.o test.x
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
