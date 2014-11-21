#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# COMPILER AND FLAGS
#--------------------------------------------------------------------#
FC     = ifort
FFLAGS = -O

#
# DEPENDENCIES
#--------------------------------------------------------------------#
INCLUDES = \
gaussbell.f \
gaussbell_h.f \

all : auxsubs.o
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
