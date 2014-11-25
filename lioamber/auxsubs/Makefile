#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# README
#--------------------------------------------------------------------#
#
#   This is the makefile for the auxsubs module, which
#   contains different kinds of miscelaneous auxilliary
#   subroutines and functions used throughout the main
#   program lio.
#
#
#   Steps to add a new file containing procedures:
#--------------------------------------------------------------------#
#
#   (1) Create a file with the code of the procedures and
#       save it in this containing folder.
#
#   (2) The file "auxsubs.f" has a section for including
#       all procedure-containing files. Add it there using
#       the following syntax:
#
#      >       include 'file.f'
#
#   (3) Add the new file as a dependency of auxsubs by
#       concatenating it in the INCLUDES variable inside
#       this file.
#
#
#   Steps to group procedures inside an interface:
#--------------------------------------------------------------------#
#
#   If one has many variations of one procedure, one can
#   group all of them together inside an interface by 
#   creating another file where the interface is specified
#   and then including this file in the section for headers
#   using the same syntax as before. This file also needs
#   to be incorporated as a dependency by concatenating
#   it to the INCLUDES variable.
#
#   One has to remember to be specially carefull when 
#   introducing the parameters of a subroutine that is
#   overloaded (part of an interface). The most advisable
#   thing to do is to always use variables to pass the
#   arguments instead of doing so directly.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#
# DEPENDENCIES
#--------------------------------------------------------------------#
INCLUDES = \
gaussbell.f \
gaussbell_h.f \
commutate.f \
commutate_h.f
#
#
# COMPILER AND FLAGS
#--------------------------------------------------------------------#
FC     = ifort
FFLAGS = -O

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
