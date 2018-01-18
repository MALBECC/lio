.DEFAULT_GOAL = all
CC = ifort

SOURCES := $(wildcard *.f90)
TARGETS := $(SOURCES:%.f90=%)
CLTARGETS := $(TARGETS:%=../bin/%)


all: $(TARGETS)
%: %.f90 
	$(CC) -o ../bin/$@ $<

clean:
	rm -r $(CLTARGETS)

