ifeq ($(gcc),1)
  CXX = gfortran
else
  CXX = ifort
endif 

ifeq ($(magma),1)
        DEFINE   += -Dmagma
        LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif

 liosolo:
	$(CXX) main.f -L ../lioamber/ -llio-g2g -L ../g2g/ -lg2g -I../lioamber/ $(LIBS) -o liosolo
 clean:
	rm -f liosolo
