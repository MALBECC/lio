#
#  Makefile for fdf example and fdf module.
#
include ../../Makefile.arch
#
# module must include io.f compilation for hybrid
#
module:  libfdf.a
libfdf.a:  fdf.o fdf_mod.o parse.o io.o 
	ar qv libfdf.a fdf.o fdf_mod.o parse.o io.o 
	-$(RANLIB) libfdf.a
	cp libfdf.a ..
	@cp -p *.*d ..
	@echo "fdf/work.pc" >> ../work.pcl



#
# Note trick for Intel Fortran Compiler
#
fdf_mod.o: fdf_mod.f parse.o
fdf.o:		fdf.h
#
#---------------------------------------------------
OBJS=	sample.o fdf.o io.o
#
sample: $(OBJS)
	$(FC) $(LDFLAGS) -o sample $(OBJS)
#
check:	sample sample.fdf
	sample < sample.fdf
#
io_sample: io_sample.o io.o
	$(FC)  $(LDFLAGS) -o io_sample io_sample.o io.o
#
io_check: io_sample
	io_sample
clean:
	rm -f *.o libfdf.a *.mod fdf_mod.o io_sample.o sample io_sample
	rm -f OUTPUT FDF.debug file2 file3 sample.out
	rm -f *.pc *.pcl
#
sample.o:	fdfdefs.h



