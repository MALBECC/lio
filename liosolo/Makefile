ifeq ($(magma),1)
        DEFINE   += -Dmagma
        LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif

liosolo: main.f
	ifort -L ../lioamber/ -llio-g2g -L ../g2g/ -lg2g  -I../lioamber/obj $(LIBS)  main.f -o liosolo
clean:
	rm -rf liosolo
