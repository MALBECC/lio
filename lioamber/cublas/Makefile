# sources and objects
MOD_SRCS=$(wildcard mods/*.f)
GARCHA_SRCS=$(wildcard *.f)

MOD_OBJ=$(MOD_SRCS:%.f=%.o) 
GARCHA_OBJ=$(MOD_OBJ) $(COMMON_SRCS:%.f=%.o) $(GARCHA_SRCS:%.f=%.o)


DEFINE   = -Dpack -fPIC   
CFLAGS   = -g -mp1 -ip -O3 -c -cpp  $(DEFINE) $(PROFILE)
CFLAGS2  = -g -O1 -c -cpp  $(DEFINE) $(PROFILE)
CFLAGS3  = -g -O3 -c -cpp -parallel   $(DEFINE) $(PROFILE)
CPPFLAGS = -P -traditional-cpp 
LIBS=-L/usr/lib -L/usr/lib64 -L$(MKLROOT)/lib/intel64 -L/lib/ -lfortran -I../garcha_mod.mod -I$(MKLROOT)/include $(PROFILE)  -lmkl_lapack95_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl 
## Targets: make ##

all: dirs liblio-g2g.so 

dirs:
	mkdir -p ../obj/garcha_g2g_obj/

liblio-g2g.so: dirs $(GARCHA_OBJ:%.o=obj/garcha_g2g_obj/%.o)
	ifort -shared -fPIC -o $@ $(GARCHA_OBJ:%.o=obj/garcha_g2g_obj/%.o) $(LIBS)


## Targets: clean ##
clean:
	rm -rf ../obj


### Objects


obj/garcha_g2g_obj/%.o: %.f
	ifort -o $@ $(CFLAGS) $(G2G_CFLAGS) $<

