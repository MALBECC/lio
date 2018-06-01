################################################################################
# UNIT TEST MAKEFILE
################################################################################

# Test for Libxc in the library folder
LIBXC ?= $(shell ldconfig -p | grep -e libxc*.so)
ifeq ("$(LIBXC)","")
    $(info >>> WARNING - Libxc not found in library folder. Please install the library to run the tests. <<<)
endif


all: libxc_proxy matrix timers

.PHONY: accumulate_point
accumulate_point:
	$(MAKE) -C accumulate_point


.PHONY: iteration
iteration:
	$(MAKE) -C iteration


.PHONY: libxc_proxy
libxc_proxy:
	$(MAKE) -C libxc_proxy


.PHONY: matrix
matrix:
	$(MAKE) -C matrix


.PHONY: timers
timers:
	$(MAKE) -C timers

.PHONY: test
test: test_gpu test_cpu

.PHONY: test_gpu
test_gpu:
	$(EXEC) libxc_proxy/main_gpu > /dev/null
	$(EXEC) echo $$?


.PHONY: test_cpu
test_cpu:
	$(EXEC) libxc_proxy/main_cpu > /dev/null
	$(EXEC) echo $$?


.PHONY: clean
clean:
	$(MAKE) clean -C libxc_proxy
	$(MAKE) clean -C matrix
	$(MAKE) clean -C timers

################################################################################
