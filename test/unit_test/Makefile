################################################################################
# UNIT TEST MAKEFILE
################################################################################

# Test for Libxc in the library folder
LIBXC = ${LIBXC_LIBS}
ifeq ("$(LIBXC)","")
    $(info >>> WARNING - Libxc not found in library folder. Please install Libxc before run the tests. <<<)
    HAS_LIBXC_INSTALLED := false
else
    HAS_LIBXC_INSTALLED := true
endif

all: libxc_proxy matrix timers

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
	$(EXEC) matrix/matrix_test > /dev/null
	$(EXEC) echo $$?
	$(EXEC) timers/timers_test > /dev/null
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
