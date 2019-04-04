################################################################################
# UNIT TEST MAKEFILE
################################################################################

# Test result
result:=1

# Test for Libxc in the library folder
LIBXC = ${LIBXC_LIBS}
ifeq ("$(LIBXC)","")
    $(info >>> WARNING - Libxc not found in library folder. Please install Libxc before run the tests. <<<)
    HAS_LIBXC_INSTALLED := false
else
    HAS_LIBXC_INSTALLED := true
endif

all: libxc_proxy matrix timers test

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
	@echo "Running gpu test..."
ifeq ($(HAS_LIBXC_INSTALLED),true)
	$(eval result:=$(shell libxc_proxy/main_gpu > /dev/null 2>&1 ; echo $$?))
	@echo "gpu unit-test01: " $(result)
	$(eval result:=$(shell matrix/matrix_test > /dev/null 2>&1 ; echo $$?))
	@echo "gpu unit-test02: " $(result)
	$(eval result:=$(shell timers/timers_test > /dev/null 2>&1 ; echo $$?))
	@echo "gpu unit-test03: " $(result)
else
	@echo ">>> WARNING: Install Libxc before run the integration with Libxc-GPU"
endif


.PHONY: test_cpu
	@echo "Running cpu test..."
test_cpu:
ifeq ($(HAS_LIBXC_INSTALLED),true)
	$(eval result:=$(shell libxc_proxy/main_cpu > /dev/null 2>&1 ; echo $$?))
	@echo "cpu integration-test: " $(result)
else
	@echo ">>> WARNING: Install Libxc before run the integration with Libxc-CPU"
endif

.PHONY: check_test
check_test:
# Check test results
ifeq ($(a),0)
	@echo "Test ok"
else
	@echo "Test fail"
endif


.PHONY: clean
clean:
	$(MAKE) clean -C libxc_proxy
	$(MAKE) clean -C matrix
	$(MAKE) clean -C timers

################################################################################
