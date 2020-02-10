################################################################################
# UNIT TEST MAKEFILE
################################################################################

#Flag for compile the tests that use Libxc
result := 999

HAS_LIBXC_INSTALLED := false

# Test for Libxc in the library folder
LIBXC = ${LIBXC_LIBS}
ifeq ("$(LIBXC)","")
    $(info >>> WARNING - Libxc not found in library folder. Please install Libxc before run the tests. <<<)
    HAS_LIBXC_INSTALLED := false
else
    HAS_LIBXC_INSTALLED := true
endif

all: accumulate_point iteration test

.PHONY: accumulate_point
accumulate_point:
	$(MAKE) -C accumulate_point


.PHONY: iteration
iteration:
	$(MAKE) -C iteration


.PHONY: test
test: test_gpu test_cpu


.PHONY: test_gpu
test_gpu:
	@echo "Running gpu test..."
ifeq ($(HAS_LIBXC_INSTALLED),true)
	$(eval result:=$(shell accumulate_point/main_gpu > /dev/null 2>&1 ; echo $$?))
	@echo "gpu integration-test01: " $(result)
	$(eval result:=$(shell accumulate_point/main_gpu > /dev/null 2>&1 ; echo $$?))
	@echo "gpu integration-test02: " $(result)
	$(eval result:=$(shell iteration/main_gpu > /dev/null 2>&1 ; echo $$?))
	@echo "gpu integration-test03: " $(result)
else
	@echo ">>> WARNING: Install Libxc before run the integration with Libxc-GPU"
endif


.PHONY: test_cpu
test_cpu:
	@echo "Running cpu test..."
ifeq ($(HAS_LIBXC_INSTALLED),true)
	$(eval result=$(shell accumulate_point/main_cpu > /dev/null 2>&1 ; echo $$?))
	@echo "cpu integration-test: " $(result)
else
	@echo ">>> WARNING: Install Libxc before run the integration with Libxc-CPU"
endif


.PHONY: clean
clean:
	$(MAKE) clean -C accumulate_point
	$(MAKE) clean -C iteration

################################################################################
