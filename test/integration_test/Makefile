################################################################################
# UNIT TEST MAKEFILE
################################################################################

#Flag for compile the tests that use Libxc
HAS_LIBXC_INSTALLED := false

# Test for Libxc in the library folder
LIBXC = ${LIBXC_LIBS}
ifeq ("$(LIBXC)","")
    $(info >>> WARNING - Libxc not found in library folder. Please install Libxc before run the tests. <<<)
    HAS_LIBXC_INSTALLED := false
else
    HAS_LIBXC_INSTALLED := true
endif


all: accumulate_point iteration

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
ifeq ($(HAS_LIBXC_INSTALLED),true)
	$(EXEC) accumulate_point/main_gpu > /dev/null
	$(EXEC) echo $$?
	$(EXEC) iteration/main_gpu > /dev/null
	$(EXEC) echo $$?
else
	@echo ">>> WARNING: Install Libxc before run the integration with Libxc-GPU"
endif

.PHONY: test_cpu
test_cpu:
ifeq ($(HAS_LIBXC_INSTALLED),true)
	$(EXEC) accumulate_point/main_cpu > /dev/null
	$(EXEC) echo $$?
else
	@echo ">>> WARNING: Install Libxc before run the integration with Libxc-CPU"
endif


.PHONY: clean
clean:
	$(MAKE) clean -C accumulate_point
	$(MAKE) clean -C iteration

################################################################################
