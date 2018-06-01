################################################################################
# UNIT TEST MAKEFILE
################################################################################

# Test for Libxc in the library folder
LIBXC ?= $(shell ldconfig -p | grep -e libxc*.so)
ifeq ("$(LIBXC)","")
    $(info >>> WARNING - Libxc not found in library folder. Please install the library to run the tests. <<<)
endif


all: accumulate_point iteration

.PHONY: accumulate_point
accumulate_point:
	$(MAKE) -C accumulate_point


.PHONY: iteration
iteration:
	$(MAKE) -C iteration


.PHONY: test
test: test_gpu


.PHONY: test_gpu
test_gpu:
	$(EXEC) accumulate_point/main_gpu > /dev/null
	$(EXEC) echo $$?
	$(EXEC) iteration/main_gpu > /dev/null
	$(EXEC) echo $$?


.PHONY: test_cpu
test_cpu:
	$(EXEC) accumulate_point/main_cpu > /dev/null
	$(EXEC) echo $$?


.PHONY: clean
clean:
	$(MAKE) clean -C accumulate_point
	$(MAKE) clean -C iteration

################################################################################
