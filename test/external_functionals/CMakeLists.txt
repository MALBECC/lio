# test/CMakeLists.txt
#include(CMakeParseArguments) # For cmake_parse_arguments

set(TEST_NAME "pbe0")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i agua.in -b basis -c agua.xyz -v ")
set(TEST_ENV "GFORTRAN_UNBUFFERED_ALL=1")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

if(USE_LIBXC_CPU)
    set(TEST_NAME "cam_b3lyp")
    set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
    set(LIO_ARGS "-i agua.in -b basis -c agua.xyz -v ")
    add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})
endif()