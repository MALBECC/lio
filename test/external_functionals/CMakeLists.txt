# test/CMakeLists.txt
include(CMakeParseArguments) # For cmake_parse_arguments


find_package(Python)
set(LIOSOLO_EXECUTABLE $<TARGET_FILE:liosolo>)
set(G2G_PATH $<TARGET_FILE:g2g>)
set(LIO-G2G_PATH $<TARGET_FILE:lio-g2g>)

function(add_lio_test _test_name _workdir _lio_args _env)

    add_test(NAME "${_test_name}" 
        COMMAND bash -c "
            ${LIOSOLO_EXECUTABLE} ${_lio_args} > output && \
            ${Python_EXECUTABLE} check_test.py > py_test_result.txt && \
            ${CMAKE_COMMAND} -E compare_files --ignore-eol py_test_result.txt py_test_result.ok
        "
        WORKING_DIRECTORY ${_workdir}
    )
    set_tests_properties("${_test_name}" PROPERTIES 
        FAIL_REGULAR_EXPRESSION "ERROR"
        DEPENDS ${LIOSOLO_EXECUTABLE} ${G2G_PATH} ${LIO-G2G_PATH}
        ENVIRONMENT "LIOHOME=${CMAKE_SOURCE_DIR}; ${_env}"
    )
endfunction()

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