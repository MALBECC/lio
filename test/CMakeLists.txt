# test/CMakeLists.txt

find_package(Python)

# DEFINE FUNCTION TO RUN TESTS
set(LIOSOLO_EXECUTABLE $<TARGET_FILE:liosolo>)
set(G2G_PATH $<TARGET_FILE:g2g>)
set(LIO-G2G_PATH $<TARGET_FILE:lio-g2g>)
function(add_lio_test _test_name _workdir _lio_args _env)

    if(TEST_CHECKS_ONLY)
        add_test(NAME "${_test_name}" 
            COMMAND bash -c "
                ${Python_EXECUTABLE} check_test.py > py_test_result.out && \
                ${CMAKE_COMMAND} -E compare_files --ignore-eol py_test_result.out py_test_result.ok
            "
            WORKING_DIRECTORY ${_workdir}
        )
    else()
        add_test(NAME "${_test_name}" 
            COMMAND bash -c "
                ${LIOSOLO_EXECUTABLE} ${_lio_args} > output && \
                ${Python_EXECUTABLE} check_test.py > py_test_result.out && \
                ${CMAKE_COMMAND} -E compare_files --ignore-eol py_test_result.out py_test_result.ok
            "
            WORKING_DIRECTORY ${_workdir}
        )
    endif()
    set_tests_properties("${_test_name}" PROPERTIES 
        FAIL_REGULAR_EXPRESSION "ERROR"
        DEPENDS ${LIOSOLO_EXECUTABLE} ${G2G_PATH} ${LIO-G2G_PATH}
        ENVIRONMENT "LIOHOME=${CMAKE_SOURCE_DIR}; ${_env}"
    )
endfunction()

# RUN TESTS
add_subdirectory(base_lio)
if(USE_LIBINT)
    add_subdirectory(external_functionals)
    if(USE_LIBXC_CPU OR USE_LIBXC_GPU)
        add_subdirectory(linear_response)
    endif()
endif()