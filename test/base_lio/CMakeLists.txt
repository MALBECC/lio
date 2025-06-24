# test/CMakeLists.txt
include(CMakeParseArguments) # For cmake_parse_arguments


find_package(Python)
set(LIOSOLO_EXECUTABLE $<TARGET_FILE:liosolo>)
set(G2G_PATH $<TARGET_FILE:g2g>)
set(LIO-G2G_PATH $<TARGET_FILE:lio-g2g>)

function(add_lio_test _test_name _workdir _lio_args)

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
        ENVIRONMENT "LIOHOME=${CMAKE_SOURCE_DIR}"
    )
endfunction()

set(TEST_NAME "agua")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i agua.in -b basis -c agua.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "OxyMol")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i Oxy-mol.in -c Oxy-mol.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "Fe3H2O6")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i fe3h2o.in -c fe3h2o.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "fosfatoQMMM")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i fos.in -c fos.xyz -b basis -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "ECP")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i cl2.in -c cl2.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "ECP2")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i fe3h2o.in -c fe3h2o.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "TDDFTField")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i chloride.in -c chloride.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "QMinPcharges")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i agua.in -b basis -c agua.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "TDDFTHCL")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i chloride.in -c chloride.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

# set(TEST_NAME "CDFT")
# set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
# set(LIO_ARGS "-i cdft.in -c poliac.xyz -v ")
# add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "TBDFT")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS " -i tbdft.in -c 10H.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "TB-DLVN")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS " -i tbdft.in -c 10H.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "LinearSearchRho")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i fe3h2o.in -b DZVP -c fe3h2o.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

# set(TEST_NAME "Zn_timers")
# set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
# set(LIO_ARGS "-i multiZn.in -b basis -c multiZn.xyz -v ")
# add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "Heme")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i heme.in -b DZVP -c heme.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "TrpDFTD")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i trp.in -c trp.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

set(TEST_NAME "AoD_test")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i aod.in -c HSOH.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS})

