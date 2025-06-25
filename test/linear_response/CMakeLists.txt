# test/CMakeLists.txt

set(TEST_NAME "agua")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i agua.in -c agua.xyz -v ")
set(TEST_ENV "GFORTRAN_UNBUFFERED_ALL=1")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "ESAs_ESAfosc")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i agua.in -c agua.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "indol")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i indol.in -c indol.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "truncatedLR_estda")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i co.in -c co.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

if(USE_LIBXC_CPU)
    set(TEST_NAME "ESAs_secondLR")
    set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
    set(LIO_ARGS "-i p3ht-2.in -c p3ht-2.xyz -v ")
    add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

    set(TEST_NAME "truncatedLR_delete-mos")
    set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
    set(LIO_ARGS "-i indol.in -c indol.xyz -v ")
    add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

    set(TEST_NAME "truncatedLR_reduced-space")
    set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
    set(LIO_ARGS "-i form_wat.in -c form_wat.xyz -v ")
    add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})
endif()
