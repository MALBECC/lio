# test/CMakeLists.txt
# include(CMakeParseArguments) # For cmake_parse_arguments

set(TEST_ENV "DUMMY=1")

set(TEST_NAME "agua")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i agua.in -b basis -c agua.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "OxyMol")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i Oxy-mol.in -c Oxy-mol.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "Fe3H2O6")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i fe3h2o.in -c fe3h2o.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "fosfatoQMMM")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i fos.in -c fos.xyz -b basis -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "ECP")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i cl2.in -c cl2.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "ECP2")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i fe3h2o.in -c fe3h2o.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "TDDFTField")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i chloride.in -c chloride.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "QMinPcharges")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i agua.in -b basis -c agua.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "TDDFTHCL")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i chloride.in -c chloride.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

# set(TEST_NAME "CDFT")
# set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
# set(LIO_ARGS "-i cdft.in -c poliac.xyz -v ")
# add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARG ${TEST_ENV}S})

set(TEST_NAME "TBDFT")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS " -i tbdft.in -c 10H.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "TB-DLVN")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS " -i tbdft.in -c 10H.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "LinearSearchRho")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i fe3h2o.in -b DZVP -c fe3h2o.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

# set(TEST_NAME "Zn_timers")
# set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
# set(LIO_ARGS "-i multiZn.in -b basis -c multiZn.xyz -v ")
# add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARG ${TEST_ENV}S})

set(TEST_NAME "Heme")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i heme.in -b DZVP -c heme.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "TrpDFTD")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i trp.in -c trp.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

set(TEST_NAME "AoD_test")
set(TEST_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}")
set(LIO_ARGS "-i aod.in -c HSOH.xyz -v ")
add_lio_test(${TEST_NAME} ${TEST_WORKING_DIR} ${LIO_ARGS} ${TEST_ENV})

