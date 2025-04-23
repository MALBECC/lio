cmake_minimum_required(VERSION 3.10)

# AVOID BUILD IN SOURCE DIRECTORY
if ( ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR} )
   message(FATAL_ERROR "You shouldn't build LIO in the source directory, please create a separate build directory and run cmake from there.")
endif()

# MAIN PROJECT - LIO
project(LIO VERSION 1.0 LANGUAGES FORTRAN CXX CUDA)

# CHECK OPTIONS COMPILATIONS
include(cmakefiles/compileOptions.cmake)

# ADD LIBRARIES
add_subdirectory(g2g)
add_subdirectory(lioamber)

# ENABLE FORTRAN
enable_language(Fortran)
add_compile_options(-O3)

# WARNING VARS
if (WARNING)
   add_compile_options(-Wall)
endif()

# LIOSOLO
add_executable(liosolo liosolo/liosolo.f90)
target_link_libraries(liosolo PUBLIC g2g lioamber)
target_include_directories(liosolo PUBLIC ${PROJECT_BINARY_DIR}/lioamber)

# LIOMD
add_executable(liomd.x liosolo/liomd.f90)
target_link_libraries(liomd.x PUBLIC g2g lioamber)
target_include_directories(liomd.x PUBLIC ${PROJECT_BINARY_DIR}/lioamber)

# INSTALL FOLDERS
install(TARGETS liosolo.x liomd.x DESTINATION bin)
install(TARGETS g2g lioamber LIBRARY DESTINATION lib)

