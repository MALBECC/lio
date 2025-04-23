include(../cmakefiles/g2gUtils.cmake)

# Get the source files for G2G lib => SRC_FILES
source_files_g2g()

# Compilation of the library
add_compile_options(-std=c++17 -O3)
add_library(g2g SHARED ${SRC_FILES})

############################# LIBRARIES ################################
check_openmp()

# EXTERNAL
if (EXTERNAL)
   message(STATUS "Using libxc and libint libraries")
   if (COMPILE_EXTERNAL)
      message(STATUS "Download and install external libraries")
      # This function will download and install libxc and libint for cpu and then
      # will be linked in the compilation of lio
      download_external_libs()
   else()
      # This function will only look in your system and the environment variables
      # for libxc and libint to link with lio
      link_external_libs()
   endif()
else()
   message(STATUS "Don't use external libraries")
   target_compile_definitions(g2g PRIVATE "USE_LIBXC=0" "LIBXC_CPU=1" "USE_LIBINT=0")
endif()
############################# END LIBRARIES ################################

# CPU compilations options
if (CPU)
  target_compile_definitions(g2g PRIVATE "CPU_KERNELS=1" "CPU_RECOMPUTE=0" "INTEL_COMP=0")
endif()

# PRECISION FLAGS
if (DOUBLE)
   message(STATUS "Compiling LIO with double precision")
   target_compile_definitions(g2g PRIVATE "FULL_DOUBLE=1")
else()
   message(STATUS "Compiling LIO with simple precision")
   target_compile_definitions(g2g PRIVATE "FULL_DOUBLE=0")
endif()

# WARNING FLAGS
if (WARNING)
  add_compile_options(-Wall -Wextra -Wshadow -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare)
endif()

# HEADERS
target_include_directories(g2g PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} cpu analytic_integral cuda)
