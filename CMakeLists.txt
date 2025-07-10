cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# AVOID BUILD IN SOURCE DIRECTORY
if ( ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR} )
   message(FATAL_ERROR "You shouldn't build LIO in the source directory, please create a separate build directory and run cmake from there.")
endif()

# Set default installation directory to be inside the build directory
# TODO: Change this hard coded definition of CMAKE_INSTALL_PREFIX. 
message(STATUS "Checking for default install path")
if(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT CMAKE_INSTALL_PREFIX)
  message(STATUS "Setting default install path to ${CMAKE_BINARY_DIR}/install")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/" CACHE PATH "Default install path" )
endif()
# In your top-level CMakeLists.txt, before defining targets

# 1. Set the RPATH that will be used for INSTALLED targets.
#    '$ORIGIN' is a special linker variable meaning "the directory of the executable".
#    This assumes you install libraries to <prefix>/lib and executables to <prefix>/bin.
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

# 2. Instruct CMake to use the install RPATH settings for executables in the
#    BUILD tree. This is the magic that makes tests "just work".
#    CMake will automatically calculate the relative paths from your test executable
#    to the libraries in the build tree.
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)

# MAIN PROJECT - LIO
project(Lio 
   VERSION 1.0.0
   DESCRIPTION "LIO is a Quantum Chemistry software based on Density Functional Theory and real time Time Dependent Density Functional Theory."
)

# Global Settings 
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include custom CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# CHECK OPTIONS 
include(ConfigureOptions)

# ADD SUBDIRECTORIES.
add_subdirectory(g2g)
add_subdirectory(lioamber)
add_subdirectory(liosolo)
add_subdirectory(tools)

# Testing
if (BUILD_TESTING)
  message(STATUS "Configuring tests")
  enable_testing()
  add_subdirectory(test)
endif()

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dat" # Source directory (note the trailing slash)
  DESTINATION "${CMAKE_INSTALL_PREFIX}/"    # Destination in the install tree
)

install(
  FILES "${CMAKE_CURRENT_SOURCE_DIR}/lio_setenv_cmake.sh" # Source file
  DESTINATION "${CMAKE_INSTALL_PREFIX}/"    # Destination in the install tree
)

# Alternative installation of the libraries due to compatibility issues with Amber.
# This should be removed once Amber installation is updated to search for the libraries
# in the standard directories ${LIOHOME}/lib instead of ${LIOHOME}/g2g and ${LIOHOME}/lioamber.
install(CODE "
  file(MAKE_DIRECTORY \"\$<INSTALL_PREFIX>/g2g\")
  file(CREATE_LINK 
    \"\$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/\$<TARGET_FILE_NAME:g2g>\"
    \"\$<INSTALL_PREFIX>/g2g/\$<TARGET_LINKER_FILE_NAME:g2g>\"
  SYMBOLIC) 
  file(CREATE_LINK 
    \"\$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/\$<TARGET_FILE_NAME:g2g>\"
    \"\$<INSTALL_PREFIX>/g2g/\$<TARGET_SONAME_FILE_NAME:g2g>\"
  SYMBOLIC) 
" COMPONENT libraries) 

install(CODE "
  file(MAKE_DIRECTORY \"\$<INSTALL_PREFIX>/lioamber\")
  file(CREATE_LINK 
    \"\$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/\$<TARGET_FILE_NAME:lio-g2g>\"
    \"\$<INSTALL_PREFIX>/lioamber/\$<TARGET_LINKER_FILE_NAME:lio-g2g>\"
  SYMBOLIC) 
  file(CREATE_LINK 
    \"\$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/\$<TARGET_FILE_NAME:lio-g2g>\"
    \"\$<INSTALL_PREFIX>/lioamber/\$<TARGET_SONAME_FILE_NAME:lio-g2g>\"
  SYMBOLIC) 
" COMPONENT libraries) 

# Print build configuration summary
include(PrintBuildInfo)
message(STATUS "Build ORIGIN variable:${ORIGIN}")


