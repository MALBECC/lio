# FindMAGMA.cmake - Find MAGMA libraries
#
# This module finds the MAGMA (Matrix Algebra on GPU and Multicore Architectures) library
#
# This module sets the following variables:
#  MAGMA_FOUND - True if MAGMA is found
#  MAGMA_INCLUDE_DIR - The MAGMA include directory
#  MAGMA_LIBRARIES - The MAGMA libraries
#  MAGMA_VERSION - The MAGMA version

# Find MAGMA root
if(NOT MAGMA_ROOT)
  set(MAGMA_ROOT $ENV{MAGMAROOT} CACHE PATH "MAGMA installation directory")
endif()

# Find include directory
find_path(MAGMA_INCLUDE_DIR
  NAMES magma.h
  PATHS
    ${MAGMA_ROOT}/include
    /usr/include/magma
    /usr/local/include/magma
)

# Find MAGMA library
find_library(MAGMA_LIBRARY
  NAMES magma
  PATHS
    ${MAGMA_ROOT}/lib
    /usr/lib
    /usr/local/lib
)

# Find MAGMA sparse library (optional)
find_library(MAGMA_SPARSE_LIBRARY
  NAMES magma_sparse
  PATHS
    ${MAGMA_ROOT}/lib
    /usr/lib
    /usr/local/lib
)

# Set libraries
set(MAGMA_LIBRARIES ${MAGMA_LIBRARY})
if(MAGMA_SPARSE_LIBRARY)
  list(APPEND MAGMA_LIBRARIES ${MAGMA_SPARSE_LIBRARY})
endif()

# Try to get version
if(MAGMA_INCLUDE_DIR)
  file(STRINGS "${MAGMA_INCLUDE_DIR}/magma_types.h" MAGMA_VERSION_MAJOR
    REGEX "^#define[ \t]+MAGMA_VERSION_MAJOR[ \t]+[0-9]+$")
  file(STRINGS "${MAGMA_INCLUDE_DIR}/magma_types.h" MAGMA_VERSION_MINOR
    REGEX "^#define[ \t]+MAGMA_VERSION_MINOR[ \t]+[0-9]+$")
  file(STRINGS "${MAGMA_INCLUDE_DIR}/magma_types.h" MAGMA_VERSION_MICRO
    REGEX "^#define[ \t]+MAGMA_VERSION_MICRO[ \t]+[0-9]+$")
  
  string(REGEX REPLACE "^#define[ \t]+MAGMA_VERSION_MAJOR[ \t]+([0-9]+)$" "\\1"
    MAGMA_VERSION_MAJOR "${MAGMA_VERSION_MAJOR}")
  string(REGEX REPLACE "^#define[ \t]+MAGMA_VERSION_MINOR[ \t]+([0-9]+)$" "\\1"
    MAGMA_VERSION_MINOR "${MAGMA_VERSION_MINOR}")
  string(REGEX REPLACE "^#define[ \t]+MAGMA_VERSION_MICRO[ \t]+([0-9]+)$" "\\1"
    MAGMA_VERSION_MICRO "${MAGMA_VERSION_MICRO}")
  
  set(MAGMA_VERSION "${MAGMA_VERSION_MAJOR}.${MAGMA_VERSION_MINOR}.${MAGMA_VERSION_MICRO}")
endif()

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MAGMA
  REQUIRED_VARS MAGMA_LIBRARY MAGMA_INCLUDE_DIR
  VERSION_VAR MAGMA_VERSION
)

mark_as_advanced(
  MAGMA_INCLUDE_DIR
  MAGMA_LIBRARY
  MAGMA_SPARSE_LIBRARY
)