# FindLibXC-gpu.cmake
# ---------------
#
# Find LibXC library
#
# This module finds the LibXC library and defines:
#
#   LibXC-gpu_FOUND        - True if LibXC was found
#   LibXC-gpu_INCLUDE_DIRS - The LibXC include directories
#   LibXC-gpu_LIBRARIES    - The LibXC library
#   LibXC-gpu_VERSION      - Version string of LibXC
#
# The following imported targets are also defined:
#
#   LibXC::xc          - The LibXC library
#
# The following variables can be set to guide the search:
#
#   LibXC-gpu_ROOT         - Root directory of LibXC installation
#   LIBXC_ROOT         - Alternative root directory variable
#   LIBXC_HOME_GPU     - LIO-specific root directory variable

# Handle different possible root directory variables
if(NOT LibXC-gpu_ROOT)
    if(DEFINED ENV{LIBXC_ROOT})
        set(LibXC-gpu_ROOT $ENV{LIBXC_ROOT})
    elseif(DEFINED ENV{LIBXC_HOME_GPU})
        set(LibXC-gpu_ROOT $ENV{LIBXC_HOME_GPU})
    elseif(LIBXC_ROOT)
        set(LibXC-gpu_ROOT ${LIBXC_ROOT})
    elseif(LIBXC_HOME_GPU)
        set(LibXC-gpu_ROOT ${LIBXC_HOME_GPU})
    endif()
endif()
message(STATUS "LibXC-gpu_ROOT is: ${LibXC-gpu_ROOT}")

# Find the header file
find_path(LibXC-gpu_INCLUDE_DIR
    NAMES xc_cuda.h
    HINTS
        ${LibXC-gpu_ROOT}
        ${LIBXC_GPU_INCLUDES}  # From LIO variables
    PATH_SUFFIXES
        include
    DOC "LibXC include directory"
)

# Find the library
find_library(LibXC-gpu_LIBRARY
    NAMES xc_cuda
    HINTS
        ${LibXC-gpu_ROOT}
        ${LIBXC_GPU_LIBS}  # From LIO variables
    PATH_SUFFIXES
        lib
        lib64
    DOC "LibXC library"
)

# Extract version information if header is found
if(LibXC-gpu_INCLUDE_DIR AND EXISTS "${LibXC-gpu_INCLUDE_DIR}/xc_version.h")
    file(STRINGS "${LibXC-gpu_INCLUDE_DIR}/xc_version.h" LibXC-gpu_VERSION_LINE
        REGEX "^#define[\t ]+XC_VERSION[\t ]+\".*\"")
    if(LibXC-gpu_VERSION_LINE)
        string(REGEX REPLACE "^#define[\t ]+XC_VERSION[\t ]+\"(.*)\"" "\\1"
            LibXC-gpu_VERSION "${LibXC-gpu_VERSION_LINE}")
    endif()
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibXC-gpu
    REQUIRED_VARS
        LibXC-gpu_LIBRARY
        LibXC-gpu_INCLUDE_DIR
    VERSION_VAR
        LibXC-gpu_VERSION
)

# Set output variables
if(LibXC-gpu_FOUND)
    set(LibXC-gpu_LIBRARIES ${LibXC-gpu_LIBRARY})
    set(LibXC-gpu_INCLUDE_DIRS ${LibXC-gpu_INCLUDE_DIR})
    
    # Create imported target
    if(NOT TARGET LibXC-gpu::xc-gpu)
        add_library(LibXC-gpu::xc-gpu UNKNOWN IMPORTED)
        set_target_properties(LibXC-gpu::xc-gpu PROPERTIES
            IMPORTED_LOCATION "${LibXC-gpu_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${LibXC-gpu_INCLUDE_DIR}"
        )
    endif()
endif()

# Hide internal variables
mark_as_advanced(LibXC-gpu_INCLUDE_DIR LibXC-gpu_LIBRARY)
