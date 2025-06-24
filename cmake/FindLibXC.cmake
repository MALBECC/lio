# FindLibXC.cmake
# ---------------
#
# Find LibXC library
#
# This module finds the LibXC library and defines:
#
#   LibXC_FOUND        - True if LibXC was found
#   LibXC_INCLUDE_DIRS - The LibXC include directories
#   LibXC_LIBRARIES    - The LibXC library
#   LibXC_VERSION      - Version string of LibXC
#
# The following imported targets are also defined:
#
#   LibXC::xc          - The LibXC library
#
# The following variables can be set to guide the search:
#
#   LibXC_ROOT         - Root directory of LibXC installation
#   LIBXC_ROOT         - Alternative root directory variable
#   LIBXC_HOME_CPU     - LIO-specific root directory variable

# Handle different possible root directory variables
if(NOT LibXC_ROOT)
    if(DEFINED ENV{LIBXC_ROOT})
        set(LibXC_ROOT $ENV{LIBXC_ROOT})
    elseif(DEFINED ENV{LIBXC_HOME_CPU})
        set(LibXC_ROOT $ENV{LIBXC_HOME_CPU})
    elseif(LIBXC_ROOT)
        set(LibXC_ROOT ${LIBXC_ROOT})
    elseif(LIBXC_HOME_CPU)
        set(LibXC_ROOT ${LIBXC_HOME_CPU})
    endif()
endif()

# Find the header file
find_path(LibXC_INCLUDE_DIR
    NAMES xc.h
    HINTS
        ${LibXC_ROOT}
        ${LIBXC_CPU_INCLUDES}  # From LIO variables
    PATH_SUFFIXES
        include
    DOC "LibXC include directory"
)

# Find the library
find_library(LibXC_LIBRARY
    NAMES xc
    HINTS
        ${LibXC_ROOT}
        ${LIBXC_CPU_LIBS}  # From LIO variables
    PATH_SUFFIXES
        lib
        lib64
    DOC "LibXC library"
)

# Extract version information if header is found
if(LibXC_INCLUDE_DIR AND EXISTS "${LibXC_INCLUDE_DIR}/xc_version.h")
    file(STRINGS "${LibXC_INCLUDE_DIR}/xc_version.h" LibXC_VERSION_LINE
        REGEX "^#define[\t ]+XC_VERSION[\t ]+\".*\"")
    if(LibXC_VERSION_LINE)
        string(REGEX REPLACE "^#define[\t ]+XC_VERSION[\t ]+\"(.*)\"" "\\1"
            LibXC_VERSION "${LibXC_VERSION_LINE}")
    endif()
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibXC
    REQUIRED_VARS
        LibXC_LIBRARY
        LibXC_INCLUDE_DIR
    VERSION_VAR
        LibXC_VERSION
)

# Set output variables
if(LibXC_FOUND)
    set(LibXC_LIBRARIES ${LibXC_LIBRARY})
    set(LibXC_INCLUDE_DIRS ${LibXC_INCLUDE_DIR})
    
    # Create imported target
    if(NOT TARGET LibXC::xc)
        add_library(LibXC::xc UNKNOWN IMPORTED)
        set_target_properties(LibXC::xc PROPERTIES
            IMPORTED_LOCATION "${LibXC_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${LibXC_INCLUDE_DIR}"
        )
    endif()
endif()

# Hide internal variables
mark_as_advanced(LibXC_INCLUDE_DIR LibXC_LIBRARY)