# FindMKL.cmake - Find Intel MKL libraries
#
# This module finds the Intel Math Kernel Library (MKL)
#
# This module sets the following variables:
#  MKL_FOUND - True if MKL is found
#  MKL_INCLUDE_DIR - The MKL include directory
#  MKL_LIBRARIES - The MKL libraries
#  MKL_ROOT - The MKL root directory

# Find MKL root
find_path(MKL_ROOT
  NAMES include/mkl.h include/mkl_cblas.h
  PATHS
    $ENV{MKLROOT}
    /opt/intel/mkl
    /opt/intel/oneapi/mkl/latest
    /opt/intel/compilers_and_libraries/linux/mkl
)

if(MKL_ROOT)
  # Find include directory
  find_path(MKL_INCLUDE_DIR
    NAMES mkl.h
    PATHS ${MKL_ROOT}/include
  )
  
  # Set library paths
  set(MKL_LIBRARY_DIR ${MKL_ROOT}/lib/intel64)
  
  # Set libraries based on architecture and threading
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(MKL_INTERFACE_LIBRARY "mkl_intel_lp64")
  else()
    set(MKL_INTERFACE_LIBRARY "mkl_intel")
  endif()
  
  # Set threading library
  if(OpenMP_FOUND)
    set(MKL_THREADING_LIBRARY "mkl_intel_thread")
  else()
    set(MKL_THREADING_LIBRARY "mkl_sequential")
  endif()
  
  # Set core library
  set(MKL_CORE_LIBRARY "mkl_core")
  
  # Add Fortran interfaces if Fortran is enabled
  if(CMAKE_Fortran_COMPILER_LOADED)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(MKL_LAPACK_LIBRARY "mkl_lapack95_lp64")
      set(MKL_BLAS_LIBRARY "mkl_blas95_lp64")
    else()
      set(MKL_LAPACK_LIBRARY "mkl_lapack95")
      set(MKL_BLAS_LIBRARY "mkl_blas95")
    endif()
    list(APPEND MKL_LIBRARIES ${MKL_BLAS_LIBRARY} ${MKL_LAPACK_LIBRARY})
  endif()
  
  # Set libraries
  list(APPEND MKL_LIBRARIES
    ${MKL_INTERFACE_LIBRARY}
    ${MKL_THREADING_LIBRARY}
    ${MKL_CORE_LIBRARY}
    iomp5 pthread m dl
  )
  
  # Convert to full paths
  set(MKL_LIBRARIES_FULL)
  foreach(LIB ${MKL_LIBRARIES})
    find_library(MKL_${LIB}_LIBRARY
      NAMES ${LIB}
      PATHS ${MKL_LIBRARY_DIR}
    )
    if(MKL_${LIB}_LIBRARY)
      list(APPEND MKL_LIBRARIES_FULL ${MKL_${LIB}_LIBRARY})
    endif()
    unset(MKL_${LIB}_LIBRARY CACHE)
  endforeach()
  
  set(MKL_LIBRARIES ${MKL_LIBRARIES_FULL})
endif()

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL
  REQUIRED_VARS MKL_ROOT MKL_INCLUDE_DIR MKL_LIBRARIES
)

mark_as_advanced(
  MKL_ROOT
  MKL_INCLUDE_DIR
  MKL_LIBRARIES
)