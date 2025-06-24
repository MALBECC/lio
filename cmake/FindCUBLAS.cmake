# FindCUBLAS.cmake - Find NVIDIA CUBLAS library
#
# This module finds the CUBLAS library that is part of the NVIDIA CUDA toolkit
#
# This module sets the following variables:
#  CUBLAS_FOUND - True if CUBLAS is found
#  CUBLAS_LIBRARIES - The CUBLAS libraries
#  CUBLAS_INCLUDE_DIRS - The CUBLAS include directories

# Find CUDA
find_package(CUDA REQUIRED)

# Find CUBLAS library
find_library(CUBLAS_LIBRARY
  NAMES cublas
  PATHS
    ${CUDA_TOOLKIT_ROOT_DIR}/lib64
    ${CUDA_TOOLKIT_ROOT_DIR}/lib
    /usr/local/cuda/lib64
    /usr/local/cuda/lib
)

# Find CUBLASLT library (part of newer CUDA releases)
find_library(CUBLASLT_LIBRARY
  NAMES cublasLt
  PATHS
    ${CUDA_TOOLKIT_ROOT_DIR}/lib64
    ${CUDA_TOOLKIT_ROOT_DIR}/lib
    /usr/local/cuda/lib64
    /usr/local/cuda/lib
)

# Set libraries
set(CUBLAS_LIBRARIES ${CUBLAS_LIBRARY})
if(CUBLASLT_LIBRARY)
  list(APPEND CUBLAS_LIBRARIES ${CUBLASLT_LIBRARY})
endif()

# Set include directories
set(CUBLAS_INCLUDE_DIRS ${CUDA_INCLUDE_DIRS})

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CUBLAS
  REQUIRED_VARS CUBLAS_LIBRARY CUBLAS_INCLUDE_DIRS
)

mark_as_advanced(
  CUBLAS_LIBRARY
  CUBLASLT_LIBRARY
)