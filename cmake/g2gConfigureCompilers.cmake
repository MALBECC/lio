# Configure compilers flags based on options.

# Enable Fortran language support.

# Set default compiler flags.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Handle optional use of Intel compilers
if(USE_INTEL_COMPILER)
    set(CMAKE_CXX_COMPILER "icpc" CACHE STRING "C++ compiler." FORCE)
    set(CMAKE_C_COMPILER "icc" CACHE STRING "C compiler." FORCE)
    message(STATUS "Attempting to use Intel compilers icpc/icc.")
  
    # Handle optional use of MKL
    if(USE_MKL)
        find_package(MKL REQUIRED)
        include_directories(${MKL_INCLUDE_DIR})
        list(APPEND EXTRA_LIBS ${MKL_LIBRARIES})
    else()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        list(APPEND EXTRA_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()
    
    # Configure Intel Fortran compiler flags.
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp ")

    # Optimization Flags.
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -ip -mp1")
    if (USE_PARALLEL)
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -parallel")
    endif()

    # Configure Intel C++ compiler flags.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -xHost -no-prec-div") 
    # Check for Intel compiler version >= 16.0 to set appropriate optimization report flag
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 16.0)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report 3")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -opt-report 3")
    endif()
else()
    # Otherwise use GNU compilers and standard BLAS/LAPACK
    message(STATUS "Using GNU Fortran compiler")
    set(CMAKE_Fortran_COMPILER "gfortran")
    set(CMAKE_CXX_COMPILER "g++")
    set(CMAKE_C_COMPILER "gcc")

    # Handle linear algebra libraries.
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    list(APPEND EXTRA_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    
    # Configure Fortran compiler flags.
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp ")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "10.0")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch ")
    endif()

    # Configure CXX compiler flags.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare")
endif()


# Debug and profiling configurations
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    if(${USE_ANALYTICS} STREQUAL 1) # Basic debug flags.
        add_compile_options(-pg) # !!! Is this redundant?
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -pg")
        set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -pg")
    elseif(USE_ANALYTICS STREQUAL 2) # Intermediate debug flags.
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -Wtabs -fcheck=all")
    elseif(USE_ANALYTICS STREQUAL 3) # Advanced debug flags.
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -Wextra -fcheck=all -fbacktrace -pedantic -ffpe-trap=invalid,zero,overflow,underflow")
    elseif(USE_ANALYTICS STREQUAL 4) # Advanced Advanced debug flags.
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall -Wextra -fcheck=all -fbacktrace -pedantic -ffpe-trap=invalid,zero,overflow,underflow")
    else()
        message(FATAL_ERROR "Unknown analytics level ${USE_ANALYTICS}. Alowed values are 1, 2, 3 and 4.")
    endif()
endif()


# CUDA Configuration
if(USE_CUDA)
    # Find CUDA package.
    find_package(CUDA REQUIRED)

    # Version-dependent architecture support
    if(${CUDA_VERSION_MAJOR} LESS 8)
        # For CUDA versions < 8, prioritize older architectures
        if(CUDA_SM20_SUPPORT)
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_20,code=compute_20")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_20,code=sm_20")
        endif()
        if(CUDA_SM30_SUPPORT)
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_30,code=compute_30")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_30,code=sm_30")
        endif()
        # Default architecture set for CUDA < 8
        if(NOT CUDA_SM20_SUPPORT AND NOT CUDA_SM30_SUPPORT AND NOT CUDA_SM35_SUPPORT)
            message(STATUS "No specific SM architecture selected for CUDA < 8, enabling SM30 by default")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_30,code=compute_30")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_30,code=sm_30")
        endif()

    # Add elseif (${CUDA_VERSION_MAJOR} GREATER_EQUAL 12) # CUDA 12 or greater not supported yet. Abort compilation.
    else()
        # For CUDA versions >= 8, use newer architectures
        # Handle user-selected architectures from CUDA_ARCH
        foreach(ARCH ${CUDA_ARCH})
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_${ARCH},code=compute_${ARCH}")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_${ARCH},code=sm_${ARCH}")
        endforeach()
        
        # Add specific architecture support based on options
        if(CUDA_SM30_SUPPORT)
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_30,code=compute_30")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_30,code=sm_30")
        endif()
        if(CUDA_SM52_SUPPORT)
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_52,code=compute_52")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_52,code=sm_52")
        endif()
        if(CUDA_SM61_SUPPORT)
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_61,code=compute_61")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_61,code=sm_61")
        endif()
        if(CUDA_SM75_SUPPORT)
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_75,code=compute_75")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_75,code=sm_75")
        endif()
        
        # Default architecture set for CUDA >= 8 (if no architectures selected)
        if(NOT CUDA_ARCH AND NOT CUDA_SM30_SUPPORT AND NOT CUDA_SM35_SUPPORT AND 
           NOT CUDA_SM50_SUPPORT AND NOT CUDA_SM52_SUPPORT AND
           NOT CUDA_SM60_SUPPORT AND NOT CUDA_SM61_SUPPORT AND NOT CUDA_SM75_SUPPORT)
            message(STATUS "No specific SM architecture selected for CUDA >= 8, enabling SM30, SM52, SM61 by default")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_30,code=compute_30")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_30,code=sm_30")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_52,code=compute_52")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_52,code=sm_52")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_61,code=compute_61")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_61,code=sm_61")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_75,code=compute_75")
            set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_75,code=sm_75")
        endif()
    endif()
    # Set CUDA flags.
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --compiler-options -fPIC -Wall")

    # Optimization and debugging flags.
    if (CMAKE_BUILD_TYPE NOT MATCHES Debug)
        # Standard CUDA optimization flags
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --compiler-options -fno-strict-aliasing")
    else()
        # Debugging flags
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -O1 --debug --device-debug")
    endif() 

    if (CUDA_USE_FAST_MATH)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -use_fast_math")
    endif()

    # Use cublas if requested.
    if (USE_CUBLAS)
        find_package(CUBLAS REQUIRED)
        list(APPEND EXTRA_LIBS ${CUBLAS_LIBRARIES})
        if (USE_MAGMA)
            find_package(MAGMA REQUIRED)
            include_directories(${MAGMA_INCLUDE_DIR})
            list(APPEND EXTRA_LIBS ${MAGMA_LIBRARIES})
        endif()
    endif()
endif()


# Check mutually exclusive LIBXC options.
if (USE_LIBXC_GPU AND USE_LIBXC_CPU)
    message(FATAL_ERROR "Cannot use both LIBXC_CPU and LIBXC_GPU. Choose one.")
endif()

# GPU LIBXC (IN-HOUSE VERSION)
if (USE_LIBXC_GPU)
    if (NOT LIBXC_GPU_DIR) 
        message(FATAL_ERROR "LIBXC_GPU_DIR is not set")
    endif()
    include_directories(${LIBXC_GPU_DIR}/include ${LIBXC_CPU_DIR}/include)
    link_directories(${LIBXC_GPU_DIR}/lib ${LIBXC_CPU_DIR}/lib)
    list(APPEND EXTRA_LIBS xc_cuda xc)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DUSE_LIBXC=1 ")
endif()

# CPU LIBXC
if (USE_LIBXC_CPU)
    if (NOT LIBXC_CPU_DIR) 
        message(FATAL_ERROR "LIBXC_CPU_DIR is not set")
    endif()
    include_directories(${LIBXC_CPU_DIR}/include)
    link_directories(${LIBXC_CPU_DIR}/lib)
    list(APPEND EXTRA_LIBS xc)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DUSE_LIBXC=1 ")
endif()

# Configure LIBINT if enabled 
if (USE_LIBINT)
    if (NOT LIBINT_DIR)
        message(FATAL_ERROR "LIBINT_DIR is not set")
    endif()
    if (NOT EIGEN_DIR)
        message(FATAL_ERROR "EIGEN_DIR is not set")
    endif()
    include_directories(
        ${LIBINT_DIR}/include
        ${LIBINT_DIR}/include/libint2
        ${EIGEN_DIR}
    )
    link_directories(${LIBINT_DIR}/lib)
    list(APPEND EXTRA_LIBS int2)
    # Use C++11 for LIBINT support.
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# Add library paths to runtime loading.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "$ORIGIN;$ORIGIN/../g2g")

