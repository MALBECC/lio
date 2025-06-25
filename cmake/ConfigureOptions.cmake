# Default Build Type.
#if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE Release CACHE STRING "Default build type (Release, Debug, etc.)" FORCE) # FORCE is okay for a default if user specifies nothing
#endif()

# Options for CUDA support.
option(USE_CUDA "Enable CUDA support" ON)            # Old cuda=1 (default)
option(USE_CUBLAS "Enable CUBLAS library" OFF)       # Old cuda=2 (implies USE_CUDA=ON)
option(USE_MAGMA "Enable MAGMA library" OFF)         # Old cuda=3 (implies USE_CUDA=ON and USE_CUBLAS=ON)

# Options for intel compilers and MKL
option(USE_INTEL_COMPILER "Use Intel compilers instead of GNU" OFF)
option(USE_MKL "Use Intel MKL library" OFF)

# Options for parallel processing (OpenMP)
option(USE_PARALLEL "Enable parallel processing" OFF)

# Option to use full double precision
option(USE_FULL_DOUBLE "Use double precision for calculations" OFF)
option(USE_CPU_DOUBLE "Use double precision for CPU calculations" ON)
option(AINT_MP "Enable multi-precision in analytic integrals" OFF)

# Options for external libraries
option(USE_LIBXC_CPU "Enable LIBXC CPU library" OFF)
option(USE_LIBXC_GPU "Use GPU version of LIBXC" OFF)
option(USE_LIBINT "Enable LIBINT library" OFF)

# Options for DFTD3 support.
option(USE_DFTD3 "Enable DFTD3 support" ON)

# Options for profiling and debugging.
set(USE_ANALYTICS 0 CACHE STRING "Debug Level for analytics" FORCE)

# CUDA-specific options
set(USE_CUDA_ARCH "35;52;61;75" CACHE STRING "CUDA architectures to compile for")
option(CUDA_USE_FAST_MATH "Use fast math in CUDA" ON)
option(CPU_RECOMPUTE "Recompute CPU energy" OFF)
option(CUDA_VERBOSE "Verbose CUDA compilation" OFF)
option(CUDA_PTX "Verbose PTX compilation" OFF)
option(CUDA_REGCOUNT "Verbose CUDA register compilation" OFF)
option(BUILD_TESTING "Build testing" OFF)
option(TEST_CHECKS_ONLY "Only perform checks for testing" OFF)

# Path configuration (can be overridden in cmake command line)
# TODO: Ckeck existence of environment variables LIBXC_HOME_GPU, LIBXC_HOME_CPU, etc.
#       Use these variables instead of the following.
#       Implement a module to find installation paths for lixc, libint, etc.
set(LIBXC_HOME_CPU "" CACHE PATH "Path to LIBXC CPU installation")
set(LIBXC_HOME_GPU "" CACHE PATH "Path to LIBXC GPU installation")
set(LIBINT_HOME "" CACHE PATH "Path to LIBINT installation")
set(EIGEN_HOME "" CACHE PATH "Path to Eigen installation")
set(MAGMA_ROOT "" CACHE PATH "Path to MAGMA installation")

