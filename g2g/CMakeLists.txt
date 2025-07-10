# CMakeLists.txt for g2g library
include(GNUInstallDirs)
# Include custom CMake modules
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Create the shared library
add_library(g2g SHARED)
set_target_properties(g2g PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    VERSION 1.0.0
    SOVERSION 1
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Set base include directories
target_include_directories(g2g 
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/datatypes
)

# Base compile definitions
target_compile_definitions(g2g PRIVATE 
    CPU_KERNELS=1
    UNIX  # From Makefile.cuda
)


# Handle OpenMP parallel options
if(USE_PARALLEL)
    find_package(OpenMP QUIET)
    if(OpenMP_FOUND)
        message(STATUS "OpenMP found, enabling for g2g library.")
        target_compile_definitions(g2g PUBLIC G2G_USE_OPENMP) # Bugfix: Allows g2g to be compiled without OpenMP (init.cpp, regenerate_partition.cpp).
        target_link_libraries(g2g PUBLIC OpenMP::OpenMP_CXX)
    else()
        message(WARNING "OpenMP requested for g2g but not found. Building without OpenMP.")
    endif()
else()
    message(STATUS "g2g: USE_PARALLEL is OFF. Not enabling OpenMP for g2g library.")
endif()

# Handle linear algebra library dependencies
if(USE_MKL AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    find_package(MKL REQUIRED)
    message(STATUS "MKL found, using MKL for g2g library.")
    target_link_libraries(g2g PRIVATE ${MKL_LIBRARIES})
    if(MKL_INCLUDE_DIRS)
        target_include_directories(g2g PRIVATE ${MKL_INCLUDE_DIRS})
    elseif(MKL_INCLUDE_DIR)
        target_include_directories(g2g PRIVATE ${MKL_INCLUDE_DIR})
    endif()
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    target_link_libraries(g2g PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    if(BLAS_INCLUDE_DIRS)
        target_include_directories(g2g PRIVATE ${BLAS_INCLUDE_DIRS})
    endif()
    if(LAPACK_INCLUDE_DIRS)
        target_include_directories(g2g PRIVATE ${LAPACK_INCLUDE_DIRS})
    endif()
endif()

# Compiler-specific options
target_compile_options(g2g PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wshadow>
    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:
        $<$<CONFIG:Debug>:-g -ggdb -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare>
        $<$<CONFIG:Release>:-O3 -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare>
    >
    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<OR:$<CXX_COMPILER_ID:IntelLLVM>,$<CXX_COMPILER_ID:Intel>>>:
        $<$<CONFIG:Debug>:-g -diag-enable=all -qopt-report=3>
        $<$<CONFIG:Release>:-O3 -qopt-report=3 -xHost -no-prec-div -ipo>
    >
)

# Debug/Release specific definitions
target_compile_definitions(g2g PRIVATE
    $<$<CONFIG:Debug>:_DEBUG>
)

# Intel compiler flag
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_definitions(g2g PRIVATE INTEL_COMP=1)
else()
    target_compile_definitions(g2g PRIVATE INTEL_COMP=0)
endif()

# CUDA Configuration
if(USE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)

    if(${CUDAToolkit_VERSION_MAJOR} GREATER_EQUAL 12)
        message(WARNING "CUDA ${CUDAToolkit_VERSION} may not be fully supported. Max tested is 11.x.")
    endif()
    
    # IMPORTANT: Enable CUDA separable compilation for proper device linking
    set_target_properties(g2g PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )    

    # Use modern CMake CUDA imported targets
    target_link_libraries(g2g PUBLIC CUDA::cudart CUDA::cuda_driver)
    target_compile_definitions(g2g PRIVATE GPU_KERNELS=1 _FORCE_INLINES)
    
    # Set CUDA architectures
    # Defaults are set in cmake/ConfigureOptions.cmake.
    # It would be better to set defaults according to CUDA version.
    message(STATUS "Setting CUDA_ARCHITECTURES: ${USE_CUDA_ARCH}")
    set_property(TARGET g2g PROPERTY CUDA_ARCHITECTURES "${USE_CUDA_ARCH}")
    
    # CUDA compile options matching Makefile.cuda
    target_compile_options(g2g PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:
            --compiler-options=-fPIC; # Semicolon
            --compiler-options=-Wall; # Semicolon
            --linker-options=-fPIC;   # Semicolon
            -Xcudafe=--diag_suppress=subscript_out_of_range
        >
        $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Debug>>:
            -g;                   # Semicolon
            -O1;                  # Semicolon
            --debug;              # Semicolon
            --generate-line-info; # Semicolon
#            --device-debug        # Semicolon (optional at the end of a group)
        >
        $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Release>>:
            -O3;
            --compiler-options=-fno-strict-aliasing
        >
        $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<BOOL:${CUDA_USE_FAST_MATH}>>:
            -use_fast_math
        >
        $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<BOOL:${CUDA_VERBOSE}>>:
            --verbose
        >
        $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<BOOL:${CUDA_PTX}>>:
            --ptx
        >
        $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<BOOL:${CUDA_REGCOUNT}>>:
            --ptxas-options=-v
        >
    )
    
    # CUBLAS/MAGMA configuration
    if(USE_CUBLAS AND USE_MAGMA)
        message(FATAL_ERROR "Cannot use both CUBLAS and MAGMA. Please choose one.")
    endif()
    
    if(USE_CUBLAS OR USE_MAGMA)
        target_compile_definitions(g2g PRIVATE CUBLAS)  # This was missing!
        target_link_libraries(g2g PUBLIC CUDA::cublas)
        
        if(USE_MAGMA)
            find_package(MAGMA REQUIRED)
            target_link_libraries(g2g PUBLIC ${MAGMA_LIBRARIES})
            target_compile_definitions(g2g PRIVATE magma)
            if(MAGMA_INCLUDE_DIRS)
                target_include_directories(g2g PRIVATE ${MAGMA_INCLUDE_DIRS})
            endif()
        endif()
    endif()
endif()

# LIBXC Configuration 
if(USE_LIBXC_CPU OR USE_LIBXC_GPU)
    if(USE_LIBXC_CPU AND USE_LIBXC_GPU)
        message(FATAL_ERROR "Cannot use both CPU and GPU versions of LIBXC.")
    endif()
    if(USE_LIBXC_CPU AND USE_CUDA)
        message(FATAL_ERROR "USE_LIBXC_CPU AND USE_CUDA are mutually exclusive.")
    endif()
    
    if(USE_LIBXC_CPU)
        # CPU version - matching Makefile.cuda libxc=1
        target_compile_definitions(g2g PRIVATE USE_LIBXC=1 LIBXC_CPU=1)
        find_package(LibXC REQUIRED)
        target_link_libraries(g2g PRIVATE ${LibXC_LIBRARIES})
        if(LibXC_INCLUDE_DIRS)
            target_include_directories(g2g PRIVATE ${LibXC_INCLUDE_DIRS})
        endif()
    else()  # USE_LIBXC_GPU
        # GPU version - matching Makefile.cuda libxc=2
        target_compile_definitions(g2g PRIVATE USE_LIBXC=1 LIBXC_CPU=0)
        find_package(LibXC-gpu REQUIRED)
        target_link_libraries(g2g PRIVATE ${LibXC-gpu_LIBRARIES})
        if(LibXC-gpu_INCLUDE_DIRS)
            target_include_directories(g2g PRIVATE ${LibXC-gpu_INCLUDE_DIRS})
        endif()
        # Also need CPU LibXC
        find_package(LibXC REQUIRED)
        target_link_libraries(g2g PRIVATE ${LibXC_LIBRARIES})
        if(LibXC_INCLUDE_DIRS)
            target_include_directories(g2g PRIVATE ${LibXC_INCLUDE_DIRS})
        endif()
    endif()
else()
    # No LIBXC - matching Makefile.cuda default
    target_compile_definitions(g2g PRIVATE USE_LIBXC=0 LIBXC_CPU=1)
endif()

# LIBINT configuration
if(USE_LIBINT)
    find_package(Libint2 2.6.0 EXACT REQUIRED)
    find_package(Eigen3 REQUIRED)

    # --- Start of Manual Property Extraction and Application ---
    # We do not link to Libint2::cxx or Eigen3::Eigen directly.
    # Instead, we extract their properties and apply them manually
    # with language-specific guards to prevent them from leaking into CUDA compilation.
    # The reason is that libint include files should only be linked to CXX sources,
    # not CUDA sources. If not done this way, the nvcc compiler will choke on 
    # Libint's headers.
    #

    # 1. Extract properties from Libint2::cxx
    get_target_property(LIBINT2_LIBRARIES Libint2::cxx INTERFACE_LINK_LIBRARIES)
    get_target_property(LIBINT2_INCLUDE_DIRS Libint2::cxx INTERFACE_INCLUDE_DIRECTORIES)

    # 2. Extract properties from Eigen3::Eigen
    get_target_property(EIGEN3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(EIGEN3_COMPILE_DEFS Eigen3::Eigen INTERFACE_COMPILE_DEFINITIONS)

    # 3. Apply the LIBRARIES to the linker (this is language-agnostic and safe)
    target_link_libraries(g2g PRIVATE ${LIBINT2_LIBRARIES})

    # 4. Apply INCLUDE_DIRECTORIES only to CXX sources
    target_include_directories(g2g PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:${LIBINT2_INCLUDE_DIRS}>
        $<$<COMPILE_LANGUAGE:CXX>:${EIGEN3_INCLUDE_DIRS}>
    )
    # --- End of Manual Property Extraction and Application ---

    # Definition -DUSE_LIBINT=1 applies ONLY to CXX.
    target_compile_definitions(g2g PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:USE_LIBINT=1>
    )
    target_compile_features(g2g PRIVATE cxx_std_11)
else()
    target_compile_definitions(g2g PRIVATE USE_LIBINT=0)
endif()
#--------------------------------------------------------

# Full double precision
if(USE_FULL_DOUBLE)
    target_compile_definitions(g2g PRIVATE FULL_DOUBLE=1)
else()
    target_compile_definitions(g2g PRIVATE FULL_DOUBLE=0)
endif()

# CPU recompute
if(CPU_RECOMPUTE)
    target_compile_definitions(g2g PRIVATE CPU_RECOMPUTE=1)
else()
    target_compile_definitions(g2g PRIVATE CPU_RECOMPUTE=0)
endif()

# AINT_MP
if(AINT_MP)
    target_compile_definitions(g2g PRIVATE AINT_MP=1)
else()
    target_compile_definitions(g2g PRIVATE AINT_MP=0)
endif()

# Source files - Base sources
set(G2G_BASE_SOURCES
    classify_functions.cpp
    global_memory_pool.cpp
    init.cpp
    matrix.cpp
    partition.cpp
    regenerate_partition.cpp
    timer.cpp
    cdft.cpp
    extern_functional.cpp
)

# CPU sources
set(G2G_CPU_SOURCES
    cpu/functions.cpp
    cpu/iteration.cpp
    cpu/weight.cpp
    cpu/compute_Wmat.cpp
)

# Analytic integral sources
set(G2G_AINT_SOURCES
    analytic_integral/aint_init.cpp
    analytic_integral/os_cutoff.cpp
    analytic_integral/qmmm_nuc.cpp
)

# PointXC sources
file(GLOB G2G_POINTXC_SOURCES pointxc/*.cpp)

# Excited state sources
set(G2G_EXCITED_SOURCES excited/saverho.cpp)
if(USE_LIBXC_CPU OR USE_LIBXC_GPU)
    list(APPEND G2G_EXCITED_SOURCES
        excited/g2g_calculateXC.cpp
        excited/g2g_calculateG.cpp
        excited/g2g_calcgradXC.cpp
        excited/calc_gradients.cpp
        excited/calc_FXC.cpp
        excited/calc_VXC.cpp
    )
else()
    list(APPEND G2G_EXCITED_SOURCES excited/fake_subs.cpp)
endif()

# LIBINT sources
if(USE_LIBINT)
    set(G2G_LIBINT_SOURCES
        libint/libintproxy.cpp
        libint/g2g_libint.cpp
    )
else()
    set(G2G_LIBINT_SOURCES libint/g2g_libint_fake.cpp)
endif()

# CIOVERLAP sources
if(BUILD_SURFACE_HOPPING)
    set(G2G_CIOVERLAP_SOURCES
        cioverlap/cioverlap.cpp
        cioverlap/g2g_cioverlap.cpp
    )
    find_package(GSL REQUIRED)
    target_link_libraries(g2g PRIVATE GSL::gsl)
    if(LA_HOME)
        target_include_directories(g2g PRIVATE ${LA_HOME}/include)
        target_link_directories(g2g PRIVATE ${LA_HOME}/lib)
        target_link_libraries(g2g PRIVATE la)
    endif()
else()
    set(G2G_CIOVERLAP_SOURCES cioverlap/cioverlap_fake.cpp)
endif()

# Handle source compilation based on CUDA availability
if(USE_CUDA)
    # When CUDA is enabled, pointxc files are compiled with nvcc (from Makefile.cuda)
    set_source_files_properties(${G2G_POINTXC_SOURCES} PROPERTIES LANGUAGE CUDA)
    
    # Add all non-CUDA C++ sources
    target_sources(g2g PRIVATE
        ${G2G_BASE_SOURCES}
        ${G2G_CPU_SOURCES}
        ${G2G_AINT_SOURCES}
        ${G2G_POINTXC_SOURCES}  # Now compiled as CUDA
        ${G2G_EXCITED_SOURCES}
        ${G2G_LIBINT_SOURCES}
        ${G2G_CIOVERLAP_SOURCES}
    )
    
    # CUDA sources
    set(G2G_CUDA_SOURCES
        cuda/iteration.cu
        cuda/compute_Wmat.cu
    )
    
    # Analytic integral CUDA sources
    list(APPEND G2G_CUDA_SOURCES
        analytic_integral/cuda/coulomb.cu
        analytic_integral/cuda/os_common.cu
        analytic_integral/cuda/qmmm.cu
    )
    
    # Excited state CUDA sources
    if(USE_LIBXC_GPU)
        list(APPEND G2G_CUDA_SOURCES
            excited/cuda/g2g_calculateXC.cu
            excited/cuda/g2g_calculateG.cu
            excited/cuda/g2g_calcgradXC.cu
        )
    else()
        list(APPEND G2G_CUDA_SOURCES excited/cuda/fake_cuda.cu)
    endif()
    
    # Add CUDA sources directly to target
    # With CUDA_SEPARABLE_COMPILATION and CUDA_RESOLVE_DEVICE_SYMBOLS set,
    # CMake will handle the device linking automatically
    target_sources(g2g PRIVATE ${G2G_CUDA_SOURCES})
else()
    # Without CUDA, compile everything as regular C++
    target_sources(g2g PRIVATE
        ${G2G_BASE_SOURCES}
        ${G2G_CPU_SOURCES}
        ${G2G_AINT_SOURCES}
        ${G2G_POINTXC_SOURCES}
        ${G2G_EXCITED_SOURCES}
        ${G2G_LIBINT_SOURCES}
        ${G2G_CIOVERLAP_SOURCES}
    )
endif()

# Additional system libraries (from original Makefile)
target_link_libraries(g2g PRIVATE rt m )

# Installation and Export Configuration
# -------------------------------------

# Export targets (only define once!)
install(TARGETS g2g
    EXPORT g2gTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install the export file
install(EXPORT g2gTargets
    FILE g2gTargets.cmake
    NAMESPACE g2g::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/g2g
)

# Generate and install config files only if the template exists
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/g2gConfig.cmake.in")
    # Generate version file
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/g2gConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    # Configure the package config file
    configure_package_config_file(
        "${CMAKE_SOURCE_DIR}/cmake/g2gConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/g2gConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/g2g
        PATH_VARS CMAKE_INSTALL_LIBDIR CMAKE_INSTALL_INCLUDEDIR
    )

    # Install config files
    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/g2gConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/g2gConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/g2g
    )

    # Optionally install custom Find modules for dependencies
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
        set(G2G_CMAKE_MODULE_PATH ${CMAKE_INSTALL_LIBDIR}/cmake/g2g/modules)
        install(DIRECTORY cmake/modules/
            DESTINATION ${G2G_CMAKE_MODULE_PATH}
            FILES_MATCHING PATTERN "Find*.cmake"
        )
    endif()
else()
    message(WARNING "g2gConfig.cmake.in not found at ${CMAKE_SOURCE_DIR}/cmake/g2gConfig.cmake.in")
    message(WARNING "Package configuration files will not be generated.")
endif()