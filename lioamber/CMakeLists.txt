# CMakeLists.txt for lioamber library
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# Enable Fortran language
enable_language(Fortran)

# Include GNU install directories
include(GNUInstallDirs)

# Create the shared library
add_library(lio-g2g SHARED)

# Set library properties
set_target_properties(lio-g2g PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION 1.0.0
    SOVERSION 1
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/obj
)
message(STATUS "Current binary directory is: ${CMAKE_CURRENT_BINARY_DIR}")

# Include directories
target_include_directories(lio-g2g
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/obj>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/datatypes
)

# Find and link to g2g library (built previously)
#find_package(g2g REQUIRED PATHS ${CMAKE_BINARY_DIR}/g2g)
target_link_libraries(lio-g2g PUBLIC g2g)

# Basic compile definitions. Flaged in original Makefile as probably not necessary.
#target_compile_definitions(lio-g2g PRIVATE
#    pack
#    G2G
#)

# Handle TD precision option
if(NOT USE_FULL_DOUBLE)
    # Simple precision for TD (precision=0 in Makefile)
    target_compile_definitions(lio-g2g PRIVATE TD_SIMPLE)
endif()

# Handle CPU precision
if(USE_CPU_DOUBLE)
    target_compile_definitions(lio-g2g PRIVATE CPU_DOUBLE)
else()
    target_compile_definitions(lio-g2g PRIVATE CPU_SIMPLE)
endif()

# Apply analytics settings based on level
if(USE_ANALYTICS EQUAL 1)
    target_compile_definitions(lio-g2g PRIVATE 
        $<$<CONFIG:Debug>:PROFILING>
    )
    target_compile_options(lio-g2g PRIVATE 
        $<$<CONFIG:Debug>:-pg>
    )
    target_link_options(lio-g2g PRIVATE 
        $<$<CONFIG:Debug>:-pg>
    )
elseif(USE_ANALYTICS EQUAL 2)
    target_compile_definitions(lio-g2g PRIVATE 
        $<$<CONFIG:Debug>:DEBUGGING>
        $<$<CONFIG:Debug>:PRINT_MATRICES>
    )
elseif(USE_ANALYTICS EQUAL 3)
    target_compile_definitions(lio-g2g PRIVATE 
        $<$<CONFIG:Debug>:DEBUGGING>
        $<$<CONFIG:Debug>:PRINT_MATRICES>
    )
elseif(USE_ANALYTICS EQUAL 4)
    target_compile_definitions(lio-g2g PRIVATE 
        $<$<CONFIG:Debug>:DEBUGGING>
        $<$<CONFIG:Debug>:PRINT_MATRICES>
        $<$<CONFIG:Debug>:FULL_CHECKS>
    )
endif()

# Compiler-specific configurations
# Define generator expressions as variables to improve readability.
set(is_intel [[$<OR:$<Fortran_COMPILER_ID:IntelLLVM>,$<Fortran_COMPILER_ID:Intel>>]])
set(is_config_debug [[$<CONFIG:Debug>]])
set(is_config_release [[$<CONFIG:Release>]])
set(is_parallel [[$<BOOL:${USE_PARALLEL}>]])
# For GNU Fortran version check
set(is_gnu "$<Fortran_COMPILER_ID:GNU>")
set(is_gnu_ge10 "$<AND:${is_gnu},$<VERSION_GREATER_EQUAL:${CMAKE_Fortran_COMPILER_VERSION},10.0>>")

target_compile_options(lio-g2g PRIVATE
    # Intel Fortran compiler settings
    $<${is_intel}:-fpp > #-module ${CMAKE_CURRENT_BINARY_DIR}/obj>
    $<$<AND:${is_intel},${is_config_debug}>:
       -g;
       -check bounds;
       -traceback
    >
    $<$<AND:${is_intel},${is_config_release}>:
        -O3;
        -ip;
        -mp1
    >
    $<$<AND:${is_intel},${is_parallel}>:
        -O3;
        -ip;
        -mp1
    >

    # GNU Fortran compiler settings
    $<${is_gnu}: # List of options if GNU
        -cpp;
        #-J${CMAKE_CURRENT_BINARY_DIR}/obj;

        # Version-dependent flag, now nested inside the is_gnu check
        $<${is_gnu_ge10}:-fallow-argument-mismatch>;

        # Config-dependent flags
        $<$<CONFIG:Release>:-O3>;                # Semicolon for list (or not needed if last)

        $<$<AND:$<CONFIG:Debug>,$<EQUAL:${USE_ANALYTICS},1>>: # This implies GNU because it's in the GNU block
            -g;-Wtabs;-fcheck=all
        >; # Semicolon for list
        # Analytics flags (these conditions already imply CONFIG:Debug and is_gnu)
        $<$<AND:$<CONFIG:Debug>,$<EQUAL:${USE_ANALYTICS},2>>: # This implies GNU because it's in the GNU block
            -g;-Wtabs;-fcheck=all # Consider if these are different from the general debug flags
        >;
        $<$<AND:$<CONFIG:Debug>,$<EQUAL:${USE_ANALYTICS},3>>:
            -g;-Wall;-Wextra;-fcheck=all;-fbacktrace;-pedantic;-ffpe-trap=invalid,zero,overflow,underflow
        >;
        $<$<AND:$<CONFIG:Debug>,$<EQUAL:${USE_ANALYTICS},4>>:
            -g;-Wall;-Wextra;-fcheck=all;-fbacktrace;-pedantic;-ffpe-trap=invalid,zero,overflow,underflow
        >
    > # End of is_gnu block
)

# Intel runtime libraries
#set(is_intel_fortran_compiler "$<OR:$<Fortran_COMPILER_ID:IntelLLVM>,$<Fortran_COMPILER_ID:Intel>>")
#set(intel_runtime_libs "iomp5;pthread") # Use semicolon for list of libraries
#target_link_libraries(lio-g2g PRIVATE
#    $<$<CONDITION:${is_intel_fortran_compiler}>:${intel_runtime_libs}>
#)
# Intel runtime libraries
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    target_link_libraries(lio-g2g PRIVATE iomp5 pthread)
endif()
    
# Intel MKL or standard BLAS/LAPACK
if(USE_MKL)
    find_package(MKL REQUIRED)
    target_link_libraries(lio-g2g PRIVATE MKL::mkl)
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    target_link_libraries(lio-g2g PRIVATE BLAS::BLAS LAPACK::LAPACK)
endif()
    
# Handle OpenMP
if(USE_PARALLEL)
    find_package(OpenMP)
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(lio-g2g PRIVATE OpenMP::OpenMP_Fortran)
    endif()
endif()

# CUDA/CUBLAS support
if(USE_CUDA AND (USE_CUBLAS OR USE_MAGMA))

    find_package(CUDAToolkit REQUIRED)
    target_compile_definitions(lio-g2g PRIVATE CUBLAS)
    target_link_libraries(lio-g2g PUBLIC CUDA::cudart CUDA::cublas)
    
    if(USE_MAGMA)
        find_package(MAGMA REQUIRED)
        target_compile_definitions(lio-g2g PRIVATE magma)
        target_link_libraries(lio-g2g PUBLIC ${MAGMA_LIBRARIES})
        if(MAGMA_INCLUDE_DIRS)
            target_include_directories(lio-g2g PRIVATE ${MAGMA_INCLUDE_DIRS})
        endif()
    endif()
    
    # Need fortran.c for CUBLAS interface
    if(EXISTS "/usr/share/doc/nvidia-cuda-toolkit/examples/fortran.c")
        add_library(fortran_cublas OBJECT /usr/share/doc/nvidia-cuda-toolkit/examples/fortran.c)
        target_compile_definitions(fortran_cublas PRIVATE CUBLAS_GFORTRAN)
        target_include_directories(fortran_cublas PRIVATE ${CUDAToolkit_INCLUDE_DIRS})
        target_sources(lio-g2g PRIVATE $<TARGET_OBJECTS:fortran_cublas>)
    endif()
endif()

# Standard libraries
target_link_libraries(lio-g2g PRIVATE m stdc++)

# Define source lists based on Makefile.depends
set(LIO_SOURCES
    # Main sources
    SCF_exter.f90
    TD.f90
    liomain.f90
    dft_get_mm_forces.f90
    dft_get_qm_forces.f90
    get_restrain_energy_forces.f90
    fock_commuts.f90
    init_lio.f90
    lio_finalize.f90
    input_read.f90
    drive.f90
    grid.f90
    packed_storage.f90
    readECP.f90
    generalECP.f90
    errlog.f90
    field_subs.f90
    propagators.f90
    convert_types.f90
    
    # Additional modules
    lionml_data.f90
    lionml_subs.f90
    build_info_cmake.f90
    ghost_atoms.f90
    geometry_optim.f90
    tbdft_subs.f90
    tbdft_data.f90
    dos_subs.f90
    dos_data.f90
    ceed_data.f90
    ceed_subs.f90
    converger_data.f90
    ehrendata.f90
    ehrensubs/ehrensubs.f90
    SCF.f90
    transport.f90
    fockbias_data.f90
    fockbias_subs.f90
    extern_functional.f90
    basis_data.f90
    initial_guess.f90
    shape_data.f90
    mask_ecp.f90
    cubegen.f90
    excited_data.f90
    fstsh_data.f90
    SCF_aux.f90
    translation.f90
    constants_mod.f90
    
)
# Add liomods subdirectories
file(GLOB LIOMODS_SOURCES liomods/*.f90)
list(APPEND LIO_SOURCES ${LIOMODS_SOURCES})

# Add subdirectory sources
# CDFT sources
list(APPEND LIO_SOURCES 
    cdft/cdft_data.f90
    cdft/cdft_subs.f90
)

# Converger sources
list(APPEND LIO_SOURCES 
    converger_subs/converger_subs.f90
)

# Debug tools sources
list(APPEND LIO_SOURCES 
    debug_tools/debug_tools.f90
)

# DFTD3 sources
if(USE_DFTD3)
    list(APPEND LIO_SOURCES 
        dftd3/dftd3.f90
    )
endif()

# Ehren sources
file(GLOB EHREN_SOURCES ehrensubs/*.f90)
list(APPEND LIO_SOURCES 
    ehrensubs/ehrensubs.f90
)

# Excited sources
list(APPEND LIO_SOURCES 
    excitedsubs/excitedsubs.f90
)

# Faint CPU sources
list(APPEND LIO_SOURCES 
    faint_cpu/subm_int1.f90
    faint_cpu/subm_int2.f90
    faint_cpu/subm_int3lu.f90
    faint_cpu/subm_int3mem.f90
    faint_cpu/subm_intsol.f90 
    faint_cpu/subm_int1G.f90
    faint_cpu/subm_intSG.f90
    faint_cpu/subm_intsolG.f90
    faint_cpu/subm_int2G.f90
    faint_cpu/subm_int3G.f90 
    faint_cpu/subm_intECPG.f90
    faint_cpu/subm_intECP.f90
    faint_cpu/subm_intdip.f90
    faint_cpu/subm_intfld.f90
    faint_cpu/faint_cpu.f90
)

# FileIO sources
list(APPEND LIO_SOURCES 
#    fileio/restart_commons.f90
#    fileio/restart_coef.f90
#    fileio/restart_rho.f90
#    fileio/restart_fock.f90
#    fileio/restart_td.f90
    fileio/fileio.f90
)

# FSTSH sources
list(APPEND LIO_SOURCES 
    fstshsubs/fstshsubs.f90
)

# Lennard Jones sources
list(APPEND LIO_SOURCES 
    lennard_jones/lj_switch.f90
    lennard_jones/lj_switch_data.f90
)

# Linear algebra sources
list(APPEND LIO_SOURCES 
    linear_algebra/linear_algebra.f90
)

# Liosubs sources
list(APPEND LIO_SOURCES 
    liosubs/liosubs.f90
)

# Liosubs math sources
list(APPEND LIO_SOURCES 
    liosubs_math/liosubs_math.f90
)

# Mathsubs sources
list(APPEND LIO_SOURCES 
    mathsubs/mathsubs.f90
)

# Properties sources
list(APPEND LIO_SOURCES 
    properties/properties.f90
    properties/properties_data.f90
)

# Typedef sources
list(APPEND LIO_SOURCES 
    typedef_cumat/typedef_cumat.f90
    typedef_operator/typedef_operator.f90
    typedef_sop/typedef_sop.f90
)

# CUBLAS sources (if enabled)
if(USE_CUDA AND (USE_CUBLAS OR USE_MAGMA))
    list(APPEND LIO_SOURCES cublasmath/cublasmath.f90)
endif()

# Python interface module (added last to resolve dependencies)
list(APPEND LIO_SOURCES python_interface.f90)

# Add all sources to the target
target_sources(lio-g2g PRIVATE ${LIO_SOURCES})

# Handle source file properties for optimization levels
# Based on Makefile.options private_flag settings
set_source_files_properties(
    dip.f90 SCFop.f90
    faint_cpu/subm_intfld.f90 faint_cpu/subm_int1G.f90 faint_cpu/subm_int3G.f90
    faint_cpu/subm_intsol.f90 faint_cpu/subm_intsolG.f90
    PROPERTIES COMPILE_OPTIONS "-O1"
)

# Special compile definitions for build_info_cmake.f90
if(TARGET git_version)
    add_dependencies(lio-g2g git_version)
endif()

# Get git version info
execute_process(
    COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT GIT_VERSION)
    set(GIT_VERSION "1.0")
endif()

execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT GIT_HASH)
    set(GIT_HASH "unknown")
endif()

# Pass build information to build_info.f90
message(STATUS "build_info compilation")
message(STATUS "Git version: ${GIT_VERSION}")
message(STATUS "Use CUDA: ${USE_CUDA}")
set_source_files_properties(build_info_cmake.f90 PROPERTIES COMPILE_DEFINITIONS
    "BPASS_VERSION=\"${GIT_VERSION}\";BPASS_GITHASH=\"${GIT_HASH}\";BPASS_CUDA=\"${USE_CUDA}\";BPASS_INTEL=${is_intel};BPASS_PARALLEL=\"${USE_PARALLEL}\";BPASS_PRECISION=\"${USE_FULL_DOUBLE}\";BPASS_ANALYTICS=${USE_ANALYTICS}"
)

# Installation rules
message(STATUS "Installation of lio-g2g")
message(STATUS "Install libraries to: ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "Install executables to: ${CMAKE_INSTALL_BINDIR}")
install(TARGETS lio-g2g
    EXPORT lioTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install Fortran module files
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/obj/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lio
    FILES_MATCHING PATTERN "*.mod"
)

# Export targets
install(EXPORT lioTargets
    FILE lioTargets.cmake
    NAMESPACE lio::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lio
)

# Generate package configuration files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/lioConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/lioConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/lioConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lio
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/lioConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/lioConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lio
)