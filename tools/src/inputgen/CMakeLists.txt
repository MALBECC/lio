# CMakeLists.txt for the inputgen tool
# Located in: tools/src/inputgen/CMakeLists.txt

cmake_minimum_required(VERSION 3.18)

# By being in a subdirectory, project() is optional unless you want to version it separately
# or need specific project-level settings here. For a simple tool, it's often omitted.
# project(InputGenTool)

# Enable Fortran language if not already enabled by a parent CMakeLists.txt
# If tools/CMakeLists.txt or the root CMakeLists.txt already did enable_language(Fortran),
# this is not strictly necessary but harmless.
if(NOT CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
endif()

# --- Fortran Module Handling ---
# Modern CMake handles Fortran module dependencies automatically for sources
# compiled within the SAME target. The Makefile's L0, L1, L2 logic implies
# a specific compilation order because modules from L0 are needed by L1, etc.
# CMake usually figures this out if all sources are part of one target.

# If these source files generate modules that are ONLY used by other files
# within this same inputgen.x executable, then simply adding them all to
# one add_executable command is often sufficient.

# However, the Makefile structure suggests a clear layering.
# If you want to be explicit or if these modules were ever intended to be
# part of a separate library (though not the case here), you might create
# INTERFACE libraries. For a single executable, it's usually simpler.

# Let's define the sources based on the Makefile's grouping.
# All sources are in the 'src' subdirectory relative to this CMakeLists.txt
# The vpath %.f90 src in the Makefile tells us this.
# We assume this CMakeLists.txt is in tools/src/inputgen/
# and the .f90 files are in tools/src/inputgen/src/

set(INPUTGEN_SOURCES_L0
    src/keywords.f90
    src/auxmod_subs.f90
)

set(INPUTGEN_SOURCES_L1
    src/parser_gaussian.f90
    src/parser_lio.f90
    # Potentially src/parser_orca.f90 if it's part of L1
)

set(INPUTGEN_SOURCES_L2
    src/data_control.f90
)

set(INPUTGEN_MAIN_SOURCE
    src/inputgen.f90
)

# Combine all sources for the executable
set(INPUTGEN_ALL_SOURCES
    ${INPUTGEN_SOURCES_L0}
    ${INPUTGEN_SOURCES_L1}
    ${INPUTGEN_SOURCES_L2}
    ${INPUTGEN_MAIN_SOURCE}
)

# --- Executable Target ---
add_executable(inputgen.x ${INPUTGEN_ALL_SOURCES})

# --- Compiler and Linker Flags (Example) ---
# These would typically be inherited from the parent project (LIO)
# or set globally. If you need specific flags ONLY for inputgen.x:

# target_compile_options(inputgen.x PRIVATE
#     # Add any specific Fortran compile flags for inputgen if different from global
#     # e.g., $<$<COMPILE_LANGUAGE:Fortran>:-O2>
# )

# target_link_libraries(inputgen.x PRIVATE
#     # Add any specific libraries inputgen.x needs to link against
#     # e.g., some_utility_library
# )

# --- Fortran Module Directory ---
# CMake handles .mod file output and dependencies automatically.
# By default, .mod files for a target are placed in a subdirectory of
# CMAKE_CURRENT_BINARY_DIR related to the target (e.g., CMakeFiles/inputgen.x.dir/)
# and this location is automatically added to the include path for other
# sources within the same target.
# You can set Fortran_MODULE_DIRECTORY if you want to explicitly control it,
# but it's often not needed for self-contained executables.
#
# For example, if you wanted all .mod files from this target in 'tools/build/inputgen_mods':
# set_target_properties(inputgen.x PROPERTIES
# Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/inputgen_mods"
# )
# And then ensure that directory is searched:
# target_include_directories(inputgen.x PRIVATE "${CMAKE_BINARY_DIR}/inputgen_mods")
# But again, for a single target, CMake's default usually works fine.

# --- Installation ---
# The Makefile copies inputgen.x to ../../inputgen
# Assuming this CMakeLists.txt is in tools/src/inputgen, then ../../inputgen
# would be tools/inputgen.
# If the main LIO project handles installation, this might be done differently.
# Let's assume we want to install it to the main project's binary directory for now,
# or to a standard install location if the LIO project defines CMAKE_INSTALL_BINDIR.

# Option A: Copy to a specific location relative to the top-level build/source.
# This mimics the `cp inputgen.x ../../inputgen` behavior if inputgen.x
# should be available in the `tools` directory of the build or source.
# This is more like a build step utility rather than a formal installation.

# If the main project (LIO) will install tools, then a proper `install(TARGETS ...)`
# command should be used, possibly in `tools/CMakeLists.txt`.

# For this specific tools/src/inputgen/CMakeLists.txt, let's focus on building the executable.
# The installation or copying can be handled by the parent `tools/CMakeLists.txt`.
# If this `inputgen.x` is meant to be an installed tool:
include(GNUInstallDirs)
install(TARGETS inputgen.x
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # e.g., <prefix>/bin
)

# --- Custom Command for Copying (if not using install()) ---
# If you strictly want to replicate the `cp inputgen.x ../../inputgen`
# and this `../../inputgen` refers to `<LIO_ROOT>/tools/inputgen.x`,
# you'd do this in the `tools/CMakeLists.txt` file, not here.
# For example, in `tools/CMakeLists.txt`:
# add_custom_command(TARGET inputgen.x POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy
#             $<TARGET_FILE:inputgen.x>
#             ${CMAKE_SOURCE_DIR}/tools/inputgen # Or wherever ../../inputgen resolves to
#     COMMENT "Copying inputgen.x to tools directory"
# )
# However, standard `install()` is preferred for actual deployment.

# --- Clean command ---
# `make clean` is handled automatically by CMake-generated Makefiles.
# The `rm -f ../../inputgen` part of the clean would be handled if the copy
# was part of a build step that CMake knows about (e.g., via add_custom_target
# if the copy was more complex or needed to be a separate build target).
# If you're just using `install()`, the uninstall target would handle removal from
# the install prefix.
