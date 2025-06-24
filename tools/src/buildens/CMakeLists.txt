# CMakeLists.txt for buildens tool

# Include GNU install directories
include(GNUInstallDirs)

# Add the executable
add_executable(buildens buildens.f90)

# Set output name and directory
set_target_properties(buildens PROPERTIES
  OUTPUT_NAME buildens.x
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
# Installation rules
install(TARGETS buildens
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)