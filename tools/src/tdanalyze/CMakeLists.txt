# CMakeLists.txt for tdanalyze tools

# Add the executables
add_executable(tdanalyze tdanalyze.f90)
add_executable(average_spectra average_spectra.f90)

# Set output properties
set_target_properties(tdanalyze PROPERTIES
  OUTPUT_NAME tdanalyze.x
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(average_spectra PROPERTIES
  OUTPUT_NAME average_spectra.x
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Installation rules
install(TARGETS tdanalyze average_spectra
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)