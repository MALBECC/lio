# CTestCustom.cmake - This file is automatically included by CTest.

# Get the locations of the library output directories.
# Note: This uses variables, not generator expressions, because this script
# is evaluated by CTest using a pre-generated configuration file.
set(G2G_LIB_DIR "${CMAKE_BINARY_DIR}/g2g") # Adjust if you set CMAKE_LIBRARY_OUTPUT_DIRECTORY
set(LIO_G2G_LIB_DIR "${CMAKE_BINARY_DIR}/lioamber")

# Build the full path
set(FULL_LIB_PATH "${G2G_LIB_DIR}:${LIO_G2G_LIB_DIR}")

# Get the old path to append to it
set(OLD_LD_PATH "$ENV{LD_LIBRARY_PATH}")
if(OLD_LD_PATH)
  set(FULL_LIB_PATH "${FULL_LIB_PATH}:${OLD_LD_PATH}")
endif()

# Set the environment variable for the ctest run
set(ENV{LD_LIBRARY_PATH} "${FULL_LIB_PATH}")

# You can print a message to verify
message(STATUS "CTestCustom: Setting LD_LIBRARY_PATH to ${FULL_LIB_PATH}")
