# CMakeLists.txt for tools directory
# Enable Fortran language
enable_language(Fortran)

# Add subdirectories for specific tools
add_subdirectory(src/buildens)
add_subdirectory(src/inputgen)
add_subdirectory(src/tdanalyze)