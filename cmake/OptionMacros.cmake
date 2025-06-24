# This file contains macros for setting options for LIO project.
#Macro for printing an option in a consistent manner
#
#Syntax: print_option(<option to print> <was specified>)
#
macro(print_option variable default)
if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    message(STATUS "Setting (unspecified) option ${variable}: ${default}")
else()
    message(STATUS "Setting option ${variable}: ${${variable}}")
endif()
endmacro()


macro(option_with_default variable msge default)
print_option(${variable} ${default})
if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
   set(${variable} ${default} CACHE STRING ${msge} FORCE)
endif()
endmacro(option_with_default)


# Stops the code if the source and the build folders are identical.
#
function(lio_ensure_out_of_source_build)

  get_filename_component(srcdir "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)

  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR
      "It is not allowed to configure and build DFTB+ from its source folder. Please, create a \
separate build directory and invoke CMake from that directory. See the INSTALL.rst file for \
detailed build instructions.")
  endif()

endfunction()
