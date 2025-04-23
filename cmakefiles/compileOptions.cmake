# Compile Options
option(CPU "LIO will be compiled to run on CPUs only" ON)
option(CUDA "LIO will be compiled with CUDA support to run on GPUs and CPUs" OFF)
option(EXTERNAL "LIO will be compiled with support for Libxc and Libint libraries" OFF)
option(COMPILE_EXTERNAL "LIO will compile libxc and libint libraries and then 
                         link them to LIO libraries and executables" OFF)
option(DOUBLE "LIO will be compiled using Double Precision" OFF)
option(WARNING "LIO will be compiled using Warning Flags" OFF)

if (CUDA) 
   message(FATAL_ERROR "For the moment no CUDA support is provided with CMake.")
endif()