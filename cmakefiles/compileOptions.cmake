# Compile Options
option(CPU "LIO will be compiled using CPUs only" ON)
option(GPU "LIO will be compiled using GPUs and CPUs" OFF)
option(EXTERNAL "LIO will be compiled using Libxc and Libint libraries" OFF)
option(COMPILE_EXTERNAL "LIO will compile libxc and libint and then 
                         will make the link with LIO" OFF)
option(DOUBLE "LIO will be compiled using Double Precision" OFF)
option(WARNING "LIO will be compiled using Warning Flags" OFF)

if (GPU) 
   message(FATAL_ERROR "For the moment no GPU compilations is provided in cmake")
endif()