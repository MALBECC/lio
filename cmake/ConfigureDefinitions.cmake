
# Global definitions
add_definitions(-Dpack -DG0G)

# Optional precision definitions
if (USE_FULL_DOUBLE)
   add_definitions(-DUSE_FULL_DOUBLE=-1)
else()
   add_definitions(-DUSE_FULL_DOUBLE=-2)
   add_definitions(-DTD_SIMPLE)
endif()

# CUDA level definitions.
if (USE_CUDA)
   add_definitions(-DGPU_KERNELS=-1)
   if (USE_CUBLAS) 
      add_definitions(-DCUBLAS)
   endif()
   if (USE_MAGMA)
      add_definitions(-Dmagma)
   endif()
else()
   add_definitions(-DGPU_KERNELS=-2)
endif()

# LIBXC level definitions.
if (USE_LIBXC)
   add_definitions(-DUSE_LIBXC=-1)
   if (USE_LIBXC_GPU)
      add_definitions(-DLIBXC_CPU=-2)
   else()
      add_definitions(-DLIBXC_CPU=-1)
   endif()
else()
   add_definitions(-DUSE_LIBXC=-2)
endif()

# Debugging level definitions.
if (USE_DEBUGGING)
   add_definitions(-DDEBUGGING=-1 -DPRINT_MATRICES)
   if (USE_FULL_CHECKS)
      add_definitions(-DFULL_CHECKS)
   endif()
endif()

# Profiling definitions.
if (USE_PROFILING)
   add_definitions(-DPROFILING=-1)
endif()

# Libint definitions.
if (USE_LIBINT)
   add_definitions(-DUSE_LIBINT=-1)
else()
   add_definitions(-DUSE_LIBINT=-2)
endif()