LIO Project

COPY SOMETHING HERE.

REQUIREMENTS
============

MKL Intel Library
Intel C++ and Fortran Compiler (can be obtained a non-commercial license).
NVIDIA CUDA (even for compiling CPU version).

COMPILATION
===========

Enter g2g folder, and makefile there.
G2G Compilation options, call following to make command, for example: make time=1

-gcc: use the GNU C Compiler suite (default is the intel suite).

-dbg: enable debugging information

-cuda: enable the compilation for NVidia GPU cards

-time: enables the timers to obtain detailed timing information from different parts of the program.

-nosync: VERIFY IN CODE

-profile: enabling gprof profiling information.

-cpu_recompute: recomputes=0 mantain in memory the value of the functions in each point (more memory is used and a 10% less execution time).

-static: generate the static library.

-full_double: generate the application using full double precision instead of mixed precision (which is the default).


INSTALLATION
