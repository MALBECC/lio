LIO Project

Lio is a Quantum Mechanical software based on Density Functional Theory (DFT) and real time-time dependent Density Functiona Theory (TD-DFT).
The most computer consumer parts are ported to use the graphical prossesors who suport CUDA.

REQUIREMENTS
============

MKL Intel Library
Intel C++ and Fortran Compiler (can be obtained a non-commercial license).
NVIDIA CUDA (even for compiling CPU version).

COMPILATION
===========

Enter g2g folder, and makefile there.
G2G Compilation options, call following to make command, for example: make cuda=1

-dbg: enable debugging information

-cuda: enable the compilation for NVidia GPU cards

-time: enables the timers to obtain detailed timing information from different parts of the program.

-nosync: VERIFY IN CODE

-profile: enabling gprof profiling information.

-cpu_recompute: recomputes=0 mantain in memory the value of the functions in each point (more memory is used and a 10% less execution time).

-static: generate the static library.

-full_double: generate the application using full double precision instead of mixed precision (which is the default).


INSTALLATION
