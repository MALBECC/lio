LIO Project
============

LIO is a Quantum Mechanical software based on _Density Functional Theory_ (DFT) and real time _Time Dependent Density Functional Theory_ (TD-DFT).

The most computationally intensive calculations are ported to use graphical processors which support CUDA
(e.g. Nvidia Maxwell, Fermi, Kepler and Tesla families).

REQUIREMENTS
------------

* Intel MKL (Math Kernel Library).
* Intel C++ and Fortran Compiler (can be obtained with a non-commercial license).
* NVIDIA CUDA (for compiling GPU version).
* GNU Make.
* Python 2.7 (for test scripts)

COMPILATION
------------

The program can be compiled using the make command. The following options can be used to modify
compilation. For example, the following compiles the GPU kernels:

```
make cuda=1 gpu=1
```

When using Intel's ICC/MKL or NVIDIA's CUDA libraries, add them to LD\_LIBRARY\_PATH environment variable before compilation. Available options for compilation include:

* _cpu_: compile CPU version.

* _gpu_: compile GPU version.

* _dbg_: enable debugging information. It also enables several asserts which degrade performance.

* _cpu_: utilized non-optimized kernels in CPU instead of GPU (not recommended, only to compare running the same code in CPU and GPU).

* _time_: enables the timers to obtain detailed timing information from different parts of the program.

* _profile_: enabling gprof profiling information.

* *cpu_recompute*: recomputes=0 mantains in memory the value of the functions for each point (more memory is used but execution time goes down by around 10%). Only used for the CPU kernels.

* *full_double*: generate the application using full double precision instead of mixed precision (which is the default).

INSTALLATION
------------

Compilation will produce two dynamic libraries, which should be added to LD\_LIBRARY\_PATH and PATH environment variables.

  1. g2g/libg2g.so
  2. lioamber/liblio-g2g.so

Then set LIOHOME environment variable, pointing to LIO location.

INSTALLATION WITH AMBER
-----------------------

  1. Compile LIO as indicated above. 
  2. Be sure to check (or edit if needed) the /src/configure2 file in AMBER so that liolibs variable correctly points to LIO library folders.
  3. Configure and compile AMBER with the -lio option (see Amber installation instructions).
  4. Done!

INSTALLATION WITH GROMACS
-------------------------

NOTE: GROMACS is not yet officially supported on the other side, but we have our own up-to-date Gromacs repository with the files needed.
  1. Compile LIO as indicated above.
  2. Compile GROMACS as usual, but changing compilation flags (see Gromacs installation instructions):
```
cmake -DGMX_GPU=0 -DGMX_THREAD_MPI=0 -DGMX_QMMM_PROGRAM="lio" -DLIO_LINK_FLAGS="-L/usr/lib -L/usr/lib64 -L/PATHTOLIOLIBRARIES -lg2g -llio-g2g"
```
  3. Done!

TESTS
-----

To run the test suite, you need to install python 2. It is already present in most Linux distributions.

The test suite can be ran from the tests directory and running

```
  ./run_tests.py
```

The first argument to the run\_tests.py program is a regular expression ([Python Documentation](https://docs.python.org/2/howto/regex.html)) which matches the folder names in the test
directory. For example, running from the tests directory:

```
  ./run_tests.py --filter_rx GPU
```

runs only the tests for GPU. For more options run

```
  ./run_tests.py --help
```

CONTRIBUTING
------------

Before contributing, make sure you have set up the git hooks for the project. That
can be done either by running a clean compile with make, or by executing

```
  make hooks
```

PUBLICATIONS
------------

1. Matías A. Nitsche, Manuel Ferreria, Esteban E. Mocskos and Mariano C. González Lebrero, _GPU Accelerated Implementation of Density Functional Theory for Hybrid QM/MM Simulations_. J. Chem. Theory Comput., 2014, 10 (3), pp 959–967.

2.  Uriel N. Morzan, Francisco F. Ramírez, M. Belén Oviedo, Cristián G. Sánchez, Damián A. Scherlis and Mariano C. González Lebrero, _Electron dynamics in complex environments with real-time time dependent density functional theory in a QM-MM framework_. J. Chem. Phys. 140, 164105 (2014).
