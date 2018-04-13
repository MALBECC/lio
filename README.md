LIO Project
============

LIO is a Quantum Mechanical software based on _Density Functional Theory_ (DFT) and real time _Time Dependent Density Functional Theory_ (TD-DFT).

The most computationally intensive calculations are ported to use graphical processors which support CUDA
(e.g. Nvidia Maxwell, Fermi, Kepler and Tesla families).

REQUIREMENTS
------------

* LAPACK or INTEL MKL.
* GNU or INTEL C++ and Fortran Compiler.
* NVIDIA CUDA (for the GPU kernels).
* GNU Make.
* Libxc for LIO (optional).

COMPILATION
------------

The program can be compiled using the make command. The following options can be used to modify
compilation. For example, the following compiles the GPU kernels:

```
make cuda=1 cpu=0
```

When using Intel's ICC/MKL or NVIDIA's CUDA libraries, add them to LD\_LIBRARY\_PATH environment variable before compilation. Available options for compilation include:

* _cpu_: compile CPU kernels (default = 0).

* _cuda_: compile GPU kernels (when = 1) and CUBLAS subroutines (when = 2). Used by default (=1).

* _intel_: use INTEL compilers (when = 1) and INTEL MKL (when = 2). Not used by default (= 0).

* _analytics_: Enables diferent levels of debug information and profiling (default = 0, max = 3).

* _precision_: When precision = 1, compile everything in double precision (default = 0, hybrid precision).

* _libxc_: compile the application to use libxc library. Requires libxc for lio installed.

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

INSTALLATION WITH LIBXC
-----------------------

  1. Download the libxc library for lio from gitlab (account needed) ([here](https://gitlab.com/eduarditoperez/libxc.git)).
  2. Compile and Install the library (follow the Readme.me instructions).
  3. In order for lio to compile with libxc, you'll need to create two variables in LD_LIBRARY_PATH
```
LIBXC_LIBXS=points to the path where the libxc libaries where installed in the file system.
LIBXC_INCLUDES=points to the path where the libxc include files where installed in the file system.
```
  4. Libxc has 3 compilation options for lio, these are
```
libxc=0 - No libxc (for backwards compatibility with lio)
libxc=1 - Use libxc in GPU mode
libxc=2 - Use libxc in CPU mode
```
  5.If you want to compile lio with libxc in GPU mode, type the following command
```
make cuda=1 libxc=1
```
  6. To run the simulations using the functionals from libxc you'll have to add the following variables in the *****.in files:
```
...
use_libxc=t
ex_functional_id=XXX
ec_functional_id=XXX
...
```
where ex_functional_id is the id for the energy-exchange functional from libxc and ec_funcional_id is the id
for the energy-correlation functional from libxc. You can see the list of available functionals ([here](https://gitlab.com/libxc/libxc/wikis/Functional-list-4.0.4))
or in see the funcs_key.c file in the src folder of libxc. Bare in mind that only the GGA functional's family are supported in
this version of libxc for lio.

TESTS
-----

The test suite can be ran from the tests directory, each subfolder contains a "correr.sh" script which performs the test.


CONTRIBUTING
------------

Before contributing, make sure you have set up the git hooks for the project, and do read the wiki and workflow of the project.

PUBLICATIONS
------------

1. Matías A. Nitsche, Manuel Ferreria, Esteban E. Mocskos and Mariano C. González Lebrero, _GPU Accelerated Implementation of Density Functional Theory for Hybrid QM/MM Simulations_. J. Chem. Theory Comput., 2014, 10 (3), pp 959–967.

2.  Uriel N. Morzan, Francisco F. Ramírez, M. Belén Oviedo, Cristián G. Sánchez, Damián A. Scherlis and Mariano C. González Lebrero, _Electron dynamics in complex environments with real-time time dependent density functional theory in a QM-MM framework_. J. Chem. Phys. 140, 164105 (2014).
