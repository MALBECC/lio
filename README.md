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
Currently there are two different implementations of the __Libxc__ library, the original version ([here](http://www.tddft.org/programs/libxc/download/)) who
runs all the functionals in __cpu__ and the modified version ([here](https://gitlab.com/eduarditoperez/libxc/tree/cuda-integration)) who
runs some functionals in __gpu__.
This version of Lio has support for both libraries depending on wich one you'll choose to use at runtime.
You can use this version of Lio with both, one of none of the __Libxc__ libraries. We'll recommend to have both installed.

In order to compile lio with libxc, follow the steps detailed below:

  1. Download the libxc [cpu](http://www.tddft.org/programs/libxc/download/) library. We recomend you to use libxc version 4.0.4.
  2. Compile and Install the __cpu__ library (follow Libxc install instructions).
  3. Download the libxc [gpu](https://gitlab.com/eduarditoperez/libxc/tree/cuda-integration) library.
  4. Compile and Install the library (follow the Libxc-gpu install [instructions](https://github.com/MALBECC/lio/wiki/Libxc-with-CUDA-support-installation-guide#instalation-guide)). Both the __gpu__ and __cpu__ libraries had to be installed in the same folder.
  5. In order for lio to compile with libxc, you'll need to create and export two environment variables in `LD_LIBRARY_PATH`
```
LIBXC_LIBS=points to the path where the libxc libaries where installed in the file system.
LIBXC_INCLUDES=points to the path where the libxc include files where installed in the file system.
```
  6. Go to the `Libxc` (gpu) installation folder and copy the next files into the `LIBXC_INCLUDES` folder defined in the step 5.
```
src/xc.h
src/xc_funcs.h
src/xc_funcs_removed.h
xc_version.h
```
  7. Libxc has 3 compilation options for lio, those are
```
libxc=0 - No libxc (DEFAULT - for backwards compatibility with lio)
libxc=1 - Use libxc in CPU mode
libxc=2 - Use libxc in GPU mode
```
  8. If you want to compile lio with libxc in GPU mode, type the following command
```
make cuda=1 libxc=2
```
  9. To validate the instalation, go to the `integration-test` folder located in `lio/test/` and run the command `make`, this will
compile and execute the integration test. After the execution of the test phase, you should see in the console:
```
Running gpu test...
gpu integration-test01:  0
gpu integration-test02:  0
gpu integration-test03:  0
Running cpu test...
cpu integration-test:  0
```
The `0` after each test means that the test ran without errors.

  9. To run the simulations using the functionals from libxc you'll have to add the following variables in the `*****.in` files:
```
file: agua.in

...
use_libxc=t
ex_functional_id=XXX
ec_functional_id=XXX
...
```
where `ex_functional_id` is the id for the exchange functional from libxc and `ec_funcional_id` is the id
for the correlation functional from libxc. You can see the list of available functionals for [gpu](https://github.com/MALBECC/lio/wiki/Libxc-available-functionals-for-GPU-version#functionals-for-gpu-version)
and the list of available functionals for [cpu](https://github.com/MALBECC/lio/wiki/Libxc-available-functionals-for-CPU#functionals-for-cpu-version).
Bare in mind that only the GGA functional's family are supported in this version of libxc for lio.


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
