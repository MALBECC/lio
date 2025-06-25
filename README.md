LIO Project
============

LIO is a Quantum Mechanical software based on _Density Functional Theory_ (DFT) and real time _Time Dependent Density Functional Theory_ (TD-DFT).

The most computationally intensive calculations are ported to use graphical processors which support CUDA
(e.g. Nvidia Maxwell, Fermi, Kepler and Tesla families).

# Compiling Project LIO with CMake

This document outlines the process for compiling Project LIO using CMake. It covers basic compilation, prerequisites, and a detailed explanation of all available build options to customize the compilation for different hardware and library configurations.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Basic Compilation](#basic-compilation)
- [Customizing the Build](#customizing-the-build)
- [Available CMake Options](#available-cmake-options)
  - [CUDA Support](#cuda-support)
  - [Intel Compilers and MKL](#intel-compilers-and-mkl)
  - [Parallelism and Precision](#parallelism-and-precision)
  - [External Libraries](#external-libraries)
  - [DFTD3 Support](#dftd3-support)
  - [Profiling and Debugging](#profiling-and-debugging)
  - [Advanced CUDA Options](#advanced-cuda-options)
  - [Testing](#testing)
  - [Library Paths](#library-paths)
- [Compilation Examples](#compilation-examples)

## Prerequisites

Before you begin, ensure you have the following software installed on your system:

- **CMake** (version 3.10 or higher recommended)
- A C++ compiler (e.g., GCC, Clang, or Intel C++ Compiler)
- **(Optional)** NVIDIA CUDA Toolkit if you plan to enable GPU support (`USE_CUDA`).
- **(Optional)** External libraries like Libxc, Libint, Eigen, or MAGMA if you wish to use them.

## Basic Compilation

This will compile the project with the default options (CUDA enabled, DFTD3 enabled).

1.  **Create a build directory:** It's best practice to perform an out-of-source build.
    ```bash
    mkdir build
    cd build
    ```

2.  **Run CMake:** This step configures the project and generates the Makefiles.
    ```bash
    cmake ..
    ```

3.  **Compile the code:** Use `make` to start the compilation. You can use the `-j` flag to specify the number of parallel jobs.
    ```bash
    # Use 4 parallel jobs for compilation
    make -j4
    ```

The compiled binaries will be located in the `build` directory.

## Customizing the Build

You can customize the build by passing options to the `cmake` command using the `-D<OPTION_NAME>=<VALUE>` syntax.

For example, to disable CUDA support and enable OpenMP parallelism, you would run:
```bash
cmake -DUSE_CUDA=OFF -DUSE_PARALLEL=ON ..
```

The following sections detail all available options.

## Available CMake Options

Options are grouped by functionality. Boolean options can be set to `ON` or `OFF`.

### CUDA Support

Options to control GPU acceleration with NVIDIA CUDA.

| Option | Description | Default Value |
| :--- | :--- | :--- |
| `USE_CUDA` | Enable general CUDA support. | `ON` |
| `USE_CUBLAS` | Enable the CUBLAS library. Implies `USE_CUDA=ON`. | `OFF` |
| `USE_MAGMA` | Enable the MAGMA library. Implies `USE_CUDA=ON` and `USE_CUBLAS=ON`. | `OFF` |

### Intel Compilers and MKL

Options for using the Intel toolchain.

| Option | Description | Default Value |
| :--- | :--- | :--- |
| `USE_INTEL_COMPILER` | Use Intel C++ Compiler (`icpc`) instead of the default GNU/Clang compiler. | `OFF` |
| `USE_MKL` | Link against the Intel Math Kernel Library (MKL). | `OFF` |

### Parallelism and Precision

Control multi-threading and numerical precision.

| Option | Description | Default Value |
| :--- | :--- | :--- |
| `USE_PARALLEL` | Enable parallel processing using OpenMP. | `OFF` |
| `USE_FULL_DOUBLE` | Use double precision for all calculations (CPU and GPU). | `OFF` |
| `USE_CPU_DOUBLE` | Use double precision for CPU-only calculations. | `ON` |
| `AINT_MP` | Enable multi-precision in analytic integrals. | `OFF` |

### External Libraries

Enable support for external scientific libraries. You may need to provide paths to these libraries (see [Library Paths](#library-paths)).

| Option | Description | Default Value |
| :--- | :--- | :--- |
| `USE_LIBXC_CPU` | Enable the CPU version of the Libxc library. | `OFF` |
| `USE_LIBXC_GPU` | Enable the GPU version of the Libxc library. | `OFF` |
| `USE_LIBINT` | Enable the Libint library. | `OFF` |

### DFTD3 Support

| Option | Description | Default Value |
| :--- | :--- | :--- |
| `USE_DFTD3` | Enable support for Grimme's D3 dispersion correction. | `ON` |

### Profiling and Debugging

| Option | Description | Default Value |
| :--- | :--- | :--- |
| `USE_ANALYTICS` | Set the debug level for analytics. `0` = off. Higher values increase verbosity. | `0` |

### Advanced CUDA Options

Fine-tune the CUDA compilation process.

| Option | Description | Default Value |
| :--- | :--- | :--- |
| `USE_CUDA_ARCH` | A semicolon-separated list of CUDA architectures to compile for (e.g., "75;86"). | `"35;52;61;75"` |
| `CUDA_USE_FAST_MATH` | Enable the `--use_fast_math` flag in `nvcc` for potentially faster but less precise math operations. | `ON` |
| `CPU_RECOMPUTE` | Recompute CPU energy (useful for debugging and verification). | `OFF` |
| `CUDA_VERBOSE` | Enable verbose output during CUDA compilation. | `OFF` |
| `CUDA_PTX` | Enable verbose PTX compilation output. | `OFF` |
| `CUDA_REGCOUNT` | Enable verbose CUDA register count output. | `OFF` |

### Testing

| Option | Description | Default Value |
| :--- | :--- | :--- |
| `BUILD_TESTING` | Build the project's test suite. | `OFF` |
| `TEST_CHECKS_ONLY`| If testing is enabled, only perform checks without running full tests. | `OFF` |

### Library Paths

Use these variables to specify the installation paths for external libraries if they are not in a standard system location.

| Variable | Description | Default Value |
| :--- | :--- | :--- |
| `LIBXC_HOME_CPU` | Path to the root of the LIBXC CPU installation. | `""` |
| `LIBXC_HOME_GPU` | Path to the root of the LIBXC GPU installation. | `""` |
| `LIBINT_HOME` | Path to the root of the LIBINT installation. | `""` |
| `EIGEN_HOME` | Path to the Eigen header library installation. | `""` |
| `MAGMA_ROOT` | Path to the root of the MAGMA installation. | `""` |

## Compilation Examples

Here are a few common compilation scenarios. Always run these commands from a clean `build` directory.

### Example 1: Standard CPU-only Build with OpenMP 

This configuration disables all GPU features and enables OpenMP for parallel CPU execution using GNU compilers.

```bash
cd build
export LIO_PREFIX=/desired/lio/install/directory # this is different from the source and build directories.
CC=gcc CXX=g++ FC=gfortran cmake -DUSE_CUDA=OFF \
      -DUSE_PARALLEL=ON -DCMAKE_INSTALL_PREFIX=${LIO_PREFIX}\
      ..
make -j8
ctest # Testing can take several hours.
make install # If provided CMAKE_INSTALL_PREFIX=/path/to/install
```

### Example 2: High-Performance GPU Build for Modern NVIDIA GPUs

This build targets modern NVIDIA architectures (Turing, Ampere), enables MAGMA for optimized linear algebra, and uses the GPU version of Libxc.
This is the default build, so appart from compilers and installation path no extra options are necessary.

```bash
cd build
export LIO_PREFIX=/desired/lio/install/directory # this is different from the source and build directories.
CC=gcc CXX=g++ FC=gfortran cmake .. -DCMAKE_INSTALL_PREFIX=${LIO_PREFIX} 
make -j8
ctest # Testing can take several hours.
make install # If provided CMAKE_INSTALL_PREFIX=/path/to/install
```

### Example 3: High-Performance GPU build with PBE0 support via LibINT library.

If you were successful in compiling the basic (default) version, you may want 
to try compiling lio with LibINT support (which allows for the use of PBE0
DFT functional) or with LibXC (standard CPU version) or its in-house GPU version.

CUDA support, parallel CPU via OpenMP, LibINT and testing activated. 
This build has support for PBE0 DFT functional using LIO's very fast DFT engine
and exact exchange integrals calculated with LibINT. 

This build therefore requires LibINT 2.6.0 which should be compiled from source previous to
this step. Its source code can be downloaded from the following website:
https://github.com/evaleev/libint/releases/tag/v2.6.0
LibInt requires EIGEN package as prerequisite (can be installed from
your linux distro repositories). In Ubuntu Eigen can be installed with:
sudo apt install libeigen3-dev

```bash
export LIBINT_HOME=/path/to/libint/installation/directory
CC=gcc CXX=g++ FC=gfortran cmake .. -DCMAKE_INSTALL_PREFIX=${LIO_PREFIX} \
                                    -DUSE_LIBINT=ON \
                                    -DUSE_FULL_DOUBLE=ON 

#If no errors are found during the cmake configure step run make:
make

If lio compiles without errors we test the compilation (the testing phase can 
take up to several hours depending on the hardware available.
ctest

If everything goes well we can install into the defined installation path:
make install
```

### Example 4: Lio with LibXC and LibINT support
LibXC external library allows more flexibility in the DFT functional selection
at the expense of a heavy performance penalty (compared to Lio's native DFT engine
which only supports the PBE functional).
LibXC support requires LibINT, so we recomend trying this build after being sure
the LibINT only build compiles and works fine.
Just as before you should download an compile libxc version 5.0.0 from the web:
https://gitlab.com/libxc/libxc/-/archive/5.0.0/libxc-5.0.0.tar.bz2
For more information visit Lio's wiki page:
https://github.com/MALBECC/lio/wiki/LIO-installation
Once LibXC and LibINT are correctly compiled you can try to compile lio:

```bash
export LIBXC_HOME_CPU=/path/to/libxc/installation # Define path to libxc in your system
export LIBINT_HOME=/path/to/libint/installation   # Define path to libint in your system
CC=gcc CXX=g++ FC=gfortran cmake ..  -DCMAKE_INSTALL_PREFIX=${LIO_PREFIX} \
                                     -DUSE_LIBINT=ON \
                                     -DUSE_FULL_DOUBLE=ON \
                                     -DUSE_LIBXC_CPU=ON 
#If no errors are found during the cmake configure step run make:
make
#If lio compiles without errors we test the compilation (the testing phase can 
#take up to several hours depending on the hardware available.
ctest
#If everything goes well we can install into the defined installation path:
make install
```

### Example 5: Debug Build with Testing

This configuration builds for debugging and includes the test suite.

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DBUILD_TESTING=ON \
      -DUSE_ANALYTICS=2 \
      ..
make
# After building, run the tests
ctest
```

OLDSTYLE COMPILATION 
---------------------

NOTICE: This is the old way of compiling Lio. Although it will work fine, it is deprecated and will go unsupported in future releases.
The program can be compiled using the make command. The following options can be used to modify
compilation. For example, the following compiles the GPU kernels:

```
make cuda=1 cpu=0
```

When using Intel's ICC/MKL or NVIDIA's CUDA libraries, add them to LD\_LIBRARY\_PATH environment variable before compilation. Available options for compilation include:

* _cpu_: compile CPU kernels (default = 0).

* _cuda_: compile GPU kernels (when = 1) and CUBLAS subroutines (when = 2). Used by default (=1).

* _intel_: use INTEL compilers (when = 1) and INTEL MKL (when = 2). Not used by default (= 0).

* _analytics_: Enables diferent levels of debug information and profiling (default = 0, max = 4).

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
  5. In order for lio to compile with libxc, you'll need to create and export the following variable
```
LIBXC_HOME_CPU=points to the path where the libxc lib and include folders where installed in the file system ( on cpu mode )
LIBXC_HOME_GPU=points to the path where the libxc lib and include folders where installed in the file system ( on gpu mode )
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
