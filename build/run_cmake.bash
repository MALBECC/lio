#!/bin/bash

# This script provides the user with some example CMake usage.
# It is meant as an example of CMake compilation.
# Provided options are meant to be edited as necessary.
#
# For more information please visit lio's wiki:
# https://github.com/MALBECC/lio/wiki/LIO-installation
# 
# Edit the following two variable for the installation
# directory (LIO_INSTALL) and the source code 
# directory (LIO_SOURCE). Replace the example paths
# by the correct ones for your system.
#
# After successfull (no errors) run of CMake do
# make # Actual compilation
# ctest # Testing (can take a few hours)
# make install # Installation of the executables and libraries


export LIO_INSTALL=/opt/lio
export LIO_SOURCE=${HOME}/src/lio

# By default Lio is configured with CUDA support, 
# parallel CPU via OpenMP, and testing activated.
# This command also uses GNU compilers.
CC=gcc CXX=g++ FC=gfortran cmake ${LIO_SOURCE} -DCMAKE_INSTALL_PREFIX=${LIO_PREFIX} 2&>1 tee cmake.log

# ALTRNATIVE (MORE ADVANCED) BUILDS
#------------------------------------------------------------------------------
# If you were successful in compiling the basic (default) version, you may want 
# to try compiling lio with LibINT support (which allows for the use of PBE0
# DFT functional) or with LibXC (standard CPU version) or its in-house GPU version.

# LIBINT SUPPORT
#------------------------------------------------------------------------------
# CUDA support, parallel CPU via OpenMP, LibINT and testing activated. 
# This build has support for PBE0 DFT functional using LIO's very fast DFT engine
# and exact exchange integrals calculated with LibINT. 
#
# This build therefore requires LibINT 2.6.0 which should be compiled from source previous to
# this step. Its source code can be downloaded from the following website:
# https://github.com/evaleev/libint/releases/tag/v2.6.0
# LibInt requires EIGEN package as prerequisite (can be installed from
# your linux distro repositories). In Ubuntu Eigen can be installed with:
# sudo apt install libeigen3-dev
#
# Uncomment the following lines to set up this build
# export LIBINT_HOME=/path/to/libint/installation/directory
# export LD_LIBRARY_PATH=$LIBINT_HOME/lib:$LD_LIBRARY_PATH
# CC=gcc CXX=g++ FC=gfortran cmake ${LIO_SOURCE} 
#                                  -DCMAKE_INSTALL_PREFIX=${LIO_PREFIX} \
#                                  -DUSE_LIBINT=ON \
#                                  -DUSE_FULL_DOUBLE=ON \
#                                  -DBUILD_TESTING=ON \
#                                  -DUSE_PARALLEL=ON 2&>1 tee cmake.log
#
# If no errors are found during the cmake configure step run make:
# make
#
# If lio compiles without errors we test the compilation (the testing phase can 
# take up to several hours depending on the hardware available.
# ctest
#
# If everything goes well we can install into the defined installation path:
# make install
#
#
# LIBXC (CPU VERSION) SUPPORT
#------------------------------------------------------------------------------
# LibXC external library allows more flexibility in the DFT functional selection
# at the expense of a heavy performance penalty (compared to Lio's native DFT engine
# which only supports the PBE functional).
# LibXC support requires LibINT, so we recomend trying this build after being sure
# the LibINT only build compiles and works fine.
# Just as before you should download an compile libxc version 5.0.0 from the web:
# https://gitlab.com/libxc/libxc/-/archive/5.0.0/libxc-5.0.0.tar.bz2
# For more information visit Lio's wiki page:
# https://github.com/MALBECC/lio/wiki/LIO-installation
# Once LibXC and LibINT are correctly compiled you can try to compile lio:
# 
# export LIBXC_HOME_CPU=/path/to/libxc/installation # Define path to libxc in your system
# export LIBINT_HOME=/path/to/libint/installation   # Define path to libint in your system
# CC=gcc CXX=g++ FC=gfortran cmake ${LIO_SOURCE} 
#                                  -DCMAKE_INSTALL_PREFIX=${LIO_PREFIX} \
#                                  -DUSE_LIBINT=ON \
#                                  -DUSE_FULL_DOUBLE=ON \
#                                  -DUSE_LIBXC_CPU=ON 2&>1 tee cmake.log
#
# If no errors are found during the cmake configure step run make:
# make
#
# If lio compiles without errors we test the compilation (the testing phase can 
# take up to several hours depending on the hardware available.
# ctest
#
# If everything goes well we can install into the defined installation path:
# make install
#

