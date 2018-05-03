# LAMA (Version 3.0.0 Chimborazo)

For building LAMA on your system please consider the installation requirements below and install required packages. 

A description on how to build LAMA is given afterwards or refer to our homepage www.libama.org --> Documentation --> Installation

## Requirements for building LAMA

Required:
 * CMake >= 2.8.8
 * C/C++ compiler (optionally: OpenMP 3.0, recommended: capable of C++11)
 * Boost-Library (header, optional for tests: unit-test-framework )

Recommended:
 * BLAS- and LAPACK-Library (Intel MKL, BLAS)
 * Nvidia GPU: CUDA >= 7.0 (with cuBLAS, cuSPARSE)
 * Intel Xeon Phi: Intel MKL
 * MPI

Optional:
 * Metis/ParMetis
 * Java
 * Documentation:
   - Sphinx for the user documentation
   - Doxygen for the system documentation

## How to build and install LAMA

extract tar.gz and change into folder
 $ tar -xzvf libama-x.x.x.tar.gz
 $ cd libama-x.x.x

create a build directory and change to it
 $ mkdir <build>
 $ cd <build>

configure cmake by giving the install prefix and pointing to the LAMA-src dir:
 $ cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> [optional options] ../scai

For a Release build be sure to define -DCMAKE_BUILD_TYPE=Release.

start the build and installation process by running make (optionally in parallel):

 $ make [-j <num_processes>]

you optionally can build the system doc (beneath the already installed userdoc) by:

 $ make doxygendoc

make install is not(!) needed (and not possible).
