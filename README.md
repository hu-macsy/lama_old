# LAMA (Version 2.0.0 Barrancas Blancas)

For building LAMA on your system please consider the installation requirements below and install required packages. 

A description on how to build LAMA is given afterwards or refer to our homepage www.libama.org --> Documentation --> Installation

## Requirements for building LAMA

Required:
 * CMake >= 2.8
 * C/C++ compiler (optionally: OpenMP 2.0, recommendent: capable of C++11)
 * Boost-Library (thread, optional for tests: test and regex )

Recommended:
 * BLAS- and LAPack-Library (Intel MKL, BLAS)
 * CUDA >= 4.0
 * Intel Xeon Phi: Intel MKL
 * MPI
 * GPI-2

Optional:
 * Metis/ParMetis
 * Documentation (optional):
   - doxygen for generating the documentation
   - sphinx for user documentation

## How to build and install LAMA

extract tar.gz and change into folder
 $ tar -xzvf libama-x.x.x.tar.gz
 $ cd libama-x.x.x

create a build directory and change to it
 $ mkdir <build>
 $ cd <build>

configure cmake by giving the install prefix and pointing to the LAMA-src dir:
 $ cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/dir> [optional options] ../scai

For a Release build define -DCMAKE_BUILD_TYPE=Release.

start the build and installation process by running make (optionally in parallel):

 $ make [-j <num_processes>]

you optionally can build the system doc (beneath the already installed userdoc) by:

 $ make doxygendoc
