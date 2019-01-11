# LAMA (Version 3.0.0 Chimborazo)

For building LAMA on your system please consider the installation requirements below and install required packages. 

A description on how to build LAMA is given afterwards or refer to our homepage www.libama.org --> Documentation --> Installation

## Requirements for building LAMA

Required:
 * CMake >= 2.8.8
 * C/C++ compiler (OpenMP 3.0, support of C++11 mandatory)

Recommended:
 * BLAS- and LAPACK-Library (Intel MKL, BLAS)
 * Nvidia GPU: CUDA >= 7.0 (with cuBLAS, cuSPARSE)
 * MPI
 * Boost-Library tests (unit-test-framework), Version 1.61 or higher

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

You can give cmake some hints where to find recommended/optional software either
by setting environment variables or by passing variables to cmake.

  $ export METIS=ROOT=<metis_installation_directory>
  $ export PARMETIS_ROOT=<parmetis_installation_directory>
  $ export BOOST_ROOT=<...>
  $ export MKL_ROOT=<...>
  $ export MPI_HOME=<...>
  $ cmake [ -DMETIS_ROOT=<metis_installation_director>  -DPARMETIS_ROOT=<..> ...] ....

start the build and installation process by running make (optionally in parallel):

 $ make [-j <num_processes>]

If the Boost unit-test framework is available you can run the tests as follows:

 $ make check

you can build the system doc (doxygen documentation) and the user doc (sphinx) as follows:

 $ make doc

You can view the documentation by calling the browser with the start file LAMA.html.

 $ <browser> doc/LAMA.html

If the build was successful, you can install the LAMA software. You have to 
make sure that you have access rights for the installation directory.

 $ [sudo] make install
