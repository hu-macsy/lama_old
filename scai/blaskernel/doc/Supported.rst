.. _supported:

Supported Implementations
=========================

The BLASKernel library provides the following implementations:

* The internal implementation provides the supported BLAS and LAPACK routines as own OpenMP 
  versions. The classes ``OpenMPBLAS1``, ``OpenMPBLAS2``, ``OpenMPBLAS3``, and 
  ``OpenMPLAPACK`` contain the corresponding implementations and register these routines
  for the Host context in the kernel registry.

* An external wrapper library provides the supported BLAS and LAPACK routines from an optimized
  vendor library (e.g. Intel MKL, OpenBLAS, ACML, GotoBLAS) for the host device. This implementation
  is rather generic to deal with different naming conventions and with different solutions 
  of parameter passing for the different vendor libraries.

* A CUDA wrapper library provides the supported BLAS and LAPACK routines on a CUDA
  device by wrapping the corresponding routines from the cuBLAS library.

