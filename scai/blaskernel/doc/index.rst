.. _main-page:

SCAI BLASKernel 
===============

* Provides wrapper for different BLAS-libraries
* Provides wrapper for different LAPACK-libraries
* Supports different backends (currently Host, CUDA, MIC)

BLAS-libraries
--------------

==============    ================   ===========================================================================   ======================================
BLAS-library      Backend            Supported ValueTypes                                                          Status 
==============    ================   ===========================================================================   ======================================
Fortran-BLAS      Host               float, double, ComplexFloat, ComplexDouble                                    Fix needed: cblas objects not found
MKL               Host/MIC           float, double, ComplexFloat, ComplexDouble                                    successfully tested
OpenMP-BLAS       Host               float, double, long double, ComplexFloat, ComplexDouble, ComplexLongDouble    successfully tested
cuBLAS            CUDA               float, double, ComplexFloat, ComplexDouble                                    successfully tested
==============    ================   ===========================================================================   ======================================

LAPACK-libraries
----------------

===============    ================   ===========================================================================   ======================================
LAPACK-library      Backend            Supported ValueTypes                                                          Status 
===============    ================   ===========================================================================   ======================================
Fortran-LAPACK      Host               float, double, ComplexFloat, ComplexDouble                                    successfully tested
MKL                 Host               float, double, ComplexFloat, ComplexDouble                                    successfully tested
OpenMP-LAPACK       Host               float, double, long double, ComplexFloat, ComplexDouble, ComplexLongDouble    successfully tested
===============    ================   ===========================================================================   ======================================


work in progress
