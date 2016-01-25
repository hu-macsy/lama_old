.. _blaskernel-tested-versions:

***************
Tested Versions
***************

BLAS-libraries
--------------

The following blas libraries has been tested. The blaskernel library contains a own implementation of the blas function
which is called OpenMP-BLAS. 

==============    ================   ===========================================================================   ======================================
BLAS-library      Backend            Supported ValueTypes                                                          Status 
==============    ================   ===========================================================================   ======================================
Fortran-BLAS      Host               float, double, ComplexFloat, ComplexDouble                                    successfully tested
MKL               Host/MIC           float, double, ComplexFloat, ComplexDouble                                    successfully tested
OpenMP-BLAS       Host               float, double, long double, ComplexFloat, ComplexDouble, ComplexLongDouble    successfully tested
cuBLAS            CUDA               float, double, ComplexFloat, ComplexDouble                                    successfully tested
==============    ================   ===========================================================================   ======================================

LAPACK-libraries
----------------

The following lapack libraries has been tested. The blaskernel library contains a own implementation of the lapack function
which is called OpenMP-LAPACK. 

===============    ================   ===========================================================================   ======================================
LAPACK-library      Backend            Supported ValueTypes                                                          Status 
===============    ================   ===========================================================================   ======================================
Fortran-LAPACK      Host               float, double, ComplexFloat, ComplexDouble                                    successfully tested
MKL                 Host               float, double, ComplexFloat, ComplexDouble                                    successfully tested
OpenMP-LAPACK       Host               float, double, long double, ComplexFloat, ComplexDouble, ComplexLongDouble    successfully tested
===============    ================   ===========================================================================   ======================================


Library blaskernel
------------------