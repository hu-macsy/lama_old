.. _blaskernel-namingscheme-wrapping:

Wrapping
========

A wrapper consists of two files:

- The trait, which provides information, e.g.
   - for the naming scheme of a implementation
   - used type for indexiation, ordering or transposition
   - for the signatures of calls
- The wrapper
   - uses the BLASKernelTrait to register its functions
   - access the corresponding calls by using the naming macro

The naming scheme for blas and lapack calls differs between implementations. 
Therefore the used naming scheme can be changed by altering the macro in the corresponding Trait. 
One possible naming for blas calls is: a identifier for the used type follew by the functionname and an underscore.
The scal operation with single precision would be named sscal\_. To create this interface the
following macro is used: 

.. code-block:: c++

    #define FORTRAN_BLAS_NAME( name, prefix ) prefix##name##_
    
    FORTRAN_BLAS_NAME( scal, s ) //will be evaluated as --> scal_ 

Currently the wrappers are implemented as a macro to create a specialization of a template. The parameters
can vary between wrappers, e.g. cublas uses its own complex type, thus the pointers need to be casted. This
is possible through their binary compatiblity.

In :ref:`blaskernel-extension` it is shown how to extend the blas-/lapackinterface. 