.. _blaskernel-namingscheme-wrapping:

Wrapping
========

Except from the internal implementation, all other implementations are wrappers as they
wrap the routines from a vendor library.

Each wrapper consists of two files:

- The trait file (e.g. CUBLASTrait.hpp) provides information, e.g.
   - for the naming scheme of a implementation
   - used type for indexiation, ordering or transposition
   - for the signatures of calls
- The wrapper file (e.g. CUBLASWrapper.hpp) 
   - uses the BLASKernelTrait to register its functions
   - access the corresponding calls by using the naming macro

The naming scheme for BLAS and LAPACK calls differs between implementations. 
Therefore the used naming scheme can be changed by altering the macro in the corresponding trait file. 
One possible naming for BLAS calls is: an identifier for the used type followed by the function name and an underscore.
The scal operation with single precision would be named ``sscal_``. To create this interface the
following macro is used: 

.. code-block:: c++

    #define FORTRAN_BLAS_NAME( name, prefix ) prefix##name##_
    
    FORTRAN_BLAS_NAME( scal, s ) //will be evaluated as --> scal_ 

Currently the wrappers are implemented as a macro to create a specialization of a template. The parameters
can vary between wrappers, e.g. cuBLAS uses its own complex type, thus the pointers need to be casted. This
is possible through their binary compatiblity.

In :ref:`blaskernel-extension` it is shown how to extend the blas-/lapackinterface. 
