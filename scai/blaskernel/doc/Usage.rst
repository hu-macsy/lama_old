.. _blaskernel-usage:

Usage
=====

Vector scaling
--------------

To scale a vector blas provides the function scal.
This examples will scale a vector with a factor of 2.
The kernel will be executed on the host.

.. code-block:: c++

    IndexType n = 128;
    
    // Create array with size n
    common::scoped_array<double> x( new double[ n ] );
    
    /*
     ... fill array with data
     */
    
    // Get function from Registry
    kregistry::KernelTraitContextFunction<BLASKernelTrait::scal<double> > scal;

    // Function call on host    
    scal[context::Host]( n, 2.0, x.get(), 1 );
    
    
Matrix times Vector
-------------------

To calculate the matrix-vector multiplication the function gemv is used. 
This function calculates:

.. math::

   y = \alpha * A * x + \beta * y
   
With given input data we just had to get the function from the Registry and execute it 
on the corresponding backend. If we asume that we have pointers to the data on the gpu, 
the call would look like this:

.. code-block:: c++

   // Get function from Registry
   kregistry::KernelTraitContextFunction<BLASKernelTrait::gemv<float> > gemv;
   
   // Function call on host
   gemv[context::CUDA]( order, trans, m, n, alpha, A, lda, x, incX, beta, y, incY );