Kernel Functions
================

A kernel function is a function that is executed on a certain device and operates on data that has been allocated on it.

.. code-block:: c++

   void kernel_add( double* a, const double* b, const double* c, int n )

   template<typename T>
   void kernel_add( T* a, const T* b, const T* c, int n );

The routine might be executed either on the CPU host or on a CUDA device or on any other context.
Each pointer argument points to data that has been allocated on the same device where the function 
is called.

A kernel function has the following properties:

- it takes full advantage of available parallelism on the device where they are implemented
  
  - OpenMP on CPU
  
  - CUDA on GPU device

- it can assume that all data is available

- is not synonymous with CUDA kernel routines as here the kernel routine is responsible for launching the corresponding kernel

- it can be implemented by using devie specific libraries like 

  - BLAS, MKL on the CPU
  
  - Thrust library or cuBLAS/cuSPARSE on the GPU

- it should never use communication in kernel routines (possible but makes encapsulation and porting difficult, e.g. SAMG setup, ScaL, ...)

- it might be executed asynchronously

Usually, a kernel function should be available for all supported context types.
