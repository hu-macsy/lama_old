Setting a Context
=================

A *Context* is a concept used by LAMA to define where to execute a calculation and therefor where to store the data.
The default context is always the **Host** (CPU). Beneath the host context there exists two other contexts in LAMA:
**CUDA** and **OpenCL** (in progress - coming with the next release), which are located on the GPU (for CUDA) and other
OpenCL supported accelerators (for OpenCL). Further backends for other accelerators are planned for further releases.

Contexts are registered by a factory (see C++ factory pattern) and you only can get an instance of a context by the
factory. If the requested context is not available, because it is not supported on the recent LAMA installation e.g. no
CUDA is found, you receive the default context. So your program using the CUDA context can be compiled and executed on 
another machine without CUDA running on the host without changing code.
 
.. code-block:: c++

   ContextPtr hostCtx = ContextFactor::getContext( Context::Host );
   ContextPtr cudaCtx = ContextFactor::getContext( Context::CUDA, 0 ); // 0 defines the CUDA device used
 
You need a *ContextPtr* to pass it to a matrix and/or vector to set the compute location for calculations with the
matrix/vector. 

.. code-block:: c++

   CSRSparseMatrix<double> csrMatrix( ... );
   csrMatrix.setContext( cudaCtx );
   
   DenseVector<double> vec( ... );
   vec.setContext( cudaCtx );
   
   DenseVector<double> res = 2.0 * vec * csrMatrix;
   
In this example the *csrMatrix* and *vec* are located at the CUDA context. Therefor the sparse matrix times vector
multiplication is executed on the GPU. The new created result vector *res* is also given the CUDA context, because the
calculation has taken place there.

Data movements between different contexts (devices) is done on demand if not declared in another way. So the initial data
of the *csrMatrix* and *vec* are copied to the GPU with the invocation of the multiplication and not with the *setContext*
command. If you want to time the calculation only the data movements have to be attached before by *prefetch* otherwise
the movements is timed, too.

.. code-block:: c++

   CSRSparseMatrix<double> csrMatrix( ... );
   csrMatrix.prefetch( cudaCtx );
   
   DenseVector<double> vec( ... );
   vec.prefetch( cudaCtx );

   // startTime   
   DenseVector<double> res = 2.0 * vec * csrMatrix;
   // endTime
   
As mentioned before the result vector now has also a CUDA context. If you want to perform other calculation with *res*
on the host you can set it to *hostCtx* or prefetch on the host explicitly.

NOTE: the location of the calculation and the context of the result is defined by the the context of the operands and
the rule to copy less data to another context. E.g.: in the given example csrMatrix has a CUDA context but vec has a host
context, the calculation and res is on the GPU; if csrMatrix has a host context and vec a CUDA context, the calculation
and res is on the host.
 